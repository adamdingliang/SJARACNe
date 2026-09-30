#!/usr/bin/env python3
"""Validate and summarize the completed BRCA100 hub-size DPI pilot."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path

from pilot_common import (
    DRIVERS,
    K_MINIMUM_RECURRENCE,
    SEEDS,
    arm_key,
    atomic_json,
    load_json,
    quartiles,
    sha256_bytes,
    sha256_file,
)


def read_run_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "arm",
        "driver",
        "hub_count",
        "seed",
        "source_commit",
        "binary_sha256",
        "pre_edges",
        "pruned_edges",
        "post_edges",
        "pruned_fraction",
        "sampling_indices",
        "sampling_sha256",
        "adjacency_sha256",
        "data_sha256",
        "anchor_data_match",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Missing required run-manifest fields in {path}")
    return rows


def gate(name: str, passed: bool, detail: str) -> dict[str, object]:
    return {"name": name, "status": "pass" if passed else "fail", "detail": detail}


def analyze(work_root: Path, *, require_anchor: bool) -> tuple[list[dict], dict]:
    design_path = work_root / "pilot_design.json"
    design = load_json(design_path)
    panel_manifest_path = work_root / "panels" / "panel_manifest.json"
    panel_manifest = load_json(panel_manifest_path)
    rows = read_run_manifest(work_root / "results" / "run_manifest.tsv")
    expected_arms = {
        arm_key(driver, count): (driver, count)
        for driver in DRIVERS
        for count in driver.counts
    }
    grouped: dict[str, list[dict[str, str]]] = {arm: [] for arm in expected_arms}
    unexpected: list[str] = []
    for row in rows:
        if row["arm"] not in grouped:
            unexpected.append(row["arm"])
        else:
            grouped[row["arm"]].append(row)

    gates: list[dict[str, object]] = []
    exact_seeds = all(
        sorted(int(row["seed"]) for row in arm_rows) == list(SEEDS)
        for arm_rows in grouped.values()
    )
    gates.append(
        gate(
            "complete_matched_runs",
            not unexpected and len(rows) == len(expected_arms) * len(SEEDS) and exact_seeds,
            f"rows={len(rows)} expected={len(expected_arms) * len(SEEDS)}; "
            "each arm must contain seeds 1..100 exactly once",
        )
    )

    sampling_errors = 0
    sampling_mismatches = 0
    for seed in SEEDS:
        seed_rows = [row for row in rows if int(row["seed"]) == seed]
        parsed: list[tuple[int, ...]] = []
        hashes: set[str] = set()
        for row in seed_rows:
            try:
                indices = tuple(int(value) for value in row["sampling_indices"].split(","))
            except ValueError:
                sampling_errors += 1
                continue
            if (
                len(indices) != 80
                or len(set(indices)) != 80
                or indices != tuple(sorted(indices))
                or any(index < 0 or index >= 100 for index in indices)
            ):
                sampling_errors += 1
            canonical = (" ".join(str(index) for index in indices) + "\n").encode(
                "ascii"
            )
            if row["sampling_sha256"] != sha256_bytes(canonical):
                sampling_errors += 1
            parsed.append(indices)
            hashes.add(row["sampling_sha256"])
        if len(seed_rows) != 6 or len(set(parsed)) != 1 or len(hashes) != 1:
            sampling_mismatches += 1
    gates.append(
        gate(
            "empirically_matched_sampling",
            sampling_errors == 0 and sampling_mismatches == 0 and len(rows) == 600,
            f"invalid_records={sampling_errors}; mismatched_seeds={sampling_mismatches}; "
            "each seed must report the same 80 original indices in all six arms",
        )
    )

    source_commits = {row["source_commit"] for row in rows}
    binary_hashes = {row["binary_sha256"] for row in rows}
    gates.append(
        gate(
            "frozen_source_and_binary",
            len(source_commits) == 1
            and source_commits == {str(design["source_commit"])}
            and len(binary_hashes) == 1
            and binary_hashes == {str(design["build"]["binary_sha256"])}
            and sha256_file(panel_manifest_path) == design["panel_manifest_sha256"],
            "all 600 runs must share the design commit, source snapshot, binary, and panel manifest",
        )
    )

    accounting_errors = 0
    fraction_errors = 0
    for row in rows:
        pre = int(row["pre_edges"])
        pruned = int(row["pruned_edges"])
        post = int(row["post_edges"])
        fraction = float(row["pruned_fraction"])
        accounting_errors += int(pre != pruned + post)
        expected_fraction = 0.0 if pre == 0 else pruned / pre
        fraction_errors += int(not math.isclose(fraction, expected_fraction, abs_tol=1e-15))
    gates.append(
        gate(
            "exact_dpi_accounting",
            accounting_errors == 0 and fraction_errors == 0 and len(rows) == 600,
            f"accounting_errors={accounting_errors}; fraction_errors={fraction_errors}",
        )
    )

    full_rows = [
        row
        for row in rows
        if int(row["hub_count"]) == expected_arms[row["arm"]][0].full_count
    ]
    anchor_evaluated = len(full_rows) == 200 and all(
        row["anchor_data_match"] != "" for row in full_rows
    )
    anchor_pass = anchor_evaluated and all(
        row["anchor_data_match"].lower() == "true" for row in full_rows
    )
    if anchor_evaluated or require_anchor:
        gates.append(
            gate(
                "full_size_anchor",
                anchor_pass,
                f"evaluated={len(full_rows) if anchor_evaluated else 0}/200; "
                "full-size seed adjacency data must match the prior operating-point sweep",
            )
        )
    else:
        gates.append(
            {
                "name": "full_size_anchor",
                "status": "not_evaluated",
                "detail": "rerun full-size arms with --anchor-root before final interpretation",
            }
        )

    summaries: list[dict[str, object]] = []
    k6_ok = True
    k6_details: list[str] = []
    for arm, (driver, count) in sorted(
        expected_arms.items(), key=lambda item: (item[1][0].key, item[1][1])
    ):
        arm_rows = grouped[arm]
        manifest_path = work_root / "results" / arm / "provisional_k6" / "manifest.json"
        output_path = (
            work_root
            / "results"
            / arm
            / "provisional_k6"
            / "consensus_support_ge6.tsv"
        )
        if not manifest_path.is_file() or not output_path.is_file():
            k6_ok = False
            k6_details.append(f"missing {arm}")
            k6 = {
                "k6_edges": None,
                "zero_filled_median_target_count": None,
                "active_hubs": None,
                "active_hub_fraction": None,
            }
        else:
            k6 = load_json(manifest_path)
            valid = (
                int(k6.get("k", -1)) == K_MINIMUM_RECURRENCE
                and int(k6.get("seed_count", -1)) == len(SEEDS)
                and int(k6.get("hub_count", -1)) == count
                and k6.get("output_sha256") == sha256_file(output_path)
            )
            if not valid:
                k6_ok = False
                k6_details.append(f"invalid {arm}")
        if len(arm_rows) == len(SEEDS):
            pre_q1, pre_median, pre_q3 = quartiles(
                int(row["pre_edges"]) for row in arm_rows
            )
            pruned_q1, pruned_median, pruned_q3 = quartiles(
                int(row["pruned_edges"]) for row in arm_rows
            )
            post_q1, post_median, post_q3 = quartiles(
                int(row["post_edges"]) for row in arm_rows
            )
            frac_q1, frac_median, frac_q3 = quartiles(
                float(row["pruned_fraction"]) for row in arm_rows
            )
        else:
            pre_q1 = pre_median = pre_q3 = None
            pruned_q1 = pruned_median = pruned_q3 = None
            post_q1 = post_median = post_q3 = None
            frac_q1 = frac_median = frac_q3 = None
        summaries.append(
            {
                "arm": arm,
                "driver": driver.key,
                "hub_count": count,
                "hub_fraction": count / driver.full_count,
                "completed_seeds": len(arm_rows),
                "pre_edges_q1": pre_q1,
                "pre_edges_median": pre_median,
                "pre_edges_q3": pre_q3,
                "pruned_edges_q1": pruned_q1,
                "pruned_edges_median": pruned_median,
                "pruned_edges_q3": pruned_q3,
                "post_edges_q1": post_q1,
                "post_edges_median": post_median,
                "post_edges_q3": post_q3,
                "pruned_fraction_q1": frac_q1,
                "pruned_fraction_median": frac_median,
                "pruned_fraction_q3": frac_q3,
                "k6_edges": k6["k6_edges"],
                "zero_filled_median_target_count": k6[
                    "zero_filled_median_target_count"
                ],
                "active_hubs": k6["active_hubs"],
                "active_hub_fraction": k6["active_hub_fraction"],
            }
        )
    gates.append(
        gate(
            "independent_k6_aggregation",
            k6_ok and len(summaries) == 6,
            "six benchmark-only direct K>=6 outputs verified"
            + ("" if not k6_details else f"; {'; '.join(k6_details)}"),
        )
    )
    expected_full_k6 = {"tf": 416408, "sig": 739958}
    observed_full_k6 = {
        str(item["driver"]): item["k6_edges"]
        for item in summaries
        if int(item["hub_count"])
        == next(driver.full_count for driver in DRIVERS if driver.key == item["driver"])
    }
    gates.append(
        gate(
            "full_size_k6_anchor",
            observed_full_k6 == expected_full_k6,
            f"expected={expected_full_k6}; observed={observed_full_k6}; "
            "must reproduce the independent prior direct-recurrence sweep",
        )
    )

    trend: dict[str, object] = {}
    for driver in DRIVERS:
        driver_summaries = [item for item in summaries if item["driver"] == driver.key]
        medians = [item["pruned_fraction_median"] for item in driver_summaries]
        if all(value is not None for value in medians):
            numeric = [float(value) for value in medians]
            trend[driver.key] = {
                "hub_counts": [item["hub_count"] for item in driver_summaries],
                "median_pruned_fractions": numeric,
                "nondecreasing": all(
                    right >= left for left, right in zip(numeric, numeric[1:])
                ),
                "full_minus_12p5": numeric[-1] - numeric[0],
                "interpretation": (
                    "descriptive one-panel screen only; no composition-robust or "
                    "Xenium conclusion"
                ),
            }
    validation = {
        "schema": "sjaracne-brca100-hub-size-dpi-validation-v1",
        "design_sha256": sha256_file(design_path),
        "run_manifest_sha256": sha256_file(work_root / "results" / "run_manifest.tsv"),
        "gates": gates,
        "all_required_gates_pass": all(
            item["status"] == "pass" for item in gates if item["status"] != "not_evaluated"
        )
        and (not require_anchor or anchor_pass),
        "trend_screen": trend,
        "limitations": [
            "One deterministic hub panel confounds hub count with panel membership.",
            "The workflow-faithful pilot changes -s and -l together; it does not isolate the DPI mechanism.",
            "Topology and pruning fraction do not establish biological network accuracy.",
            "The full BRCA100 target universe is fixed; results do not validate Xenium 5k.",
            "K>=6 output is an independent benchmark calculation, not production consensus code.",
        ],
    }
    return summaries, validation


def write_summary(path: Path, summaries: list[dict]) -> None:
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    os.replace(temporary, path)


def render_report(summaries: list[dict], validation: dict) -> str:
    lines = [
        "# BRCA100 hub-size DPI pilot",
        "",
        "Primary question: with the whole-transcriptome BRCA100 target universe and all other inference settings fixed, how does hub-list size change DPI pruning?",
        "",
        "| Network | Hubs | Fraction | Median pre-DPI | Median pruned | Median post-DPI | Median pruned fraction (IQR) | K>=6 edges | Median targets, zero-filled | Active hubs |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for item in summaries:
        fraction = item["pruned_fraction_median"]
        if fraction is None:
            fraction_text = "NA"
        else:
            fraction_text = (
                f"{float(fraction):.4f} "
                f"({float(item['pruned_fraction_q1']):.4f}-{float(item['pruned_fraction_q3']):.4f})"
            )
        lines.append(
            "| {driver} | {hub_count:,} | {hub_fraction:.1%} | {pre} | {pruned} | "
            "{post} | {fraction} | {k6} | {targets} | {active} |".format(
                driver=str(item["driver"]).upper(),
                hub_count=int(item["hub_count"]),
                hub_fraction=float(item["hub_fraction"]),
                pre="NA" if item["pre_edges_median"] is None else f"{item['pre_edges_median']:,.0f}",
                pruned="NA" if item["pruned_edges_median"] is None else f"{item['pruned_edges_median']:,.0f}",
                post="NA" if item["post_edges_median"] is None else f"{item['post_edges_median']:,.0f}",
                fraction=fraction_text,
                k6="NA" if item["k6_edges"] is None else f"{int(item['k6_edges']):,}",
                targets=(
                    "NA"
                    if item["zero_filled_median_target_count"] is None
                    else f"{float(item['zero_filled_median_target_count']):.1f}"
                ),
                active=(
                    "NA"
                    if item["active_hub_fraction"] is None
                    else f"{float(item['active_hub_fraction']):.1%}"
                ),
            )
        )
    lines.extend(["", "## Validation gates", ""])
    for item in validation["gates"]:
        lines.append(f"- **{item['status'].upper()} -- {item['name']}**: {item['detail']}")
    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "This is a one-panel, workflow-faithful screen. A rising pruning fraction would justify the fixed-annotation control and replicated panels; it would not by itself identify a biologically optimal hub count or validate Xenium 5k.",
            "",
            "The K>=6 results are independently counted benchmark outputs and must not be represented as the production minimum-recurrence implementation.",
        ]
    )
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--work-root",
        type=Path,
        default=Path.home() / "sjaracne-benchmarks" / "brca100-hub-size-dpi-pilot",
    )
    parser.add_argument(
        "--require-anchor",
        action="store_true",
        help="Fail the full-size anchor gate when --anchor-root was not used",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    summaries, validation = analyze(args.work_root, require_anchor=args.require_anchor)
    results_root = args.work_root / "results"
    write_summary(results_root / "dpi_summary.tsv", summaries)
    atomic_json(results_root / "validation_gates.json", validation)
    report = render_report(summaries, validation)
    temporary = results_root / "REPORT.md.partial"
    temporary.write_text(report, encoding="utf-8", newline="\n")
    os.replace(temporary, results_root / "REPORT.md")
    print(results_root / "REPORT.md")
    return 0 if validation["all_required_gates_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
