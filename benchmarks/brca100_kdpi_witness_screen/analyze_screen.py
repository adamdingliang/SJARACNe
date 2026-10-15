#!/usr/bin/env python3
"""Analyze the completed BRCA100 SIG K_DPI witness-sidecar screen."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import sys

import run_screen
from screen_common import (
    BASELINE_DEFAULT,
    HUB_COUNTS,
    K_VALUES,
    SEEDS,
    WORK_ROOT_DEFAULT,
    aggregate_sidecar_group,
    arm_key,
    atomic_json,
    expression_index,
    load_json,
    panel_source_indices,
    parse_witness_sidecar,
    quartiles,
    read_unique_ids,
    result_root,
    sha256_file,
)


SOURCE_GROUPS = ("all_sources", "common_1335_sources")


def gate(name: str, passed: bool, detail: str) -> dict[str, object]:
    return {"name": name, "status": "pass" if passed else "fail", "detail": detail}


def compare_common_source_pre_edges(
    pre_edges_by_panel_seed: dict[tuple[int, int], dict[int, int]],
    common_source_indices: set[int],
    *,
    hub_counts: tuple[int, ...] = HUB_COUNTS,
    seeds: tuple[int, ...] = SEEDS,
) -> dict[str, int]:
    """Compare every common source row to the small-panel pre-DPI reference."""

    expected = len(common_source_indices) * (len(hub_counts) - 1) * len(seeds)
    comparisons = 0
    mismatches = 0
    missing_panel_pairs = 0
    small = hub_counts[0]
    for seed in seeds:
        reference = pre_edges_by_panel_seed.get((small, seed))
        for count in hub_counts[1:]:
            candidate = pre_edges_by_panel_seed.get((count, seed))
            if reference is None or candidate is None:
                missing_panel_pairs += 1
                continue
            for source_index in common_source_indices:
                if source_index not in reference or source_index not in candidate:
                    mismatches += 1
                    continue
                comparisons += 1
                mismatches += int(reference[source_index] != candidate[source_index])
    return {
        "common_pre_edge_expected_comparisons": expected,
        "common_pre_edge_comparisons": comparisons,
        "common_pre_edge_mismatches": mismatches,
        "common_pre_edge_missing_panel_pairs": missing_panel_pairs,
    }


def marker_paths(work_root: Path, count: int, seed: int) -> dict[str, Path]:
    root = result_root(work_root, count)
    stem = f"TF_run_{seed:03d}"
    return {
        "marker": root / "seed_metadata" / f"{stem}.json",
        "adjacency": root / "adjacency" / f"{stem}.adj",
        "sidecar": root / "witness_sidecars" / f"{stem}.tsv",
        "stdout": root / "logs" / f"{stem}.stdout.log",
    }


def collect_seed_metrics(
    work_root: Path, baseline_root: Path, design: dict
) -> tuple[list[dict], list[str], dict[str, int]]:
    baseline = run_screen.inspect_baseline(baseline_root)
    expected_baseline = design["baseline_import"]
    for key in ("source_commit", "source_tree_fingerprint", "config_sha256"):
        if baseline[key] != expected_baseline[key]:
            raise RuntimeError(f"Baseline {key} differs from frozen screen design")
    for key in ("pilot_design", "panel_manifest", "expression", "null_model"):
        if baseline[key]["sha256"] != expected_baseline[key]["sha256"]:
            raise RuntimeError(f"Baseline {key} hash differs from frozen screen design")

    expression_mapping, all_expression_ids = expression_index(
        Path(baseline["expression"]["path"])
    )
    panel_ids: dict[int, set[str]] = {}
    source_indices: dict[int, set[int]] = {}
    for count in HUB_COUNTS:
        ids = set(read_unique_ids(Path(baseline["panels"][str(count)]["path"]), count))
        panel_ids[count] = ids
        source_indices[count] = panel_source_indices(expression_mapping, ids)
    common = source_indices[HUB_COUNTS[0]]
    nested_failures = sum(
        not common.issubset(source_indices[count]) for count in HUB_COUNTS
    )

    metrics: list[dict] = []
    errors: list[str] = []
    counters = {
        "completed": 0,
        "k1_reproduction_failures": 0,
        "sidecar_accounting_failures": 0,
        "provenance_failures": 0,
        "nested_failures": nested_failures,
        "monotonic_failures": 0,
    }
    common_pre_edges: dict[tuple[int, int], dict[int, int]] = {}
    candidate = design["candidate"]
    design_sha256 = sha256_file(work_root / "screen_design.json")
    for count in HUB_COUNTS:
        for seed in SEEDS:
            paths = marker_paths(work_root, count, seed)
            executed_output = (
                result_root(work_root, count)
                / "work"
                / f"TF_run_{seed:03d}.adj.partial"
            )
            expected_sidecar_paths = run_screen.expected_sidecar_provenance(
                baseline=baseline, count=count, executed_output=executed_output
            )
            missing = [name for name, path in paths.items() if not path.is_file()]
            if missing:
                errors.append(f"{arm_key(count)}/{seed:03d}: missing {missing}")
                continue
            try:
                marker = load_json(paths["marker"])
                if (
                    marker.get("schema") != run_screen.SEED_SCHEMA
                    or marker.get("design_sha256") != design_sha256
                    or marker.get("driver") != "sig"
                    or int(marker.get("hub_count", -1)) != count
                    or int(marker.get("seed", -1)) != seed
                ):
                    raise ValueError("seed marker identity/design mismatch")
                provenance_ok = (
                    marker.get("candidate_source_tree_fingerprint")
                    == candidate["source_tree_fingerprint"]
                    and marker.get("candidate_binary_sha256") == candidate["binary_sha256"]
                    and marker.get("config_sha256") == candidate["config_sha256"]
                    and marker.get("null_model_sha256") == candidate["null_model_sha256"]
                )
                counters["provenance_failures"] += int(not provenance_ok)
                baseline_marker, _ = run_screen.load_baseline_seed(
                    baseline_root, count, seed
                )
                adjacency, dpi, sampling, sidecar_summary, comparison = (
                    run_screen.validate_candidate(
                        output=paths["adjacency"],
                        sidecar_path=paths["sidecar"],
                        stdout_path=paths["stdout"],
                        panel_ids=panel_ids[count],
                        all_expression_ids=all_expression_ids,
                        source_indices=source_indices[count],
                        baseline_marker=baseline_marker,
                        expected_sidecar_paths=expected_sidecar_paths,
                    )
                )
                reproduction_ok = all(
                    comparison[key]
                    for key in (
                        "sampling_match",
                        "dpi_stats_match",
                        "adjacency_data_match",
                        "adjacency_edge_count_match",
                    )
                )
                counters["k1_reproduction_failures"] += int(not reproduction_ok)
                if marker.get("adjacency") != adjacency or marker.get("dpi") != dpi or marker.get("sampling") != sampling:
                    raise ValueError("raw candidate outputs differ from seed marker")
                if marker.get("sidecar") != sidecar_summary:
                    raise ValueError("raw sidecar differs from seed marker")
                parsed = parse_witness_sidecar(
                    paths["sidecar"],
                    expected_source_indices=source_indices[count],
                    expected_provenance=expected_sidecar_paths,
                )
                totals = parsed["totals"]
                accounting_ok = (
                    int(totals["pre_edges"]) == int(dpi["pre_edges"])
                    and int(totals["witnesses_ge_1"]) == int(dpi["pruned_edges"])
                )
                counters["sidecar_accounting_failures"] += int(not accounting_ok)
                common_pre_edges[(count, seed)] = {
                    source_index: int(parsed["rows"][source_index]["pre_edges"])
                    for source_index in common
                }

                group_indices = {
                    "all_sources": source_indices[count],
                    "common_1335_sources": common,
                }
                for group, indices in group_indices.items():
                    by_k = {
                        k: aggregate_sidecar_group(parsed, indices, k) for k in K_VALUES
                    }
                    pruned_counts = [int(by_k[k]["pruned_edges"]) for k in K_VALUES]
                    counters["monotonic_failures"] += int(
                        any(right > left for left, right in zip(pruned_counts, pruned_counts[1:]))
                    )
                    k1_pruned = pruned_counts[0]
                    for k in K_VALUES:
                        item = by_k[k]
                        metrics.append(
                            {
                                "source_group": group,
                                "hub_count": count,
                                "seed": seed,
                                "k_dpi": k,
                                "source_rows": len(indices),
                                "pre_edges": item["pre_edges"],
                                "pruned_edges": item["pruned_edges"],
                                "post_edges": item["post_edges"],
                                "pruned_fraction": item["pruned_fraction"],
                                "pruning_retained_vs_k1": (
                                    None
                                    if k1_pruned == 0
                                    else int(item["pruned_edges"]) / k1_pruned
                                ),
                            }
                        )
                counters["completed"] += 1
            except Exception as error:
                errors.append(f"{arm_key(count)}/{seed:03d}: {error}")
    counters.update(
        compare_common_source_pre_edges(common_pre_edges, common)
    )
    return metrics, errors, counters


def summarize_groups(metrics: list[dict]) -> list[dict]:
    output: list[dict] = []
    for group in SOURCE_GROUPS:
        for count in HUB_COUNTS:
            for k in K_VALUES:
                rows = [
                    row
                    for row in metrics
                    if row["source_group"] == group
                    and row["hub_count"] == count
                    and row["k_dpi"] == k
                ]

                def qs(field: str) -> tuple[float | None, float | None, float | None]:
                    values = [row[field] for row in rows if row[field] is not None]
                    return quartiles(values) if len(values) >= 2 else (None, None, None)

                pre = qs("pre_edges")
                pruned = qs("pruned_edges")
                fraction = qs("pruned_fraction")
                retained = qs("pruning_retained_vs_k1")
                output.append(
                    {
                        "source_group": group,
                        "hub_count": count,
                        "k_dpi": k,
                        "completed_seeds": len(rows),
                        "source_rows": None if not rows else rows[0]["source_rows"],
                        "pre_edges_q1": pre[0],
                        "pre_edges_median": pre[1],
                        "pre_edges_q3": pre[2],
                        "pruned_edges_q1": pruned[0],
                        "pruned_edges_median": pruned[1],
                        "pruned_edges_q3": pruned[2],
                        "pruned_fraction_q1": fraction[0],
                        "pruned_fraction_median": fraction[1],
                        "pruned_fraction_q3": fraction[2],
                        "pruning_retained_vs_k1_q1": retained[0],
                        "pruning_retained_vs_k1_median": retained[1],
                        "pruning_retained_vs_k1_q3": retained[2],
                    }
                )
    return output


def paired_effects(metrics: list[dict]) -> list[dict]:
    lookup = {
        (row["source_group"], row["hub_count"], row["seed"], row["k_dpi"]): row
        for row in metrics
    }
    output: list[dict] = []
    small, full = HUB_COUNTS[0], HUB_COUNTS[-1]
    span_thousands = (full - small) / 1000.0
    for group in SOURCE_GROUPS:
        for k in K_VALUES:
            seed_rows: list[dict] = []
            for seed in SEEDS:
                small_row = lookup.get((group, small, seed, k))
                full_row = lookup.get((group, full, seed, k))
                small_k1 = lookup.get((group, small, seed, 1))
                full_k1 = lookup.get((group, full, seed, 1))
                if any(item is None for item in (small_row, full_row, small_k1, full_k1)):
                    continue
                small_fraction = float(small_row["pruned_fraction"])
                full_fraction = float(full_row["pruned_fraction"])
                delta = full_fraction - small_fraction
                delta_k1 = float(full_k1["pruned_fraction"]) - float(
                    small_k1["pruned_fraction"]
                )
                seed_rows.append(
                    {
                        "delta": delta,
                        "absolute_delta": abs(delta),
                        "slope": delta / span_thousands,
                        "absolute_slope": abs(delta) / span_thousands,
                        "ratio": None if small_fraction == 0.0 else full_fraction / small_fraction,
                        "normalized_gap": None if delta_k1 == 0.0 else delta / delta_k1,
                        "absolute_normalized_gap": (
                            None if delta_k1 == 0.0 else abs(delta) / abs(delta_k1)
                        ),
                        "full_gt_small": delta > 0.0,
                        "full_lt_small": delta < 0.0,
                        "gap_reduced": abs(delta) < abs(delta_k1),
                        "same_direction": delta == 0.0 or delta_k1 == 0.0 or (delta > 0) == (delta_k1 > 0),
                        "full_retained": full_row["pruning_retained_vs_k1"],
                        "small_retained": small_row["pruning_retained_vs_k1"],
                    }
                )

            def qs(field: str) -> tuple[float | None, float | None, float | None]:
                values = [row[field] for row in seed_rows if row[field] is not None]
                return quartiles(values) if len(values) >= 2 else (None, None, None)

            delta_q = qs("delta")
            absolute_delta_q = qs("absolute_delta")
            slope_q = qs("slope")
            absolute_slope_q = qs("absolute_slope")
            ratio_q = qs("ratio")
            normalized_q = qs("normalized_gap")
            absolute_normalized_q = qs("absolute_normalized_gap")
            full_retained_q = qs("full_retained")
            small_retained_q = qs("small_retained")
            output.append(
                {
                    "source_group": group,
                    "k_dpi": k,
                    "paired_seeds": len(seed_rows),
                    "delta_pruned_fraction_q1": delta_q[0],
                    "delta_pruned_fraction_median": delta_q[1],
                    "delta_pruned_fraction_q3": delta_q[2],
                    "absolute_delta_pruned_fraction_q1": absolute_delta_q[0],
                    "absolute_delta_pruned_fraction_median": absolute_delta_q[1],
                    "absolute_delta_pruned_fraction_q3": absolute_delta_q[2],
                    "slope_per_1000_hubs_q1": slope_q[0],
                    "slope_per_1000_hubs_median": slope_q[1],
                    "slope_per_1000_hubs_q3": slope_q[2],
                    "absolute_slope_per_1000_hubs_q1": absolute_slope_q[0],
                    "absolute_slope_per_1000_hubs_median": absolute_slope_q[1],
                    "absolute_slope_per_1000_hubs_q3": absolute_slope_q[2],
                    "full_over_small_ratio_q1": ratio_q[0],
                    "full_over_small_ratio_median": ratio_q[1],
                    "full_over_small_ratio_q3": ratio_q[2],
                    "normalized_gap_vs_k1_q1": normalized_q[0],
                    "normalized_gap_vs_k1_median": normalized_q[1],
                    "normalized_gap_vs_k1_q3": normalized_q[2],
                    "absolute_normalized_gap_vs_k1_q1": absolute_normalized_q[0],
                    "absolute_normalized_gap_vs_k1_median": absolute_normalized_q[1],
                    "absolute_normalized_gap_vs_k1_q3": absolute_normalized_q[2],
                    "small_pruning_retained_vs_k1_median": small_retained_q[1],
                    "full_pruning_retained_vs_k1_median": full_retained_q[1],
                    "full_gt_small_seeds": sum(row["full_gt_small"] for row in seed_rows),
                    "full_lt_small_seeds": sum(row["full_lt_small"] for row in seed_rows),
                    "equal_seeds": sum(
                        not row["full_gt_small"] and not row["full_lt_small"]
                        for row in seed_rows
                    ),
                    "abs_gap_reduced_vs_k1_seeds": sum(row["gap_reduced"] for row in seed_rows),
                    "same_direction_as_k1_seeds": sum(row["same_direction"] for row in seed_rows),
                }
            )
    return output


def validate(
    metrics: list[dict], errors: list[str], counters: dict[str, int], paired: list[dict]
) -> dict[str, object]:
    expected_runs = len(HUB_COUNTS) * len(SEEDS)
    expected_metrics = expected_runs * len(SOURCE_GROUPS) * len(K_VALUES)
    gates = [
        gate(
            "complete_matched_runs",
            counters["completed"] == expected_runs and not errors and len(metrics) == expected_metrics,
            f"runs={counters['completed']}/{expected_runs}; metrics={len(metrics)}/{expected_metrics}; errors={len(errors)}",
        ),
        gate(
            "exact_k1_reproduction",
            counters["k1_reproduction_failures"] == 0 and counters["completed"] == expected_runs,
            f"sample/DPI/adjacency-data failures={counters['k1_reproduction_failures']}",
        ),
        gate(
            "exact_sidecar_accounting",
            counters["sidecar_accounting_failures"] == 0 and counters["completed"] == expected_runs,
            f"sum(pre) or sum(witnesses_ge_1) failures={counters['sidecar_accounting_failures']}",
        ),
        gate(
            "frozen_candidate_provenance",
            counters["provenance_failures"] == 0 and counters["completed"] == expected_runs,
            f"marker provenance failures={counters['provenance_failures']}",
        ),
        gate(
            "nested_common_source_coverage",
            counters["nested_failures"] == 0,
            f"nested-panel failures={counters['nested_failures']}",
        ),
        gate(
            "identical_common_source_pre_dpi_edge_counts",
            counters["common_pre_edge_comparisons"]
            == counters["common_pre_edge_expected_comparisons"]
            and counters["common_pre_edge_mismatches"] == 0
            and counters["common_pre_edge_missing_panel_pairs"] == 0,
            "source-row comparisons="
            f"{counters['common_pre_edge_comparisons']}/"
            f"{counters['common_pre_edge_expected_comparisons']}; "
            f"mismatches={counters['common_pre_edge_mismatches']}; "
            "missing panel/seed pairs="
            f"{counters['common_pre_edge_missing_panel_pairs']}",
        ),
        gate(
            "per_seed_pruning_monotone_with_kdpi",
            counters["monotonic_failures"] == 0 and counters["completed"] == expected_runs,
            f"source-group trajectories with violations={counters['monotonic_failures']}",
        ),
    ]
    lookup = {
        (row["source_group"], row["hub_count"], row["seed"], row["k_dpi"]): row
        for row in metrics
    }
    equality_failures = 0
    for seed in SEEDS:
        for k in K_VALUES:
            all_row = lookup.get(("all_sources", HUB_COUNTS[0], seed, k))
            common_row = lookup.get(("common_1335_sources", HUB_COUNTS[0], seed, k))
            if all_row is None or common_row is None or all_row != {
                **common_row,
                "source_group": "all_sources",
            }:
                equality_failures += 1
    gates.append(
        gate(
            "small_panel_all_equals_common",
            equality_failures == 0,
            f"mismatched seed/K rows={equality_failures}",
        )
    )
    paired_complete = all(row["paired_seeds"] == len(SEEDS) for row in paired)
    gates.append(
        gate(
            "complete_paired_full_minus_small_effects",
            paired_complete and len(paired) == len(SOURCE_GROUPS) * len(K_VALUES),
            f"paired rows={len(paired)}/{len(SOURCE_GROUPS) * len(K_VALUES)}",
        )
    )
    return {
        "schema": "sjaracne-brca100-kdpi-witness-validation-v1",
        "all_required_gates_pass": all(item["status"] == "pass" for item in gates),
        "gates": gates,
        "collection_errors": errors,
        "limitations": [
            "Ten seeds and one deterministic nested SIG panel trajectory are only a screening design.",
            "The sidecar infers K_DPI outputs from exact per-edge witness multiplicities; it does not emit separate adjacency files for K_DPI>1.",
            "No K_edge consensus aggregation or downstream activity robustness is evaluated.",
            "The screen does not establish biological correctness for removed or retained edges.",
            "A reduced size gap is not useful if it is achieved by making DPI nearly inert.",
        ],
    }


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    fields = list(rows[0]) if rows else ["empty"]
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def fmt(value: object, digits: int = 4) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    return str(value)


def render_report(
    summaries: list[dict], paired: list[dict], validation: dict[str, object]
) -> str:
    lines = [
        "# BRCA100 SIG `K_DPI` witness screen",
        "",
        "This 10-seed matched screen changes only the number of qualifying DPI witnesses "
        "required to remove an edge. It does not retune `p_b`, AP-MI, epsilon, panels, or samples.",
        "",
        "## Primary full-minus-small result",
        "",
        "`common_1335_sources` evaluates the same 1,335 source rows at every H and is the primary "
        "technical-bias readout. `all_sources` changes the evaluated source population with H.",
        "",
        "| source group | K_DPI | paired seeds | median signed full-small gap | median absolute gap | absolute slope / 1,000 hubs | full/small ratio | normalized absolute gap vs K=1 | small-panel pruning retained | full-panel pruning retained | seeds gap reduced |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in paired:
        lines.append(
            f"| {row['source_group']} | {row['k_dpi']} | {row['paired_seeds']} | "
            f"{fmt(row['delta_pruned_fraction_median'], 5)} | "
            f"{fmt(row['absolute_delta_pruned_fraction_median'], 5)} | "
            f"{fmt(row['absolute_slope_per_1000_hubs_median'], 5)} | "
            f"{fmt(row['full_over_small_ratio_median'], 5)} | "
            f"{fmt(row['absolute_normalized_gap_vs_k1_median'], 5)} | "
            f"{fmt(row['small_pruning_retained_vs_k1_median'], 5)} | "
            f"{fmt(row['full_pruning_retained_vs_k1_median'], 5)} | "
            f"{row['abs_gap_reduced_vs_k1_seeds']}/{row['paired_seeds']} |"
        )
    lines.extend(
        [
            "",
            "## Pruning trajectories",
            "",
            "The pruning fraction is the number of source-row pre-DPI edges with at least K "
            "qualifying witnesses divided by all pre-DPI edges in the stated source population.",
            "",
            "| source group | H | K_DPI | seeds | median pruning fraction [IQR] | median pruning retained vs K=1 |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in summaries:
        if row["pruned_fraction_median"] is None:
            fraction = "NA"
        else:
            fraction = (
                f"{row['pruned_fraction_median']:.4f} "
                f"[{row['pruned_fraction_q1']:.4f}, {row['pruned_fraction_q3']:.4f}]"
            )
        lines.append(
            f"| {row['source_group']} | {row['hub_count']:,} | {row['k_dpi']} | "
            f"{row['completed_seeds']} | {fraction} | "
            f"{fmt(row['pruning_retained_vs_k1_median'], 5)} |"
        )
    lines.extend(["", "## Validation", ""])
    for item in validation["gates"]:
        lines.append(f"- **{item['status'].upper()} -- {item['name']}**: {item['detail']}")
    lines.extend(
        [
            "",
            "## Defensible interpretation",
            "",
            "A higher K_DPI is promising only if it materially lowers the common-source paired "
            "gap and does so consistently across seeds while retaining a nontrivial fraction of "
            "K=1 pruning. A flatter trajectory caused by eliminating nearly all DPI pruning is "
            "not a useful correction. This screen cannot choose a biological default.",
            "",
        ]
    )
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-root", type=Path, default=WORK_ROOT_DEFAULT)
    parser.add_argument("--baseline-root", type=Path, default=BASELINE_DEFAULT)
    parser.add_argument(
        "--require-complete",
        action="store_true",
        help="Return nonzero unless every validation gate passes",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    design_path = args.work_root / "screen_design.json"
    if not design_path.is_file():
        raise FileNotFoundError(f"Run prepare/infer first: {design_path}")
    design = load_json(design_path)
    if design.get("schema") != run_screen.SCHEMA:
        raise RuntimeError("Unexpected screen design schema")
    run_screen.verify_frozen_harness_hashes(design.get("harness_sha256", {}))
    metrics, errors, counters = collect_seed_metrics(
        args.work_root, args.baseline_root, design
    )
    summaries = summarize_groups(metrics)
    paired = paired_effects(metrics)
    validation = validate(metrics, errors, counters, paired)
    results = args.work_root / "results"
    write_tsv(results / "seed_level_metrics.tsv", metrics)
    write_tsv(results / "source_group_summary.tsv", summaries)
    write_tsv(results / "paired_hub_size_effects.tsv", paired)
    atomic_json(results / "validation_gates.json", validation)
    report = render_report(summaries, paired, validation)
    temporary = results / "REPORT.md.partial"
    temporary.write_text(report, encoding="utf-8", newline="\n")
    os.replace(temporary, results / "REPORT.md")
    print(f"[ANALYZE] wrote {results / 'REPORT.md'}")
    if args.require_complete and not validation["all_required_gates_pass"]:
        return 2
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
