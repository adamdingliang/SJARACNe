#!/usr/bin/env python3
"""Run the matched BRCA100 SIG K_DPI witness-sidecar screen."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from screen_common import (
    BASELINE_DEFAULT,
    DRIVER_FILENAME,
    EXPECTED_BASELINE,
    EXPRESSION_FILENAME,
    HUB_COUNTS,
    K_VALUES,
    MI_CUTOFF,
    N_PAR,
    P_TOKEN,
    P_VALUE,
    REPO_DEFAULT,
    SEEDS,
    WORK_ROOT_DEFAULT,
    arm_key,
    atomic_json,
    baseline_arm_root,
    ensure_exact_json,
    expression_index,
    fingerprint,
    load_json,
    panel_path,
    panel_source_indices,
    parse_dpi_stats,
    parse_sampling,
    parse_witness_sidecar,
    read_unique_ids,
    result_root,
    sha256_file,
    utc_now,
    validate_adjacency,
)


THIS_DIR = Path(__file__).resolve().parent
PILOT_HARNESS_DIR = THIS_DIR.parent / "brca100_hub_size_dpi_pilot"
if str(PILOT_HARNESS_DIR) not in sys.path:
    sys.path.insert(0, str(PILOT_HARNESS_DIR))
import run_pilot as pilot_build  # noqa: E402


SCHEMA = "sjaracne-brca100-kdpi-witness-screen-v1"
SEED_SCHEMA = "sjaracne-brca100-kdpi-witness-seed-v1"


def harness_file_paths() -> dict[str, Path]:
    paths = (
        THIS_DIR / "screen_common.py",
        THIS_DIR / "run_screen.py",
        THIS_DIR / "analyze_screen.py",
        PILOT_HARNESS_DIR / "pilot_common.py",
        PILOT_HARNESS_DIR / "run_pilot.py",
    )
    return {str(path.relative_to(THIS_DIR.parent)): path for path in paths}


def harness_hashes(paths: dict[str, Path] | None = None) -> dict[str, str]:
    selected = harness_file_paths() if paths is None else paths
    missing = [label for label, path in selected.items() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing harness files: {missing}")
    return {label: sha256_file(path) for label, path in selected.items()}


def verify_frozen_harness_hashes(
    expected: dict[str, str], paths: dict[str, Path] | None = None
) -> dict[str, str]:
    observed = harness_hashes(paths)
    if observed != expected:
        missing = sorted(set(expected) - set(observed))
        extra = sorted(set(observed) - set(expected))
        changed = sorted(
            label
            for label in set(expected) & set(observed)
            if expected[label] != observed[label]
        )
        raise RuntimeError(
            "Frozen harness hashes do not match current files: "
            f"missing={missing}; extra={extra}; changed={changed}"
        )
    return observed


def selected_counts(specification: str) -> list[int]:
    if specification == "all":
        return list(HUB_COUNTS)
    requested = [int(item.strip()) for item in specification.split(",") if item.strip()]
    if (
        not requested
        or len(requested) != len(set(requested))
        or set(requested) - set(HUB_COUNTS)
    ):
        raise ValueError(f"--hub-counts must be unique values from {HUB_COUNTS}")
    return requested


def parse_gnu_time(path: Path) -> dict[str, float | int]:
    values: dict[str, float | int] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        key, value = line.split("=", 1)
        values[key] = int(value) if key == "max_rss_kib" else float(value)
    required = {"elapsed_s", "user_s", "system_s", "max_rss_kib"}
    if set(values) != required:
        raise ValueError(f"Incomplete GNU time record in {path}")
    return values


def checked_run(
    command: list[str], *, cwd: Path, stdout_path: Path, stderr_path: Path
) -> None:
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        subprocess.run(command, cwd=cwd, stdout=stdout, stderr=stderr, check=True)


def _verified(path: Path, expected_sha256: str, label: str) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing frozen {label}: {path}")
    observed = sha256_file(path)
    if observed != expected_sha256:
        raise RuntimeError(
            f"Frozen {label} hash mismatch: expected {expected_sha256}, got {observed}"
        )
    return {"path": str(path.resolve()), "sha256": observed, "bytes": path.stat().st_size}


def inspect_baseline(baseline_root: Path) -> dict[str, object]:
    """Resolve and hash every frozen pilot artifact reused by this screen."""

    design_path = baseline_root / "pilot_design.json"
    panel_manifest_path = baseline_root / "panels" / "panel_manifest.json"
    design = load_json(design_path)
    panel_manifest = load_json(panel_manifest_path)
    if design.get("schema") != EXPECTED_BASELINE["schema"]:
        raise RuntimeError("Unexpected baseline pilot schema")
    for key in ("source_commit", "source_tree_fingerprint"):
        if design.get(key) != EXPECTED_BASELINE[key]:
            raise RuntimeError(f"Unexpected frozen baseline {key}")
    build_stage = str(design["build"]["stage"])
    build_manifest_path = baseline_root / "builds" / build_stage / "build_manifest.json"
    build = load_json(build_manifest_path)
    for key in ("binary_sha256", "config_sha256", "null_model_sha256"):
        if design["build"].get(key) != EXPECTED_BASELINE[key] or build.get(key) != EXPECTED_BASELINE[key]:
            raise RuntimeError(f"Unexpected frozen baseline {key}")

    expression = baseline_root / "inputs" / EXPRESSION_FILENAME
    artifacts: dict[str, object] = {
        "baseline_root": str(baseline_root.resolve()),
        "pilot_design": _verified(design_path, sha256_file(design_path), "pilot design"),
        "panel_manifest": _verified(
            panel_manifest_path, sha256_file(panel_manifest_path), "panel manifest"
        ),
        "build_manifest": _verified(
            build_manifest_path, sha256_file(build_manifest_path), "build manifest"
        ),
        "source_commit": design["source_commit"],
        "source_tree_fingerprint": design["source_tree_fingerprint"],
        "binary": _verified(
            Path(build["binary"]), EXPECTED_BASELINE["binary_sha256"], "binary"
        ),
        "config_directory": str(Path(build["config_directory"]).resolve()),
        "config_sha256": build["config_sha256"],
        "null_model": _verified(
            Path(build["null_model"]), EXPECTED_BASELINE["null_model_sha256"], "null model"
        ),
        "expression": _verified(
            expression, EXPECTED_BASELINE["expression_sha256"], "expression matrix"
        ),
        "panels": {},
    }
    previous_ids: set[str] = set()
    for count in HUB_COUNTS:
        path = panel_path(baseline_root, count)
        ids = set(read_unique_ids(path, count))
        if not previous_ids.issubset(ids):
            raise RuntimeError(f"Frozen SIG panel is not nested at H={count}")
        previous_ids = ids
        panel_record = next(
            record
            for record in panel_manifest["drivers"]["sig"]["panels"]
            if int(record["hub_count"]) == count
        )
        expected = EXPECTED_BASELINE["panel_sha256"][count]
        if panel_record["sha256"] != expected:
            raise RuntimeError(f"Baseline panel manifest drift for H={count}")
        artifacts["panels"][str(count)] = _verified(path, expected, f"SIG panel H={count}")
    return artifacts


def prepare(
    *, repo: Path, work_root: Path, baseline_root: Path, source_commit: str
) -> tuple[dict, dict, dict]:
    baseline = inspect_baseline(baseline_root)
    resolved_commit = pilot_build.resolve_commit(repo, source_commit)
    if resolved_commit != EXPECTED_BASELINE["source_commit"]:
        raise RuntimeError(
            "Candidate snapshot must be overlaid on the frozen pilot commit "
            f"{EXPECTED_BASELINE['source_commit']}; got {resolved_commit}"
        )
    source_files = list((repo / "SJARACNe" / "src").glob("*.cpp")) + list(
        (repo / "SJARACNe" / "src").glob("*.h")
    )
    if not any(b"[DPI_WITNESS]" in path.read_bytes() for path in source_files):
        raise RuntimeError(
            "Candidate source tree lacks the required -W/[DPI_WITNESS] diagnostic interface"
        )
    candidate = pilot_build.extract_or_snapshot_and_build(repo, work_root, resolved_commit)
    if candidate.get("commit") != resolved_commit:
        raise RuntimeError("Candidate build manifest commit mismatch")
    for key in ("config_sha256", "null_model_sha256"):
        if candidate.get(key) != EXPECTED_BASELINE[key]:
            raise RuntimeError(
                f"Candidate {key} differs from the frozen baseline; this is not a matched screen"
            )

    harness = harness_hashes()
    design = {
        "schema": SCHEMA,
        "question": (
            "Does requiring multiple qualifying DPI witnesses reduce SIG hub-list-size-dependent "
            "pruning without merely disabling most DPI removals?"
        ),
        "baseline_import": baseline,
        "candidate": {
            "base_commit": resolved_commit,
            "source_mode": candidate["source_mode"],
            "source_tree_fingerprint": candidate["source_provenance"]["source_tree_fingerprint"],
            "source_provenance": candidate["source_provenance"],
            "build_stage": candidate["stage"],
            "build_manifest": str(
                work_root / "builds" / candidate["stage"] / "build_manifest.json"
            ),
            "binary": candidate["binary"],
            "binary_sha256": candidate["binary_sha256"],
            "config_directory": candidate["config_directory"],
            "config_sha256": candidate["config_sha256"],
            "null_model": candidate["null_model"],
            "null_model_sha256": candidate["null_model_sha256"],
            "compiler": candidate["compiler"],
        },
        "harness_sha256": harness,
        "grid": {
            "driver": "sig",
            "hub_counts": list(HUB_COUNTS),
            "seeds": list(SEEDS),
            "k_dpi_values": list(K_VALUES),
            "p_value": P_VALUE,
            "p_token": P_TOKEN,
            "mi_cutoff": MI_CUTOFF,
            "npar": N_PAR,
            "dpi_epsilon": 0,
            "sampling": "fixed 80% without replacement",
            "source_and_annotation_lists": "same frozen panel (-s=panel, -l=panel)",
            "common_source_definition": "zero-based expression-row indices in H=1335 panel",
        },
        "execution": {
            "direct_runs": len(HUB_COUNTS) * len(SEEDS),
            "one_run_per_panel_seed": True,
            "candidate_adjacency_k_dpi": 1,
            "sidecar_option": "-W",
            "consensus_k_edge": "not run",
        },
        "baseline_policy": (
            "read-only; candidate K_DPI=1 adjacency data, DPI accounting, and sample schedule "
            "must exactly reproduce it"
        ),
    }
    work_root.mkdir(parents=True, exist_ok=True)
    ensure_exact_json(work_root / "screen_design.json", design)
    verify_frozen_harness_hashes(design["harness_sha256"])
    return baseline, candidate, design


def command_for_job(
    *, baseline: dict, candidate: dict, count: int, seed: int, output: Path, sidecar: Path
) -> list[str]:
    panel = Path(baseline["panels"][str(count)]["path"])
    return [
        str(candidate["binary"]),
        "-i",
        str(baseline["expression"]["path"]),
        "-l",
        str(panel),
        "-s",
        str(panel),
        "-p",
        P_TOKEN,
        "-e",
        "0",
        "-a",
        "adaptive_partitioning",
        "-H",
        str(candidate["config_directory"]).rstrip("/") + "/",
        "-N",
        str(N_PAR),
        "-S",
        str(seed),
        "-v",
        "on",
        "-o",
        str(output),
        "-u",
        "80%",
        "-M",
        str(candidate["null_model"]),
        "-W",
        str(sidecar),
    ]


def expected_sidecar_provenance(
    *, baseline: dict, count: int, executed_output: Path
) -> dict[str, str]:
    """Return path provenance emitted for the command's temporary output.

    The network and sidecar are atomically published after validation. The
    sidecar correctly retains the executed temporary network path, so resume
    revalidates that deterministic former path rather than rewriting provenance.
    """

    panel = str(baseline["panels"][str(count)]["path"])
    return {
        "input_file": str(baseline["expression"]["path"]),
        "input_adjacency_file": "",
        "network_output_file": str(executed_output),
        "subnetwork_file": panel,
        "annotation_file": panel,
    }


def baseline_seed_paths(baseline_root: Path, count: int, seed: int) -> dict[str, Path]:
    root = baseline_arm_root(baseline_root, count)
    stem = f"TF_run_{seed:03d}"
    return {
        "marker": root / "seed_metadata" / f"{stem}.json",
        "adjacency": root / "adjacency" / f"{stem}.adj",
        "stdout": root / "logs" / f"{stem}.stdout.log",
    }


def load_baseline_seed(baseline_root: Path, count: int, seed: int) -> tuple[dict, dict]:
    paths = baseline_seed_paths(baseline_root, count, seed)
    if any(not path.is_file() for path in paths.values()):
        missing = [name for name, path in paths.items() if not path.is_file()]
        raise FileNotFoundError(f"Missing baseline seed artifacts H={count}/seed={seed}: {missing}")
    marker = load_json(paths["marker"])
    if (
        marker.get("driver") != "sig"
        or int(marker.get("hub_count", -1)) != count
        or int(marker.get("seed", -1)) != seed
        or marker.get("binary_sha256") != EXPECTED_BASELINE["binary_sha256"]
        or marker.get("null_model_sha256") != EXPECTED_BASELINE["null_model_sha256"]
        or marker.get("source_tree_fingerprint") != EXPECTED_BASELINE["source_tree_fingerprint"]
        or marker.get("adjacency", {}).get("full_sha256") != sha256_file(paths["adjacency"])
        or marker.get("stdout_sha256") != sha256_file(paths["stdout"])
    ):
        raise RuntimeError(f"Frozen baseline seed is inconsistent: H={count}/seed={seed}")
    dpi = marker["dpi"]
    if int(dpi["pre_edges"]) != int(dpi["pruned_edges"]) + int(dpi["post_edges"]):
        raise RuntimeError(f"Baseline DPI accounting mismatch: H={count}/seed={seed}")
    return marker, paths


def validate_candidate(
    *,
    output: Path,
    sidecar_path: Path,
    stdout_path: Path,
    panel_ids: set[str],
    all_expression_ids: set[str],
    source_indices: set[int],
    baseline_marker: dict,
    expected_sidecar_paths: dict[str, str],
) -> tuple[dict, dict, dict, dict, dict]:
    adjacency = validate_adjacency(
        output, panel_ids=panel_ids, all_expression_ids=all_expression_ids
    )
    dpi = parse_dpi_stats(stdout_path)
    sampling = parse_sampling(stdout_path)
    sidecar = parse_witness_sidecar(
        sidecar_path,
        expected_source_indices=source_indices,
        expected_provenance=expected_sidecar_paths,
    )
    baseline_dpi = baseline_marker["dpi"]
    dpi_fields = ("pre_edges", "pruned_edges", "post_edges", "dpi_applied")
    dpi_match = all(dpi[field] == baseline_dpi[field] for field in dpi_fields)
    sampling_match = sampling["sha256"] == baseline_marker["sampling"]["sha256"]
    data_match = adjacency["data_sha256"] == baseline_marker["adjacency"]["data_sha256"]
    edge_count_match = adjacency["edges"] == baseline_marker["adjacency"]["edges"]
    if adjacency["edges"] != dpi["post_edges"]:
        raise ValueError(f"Adjacency/DPI post-edge mismatch for {output}")
    if not dpi_match:
        raise ValueError(f"Candidate K_DPI=1 DPI counts differ from frozen baseline for {output}")
    if not sampling_match:
        raise ValueError(f"Candidate sample differs from frozen baseline for {output}")
    if not data_match or not edge_count_match:
        raise ValueError(f"Candidate K_DPI=1 adjacency data differ from frozen baseline for {output}")
    totals = sidecar["totals"]
    if int(totals["pre_edges"]) != int(dpi["pre_edges"]):
        raise ValueError(f"Sidecar pre-edge total differs from DPI stats for {output}")
    if int(totals["witnesses_ge_1"]) != int(dpi["pruned_edges"]):
        raise ValueError(f"Sidecar K=1 total differs from DPI stats for {output}")
    if int(sidecar["provenance"]["k1_pruned_edges"]) != int(dpi["pruned_edges"]):
        raise ValueError(f"Sidecar provenance K=1 total differs from DPI stats for {output}")
    baseline_comparison = {
        "sampling_match": sampling_match,
        "dpi_stats_match": dpi_match,
        "adjacency_data_match": data_match,
        "adjacency_edge_count_match": edge_count_match,
        "baseline_adjacency_data_sha256": baseline_marker["adjacency"]["data_sha256"],
    }
    sidecar_summary = {
        "sha256": sidecar["sha256"],
        "bytes": sidecar["bytes"],
        "source_rows": sidecar["source_rows"],
        "provenance": sidecar["provenance"],
        "totals": totals,
    }
    return adjacency, dpi, sampling, sidecar_summary, baseline_comparison


def run_seed_job(
    *,
    work_root: Path,
    baseline_root: Path,
    baseline: dict,
    candidate: dict,
    design_sha256: str,
    count: int,
    seed: int,
    all_expression_ids: set[str],
    panel_ids: set[str],
    source_indices: set[int],
) -> tuple[str, bool, dict]:
    root = result_root(work_root, count)
    for relative in ("adjacency", "witness_sidecars", "logs", "seed_metadata", "work"):
        (root / relative).mkdir(parents=True, exist_ok=True)
    stem = f"TF_run_{seed:03d}"
    output = root / "adjacency" / f"{stem}.adj"
    sidecar_path = root / "witness_sidecars" / f"{stem}.tsv"
    marker_path = root / "seed_metadata" / f"{stem}.json"
    stdout_path = root / "logs" / f"{stem}.stdout.log"
    stderr_path = root / "logs" / f"{stem}.stderr.log"
    time_path = root / "logs" / f"{stem}.time.txt"
    partial_output = root / "work" / f"{stem}.adj.partial"
    partial_sidecar = root / "work" / f"{stem}.witness.tsv.partial"
    expected_sidecar_paths = expected_sidecar_provenance(
        baseline=baseline, count=count, executed_output=partial_output
    )
    baseline_marker, baseline_paths = load_baseline_seed(baseline_root, count, seed)
    preview = command_for_job(
        baseline=baseline,
        candidate=candidate,
        count=count,
        seed=seed,
        output=output,
        sidecar=sidecar_path,
    )
    fingerprint_payload = {
        "schema": SEED_SCHEMA,
        "design_sha256": design_sha256,
        "hub_count": count,
        "seed": seed,
        "candidate_binary_sha256": candidate["binary_sha256"],
        "expression_sha256": baseline["expression"]["sha256"],
        "panel_sha256": baseline["panels"][str(count)]["sha256"],
        "baseline_marker_sha256": sha256_file(baseline_paths["marker"]),
        "baseline_adjacency_sha256": baseline_marker["adjacency"]["full_sha256"],
        "command_without_outputs": [
            "<OUTPUT>" if value == str(output) else "<SIDECAR>" if value == str(sidecar_path) else value
            for value in preview
        ],
    }
    run_fingerprint = fingerprint(fingerprint_payload)

    completed = (marker_path.is_file(), output.is_file(), sidecar_path.is_file())
    if all(completed):
        existing = load_json(marker_path)
        adjacency, dpi, sampling, sidecar, comparison = validate_candidate(
            output=output,
            sidecar_path=sidecar_path,
            stdout_path=stdout_path,
            panel_ids=panel_ids,
            all_expression_ids=all_expression_ids,
            source_indices=source_indices,
            baseline_marker=baseline_marker,
            expected_sidecar_paths=expected_sidecar_paths,
        )
        if (
            existing.get("fingerprint") == run_fingerprint
            and existing.get("adjacency") == adjacency
            and existing.get("dpi") == dpi
            and existing.get("sampling") == sampling
            and existing.get("sidecar") == sidecar
            and existing.get("baseline_comparison") == comparison
            and existing.get("stdout_sha256") == sha256_file(stdout_path)
            and existing.get("stderr_sha256") == sha256_file(stderr_path)
        ):
            return f"{arm_key(count)}/{seed:03d}", True, existing
        raise RuntimeError(f"Stale or inconsistent completed seed: {marker_path}")
    if any(completed):
        raise RuntimeError(
            f"Orphan marker/adjacency/sidecar requires manual inspection: {arm_key(count)}/{seed:03d}"
        )

    for partial in (partial_output, partial_sidecar):
        if partial.exists():
            partial.unlink()
    command = command_for_job(
        baseline=baseline,
        candidate=candidate,
        count=count,
        seed=seed,
        output=partial_output,
        sidecar=partial_sidecar,
    )
    timed = [
        "/usr/bin/time",
        "-f",
        "elapsed_s=%e\nuser_s=%U\nsystem_s=%S\nmax_rss_kib=%M",
        "-o",
        str(time_path),
        *command,
    ]
    started = utc_now()
    clock = time.perf_counter()
    try:
        checked_run(timed, cwd=root, stdout_path=stdout_path, stderr_path=stderr_path)
    except Exception as error:
        raise RuntimeError(
            f"Seed failed ({arm_key(count)}/{seed}); see {stderr_path}"
        ) from error
    wall = time.perf_counter() - clock
    adjacency, dpi, sampling, sidecar, comparison = validate_candidate(
        output=partial_output,
        sidecar_path=partial_sidecar,
        stdout_path=stdout_path,
        panel_ids=panel_ids,
        all_expression_ids=all_expression_ids,
        source_indices=source_indices,
        baseline_marker=baseline_marker,
        expected_sidecar_paths=expected_sidecar_paths,
    )
    os.replace(partial_output, output)
    os.replace(partial_sidecar, sidecar_path)
    adjacency["full_sha256"] = sha256_file(output)
    sidecar["sha256"] = sha256_file(sidecar_path)
    record = {
        "schema": SEED_SCHEMA,
        "fingerprint": run_fingerprint,
        "design_sha256": design_sha256,
        "arm": arm_key(count),
        "driver": "sig",
        "hub_count": count,
        "seed": seed,
        "k_dpi_values": list(K_VALUES),
        "candidate_adjacency_k_dpi": 1,
        "executed_command": command,
        "canonical_command_after_atomic_publish": command_for_job(
            baseline=baseline,
            candidate=candidate,
            count=count,
            seed=seed,
            output=output,
            sidecar=sidecar_path,
        ),
        "candidate_source_tree_fingerprint": candidate["source_provenance"]["source_tree_fingerprint"],
        "candidate_binary_sha256": candidate["binary_sha256"],
        "config_sha256": candidate["config_sha256"],
        "null_model_sha256": candidate["null_model_sha256"],
        "expression_sha256": baseline["expression"]["sha256"],
        "panel_sha256": baseline["panels"][str(count)]["sha256"],
        "started_at_utc": started,
        "finished_at_utc": utc_now(),
        "high_resolution_wall_s": wall,
        "gnu_time": parse_gnu_time(time_path),
        "stdout_sha256": sha256_file(stdout_path),
        "stderr_sha256": sha256_file(stderr_path),
        "stderr_bytes": stderr_path.stat().st_size,
        "adjacency": adjacency,
        "dpi": dpi,
        "sampling": sampling,
        "sidecar": sidecar,
        "baseline": {
            "root": str(baseline_root.resolve()),
            "marker_path": str(baseline_paths["marker"]),
            "marker_sha256": sha256_file(baseline_paths["marker"]),
            "adjacency_path": str(baseline_paths["adjacency"]),
            "adjacency_sha256": baseline_marker["adjacency"]["full_sha256"],
        },
        "baseline_comparison": comparison,
    }
    atomic_json(marker_path, record)
    return f"{arm_key(count)}/{seed:03d}", False, record


def write_run_manifest(work_root: Path) -> None:
    rows: list[dict[str, object]] = []
    for marker_path in sorted(
        (work_root / "results").glob("sig_h*/seed_metadata/TF_run_*.json")
    ):
        record = load_json(marker_path)
        rows.append(
            {
                "arm": record["arm"],
                "hub_count": record["hub_count"],
                "seed": record["seed"],
                "candidate_binary_sha256": record["candidate_binary_sha256"],
                "elapsed_s": record["gnu_time"]["elapsed_s"],
                "max_rss_kib": record["gnu_time"]["max_rss_kib"],
                "pre_edges": record["dpi"]["pre_edges"],
                "k1_pruned_edges": record["dpi"]["pruned_edges"],
                "k1_post_edges": record["dpi"]["post_edges"],
                "adjacency_data_match": record["baseline_comparison"]["adjacency_data_match"],
                "dpi_stats_match": record["baseline_comparison"]["dpi_stats_match"],
                "sampling_match": record["baseline_comparison"]["sampling_match"],
                "sidecar_sha256": record["sidecar"]["sha256"],
            }
        )
    path = work_root / "run_manifest.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    fields = list(rows[0]) if rows else ["empty"]
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=("prepare", "infer", "all"), default="all")
    parser.add_argument("--repo", type=Path, default=REPO_DEFAULT)
    parser.add_argument("--work-root", type=Path, default=WORK_ROOT_DEFAULT)
    parser.add_argument("--baseline-root", type=Path, default=BASELINE_DEFAULT)
    parser.add_argument("--source-commit", default=EXPECTED_BASELINE["source_commit"])
    parser.add_argument("--hub-counts", default="all")
    parser.add_argument("--seed-start", type=int, default=SEEDS[0])
    parser.add_argument("--seed-end", type=int, default=SEEDS[-1])
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.workers < 1:
        raise ValueError("--workers must be positive")
    if args.seed_start < SEEDS[0] or args.seed_end > SEEDS[-1] or args.seed_start > args.seed_end:
        raise ValueError(f"Seed range must be within {SEEDS[0]}..{SEEDS[-1]}")
    counts = selected_counts(args.hub_counts)
    baseline, candidate, design = prepare(
        repo=args.repo.resolve(),
        work_root=args.work_root,
        baseline_root=args.baseline_root,
        source_commit=args.source_commit,
    )
    design_sha256 = sha256_file(args.work_root / "screen_design.json")
    print(
        f"[PREPARE] candidate={candidate['binary_sha256']} design={design_sha256} "
        f"work_root={args.work_root}"
    )
    if args.phase == "prepare":
        return 0

    expression_mapping, all_expression_ids = expression_index(
        Path(baseline["expression"]["path"])
    )
    panels: dict[int, set[str]] = {}
    source_indices: dict[int, set[int]] = {}
    for count in HUB_COUNTS:
        ids = set(read_unique_ids(Path(baseline["panels"][str(count)]["path"]), count))
        panels[count] = ids
        source_indices[count] = panel_source_indices(expression_mapping, ids)
    common = source_indices[HUB_COUNTS[0]]
    if any(not common.issubset(source_indices[count]) for count in HUB_COUNTS):
        raise RuntimeError("The frozen H=1335 common-source indices are not nested")

    jobs = [
        (count, seed)
        for count in counts
        for seed in range(args.seed_start, args.seed_end + 1)
    ]
    if args.dry_run:
        for count, seed in jobs:
            root = result_root(args.work_root, count)
            command = command_for_job(
                baseline=baseline,
                candidate=candidate,
                count=count,
                seed=seed,
                output=root / "adjacency" / f"TF_run_{seed:03d}.adj",
                sidecar=root / "witness_sidecars" / f"TF_run_{seed:03d}.tsv",
            )
            print("[DRY-RUN] " + " ".join(command))
        return 0

    completed = 0
    reused = 0
    failures: list[str] = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                run_seed_job,
                work_root=args.work_root,
                baseline_root=args.baseline_root,
                baseline=baseline,
                candidate=candidate,
                design_sha256=design_sha256,
                count=count,
                seed=seed,
                all_expression_ids=all_expression_ids,
                panel_ids=panels[count],
                source_indices=source_indices[count],
            ): (count, seed)
            for count, seed in jobs
        }
        for future in as_completed(futures):
            count, seed = futures[future]
            try:
                label, was_reused, _ = future.result()
                completed += 1
                reused += int(was_reused)
                print(f"[INFER] {completed}/{len(jobs)} {label} {'reused' if was_reused else 'completed'}")
            except Exception as error:
                message = f"{arm_key(count)}/{seed:03d}: {error}"
                failures.append(message)
                print(f"[FAIL] {message}", file=sys.stderr)
    write_run_manifest(args.work_root)
    print(f"[DONE] completed={completed} reused={reused} failures={len(failures)}")
    if failures:
        return 2
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
