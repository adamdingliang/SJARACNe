#!/usr/bin/env python3
"""Run a one-pass recurrence-threshold sweep on fixed SJARACNe subnetworks.

The expensive AP-MI calculation and DPI pruning are inputs to this benchmark.
For each driver class, the 100 immutable post-DPI adjacency files are aggregated
once.  The resulting recurrence table is then materialized at every integer
minimum support count from 6 through 20.

The Poisson-binomial probabilities reported here are exact for SJARACNe's
current plug-in edge-occupancy model after freezing its estimated inputs:

    q_i = E_i / U,

where E_i is the edge count in bootstrap run i and U is the observed union of
ordered edges.  This script deliberately calls those probabilities
"plug-in model tails" and never labels them FDR estimates.  They are not
tails conditioned on an edge having appeared in the observed union.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from typing import Any, Iterable, Iterator


SCHEMA = "sjaracne-brca100-consensus-recurrence-sweep-v1"
ARM_SCHEMA = "sjaracne-brca100-consensus-recurrence-arm-v1"
MATERIALIZED_SCHEMA = "sjaracne-brca100-consensus-recurrence-network-v1"
ANALYSIS_SCHEMA = "sjaracne-brca100-consensus-recurrence-analysis-v1"
BOOTSTRAP_RUNS = 100
MIN_SUPPORTS = tuple(range(6, 21))
CURRENT_LEGACY_SUPPORT = 9
EXPECTED_COMMIT = "7633ebb4a0d966dbda15a4e32d0efa492fb71aeb"
EXPECTED_BINARY_SHA256 = (
    "61180915be61455887ac261499455ac933cf976b292cab1a8180d759acc3ac2d"
)
EXPECTED_INPUT_SHA256 = {
    "BRCA100.exp": "ad8a334f5f8cdf46a1000d3ee259b35258a18b3da2e314bb3a0cf7a421d98bc8",
    "BRCA100_TF.txt": "9b1219a489b99432175e4c4ad46add7b06f25aae388ee8dd3261fa91e4c43ffd",
    "BRCA100_SIG.txt": "16ca27df655f16684f880a4ad719c4e2ae3f8dc0d7e6b9eccdd24cd97c40797c",
}
EXPECTED_NETBID_ENVIRONMENT = {
    "R": "R version 4.4.3 (2025-02-28)",
    "NetBID2": "2.2.0",
    "NetBID2_remote_sha": "5defa454d600b94f5dd6d1f9f4428f99759a6821",
    "igraph": "2.3.3",
}
ARMS = {
    "tf": {
        "point": "p1e-03",
        "p_value": 1e-3,
        "driver_file": "BRCA100_TF.txt",
        "prefix": "TF_",
        "expected_k9_edges": 224608,
    },
    "sig": {
        "point": "p5e-04",
        "p_value": 5e-4,
        "driver_file": "BRCA100_SIG.txt",
        "prefix": "SIG_",
        "expected_k9_edges": 379053,
    },
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def serialized_json(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def atomic_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("wb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def atomic_json(path: Path, payload: object) -> None:
    atomic_bytes(path, serialized_json(payload))


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def fingerprint(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def ensure_exact_json(path: Path, payload: object, description: str) -> None:
    expected = serialized_json(payload)
    if path.exists():
        if path.read_bytes() != expected:
            raise ValueError(f"Existing {description} is incompatible: {path}")
        return
    atomic_bytes(path, expected)


def poisson_binomial_tails(probabilities: Iterable[float]) -> list[float]:
    """Return tails[k] = P(S >= k) for independent Bernoulli probabilities."""

    values = [float(value) for value in probabilities]
    if any(not math.isfinite(value) or value < 0.0 or value > 1.0 for value in values):
        raise ValueError("Poisson-binomial probabilities must be finite and in [0, 1]")
    probability_mass = [1.0] + [0.0] * len(values)
    processed = 0
    for probability in values:
        for successes in range(processed + 1, 0, -1):
            probability_mass[successes] = (
                probability_mass[successes] * (1.0 - probability)
                + probability_mass[successes - 1] * probability
            )
        probability_mass[0] *= 1.0 - probability
        processed += 1
    tails = [0.0] * (len(values) + 1)
    running = 0.0
    for successes in range(len(values), -1, -1):
        running += probability_mass[successes]
        tails[successes] = min(1.0, max(0.0, running))
    tails[0] = 1.0
    return tails


def legacy_uprob(z_score: float) -> float:
    """Reproduce the upper-normal-tail approximation used by SJARACNe."""

    probability = 0.0
    absolute = abs(z_score)
    if absolute < 1.9:
        probability = (
            1
            + absolute
            * (
                0.049867347
                + absolute
                * (
                    0.0211410061
                    + absolute * 0.0032776263
                    + absolute
                    * (
                        0.0000380036
                        + absolute * (0.0000488906 + absolute * 0.000005383)
                    )
                )
            )
        ) ** (-16) / 2
    elif absolute <= 100:
        for index in range(18, 0, -1):
            probability = index / (absolute + probability)
        probability = (
            math.exp(-0.5 * absolute * absolute)
            / math.sqrt(2 * math.pi)
            / (absolute + probability)
        )
    if z_score < 0:
        probability = 1 - probability
    return probability


def parse_metric_tsv(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["metric", "value"]:
            raise ValueError(f"Unexpected metric table header: {path}")
        result: dict[str, str] = {}
        for row in reader:
            metric = row["metric"]
            if metric in result:
                raise ValueError(f"Duplicate metric {metric!r} in {path}")
            result[metric] = row["value"]
    return result


def parse_bootstrap_info(path: Path) -> dict[str, float]:
    result: dict[str, float] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if ":" not in raw_line:
            continue
        key, value = raw_line.split(":", 1)
        value = value.strip()
        if value.startswith("N/A"):
            continue
        try:
            result[key.strip()] = float(value)
        except ValueError:
            continue
    return result


def load_source_manifest(source_root: Path) -> list[dict[str, str]]:
    path = source_root / "results" / "run_manifest.tsv"
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 2600:
        raise ValueError(f"Expected 2,600 source inference rows, got {len(rows)}")
    return rows


def adjacency_inventory(
    source_root: Path,
    driver: str,
    source_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    arm = ARMS[driver]
    adjacency_root = source_root / "results" / str(arm["point"]) / driver / "adjacency"
    files = sorted(adjacency_root.glob("*.adj"))
    expected_names = [f"TF_run_{seed:03d}.adj" for seed in range(1, BOOTSTRAP_RUNS + 1)]
    if [path.name for path in files] != expected_names:
        raise ValueError(f"Adjacency inventory is not the exact seed 1..100 set: {adjacency_root}")
    source_index = {
        int(row["seed"]): row
        for row in source_rows
        if row["point"] == arm["point"] and row["driver"] == driver
    }
    if set(source_index) != set(range(1, BOOTSTRAP_RUNS + 1)):
        raise ValueError(f"Source manifest is incomplete for {driver}")
    records: list[dict[str, object]] = []
    for seed, path in enumerate(files, start=1):
        row = source_index[seed]
        actual_hash = sha256_file(path)
        actual_bytes = path.stat().st_size
        if actual_hash != row["adjacency_sha256"]:
            raise ValueError(f"Adjacency SHA-256 mismatch: {path}")
        if actual_bytes != int(row["adjacency_bytes"]):
            raise ValueError(f"Adjacency size mismatch: {path}")
        records.append(
            {
                "seed": seed,
                "path": str(path),
                "bytes": actual_bytes,
                "sha256": actual_hash,
                "data_sha256": row["data_sha256"],
                "edges": int(row["edges"]),
            }
        )
    return records


def validate_point_manifest(source_root: Path, driver: str) -> dict[str, Any]:
    arm = ARMS[driver]
    path = source_root / "results" / str(arm["point"]) / "point_manifest.json"
    point = read_json(path)
    required = {
        "key": arm["point"],
        "commit": EXPECTED_COMMIT,
        "binary_sha256": EXPECTED_BINARY_SHA256,
        "m": 80,
        "npar": 40,
        "dpi_epsilon": 0,
        "sampling": "fixed 80% without replacement",
    }
    for field, expected in required.items():
        if point.get(field) != expected:
            raise ValueError(f"Unexpected {field} in {path}: {point.get(field)!r}")
    if not math.isclose(float(point["p_value"]), float(arm["p_value"]), rel_tol=0, abs_tol=1e-15):
        raise ValueError(f"Unexpected per-subsample p-value in {path}")
    if point.get("tail_extrapolated") is not False:
        raise ValueError(f"The selected point must be inside the held-out range: {path}")
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "payload": point,
    }


def prepare_design(repo_root: Path, source_root: Path, work_root: Path) -> dict[str, Any]:
    if work_root.resolve() == source_root.resolve():
        raise ValueError("The recurrence work root must differ from the source sweep root")
    source_design_path = source_root / "sweep_design.json"
    source_run_manifest_path = source_root / "results" / "run_manifest.tsv"
    source_rows = load_source_manifest(source_root)
    inputs: dict[str, object] = {}
    for filename, expected_hash in EXPECTED_INPUT_SHA256.items():
        path = source_root / "inputs" / filename
        actual_hash = sha256_file(path)
        if actual_hash != expected_hash:
            raise ValueError(f"Pinned input SHA-256 mismatch: {path}")
        inputs[filename] = {
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": actual_hash,
        }
    arm_records: dict[str, object] = {}
    for driver in ("tf", "sig"):
        arm_records[driver] = {
            "source_point": validate_point_manifest(source_root, driver),
            "adjacencies": adjacency_inventory(source_root, driver, source_rows),
            "driver_file": inputs[str(ARMS[driver]["driver_file"])],
            "per_subsample_p": ARMS[driver]["p_value"],
        }
    source_design = read_json(source_design_path)
    if source_design.get("commit") != EXPECTED_COMMIT:
        raise ValueError("Unexpected source sweep commit")
    design_core = {
        "schema": SCHEMA,
        "purpose": (
            "Sweep explicit recurrence counts while freezing AP-MI, m=80, Npar=40, "
            "DPI=0, seed networks, and BRCA100 inputs"
        ),
        "source_work_root": str(source_root),
        "source_sweep_design": {
            "path": str(source_design_path),
            "sha256": sha256_file(source_design_path),
        },
        "source_run_manifest": {
            "path": str(source_run_manifest_path),
            "sha256": sha256_file(source_run_manifest_path),
        },
        "repo_commit": EXPECTED_COMMIT,
        "benchmark_scripts": {
            "runner_sha256": sha256_file(Path(__file__).resolve()),
            "aggregator_source_sha256": sha256_file(
                repo_root
                / "benchmarks"
                / "brca100_consensus_recurrence_sweep"
                / "aggregate_recurrence.cpp"
            ),
            "netbid_r_sha256": sha256_file(
                repo_root
                / "benchmarks"
                / "brca100_pr67_threshold_sweep"
                / "run_netbid_qc.R"
            ),
        },
        "inputs": inputs,
        "arms": arm_records,
        "bootstrap_runs": BOOTSTRAP_RUNS,
        "minimum_supports": list(MIN_SUPPORTS),
        "legacy_consensus_p": 1e-5,
        "legacy_minimum_support": CURRENT_LEGACY_SUPPORT,
        "occupancy_null": {
            "definition": "q_i = E_i / U; independent Bernoulli occurrences across runs",
            "frozen_plugin_inputs": "observed ordered-edge union U and run sizes E_i",
            "tail": "untruncated Poisson-binomial P(S >= K); not conditioned on S >= 1",
            "warning": "Plug-in stability model; not truth, FDR, or external validation",
        },
    }
    design = dict(design_core)
    design["fingerprint"] = fingerprint(design_core)
    work_root.mkdir(parents=True, exist_ok=True)
    ensure_exact_json(work_root / "design.json", design, "recurrence sweep design")
    return design


def build_aggregator(repo_root: Path, work_root: Path) -> dict[str, Any]:
    source = repo_root / "benchmarks" / "brca100_consensus_recurrence_sweep" / "aggregate_recurrence.cpp"
    if not source.is_file():
        raise ValueError(f"Missing recurrence aggregator source: {source}")
    binary_root = work_root / "bin"
    binary_root.mkdir(parents=True, exist_ok=True)
    binary = binary_root / "aggregate_recurrence"
    flags = ["-O3", "-std=c++11", "-Wall", "-Wextra", "-pedantic"]
    compiler_version = subprocess.check_output(["g++", "--version"], text=True).splitlines()[0]
    build_core = {
        "source": str(source),
        "source_sha256": sha256_file(source),
        "compiler": compiler_version,
        "flags": flags,
    }
    expected_fingerprint = fingerprint(build_core)
    manifest_path = binary_root / "build_manifest.json"
    if manifest_path.exists() and binary.exists():
        manifest = read_json(manifest_path)
        if (
            manifest.get("fingerprint") == expected_fingerprint
            and manifest.get("binary_sha256") == sha256_file(binary)
        ):
            return manifest
    partial = binary.with_name(binary.name + ".partial")
    command = ["g++", *flags, str(source), "-o", str(partial)]
    subprocess.run(command, check=True)
    os.replace(partial, binary)
    binary.chmod(0o755)
    manifest = dict(build_core)
    manifest.update(
        {
            "fingerprint": expected_fingerprint,
            "command": command,
            "binary": str(binary),
            "binary_sha256": sha256_file(binary),
        }
    )
    atomic_json(manifest_path, manifest)
    return manifest


def read_run_counts(path: Path) -> list[dict[str, object]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != BOOTSTRAP_RUNS:
        raise ValueError(f"Expected 100 run counts in {path}")
    result: list[dict[str, object]] = []
    for expected, row in enumerate(rows, start=1):
        if int(row["run_ordinal"]) != expected:
            raise ValueError(f"Non-canonical run order in {path}")
        result.append(
            {
                "run_ordinal": expected,
                "adjacency_file": row["adjacency_file"],
                "edge_count": int(row["edge_count"]),
            }
        )
    return result


def parse_recurrence_rows(path: Path) -> Iterator[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = [
            "source",
            "target",
            "mean_observed_MI",
            "consensus_MI",
            "support_count",
            "support_fraction",
        ]
        if reader.fieldnames != expected:
            raise ValueError(f"Unexpected recurrence header: {path}")
        for row in reader:
            yield row


def validate_recurrence_table(path: Path) -> dict[int, int]:
    counts = {support: 0 for support in MIN_SUPPORTS}
    previous: tuple[str, str] | None = None
    rows = 0
    for row in parse_recurrence_rows(path):
        key = (row["source"], row["target"])
        if previous is not None and key <= previous:
            raise ValueError(f"Recurrence rows are not unique and sorted: {path}")
        previous = key
        support = int(row["support_count"])
        fraction = float(row["support_fraction"])
        mean_mi = float(row["mean_observed_MI"])
        consensus_mi = float(row["consensus_MI"])
        if support < MIN_SUPPORTS[0] or support > BOOTSTRAP_RUNS:
            raise ValueError(f"Invalid support count in {path}: {support}")
        if not math.isclose(fraction, support / BOOTSTRAP_RUNS, rel_tol=0, abs_tol=1e-15):
            raise ValueError(f"Invalid support fraction in {path}")
        if not math.isfinite(mean_mi) or not math.isfinite(consensus_mi) or mean_mi <= 0:
            raise ValueError(f"Invalid MI in {path}")
        if f"{mean_mi:.4f}" != row["consensus_MI"]:
            raise ValueError(f"Consensus MI is not the required four-decimal value in {path}")
        for threshold in MIN_SUPPORTS:
            if support >= threshold:
                counts[threshold] += 1
        rows += 1
    if rows == 0:
        raise ValueError(f"Recurrence table is empty: {path}")
    if any(counts[left] < counts[right] for left, right in zip(MIN_SUPPORTS, MIN_SUPPORTS[1:])):
        raise ValueError("Recurrence counts must be nonincreasing")
    return counts


def validate_k9_anchor(source_root: Path, driver: str, recurrence_path: Path) -> None:
    arm = ARMS[driver]
    existing = (
        source_root
        / "results"
        / str(arm["point"])
        / driver
        / "consensus"
        / "consensus_support.tsv"
    )
    with existing.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        source_rows = list(reader)
    recurrence_rows = [
        row for row in parse_recurrence_rows(recurrence_path)
        if int(row["support_count"]) >= CURRENT_LEGACY_SUPPORT
    ]
    if len(recurrence_rows) != int(arm["expected_k9_edges"]):
        raise ValueError(f"K=9 recurrence count mismatch for {driver}")
    if len(source_rows) != len(recurrence_rows):
        raise ValueError(f"Published consensus size mismatch for {driver}")
    for index, (current, prior) in enumerate(zip(recurrence_rows, source_rows)):
        if (current["source"], current["target"]) != (prior["source"], prior["target"]):
            raise ValueError(f"K=9 edge mismatch for {driver} at row {index + 2}")
        if int(current["support_count"]) != int(prior["support_count"]):
            raise ValueError(f"K=9 support mismatch for {driver} at row {index + 2}")
        if current["consensus_MI"] != f"{float(prior['consensus_MI']):.4f}":
            raise ValueError(f"K=9 MI mismatch for {driver} at row {index + 2}")


def aggregate_arm(
    repo_root: Path,
    source_root: Path,
    work_root: Path,
    design: dict[str, Any],
    build: dict[str, Any],
    driver: str,
) -> dict[str, Any]:
    arm_root = work_root / "aggregate" / driver
    manifest_path = arm_root / "aggregate_manifest.json"
    edge_path = arm_root / "recurrence_edges.tsv"
    run_counts_path = arm_root / "run_counts.tsv"
    summary_path = arm_root / "aggregate_summary.tsv"
    expected_input = {
        "design_fingerprint": design["fingerprint"],
        "build_fingerprint": build["fingerprint"],
        "driver": driver,
        "source_point": ARMS[driver]["point"],
    }
    if manifest_path.exists():
        manifest = read_json(manifest_path)
        if manifest.get("input") != expected_input:
            raise ValueError(f"Incompatible aggregate manifest: {manifest_path}")
        for field, path in (
            ("edge_sha256", edge_path),
            ("run_counts_sha256", run_counts_path),
            ("summary_sha256", summary_path),
        ):
            if manifest.get(field) != sha256_file(path):
                raise ValueError(f"Aggregate output hash mismatch: {path}")
        for field, path in (
            ("plugin_tail_sha256", arm_root / "plugin_tail.tsv"),
            ("stdout_sha256", arm_root / "aggregate.stdout.log"),
            ("stderr_sha256", arm_root / "aggregate.stderr.log"),
        ):
            if manifest.get(field) != sha256_file(path):
                raise ValueError(f"Aggregate provenance hash mismatch: {path}")
        return manifest
    if arm_root.exists():
        shutil.rmtree(arm_root)
    partial_root = arm_root.with_name(arm_root.name + ".partial")
    if partial_root.exists():
        shutil.rmtree(partial_root)
    partial_root.mkdir(parents=True)
    command = [
        str(build["binary"]),
        str(source_root / "results" / str(ARMS[driver]["point"]) / driver / "adjacency"),
        str(partial_root / edge_path.name),
        str(partial_root / run_counts_path.name),
        str(partial_root / summary_path.name),
    ]
    started = time.monotonic()
    completed = subprocess.run(command, text=True, capture_output=True)
    elapsed = time.monotonic() - started
    (partial_root / "aggregate.stdout.log").write_text(completed.stdout, encoding="utf-8")
    (partial_root / "aggregate.stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"Recurrence aggregation failed for {driver}:\n{completed.stderr}"
        )
    os.replace(partial_root, arm_root)
    counts = validate_recurrence_table(edge_path)
    validate_k9_anchor(source_root, driver, edge_path)
    run_counts = read_run_counts(run_counts_path)
    summary = parse_metric_tsv(summary_path)
    union_edges = int(summary["union_edges"])
    if int(summary["bootstrap_runs"]) != BOOTSTRAP_RUNS:
        raise ValueError(f"Unexpected bootstrap count for {driver}")
    if int(summary["minimum_support"]) != MIN_SUPPORTS[0]:
        raise ValueError(f"Unexpected minimum support for {driver}")
    if int(summary["retained_edges"]) != counts[MIN_SUPPORTS[0]]:
        raise ValueError(f"Retained-edge summary mismatch for {driver}")
    probabilities = [int(row["edge_count"]) / union_edges for row in run_counts]
    tails = poisson_binomial_tails(probabilities)
    mu = sum(probabilities)
    variance = sum(value * (1.0 - value) for value in probabilities)
    sigma = math.sqrt(variance)
    prior_bootstrap = parse_bootstrap_info(
        source_root
        / "results"
        / str(ARMS[driver]["point"])
        / driver
        / "consensus"
        / "bootstrap_info_.txt"
    )
    if int(prior_bootstrap["Total edge tested"]) != union_edges:
        raise ValueError(f"Union-edge mismatch versus published consensus for {driver}")
    if not math.isclose(prior_bootstrap["mu"], mu, rel_tol=0, abs_tol=5e-13):
        raise ValueError(f"mu mismatch versus published consensus for {driver}")
    if not math.isclose(prior_bootstrap["sigma"], sigma, rel_tol=0, abs_tol=5e-13):
        raise ValueError(f"sigma mismatch versus published consensus for {driver}")
    tail_path = arm_root / "plugin_tail.tsv"
    with tail_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "minimum_support",
                "support_fraction",
                "plugin_poisson_binomial_tail",
                "legacy_normal_tail",
                "plugin_null_edge_burden_proxy",
                "observed_edges",
            ]
        )
        for threshold in MIN_SUPPORTS:
            normal_tail = legacy_uprob((threshold - mu) / sigma) if sigma else 0.0
            writer.writerow(
                [
                    threshold,
                    threshold / BOOTSTRAP_RUNS,
                    f"{tails[threshold]:.17g}",
                    f"{normal_tail:.17g}",
                    f"{union_edges * tails[threshold]:.17g}",
                    counts[threshold],
                ]
            )
    manifest = {
        "schema": ARM_SCHEMA,
        "input": expected_input,
        "command": command,
        "elapsed_seconds": elapsed,
        "union_edges": union_edges,
        "run_edge_counts": [int(row["edge_count"]) for row in run_counts],
        "mu": mu,
        "sigma": sigma,
        "edge_counts_by_minimum_support": {str(key): value for key, value in counts.items()},
        "k9_anchor_reproduced": True,
        "edge_path": str(edge_path),
        "edge_bytes": edge_path.stat().st_size,
        "edge_sha256": sha256_file(edge_path),
        "run_counts_sha256": sha256_file(run_counts_path),
        "summary_sha256": sha256_file(summary_path),
        "plugin_tail_sha256": sha256_file(tail_path),
        "stdout_sha256": sha256_file(arm_root / "aggregate.stdout.log"),
        "stderr_sha256": sha256_file(arm_root / "aggregate.stderr.log"),
    }
    atomic_json(manifest_path, manifest)
    return manifest


def write_base_three_column(recurrence_path: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["source", "target", "MI"])
        for row in parse_recurrence_rows(recurrence_path):
            writer.writerow([row["source"], row["target"], row["consensus_MI"]])


def enhance_base_network(repo_root: Path, source_root: Path, work_root: Path, driver: str) -> Path:
    aggregate_root = work_root / "aggregate" / driver
    recurrence_path = aggregate_root / "recurrence_edges.tsv"
    enhanced_root = aggregate_root / "enhanced_minimum_k006"
    ncol_path = enhanced_root / "consensus_network_ncol_.txt"
    manifest_path = enhanced_root / "enhanced_manifest.json"
    input_core = {
        "recurrence_sha256": sha256_file(recurrence_path),
        "expression_sha256": EXPECTED_INPUT_SHA256["BRCA100.exp"],
        "enhancer_sha256": sha256_file(
            repo_root / "SJARACNe" / "bin" / "create_consensus_network.py"
        ),
        "minimum_support": MIN_SUPPORTS[0],
    }
    if manifest_path.exists():
        manifest = read_json(manifest_path)
        if manifest.get("input") != input_core:
            raise ValueError(f"Incompatible enhanced-network manifest: {manifest_path}")
        if manifest.get("ncol_sha256") != sha256_file(ncol_path):
            raise ValueError(f"Enhanced network hash mismatch: {ncol_path}")
        return ncol_path
    if enhanced_root.exists():
        shutil.rmtree(enhanced_root)
    partial_root = enhanced_root.with_name(enhanced_root.name + ".partial")
    if partial_root.exists():
        shutil.rmtree(partial_root)
    partial_root.mkdir(parents=True)
    three_col = partial_root / "consensus_network_3col_.txt"
    write_base_three_column(recurrence_path, three_col)
    sys.path.insert(0, str(repo_root))
    try:
        from SJARACNe.bin.create_consensus_network import create_enhanced_consensus_network
        started = time.monotonic()
        create_enhanced_consensus_network(
            str(source_root / "inputs" / "BRCA100.exp"),
            str(three_col),
            str(partial_root),
        )
        elapsed = time.monotonic() - started
    finally:
        sys.path.pop(0)
    os.replace(partial_root, enhanced_root)
    manifest = {
        "schema": "sjaracne-brca100-consensus-recurrence-enhanced-base-v1",
        "input": input_core,
        "elapsed_seconds": elapsed,
        "three_col_sha256": sha256_file(enhanced_root / three_col.name),
        "ncol_sha256": sha256_file(ncol_path),
        "rows": sum(1 for _ in ncol_path.open("r", encoding="utf-8")) - 1,
    }
    atomic_json(manifest_path, manifest)
    return ncol_path


def materialize_arm(
    source_root: Path,
    work_root: Path,
    driver: str,
    aggregate_manifest: dict[str, Any],
    base_ncol: Path,
) -> list[dict[str, Any]]:
    results_root = work_root / "results" / driver
    completed_manifests = [
        results_root / f"k{threshold:03d}" / "network_manifest.json"
        for threshold in MIN_SUPPORTS
    ]
    if all(path.exists() for path in completed_manifests):
        manifests = [read_json(path) for path in completed_manifests]
        for manifest, threshold in zip(manifests, MIN_SUPPORTS):
            if manifest.get("minimum_support") != threshold:
                raise ValueError(f"Invalid materialized threshold for {driver}")
            for field, filename in (
                ("three_col_sha256", "consensus_network_3col_.txt"),
                ("ncol_sha256", "consensus_network_ncol_.txt"),
                ("support_sha256", "consensus_support.tsv"),
                ("bootstrap_info_sha256", "bootstrap_info_.txt"),
                ("parameter_info_sha256", "parameter_info_.txt"),
            ):
                path = results_root / f"k{threshold:03d}" / filename
                if manifest.get(field) != sha256_file(path):
                    raise ValueError(f"Materialized output hash mismatch: {path}")
            expected_aggregate = sha256_file(
                work_root / "aggregate" / driver / "aggregate_manifest.json"
            )
            if manifest.get("aggregate_manifest_sha256") != expected_aggregate:
                raise ValueError(
                    f"Materialized aggregate linkage mismatch for {driver} K={threshold}"
                )
        return manifests
    if results_root.exists():
        shutil.rmtree(results_root)
    partial_root = results_root.with_name(results_root.name + ".partial")
    if partial_root.exists():
        shutil.rmtree(partial_root)
    handles: dict[int, dict[str, Any]] = {}
    try:
        for threshold in MIN_SUPPORTS:
            threshold_root = partial_root / f"k{threshold:03d}"
            threshold_root.mkdir(parents=True)
            three = (threshold_root / "consensus_network_3col_.txt").open(
                "w", encoding="utf-8", newline=""
            )
            ncol = (threshold_root / "consensus_network_ncol_.txt").open(
                "w", encoding="utf-8", newline=""
            )
            support = (threshold_root / "consensus_support.tsv").open(
                "w", encoding="utf-8", newline=""
            )
            three.write("source\ttarget\tMI\n")
            ncol.write(
                "source\ttarget\tsource.symbol\ttarget.symbol\tMI\tpearson\t"
                "spearman\tslope\tp-value\n"
            )
            support.write(
                "source\ttarget\tconsensus_MI\tsupport_count\tsupport_fraction\t"
                "mean_observed_MI\n"
            )
            handles[threshold] = {
                "root": threshold_root,
                "three": three,
                "ncol": ncol,
                "support": support,
                "rows": 0,
            }
        recurrence_iterator = parse_recurrence_rows(
            work_root / "aggregate" / driver / "recurrence_edges.tsv"
        )
        with base_ncol.open("r", encoding="utf-8", newline="") as ncol_input:
            ncol_header = ncol_input.readline()
            expected_header = (
                "source\ttarget\tsource.symbol\ttarget.symbol\tMI\tpearson\t"
                "spearman\tslope\tp-value\n"
            )
            if ncol_header != expected_header:
                raise ValueError(f"Unexpected enhanced-network header: {base_ncol}")
            recurrence_count = 0
            for recurrence, ncol_line in zip(recurrence_iterator, ncol_input):
                fields = ncol_line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise ValueError(f"Malformed enhanced-network row: {base_ncol}")
                if (recurrence["source"], recurrence["target"]) != (fields[0], fields[1]):
                    raise ValueError("Enhanced and recurrence edge order differs")
                if recurrence["consensus_MI"] != fields[4]:
                    raise ValueError("Enhanced and recurrence MI differs")
                support_count = int(recurrence["support_count"])
                for threshold in MIN_SUPPORTS:
                    if support_count < threshold:
                        continue
                    output = handles[threshold]
                    output["three"].write(
                        f"{fields[0]}\t{fields[1]}\t{fields[4]}\n"
                    )
                    output["ncol"].write(ncol_line)
                    output["support"].write(
                        "\t".join(
                            [
                                fields[0],
                                fields[1],
                                fields[4],
                                recurrence["support_count"],
                                recurrence["support_fraction"],
                                recurrence["mean_observed_MI"],
                            ]
                        )
                        + "\n"
                    )
                    output["rows"] += 1
                recurrence_count += 1
            if ncol_input.readline() != "":
                raise ValueError("Enhanced network has more rows than recurrence table")
            try:
                next(recurrence_iterator)
                raise ValueError("Recurrence table has more rows than enhanced network")
            except StopIteration:
                pass
            if recurrence_count != aggregate_manifest["edge_counts_by_minimum_support"]["6"]:
                raise ValueError("Enhanced network row count mismatch")
    finally:
        for output in handles.values():
            for name in ("three", "ncol", "support"):
                output[name].close()
    tail_rows: dict[int, dict[str, str]] = {}
    tail_path = work_root / "aggregate" / driver / "plugin_tail.tsv"
    with tail_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            tail_rows[int(row["minimum_support"])] = row
    manifests: list[dict[str, Any]] = []
    for threshold in MIN_SUPPORTS:
        output = handles[threshold]
        threshold_root = Path(output["root"])
        row_count = int(output["rows"])
        expected_rows = int(
            aggregate_manifest["edge_counts_by_minimum_support"][str(threshold)]
        )
        if row_count != expected_rows:
            raise ValueError(f"Materialized row mismatch for {driver} K={threshold}")
        tail = tail_rows[threshold]
        bootstrap_text = (
            f"Total edge tested: {aggregate_manifest['union_edges']}\n"
            f"Bootstrap No: {BOOTSTRAP_RUNS}\n"
            f"mu: {aggregate_manifest['mu']:.17g}\n"
            f"sigma: {aggregate_manifest['sigma']:.17g}\n"
            f"Minimum support count: {threshold}\n"
            f"Minimum support fraction: {threshold / BOOTSTRAP_RUNS:.17g}\n"
            f"Plug-in Poisson-binomial tail: {tail['plugin_poisson_binomial_tail']}\n"
            f"Legacy normal-approximation tail: {tail['legacy_normal_tail']}\n"
            "Warning: untruncated plug-in occupancy model; not FDR or biological truth\n"
        )
        (threshold_root / "bootstrap_info_.txt").write_text(
            bootstrap_text, encoding="utf-8"
        )
        parameter_text = (
            f">  Source per-subsample p: {ARMS[driver]['p_value']}\n"
            f">  Source point: {ARMS[driver]['point']}\n"
            ">  Sampling: fixed 80% without replacement\n"
            ">  Npar: 40\n"
            ">  DPI epsilon: 0\n"
            f">  Consensus recurrence rule: support_count >= {threshold}\n"
            ">  Consensus tail method: plug-in Poisson-binomial using frozen U and E_i\n"
        )
        (threshold_root / "parameter_info_.txt").write_text(
            parameter_text, encoding="utf-8"
        )
        manifest = {
            "schema": MATERIALIZED_SCHEMA,
            "driver": driver,
            "source_point": ARMS[driver]["point"],
            "per_subsample_p": ARMS[driver]["p_value"],
            "minimum_support": threshold,
            "support_fraction": threshold / BOOTSTRAP_RUNS,
            "edges": row_count,
            "plugin_poisson_binomial_tail": float(
                tail["plugin_poisson_binomial_tail"]
            ),
            "legacy_normal_tail": float(tail["legacy_normal_tail"]),
            "plugin_null_edge_burden_proxy": float(
                tail["plugin_null_edge_burden_proxy"]
            ),
            "aggregate_manifest_sha256": sha256_file(
                work_root / "aggregate" / driver / "aggregate_manifest.json"
            ),
            "three_col_sha256": sha256_file(
                threshold_root / "consensus_network_3col_.txt"
            ),
            "ncol_sha256": sha256_file(
                threshold_root / "consensus_network_ncol_.txt"
            ),
            "support_sha256": sha256_file(threshold_root / "consensus_support.tsv"),
            "bootstrap_info_sha256": sha256_file(
                threshold_root / "bootstrap_info_.txt"
            ),
            "parameter_info_sha256": sha256_file(
                threshold_root / "parameter_info_.txt"
            ),
        }
        atomic_json(threshold_root / "network_manifest.json", manifest)
        manifests.append(manifest)
    os.replace(partial_root, results_root)
    return manifests


def probe_netbid_environment(wrapper: Path, r_script: Path) -> dict[str, str]:
    completed = subprocess.run(
        [str(wrapper), "Rscript", str(r_script), "--probe"],
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"NetBID2 environment probe failed:\n{completed.stderr}")
    rows = list(csv.DictReader(completed.stdout.splitlines(), delimiter="\t"))
    environment: dict[str, str] = {}
    for row in rows:
        component = row.get("component", "")
        if not component or component in environment:
            raise ValueError("Malformed or duplicate NetBID2 environment component")
        environment[component] = row.get("version", "")
    if environment != EXPECTED_NETBID_ENVIRONMENT:
        raise ValueError(
            f"Pinned NetBID2 environment mismatch: {environment!r}"
        )
    return environment


def run_netbid(repo_root: Path, source_root: Path, work_root: Path) -> None:
    wrapper = repo_root / "benchmarks" / "brca100_netbid_qc" / "netbid2-r"
    r_script = repo_root / "benchmarks" / "brca100_pr67_threshold_sweep" / "run_netbid_qc.R"
    environment = probe_netbid_environment(wrapper, r_script)
    for driver in ("tf", "sig"):
        driver_file = source_root / "inputs" / str(ARMS[driver]["driver_file"])
        for threshold in MIN_SUPPORTS:
            threshold_root = work_root / "results" / driver / f"k{threshold:03d}"
            network_manifest = read_json(threshold_root / "network_manifest.json")
            final_root = threshold_root / "netbid2_qc"
            manifest_path = threshold_root / "netbid2_manifest.json"
            input_core = {
                "network_manifest_sha256": sha256_file(
                    threshold_root / "network_manifest.json"
                ),
                "network_sha256": network_manifest["ncol_sha256"],
                "driver_sha256": EXPECTED_INPUT_SHA256[str(ARMS[driver]["driver_file"])],
                "r_script_sha256": sha256_file(r_script),
                "wrapper_sha256": sha256_file(wrapper),
                "environment": environment,
                "generate_html": False,
            }
            if manifest_path.exists():
                manifest = read_json(manifest_path)
                if manifest.get("input") != input_core:
                    raise ValueError(f"Incompatible NetBID2 manifest: {manifest_path}")
                for filename, expected_hash in manifest["outputs"].items():
                    if sha256_file(final_root / filename) != expected_hash:
                        raise ValueError(f"NetBID2 output hash mismatch: {final_root / filename}")
                continue
            if final_root.exists():
                shutil.rmtree(final_root)
            partial_root = final_root.with_name(final_root.name + ".partial")
            if partial_root.exists():
                shutil.rmtree(partial_root)
            command = [
                str(wrapper),
                "Rscript",
                str(r_script),
                str(threshold_root / "consensus_network_ncol_.txt"),
                str(driver_file),
                str(partial_root),
                str(ARMS[driver]["prefix"]),
                "false",
            ]
            completed = subprocess.run(command, text=True, capture_output=True)
            (threshold_root / "netbid2.stdout.log").write_text(
                completed.stdout, encoding="utf-8"
            )
            (threshold_root / "netbid2.stderr.log").write_text(
                completed.stderr, encoding="utf-8"
            )
            if completed.returncode != 0:
                raise RuntimeError(
                    f"NetBID2 summary failed for {driver} K={threshold}:\n{completed.stderr}"
                )
            os.replace(partial_root, final_root)
            outputs = {
                path.name: sha256_file(path)
                for path in sorted(final_root.iterdir())
                if path.is_file()
            }
            required = {
                "network_summary.tsv",
                "driver_target_sizes.tsv",
                "netbid_environment.tsv",
            }
            if set(outputs) != required:
                raise ValueError(f"Unexpected NetBID2 output inventory: {final_root}")
            manifest = {
                "schema": "sjaracne-brca100-consensus-recurrence-netbid2-v1",
                "driver": driver,
                "minimum_support": threshold,
                "input": input_core,
                "command": command,
                "outputs": outputs,
                "stdout_sha256": sha256_file(threshold_root / "netbid2.stdout.log"),
                "stderr_sha256": sha256_file(threshold_root / "netbid2.stderr.log"),
            }
            atomic_json(manifest_path, manifest)


def read_netbid_metrics(path: Path) -> dict[str, float]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["metric", "value"]:
            raise ValueError(f"Unexpected NetBID2 summary header: {path}")
        result: dict[str, float] = {}
        for row in reader:
            result[row["metric"]] = float(row["value"])
    return result


def quantile(sorted_values: list[float], probability: float) -> float:
    if not sorted_values:
        return math.nan
    position = (len(sorted_values) - 1) * probability
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return sorted_values[lower]
    return sorted_values[lower] + (position - lower) * (
        sorted_values[upper] - sorted_values[lower]
    )


def analyze(work_root: Path) -> None:
    analysis_root = work_root / "analysis"
    if analysis_root.exists():
        manifest_path = analysis_root / "analysis_manifest.json"
        if not manifest_path.is_file():
            raise ValueError(f"Incomplete existing analysis root: {analysis_root}")
        manifest = read_json(manifest_path)
        expected_outputs = {
            "network_summary_sha256": analysis_root / "network_summary.tsv",
            "driver_target_coverage_sha256": analysis_root / "driver_target_coverage.tsv",
            "plot_png_sha256": analysis_root / "plots" / "recurrence_density_coverage.png",
            "plot_svg_sha256": analysis_root / "plots" / "recurrence_density_coverage.svg",
        }
        for field, path in expected_outputs.items():
            if manifest.get(field) != sha256_file(path):
                raise ValueError(f"Existing analysis output hash mismatch: {path}")
        print(f"Validated existing analysis: {analysis_root}")
        return
    partial_root = analysis_root.with_name(analysis_root.name + ".partial")
    if partial_root.exists():
        shutil.rmtree(partial_root)
    plots_root = partial_root / "plots"
    plots_root.mkdir(parents=True)
    rows: list[dict[str, object]] = []
    target_rows: list[dict[str, object]] = []
    for driver in ("tf", "sig"):
        for threshold in MIN_SUPPORTS:
            threshold_root = work_root / "results" / driver / f"k{threshold:03d}"
            manifest = read_json(threshold_root / "network_manifest.json")
            netbid = read_netbid_metrics(
                threshold_root / "netbid2_qc" / "network_summary.tsv"
            )
            with (
                threshold_root / "netbid2_qc" / "driver_target_sizes.tsv"
            ).open("r", encoding="utf-8", newline="") as handle:
                target_sizes = [
                    int(row["target_count"])
                    for row in csv.DictReader(handle, delimiter="\t")
                ]
            if int(netbid["edges"]) != int(manifest["edges"]):
                raise ValueError(
                    f"NetBID2/network edge count mismatch for {driver} K={threshold}"
                )
            if len(target_sizes) != int(netbid["candidate_drivers"]):
                raise ValueError(
                    f"NetBID2 candidate-driver count mismatch for {driver} K={threshold}"
                )
            if sum(target_sizes) != int(manifest["edges"]):
                raise ValueError(
                    f"Target-size sum does not equal edge count for {driver} K={threshold}"
                )
            support_values: list[int] = []
            mi_values: list[float] = []
            with (threshold_root / "consensus_support.tsv").open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    support_values.append(int(row["support_count"]))
                    mi_values.append(float(row["mean_observed_MI"]))
            support_values.sort()
            mi_values.sort()
            if len(support_values) != manifest["edges"]:
                raise ValueError(f"Support/network edge count mismatch for {driver} K={threshold}")
            target_record = {
                "driver": driver,
                "minimum_support": threshold,
                "candidate_drivers": len(target_sizes),
            }
            for minimum_targets in (1, 10, 20, 50, 100):
                count = sum(value >= minimum_targets for value in target_sizes)
                target_record[f"drivers_ge_{minimum_targets}_targets"] = count
                target_record[f"fraction_ge_{minimum_targets}_targets"] = count / len(target_sizes)
            target_rows.append(target_record)
            rows.append(
                {
                    "driver": driver,
                    "per_subsample_p": ARMS[driver]["p_value"],
                    "minimum_support": threshold,
                    "minimum_support_fraction": threshold / BOOTSTRAP_RUNS,
                    "plugin_poisson_binomial_tail": manifest[
                        "plugin_poisson_binomial_tail"
                    ],
                    "legacy_normal_tail": manifest["legacy_normal_tail"],
                    "plugin_null_edge_burden_proxy": manifest[
                        "plugin_null_edge_burden_proxy"
                    ],
                    "consensus_edges": int(netbid["edges"]),
                    "active_drivers": int(netbid["active_drivers"]),
                    "active_driver_fraction": netbid["active_driver_fraction"],
                    "incident_nodes": int(netbid["incident_nodes"]),
                    "weak_components": int(netbid["weak_components"]),
                    "largest_component_fraction": netbid[
                        "largest_weak_component_fraction"
                    ],
                    "median_targets_zero_filled": netbid["target_size_zero_median"],
                    "q25_targets_zero_filled": netbid["target_size_zero_q25"],
                    "q75_targets_zero_filled": netbid["target_size_zero_q75"],
                    "scale_free_adjusted_r2": netbid["scale_free_adjusted_r2"],
                    "support_mean": sum(support_values) / len(support_values),
                    "support_median": quantile(support_values, 0.5),
                    "support_q25": quantile(support_values, 0.25),
                    "support_q75": quantile(support_values, 0.75),
                    "edges_support_ge_20": sum(value >= 20 for value in support_values),
                    "edges_support_ge_50": sum(value >= 50 for value in support_values),
                    "edges_support_ge_80": sum(value >= 80 for value in support_values),
                    "mi_mean": sum(mi_values) / len(mi_values),
                    "mi_median": quantile(mi_values, 0.5),
                }
            )
    for driver in ("tf", "sig"):
        driver_rows = [row for row in rows if row["driver"] == driver]
        edge_counts = [int(row["consensus_edges"]) for row in driver_rows]
        if any(left < right for left, right in zip(edge_counts, edge_counts[1:])):
            raise ValueError(f"Consensus edge counts are not nested for {driver}")
        current = next(
            row for row in driver_rows
            if row["minimum_support"] == CURRENT_LEGACY_SUPPORT
        )
        if current["consensus_edges"] != ARMS[driver]["expected_k9_edges"]:
            raise ValueError(f"K=9 analysis anchor failed for {driver}")
    summary_path = partial_root / "network_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    target_path = partial_root / "driver_target_coverage.tsv"
    with target_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(target_rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(target_rows)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matplotlib.rcParams["svg.hashsalt"] = "sjaracne-recurrence-sweep-v1"
    figure, axes = plt.subplots(3, 2, figsize=(11, 10), sharex=True)
    metrics = [
        ("consensus_edges", "Consensus edges", "log"),
        ("median_targets_zero_filled", "Median targets", "linear"),
        ("active_driver_fraction", "Active-driver fraction", "linear"),
    ]
    for column, driver in enumerate(("tf", "sig")):
        driver_rows = [row for row in rows if row["driver"] == driver]
        x = [int(row["minimum_support"]) for row in driver_rows]
        for row_index, (field, ylabel, scale) in enumerate(metrics):
            axis = axes[row_index, column]
            axis.plot(x, [float(row[field]) for row in driver_rows], marker="o")
            axis.axvline(CURRENT_LEGACY_SUPPORT, color="black", linestyle="--", linewidth=1)
            axis.set_title(driver.upper())
            axis.set_ylabel(ylabel)
            axis.set_yscale(scale)
            axis.grid(alpha=0.25)
            if row_index == len(metrics) - 1:
                axis.set_xlabel("Minimum recurrence count K (out of 100)")
    figure.suptitle(
        "BRCA100 recurrence sweep; dashed line is legacy consensus boundary K=9"
    )
    figure.tight_layout()
    for suffix in ("png", "svg"):
        figure.savefig(
            plots_root / f"recurrence_density_coverage.{suffix}",
            dpi=180,
            metadata={"Date": None},
        )
    plt.close(figure)
    manifest_core = {
        "schema": ANALYSIS_SCHEMA,
        "design_sha256": sha256_file(work_root / "design.json"),
        "minimum_supports": list(MIN_SUPPORTS),
        "drivers": ["tf", "sig"],
        "network_summary_sha256": sha256_file(summary_path),
        "driver_target_coverage_sha256": sha256_file(target_path),
        "plot_png_sha256": sha256_file(
            plots_root / "recurrence_density_coverage.png"
        ),
        "plot_svg_sha256": sha256_file(
            plots_root / "recurrence_density_coverage.svg"
        ),
        "selection": None,
        "selection_note": (
            "No operating K is selected without downstream activity or biological validation"
        ),
    }
    manifest = dict(manifest_core)
    manifest["fingerprint"] = fingerprint(manifest_core)
    atomic_json(partial_root / "analysis_manifest.json", manifest)
    os.replace(partial_root, analysis_root)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    parser.add_argument("--source-work-root", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument(
        "--phase",
        choices=("prepare", "aggregate", "materialize", "netbid", "analyze", "all"),
        default="all",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    source_root = args.source_work_root.resolve()
    work_root = args.work_root.resolve()
    design = prepare_design(repo_root, source_root, work_root)
    if args.phase == "prepare":
        print(f"Prepared immutable design: {work_root / 'design.json'}")
        return 0
    build = build_aggregator(repo_root, work_root)
    aggregates: dict[str, dict[str, Any]] = {}
    for driver in ("tf", "sig"):
        aggregates[driver] = aggregate_arm(
            repo_root, source_root, work_root, design, build, driver
        )
    if args.phase == "aggregate":
        print("Aggregated both driver classes")
        return 0
    for driver in ("tf", "sig"):
        base_ncol = enhance_base_network(repo_root, source_root, work_root, driver)
        materialize_arm(
            source_root, work_root, driver, aggregates[driver], base_ncol
        )
    if args.phase == "materialize":
        print("Materialized K=6..20 for both driver classes")
        return 0
    run_netbid(repo_root, source_root, work_root)
    if args.phase == "netbid":
        print("Completed NetBID2 summaries")
        return 0
    analyze(work_root)
    print(f"Completed recurrence sweep: {work_root / 'analysis'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
