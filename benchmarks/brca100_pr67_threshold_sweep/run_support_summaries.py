#!/usr/bin/env python3
"""Reconstruct retained-edge support for every completed PR67 sweep arm."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess


DRIVERS = ("tf", "sig")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    os.replace(partial, path)


def load_json(path: Path) -> object:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def selected_values(specification: str, available: list[str], label: str) -> list[str]:
    if specification == "all":
        return available
    requested = [item.strip() for item in specification.split(",") if item.strip()]
    unknown = sorted(set(requested) - set(available))
    if unknown:
        raise ValueError(f"Unknown {label}(s): {', '.join(unknown)}")
    if not requested or len(set(requested)) != len(requested):
        raise ValueError(f"Empty or duplicate {label} selection")
    return requested


def validate_support(path: Path, expected_edges: int) -> int:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing/empty support output: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        expected_header = [
            "source",
            "target",
            "consensus_MI",
            "support_count",
            "support_fraction",
            "mean_observed_MI",
            "consensus_MI_roundtrip_match",
        ]
        if header != expected_header:
            raise ValueError(f"Unexpected support header in {path}: {header}")
        rows = 0
        edges: set[tuple[str, str]] = set()
        for line_number, line in enumerate(handle, 2):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(expected_header) or not fields[0] or not fields[1]:
                raise ValueError(f"Malformed support row {path}:{line_number}")
            consensus_mi = float(fields[2])
            count = int(fields[3])
            fraction = float(fields[4])
            mean_mi = float(fields[5])
            roundtrip_match = fields[6]
            edge = (fields[0], fields[1])
            if edge in edges or edge[0] == edge[1]:
                raise ValueError(f"Duplicate/self edge at {path}:{line_number}")
            edges.add(edge)
            if (
                not math.isfinite(fraction)
                or not 1 <= count <= 100
                or abs(fraction - count / 100.0) > 1e-12
            ):
                raise ValueError(f"Invalid support at {path}:{line_number}")
            if not (
                math.isfinite(consensus_mi)
                and math.isfinite(mean_mi)
                and consensus_mi > 0.0
                and mean_mi > 0.0
            ):
                raise ValueError(f"Invalid MI at {path}:{line_number}")
            if roundtrip_match.lower() not in {"true", "1", "yes"}:
                raise ValueError(
                    f"Consensus/mean MI round-trip mismatch at {path}:{line_number}"
                )
            rows += 1
    if rows != expected_edges:
        raise ValueError(
            f"Support rows disagree with consensus edges for {path}: "
            f"{rows} != {expected_edges}"
        )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--benchmark-repo",
        type=Path,
        default=Path("/mnt/d/GitHub/SJARACNe-brca100-netbid-qc"),
    )
    parser.add_argument(
        "--work-root",
        type=Path,
        default=(
            Path.home()
            / "sjaracne-benchmarks"
            / "brca100-pr67-threshold-sweep-20260818"
        ),
    )
    parser.add_argument("--points", default="all")
    parser.add_argument("--drivers", default="all")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results_root = args.work_root / "results"
    sweep_design_path = args.work_root / "sweep_design.json"
    sweep_design = load_json(sweep_design_path)
    if not isinstance(sweep_design, dict):
        raise ValueError(f"Expected a JSON object in {sweep_design_path}")
    expected_points = [str(item["key"]) for item in sweep_design["all_points"]]
    if not expected_points or len(set(expected_points)) != len(expected_points):
        raise ValueError(f"Invalid point list in {sweep_design_path}")
    point_roots = sorted(
        path.parent for path in results_root.glob("*/point_manifest.json")
    )
    point_map = {path.name: path for path in point_roots}
    if set(point_map) != set(expected_points):
        raise ValueError(
            "Point manifests do not exactly match sweep_design.json: "
            f"expected={sorted(expected_points)}, found={sorted(point_map)}"
        )
    points = selected_values(args.points, expected_points, "point")
    drivers = selected_values(args.drivers, list(DRIVERS), "driver")

    source = (
        args.benchmark_repo
        / "benchmarks/brca100_netbid_qc/summarize_consensus_support.cpp"
    )
    tool_root = args.work_root / "tools"
    tool_root.mkdir(exist_ok=True)
    binary = tool_root / "summarize_consensus_support"
    build_log = tool_root / "summarize_consensus_support.build.log"
    command = [
        "g++", "-O3", "-std=c++11", "-Wall", "-Wextra", "-Wpedantic",
        "-Wconversion", "-Wshadow", "-o", str(binary), str(source),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    build_log.write_text(
        completed.stdout + completed.stderr, encoding="utf-8", newline="\n"
    )
    if completed.returncode != 0:
        raise RuntimeError(f"Support helper build failed; see {build_log}")
    compiler = subprocess.run(
        ["g++", "--version"], text=True, capture_output=True, check=True
    ).stdout.splitlines()[0]
    binary_hash = sha256_file(binary)
    source_hash = sha256_file(source)
    consensus_script_hash = sha256_file(
        args.benchmark_repo / "SJARACNe/bin/create_consensus_network.py"
    )

    processed_records = []
    for point in points:
        point_manifest_path = point_map[point] / "point_manifest.json"
        point_manifest = load_json(point_manifest_path)
        if not isinstance(point_manifest, dict) or point_manifest.get("key") != point:
            raise ValueError(f"Invalid point manifest identity: {point_manifest_path}")
        for driver in drivers:
            arm_root = point_map[point] / driver
            adjacency = arm_root / "adjacency"
            consensus = arm_root / "consensus/consensus_network_ncol_.txt"
            consensus_manifest = arm_root / "consensus_manifest.json"
            output = arm_root / "consensus/consensus_support.tsv"
            manifest = arm_root / "support_summary_manifest.json"
            pending = arm_root / "support_summary_manifest.pending.json"
            adjacency_files = sorted(adjacency.glob("TF_run_*.adj"))
            if len(adjacency_files) != 100:
                raise RuntimeError(
                    f"Expected 100 adjacency inputs for {point}/{driver}, "
                    f"got {len(adjacency_files)}"
                )
            adjacency_hashes = [sha256_file(path) for path in adjacency_files]
            consensus_record = load_json(consensus_manifest)
            consensus_hash = sha256_file(consensus)
            expected_consensus_fingerprint = fingerprint(
                {
                    "stage": point,
                    "driver": driver,
                    "adjacency_hashes": adjacency_hashes,
                    "consensus_p": 1e-5,
                    "consensus_script_sha256": consensus_script_hash,
                }
            )
            if (
                consensus_record.get("fingerprint") != expected_consensus_fingerprint
                or consensus_record.get("ncol", {}).get("sha256") != consensus_hash
            ):
                raise RuntimeError(
                    f"Consensus provenance does not match inputs: {point}/{driver}"
                )
            expected_edges = int(consensus_record["ncol"]["edges"])
            consensus_manifest_hash = sha256_file(consensus_manifest)
            point_manifest_hash = sha256_file(point_manifest_path)
            fingerprint_payload = {
                "schema": "sjaracne-brca100-pr67-p-sweep-support-v1",
                "point": point,
                "p_value": point_manifest["p_value"],
                "driver": driver,
                "consensus_sha256": consensus_hash,
                "consensus_fingerprint": expected_consensus_fingerprint,
                "consensus_manifest_sha256": consensus_manifest_hash,
                "point_manifest_sha256": point_manifest_hash,
                "source_sha256": source_hash,
                "binary_sha256": binary_hash,
                "adjacency_sha256": adjacency_hashes,
            }
            run_fingerprint = fingerprint(fingerprint_payload)
            if manifest.is_file() and output.is_file():
                existing = load_json(manifest)
                rows = validate_support(output, expected_edges)
                if (
                    existing.get("fingerprint") == run_fingerprint
                    and existing.get("output_sha256") == sha256_file(output)
                    and existing.get("retained_edges") == rows
                ):
                    if pending.is_file():
                        pending_record = load_json(pending)
                        if pending_record != existing:
                            raise RuntimeError(
                                f"Final and pending support manifests disagree: {arm_root}"
                            )
                        pending.unlink()
                    processed_records.append(existing)
                    print(f"[SUPPORT] {point}/{driver} resume", flush=True)
                    continue
                raise RuntimeError(f"Stale support output: {output}")

            if pending.is_file() and output.is_file():
                recovery = load_json(pending)
                rows = validate_support(output, expected_edges)
                if (
                    recovery.get("fingerprint") != run_fingerprint
                    or recovery.get("output_sha256") != sha256_file(output)
                    or recovery.get("retained_edges") != rows
                ):
                    raise RuntimeError(f"Invalid pending support recovery: {pending}")
                atomic_json(manifest, recovery)
                pending.unlink()
                processed_records.append(recovery)
                print(f"[SUPPORT] {point}/{driver} recovered", flush=True)
                continue

            if output.exists() or manifest.exists():
                raise RuntimeError(f"Unverifiable partial support state: {arm_root}")
            temporary = output.with_name(output.name + ".partial")
            if temporary.exists():
                temporary.unlink()
            stdout_path = arm_root / "logs/support_summary.stdout.log"
            stderr_path = arm_root / "logs/support_summary.stderr.log"
            time_path = arm_root / "logs/support_summary.time.txt"
            run_command = [
                "/usr/bin/time",
                "-f", "elapsed_s=%e\nuser_s=%U\nsystem_s=%S\nmax_rss_kib=%M",
                "-o", str(time_path),
                str(binary), str(consensus), str(adjacency), str(temporary),
            ]
            print(f"[SUPPORT] {point}/{driver}", flush=True)
            with stdout_path.open("w", encoding="utf-8", newline="\n") as stdout, \
                    stderr_path.open("w", encoding="utf-8", newline="\n") as stderr:
                subprocess.run(run_command, stdout=stdout, stderr=stderr, check=True)
            rows = validate_support(temporary, expected_edges)
            record = {
                "schema": "sjaracne-brca100-pr67-p-sweep-support-v1",
                "fingerprint": run_fingerprint,
                "point": point,
                "p_value": point_manifest["p_value"],
                "driver": driver,
                "command": run_command,
                "build_command": command,
                "compiler": compiler,
                "source_sha256": source_hash,
                "binary_sha256": binary_hash,
                "point_manifest_sha256": point_manifest_hash,
                "consensus_sha256": consensus_hash,
                "consensus_fingerprint": expected_consensus_fingerprint,
                "consensus_manifest_sha256": consensus_manifest_hash,
                "adjacency_sha256": adjacency_hashes,
                "output": str(output),
                "output_sha256": sha256_file(temporary),
                "retained_edges": rows,
            }
            atomic_json(pending, record)
            os.replace(temporary, output)
            atomic_json(manifest, record)
            pending.unlink()
            processed_records.append(record)

    atomic_json(
        results_root / "support_summary_manifest.json",
        {
            "schema": "sjaracne-brca100-pr67-p-sweep-support-aggregate-v1",
            "sweep_design_sha256": sha256_file(sweep_design_path),
            "points": points,
            "drivers": drivers,
            "records": processed_records,
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
