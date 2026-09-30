#!/usr/bin/env python3
"""Reconstruct retained-edge support across the 100 adjacency networks."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess


STAGES = ("baseline_12113fb", "pr66_5809183", "pr67_7633ebb")
DRIVERS = ("tf", "sig")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    os.replace(partial, path)


def main() -> int:
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
            / "brca100-netbid-qc-20260817-rerun"
        ),
    )
    args = parser.parse_args()

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

    binary_hash = sha256_file(binary)
    records = []
    for stage in STAGES:
        for driver in DRIVERS:
            arm_root = args.work_root / "results" / stage / driver
            consensus = arm_root / "consensus/consensus_network_ncol_.txt"
            adjacency = arm_root / "adjacency"
            consensus_manifest = arm_root / "consensus_manifest.json"
            output = arm_root / "consensus/consensus_support.tsv"
            manifest = arm_root / "support_summary_manifest.json"
            adjacency_files = sorted(adjacency.glob("TF_run_*.adj"))
            if len(adjacency_files) != 100:
                raise RuntimeError(
                    f"Expected 100 adjacency inputs for {stage}/{driver}, "
                    f"got {len(adjacency_files)}"
                )
            adjacency_hashes = [sha256_file(path) for path in adjacency_files]
            consensus_sha256 = sha256_file(consensus)
            consensus_record = json.loads(
                consensus_manifest.read_text(encoding="utf-8")
            )
            expected_consensus_fingerprint = fingerprint(
                {
                    "stage": stage,
                    "driver": driver,
                    "adjacency_hashes": adjacency_hashes,
                    "consensus_p": 1e-5,
                    "consensus_script_sha256": sha256_file(
                        args.benchmark_repo
                        / "SJARACNe/bin/create_consensus_network.py"
                    ),
                }
            )
            if (
                consensus_record.get("fingerprint")
                != expected_consensus_fingerprint
                or consensus_record.get("ncol", {}).get("sha256")
                != consensus_sha256
            ):
                raise RuntimeError(
                    f"Consensus provenance does not match adjacency inputs: "
                    f"{stage}/{driver}"
                )
            run_fingerprint = fingerprint(
                {
                    "stage": stage,
                    "driver": driver,
                    "consensus_sha256": consensus_sha256,
                    "consensus_fingerprint": expected_consensus_fingerprint,
                    "consensus_manifest_sha256": sha256_file(consensus_manifest),
                    "source_sha256": sha256_file(source),
                    "binary_sha256": binary_hash,
                    "adjacency_sha256": adjacency_hashes,
                }
            )
            if manifest.is_file() and output.is_file():
                existing = json.loads(manifest.read_text(encoding="utf-8"))
                if (
                    existing.get("fingerprint") == run_fingerprint
                    and existing.get("output_sha256") == sha256_file(output)
                ):
                    print(f"[SUPPORT] {stage}/{driver} resume", flush=True)
                    records.append(existing)
                    continue
                raise RuntimeError(f"Stale support output: {output}")
            stdout_path = arm_root / "logs/support_summary.stdout.log"
            stderr_path = arm_root / "logs/support_summary.stderr.log"
            time_path = arm_root / "logs/support_summary.time.txt"
            run_command = [
                "/usr/bin/time",
                "-f", "elapsed_s=%e\nuser_s=%U\nsystem_s=%S\nmax_rss_kib=%M",
                "-o", str(time_path),
                str(binary), str(consensus), str(adjacency), str(output),
            ]
            print(f"[SUPPORT] {stage}/{driver}", flush=True)
            with stdout_path.open("w", encoding="utf-8", newline="\n") as stdout, \
                    stderr_path.open("w", encoding="utf-8", newline="\n") as stderr:
                subprocess.run(run_command, stdout=stdout, stderr=stderr, check=True)
            if not output.is_file() or output.stat().st_size == 0:
                raise RuntimeError(f"Missing support output: {output}")
            with output.open("r", encoding="utf-8", newline="") as handle:
                rows = sum(1 for _ in handle) - 1
            record = {
                "fingerprint": run_fingerprint,
                "stage": stage,
                "driver": driver,
                "command": run_command,
                "binary_sha256": binary_hash,
                "output": str(output),
                "output_sha256": sha256_file(output),
                "retained_edges": rows,
            }
            atomic_json(manifest, record)
            records.append(record)

    atomic_json(args.work_root / "results/support_summary_manifest.json", records)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
