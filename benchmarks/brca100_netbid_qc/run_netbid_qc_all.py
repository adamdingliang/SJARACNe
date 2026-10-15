#!/usr/bin/env python3
"""Run the six matched NetBID2 QC reports after consensus generation."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess


STAGES = ("baseline_12113fb", "pr66_5809183", "pr67_7633ebb")
DRIVERS = {
    "tf": ("BRCA100_TF.txt", "TF_"),
    "sig": ("BRCA100_SIG.txt", "SIG_"),
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    os.replace(temporary, path)


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

    wrapper = args.benchmark_repo / "benchmarks/brca100_netbid_qc/netbid2-r"
    qc_script = args.benchmark_repo / "benchmarks/brca100_netbid_qc/run_netbid_qc.R"
    records = []

    for stage in STAGES:
        for driver, (driver_filename, prefix) in DRIVERS.items():
            arm_root = args.work_root / "results" / stage / driver
            consensus = arm_root / "consensus" / "consensus_network_ncol_.txt"
            driver_file = args.work_root / "inputs" / driver_filename
            qc_root = arm_root / "netbid2_qc"
            partial_qc_root = arm_root / "netbid2_qc.partial"
            manifest = arm_root / "netbid2_qc_manifest.json"
            pending_manifest = arm_root / "netbid2_qc_manifest.pending.json"
            expected_html = qc_root / f"{prefix}netQC.html"
            expected_summary = qc_root / "network_summary.tsv"
            expected_target_sizes = qc_root / "driver_target_sizes.tsv"
            stdout_path = arm_root / "logs" / "netbid2_qc.stdout.log"
            stderr_path = arm_root / "logs" / "netbid2_qc.stderr.log"
            command = [
                str(wrapper),
                "Rscript",
                str(qc_script),
                str(consensus),
                str(driver_file),
                str(partial_qc_root),
                prefix,
            ]
            fingerprint = hashlib.sha256(
                (
                    sha256_file(consensus)
                    + sha256_file(driver_file)
                    + sha256_file(qc_script)
                    + "NetBID2-2.2.0-5defa454"
                ).encode("ascii")
            ).hexdigest()

            if (
                manifest.is_file()
                and expected_html.is_file()
                and expected_summary.is_file()
                and expected_target_sizes.is_file()
            ):
                existing = json.loads(manifest.read_text(encoding="utf-8"))
                if (
                    existing.get("fingerprint") == fingerprint
                    and existing.get("html_sha256") == sha256_file(expected_html)
                    and existing.get("network_summary_sha256")
                    == sha256_file(expected_summary)
                    and existing.get("driver_target_sizes_sha256")
                    == sha256_file(expected_target_sizes)
                ):
                    print(f"[QC] {stage}/{driver} resume", flush=True)
                    if pending_manifest.exists():
                        pending_manifest.unlink()
                    records.append(existing)
                    continue
                raise RuntimeError(f"Stale QC output at {qc_root}")
            if qc_root.exists() and not manifest.exists():
                if not pending_manifest.is_file():
                    raise RuntimeError(f"Unverifiable orphan QC directory: {qc_root}")
                for required in (expected_html, expected_summary, expected_target_sizes):
                    if not required.is_file() or required.stat().st_size == 0:
                        raise RuntimeError(f"Incomplete recovered QC output: {required}")
                record = json.loads(pending_manifest.read_text(encoding="utf-8"))
                expected_hashes = {
                    "fingerprint": fingerprint,
                    "html_sha256": sha256_file(expected_html),
                    "network_summary_sha256": sha256_file(expected_summary),
                    "driver_target_sizes_sha256": sha256_file(
                        expected_target_sizes
                    ),
                }
                if any(record.get(key) != value for key, value in expected_hashes.items()):
                    raise RuntimeError(f"Orphan QC fails pending manifest: {qc_root}")
                record["recovered_after_interrupted_manifest"] = True
                atomic_json(manifest, record)
                pending_manifest.unlink()
                records.append(record)
                print(f"[QC] {stage}/{driver} recovered", flush=True)
                continue
            if qc_root.exists():
                raise RuntimeError(f"Unexpected existing QC directory {qc_root}")
            if partial_qc_root.exists():
                shutil.rmtree(partial_qc_root)
            if pending_manifest.exists():
                pending_manifest.unlink()

            print(f"[QC] {stage}/{driver}", flush=True)
            with stdout_path.open("w", encoding="utf-8", newline="\n") as stdout, \
                    stderr_path.open("w", encoding="utf-8", newline="\n") as stderr:
                subprocess.run(command, stdout=stdout, stderr=stderr, check=True)
            partial_outputs = (
                partial_qc_root / f"{prefix}netQC.html",
                partial_qc_root / "network_summary.tsv",
                partial_qc_root / "driver_target_sizes.tsv",
            )
            for required in partial_outputs:
                if not required.is_file() or required.stat().st_size == 0:
                    raise RuntimeError(f"Missing NetBID2 output {required}")
            record = {
                "fingerprint": fingerprint,
                "stage": stage,
                "driver": driver,
                "command": command,
                "consensus_sha256": sha256_file(consensus),
                "html": str(expected_html),
                "html_sha256": sha256_file(partial_outputs[0]),
                "html_bytes": partial_outputs[0].stat().st_size,
                "network_summary_sha256": sha256_file(
                    partial_outputs[1]
                ),
                "driver_target_sizes_sha256": sha256_file(
                    partial_outputs[2]
                ),
                "stdout_sha256": sha256_file(stdout_path),
                "stderr_sha256": sha256_file(stderr_path),
                "stderr_bytes": stderr_path.stat().st_size,
            }
            atomic_json(pending_manifest, record)
            os.replace(partial_qc_root, qc_root)
            atomic_json(manifest, record)
            pending_manifest.unlink()
            records.append(record)

    atomic_json(args.work_root / "results" / "netbid2_qc_manifest.json", {"runs": records})
    print("[QC] all six reports complete", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
