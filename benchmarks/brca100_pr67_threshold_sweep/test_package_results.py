#!/usr/bin/env python3
"""Synthetic contract tests for package_results.py."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import package_results as pkg


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def inventory(root: Path) -> list[dict[str, object]]:
    return [
        {
            "path": path.relative_to(root).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": pkg.sha256_file(path),
        }
        for path in sorted(
            (item for item in root.rglob("*") if item.is_file()),
            key=lambda item: item.relative_to(root).as_posix(),
        )
    ]


def finalized(value: dict[str, object]) -> dict[str, object]:
    result = dict(value)
    result["fingerprint"] = pkg.json_fingerprint(result)
    return result


def make_fixture(root: Path) -> tuple[dict[str, str], str]:
    inputs: dict[str, dict[str, object]] = {}
    expected_input_hashes: dict[str, str] = {}
    for filename, content in (
        ("BRCA100.exp", b"expression\n"),
        ("BRCA100_TF.txt", b"tf\n"),
        ("BRCA100_SIG.txt", b"sig\n"),
    ):
        path = root / "inputs" / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
        checksum = pkg.sha256_file(path)
        expected_input_hashes[filename] = checksum
        inputs[filename] = {"path": str(path), "sha256": checksum, "bytes": len(content)}
    inputs["expression_id_count"] = {"count": 28278}

    build_root = root / "builds" / "pr67_7633ebb"
    binary = build_root / "bin" / "sjaracne.exe"
    config = build_root / "source" / "SJARACNe" / "config"
    model = config / "apmi_null" / "apmi_null_m00080_npar040.model"
    binary.parent.mkdir(parents=True)
    model.parent.mkdir(parents=True)
    binary.write_bytes(b"binary\n")
    model.write_bytes(b"model\n")
    model_hash = pkg.sha256_file(model)
    build = {
        "stage": "pr67_7633ebb",
        "commit": pkg.PR67_COMMIT,
        "binary": str(binary),
        "binary_sha256": pkg.sha256_file(binary),
        "compiler": "synthetic",
        "config_directory": str(config),
        "config_sha256": pkg.sha256_directory(config),
        "null_model": str(model),
        "null_model_sha256": model_hash,
    }
    build_path = build_root / "build_manifest.json"
    write_json(build_path, build)
    support_binary = root / "tools" / "summarize_consensus_support"
    support_binary.parent.mkdir(parents=True)
    support_binary.write_bytes(b"helper\n")
    support_binary_hash = pkg.sha256_file(support_binary)
    support_source = (
        Path(pkg.__file__).resolve().parents[1]
        / "brca100_netbid_qc"
        / "summarize_consensus_support.cpp"
    )
    support_source_hash = pkg.sha256_file(support_source)
    benchmark_root = Path(pkg.__file__).resolve().parents[1]
    wrapper = benchmark_root / "brca100_netbid_qc" / "netbid2-r"
    r_script = Path(pkg.__file__).resolve().with_name("run_netbid_qc.R")
    wrapper_hash = pkg.sha256_file(wrapper)
    r_script_hash = pkg.sha256_file(r_script)
    netbid_environment = {
        "R": "synthetic-R",
        "NetBID2": "synthetic-NetBID2",
        "NetBID2_remote_sha": "synthetic-remote-sha",
        "igraph": "synthetic-igraph",
    }

    design_points: list[dict[str, object]] = []
    for index, (key, probability) in enumerate(pkg.POINTS, 1):
        design_points.append(
            {
                "key": key,
                "label": key,
                "role": "synthetic",
                "p_token": format(probability, ".17g"),
                "p_value": probability,
                "mi_cutoff": 1.0 / index,
                "mi_cutoff_header": format(1.0 / index, ".6g"),
                "p_header": format(probability, ".6g"),
                "tail_extrapolated": False,
                "validation_class": "synthetic-validation",
                "validated_p_min": 1e-8,
                "validated_p_max": 1e-2,
            }
        )
    design = {
        "schema": "sjaracne-brca100-pr67-p-sweep-v2",
        "commit": pkg.PR67_COMMIT,
        "binary_sha256": build["binary_sha256"],
        "config_sha256": build["config_sha256"],
        "null_model_sha256": model_hash,
        "all_points": design_points,
        "fixed_parameters": {
            "sampling": "fixed 80% without replacement",
            "m": 80,
            "npar": 40,
            "dpi_epsilon": 0,
            "consensus_p": 1e-5,
            "seeds": list(pkg.SEEDS),
        },
        "inputs": inputs,
    }
    design_path = root / "sweep_design.json"
    write_json(design_path, design)
    legacy_design = {
        **design,
        "schema": "sjaracne-brca100-pr67-p-sweep-v1",
        "all_points": design_points[: len(pkg.LEGACY_POINT_KEYS)],
    }
    history_root = root / "sweep_design_history"
    legacy_placeholder = history_root / "legacy.sweep_design.json"
    write_json(legacy_placeholder, legacy_design)
    legacy_hash = pkg.sha256_file(legacy_placeholder)
    archive_path = history_root / f"{legacy_hash}.sweep_design.json"
    legacy_placeholder.replace(archive_path)
    design_hash = pkg.sha256_file(design_path)
    migration_path = (
        history_root / f"{legacy_hash}_to_{design_hash}.migration.json"
    )
    migration = {
        "schema": "sjaracne-brca100-pr67-p-sweep-design-migration-v1",
        "operation": "append-only-point-extension",
        "from": {
            "schema": "sjaracne-brca100-pr67-p-sweep-v1",
            "sha256": legacy_hash,
            "archived_path": archive_path.relative_to(root).as_posix(),
            "point_keys": list(pkg.LEGACY_POINT_KEYS),
        },
        "to": {
            "schema": "sjaracne-brca100-pr67-p-sweep-v2",
            "sha256": design_hash,
            "active_path": "sweep_design.json",
            "point_keys": list(pkg.POINT_KEYS),
        },
        "appended_points": design_points[len(pkg.LEGACY_POINT_KEYS) :],
        "manifest_path": migration_path.relative_to(root).as_posix(),
    }
    migration["fingerprint"] = pkg.json_fingerprint(migration)
    write_json(migration_path, migration)

    point_by_key: dict[str, dict[str, object]] = {}
    point_hashes: dict[str, str] = {}
    for point in design_points:
        key = str(point["key"])
        record = {
            "schema": "sjaracne-brca100-pr67-p-sweep-point-v1",
            **point,
            "commit": pkg.PR67_COMMIT,
            "binary_sha256": build["binary_sha256"],
            "config_sha256": build["config_sha256"],
            "null_model_sha256": model_hash,
            "sampling": "fixed 80% without replacement",
            "m": 80,
            "npar": 40,
            "dpi_epsilon": 0,
            "consensus_p": 1e-5,
            "seeds": list(pkg.SEEDS),
            "inputs": inputs,
        }
        path = root / "results" / key / "point_manifest.json"
        write_json(path, record)
        point_by_key[key] = record
        point_hashes[key] = pkg.sha256_file(path)

    run_rows: list[dict[str, object]] = []
    arm_partial: dict[tuple[str, str], dict[str, str]] = {}
    support_records: list[dict[str, object]] = []
    netbid_records: list[dict[str, object]] = []
    for key in pkg.POINT_KEYS:
        point = point_by_key[key]
        for driver in pkg.DRIVERS:
            arm = root / "results" / key / driver
            adjacency_hashes: list[str] = []
            adjacency_set: list[tuple[str, str]] = []
            metadata_set: list[tuple[str, str]] = []
            for seed in pkg.SEEDS:
                stem = f"TF_run_{seed:03d}"
                adjacency = arm / "adjacency" / f"{stem}.adj"
                adjacency.parent.mkdir(parents=True, exist_ok=True)
                adjacency.write_bytes(b"a")
                adjacency_hash = pkg.sha256_file(adjacency)
                adjacency_hashes.append(adjacency_hash)
                adjacency_set.append(
                    (adjacency.relative_to(arm).as_posix(), adjacency_hash)
                )
                log_root = arm / "logs"
                log_root.mkdir(parents=True, exist_ok=True)
                stdout = log_root / f"{stem}.stdout.log"
                stderr = log_root / f"{stem}.stderr.log"
                timing = log_root / f"{stem}.time.txt"
                stdout.write_bytes(b"")
                stderr.write_bytes(b"")
                timing.write_bytes(b"elapsed_s=0\n")
                metadata = {
                    "schema": "sjaracne-brca100-pr67-p-sweep-seed-v1",
                    "point": {"key": key},
                    "driver": driver,
                    "seed": seed,
                    "adjacency": {"full_sha256": adjacency_hash, "bytes": 1},
                    "stdout_sha256": pkg.sha256_file(stdout),
                    "stderr_sha256": pkg.sha256_file(stderr),
                }
                metadata_path = arm / "seed_metadata" / f"{stem}.json"
                write_json(metadata_path, metadata)
                metadata_set.append(
                    (metadata_path.relative_to(arm).as_posix(), pkg.sha256_file(metadata_path))
                )
                run_rows.append(
                    {
                        "point": key,
                        "p_value": point["p_value"],
                        "mi_cutoff": point["mi_cutoff"],
                        "validation_class": point["validation_class"],
                        "commit": pkg.PR67_COMMIT,
                        "driver": driver,
                        "seed": seed,
                        "binary_sha256": build["binary_sha256"],
                        "edges": 1,
                        "source_rows": 1,
                        "adjacency_bytes": 1,
                        "adjacency_sha256": adjacency_hash,
                        "data_sha256": hashlib.sha256(b"data").hexdigest(),
                        "stderr_bytes": 0,
                    }
                )

            consensus_root = arm / "consensus"
            consensus_root.mkdir(parents=True, exist_ok=True)
            consensus_files: dict[str, Path] = {}
            for filename in (
                "consensus_network_ncol_.txt",
                "consensus_network_3col_.txt",
                "parameter_info_.txt",
                "bootstrap_info_.txt",
            ):
                path = consensus_root / filename
                path.write_bytes(filename.encode("ascii"))
                consensus_files[filename] = path
            for filename in ("consensus.stdout.log", "consensus.stderr.log", "consensus.time.txt"):
                (arm / "logs" / filename).write_bytes(b"log\n")
            consensus_fingerprint = hashlib.sha256(f"consensus:{key}:{driver}".encode()).hexdigest()
            consensus_manifest = {
                "fingerprint": consensus_fingerprint,
                "stage": key,
                "driver": driver,
                "ncol": {
                    "sha256": pkg.sha256_file(consensus_files["consensus_network_ncol_.txt"]),
                    "bytes": consensus_files["consensus_network_ncol_.txt"].stat().st_size,
                    "edges": 1,
                },
                "consensus_3col_sha256": pkg.sha256_file(
                    consensus_files["consensus_network_3col_.txt"]
                ),
                "parameter_info_sha256": pkg.sha256_file(consensus_files["parameter_info_.txt"]),
                "bootstrap_info_sha256": pkg.sha256_file(consensus_files["bootstrap_info_.txt"]),
            }
            consensus_manifest_path = arm / "consensus_manifest.json"
            write_json(consensus_manifest_path, consensus_manifest)

            support_output = consensus_root / "consensus_support.tsv"
            support_output.write_bytes(b"support\n")
            for filename in (
                "support_summary.stdout.log",
                "support_summary.stderr.log",
                "support_summary.time.txt",
            ):
                (arm / "logs" / filename).write_bytes(b"log\n")
            support_fingerprint = hashlib.sha256(f"support:{key}:{driver}".encode()).hexdigest()
            support_manifest = {
                "schema": "sjaracne-brca100-pr67-p-sweep-support-v1",
                "fingerprint": support_fingerprint,
                "point": key,
                "p_value": point["p_value"],
                "driver": driver,
                "consensus_sha256": consensus_manifest["ncol"]["sha256"],
                "consensus_manifest_sha256": pkg.sha256_file(consensus_manifest_path),
                "point_manifest_sha256": point_hashes[key],
                "source_sha256": support_source_hash,
                "binary_sha256": support_binary_hash,
                "adjacency_sha256": adjacency_hashes,
                "output": str(support_output),
                "output_sha256": pkg.sha256_file(support_output),
            }
            support_manifest_path = arm / "support_summary_manifest.json"
            write_json(support_manifest_path, support_manifest)
            support_records.append(support_manifest)

            netbid_root = arm / "netbid2_qc"
            netbid_root.mkdir(parents=True)
            for filename in (
                "driver_target_sizes.tsv",
                "netbid_environment.tsv",
                "network_summary.tsv",
            ):
                (netbid_root / filename).write_bytes((filename + "\n").encode("ascii"))
            netbid_stdout = arm / "logs" / "netbid2_qc.stdout.log"
            netbid_stderr = arm / "logs" / "netbid2_qc.stderr.log"
            netbid_stdout.write_bytes(b"")
            netbid_stderr.write_bytes(b"")
            driver_filename, prefix = pkg.NETBID_DRIVERS[driver]
            driver_path = root / "inputs" / driver_filename
            netbid_payload = {
                "schema": "sjaracne-brca100-pr67-p-sweep-netbid2-v1",
                "mode": "summary",
                "point": key,
                "p_value": point["p_value"],
                "mi_cutoff": point["mi_cutoff"],
                "point_manifest_sha256": point_hashes[key],
                "sweep_design_sha256": pkg.sha256_file(design_path),
                "driver": driver,
                "driver_sha256": pkg.sha256_file(driver_path),
                "consensus_sha256": consensus_manifest["ncol"]["sha256"],
                "consensus_manifest_sha256": pkg.sha256_file(consensus_manifest_path),
                "r_script_sha256": r_script_hash,
                "wrapper_sha256": wrapper_hash,
                "environment": netbid_environment,
                "prefix": prefix,
            }
            netbid_fingerprint = pkg.json_fingerprint(netbid_payload)
            netbid_manifest = {
                **netbid_payload,
                "fingerprint": netbid_fingerprint,
                "command": [
                    str(wrapper), "Rscript", str(r_script),
                    str(consensus_files["consensus_network_ncol_.txt"]),
                    str(driver_path), str(arm / "netbid2_qc.partial"), prefix, "false",
                ],
                "finished_at_utc": "synthetic",
                "output": str(netbid_root),
                "output_inventory": inventory(netbid_root),
                "stdout_sha256": pkg.sha256_file(netbid_stdout),
                "stderr_sha256": pkg.sha256_file(netbid_stderr),
                "stderr_bytes": 0,
            }
            netbid_manifest_path = arm / "netbid2_qc_manifest.json"
            write_json(netbid_manifest_path, netbid_manifest)
            netbid_records.append(netbid_manifest)

            netbid_set = pkg.set_digest(
                [(item["path"], item["sha256"]) for item in netbid_manifest["output_inventory"]]
            )
            arm_partial[(key, driver)] = {
                "adjacency_set_sha256": pkg.set_digest(adjacency_set),
                "seed_metadata_set_sha256": pkg.set_digest(metadata_set),
                "run_manifest_rows_sha256": hashlib.sha256(f"rows:{key}:{driver}".encode()).hexdigest(),
                "point_manifest_sha256": point_hashes[key],
                "consensus_sha256": consensus_manifest["ncol"]["sha256"],
                "consensus_manifest_sha256": pkg.sha256_file(consensus_manifest_path),
                "consensus_fingerprint": consensus_fingerprint,
                "support_sha256": support_manifest["output_sha256"],
                "support_manifest_sha256": pkg.sha256_file(support_manifest_path),
                "support_fingerprint": support_fingerprint,
                "support_source_sha256": support_source_hash,
                "support_binary_sha256": support_binary_hash,
                "netbid2_qc_set_sha256": netbid_set,
                "netbid2_network_summary_sha256": pkg.sha256_file(
                    netbid_root / "network_summary.tsv"
                ),
                "netbid2_manifest_sha256": pkg.sha256_file(netbid_manifest_path),
                "netbid2_fingerprint": netbid_fingerprint,
            }

    run_manifest = root / "results" / "run_manifest.tsv"
    with run_manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(run_rows[0]))
        writer.writeheader()
        writer.writerows(run_rows)

    support_aggregate = {
        "schema": "sjaracne-brca100-pr67-p-sweep-support-aggregate-v1",
        "sweep_design_sha256": pkg.sha256_file(design_path),
        "points": list(pkg.POINT_KEYS),
        "drivers": list(pkg.DRIVERS),
        "records": support_records,
    }
    support_aggregate_path = root / "results" / "support_summary_manifest.json"
    write_json(support_aggregate_path, support_aggregate)
    netbid_aggregate = finalized(
        {
            "schema": "sjaracne-brca100-pr67-p-sweep-netbid2-summary-aggregate-v1",
            "environment": netbid_environment,
            "sweep_design_sha256": pkg.sha256_file(design_path),
            "all_sweep_points": list(pkg.POINT_KEYS),
            "selection": {"points": list(pkg.POINT_KEYS), "drivers": list(pkg.DRIVERS)},
            "summary_runs": netbid_records,
        }
    )
    netbid_aggregate_path = root / "results" / "netbid2_qc_manifest.json"
    write_json(netbid_aggregate_path, netbid_aggregate)
    html_point = "p5e-05"
    html_records: list[dict[str, object]] = []
    for driver in pkg.DRIVERS:
        arm = root / "results" / html_point / driver
        point = point_by_key[html_point]
        driver_filename, prefix = pkg.NETBID_DRIVERS[driver]
        driver_path = root / "inputs" / driver_filename
        consensus_path = arm / "consensus" / "consensus_network_ncol_.txt"
        consensus_manifest_path = arm / "consensus_manifest.json"
        html_root = arm / "netbid2_qc_html"
        html_root.mkdir()
        for filename in (
            "driver_target_sizes.tsv",
            "netbid_environment.tsv",
            "network_summary.tsv",
            f"{prefix}netQC.Rmd",
            f"{prefix}netQC.html",
        ):
            (html_root / filename).write_bytes((filename + "\n").encode("ascii"))
        stdout = arm / "logs" / "netbid2_qc_html.stdout.log"
        stderr = arm / "logs" / "netbid2_qc_html.stderr.log"
        stdout.write_bytes(b"")
        stderr.write_bytes(b"")
        html_payload = {
            "schema": "sjaracne-brca100-pr67-p-sweep-netbid2-v1",
            "mode": "html",
            "point": html_point,
            "p_value": point["p_value"],
            "mi_cutoff": point["mi_cutoff"],
            "point_manifest_sha256": point_hashes[html_point],
            "sweep_design_sha256": pkg.sha256_file(design_path),
            "driver": driver,
            "driver_sha256": pkg.sha256_file(driver_path),
            "consensus_sha256": pkg.sha256_file(consensus_path),
            "consensus_manifest_sha256": pkg.sha256_file(consensus_manifest_path),
            "r_script_sha256": r_script_hash,
            "wrapper_sha256": wrapper_hash,
            "environment": netbid_environment,
            "prefix": prefix,
        }
        record = {
            **html_payload,
            "fingerprint": pkg.json_fingerprint(html_payload),
            "command": [
                str(wrapper), "Rscript", str(r_script), str(consensus_path),
                str(driver_path), str(arm / "netbid2_qc_html.partial"), prefix, "true",
            ],
            "finished_at_utc": "synthetic",
            "output": str(html_root),
            "output_inventory": inventory(html_root),
            "stdout_sha256": pkg.sha256_file(stdout),
            "stderr_sha256": pkg.sha256_file(stderr),
            "stderr_bytes": 0,
        }
        write_json(arm / "netbid2_qc_html_manifest.json", record)
        html_records.append(record)
    html_aggregate = finalized(
        {
            "schema": "sjaracne-brca100-pr67-p-sweep-netbid2-html-aggregate-v1",
            "environment": netbid_environment,
            "sweep_design_sha256": pkg.sha256_file(design_path),
            "all_sweep_points": list(pkg.POINT_KEYS),
            "selection": {
                "points": list(pkg.POINT_KEYS),
                "drivers": list(pkg.DRIVERS),
                "html_points": [html_point],
            },
            "html_runs": html_records,
        }
    )
    write_json(root / "results" / "netbid2_qc_html_manifest.json", html_aggregate)

    invocations = {
        "schema": "sjaracne-brca100-pr67-p-sweep-invocations-v1",
        "invocations": [
            {
                "status": "complete", "phase": "infer", "points": list(pkg.POINT_KEYS),
                "drivers": list(pkg.DRIVERS), "seed_start": 1, "seed_end": 100,
                "workers": 2, "finished_at_utc": "synthetic",
                "inference_jobs": len(pkg.POINT_KEYS) * len(pkg.DRIVERS) * len(pkg.SEEDS),
                "inference_new_jobs": len(pkg.POINT_KEYS) * len(pkg.DRIVERS) * len(pkg.SEEDS),
                "inference_resumed_jobs": 0,
            },
            {
                "status": "complete", "phase": "consensus", "points": list(pkg.POINT_KEYS),
                "drivers": list(pkg.DRIVERS), "seed_start": 1, "seed_end": 100,
                "workers": 1, "finished_at_utc": "synthetic",
            },
        ],
    }
    write_json(root / "invocations.json", invocations)

    validation = root / "results" / "validation"
    table = validation / "anchor_seed_equivalence.tsv"
    validation.mkdir(parents=True)
    fields = [
        "sweep_point", "prior_stage", "driver", "seed", "data_sha256", "edges",
        "source_rows", "sweep_metadata_sha256", "prior_metadata_sha256",
    ]
    with table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for point, stage in (
            ("p1e-07", "pr67_7633ebb"),
            ("p_pr66_cutoff_match", "pr66_5809183"),
        ):
            for driver in pkg.DRIVERS:
                for seed in pkg.SEEDS:
                    writer.writerow(
                        {
                            "sweep_point": point,
                            "prior_stage": stage,
                            "driver": driver,
                            "seed": seed,
                            "data_sha256": hashlib.sha256(b"data").hexdigest(),
                            "edges": 1,
                            "source_rows": 1,
                            "sweep_metadata_sha256": hashlib.sha256(b"sweep").hexdigest(),
                            "prior_metadata_sha256": hashlib.sha256(b"prior").hexdigest(),
                        }
                    )
    anchor_manifest = {
        "schema": "sjaracne-brca100-pr67-p-sweep-anchor-equivalence-v1",
        "drivers": list(pkg.DRIVERS),
        "seeds": list(pkg.SEEDS),
        "comparisons": 2 * len(pkg.DRIVERS) * len(pkg.SEEDS),
        "all_data_sections_equal": True,
        "table": str(table),
        "table_sha256": pkg.sha256_file(table),
        "sweep_design_sha256": pkg.sha256_file(design_path),
    }
    anchor_manifest_path = validation / "anchor_seed_equivalence_manifest.json"
    write_json(anchor_manifest_path, anchor_manifest)

    analysis_root = root / "results" / "analysis"
    for relative in sorted(pkg.ANALYSIS_OUTPUTS):
        path = analysis_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"{}\n" if relative == "selection.json" else b"analysis\n")
    output_files = {
        relative: {
            "bytes": (analysis_root / relative).stat().st_size,
            "sha256": pkg.sha256_file(analysis_root / relative),
        }
        for relative in sorted(pkg.ANALYSIS_OUTPUTS)
    }
    analysis_manifest = {
        "schema": "sjaracne-brca100-pr67-threshold-sweep-analysis-v1",
        "work_root": str(root),
        "output_root": str(analysis_root),
        "p_keys_in_increasing_p_order": list(pkg.POINT_KEYS),
        "sweep_design_sha256": pkg.sha256_file(design_path),
        "run_manifest_sha256": pkg.sha256_file(run_manifest),
        "support_aggregate_manifest_sha256": pkg.sha256_file(support_aggregate_path),
        "netbid2_aggregate_manifest_sha256": pkg.sha256_file(netbid_aggregate_path),
        "build": {
            "commit": pkg.PR67_COMMIT,
            "binary_sha256": build["binary_sha256"],
            "config_sha256": build["config_sha256"],
            "null_model_sha256": model_hash,
            "build_manifest_sha256": pkg.sha256_file(build_path),
        },
        "matched_seed_count": len(pkg.SEEDS),
        "arms": {f"{key}/{driver}": arm_partial[(key, driver)] for key in pkg.POINT_KEYS for driver in pkg.DRIVERS},
        "output_files": output_files,
        "operating_point_selection": {},
        "pr66_context": {
            "cutoff_match_consensus_exact": True,
            "anchor_seed_equivalence": {
                "manifest_sha256": pkg.sha256_file(anchor_manifest_path),
                "table_sha256": pkg.sha256_file(table),
                "comparisons": 2 * len(pkg.DRIVERS) * len(pkg.SEEDS),
                "all_data_sections_equal": True,
            },
        },
    }
    write_json(analysis_root / "analysis_manifest.json", analysis_manifest)
    return expected_input_hashes, model_hash


def replace_html_record(
    root: Path, key: str, driver: str, record: dict[str, object]
) -> None:
    """Keep synthetic per-arm/root HTML records coherent after a test mutation."""
    manifest_path = root / "results" / key / driver / "netbid2_qc_html_manifest.json"
    write_json(manifest_path, record)
    aggregate_path = root / "results" / "netbid2_qc_html_manifest.json"
    aggregate = json.loads(aggregate_path.read_text(encoding="utf-8"))
    replaced = 0
    for index, candidate in enumerate(aggregate["html_runs"]):
        if candidate.get("point") == key and candidate.get("driver") == driver:
            aggregate["html_runs"][index] = record
            replaced += 1
    if replaced != 1:
        raise AssertionError("Synthetic HTML aggregate did not contain exactly one arm")
    aggregate.pop("fingerprint", None)
    aggregate["fingerprint"] = pkg.json_fingerprint(aggregate)
    write_json(aggregate_path, aggregate)


class PackageResultsTest(unittest.TestCase):
    def test_full_synthetic_package_and_tamper_rejection(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            work = base / "work"
            work.mkdir()
            output = base / "package"
            with mock.patch.object(pkg, "SEEDS", (1, 2)):
                input_hashes, model_hash = make_fixture(work)
                patches = (
                    mock.patch.object(pkg, "EXPECTED_INPUT_SHA256", input_hashes),
                    mock.patch.object(pkg, "NULL_MODEL_SHA256", model_hash),
                )
                with patches[0], patches[1]:
                    pkg.package_results(work, output)
                    self.assertTrue((output / "package_manifest.json").is_file())
                    self.assertEqual(
                        (output / ".gitattributes").read_bytes(), b"* -text -whitespace\n"
                    )
                    self.assertFalse(any(output.rglob("*.adj")))
                    self.assertFalse(any(path.name.endswith("netQC.html") for path in output.rglob("*")))
                    self.assertFalse(any(path.name.endswith("netQC.Rmd") for path in output.rglob("*")))
                    self.assertTrue(
                        (
                            output
                            / "provenance"
                            / "aggregates"
                            / "netbid2_qc_html_manifest.json"
                        ).is_file()
                    )
                    self.assertEqual(
                        len(
                            list(
                                (
                                    output
                                    / "provenance"
                                    / "sweep_design_history"
                                ).glob("*.migration.json")
                            )
                        ),
                        1,
                    )
                    checksum_lines = (output / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
                    nonself = [
                        path for path in output.rglob("*")
                        if path.is_file() and path.name != "SHA256SUMS"
                    ]
                    self.assertEqual(len(checksum_lines), len(nonself))
                    omitted = json.loads((output / "omitted_artifacts.json").read_text(encoding="utf-8"))
                    self.assertEqual(
                        omitted["category_counts"]["seed-adjacency"],
                        len(pkg.POINT_KEYS) * len(pkg.DRIVERS) * len(pkg.SEEDS),
                    )

                    archive = next(
                        (work / "sweep_design_history").glob("*.sweep_design.json")
                    )
                    archive_bytes = archive.read_bytes()
                    archive.write_bytes(archive_bytes + b" ")
                    rejected_history = base / "rejected-history"
                    with self.assertRaisesRegex(
                        ValueError, "append-only sweep-design migration"
                    ):
                        pkg.package_results(work, rejected_history)
                    self.assertFalse(rejected_history.exists())
                    archive.write_bytes(archive_bytes)

                    (work / "results" / "analysis" / "network_summary.tsv").write_bytes(b"tampered\n")
                    rejected = base / "rejected"
                    with self.assertRaisesRegex(ValueError, "Analysis output"):
                        pkg.package_results(work, rejected)
                    self.assertFalse(rejected.exists())
                    self.assertFalse(rejected.with_name("rejected.partial").exists())

    def test_same_size_adjacency_and_consensus_support_tampers_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            work = base / "work"
            work.mkdir()
            with mock.patch.object(pkg, "SEEDS", (1, 2)):
                input_hashes, model_hash = make_fixture(work)
                with (
                    mock.patch.object(pkg, "EXPECTED_INPUT_SHA256", input_hashes),
                    mock.patch.object(pkg, "NULL_MODEL_SHA256", model_hash),
                ):
                    targets = (
                        work
                        / "results"
                        / "p1e-07"
                        / "tf"
                        / "adjacency"
                        / "TF_run_001.adj",
                        work
                        / "results"
                        / "p1e-07"
                        / "tf"
                        / "consensus"
                        / "consensus_network_3col_.txt",
                        work
                        / "results"
                        / "p1e-07"
                        / "tf"
                        / "consensus"
                        / "consensus_support.tsv",
                    )
                    for index, target in enumerate(targets):
                        with self.subTest(target=target.name):
                            original = target.read_bytes()
                            replacement = bytes([original[0] ^ 1]) + original[1:]
                            self.assertEqual(len(replacement), len(original))
                            target.write_bytes(replacement)
                            rejected = base / f"rejected-omitted-{index}"
                            with self.assertRaisesRegex(
                                ValueError, "Omitted artifact hash mismatch"
                            ):
                                pkg.package_results(work, rejected)
                            self.assertFalse(rejected.exists())
                            self.assertFalse(
                                rejected.with_name(rejected.name + ".partial").exists()
                            )
                            target.write_bytes(original)

    def test_optional_html_requires_complete_producer_records(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            work = base / "work"
            work.mkdir()
            key, driver = "p5e-05", "tf"
            manifest_path = (
                work / "results" / key / driver / "netbid2_qc_html_manifest.json"
            )
            aggregate_path = work / "results" / "netbid2_qc_html_manifest.json"
            with mock.patch.object(pkg, "SEEDS", (1, 2)):
                input_hashes, model_hash = make_fixture(work)
                original_manifest = manifest_path.read_bytes()
                original_aggregate = aggregate_path.read_bytes()
                with (
                    mock.patch.object(pkg, "EXPECTED_INPUT_SHA256", input_hashes),
                    mock.patch.object(pkg, "NULL_MODEL_SHA256", model_hash),
                ):
                    for label, mutate in (
                        ("arbitrary", lambda record: record.__setitem__("arbitrary", "field")),
                        ("truncated", lambda record: record.pop("driver_sha256")),
                    ):
                        with self.subTest(label=label):
                            record = json.loads(original_manifest.decode("utf-8"))
                            mutate(record)
                            replace_html_record(work, key, driver, record)
                            rejected = base / f"rejected-html-{label}"
                            with self.assertRaisesRegex(
                                ValueError, "missing/arbitrary fields"
                            ):
                                pkg.package_results(work, rejected)
                            self.assertFalse(rejected.exists())
                            manifest_path.write_bytes(original_manifest)
                            aggregate_path.write_bytes(original_aggregate)

                    aggregate = json.loads(original_aggregate.decode("utf-8"))
                    aggregate["html_runs"] = aggregate["html_runs"][:-1]
                    aggregate.pop("fingerprint")
                    aggregate["fingerprint"] = pkg.json_fingerprint(aggregate)
                    write_json(aggregate_path, aggregate)
                    rejected = base / "rejected-html-incomplete-aggregate"
                    with self.assertRaisesRegex(
                        ValueError, "arm ordering/coverage"
                    ):
                        pkg.package_results(work, rejected)
                    self.assertFalse(rejected.exists())

    def test_optional_html_shared_tsvs_must_equal_stable_summary(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            work = base / "work"
            work.mkdir()
            key, driver = "p5e-05", "tf"
            with mock.patch.object(pkg, "SEEDS", (1, 2)):
                input_hashes, model_hash = make_fixture(work)
                html_root = work / "results" / key / driver / "netbid2_qc_html"
                shared = html_root / "network_summary.tsv"
                original = shared.read_bytes()
                shared.write_bytes(bytes([original[0] ^ 1]) + original[1:])
                manifest_path = (
                    work / "results" / key / driver / "netbid2_qc_html_manifest.json"
                )
                record = json.loads(manifest_path.read_text(encoding="utf-8"))
                record["output_inventory"] = inventory(html_root)
                replace_html_record(work, key, driver, record)
                with (
                    mock.patch.object(pkg, "EXPECTED_INPUT_SHA256", input_hashes),
                    mock.patch.object(pkg, "NULL_MODEL_SHA256", model_hash),
                ):
                    rejected = base / "rejected-html-shared-tsv"
                    with self.assertRaisesRegex(ValueError, "shared TSV mismatch"):
                        pkg.package_results(work, rejected)
                    self.assertFalse(rejected.exists())


if __name__ == "__main__":
    unittest.main()
