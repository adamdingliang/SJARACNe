#!/usr/bin/env python3
"""Create a compact, checksummed package of a completed recurrence sweep."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
from typing import Any


DRIVERS = ("tf", "sig")
MINIMUM_SUPPORTS = tuple(range(6, 21))
PACKAGE_SCHEMA = "sjaracne-brca100-consensus-recurrence-package-v1"
ARMS_FOR_PACKAGE = {
    "tf": {"point": "p1e-03", "driver_file": "BRCA100_TF.txt"},
    "sig": {"point": "p5e-04", "driver_file": "BRCA100_SIG.txt"},
}
EXPECTED_NETBID_ENVIRONMENT = {
    "R": "R version 4.4.3 (2025-02-28)",
    "NetBID2": "2.2.0",
    "NetBID2_remote_sha": "5defa454d600b94f5dd6d1f9f4428f99759a6821",
    "igraph": "2.3.3",
}
SOURCE_NETBID_SCHEMA = "sjaracne-brca100-pr67-p-sweep-netbid2-v1"
SOURCE_NETBID_FINGERPRINT_FIELDS = frozenset(
    {
        "schema",
        "mode",
        "point",
        "p_value",
        "mi_cutoff",
        "point_manifest_sha256",
        "sweep_design_sha256",
        "driver",
        "driver_sha256",
        "consensus_sha256",
        "consensus_manifest_sha256",
        "r_script_sha256",
        "wrapper_sha256",
        "environment",
        "prefix",
    }
)
SOURCE_NETBID_OUTPUTS = {
    "driver_target_sizes.tsv",
    "netbid_environment.tsv",
    "network_summary.tsv",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def fingerprint(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def serialized_json(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def copy_verified(
    source: Path,
    destination_root: Path,
    relative: Path,
    copied: list[dict[str, object]],
) -> None:
    if not source.is_file():
        raise ValueError(f"Missing package source: {source}")
    destination = destination_root / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    source_hash = sha256_file(source)
    if sha256_file(destination) != source_hash:
        raise ValueError(f"Package copy mismatch: {source}")
    copied.append(
        {
            "path": relative.as_posix(),
            "bytes": source.stat().st_size,
            "sha256": source_hash,
        }
    )


def record_omitted(
    path: Path,
    work_root: Path,
    expected_hash: str,
    omitted: list[dict[str, object]],
    recorded_path: str | None = None,
) -> None:
    if not path.is_file():
        raise ValueError(f"Missing omitted artifact: {path}")
    actual_hash = sha256_file(path)
    if actual_hash != expected_hash:
        raise ValueError(f"Omitted artifact hash mismatch: {path}")
    omitted.append(
        {
            "path": recorded_path or path.relative_to(work_root).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": actual_hash,
        }
    )


def write_sha256s(root: Path) -> None:
    paths = sorted(
        path for path in root.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS"
    )
    lines = [f"{sha256_file(path)}  {path.relative_to(root).as_posix()}\n" for path in paths]
    (root / "SHA256SUMS").write_text("".join(lines), encoding="utf-8", newline="\n")


def read_netbid_environment(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["component", "version"]:
            raise ValueError(f"Unexpected NetBID2 environment header: {path}")
        result: dict[str, str] = {}
        for row in reader:
            component = row["component"]
            if component in result:
                raise ValueError(f"Duplicate NetBID2 environment component: {path}")
            result[component] = row["version"]
    return result


def package(work_root: Path, output_root: Path, results_markdown: Path | None) -> None:
    if output_root.exists():
        raise ValueError(f"Output root already exists: {output_root}")
    if work_root in output_root.parents or output_root in work_root.parents:
        raise ValueError("Package root and live work root must be separate")
    partial_root = output_root.with_name(output_root.name + ".partial")
    owner_path = output_root.with_name(output_root.name + ".partial.owner.json")
    if partial_root.exists():
        if not owner_path.is_file():
            raise ValueError(f"Unrecognized partial output root: {partial_root}")
        owner = load_json(owner_path)
        expected_owner_fields = {
            "schema": "sjaracne-brca100-consensus-recurrence-package-partial-v1",
            "work_root": str(work_root),
            "output_root": str(output_root),
        }
        if any(owner.get(key) != value for key, value in expected_owner_fields.items()):
            raise ValueError(f"Incompatible partial output root: {partial_root}")
        owner_pid = int(owner.get("pid", -1))
        if owner_pid <= 0:
            raise ValueError(f"Invalid partial-package owner PID: {owner_pid}")
        try:
            os.kill(owner_pid, 0)
        except ProcessLookupError:
            pass
        except PermissionError as error:
            raise ValueError(f"Cannot verify partial-package owner PID {owner_pid}") from error
        else:
            raise ValueError(
                f"Another package process still owns {partial_root} (PID {owner_pid})"
            )
        shutil.rmtree(partial_root)
        owner_path.unlink()
    elif owner_path.exists():
        owner = load_json(owner_path)
        expected_owner_fields = {
            "schema": "sjaracne-brca100-consensus-recurrence-package-partial-v1",
            "work_root": str(work_root),
            "output_root": str(output_root),
        }
        if any(owner.get(key) != value for key, value in expected_owner_fields.items()):
            raise ValueError(f"Incompatible orphan package-owner record: {owner_path}")
        owner_pid = int(owner.get("pid", -1))
        if owner_pid <= 0:
            raise ValueError(f"Invalid orphan package-owner PID: {owner_pid}")
        try:
            os.kill(owner_pid, 0)
        except ProcessLookupError:
            owner_path.unlink()
        except PermissionError as error:
            raise ValueError(f"Cannot verify orphan owner PID {owner_pid}") from error
        else:
            raise ValueError(
                f"Package owner PID {owner_pid} is active without its partial root"
            )
    copied: list[dict[str, object]] = []
    omitted: list[dict[str, object]] = []
    owner_payload = {
        "schema": "sjaracne-brca100-consensus-recurrence-package-partial-v1",
        "work_root": str(work_root),
        "output_root": str(output_root),
        "pid": os.getpid(),
    }
    with owner_path.open("xb") as handle:
        handle.write(serialized_json(owner_payload))
        handle.flush()
        os.fsync(handle.fileno())
    try:
        partial_root.mkdir(parents=True)
    except Exception:
        owner_path.unlink(missing_ok=True)
        raise
    (partial_root / ".gitattributes").write_text(
        "* -text -whitespace\n", encoding="utf-8", newline="\n"
    )
    design = load_json(work_root / "design.json")
    if design.get("schema") != "sjaracne-brca100-consensus-recurrence-sweep-v1":
        raise ValueError("Unexpected recurrence design schema")
    if design.get("fingerprint") != fingerprint(
        {key: value for key, value in design.items() if key != "fingerprint"}
    ):
        raise ValueError("Recurrence design fingerprint mismatch")
    repo_root = Path(__file__).resolve().parents[2]
    script_sources = {
        "run_recurrence_sweep.py": (
            repo_root
            / "benchmarks"
            / "brca100_consensus_recurrence_sweep"
            / "run_recurrence_sweep.py"
        ),
        "aggregate_recurrence.cpp": (
            repo_root
            / "benchmarks"
            / "brca100_consensus_recurrence_sweep"
            / "aggregate_recurrence.cpp"
        ),
        "package_results.py": Path(__file__).resolve(),
        "run_netbid_qc.R": (
            repo_root / "benchmarks" / "brca100_pr67_threshold_sweep" / "run_netbid_qc.R"
        ),
        "netbid2-r": repo_root / "benchmarks" / "brca100_netbid_qc" / "netbid2-r",
    }
    expected_script_hashes = design["benchmark_scripts"]
    for field, filename in (
        ("runner_sha256", "run_recurrence_sweep.py"),
        ("aggregator_source_sha256", "aggregate_recurrence.cpp"),
        ("netbid_r_sha256", "run_netbid_qc.R"),
    ):
        if sha256_file(script_sources[filename]) != expected_script_hashes[field]:
            raise ValueError(f"Benchmark script no longer matches immutable design: {filename}")
    for filename, source in script_sources.items():
        copy_verified(
            source, partial_root, Path("provenance/scripts") / filename, copied
        )
    copy_verified(
        work_root / "design.json", partial_root, Path("provenance/design.json"), copied
    )
    for field, filename in (
        ("source_sweep_design", "sweep_design.json"),
        ("source_run_manifest", "run_manifest.tsv"),
    ):
        record = design[field]
        source = Path(record["path"])
        if sha256_file(source) != record["sha256"]:
            raise ValueError(f"Upstream source provenance mismatch: {source}")
        copy_verified(
            source,
            partial_root,
            Path("provenance/source_sweep") / filename,
            copied,
        )
    build_manifest_path = work_root / "bin/build_manifest.json"
    build_manifest = load_json(build_manifest_path)
    if build_manifest.get("source_sha256") != expected_script_hashes[
        "aggregator_source_sha256"
    ]:
        raise ValueError("Aggregator build/source linkage mismatch")
    build_core = {
        "source": build_manifest["source"],
        "source_sha256": build_manifest["source_sha256"],
        "compiler": build_manifest["compiler"],
        "flags": build_manifest["flags"],
    }
    if build_manifest.get("fingerprint") != fingerprint(build_core):
        raise ValueError("Aggregator build fingerprint mismatch")
    record_omitted(
        Path(build_manifest["binary"]),
        work_root,
        build_manifest["binary_sha256"],
        omitted,
    )
    copy_verified(
        build_manifest_path,
        partial_root,
        Path("provenance/build_manifest.json"),
        copied,
    )
    analysis_manifest = load_json(work_root / "analysis/analysis_manifest.json")
    if analysis_manifest.get("design_sha256") != sha256_file(work_root / "design.json"):
        raise ValueError("Analysis/design linkage mismatch")
    if analysis_manifest.get("fingerprint") != fingerprint(
        {key: value for key, value in analysis_manifest.items() if key != "fingerprint"}
    ):
        raise ValueError("Analysis fingerprint mismatch")
    analysis_files = {
        "analysis_manifest.json": None,
        "network_summary.tsv": analysis_manifest["network_summary_sha256"],
        "driver_target_coverage.tsv": analysis_manifest[
            "driver_target_coverage_sha256"
        ],
        "plots/recurrence_density_coverage.png": analysis_manifest["plot_png_sha256"],
        "plots/recurrence_density_coverage.svg": analysis_manifest["plot_svg_sha256"],
    }
    for filename, expected_hash in analysis_files.items():
        source = work_root / "analysis" / filename
        if expected_hash is not None and sha256_file(source) != expected_hash:
            raise ValueError(f"Analysis hash mismatch: {source}")
        copy_verified(source, partial_root, Path("analysis") / filename, copied)
    for driver in DRIVERS:
        aggregate_root = work_root / "aggregate" / driver
        aggregate_manifest = load_json(aggregate_root / "aggregate_manifest.json")
        expected_aggregate_input = {
            "design_fingerprint": design["fingerprint"],
            "build_fingerprint": build_manifest["fingerprint"],
            "driver": driver,
            "source_point": ARMS_FOR_PACKAGE[driver]["point"],
        }
        if aggregate_manifest.get("schema") != "sjaracne-brca100-consensus-recurrence-arm-v1":
            raise ValueError(f"Unexpected aggregate schema for {driver}")
        if aggregate_manifest.get("input") != expected_aggregate_input:
            raise ValueError(f"Aggregate/design linkage mismatch for {driver}")
        if aggregate_manifest.get("k9_anchor_reproduced") is not True:
            raise ValueError(f"K=9 anchor was not reproduced for {driver}")
        aggregate_files = {
            "aggregate_manifest.json": None,
            "aggregate_summary.tsv": aggregate_manifest["summary_sha256"],
            "run_counts.tsv": aggregate_manifest["run_counts_sha256"],
            "plugin_tail.tsv": aggregate_manifest["plugin_tail_sha256"],
            "aggregate.stdout.log": aggregate_manifest["stdout_sha256"],
            "aggregate.stderr.log": aggregate_manifest["stderr_sha256"],
        }
        for filename, expected_hash in aggregate_files.items():
            source = aggregate_root / filename
            if expected_hash is not None and sha256_file(source) != expected_hash:
                raise ValueError(f"Aggregate manifest hash mismatch: {source}")
            copy_verified(
                source,
                partial_root,
                Path("provenance/aggregate") / driver / filename,
                copied,
            )
        record_omitted(
            aggregate_root / "recurrence_edges.tsv",
            work_root,
            aggregate_manifest["edge_sha256"],
            omitted,
        )
        enhanced_root = aggregate_root / "enhanced_minimum_k006"
        enhanced_manifest = load_json(enhanced_root / "enhanced_manifest.json")
        expected_enhanced_input = {
            "recurrence_sha256": aggregate_manifest["edge_sha256"],
            "expression_sha256": design["inputs"]["BRCA100.exp"]["sha256"],
            "enhancer_sha256": sha256_file(
                repo_root / "SJARACNe" / "bin" / "create_consensus_network.py"
            ),
            "minimum_support": 6,
        }
        if enhanced_manifest.get("schema") != (
            "sjaracne-brca100-consensus-recurrence-enhanced-base-v1"
        ):
            raise ValueError(f"Unexpected enhanced-network schema for {driver}")
        if enhanced_manifest.get("input") != expected_enhanced_input:
            raise ValueError(f"Enhanced/aggregate linkage mismatch for {driver}")
        if enhanced_manifest.get("rows") != aggregate_manifest[
            "edge_counts_by_minimum_support"
        ]["6"]:
            raise ValueError(f"Enhanced-network row count mismatch for {driver}")
        copy_verified(
            enhanced_root / "enhanced_manifest.json",
            partial_root,
            Path("provenance/aggregate") / driver / "enhanced_manifest.json",
            copied,
        )
        for filename, field in (
            ("consensus_network_3col_.txt", "three_col_sha256"),
            ("consensus_network_ncol_.txt", "ncol_sha256"),
        ):
            record_omitted(
                enhanced_root / filename,
                work_root,
                enhanced_manifest[field],
                omitted,
            )
        with (aggregate_root / "plugin_tail.tsv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            plugin_tail_rows = {
                int(row["minimum_support"]): row
                for row in csv.DictReader(handle, delimiter="\t")
            }
        if set(plugin_tail_rows) != set(MINIMUM_SUPPORTS):
            raise ValueError(f"Incomplete plug-in tail table for {driver}")
        for threshold in MINIMUM_SUPPORTS:
            arm_root = work_root / "results" / driver / f"k{threshold:03d}"
            network_manifest = load_json(arm_root / "network_manifest.json")
            if network_manifest.get("schema") != (
                "sjaracne-brca100-consensus-recurrence-network-v1"
            ):
                raise ValueError(f"Unexpected threshold-network schema: {arm_root}")
            if network_manifest.get("driver") != driver:
                raise ValueError(f"Threshold-network driver mismatch: {arm_root}")
            if network_manifest.get("source_point") != ARMS_FOR_PACKAGE[driver]["point"]:
                raise ValueError(f"Threshold-network point mismatch: {arm_root}")
            if float(network_manifest.get("per_subsample_p")) != float(
                design["arms"][driver]["per_subsample_p"]
            ):
                raise ValueError(f"Threshold-network p-value mismatch: {arm_root}")
            if network_manifest.get("minimum_support") != threshold:
                raise ValueError(f"Threshold manifest mismatch: {arm_root}")
            if float(network_manifest.get("support_fraction")) != threshold / 100:
                raise ValueError(f"Threshold support-fraction mismatch: {arm_root}")
            if network_manifest.get("edges") != aggregate_manifest[
                "edge_counts_by_minimum_support"
            ][str(threshold)]:
                raise ValueError(f"Threshold edge-count mismatch: {arm_root}")
            if network_manifest.get("aggregate_manifest_sha256") != sha256_file(
                aggregate_root / "aggregate_manifest.json"
            ):
                raise ValueError(f"Threshold/aggregate linkage mismatch: {arm_root}")
            tail_row = plugin_tail_rows[threshold]
            tail_checks = {
                "plugin_poisson_binomial_tail": "plugin_poisson_binomial_tail",
                "legacy_normal_tail": "legacy_normal_tail",
                "plugin_null_edge_burden_proxy": "plugin_null_edge_burden_proxy",
            }
            for manifest_field, table_field in tail_checks.items():
                if not math.isclose(
                    float(network_manifest[manifest_field]),
                    float(tail_row[table_field]),
                    rel_tol=0,
                    abs_tol=1e-15,
                ):
                    raise ValueError(
                        f"Threshold plug-in-tail mismatch for {driver} K={threshold}"
                    )
            netbid_manifest = load_json(arm_root / "netbid2_manifest.json")
            if netbid_manifest.get("schema") != (
                "sjaracne-brca100-consensus-recurrence-netbid2-v1"
            ):
                raise ValueError(f"Unexpected NetBID2 schema: {arm_root}")
            if netbid_manifest.get("driver") != driver:
                raise ValueError(f"NetBID2 driver mismatch: {arm_root}")
            if netbid_manifest.get("minimum_support") != threshold:
                raise ValueError(f"NetBID2 threshold mismatch: {arm_root}")
            if netbid_manifest["input"]["network_manifest_sha256"] != sha256_file(
                arm_root / "network_manifest.json"
            ):
                raise ValueError(f"NetBID2/network-manifest linkage mismatch: {arm_root}")
            if netbid_manifest["input"]["network_sha256"] != network_manifest["ncol_sha256"]:
                raise ValueError(f"NetBID2/network linkage mismatch: {arm_root}")
            driver_filename = str(ARMS_FOR_PACKAGE[driver]["driver_file"])
            expected_driver_hash = design["inputs"][driver_filename]["sha256"]
            if netbid_manifest["input"]["driver_sha256"] != expected_driver_hash:
                raise ValueError(f"NetBID2/driver linkage mismatch: {arm_root}")
            if netbid_manifest["input"]["r_script_sha256"] != sha256_file(
                script_sources["run_netbid_qc.R"]
            ):
                raise ValueError(f"NetBID2 R-script linkage mismatch: {arm_root}")
            if netbid_manifest["input"]["wrapper_sha256"] != sha256_file(
                script_sources["netbid2-r"]
            ):
                raise ValueError(f"NetBID2 wrapper linkage mismatch: {arm_root}")
            if netbid_manifest["input"].get("environment") != EXPECTED_NETBID_ENVIRONMENT:
                raise ValueError(f"NetBID2 environment mismatch: {arm_root}")
            if netbid_manifest["input"].get("generate_html") is not False:
                raise ValueError(f"Unexpected NetBID2 HTML mode: {arm_root}")
            required_outputs = {
                "network_summary.tsv",
                "driver_target_sizes.tsv",
                "netbid_environment.tsv",
            }
            if set(netbid_manifest.get("outputs", {})) != required_outputs:
                raise ValueError(f"Incomplete NetBID2 output inventory: {arm_root}")
            expected_arm_files = {
                "network_manifest.json": None,
                "bootstrap_info_.txt": network_manifest["bootstrap_info_sha256"],
                "parameter_info_.txt": network_manifest["parameter_info_sha256"],
                "netbid2_manifest.json": None,
                "netbid2.stdout.log": netbid_manifest["stdout_sha256"],
                "netbid2.stderr.log": netbid_manifest["stderr_sha256"],
            }
            for filename, expected_hash in expected_arm_files.items():
                source = arm_root / filename
                if expected_hash is not None and sha256_file(source) != expected_hash:
                    raise ValueError(f"Arm manifest hash mismatch: {source}")
                copy_verified(
                    source,
                    partial_root,
                    Path("provenance/arms") / driver / f"k{threshold:03d}" / filename,
                    copied,
                )
            for filename, expected_hash in netbid_manifest["outputs"].items():
                source = arm_root / "netbid2_qc" / filename
                if sha256_file(source) != expected_hash:
                    raise ValueError(f"NetBID2 output hash mismatch: {source}")
                copy_verified(
                    source,
                    partial_root,
                    Path("netbid2") / driver / f"k{threshold:03d}" / filename,
                    copied,
                )
            if read_netbid_environment(
                arm_root / "netbid2_qc" / "netbid_environment.tsv"
            ) != EXPECTED_NETBID_ENVIRONMENT:
                raise ValueError(f"NetBID2 environment artifact mismatch: {arm_root}")
            for filename, field in (
                ("consensus_network_3col_.txt", "three_col_sha256"),
                ("consensus_network_ncol_.txt", "ncol_sha256"),
                ("consensus_support.tsv", "support_sha256"),
            ):
                record_omitted(arm_root / filename, work_root, network_manifest[field], omitted)
        source_root = Path(design["source_work_root"])
        source_point = str(ARMS_FOR_PACKAGE[driver]["point"])
        source_arm = source_root / "results" / source_point / driver
        source_consensus_manifest = load_json(source_arm / "consensus_manifest.json")
        source_support_manifest = load_json(source_arm / "support_summary_manifest.json")
        source_point_record = design["arms"][driver]["source_point"]
        if source_support_manifest["point_manifest_sha256"] != source_point_record[
            "sha256"
        ]:
            raise ValueError(f"K=9 point/support linkage mismatch for {driver}")
        if sha256_file(Path(source_point_record["path"])) != source_point_record["sha256"]:
            raise ValueError(f"K=9 point manifest bytes changed for {driver}")
        source_point_manifest = load_json(Path(source_point_record["path"]))
        if source_support_manifest.get("driver") != driver:
            raise ValueError(f"K=9 support driver mismatch for {driver}")
        if source_support_manifest.get("point") != source_point:
            raise ValueError(f"K=9 support point mismatch for {driver}")
        if float(source_support_manifest.get("p_value")) != float(
            design["arms"][driver]["per_subsample_p"]
        ):
            raise ValueError(f"K=9 support p-value mismatch for {driver}")
        if source_support_manifest["consensus_manifest_sha256"] != sha256_file(
            source_arm / "consensus_manifest.json"
        ):
            raise ValueError(f"K=9 support/consensus SHA linkage mismatch for {driver}")
        if source_support_manifest["consensus_fingerprint"] != source_consensus_manifest[
            "fingerprint"
        ]:
            raise ValueError(f"K=9 support/consensus linkage mismatch for {driver}")
        for filename in (
            "consensus_manifest.json",
            "support_summary_manifest.json",
        ):
            source = source_arm / filename
            copy_verified(
                source,
                partial_root,
                Path("provenance/k9_anchor") / driver / filename,
                copied,
            )
        copy_verified(
            Path(source_point_record["path"]),
            partial_root,
            Path("provenance/k9_anchor") / driver / "point_manifest.json",
            copied,
        )
        source_netbid_manifest_path = source_arm / "netbid2_qc_manifest.json"
        source_netbid_manifest = load_json(source_netbid_manifest_path)
        source_netbid_fingerprint_payload = {
            field: source_netbid_manifest[field]
            for field in SOURCE_NETBID_FINGERPRINT_FIELDS
        }
        if source_netbid_manifest.get("schema") != SOURCE_NETBID_SCHEMA:
            raise ValueError(f"Unexpected source K=9 NetBID2 schema for {driver}")
        if source_netbid_manifest.get("fingerprint") != fingerprint(
            source_netbid_fingerprint_payload
        ):
            raise ValueError(f"Source K=9 NetBID2 fingerprint mismatch for {driver}")
        expected_source_netbid = {
            "mode": "summary",
            "point": source_point,
            "p_value": float(design["arms"][driver]["per_subsample_p"]),
            "mi_cutoff": float(source_point_manifest["mi_cutoff"]),
            "point_manifest_sha256": source_point_record["sha256"],
            "sweep_design_sha256": design["source_sweep_design"]["sha256"],
            "driver": driver,
            "driver_sha256": design["inputs"][
                str(ARMS_FOR_PACKAGE[driver]["driver_file"])
            ]["sha256"],
            "consensus_sha256": source_consensus_manifest["ncol"]["sha256"],
            "consensus_manifest_sha256": sha256_file(
                source_arm / "consensus_manifest.json"
            ),
            "r_script_sha256": sha256_file(script_sources["run_netbid_qc.R"]),
            "wrapper_sha256": sha256_file(script_sources["netbid2-r"]),
            "environment": EXPECTED_NETBID_ENVIRONMENT,
            "prefix": "TF_" if driver == "tf" else "SIG_",
        }
        for field, expected_value in expected_source_netbid.items():
            observed_value = source_netbid_manifest.get(field)
            if field in {"p_value", "mi_cutoff"}:
                if not math.isclose(
                    float(observed_value),
                    float(expected_value),
                    rel_tol=0,
                    abs_tol=1e-15,
                ):
                    raise ValueError(
                        f"Source K=9 NetBID2 {field} mismatch for {driver}"
                    )
            elif observed_value != expected_value:
                raise ValueError(
                    f"Source K=9 NetBID2 {field} mismatch for {driver}"
                )
        source_netbid_inventory = source_netbid_manifest.get("output_inventory")
        if not isinstance(source_netbid_inventory, list):
            raise ValueError(f"Malformed source K=9 NetBID2 inventory for {driver}")
        inventory_by_name = {
            entry["path"]: entry
            for entry in source_netbid_inventory
            if isinstance(entry, dict) and set(entry) == {"path", "bytes", "sha256"}
        }
        if set(inventory_by_name) != SOURCE_NETBID_OUTPUTS or len(
            inventory_by_name
        ) != len(source_netbid_inventory):
            raise ValueError(f"Incomplete source K=9 NetBID2 inventory for {driver}")
        current_k9_netbid = work_root / "results" / driver / "k009" / "netbid2_qc"
        for filename in sorted(SOURCE_NETBID_OUTPUTS):
            if Path(filename).name != filename:
                raise ValueError(f"Unsafe source K=9 NetBID2 output path: {filename}")
            entry = inventory_by_name[filename]
            source_output = source_arm / "netbid2_qc" / filename
            current_output = current_k9_netbid / filename
            if (
                source_output.stat().st_size != entry["bytes"]
                or current_output.stat().st_size != entry["bytes"]
                or sha256_file(source_output) != entry["sha256"]
                or sha256_file(current_output) != entry["sha256"]
            ):
                raise ValueError(
                    f"K=9 NetBID2 output was not reproduced for {driver}: {filename}"
                )
        copy_verified(
            source_netbid_manifest_path,
            partial_root,
            Path("provenance/k9_anchor") / driver / "netbid2_qc_manifest.json",
            copied,
        )
        source_consensus = source_arm / "consensus"
        if sha256_file(source_consensus / "bootstrap_info_.txt") != source_consensus_manifest[
            "bootstrap_info_sha256"
        ]:
            raise ValueError(f"K=9 bootstrap evidence mismatch for {driver}")
        copy_verified(
            source_consensus / "bootstrap_info_.txt",
            partial_root,
            Path("provenance/k9_anchor") / driver / "bootstrap_info_.txt",
            copied,
        )
        record_omitted(
            source_consensus / "consensus_support.tsv",
            work_root,
            source_support_manifest["output_sha256"],
            omitted,
            recorded_path=(
                f"source_sweep/results/{source_point}/{driver}/consensus/"
                "consensus_support.tsv"
            ),
        )
        k9_manifest = load_json(work_root / "results" / driver / "k009/network_manifest.json")
        if source_support_manifest["retained_edges"] != k9_manifest["edges"]:
            raise ValueError(f"K=9 retained-edge mismatch for {driver}")
        if source_support_manifest["consensus_sha256"] != source_consensus_manifest["ncol"][
            "sha256"
        ]:
            raise ValueError(f"K=9 source NCOL linkage mismatch for {driver}")
        if k9_manifest["ncol_sha256"] != source_consensus_manifest["ncol"]["sha256"]:
            raise ValueError(f"K=9 NCOL was not reproduced byte-for-byte for {driver}")
        if k9_manifest["three_col_sha256"] != source_consensus_manifest[
            "consensus_3col_sha256"
        ]:
            raise ValueError(f"K=9 3-column network was not reproduced for {driver}")
    omitted.sort(key=lambda row: str(row["path"]))
    omitted_payload = {
        "schema": "sjaracne-brca100-consensus-recurrence-omitted-v1",
        "count": len(omitted),
        "files": omitted,
    }
    (partial_root / "omitted_artifacts.json").write_bytes(serialized_json(omitted_payload))
    if results_markdown is not None:
        copy_verified(results_markdown, partial_root, Path("RESULTS.md"), copied)
    package_core = {
        "schema": PACKAGE_SCHEMA,
        "design_fingerprint": design["fingerprint"],
        "analysis_fingerprint": analysis_manifest["fingerprint"],
        "live_work_root": str(work_root),
        "packager_sha256": sha256_file(Path(__file__).resolve()),
        "checksum_policy": "SHA-256 over exact bytes; package-local files are Git -text",
        "copied_file_count": len(copied),
        "omitted_file_count": len(omitted),
        "copied_files": sorted(copied, key=lambda row: str(row["path"])),
        "omitted_manifest_sha256": sha256_file(partial_root / "omitted_artifacts.json"),
    }
    package_payload = dict(package_core)
    package_payload["fingerprint"] = hashlib.sha256(
        json.dumps(package_core, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    (partial_root / "package_manifest.json").write_bytes(serialized_json(package_payload))
    write_sha256s(partial_root)
    os.replace(partial_root, output_root)
    owner_path.unlink()
    print(f"Packaged recurrence sweep: {output_root}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--results-markdown", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    package(args.work_root.resolve(), args.output_root.resolve(), args.results_markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
