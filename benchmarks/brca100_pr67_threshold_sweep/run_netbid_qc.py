#!/usr/bin/env python3
"""Run resumable NetBID2 QC for dynamically discovered PR67 sweep points."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
from pathlib import PurePosixPath
import shutil
import subprocess
from typing import Any


DRIVERS = {
    "tf": ("BRCA100_TF.txt", "TF_", 2608),
    "sig": ("BRCA100_SIG.txt", "SIG_", 10680),
}
POINT_SCHEMA = "sjaracne-brca100-pr67-p-sweep-point-v1"
LEGACY_SWEEP_DESIGN_SCHEMA = "sjaracne-brca100-pr67-p-sweep-v1"
SWEEP_DESIGN_SCHEMA = "sjaracne-brca100-pr67-p-sweep-v2"
SWEEP_DESIGN_MIGRATION_SCHEMA = (
    "sjaracne-brca100-pr67-p-sweep-design-migration-v1"
)
NETBID_MANIFEST_MIGRATION_SCHEMA = (
    "sjaracne-brca100-pr67-p-sweep-netbid2-manifest-migration-v1"
)
SWEEP_DESIGN_HISTORY_DIRECTORY = "sweep_design_history"
NETBID_MANIFEST_HISTORY_DIRECTORY = "netbid2_manifest_history"
RUN_SCHEMA = "sjaracne-brca100-pr67-p-sweep-netbid2-v1"
RUN_FINGERPRINT_FIELDS = frozenset(
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
SUMMARY_AGGREGATE_SCHEMA = (
    "sjaracne-brca100-pr67-p-sweep-netbid2-summary-aggregate-v1"
)
HTML_AGGREGATE_SCHEMA = (
    "sjaracne-brca100-pr67-p-sweep-netbid2-html-aggregate-v1"
)
REQUIRED_METRICS = {
    "candidate_drivers",
    "active_drivers",
    "active_driver_fraction",
    "edges",
    "incident_nodes",
    "weak_components",
    "largest_weak_component",
    "largest_weak_component_fraction",
    "density",
    "target_size_zero_mean",
    "target_size_zero_median",
    "target_size_zero_q25",
    "target_size_zero_q75",
    "target_size_zero_max",
    "target_size_active_mean",
    "target_size_active_median",
    "target_size_active_q25",
    "target_size_active_q75",
    "target_size_active_max",
    "scale_free_adjusted_r2",
}


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


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


def serialized_json(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def atomic_bytes(path: Path, value: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("wb") as handle:
        handle.write(value)
    os.replace(partial, path)


def atomic_json(path: Path, value: object) -> None:
    atomic_bytes(path, serialized_json(value))


def ensure_exact_bytes(path: Path, value: bytes, description: str) -> None:
    if path.exists():
        if not path.is_file() or path.read_bytes() != value:
            raise RuntimeError(f"Incompatible existing {description}: {path}")
        return
    atomic_bytes(path, value)


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def discover_points(
    work_root: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any], str]:
    results_root = work_root / "results"
    sweep_design_path = work_root / "sweep_design.json"
    if not sweep_design_path.is_file():
        raise FileNotFoundError(sweep_design_path)
    sweep_design = load_json(sweep_design_path)
    if sweep_design.get("schema") not in {
        LEGACY_SWEEP_DESIGN_SCHEMA,
        SWEEP_DESIGN_SCHEMA,
    }:
        raise ValueError(f"Unexpected sweep design schema: {sweep_design_path}")
    design_points = sweep_design.get("all_points")
    if not isinstance(design_points, list) or not design_points:
        raise ValueError(f"Sweep design has no all_points list: {sweep_design_path}")
    if any(not isinstance(item, dict) or not isinstance(item.get("key"), str)
           for item in design_points):
        raise ValueError(f"Malformed all_points in {sweep_design_path}")
    design_by_key = {item["key"]: item for item in design_points}
    if len(design_by_key) != len(design_points):
        raise ValueError(f"Duplicate point key in {sweep_design_path}")

    manifest_paths = sorted(results_root.glob("*/point_manifest.json"))
    manifest_keys = {path.parent.name for path in manifest_paths}
    if manifest_keys != set(design_by_key):
        raise ValueError(
            "Point manifests do not exactly cover sweep_design.json all_points: "
            f"missing={sorted(set(design_by_key) - manifest_keys)}, "
            f"extra={sorted(manifest_keys - set(design_by_key))}"
        )
    points: list[dict[str, Any]] = []
    seen_probabilities: set[float] = set()
    for manifest_path in manifest_paths:
        manifest = load_json(manifest_path)
        key = manifest.get("key")
        if manifest.get("schema") != POINT_SCHEMA:
            raise ValueError(f"Unexpected point schema: {manifest_path}")
        if not isinstance(key, str) or key != manifest_path.parent.name:
            raise ValueError(f"Point key/path mismatch: {manifest_path}")
        design_point = design_by_key[key]
        mismatched_fields = [
            field for field, value in design_point.items()
            if manifest.get(field) != value
        ]
        if mismatched_fields:
            raise ValueError(
                f"Point manifest disagrees with sweep design at {manifest_path}: "
                + ", ".join(mismatched_fields)
            )
        probability = float(manifest["p_value"])
        cutoff = float(manifest["mi_cutoff"])
        if not math.isfinite(probability) or probability <= 0.0:
            raise ValueError(f"Invalid point probability: {manifest_path}")
        if not math.isfinite(cutoff) or cutoff <= 0.0:
            raise ValueError(f"Invalid point cutoff: {manifest_path}")
        if probability in seen_probabilities:
            raise ValueError(f"Duplicate point probability {probability:.17g}")
        seen_probabilities.add(probability)
        points.append(
            {
                "key": key,
                "p_value": probability,
                "mi_cutoff": cutoff,
                "manifest": manifest,
                "manifest_path": manifest_path,
                "root": manifest_path.parent,
            }
        )
    if not points:
        raise ValueError(f"No immutable point manifests found under {results_root}")
    return (
        sorted(points, key=lambda item: (item["p_value"], item["key"])),
        sweep_design,
        sha256_file(sweep_design_path),
    )


def require_exact_relative_path(
    value: object, expected: str, field: str, source: Path
) -> None:
    if not isinstance(value, str) or not value or "\\" in value:
        raise ValueError(f"Invalid {field} in {source}: {value!r}")
    candidate = PurePosixPath(value)
    if candidate.is_absolute() or any(part in {"", ".", ".."} for part in candidate.parts):
        raise ValueError(f"Unsafe {field} in {source}: {value!r}")
    if candidate.as_posix() != value or value != expected:
        raise ValueError(
            f"Unexpected {field} in {source}: {value!r} != {expected!r}"
        )


def load_sweep_design_migration(
    *,
    work_root: Path,
    sweep_design: dict[str, Any],
    sweep_design_hash: str,
) -> dict[str, Any] | None:
    """Validate the exact append-only design history that permits reuse.

    A freshly created v2 root has no history and returns ``None``.  A migrated
    root is accepted only when its archived v1 bytes and deterministic
    migration manifest reproduce the active v2 design exactly.
    """
    schema = sweep_design.get("schema")
    if schema == LEGACY_SWEEP_DESIGN_SCHEMA:
        return None
    if schema != SWEEP_DESIGN_SCHEMA:
        raise ValueError(f"Unexpected active sweep schema: {schema!r}")

    history_root = work_root / SWEEP_DESIGN_HISTORY_DIRECTORY
    if not history_root.exists():
        return None
    if not history_root.is_dir():
        raise RuntimeError(f"Sweep-design history is not a directory: {history_root}")
    candidates = sorted(
        history_root.glob(f"*_to_{sweep_design_hash}.migration.json")
    )
    if not candidates:
        return None
    if len(candidates) != 1:
        raise RuntimeError(
            "Expected exactly one sweep-design migration into the active design: "
            + ", ".join(str(path) for path in candidates)
        )
    migration_path = candidates[0]
    migration = load_json(migration_path)
    expected_fields = {
        "schema",
        "operation",
        "from",
        "to",
        "appended_points",
        "manifest_path",
        "fingerprint",
    }
    if set(migration) != expected_fields:
        raise ValueError(
            f"Unexpected sweep-design migration fields in {migration_path}"
        )
    payload = dict(migration)
    observed_fingerprint = payload.pop("fingerprint")
    if observed_fingerprint != fingerprint(payload):
        raise ValueError(
            f"Sweep-design migration fingerprint mismatch: {migration_path}"
        )
    if (
        migration.get("schema") != SWEEP_DESIGN_MIGRATION_SCHEMA
        or migration.get("operation") != "append-only-point-extension"
    ):
        raise ValueError(f"Unexpected sweep-design migration: {migration_path}")

    from_record = migration.get("from")
    to_record = migration.get("to")
    if (
        not isinstance(from_record, dict)
        or set(from_record) != {"schema", "sha256", "archived_path", "point_keys"}
        or not isinstance(to_record, dict)
        or set(to_record) != {"schema", "sha256", "active_path", "point_keys"}
    ):
        raise ValueError(f"Malformed design endpoints in {migration_path}")
    prior_hash = from_record.get("sha256")
    if (
        not isinstance(prior_hash, str)
        or len(prior_hash) != 64
        or any(character not in "0123456789abcdef" for character in prior_hash)
        or from_record.get("schema") != LEGACY_SWEEP_DESIGN_SCHEMA
        or to_record.get("schema") != SWEEP_DESIGN_SCHEMA
        or to_record.get("sha256") != sweep_design_hash
    ):
        raise ValueError(f"Invalid design endpoints in {migration_path}")

    archive_relative = (
        f"{SWEEP_DESIGN_HISTORY_DIRECTORY}/{prior_hash}.sweep_design.json"
    )
    migration_relative = (
        f"{SWEEP_DESIGN_HISTORY_DIRECTORY}/"
        f"{prior_hash}_to_{sweep_design_hash}.migration.json"
    )
    require_exact_relative_path(
        from_record.get("archived_path"),
        archive_relative,
        "from.archived_path",
        migration_path,
    )
    require_exact_relative_path(
        to_record.get("active_path"),
        "sweep_design.json",
        "to.active_path",
        migration_path,
    )
    require_exact_relative_path(
        migration.get("manifest_path"),
        migration_relative,
        "manifest_path",
        migration_path,
    )
    if migration_path != work_root / migration_relative:
        raise ValueError(f"Sweep-design migration path mismatch: {migration_path}")

    archive_path = work_root / archive_relative
    if not archive_path.is_file() or sha256_file(archive_path) != prior_hash:
        raise RuntimeError(f"Archived sweep-design bytes changed: {archive_path}")
    archived_design = load_json(archive_path)
    old_points = archived_design.get("all_points")
    new_points = sweep_design.get("all_points")
    old_keys = from_record.get("point_keys")
    new_keys = to_record.get("point_keys")
    if (
        archived_design.get("schema") != LEGACY_SWEEP_DESIGN_SCHEMA
        or not isinstance(old_points, list)
        or not isinstance(new_points, list)
        or not isinstance(old_keys, list)
        or not isinstance(new_keys, list)
        or any(not isinstance(key, str) for key in old_keys + new_keys)
        or old_keys != [point.get("key") for point in old_points]
        or new_keys != [point.get("key") for point in new_points]
        or len(set(old_keys)) != len(old_keys)
        or len(set(new_keys)) != len(new_keys)
        or new_points[: len(old_points)] != old_points
        or migration.get("appended_points") != new_points[len(old_points) :]
    ):
        raise ValueError(f"Design migration is not an exact point append: {migration_path}")
    old_invariants = {
        key: value
        for key, value in archived_design.items()
        if key not in {"schema", "all_points"}
    }
    new_invariants = {
        key: value
        for key, value in sweep_design.items()
        if key not in {"schema", "all_points"}
    }
    if old_invariants != new_invariants:
        raise ValueError(
            f"Design migration changes a fixed sweep invariant: {migration_path}"
        )
    return {
        "prior_sweep_design_sha256": prior_hash,
        "current_sweep_design_sha256": sweep_design_hash,
        "prior_point_keys": tuple(old_keys),
        "current_point_keys": tuple(new_keys),
        "sweep_migration_path": migration_path,
        "sweep_migration_relative": migration_relative,
        "sweep_migration_sha256": sha256_file(migration_path),
        "sweep_migration_fingerprint": observed_fingerprint,
    }


def select_values(
    specification: str, available: list[str], label: str, *, allow_none: bool = False
) -> list[str]:
    if specification == "all":
        return list(available)
    if allow_none and specification in ("", "none"):
        return []
    requested = [item.strip() for item in specification.split(",") if item.strip()]
    unknown = sorted(set(requested) - set(available))
    if unknown:
        raise ValueError(f"Unknown {label}(s): {', '.join(unknown)}")
    if not requested:
        raise ValueError(f"At least one {label} must be selected")
    if len(set(requested)) != len(requested):
        raise ValueError(f"Duplicate {label} selected")
    return requested


def probe_environment(wrapper: Path, r_script: Path) -> tuple[dict[str, str], str]:
    completed = subprocess.run(
        [str(wrapper), "Rscript", str(r_script), "--probe"],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "NetBID2 environment probe failed:\n" + completed.stderr.strip()
        )
    rows = list(csv.DictReader(completed.stdout.splitlines(), delimiter="\t"))
    if not rows or any(set(row) != {"component", "version"} for row in rows):
        raise ValueError("Unexpected NetBID2 environment probe output")
    environment = {row["component"]: row["version"] for row in rows}
    if (
        len(rows) != 4
        or len(environment) != len(rows)
        or set(environment) != {"R", "NetBID2", "NetBID2_remote_sha", "igraph"}
    ):
        raise ValueError(f"Incomplete NetBID2 environment probe: {environment}")
    return environment, completed.stderr


def inventory(root: Path) -> list[dict[str, object]]:
    if not root.is_dir():
        raise RuntimeError(f"Missing output directory: {root}")
    records = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        if path.stat().st_size == 0:
            raise RuntimeError(f"Empty NetBID2 artifact: {path}")
        records.append(
            {
                "path": path.relative_to(root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    if not records:
        raise RuntimeError(f"Empty NetBID2 output directory: {root}")
    return records


def read_network_summary(path: Path) -> dict[str, float]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or any(set(row) != {"metric", "value"} for row in rows):
        raise ValueError(f"Unexpected network summary format: {path}")
    if any(not row["metric"] for row in rows):
        raise ValueError(f"Empty metric in {path}")
    if len({row["metric"] for row in rows}) != len(rows):
        raise ValueError(f"Duplicate metric in {path}")
    metrics = {
        row["metric"]: float("nan") if row["value"] == "NA" else float(row["value"])
        for row in rows
    }
    if set(metrics) != REQUIRED_METRICS:
        raise ValueError(
            f"Unexpected metrics in {path}: missing={sorted(REQUIRED_METRICS - set(metrics))}, "
            f"extra={sorted(set(metrics) - REQUIRED_METRICS)}"
        )
    return metrics


def validate_output(
    root: Path,
    *,
    mode: str,
    prefix: str,
    driver_ids: list[str],
    expected_edges: int,
    expected_environment: dict[str, str],
) -> list[dict[str, object]]:
    required = [
        root / "network_summary.tsv",
        root / "driver_target_sizes.tsv",
        root / "netbid_environment.tsv",
    ]
    if mode == "html":
        required.append(root / f"{prefix}netQC.html")
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing/empty NetBID2 output: {path}")

    metrics = read_network_summary(root / "network_summary.tsv")
    if int(metrics["candidate_drivers"]) != len(driver_ids):
        raise ValueError(f"Candidate-driver count mismatch in {root}")
    if int(metrics["edges"]) != expected_edges:
        raise ValueError(f"Consensus/NetBID2 edge-count mismatch in {root}")

    with (root / "driver_target_sizes.tsv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if any(set(row) != {"driver", "target_count"} for row in rows):
        raise ValueError(f"Unexpected driver-target format in {root}")
    if [row["driver"] for row in rows] != driver_ids:
        raise ValueError(f"Driver order/content mismatch in {root}")
    target_counts = [int(row["target_count"]) for row in rows]
    if any(value < 0 for value in target_counts):
        raise ValueError(f"Negative driver target count in {root}")
    if sum(target_counts) != expected_edges:
        raise ValueError(f"Driver target counts do not sum to edge count in {root}")
    if sum(value > 0 for value in target_counts) != int(metrics["active_drivers"]):
        raise ValueError(f"Active-driver count mismatch in {root}")

    with (root / "netbid_environment.tsv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        environment_rows = list(csv.DictReader(handle, delimiter="\t"))
    output_environment = {
        row.get("component"): row.get("version") for row in environment_rows
    }
    if (
        len(environment_rows) != len(expected_environment)
        or len(output_environment) != len(environment_rows)
        or any(set(row) != {"component", "version"} for row in environment_rows)
        or output_environment != expected_environment
    ):
        raise ValueError(f"Unexpected environment table in {root}")
    return inventory(root)


def validate_record(
    *,
    root: Path,
    record: dict[str, Any],
    expected_fingerprint: str,
    mode: str,
    prefix: str,
    driver_ids: list[str],
    expected_edges: int,
    stdout_path: Path,
    stderr_path: Path,
    expected_environment: dict[str, str],
) -> None:
    actual_inventory = validate_output(
        root,
        mode=mode,
        prefix=prefix,
        driver_ids=driver_ids,
        expected_edges=expected_edges,
        expected_environment=expected_environment,
    )
    for log_path in (stdout_path, stderr_path):
        if not log_path.is_file():
            raise RuntimeError(f"Missing NetBID2 log: {log_path}")
    expected = {
        "fingerprint": expected_fingerprint,
        "output_inventory": actual_inventory,
        "stdout_sha256": sha256_file(stdout_path),
        "stderr_sha256": sha256_file(stderr_path),
        "stderr_bytes": stderr_path.stat().st_size,
    }
    mismatches = [key for key, value in expected.items() if record.get(key) != value]
    if mismatches:
        raise RuntimeError(
            f"NetBID2 manifest validation failed at {root}: {', '.join(mismatches)}"
        )


def validate_completed_run_record(
    *,
    source: Path,
    root: Path,
    record: dict[str, Any],
    expected_payload: dict[str, Any],
    expected_command: list[str],
    mode: str,
    prefix: str,
    driver_ids: list[str],
    expected_edges: int,
    stdout_path: Path,
    stderr_path: Path,
    expected_environment: dict[str, str],
) -> None:
    """Validate the complete immutable run contract, not just its checksum."""
    if set(expected_payload) != RUN_FINGERPRINT_FIELDS:
        raise RuntimeError("Internal NetBID2 fingerprint payload field mismatch")
    required_fields = RUN_FINGERPRINT_FIELDS | {
        "fingerprint",
        "command",
        "finished_at_utc",
        "output",
        "output_inventory",
        "stdout_sha256",
        "stderr_sha256",
        "stderr_bytes",
    }
    allowed_fields = required_fields | {"recovered_after_interrupted_manifest"}
    record_fields = set(record)
    if record_fields != required_fields and record_fields != allowed_fields:
        raise RuntimeError(
            f"NetBID2 manifest has missing/arbitrary fields at {source}: "
            f"missing={sorted(required_fields - record_fields)}, "
            f"extra={sorted(record_fields - allowed_fields)}"
        )
    if (
        "recovered_after_interrupted_manifest" in record
        and record["recovered_after_interrupted_manifest"] is not True
    ):
        raise RuntimeError(f"Invalid NetBID2 recovery marker in {source}")
    mismatches = [
        field
        for field, expected in expected_payload.items()
        if record.get(field) != expected
    ]
    if mismatches:
        raise RuntimeError(
            f"NetBID2 fingerprint inputs changed in {source}: "
            + ", ".join(mismatches)
        )
    expected_fingerprint = fingerprint(expected_payload)
    if record.get("fingerprint") != expected_fingerprint:
        raise RuntimeError(f"NetBID2 fingerprint mismatch in {source}")
    if record.get("command") != expected_command:
        raise RuntimeError(f"NetBID2 command mismatch in {source}")
    if (
        not isinstance(record.get("finished_at_utc"), str)
        or not record["finished_at_utc"]
    ):
        raise RuntimeError(f"Missing NetBID2 completion timestamp in {source}")
    output = record.get("output")
    if not isinstance(output, str) or Path(output).resolve() != root.resolve():
        raise RuntimeError(f"NetBID2 output path mismatch in {source}")
    validate_record(
        root=root,
        record=record,
        expected_fingerprint=expected_fingerprint,
        mode=mode,
        prefix=prefix,
        driver_ids=driver_ids,
        expected_edges=expected_edges,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        expected_environment=expected_environment,
    )


def manifest_fingerprint_payload(record: dict[str, Any]) -> dict[str, Any]:
    if not RUN_FINGERPRINT_FIELDS.issubset(record):
        raise ValueError("NetBID2 manifest lacks fingerprint payload fields")
    return {field: record[field] for field in RUN_FINGERPRINT_FIELDS}


def validate_netbid_migration_audit(
    *,
    work_root: Path,
    migration: dict[str, Any],
    audit_path: Path,
    audit: dict[str, Any],
) -> list[dict[str, Any]]:
    expected_fields = {
        "schema",
        "operation",
        "sweep_design_migration",
        "from_sweep_design_sha256",
        "to_sweep_design_sha256",
        "migrated_runs",
        "fingerprint",
    }
    if set(audit) != expected_fields:
        raise RuntimeError(f"Unexpected NetBID2 migration audit fields: {audit_path}")
    payload = dict(audit)
    observed_fingerprint = payload.pop("fingerprint")
    if observed_fingerprint != fingerprint(payload):
        raise RuntimeError(f"NetBID2 migration audit fingerprint mismatch: {audit_path}")
    expected_sweep_migration = {
        "manifest_path": migration["sweep_migration_relative"],
        "manifest_sha256": migration["sweep_migration_sha256"],
        "fingerprint": migration["sweep_migration_fingerprint"],
    }
    if (
        audit.get("schema") != NETBID_MANIFEST_MIGRATION_SCHEMA
        or audit.get("operation")
        != "refresh-sweep-design-fingerprint-after-append-only-extension"
        or audit.get("sweep_design_migration") != expected_sweep_migration
        or audit.get("from_sweep_design_sha256")
        != migration["prior_sweep_design_sha256"]
        or audit.get("to_sweep_design_sha256")
        != migration["current_sweep_design_sha256"]
    ):
        raise RuntimeError(f"NetBID2 migration audit provenance mismatch: {audit_path}")
    runs = audit.get("migrated_runs")
    if not isinstance(runs, list) or any(not isinstance(item, dict) for item in runs):
        raise RuntimeError(f"Malformed migrated_runs in {audit_path}")
    expected_entry_fields = {
        "point",
        "driver",
        "mode",
        "active_manifest_path",
        "archived_manifest_path",
        "archived_manifest_sha256",
        "current_manifest_sha256",
        "old_fingerprint",
        "new_fingerprint",
    }
    identities: list[tuple[str, str, str]] = []
    pair_name = (
        f"{migration['prior_sweep_design_sha256']}_to_"
        f"{migration['current_sweep_design_sha256']}"
    )
    for entry in runs:
        if set(entry) != expected_entry_fields:
            raise RuntimeError(f"Malformed run entry in {audit_path}")
        point = entry.get("point")
        driver = entry.get("driver")
        mode = entry.get("mode")
        if (
            not isinstance(point, str)
            or point not in migration["prior_point_keys"]
            or driver not in DRIVERS
            or mode not in {"summary", "html"}
        ):
            raise RuntimeError(f"Invalid migrated run identity in {audit_path}")
        suffix = "netbid2_qc" if mode == "summary" else "netbid2_qc_html"
        active_relative = f"results/{point}/{driver}/{suffix}_manifest.json"
        archived_relative = (
            f"{NETBID_MANIFEST_HISTORY_DIRECTORY}/{pair_name}/arms/"
            f"{point}/{driver}/{suffix}_manifest.json"
        )
        require_exact_relative_path(
            entry.get("active_manifest_path"),
            active_relative,
            "active_manifest_path",
            audit_path,
        )
        require_exact_relative_path(
            entry.get("archived_manifest_path"),
            archived_relative,
            "archived_manifest_path",
            audit_path,
        )
        active_path = work_root / active_relative
        archived_path = work_root / archived_relative
        for field in (
            "archived_manifest_sha256",
            "current_manifest_sha256",
            "old_fingerprint",
            "new_fingerprint",
        ):
            value = entry.get(field)
            if (
                not isinstance(value, str)
                or len(value) != 64
                or any(character not in "0123456789abcdef" for character in value)
            ):
                raise RuntimeError(f"Invalid {field} in {audit_path}")
        if (
            not active_path.is_file()
            or not archived_path.is_file()
            or sha256_file(active_path) != entry["current_manifest_sha256"]
            or sha256_file(archived_path) != entry["archived_manifest_sha256"]
        ):
            raise RuntimeError(f"Migrated NetBID2 manifest bytes changed: {point}/{driver}/{mode}")
        archived_record = load_json(archived_path)
        current_record = load_json(active_path)
        if (
            archived_record.get("sweep_design_sha256")
            != migration["prior_sweep_design_sha256"]
            or current_record.get("sweep_design_sha256")
            != migration["current_sweep_design_sha256"]
            or archived_record.get("fingerprint") != entry["old_fingerprint"]
            or current_record.get("fingerprint") != entry["new_fingerprint"]
            or fingerprint(manifest_fingerprint_payload(archived_record))
            != entry["old_fingerprint"]
            or fingerprint(manifest_fingerprint_payload(current_record))
            != entry["new_fingerprint"]
        ):
            raise RuntimeError(f"Migrated NetBID2 fingerprint changed: {point}/{driver}/{mode}")
        projected = dict(archived_record)
        projected["sweep_design_sha256"] = migration["current_sweep_design_sha256"]
        projected_payload = manifest_fingerprint_payload(projected)
        projected["fingerprint"] = fingerprint(projected_payload)
        if projected != current_record:
            raise RuntimeError(
                "Migrated NetBID2 manifests differ by more than the sweep-design "
                f"hash/fingerprint: {point}/{driver}/{mode}"
            )
        identities.append((point, driver, mode))
    if identities != sorted(identities) or len(set(identities)) != len(identities):
        raise RuntimeError(f"NetBID2 migration audit order/uniqueness mismatch: {audit_path}")
    return runs


def record_netbid_manifest_migration(
    *,
    work_root: Path,
    migration: dict[str, Any],
    entry: dict[str, Any],
) -> None:
    pair_name = (
        f"{migration['prior_sweep_design_sha256']}_to_"
        f"{migration['current_sweep_design_sha256']}"
    )
    audit_path = (
        work_root / NETBID_MANIFEST_HISTORY_DIRECTORY / pair_name / "migration.json"
    )
    if audit_path.is_file():
        existing = load_json(audit_path)
        runs = validate_netbid_migration_audit(
            work_root=work_root,
            migration=migration,
            audit_path=audit_path,
            audit=existing,
        )
    elif audit_path.exists():
        raise RuntimeError(f"NetBID2 migration audit is not a file: {audit_path}")
    else:
        runs = []
    identity = (entry["point"], entry["driver"], entry["mode"])
    existing_for_identity = [
        item
        for item in runs
        if (item["point"], item["driver"], item["mode"]) == identity
    ]
    if existing_for_identity and existing_for_identity != [entry]:
        raise RuntimeError(f"NetBID2 migration audit entry changed: {identity}")
    if not existing_for_identity:
        runs = [*runs, entry]
    runs = sorted(runs, key=lambda item: (item["point"], item["driver"], item["mode"]))
    payload: dict[str, Any] = {
        "schema": NETBID_MANIFEST_MIGRATION_SCHEMA,
        "operation": "refresh-sweep-design-fingerprint-after-append-only-extension",
        "sweep_design_migration": {
            "manifest_path": migration["sweep_migration_relative"],
            "manifest_sha256": migration["sweep_migration_sha256"],
            "fingerprint": migration["sweep_migration_fingerprint"],
        },
        "from_sweep_design_sha256": migration["prior_sweep_design_sha256"],
        "to_sweep_design_sha256": migration["current_sweep_design_sha256"],
        "migrated_runs": runs,
    }
    payload["fingerprint"] = fingerprint(payload)
    if not audit_path.is_file() or load_json(audit_path) != payload:
        atomic_json(audit_path, payload)
    validate_netbid_migration_audit(
        work_root=work_root,
        migration=migration,
        audit_path=audit_path,
        audit=load_json(audit_path),
    )


def migrate_completed_manifest_if_eligible(
    *,
    work_root: Path,
    migration: dict[str, Any] | None,
    point: str,
    driver: str,
    mode: str,
    manifest_path: Path,
    output_root: Path,
    record: dict[str, Any],
    expected_payload: dict[str, Any],
    expected_command: list[str],
    prefix: str,
    driver_ids: list[str],
    expected_edges: int,
    stdout_path: Path,
    stderr_path: Path,
    expected_environment: dict[str, str],
) -> dict[str, Any]:
    """Refresh only the design hash of a rigorously validated legacy record."""
    current_fingerprint = fingerprint(expected_payload)
    if migration is None or point not in migration["prior_point_keys"]:
        validate_completed_run_record(
            source=manifest_path,
            root=output_root,
            record=record,
            expected_payload=expected_payload,
            expected_command=expected_command,
            mode=mode,
            prefix=prefix,
            driver_ids=driver_ids,
            expected_edges=expected_edges,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            expected_environment=expected_environment,
        )
        return record
    if expected_payload["sweep_design_sha256"] != migration[
        "current_sweep_design_sha256"
    ]:
        raise RuntimeError("Active NetBID2 payload does not use the migrated design")

    prior_payload = dict(expected_payload)
    prior_payload["sweep_design_sha256"] = migration["prior_sweep_design_sha256"]
    prior_fingerprint = fingerprint(prior_payload)
    mode_suffix = "netbid2_qc" if mode == "summary" else "netbid2_qc_html"
    pair_name = (
        f"{migration['prior_sweep_design_sha256']}_to_"
        f"{migration['current_sweep_design_sha256']}"
    )
    archive_path = (
        work_root
        / NETBID_MANIFEST_HISTORY_DIRECTORY
        / pair_name
        / "arms"
        / point
        / driver
        / f"{mode_suffix}_manifest.json"
    )
    manifest_bytes = manifest_path.read_bytes()
    migrated = False
    if record.get("fingerprint") == prior_fingerprint:
        validate_completed_run_record(
            source=manifest_path,
            root=output_root,
            record=record,
            expected_payload=prior_payload,
            expected_command=expected_command,
            mode=mode,
            prefix=prefix,
            driver_ids=driver_ids,
            expected_edges=expected_edges,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            expected_environment=expected_environment,
        )
        ensure_exact_bytes(
            archive_path, manifest_bytes, "archived pre-extension NetBID2 manifest"
        )
        record = dict(record)
        record["sweep_design_sha256"] = migration["current_sweep_design_sha256"]
        record["fingerprint"] = current_fingerprint
        atomic_json(manifest_path, record)
        record = load_json(manifest_path)
        migrated = True
    elif record.get("fingerprint") != current_fingerprint:
        raise RuntimeError(
            "Completed NetBID2 manifest matches neither the active design nor the "
            f"single archived predecessor: {manifest_path}"
        )

    validate_completed_run_record(
        source=manifest_path,
        root=output_root,
        record=record,
        expected_payload=expected_payload,
        expected_command=expected_command,
        mode=mode,
        prefix=prefix,
        driver_ids=driver_ids,
        expected_edges=expected_edges,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        expected_environment=expected_environment,
    )
    if not archive_path.exists():
        if migrated:
            raise RuntimeError(f"NetBID2 manifest archive disappeared: {archive_path}")
        return record
    if not archive_path.is_file():
        raise RuntimeError(f"NetBID2 manifest archive is not a file: {archive_path}")
    archived_record = load_json(archive_path)
    validate_completed_run_record(
        source=archive_path,
        root=output_root,
        record=archived_record,
        expected_payload=prior_payload,
        expected_command=expected_command,
        mode=mode,
        prefix=prefix,
        driver_ids=driver_ids,
        expected_edges=expected_edges,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        expected_environment=expected_environment,
    )
    projected = dict(archived_record)
    projected["sweep_design_sha256"] = migration["current_sweep_design_sha256"]
    projected["fingerprint"] = current_fingerprint
    if projected != record:
        raise RuntimeError(
            "Archived/current NetBID2 manifests differ by more than the exact "
            f"design-hash migration: {manifest_path}"
        )
    active_relative = manifest_path.relative_to(work_root).as_posix()
    archived_relative = archive_path.relative_to(work_root).as_posix()
    entry = {
        "point": point,
        "driver": driver,
        "mode": mode,
        "active_manifest_path": active_relative,
        "archived_manifest_path": archived_relative,
        "archived_manifest_sha256": sha256_file(archive_path),
        "current_manifest_sha256": sha256_file(manifest_path),
        "old_fingerprint": prior_fingerprint,
        "new_fingerprint": current_fingerprint,
    }
    record_netbid_manifest_migration(
        work_root=work_root, migration=migration, entry=entry
    )
    return record


def run_mode(
    *,
    work_root: Path,
    point: dict[str, Any],
    driver: str,
    driver_file: Path,
    driver_ids: list[str],
    prefix: str,
    mode: str,
    wrapper: Path,
    wrapper_hash: str,
    r_script: Path,
    r_script_hash: str,
    environment: dict[str, str],
    sweep_design_hash: str,
    design_migration: dict[str, Any] | None,
) -> dict[str, Any]:
    arm_root = point["root"] / driver
    consensus = arm_root / "consensus/consensus_network_ncol_.txt"
    consensus_manifest_path = arm_root / "consensus_manifest.json"
    if not consensus.is_file() or not consensus_manifest_path.is_file():
        raise RuntimeError(f"Consensus is incomplete for {point['key']}/{driver}")
    consensus_manifest = load_json(consensus_manifest_path)
    consensus_hash = sha256_file(consensus)
    if (
        consensus_manifest.get("stage") != point["key"]
        or consensus_manifest.get("driver") != driver
        or consensus_manifest.get("ncol", {}).get("sha256") != consensus_hash
    ):
        raise RuntimeError(f"Consensus provenance mismatch for {point['key']}/{driver}")
    expected_edges = int(consensus_manifest["ncol"]["edges"])

    mode_suffix = "netbid2_qc" if mode == "summary" else "netbid2_qc_html"
    output_root = arm_root / mode_suffix
    partial_root = arm_root / f"{mode_suffix}.partial"
    manifest_path = arm_root / f"{mode_suffix}_manifest.json"
    pending_path = arm_root / f"{mode_suffix}_manifest.pending.json"
    log_root = arm_root / "logs"
    log_root.mkdir(parents=True, exist_ok=True)
    stdout_path = log_root / f"{mode_suffix}.stdout.log"
    stderr_path = log_root / f"{mode_suffix}.stderr.log"

    fingerprint_payload = {
        "schema": RUN_SCHEMA,
        "mode": mode,
        "point": point["key"],
        "p_value": point["p_value"],
        "mi_cutoff": point["mi_cutoff"],
        "point_manifest_sha256": sha256_file(point["manifest_path"]),
        "sweep_design_sha256": sweep_design_hash,
        "driver": driver,
        "driver_sha256": sha256_file(driver_file),
        "consensus_sha256": consensus_hash,
        "consensus_manifest_sha256": sha256_file(consensus_manifest_path),
        "r_script_sha256": r_script_hash,
        "wrapper_sha256": wrapper_hash,
        "environment": environment,
        "prefix": prefix,
    }
    run_fingerprint = fingerprint(fingerprint_payload)
    command = [
        str(wrapper),
        "Rscript",
        str(r_script),
        str(consensus),
        str(driver_file),
        str(partial_root),
        prefix,
        "true" if mode == "html" else "false",
    ]

    if manifest_path.is_file() and output_root.is_dir():
        if partial_root.exists():
            raise RuntimeError(
                f"Completed and partial NetBID2 outputs both exist: {arm_root}"
            )
        record = load_json(manifest_path)
        if pending_path.exists():
            pending_record = load_json(pending_path)
            if pending_record != record:
                raise RuntimeError(
                    f"Pending/completed NetBID2 manifests disagree: {pending_path}"
                )
        record = migrate_completed_manifest_if_eligible(
            work_root=work_root,
            migration=design_migration,
            point=point["key"],
            driver=driver,
            mode=mode,
            manifest_path=manifest_path,
            output_root=output_root,
            record=record,
            expected_payload=fingerprint_payload,
            expected_command=command,
            prefix=prefix,
            driver_ids=driver_ids,
            expected_edges=expected_edges,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            expected_environment=environment,
        )
        if pending_path.exists():
            pending_path.unlink()
        print(f"[NETBID2 {mode.upper()}] {point['key']}/{driver} resume", flush=True)
        return record

    if manifest_path.exists() != output_root.exists():
        if manifest_path.exists():
            raise RuntimeError(f"Manifest exists without output: {manifest_path}")
        if not pending_path.is_file():
            raise RuntimeError(f"Unverifiable orphan NetBID2 output: {output_root}")

    if pending_path.is_file() and (output_root.is_dir() or partial_root.is_dir()):
        if output_root.is_dir() and partial_root.exists():
            raise RuntimeError(f"Both final and partial NetBID2 outputs exist: {arm_root}")
        recovery_root = output_root if output_root.is_dir() else partial_root
        record = load_json(pending_path)
        validate_record(
            root=recovery_root,
            record=record,
            expected_fingerprint=run_fingerprint,
            mode=mode,
            prefix=prefix,
            driver_ids=driver_ids,
            expected_edges=expected_edges,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            expected_environment=environment,
        )
        if recovery_root == partial_root:
            os.replace(partial_root, output_root)
        record["recovered_after_interrupted_manifest"] = True
        # Keep pending and completed records byte-equivalent across the final
        # narrow crash window so a later M+O+P resume can verify before unlink.
        atomic_json(pending_path, record)
        atomic_json(manifest_path, record)
        pending_path.unlink()
        print(f"[NETBID2 {mode.upper()}] {point['key']}/{driver} recovered", flush=True)
        return record

    if output_root.exists() or manifest_path.exists() or pending_path.exists():
        raise RuntimeError(f"Unverifiable NetBID2 state under {arm_root}")
    if partial_root.exists():
        shutil.rmtree(partial_root)
    print(f"[NETBID2 {mode.upper()}] {point['key']}/{driver}", flush=True)
    with stdout_path.open("w", encoding="utf-8", newline="\n") as stdout, \
            stderr_path.open("w", encoding="utf-8", newline="\n") as stderr:
        subprocess.run(command, stdout=stdout, stderr=stderr, check=True)
    output_inventory = validate_output(
        partial_root,
        mode=mode,
        prefix=prefix,
        driver_ids=driver_ids,
        expected_edges=expected_edges,
        expected_environment=environment,
    )
    record = {
        **fingerprint_payload,
        "fingerprint": run_fingerprint,
        "command": command,
        "finished_at_utc": utc_now(),
        "output": str(output_root),
        "output_inventory": output_inventory,
        "stdout_sha256": sha256_file(stdout_path),
        "stderr_sha256": sha256_file(stderr_path),
        "stderr_bytes": stderr_path.stat().st_size,
    }
    atomic_json(pending_path, record)
    os.replace(partial_root, output_root)
    atomic_json(manifest_path, record)
    pending_path.unlink()
    return record


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
    parser.add_argument(
        "--html-points",
        default="none",
        help="Comma-separated selected point keys, all, or none (default)",
    )
    return parser.parse_args()


def finalized_aggregate(payload: dict[str, Any]) -> dict[str, Any]:
    result = dict(payload)
    result["fingerprint"] = fingerprint(result)
    return result


def main() -> int:
    args = parse_args()
    results_root = args.work_root / "results"
    discovered, sweep_design, sweep_design_hash = discover_points(args.work_root)
    design_migration = load_sweep_design_migration(
        work_root=args.work_root,
        sweep_design=sweep_design,
        sweep_design_hash=sweep_design_hash,
    )
    point_map = {point["key"]: point for point in discovered}
    available_points = [point["key"] for point in discovered]
    point_keys = select_values(args.points, available_points, "point")
    driver_keys = select_values(args.drivers, list(DRIVERS), "driver")
    html_points = select_values(
        args.html_points, available_points, "HTML point", allow_none=True
    )
    outside_selection = sorted(set(html_points) - set(point_keys))
    if outside_selection:
        raise ValueError(
            "--html-points must be a subset of --points: "
            + ", ".join(outside_selection)
        )
    complete_selection = (
        point_keys == available_points and driver_keys == list(DRIVERS)
    )
    if html_points and not complete_selection:
        raise ValueError(
            "HTML generation requires --points all --drivers all so every "
            "stable summary record is revalidated first"
        )

    wrapper = args.benchmark_repo / "benchmarks/brca100_netbid_qc/netbid2-r"
    r_script = (
        args.benchmark_repo
        / "benchmarks/brca100_pr67_threshold_sweep/run_netbid_qc.R"
    )
    for required in (wrapper, r_script):
        if not required.is_file():
            raise FileNotFoundError(required)
    environment, probe_stderr = probe_environment(wrapper, r_script)
    if probe_stderr:
        print("[NETBID2 PROBE STDERR] " + probe_stderr.rstrip(), flush=True)
    wrapper_hash = sha256_file(wrapper)
    r_script_hash = sha256_file(r_script)

    summary_records: list[dict[str, Any]] = []
    summary_by_arm: dict[tuple[str, str], dict[str, Any]] = {}
    arm_context: dict[tuple[str, str], tuple[Path, list[str], str]] = {}
    for key in point_keys:
        point = point_map[key]
        inputs = point["manifest"].get("inputs", {})
        for driver in driver_keys:
            filename, prefix, expected_count = DRIVERS[driver]
            driver_file = args.work_root / "inputs" / filename
            if not driver_file.is_file():
                raise FileNotFoundError(driver_file)
            driver_ids = [
                line.strip()
                for line in driver_file.read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
            if len(driver_ids) != expected_count or len(set(driver_ids)) != len(driver_ids):
                raise ValueError(f"Unexpected {driver} driver list")
            expected_input_hash = inputs.get(filename, {}).get("sha256")
            if expected_input_hash != sha256_file(driver_file):
                raise RuntimeError(
                    f"Driver input does not match point manifest: {key}/{driver}"
                )
            summary = run_mode(
                work_root=args.work_root,
                point=point,
                driver=driver,
                driver_file=driver_file,
                driver_ids=driver_ids,
                prefix=prefix,
                mode="summary",
                wrapper=wrapper,
                wrapper_hash=wrapper_hash,
                r_script=r_script,
                r_script_hash=r_script_hash,
                environment=environment,
                sweep_design_hash=sweep_design_hash,
                design_migration=design_migration,
            )
            summary_records.append(summary)
            summary_by_arm[(key, driver)] = summary
            arm_context[(key, driver)] = (driver_file, driver_ids, prefix)

    summary_aggregate = finalized_aggregate({
        "schema": SUMMARY_AGGREGATE_SCHEMA,
        "environment": environment,
        "sweep_design_sha256": sweep_design_hash,
        "all_sweep_points": [item["key"] for item in discovered],
        "selection": {
            "points": point_keys,
            "drivers": driver_keys,
        },
        "summary_runs": summary_records,
    })
    summary_aggregate_path = results_root / "netbid2_qc_manifest.json"
    if html_points:
        if not summary_aggregate_path.is_file():
            raise RuntimeError(
                "Full summary aggregate is missing; first run with "
                "--points all --drivers all --html-points none"
            )
        existing_summary_aggregate = load_json(summary_aggregate_path)
        if existing_summary_aggregate != summary_aggregate:
            raise RuntimeError(
                "Stable NetBID2 summary aggregate is stale or incomplete; "
                "rerun summary-only QC before HTML generation"
            )
    elif complete_selection:
        if (
            not summary_aggregate_path.is_file()
            or load_json(summary_aggregate_path) != summary_aggregate
        ):
            atomic_json(summary_aggregate_path, summary_aggregate)
    else:
        print(
            "[NETBID2] partial summary selection: stable root aggregate unchanged",
            flush=True,
        )

    html_records: list[dict[str, Any]] = []
    for key in html_points:
        point = point_map[key]
        for driver in driver_keys:
            driver_file, driver_ids, prefix = arm_context[(key, driver)]
            html = run_mode(
                work_root=args.work_root,
                point=point,
                driver=driver,
                driver_file=driver_file,
                driver_ids=driver_ids,
                prefix=prefix,
                mode="html",
                wrapper=wrapper,
                wrapper_hash=wrapper_hash,
                r_script=r_script,
                r_script_hash=r_script_hash,
                environment=environment,
                sweep_design_hash=sweep_design_hash,
                design_migration=design_migration,
            )
            summary = summary_by_arm[(key, driver)]
            for filename_to_match in (
                "network_summary.tsv",
                "driver_target_sizes.tsv",
                "netbid_environment.tsv",
            ):
                summary_path = Path(summary["output"]) / filename_to_match
                html_path = Path(html["output"]) / filename_to_match
                if sha256_file(summary_path) != sha256_file(html_path):
                    raise RuntimeError(
                        f"Summary/full HTML mismatch: {key}/{driver}/"
                        f"{filename_to_match}"
                    )
            html_records.append(html)

    if html_points:
        html_aggregate = finalized_aggregate({
            "schema": HTML_AGGREGATE_SCHEMA,
            "environment": environment,
            "sweep_design_sha256": sweep_design_hash,
            "all_sweep_points": [item["key"] for item in discovered],
            "selection": {
                "points": point_keys,
                "drivers": driver_keys,
                "html_points": html_points,
            },
            "html_runs": html_records,
        })
        atomic_json(
            results_root / "netbid2_qc_html_manifest.json", html_aggregate
        )
    print(
        f"[NETBID2] complete: {len(summary_records)} summary arm(s), "
        f"{len(html_records)} HTML arm(s)",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
