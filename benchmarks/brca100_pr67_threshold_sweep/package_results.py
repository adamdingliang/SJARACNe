#!/usr/bin/env python3
"""Create a compact, checksummed package from a completed PR67 p-value sweep.

The packager is deliberately fail-closed.  It accepts only the exact
thirteen-point, two-driver BRCA100 design, a completed analysis with PR66
anchor evidence, all 26 consensus/support/NetBID2 summary arms, and an
immutable full-run manifest.
Large scientific artifacts are never copied.  Their analysis-validated hashes
are rechecked against the current bytes and recorded in
``omitted_artifacts.json`` instead.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import shutil
from typing import Any


DRIVERS = ("tf", "sig")
NETBID_DRIVERS = {
    "tf": ("BRCA100_TF.txt", "TF_"),
    "sig": ("BRCA100_SIG.txt", "SIG_"),
}
NETBID_ENVIRONMENT_KEYS = frozenset(
    {"R", "NetBID2", "NetBID2_remote_sha", "igraph"}
)
NETBID_SHARED_OUTPUTS = (
    "driver_target_sizes.tsv",
    "netbid_environment.tsv",
    "network_summary.tsv",
)
POINTS = (
    ("p1e-07", 1e-7),
    ("p1e-06", 1e-6),
    ("p1e-05", 1e-5),
    ("p2e-05", 2e-5),
    ("p5e-05", 5e-5),
    ("p1e-04", 1e-4),
    ("p2e-04", 2e-4),
    ("p3e-04", 3e-4),
    ("p_pr66_cutoff_match", 0.000352804562601613),
    ("p4e-04", 4e-4),
    ("p5e-04", 5e-4),
    ("p7p5e-04", 7.5e-4),
    ("p1e-03", 1e-3),
)
POINT_KEYS = tuple(key for key, _value in POINTS)
POINT_P = dict(POINTS)
LEGACY_POINT_KEYS = POINT_KEYS[:9]
SEEDS = tuple(range(1, 101))
PR67_COMMIT = "7633ebb4a0d966dbda15a4e32d0efa492fb71aeb"
NULL_MODEL_SHA256 = (
    "e3a8522682a8ea239821aaa10b12db72d00e07bfdcad43599d8e76a06be80944"
)
EXPECTED_INPUT_SHA256 = {
    "BRCA100.exp": "ad8a334f5f8cdf46a1000d3ee259b35258a18b3da2e314bb3a0cf7a421d98bc8",
    "BRCA100_TF.txt": "9b1219a489b99432175e4c4ad46add7b06f25aae388ee8dd3261fa91e4c43ffd",
    "BRCA100_SIG.txt": "16ca27df655f16684f880a4ad719c4e2ae3f8dc0d7e6b9eccdd24cd97c40797c",
}
ANALYSIS_OUTPUTS = frozenset(
    {
        "seed_manifest_summary.tsv",
        "point_manifest_summary.tsv",
        "network_summary.tsv",
        "adjacent_overlap.tsv",
        "anchor_overlap.tsv",
        "operating_point_screen.tsv",
        "selection.json",
        "arm_provenance.tsv",
        "pr66_context_summary.tsv",
        "pr66_context_overlap.tsv",
        "plots/core_metrics_vs_log10_p.png",
        "plots/core_metrics_vs_log10_p.svg",
        "plots/edge_overlap_vs_log10_p.png",
        "plots/edge_overlap_vs_log10_p.svg",
        "plots/coverage_vs_nominal_null_burden.png",
        "plots/coverage_vs_nominal_null_burden.svg",
    }
)
SHA256_HEX = frozenset("0123456789abcdef")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_directory(path: Path) -> str:
    require_directory(path)
    digest = hashlib.sha256()
    files = sorted(item for item in path.rglob("*") if item.is_file())
    for child in files:
        if child.is_symlink():
            raise ValueError(f"Symlink is not accepted in hashed directory: {child}")
        digest.update(child.relative_to(path).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(bytes.fromhex(sha256_file(child)))
    return digest.hexdigest()


def json_fingerprint(payload: object) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def atomic_json(path: Path, payload: object) -> None:
    temporary = path.with_name(path.name + ".partial")
    data = (
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False)
        + "\n"
    ).encode("utf-8")
    temporary.write_bytes(data)
    os.replace(temporary, path)


def load_json(path: Path) -> dict[str, Any]:
    require_file(path)
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError(f"Cannot read JSON object {path}: {error}") from error
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def require_file(path: Path) -> Path:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Missing/non-regular source file: {path}")
    return path


def require_directory(path: Path) -> Path:
    if path.is_symlink() or not path.is_dir():
        raise ValueError(f"Missing/non-regular source directory: {path}")
    return path


def require_sha256(value: Any, field: str, source: Path) -> str:
    if not isinstance(value, str) or len(value) != 64 or any(
        character not in SHA256_HEX for character in value
    ):
        raise ValueError(f"Invalid {field} in {source}: {value!r}")
    return value


def safe_relative(value: Any, field: str, source: Path) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"Invalid {field} in {source}: {value!r}")
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or "." in path.parts:
        raise ValueError(f"Unsafe {field} in {source}: {value!r}")
    normalized = path.as_posix()
    if normalized != value.replace("\\", "/"):
        raise ValueError(f"Non-canonical {field} in {source}: {value!r}")
    return normalized


def same_path(recorded: Any, expected: Path, field: str, source: Path) -> None:
    if not isinstance(recorded, str) or Path(recorded).resolve() != expected.resolve():
        raise ValueError(
            f"{field} path mismatch in {source}: {recorded!r} != {str(expected)!r}"
        )


def work_relative(path: Path, work_root: Path) -> str:
    try:
        return path.resolve().relative_to(work_root).as_posix()
    except ValueError as error:
        raise ValueError(f"Artifact is outside work root: {path}") from error


def set_digest(items: list[tuple[str, str]]) -> str:
    digest = hashlib.sha256()
    for relative, checksum in items:
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(bytes.fromhex(checksum))
    return digest.hexdigest()


class PackageEvidence:
    def __init__(self, work_root: Path) -> None:
        self.work_root = work_root
        self.copy_plan: dict[str, tuple[Path, str]] = {}
        self.omitted_files: dict[str, dict[str, Any]] = {}
        self.omitted_sets: dict[tuple[str, str], dict[str, Any]] = {}

    def copy(self, source: Path, destination: str, role: str) -> None:
        require_file(source)
        destination = safe_relative(destination, "package destination", source)
        if destination in self.copy_plan:
            raise ValueError(f"Duplicate package destination: {destination}")
        self.copy_plan[destination] = (source, role)

    def omit_file(
        self,
        path: Path,
        category: str,
        checksum: str,
        *,
        verification: str,
        expected_bytes: int | None = None,
    ) -> None:
        require_file(path)
        checksum = require_sha256(checksum, "omitted artifact SHA256", path)
        size = path.stat().st_size
        if expected_bytes is not None and size != expected_bytes:
            raise ValueError(
                f"Omitted artifact size mismatch for {path}: {size} != {expected_bytes}"
            )
        actual_checksum = sha256_file(path)
        if path.stat().st_size != size:
            raise ValueError(f"Omitted artifact changed size while hashing: {path}")
        if actual_checksum != checksum:
            raise ValueError(
                f"Omitted artifact hash mismatch for {path}: "
                f"{actual_checksum} != {checksum}"
            )
        relative = work_relative(path, self.work_root)
        record = {
            "bytes": size,
            "category": category,
            "path": relative,
            "sha256": checksum,
            "verification": verification,
        }
        previous = self.omitted_files.get(relative)
        if previous is not None and previous != record:
            raise ValueError(f"Conflicting omitted-artifact evidence: {relative}")
        self.omitted_files[relative] = record

    def omit_rehashed(self, path: Path, category: str) -> None:
        self.omit_file(
            path,
            category,
            sha256_file(require_file(path)),
            verification="rehashed-by-packager",
        )

    def omit_set(
        self,
        root: Path,
        category: str,
        checksum: str,
        file_count: int,
        verification: str,
    ) -> None:
        checksum = require_sha256(checksum, "omitted set SHA256", root)
        relative = work_relative(root, self.work_root)
        key = (category, relative)
        record = {
            "category": category,
            "file_count": file_count,
            "root": relative,
            "sha256": checksum,
            "verification": verification,
        }
        previous = self.omitted_sets.get(key)
        if previous is not None and previous != record:
            raise ValueError(f"Conflicting omitted-set evidence: {relative}")
        self.omitted_sets[key] = record


def validate_analysis(work_root: Path) -> tuple[dict[str, Any], Path]:
    analysis_root = work_root / "results" / "analysis"
    require_directory(analysis_root)
    if analysis_root.with_name("analysis.partial").exists():
        raise ValueError("Unresolved analysis.partial directory")
    manifest_path = analysis_root / "analysis_manifest.json"
    manifest = load_json(manifest_path)
    if manifest.get("schema") != "sjaracne-brca100-pr67-threshold-sweep-analysis-v1":
        raise ValueError(f"Unexpected analysis schema: {manifest_path}")
    same_path(manifest.get("work_root"), work_root, "analysis work_root", manifest_path)
    same_path(manifest.get("output_root"), analysis_root, "analysis output_root", manifest_path)
    if manifest.get("p_keys_in_increasing_p_order") != list(POINT_KEYS):
        raise ValueError(
            "Analysis does not cover the exact ordered thirteen-point sweep"
        )
    if int(manifest.get("matched_seed_count", -1)) != len(SEEDS):
        raise ValueError(f"Analysis does not record {len(SEEDS)} matched seeds")
    if manifest.get("build", {}).get("commit") != PR67_COMMIT:
        raise ValueError("Analysis was not produced from the pinned PR67 commit")
    if manifest.get("build", {}).get("null_model_sha256") != NULL_MODEL_SHA256:
        raise ValueError("Analysis null-model hash is not the pinned m=80 model")

    outputs = manifest.get("output_files")
    if not isinstance(outputs, dict) or set(outputs) != ANALYSIS_OUTPUTS:
        raise ValueError(
            "Analysis output inventory is incomplete or unexpected: "
            f"missing={sorted(ANALYSIS_OUTPUTS - set(outputs or {}))}, "
            f"extra={sorted(set(outputs or {}) - ANALYSIS_OUTPUTS)}"
        )
    actual = {
        item.relative_to(analysis_root).as_posix()
        for item in analysis_root.rglob("*")
        if item.is_file()
    }
    expected_actual = set(outputs) | {"analysis_manifest.json"}
    if actual != expected_actual:
        raise ValueError(
            "Analysis directory contains missing/unlisted files: "
            f"missing={sorted(expected_actual - actual)}, extra={sorted(actual - expected_actual)}"
        )
    for relative in sorted(outputs):
        safe_relative(relative, "analysis output", manifest_path)
        path = analysis_root / relative
        require_file(path)
        record = outputs[relative]
        if not isinstance(record, dict):
            raise ValueError(f"Invalid analysis output record: {relative}")
        expected_hash = require_sha256(record.get("sha256"), relative, manifest_path)
        if path.stat().st_size != int(record.get("bytes", -1)):
            raise ValueError(f"Analysis output size changed: {path}")
        if sha256_file(path) != expected_hash:
            raise ValueError(f"Analysis output hash changed: {path}")
    if load_json(analysis_root / "selection.json") != manifest.get(
        "operating_point_selection"
    ):
        raise ValueError("selection.json disagrees with analysis_manifest.json")

    arms = manifest.get("arms")
    expected_arms = {f"{point}/{driver}" for point in POINT_KEYS for driver in DRIVERS}
    if not isinstance(arms, dict) or set(arms) != expected_arms:
        raise ValueError(
            "Analysis arm provenance does not exactly cover "
            f"{len(POINT_KEYS) * len(DRIVERS)} arms"
        )
    required_arm_hashes = {
        "adjacency_set_sha256",
        "seed_metadata_set_sha256",
        "run_manifest_rows_sha256",
        "point_manifest_sha256",
        "consensus_sha256",
        "consensus_manifest_sha256",
        "consensus_fingerprint",
        "support_sha256",
        "support_manifest_sha256",
        "support_fingerprint",
        "support_source_sha256",
        "support_binary_sha256",
        "netbid2_qc_set_sha256",
        "netbid2_network_summary_sha256",
        "netbid2_manifest_sha256",
        "netbid2_fingerprint",
    }
    for arm, record in arms.items():
        if not isinstance(record, dict) or not required_arm_hashes.issubset(record):
            raise ValueError(f"Incomplete analysis provenance for {arm}")
        for field in required_arm_hashes:
            require_sha256(record[field], f"{arm}.{field}", manifest_path)

    context = manifest.get("pr66_context")
    if not isinstance(context, dict) or context.get("cutoff_match_consensus_exact") is not True:
        raise ValueError("Completed PR66 context and exact cutoff-match evidence are required")
    anchor = context.get("anchor_seed_equivalence")
    if not isinstance(anchor, dict) or (
        int(anchor.get("comparisons", -1)) != 2 * len(DRIVERS) * len(SEEDS)
        or anchor.get("all_data_sections_equal") is not True
    ):
        raise ValueError("Analysis lacks completed full anchor evidence")
    return manifest, analysis_root


def validate_design_and_build(
    work_root: Path, analysis: dict[str, Any], evidence: PackageEvidence
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    design_path = work_root / "sweep_design.json"
    design = load_json(design_path)
    design_hash = sha256_file(design_path)
    if design_hash != analysis.get("sweep_design_sha256"):
        raise ValueError("sweep_design.json changed after analysis")
    if (
        design.get("schema") != "sjaracne-brca100-pr67-p-sweep-v2"
        or design.get("commit") != PR67_COMMIT
        or design.get("null_model_sha256") != NULL_MODEL_SHA256
    ):
        raise ValueError("Unexpected sweep design identity")
    fixed = design.get("fixed_parameters")
    expected_fixed = {
        "sampling": "fixed 80% without replacement",
        "m": 80,
        "npar": 40,
        "dpi_epsilon": 0,
        "consensus_p": 1e-5,
        "seeds": list(SEEDS),
    }
    if fixed != expected_fixed:
        raise ValueError("Sweep fixed parameters changed")
    points = design.get("all_points")
    if not isinstance(points, list) or [point.get("key") for point in points] != list(
        POINT_KEYS
    ):
        raise ValueError(
            "Sweep design does not contain the exact ordered thirteen points"
        )
    for point in points:
        key = point["key"]
        if not math.isclose(
            float(point.get("p_value")), POINT_P[key], rel_tol=1e-14, abs_tol=0.0
        ):
            raise ValueError(f"Sweep p-value changed for {key}")

    inputs = design.get("inputs")
    if not isinstance(inputs, dict):
        raise ValueError("Missing sweep input provenance")
    for filename, expected_hash in EXPECTED_INPUT_SHA256.items():
        record = inputs.get(filename)
        path = work_root / "inputs" / filename
        if not isinstance(record, dict):
            raise ValueError(f"Missing input metadata for {filename}")
        same_path(record.get("path"), path, f"input {filename}", design_path)
        if (
            record.get("sha256") != expected_hash
            or int(record.get("bytes", -1)) != require_file(path).stat().st_size
            or sha256_file(path) != expected_hash
        ):
            raise ValueError(f"Input file changed: {path}")
        evidence.omit_file(
            path,
            "pinned-input",
            expected_hash,
            verification="rehashed-by-packager",
            expected_bytes=int(record["bytes"]),
        )

    build_path = work_root / "builds" / "pr67_7633ebb" / "build_manifest.json"
    build = load_json(build_path)
    if sha256_file(build_path) != analysis.get("build", {}).get("build_manifest_sha256"):
        raise ValueError("Build manifest changed after analysis")
    binary = work_root / "builds" / "pr67_7633ebb" / "bin" / "sjaracne.exe"
    config = (
        work_root / "builds" / "pr67_7633ebb" / "source" / "SJARACNe" / "config"
    )
    model = config / "apmi_null" / "apmi_null_m00080_npar040.model"
    if build.get("stage") != "pr67_7633ebb" or build.get("commit") != PR67_COMMIT:
        raise ValueError("Unexpected build manifest identity")
    same_path(build.get("binary"), binary, "build binary", build_path)
    same_path(build.get("config_directory"), config, "build config", build_path)
    same_path(build.get("null_model"), model, "build null model", build_path)
    binary_hash = sha256_file(require_file(binary))
    model_hash = sha256_file(require_file(model))
    config_hash = sha256_directory(config)
    checks = {
        "binary_sha256": binary_hash,
        "config_sha256": config_hash,
        "null_model_sha256": model_hash,
    }
    for field, observed in checks.items():
        if build.get(field) != observed or design.get(field) != observed:
            raise ValueError(f"Build/design {field} mismatch")
    if model_hash != NULL_MODEL_SHA256:
        raise ValueError("Built null model is not the pinned m=80 model")
    evidence.omit_file(
        binary, "build-binary", binary_hash, verification="rehashed-by-packager"
    )
    evidence.omit_file(
        model, "null-model", model_hash, verification="rehashed-by-packager"
    )
    config_files = [item for item in sorted(config.rglob("*")) if item.is_file()]
    evidence.omit_set(
        config,
        "build-config-tree",
        config_hash,
        len(config_files),
        "rehashed-by-packager",
    )
    evidence.copy(design_path, "provenance/sweep_design.json", "sweep-design")

    history_root = work_root / "sweep_design_history"
    migration_paths = sorted(history_root.glob("*.migration.json"))
    archive_paths = sorted(history_root.glob("*.sweep_design.json"))
    if len(migration_paths) != 1 or len(archive_paths) != 1:
        raise ValueError(
            "Extended sweep requires exactly one archived v1 design and migration record"
        )
    migration_path = migration_paths[0]
    archive_path = archive_paths[0]
    migration = load_json(migration_path)
    archived = load_json(archive_path)
    migration_fingerprint = migration.get("fingerprint")
    migration_payload = dict(migration)
    migration_payload.pop("fingerprint", None)
    archive_hash = sha256_file(archive_path)
    archived_invariants = {
        key: value
        for key, value in archived.items()
        if key not in {"schema", "all_points"}
    }
    active_invariants = {
        key: value
        for key, value in design.items()
        if key not in {"schema", "all_points"}
    }
    expected_history = {
        "schema": "sjaracne-brca100-pr67-p-sweep-design-migration-v1",
        "operation": "append-only-point-extension",
        "from_schema": "sjaracne-brca100-pr67-p-sweep-v1",
        "to_schema": "sjaracne-brca100-pr67-p-sweep-v2",
    }
    if (
        migration.get("schema") != expected_history["schema"]
        or migration.get("operation") != expected_history["operation"]
        or not isinstance(migration.get("from"), dict)
        or not isinstance(migration.get("to"), dict)
        or migration["from"].get("schema") != expected_history["from_schema"]
        or migration["to"].get("schema") != expected_history["to_schema"]
        or migration["from"].get("sha256") != archive_hash
        or migration["to"].get("sha256") != design_hash
        or migration["from"].get("archived_path")
        != archive_path.relative_to(work_root).as_posix()
        or migration["to"].get("active_path") != "sweep_design.json"
        or migration.get("manifest_path")
        != migration_path.relative_to(work_root).as_posix()
        or migration_fingerprint != json_fingerprint(migration_payload)
        or archived.get("schema") != expected_history["from_schema"]
        or archived.get("all_points") != points[: len(LEGACY_POINT_KEYS)]
        or archived_invariants != active_invariants
        or migration["from"].get("point_keys") != list(LEGACY_POINT_KEYS)
        or migration["to"].get("point_keys") != list(POINT_KEYS)
        or migration.get("appended_points") != points[len(LEGACY_POINT_KEYS) :]
    ):
        raise ValueError("Invalid append-only sweep-design migration evidence")
    evidence.copy(
        archive_path,
        f"provenance/sweep_design_history/{archive_path.name}",
        "archived-sweep-design",
    )
    evidence.copy(
        migration_path,
        f"provenance/sweep_design_history/{migration_path.name}",
        "sweep-design-migration",
    )
    evidence.copy(
        build_path,
        "provenance/builds/pr67_7633ebb/build_manifest.json",
        "build-manifest",
    )

    manifests = sorted((work_root / "results").glob("*/point_manifest.json"))
    if [path.parent.name for path in manifests] != sorted(POINT_KEYS):
        raise ValueError(
            "Point manifests do not exactly cover the thirteen-point sweep"
        )
    point_records: list[dict[str, Any]] = []
    design_by_key = {point["key"]: point for point in points}
    for key in POINT_KEYS:
        path = work_root / "results" / key / "point_manifest.json"
        record = load_json(path)
        design_point = design_by_key[key]
        if (
            record.get("schema") != "sjaracne-brca100-pr67-p-sweep-point-v1"
            or record.get("key") != key
            or record.get("commit") != PR67_COMMIT
            or record.get("m") != 80
            or record.get("npar") != 40
            or record.get("dpi_epsilon") != 0
            or record.get("consensus_p") != 1e-5
            or record.get("seeds") != list(SEEDS)
            or record.get("sampling") != "fixed 80% without replacement"
            or record.get("p_value") != design_point.get("p_value")
            or record.get("mi_cutoff") != design_point.get("mi_cutoff")
            or record.get("inputs") != inputs
        ):
            raise ValueError(f"Point manifest disagrees with sweep design: {path}")
        manifest_hash = sha256_file(path)
        for driver in DRIVERS:
            if analysis["arms"][f"{key}/{driver}"]["point_manifest_sha256"] != manifest_hash:
                raise ValueError(f"Point manifest changed after analysis: {path}")
        evidence.copy(
            path,
            f"provenance/points/{key}/point_manifest.json",
            "point-manifest",
        )
        point_records.append(record)
    return design, point_records


def validate_invocations(work_root: Path, evidence: PackageEvidence) -> None:
    path = work_root / "invocations.json"
    value = load_json(path)
    if value.get("schema") != "sjaracne-brca100-pr67-p-sweep-invocations-v1":
        raise ValueError("Unexpected invocation-history schema")
    invocations = value.get("invocations")
    if not isinstance(invocations, list) or not invocations:
        raise ValueError("Empty invocation history")

    def complete_phase(invocation: Any, phases: set[str]) -> bool:
        return (
            isinstance(invocation, dict)
            and invocation.get("status") == "complete"
            and invocation.get("phase") in phases
            and invocation.get("points") == list(POINT_KEYS)
            and invocation.get("drivers") == list(DRIVERS)
            and invocation.get("seed_start") == 1
            and invocation.get("seed_end") == 100
            and isinstance(invocation.get("finished_at_utc"), str)
        )

    inference = [item for item in invocations if complete_phase(item, {"infer", "all"})]
    consensus = [item for item in invocations if complete_phase(item, {"consensus", "all"})]
    if not inference or not consensus:
        raise ValueError("Invocation history lacks completed full inference and consensus runs")
    if not any(
        item.get("inference_jobs") == len(POINT_KEYS) * len(DRIVERS) * len(SEEDS)
        and item.get("inference_new_jobs", 0) + item.get("inference_resumed_jobs", 0)
        == len(POINT_KEYS) * len(DRIVERS) * len(SEEDS)
        for item in inference
    ):
        raise ValueError("Invocation history lacks a completed full-grid inference record")
    evidence.copy(path, "provenance/invocations.json", "invocation-history")


def validate_run_manifest(
    work_root: Path,
    analysis: dict[str, Any],
    point_records: list[dict[str, Any]],
    evidence: PackageEvidence,
) -> dict[tuple[str, str], list[dict[str, str]]]:
    path = work_root / "results" / "run_manifest.tsv"
    require_file(path)
    if sha256_file(path) != analysis.get("run_manifest_sha256"):
        raise ValueError("run_manifest.tsv changed after analysis")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = set(reader.fieldnames or [])
    required = {
        "point", "p_value", "mi_cutoff", "validation_class", "commit", "driver",
        "seed", "binary_sha256", "edges", "source_rows", "adjacency_bytes",
        "adjacency_sha256", "data_sha256", "stderr_bytes",
    }
    expected_job_count = len(POINT_KEYS) * len(DRIVERS) * len(SEEDS)
    if not required.issubset(fields) or len(rows) != expected_job_count:
        raise ValueError(
            f"run_manifest.tsv is not the exact {expected_job_count}-row completed manifest"
        )
    expected_order = [
        (point, driver, seed)
        for point in POINT_KEYS
        for driver in DRIVERS
        for seed in SEEDS
    ]
    expected_identities = set(expected_order)
    observed_order: list[tuple[str, str, int]] = []
    by_arm: dict[tuple[str, str], list[dict[str, str]]] = {
        (point, driver): [] for point in POINT_KEYS for driver in DRIVERS
    }
    point_by_key = {record["key"]: record for record in point_records}
    for row in rows:
        try:
            key = row["point"]
            driver = row["driver"]
            seed = int(row["seed"])
            adjacency_bytes = int(row["adjacency_bytes"])
            edges = int(row["edges"])
            source_rows = int(row["source_rows"])
            stderr_bytes = int(row["stderr_bytes"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("Invalid numeric/identity field in run_manifest.tsv") from error
        identity = (key, driver, seed)
        observed_order.append(identity)
        if identity not in expected_identities:
            raise ValueError(f"Unexpected run-manifest row: {identity}")
        point = point_by_key[key]
        if (
            row["commit"] != PR67_COMMIT
            or row["binary_sha256"] != analysis["build"]["binary_sha256"]
            or row["validation_class"] != point["validation_class"]
            or not math.isclose(float(row["p_value"]), float(point["p_value"]), rel_tol=1e-14)
            or not math.isclose(float(row["mi_cutoff"]), float(point["mi_cutoff"]), rel_tol=1e-14)
            or min(adjacency_bytes, edges, source_rows, stderr_bytes) < 0
        ):
            raise ValueError(f"Run-manifest provenance mismatch: {identity}")
        adjacency_hash = require_sha256(row["adjacency_sha256"], "adjacency_sha256", path)
        require_sha256(row["data_sha256"], "data_sha256", path)
        adjacency_path = (
            work_root / "results" / key / driver / "adjacency" / f"TF_run_{seed:03d}.adj"
        )
        evidence.omit_file(
            adjacency_path,
            "seed-adjacency",
            adjacency_hash,
            verification="rehashed-by-packager-against-run-manifest",
            expected_bytes=adjacency_bytes,
        )
        by_arm[(key, driver)].append(row)
    if observed_order != expected_order:
        raise ValueError("run_manifest.tsv ordering/coverage is not canonical")

    for key in POINT_KEYS:
        for driver in DRIVERS:
            arm_root = work_root / "results" / key / driver
            rows_arm = by_arm[(key, driver)]
            adjacency_set_items: list[tuple[str, str]] = []
            metadata_set_items: list[tuple[str, str]] = []
            adjacency_root = arm_root / "adjacency"
            metadata_root = arm_root / "seed_metadata"
            expected_adj_names = [f"TF_run_{seed:03d}.adj" for seed in SEEDS]
            expected_meta_names = [f"TF_run_{seed:03d}.json" for seed in SEEDS]
            if sorted(path.name for path in adjacency_root.iterdir()) != expected_adj_names:
                raise ValueError(f"Adjacency directory is not exactly seeds 1..100: {adjacency_root}")
            if sorted(path.name for path in metadata_root.iterdir()) != expected_meta_names:
                raise ValueError(f"Metadata directory is not exactly seeds 1..100: {metadata_root}")
            for row, seed in zip(rows_arm, SEEDS):
                stem = f"TF_run_{seed:03d}"
                adjacency_path = adjacency_root / f"{stem}.adj"
                adjacency_set_items.append(
                    (adjacency_path.relative_to(arm_root).as_posix(), row["adjacency_sha256"])
                )
                metadata_path = metadata_root / f"{stem}.json"
                metadata = load_json(metadata_path)
                metadata_hash = sha256_file(metadata_path)
                if (
                    metadata.get("schema") != "sjaracne-brca100-pr67-p-sweep-seed-v1"
                    or metadata.get("point", {}).get("key") != key
                    or metadata.get("driver") != driver
                    or metadata.get("seed") != seed
                    or metadata.get("adjacency", {}).get("full_sha256")
                    != row["adjacency_sha256"]
                    or int(metadata.get("adjacency", {}).get("bytes", -1))
                    != int(row["adjacency_bytes"])
                ):
                    raise ValueError(f"Seed metadata identity mismatch: {metadata_path}")
                metadata_set_items.append(
                    (metadata_path.relative_to(arm_root).as_posix(), metadata_hash)
                )
                evidence.omit_file(
                    metadata_path,
                    "seed-metadata",
                    metadata_hash,
                    verification="rehashed-by-packager",
                )
                for suffix, field in (("stdout.log", "stdout_sha256"), ("stderr.log", "stderr_sha256")):
                    log_path = arm_root / "logs" / f"{stem}.{suffix}"
                    expected_hash = require_sha256(metadata.get(field), field, metadata_path)
                    if sha256_file(require_file(log_path)) != expected_hash:
                        raise ValueError(f"Seed log changed after analysis: {log_path}")
                    evidence.omit_file(
                        log_path,
                        "seed-log",
                        expected_hash,
                        verification="rehashed-by-packager",
                    )
                time_path = arm_root / "logs" / f"{stem}.time.txt"
                evidence.omit_rehashed(time_path, "seed-timing-log")
            arm_analysis = analysis["arms"][f"{key}/{driver}"]
            adjacency_set_hash = set_digest(adjacency_set_items)
            metadata_set_hash = set_digest(metadata_set_items)
            if adjacency_set_hash != arm_analysis["adjacency_set_sha256"]:
                raise ValueError(f"Adjacency set digest changed for {key}/{driver}")
            if metadata_set_hash != arm_analysis["seed_metadata_set_sha256"]:
                raise ValueError(f"Seed-metadata set digest changed for {key}/{driver}")
            evidence.omit_set(
                adjacency_root,
                "seed-adjacency-set",
                adjacency_set_hash,
                len(SEEDS),
                "derived-from-analysis-validated-run-manifest",
            )
            evidence.omit_set(
                metadata_root,
                "seed-metadata-set",
                metadata_set_hash,
                len(SEEDS),
                "rehashed-by-packager",
            )
    evidence.copy(path, "provenance/run_manifest.tsv", "completed-run-manifest")
    return by_arm


def validate_anchor(
    work_root: Path, analysis: dict[str, Any], evidence: PackageEvidence
) -> None:
    root = work_root / "results" / "validation"
    table_path = root / "anchor_seed_equivalence.tsv"
    manifest_path = root / "anchor_seed_equivalence_manifest.json"
    manifest = load_json(manifest_path)
    context = analysis["pr66_context"]["anchor_seed_equivalence"]
    table_hash = sha256_file(require_file(table_path))
    manifest_hash = sha256_file(manifest_path)
    if (
        manifest.get("schema") != "sjaracne-brca100-pr67-p-sweep-anchor-equivalence-v1"
        or manifest.get("drivers") != list(DRIVERS)
        or manifest.get("seeds") != list(SEEDS)
        or manifest.get("comparisons") != 2 * len(DRIVERS) * len(SEEDS)
        or manifest.get("all_data_sections_equal") is not True
        or manifest.get("table_sha256") != table_hash
        or manifest.get("sweep_design_sha256") != analysis["sweep_design_sha256"]
        or context.get("table_sha256") != table_hash
        or context.get("manifest_sha256") != manifest_hash
    ):
        raise ValueError("Anchor-equivalence manifest changed or is incomplete")
    same_path(manifest.get("table"), table_path, "anchor table", manifest_path)
    with table_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        columns = tuple(reader.fieldnames or [])
    expected_columns = (
        "sweep_point", "prior_stage", "driver", "seed", "data_sha256", "edges",
        "source_rows", "sweep_metadata_sha256", "prior_metadata_sha256",
    )
    expected_comparisons = 2 * len(DRIVERS) * len(SEEDS)
    if columns != expected_columns or len(rows) != expected_comparisons:
        raise ValueError(
            f"Anchor-equivalence table is not the exact {expected_comparisons}-row schema"
        )
    stages = {"p1e-07": "pr67_7633ebb", "p_pr66_cutoff_match": "pr66_5809183"}
    expected = {
        (point, stage, driver, seed)
        for point, stage in stages.items()
        for driver in DRIVERS
        for seed in SEEDS
    }
    observed: set[tuple[str, str, str, int]] = set()
    for row in rows:
        try:
            item = (row["sweep_point"], row["prior_stage"], row["driver"], int(row["seed"]))
            if int(row["edges"]) < 0 or int(row["source_rows"]) < 0:
                raise ValueError
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("Invalid anchor-equivalence table row") from error
        if item in observed or item not in expected:
            raise ValueError(f"Duplicate/unexpected anchor-equivalence row: {item}")
        observed.add(item)
        for field in ("data_sha256", "sweep_metadata_sha256", "prior_metadata_sha256"):
            require_sha256(row[field], field, table_path)
    if observed != expected:
        raise ValueError("Anchor-equivalence table does not cover both anchors")
    evidence.copy(
        table_path,
        "validation/anchor_seed_equivalence.tsv",
        "anchor-validation-table",
    )
    evidence.copy(
        manifest_path,
        "validation/anchor_seed_equivalence_manifest.json",
        "anchor-validation-manifest",
    )


def validate_consensus_support(
    work_root: Path,
    analysis: dict[str, Any],
    by_arm: dict[tuple[str, str], list[dict[str, str]]],
    evidence: PackageEvidence,
) -> None:
    results = work_root / "results"
    support_binary = work_root / "tools" / "summarize_consensus_support"
    support_binary_hash = sha256_file(require_file(support_binary))
    support_source = (
        Path(__file__).resolve().parents[1]
        / "brca100_netbid_qc"
        / "summarize_consensus_support.cpp"
    )
    support_source_hash = sha256_file(require_file(support_source))
    expected_helper_hashes = {
        (
            analysis["arms"][f"{key}/{driver}"]["support_source_sha256"],
            analysis["arms"][f"{key}/{driver}"]["support_binary_sha256"],
        )
        for key in POINT_KEYS
        for driver in DRIVERS
    }
    if expected_helper_hashes != {(support_source_hash, support_binary_hash)}:
        raise ValueError("Support helper source/binary changed after analysis")
    evidence.omit_file(
        support_binary,
        "support-helper-binary",
        support_binary_hash,
        verification="rehashed-by-packager",
    )
    aggregate_path = results / "support_summary_manifest.json"
    aggregate = load_json(aggregate_path)
    if sha256_file(aggregate_path) != analysis.get("support_aggregate_manifest_sha256"):
        raise ValueError("Support aggregate changed after analysis")
    if (
        aggregate.get("schema") != "sjaracne-brca100-pr67-p-sweep-support-aggregate-v1"
        or aggregate.get("sweep_design_sha256") != analysis["sweep_design_sha256"]
        or aggregate.get("points") != list(POINT_KEYS)
        or aggregate.get("drivers") != list(DRIVERS)
    ):
        raise ValueError("Support aggregate does not describe the complete sweep")
    records = aggregate.get("records")
    expected_arm_count = len(POINT_KEYS) * len(DRIVERS)
    if not isinstance(records, list) or len(records) != expected_arm_count:
        raise ValueError(
            f"Support aggregate does not contain {expected_arm_count} arms"
        )
    expected_order = [(point, driver) for point in POINT_KEYS for driver in DRIVERS]
    observed_order = [(record.get("point"), record.get("driver")) for record in records]
    if observed_order != expected_order:
        raise ValueError("Support aggregate arm ordering/coverage is not canonical")
    aggregate_by_arm = dict(zip(expected_order, records))
    evidence.copy(
        aggregate_path,
        "provenance/aggregates/support_summary_manifest.json",
        "support-aggregate-manifest",
    )

    for key, driver in expected_order:
        arm_root = results / key / driver
        arm_analysis = analysis["arms"][f"{key}/{driver}"]
        pending = (
            arm_root / "consensus_manifest.pending.json",
            arm_root / "support_summary_manifest.pending.json",
        )
        if any(path.exists() for path in pending):
            raise ValueError(f"Unresolved consensus/support pending state: {key}/{driver}")
        consensus_manifest_path = arm_root / "consensus_manifest.json"
        consensus = load_json(consensus_manifest_path)
        if sha256_file(consensus_manifest_path) != arm_analysis["consensus_manifest_sha256"]:
            raise ValueError(f"Consensus manifest changed after analysis: {key}/{driver}")
        if (
            consensus.get("stage") != key
            or consensus.get("driver") != driver
            or consensus.get("fingerprint") != arm_analysis["consensus_fingerprint"]
        ):
            raise ValueError(f"Consensus manifest identity mismatch: {key}/{driver}")
        consensus_root = arm_root / "consensus"
        ncol = consensus.get("ncol")
        if not isinstance(ncol, dict):
            raise ValueError(f"Missing consensus ncol record: {key}/{driver}")
        consensus_files = (
            ("consensus_network_ncol_.txt", ncol.get("sha256"), ncol.get("bytes")),
            ("consensus_network_3col_.txt", consensus.get("consensus_3col_sha256"), None),
            ("parameter_info_.txt", consensus.get("parameter_info_sha256"), None),
            ("bootstrap_info_.txt", consensus.get("bootstrap_info_sha256"), None),
        )
        if ncol.get("sha256") != arm_analysis["consensus_sha256"]:
            raise ValueError(f"Consensus ncol hash changed: {key}/{driver}")
        for filename, checksum, byte_count in consensus_files:
            evidence.omit_file(
                consensus_root / filename,
                "consensus-output",
                require_sha256(checksum, filename, consensus_manifest_path),
                verification="rehashed-by-packager-against-consensus-manifest",
                expected_bytes=int(byte_count) if byte_count is not None else None,
            )
        for filename in ("consensus.stdout.log", "consensus.stderr.log", "consensus.time.txt"):
            evidence.omit_rehashed(arm_root / "logs" / filename, "consensus-log")

        support_manifest_path = arm_root / "support_summary_manifest.json"
        support = load_json(support_manifest_path)
        if (
            sha256_file(support_manifest_path) != arm_analysis["support_manifest_sha256"]
            or support != aggregate_by_arm[(key, driver)]
            or support.get("schema") != "sjaracne-brca100-pr67-p-sweep-support-v1"
            or support.get("point") != key
            or support.get("driver") != driver
            or support.get("fingerprint") != arm_analysis["support_fingerprint"]
            or support.get("output_sha256") != arm_analysis["support_sha256"]
            or support.get("consensus_sha256") != arm_analysis["consensus_sha256"]
            or support.get("source_sha256") != support_source_hash
            or support.get("binary_sha256") != support_binary_hash
            or support.get("point_manifest_sha256")
            != arm_analysis["point_manifest_sha256"]
            or support.get("consensus_manifest_sha256")
            != arm_analysis["consensus_manifest_sha256"]
            or support.get("adjacency_sha256")
            != [row["adjacency_sha256"] for row in by_arm[(key, driver)]]
        ):
            raise ValueError(f"Support manifest changed or is inconsistent: {key}/{driver}")
        support_path = consensus_root / "consensus_support.tsv"
        same_path(support.get("output"), support_path, "support output", support_manifest_path)
        evidence.omit_file(
            support_path,
            "consensus-support",
            support["output_sha256"],
            verification="rehashed-by-packager-against-support-manifest",
        )
        for filename in (
            "support_summary.stdout.log",
            "support_summary.stderr.log",
            "support_summary.time.txt",
        ):
            evidence.omit_rehashed(arm_root / "logs" / filename, "support-log")
        evidence.copy(
            consensus_manifest_path,
            f"provenance/arms/{key}/{driver}/consensus_manifest.json",
            "consensus-manifest",
        )
        evidence.copy(
            support_manifest_path,
            f"provenance/arms/{key}/{driver}/support_summary_manifest.json",
            "support-manifest",
        )


def pkg_network_hash(
    manifest: dict[str, Any], filename: str, manifest_path: Path
) -> str:
    inventory = manifest.get("output_inventory")
    if not isinstance(inventory, list):
        raise ValueError(f"Missing output inventory in {manifest_path}")
    matches = [record for record in inventory if record.get("path") == filename]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one {filename} in {manifest_path}")
    return require_sha256(matches[0].get("sha256"), filename, manifest_path)


def validate_inventory(
    root: Path,
    inventory: Any,
    evidence: PackageEvidence,
    category: str,
    manifest_path: Path,
) -> tuple[list[str], str]:
    require_directory(root)
    if not isinstance(inventory, list) or not inventory:
        raise ValueError(f"Empty output inventory in {manifest_path}")
    paths: list[str] = []
    for record in inventory:
        if not isinstance(record, dict) or set(record) != {"path", "bytes", "sha256"}:
            raise ValueError(f"Invalid output inventory record in {manifest_path}")
        relative = safe_relative(record.get("path"), "output inventory path", manifest_path)
        paths.append(relative)
        path = root / relative
        checksum = require_sha256(record.get("sha256"), relative, manifest_path)
        byte_count = record.get("bytes")
        if isinstance(byte_count, bool) or not isinstance(byte_count, int) or byte_count <= 0:
            raise ValueError(f"Invalid output inventory byte count: {path}")
        if path.stat().st_size != byte_count or sha256_file(require_file(path)) != checksum:
            raise ValueError(f"NetBID2 output changed: {path}")
        evidence.omit_file(
            path,
            category,
            checksum,
            verification="rehashed-by-packager",
            expected_bytes=byte_count,
        )
    if len(paths) != len(set(paths)):
        raise ValueError(f"NetBID2 output inventory is not unique: {manifest_path}")
    actual_paths = sorted(
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file()
    )
    if paths != actual_paths:
        raise ValueError(
            f"NetBID2 output inventory is incomplete/non-canonical in {manifest_path}: "
            f"recorded={paths}, actual={actual_paths}"
        )
    return paths, set_digest(
        [(str(record["path"]), str(record["sha256"])) for record in inventory]
    )


def files_equal(first: Path, second: Path) -> bool:
    """Compare two regular files exactly without loading either one in memory."""
    require_file(first)
    require_file(second)
    if first.stat().st_size != second.stat().st_size:
        return False
    with first.open("rb") as left, second.open("rb") as right:
        while True:
            left_chunk = left.read(1024 * 1024)
            right_chunk = right.read(1024 * 1024)
            if left_chunk != right_chunk:
                return False
            if not left_chunk:
                return True


def validate_netbid_environment(value: Any, source: Path) -> dict[str, str]:
    if (
        not isinstance(value, dict)
        or set(value) != NETBID_ENVIRONMENT_KEYS
        or any(not isinstance(item, str) or not item for item in value.values())
    ):
        raise ValueError(f"Incomplete NetBID2 environment in {source}: {value!r}")
    return value


def validate_netbid_run_contract(
    *,
    record: dict[str, Any],
    manifest_path: Path,
    mode: str,
    key: str,
    driver: str,
    work_root: Path,
    analysis: dict[str, Any],
    environment: dict[str, str],
    wrapper: Path,
    wrapper_sha256: str,
    r_script: Path,
    r_script_sha256: str,
) -> None:
    """Rebuild the exact fingerprint/input contract emitted by run_netbid_qc.py."""
    if mode not in {"summary", "html"}:
        raise ValueError(f"Invalid NetBID2 mode requested by packager: {mode}")
    filename, prefix = NETBID_DRIVERS[driver]
    arm_root = work_root / "results" / key / driver
    point_manifest_path = work_root / "results" / key / "point_manifest.json"
    point_manifest = load_json(point_manifest_path)
    consensus_path = arm_root / "consensus" / "consensus_network_ncol_.txt"
    consensus_manifest_path = arm_root / "consensus_manifest.json"
    driver_path = work_root / "inputs" / filename
    suffix = "netbid2_qc" if mode == "summary" else "netbid2_qc_html"
    output_root = arm_root / suffix
    partial_root = arm_root / f"{suffix}.partial"
    arm_analysis = analysis["arms"][f"{key}/{driver}"]

    fingerprint_payload = {
        "schema": "sjaracne-brca100-pr67-p-sweep-netbid2-v1",
        "mode": mode,
        "point": key,
        "p_value": point_manifest.get("p_value"),
        "mi_cutoff": point_manifest.get("mi_cutoff"),
        "point_manifest_sha256": sha256_file(point_manifest_path),
        "sweep_design_sha256": analysis["sweep_design_sha256"],
        "driver": driver,
        "driver_sha256": sha256_file(require_file(driver_path)),
        "consensus_sha256": sha256_file(require_file(consensus_path)),
        "consensus_manifest_sha256": sha256_file(consensus_manifest_path),
        "r_script_sha256": r_script_sha256,
        "wrapper_sha256": wrapper_sha256,
        "environment": environment,
        "prefix": prefix,
    }
    required_fields = set(fingerprint_payload) | {
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
        raise ValueError(
            f"NetBID2 {mode} manifest has missing/arbitrary fields at "
            f"{manifest_path}: missing={sorted(required_fields - record_fields)}, "
            f"extra={sorted(record_fields - allowed_fields)}"
        )
    if (
        "recovered_after_interrupted_manifest" in record
        and record["recovered_after_interrupted_manifest"] is not True
    ):
        raise ValueError(f"Invalid NetBID2 recovery marker in {manifest_path}")
    mismatches = [
        field
        for field, expected in fingerprint_payload.items()
        if record.get(field) != expected
    ]
    if mismatches:
        raise ValueError(
            f"NetBID2 {mode} fingerprint inputs changed in {manifest_path}: "
            + ", ".join(mismatches)
        )
    expected_fingerprint = json_fingerprint(fingerprint_payload)
    if record.get("fingerprint") != expected_fingerprint:
        raise ValueError(f"NetBID2 {mode} fingerprint mismatch in {manifest_path}")
    expected_command = [
        str(wrapper),
        "Rscript",
        str(r_script),
        str(consensus_path),
        str(driver_path),
        str(partial_root),
        prefix,
        "false" if mode == "summary" else "true",
    ]
    if record.get("command") != expected_command:
        raise ValueError(f"NetBID2 {mode} command mismatch in {manifest_path}")
    if not isinstance(record.get("finished_at_utc"), str) or not record["finished_at_utc"]:
        raise ValueError(f"Missing NetBID2 completion timestamp in {manifest_path}")
    stderr_bytes = record.get("stderr_bytes")
    if isinstance(stderr_bytes, bool) or not isinstance(stderr_bytes, int) or stderr_bytes < 0:
        raise ValueError(f"Invalid NetBID2 stderr byte count in {manifest_path}")
    same_path(record.get("output"), output_root, f"NetBID2 {mode} output", manifest_path)
    if (
        fingerprint_payload["point_manifest_sha256"]
        != arm_analysis["point_manifest_sha256"]
        or fingerprint_payload["consensus_sha256"]
        != arm_analysis["consensus_sha256"]
        or fingerprint_payload["consensus_manifest_sha256"]
        != arm_analysis["consensus_manifest_sha256"]
    ):
        raise ValueError(f"NetBID2 {mode} inputs disagree with analysis: {key}/{driver}")


def validate_netbid(work_root: Path, analysis: dict[str, Any], evidence: PackageEvidence) -> None:
    results = work_root / "results"
    benchmark_root = Path(__file__).resolve().parents[1]
    wrapper = require_file(benchmark_root / "brca100_netbid_qc" / "netbid2-r")
    r_script = require_file(Path(__file__).resolve().with_name("run_netbid_qc.R"))
    wrapper_sha256 = sha256_file(wrapper)
    r_script_sha256 = sha256_file(r_script)
    aggregate_path = results / "netbid2_qc_manifest.json"
    aggregate = load_json(aggregate_path)
    aggregate_hash = sha256_file(aggregate_path)
    if aggregate_hash != analysis.get("netbid2_aggregate_manifest_sha256"):
        raise ValueError("NetBID2 summary aggregate changed after analysis")
    fingerprint = aggregate.get("fingerprint")
    payload = dict(aggregate)
    payload.pop("fingerprint", None)
    expected_order = [(point, driver) for point in POINT_KEYS for driver in DRIVERS]
    records = aggregate.get("summary_runs")
    summary_selection = aggregate.get("selection")
    environment = validate_netbid_environment(aggregate.get("environment"), aggregate_path)
    if (
        set(aggregate)
        != {
            "schema", "environment", "sweep_design_sha256", "all_sweep_points",
            "selection", "summary_runs", "fingerprint",
        }
        or not isinstance(summary_selection, dict)
        or set(summary_selection) != {"points", "drivers"}
        or aggregate.get("schema")
        != "sjaracne-brca100-pr67-p-sweep-netbid2-summary-aggregate-v1"
        or aggregate.get("sweep_design_sha256") != analysis["sweep_design_sha256"]
        or aggregate.get("all_sweep_points") != list(POINT_KEYS)
        or aggregate.get("selection")
        != {"points": list(POINT_KEYS), "drivers": list(DRIVERS)}
        or not isinstance(records, list)
        or any(not isinstance(record, dict) for record in records)
        or [(record.get("point"), record.get("driver")) for record in records]
        != expected_order
        or fingerprint != json_fingerprint(payload)
    ):
        raise ValueError("NetBID2 summary aggregate is incomplete or has a bad fingerprint")
    evidence.copy(
        aggregate_path,
        "provenance/aggregates/netbid2_qc_manifest.json",
        "netbid2-summary-aggregate-manifest",
    )
    summary_manifests: dict[tuple[str, str], dict[str, Any]] = {}
    summary_roots: dict[tuple[str, str], Path] = {}
    for (key, driver), aggregate_record in zip(expected_order, records):
        arm_root = results / key / driver
        manifest_path = arm_root / "netbid2_qc_manifest.json"
        manifest = load_json(manifest_path)
        arm_analysis = analysis["arms"][f"{key}/{driver}"]
        validate_netbid_run_contract(
            record=manifest,
            manifest_path=manifest_path,
            mode="summary",
            key=key,
            driver=driver,
            work_root=work_root,
            analysis=analysis,
            environment=environment,
            wrapper=wrapper,
            wrapper_sha256=wrapper_sha256,
            r_script=r_script,
            r_script_sha256=r_script_sha256,
        )
        if (
            manifest != aggregate_record
            or sha256_file(manifest_path) != arm_analysis["netbid2_manifest_sha256"]
            or manifest.get("schema") != "sjaracne-brca100-pr67-p-sweep-netbid2-v1"
            or manifest.get("mode") != "summary"
            or manifest.get("point") != key
            or manifest.get("driver") != driver
            or manifest.get("fingerprint") != arm_analysis["netbid2_fingerprint"]
        ):
            raise ValueError(f"NetBID2 summary manifest changed: {key}/{driver}")
        output_root = arm_root / "netbid2_qc"
        same_path(manifest.get("output"), output_root, "NetBID2 output", manifest_path)
        inventory_names, inventory_digest = validate_inventory(
            output_root, manifest.get("output_inventory"), evidence,
            "netbid2-summary-output", manifest_path,
        )
        if inventory_names != list(NETBID_SHARED_OUTPUTS):
            raise ValueError(f"Unexpected NetBID2 summary file set: {key}/{driver}")
        if (
            inventory_digest != arm_analysis["netbid2_qc_set_sha256"]
            or pkg_network_hash(manifest, "network_summary.tsv", manifest_path)
            != arm_analysis["netbid2_network_summary_sha256"]
        ):
            raise ValueError(f"NetBID2 output set changed after analysis: {key}/{driver}")
        for stream in ("stdout", "stderr"):
            log = arm_root / "logs" / f"netbid2_qc.{stream}.log"
            checksum = require_sha256(manifest.get(f"{stream}_sha256"), stream, manifest_path)
            if sha256_file(require_file(log)) != checksum:
                raise ValueError(f"NetBID2 log changed: {log}")
            evidence.omit_file(
                log,
                "netbid2-summary-log",
                checksum,
                verification="rehashed-by-packager",
            )
        if int(manifest.get("stderr_bytes", -1)) != (
            arm_root / "logs" / "netbid2_qc.stderr.log"
        ).stat().st_size:
            raise ValueError(f"NetBID2 stderr size changed: {key}/{driver}")
        for pending in (
            arm_root / "netbid2_qc.partial",
            arm_root / "netbid2_qc_manifest.pending.json",
        ):
            if pending.exists():
                raise ValueError(f"Unresolved NetBID2 summary state: {pending}")
        evidence.copy(
            manifest_path,
            f"provenance/arms/{key}/{driver}/netbid2_qc_manifest.json",
            "netbid2-summary-manifest",
        )
        summary_manifests[(key, driver)] = manifest
        summary_roots[(key, driver)] = output_root

    html_aggregate_path = results / "netbid2_qc_html_manifest.json"
    per_arm_html = [
        results / key / driver / "netbid2_qc_html_manifest.json"
        for key in POINT_KEYS
        for driver in DRIVERS
    ]
    if not html_aggregate_path.exists():
        unexpected: list[Path] = []
        for manifest_path in per_arm_html:
            arm_root = manifest_path.parent
            unexpected.extend(
                path
                for path in (
                    manifest_path,
                    arm_root / "netbid2_qc_html",
                    arm_root / "netbid2_qc_html.partial",
                    arm_root / "netbid2_qc_html_manifest.pending.json",
                )
                if path.exists()
            )
        if unexpected:
            raise ValueError(
                "Per-arm HTML state exists without the root aggregate: "
                f"{unexpected[:3]}"
            )
        return
    html_aggregate = load_json(html_aggregate_path)
    html_payload = dict(html_aggregate)
    html_fingerprint = html_payload.pop("fingerprint", None)
    selection = html_aggregate.get("selection")
    html_points = selection.get("html_points") if isinstance(selection, dict) else None
    html_records = html_aggregate.get("html_runs")
    html_environment = validate_netbid_environment(
        html_aggregate.get("environment"), html_aggregate_path
    )
    if (
        set(html_aggregate)
        != {
            "schema", "environment", "sweep_design_sha256", "all_sweep_points",
            "selection", "html_runs", "fingerprint",
        }
        or not isinstance(selection, dict)
        or set(selection) != {"points", "drivers", "html_points"}
        or html_aggregate.get("schema")
        != "sjaracne-brca100-pr67-p-sweep-netbid2-html-aggregate-v1"
        or html_aggregate.get("sweep_design_sha256") != analysis["sweep_design_sha256"]
        or html_aggregate.get("all_sweep_points") != list(POINT_KEYS)
        or html_environment != environment
        or not isinstance(html_points, list)
        or not html_points
        or len(html_points) != len(set(html_points))
        or any(point not in POINT_KEYS for point in html_points)
        or selection.get("points") != list(POINT_KEYS)
        or selection.get("drivers") != list(DRIVERS)
        or not isinstance(html_records, list)
        or any(not isinstance(record, dict) for record in html_records)
        or html_fingerprint != json_fingerprint(html_payload)
    ):
        raise ValueError("Optional NetBID2 HTML aggregate is incomplete or invalid")
    expected_html_order = [(point, driver) for point in html_points for driver in DRIVERS]
    if [
        (record.get("point"), record.get("driver")) for record in html_records
    ] != expected_html_order:
        raise ValueError("Optional NetBID2 HTML arm ordering/coverage is invalid")
    evidence.copy(
        html_aggregate_path,
        "provenance/aggregates/netbid2_qc_html_manifest.json",
        "netbid2-html-aggregate-manifest",
    )
    record_map = dict(zip(expected_html_order, html_records))
    for key, driver in [(point, driver) for point in POINT_KEYS for driver in DRIVERS]:
        arm_root = results / key / driver
        manifest_path = arm_root / "netbid2_qc_html_manifest.json"
        output_root = arm_root / "netbid2_qc_html"
        pending_paths = (
            arm_root / "netbid2_qc_html.partial",
            arm_root / "netbid2_qc_html_manifest.pending.json",
        )
        if key not in html_points:
            if (
                manifest_path.exists()
                or output_root.exists()
                or any(path.exists() for path in pending_paths)
            ):
                raise ValueError(f"Unaggregated optional HTML state: {key}/{driver}")
            continue
        if any(path.exists() for path in pending_paths):
            raise ValueError(f"Unresolved optional HTML state: {key}/{driver}")
        manifest = load_json(manifest_path)
        validate_netbid_run_contract(
            record=manifest,
            manifest_path=manifest_path,
            mode="html",
            key=key,
            driver=driver,
            work_root=work_root,
            analysis=analysis,
            environment=environment,
            wrapper=wrapper,
            wrapper_sha256=wrapper_sha256,
            r_script=r_script,
            r_script_sha256=r_script_sha256,
        )
        if (
            manifest != record_map[(key, driver)]
            or manifest.get("schema") != "sjaracne-brca100-pr67-p-sweep-netbid2-v1"
            or manifest.get("mode") != "html"
            or manifest.get("point") != key
            or manifest.get("driver") != driver
        ):
            raise ValueError(f"Optional HTML manifest mismatch: {key}/{driver}")
        same_path(manifest.get("output"), output_root, "NetBID2 HTML output", manifest_path)
        inventory_names, _inventory_digest = validate_inventory(
            output_root, manifest.get("output_inventory"), evidence,
            "netbid2-html-output", manifest_path,
        )
        prefix = NETBID_DRIVERS[driver][1]
        expected_inventory_names = sorted(
            [
                *NETBID_SHARED_OUTPUTS,
                f"{prefix}netQC.Rmd",
                f"{prefix}netQC.html",
            ]
        )
        if inventory_names != expected_inventory_names:
            raise ValueError(f"Unexpected HTML output inventory: {key}/{driver}")
        summary_manifest = summary_manifests[(key, driver)]
        summary_inventory = {
            item["path"]: item for item in summary_manifest["output_inventory"]
        }
        html_inventory = {item["path"]: item for item in manifest["output_inventory"]}
        for shared_name in NETBID_SHARED_OUTPUTS:
            summary_path = summary_roots[(key, driver)] / shared_name
            html_path = output_root / shared_name
            if (
                html_inventory.get(shared_name) != summary_inventory.get(shared_name)
                or not files_equal(summary_path, html_path)
            ):
                raise ValueError(
                    f"HTML/stable-summary shared TSV mismatch: "
                    f"{key}/{driver}/{shared_name}"
                )
        for stream in ("stdout", "stderr"):
            log = arm_root / "logs" / f"netbid2_qc_html.{stream}.log"
            checksum = require_sha256(manifest.get(f"{stream}_sha256"), stream, manifest_path)
            if sha256_file(require_file(log)) != checksum:
                raise ValueError(f"Optional HTML log changed: {log}")
            evidence.omit_file(
                log,
                "netbid2-html-log",
                checksum,
                verification="rehashed-by-packager",
            )
        if int(manifest.get("stderr_bytes", -1)) != (
            arm_root / "logs" / "netbid2_qc_html.stderr.log"
        ).stat().st_size:
            raise ValueError(f"Optional HTML stderr size changed: {key}/{driver}")
        evidence.copy(
            manifest_path,
            f"provenance/arms/{key}/{driver}/netbid2_qc_html_manifest.json",
            "netbid2-html-manifest",
        )


def add_analysis_copy_plan(
    analysis_root: Path, analysis: dict[str, Any], evidence: PackageEvidence
) -> None:
    for relative in sorted(analysis["output_files"]):
        evidence.copy(analysis_root / relative, f"analysis/{relative}", "analysis-output")
    evidence.copy(
        analysis_root / "analysis_manifest.json",
        "analysis/analysis_manifest.json",
        "analysis-manifest",
    )


def validate_optional_netbid_migration_history(
    work_root: Path, evidence: PackageEvidence
) -> None:
    """Validate and package any audited v1-to-v2 NetBID manifest refresh."""
    history_root = work_root / "netbid2_manifest_history"
    if not history_root.exists():
        return
    if not history_root.is_dir():
        raise ValueError(f"NetBID2 manifest history is not a directory: {history_root}")
    migration_paths = sorted(history_root.glob("*_to_*/migration.json"))
    if len(migration_paths) != 1:
        raise ValueError("Expected exactly one NetBID2 manifest migration audit")

    import run_netbid_qc as netbid_qc

    sweep_path = work_root / "sweep_design.json"
    sweep_design = load_json(sweep_path)
    sweep_hash = sha256_file(sweep_path)
    design_migration = netbid_qc.load_sweep_design_migration(
        work_root=work_root,
        sweep_design=sweep_design,
        sweep_design_hash=sweep_hash,
    )
    if design_migration is None:
        raise ValueError("NetBID2 history exists without a sweep-design migration")
    migration_path = migration_paths[0]
    runs = netbid_qc.validate_netbid_migration_audit(
        work_root=work_root,
        migration=design_migration,
        audit_path=migration_path,
        audit=load_json(migration_path),
    )
    summary_identities = {
        (str(item["point"]), str(item["driver"]))
        for item in runs
        if item.get("mode") == "summary"
    }
    expected_summaries = {
        (point, driver) for point in LEGACY_POINT_KEYS for driver in DRIVERS
    }
    if summary_identities != expected_summaries:
        raise ValueError(
            "NetBID2 migration audit does not cover every legacy summary arm"
        )
    evidence.copy(
        migration_path,
        "provenance/netbid2_manifest_history/migration.json",
        "netbid2-manifest-migration",
    )
    for item in runs:
        archive = work_root / str(item["archived_manifest_path"])
        evidence.copy(
            archive,
            "provenance/netbid2_manifest_history/"
            f"arms/{item['point']}/{item['driver']}/{item['mode']}_manifest.json",
            "archived-netbid2-manifest",
        )


def assert_no_forbidden_payloads(root: Path) -> None:
    forbidden_names = {
        "consensus_network_ncol_.txt",
        "consensus_network_3col_.txt",
        "consensus_support.tsv",
    }
    violations: list[str] = []
    for path in root.rglob("*"):
        if path.is_symlink():
            violations.append(path.relative_to(root).as_posix())
        elif path.is_file() and (
            path.suffix == ".adj"
            or path.name in forbidden_names
            or path.name.endswith("netQC.html")
            or path.name.endswith("netQC.Rmd")
        ):
            violations.append(path.relative_to(root).as_posix())
    if violations:
        raise ValueError(f"Forbidden full artifacts entered compact package: {violations}")


def write_sha256s(root: Path) -> None:
    checksum_path = root / "SHA256SUMS"
    if checksum_path.exists():
        checksum_path.unlink()
    files = sorted(
        (path for path in root.rglob("*") if path.is_file()),
        key=lambda path: path.relative_to(root).as_posix(),
    )
    lines = [
        f"{sha256_file(path)}  {path.relative_to(root).as_posix()}\n" for path in files
    ]
    checksum_path.write_bytes("".join(lines).encode("utf-8"))
    with checksum_path.open("r", encoding="utf-8", newline="") as handle:
        parsed = [line.rstrip("\n").split("  ", 1) for line in handle]
    if [relative for _digest, relative in parsed] != [
        path.relative_to(root).as_posix() for path in files
    ]:
        raise RuntimeError("SHA256SUMS ordering validation failed")
    for (digest, relative), path in zip(parsed, files):
        if digest != sha256_file(path) or relative != path.relative_to(root).as_posix():
            raise RuntimeError(f"SHA256SUMS self-validation failed: {relative}")


def package_results(work_root: Path, output_root: Path) -> Path:
    work_root = work_root.resolve(strict=True)
    require_directory(work_root)
    output_root = output_root.resolve(strict=False)
    if output_root.exists():
        raise ValueError(f"Output root must not exist: {output_root}")
    try:
        output_root.relative_to(work_root)
    except ValueError:
        pass
    else:
        raise ValueError("Output root must be outside the live sweep work root")
    try:
        work_root.relative_to(output_root)
    except ValueError:
        pass
    else:
        raise ValueError("Output root cannot contain the live sweep work root")
    staging = output_root.with_name(output_root.name + ".partial")
    if staging.exists():
        raise ValueError(f"Staging root already exists: {staging}")

    analysis, analysis_root = validate_analysis(work_root)
    evidence = PackageEvidence(work_root)
    _design, point_records = validate_design_and_build(work_root, analysis, evidence)
    validate_invocations(work_root, evidence)
    by_arm = validate_run_manifest(work_root, analysis, point_records, evidence)
    validate_anchor(work_root, analysis, evidence)
    validate_consensus_support(work_root, analysis, by_arm, evidence)
    validate_netbid(work_root, analysis, evidence)
    validate_optional_netbid_migration_history(work_root, evidence)
    add_analysis_copy_plan(analysis_root, analysis, evidence)

    output_root.parent.mkdir(parents=True, exist_ok=True)
    staging.mkdir(parents=False, exist_ok=False)
    try:
        copied: list[dict[str, Any]] = []
        for destination in sorted(evidence.copy_plan):
            source, role = evidence.copy_plan[destination]
            target = staging / destination
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
            source_hash = sha256_file(source)
            target_hash = sha256_file(target)
            if target_hash != source_hash or target.stat().st_size != source.stat().st_size:
                raise RuntimeError(f"Byte-preserving copy failed: {source} -> {target}")
            copied.append(
                {
                    "bytes": target.stat().st_size,
                    "package_path": destination,
                    "role": role,
                    "sha256": target_hash,
                    "source_work_root_relative": work_relative(source, work_root),
                }
            )

        attributes = staging / ".gitattributes"
        attributes.write_bytes(b"* -text -whitespace\n")
        omitted_payload: dict[str, Any] = {
            "schema": "sjaracne-brca100-pr67-threshold-sweep-omitted-artifacts-v1",
            "sweep_design_sha256": analysis["sweep_design_sha256"],
            "analysis_manifest_sha256": sha256_file(
                analysis_root / "analysis_manifest.json"
            ),
            "scope": (
                "Manifest-addressed inputs, build products, seed outputs/logs, "
                "consensus/support products, and NetBID2 products intentionally "
                "excluded from the compact package. Extracted source trees and "
                "transient work directories are outside this inventory."
            ),
            "files": [
                evidence.omitted_files[key] for key in sorted(evidence.omitted_files)
            ],
            "file_sets": [
                evidence.omitted_sets[key] for key in sorted(evidence.omitted_sets)
            ],
        }
        omitted_payload["file_count"] = len(omitted_payload["files"])
        omitted_payload["file_set_count"] = len(omitted_payload["file_sets"])
        omitted_payload["category_counts"] = dict(
            sorted(Counter(item["category"] for item in omitted_payload["files"]).items())
        )
        omitted_payload["fingerprint"] = json_fingerprint(omitted_payload)
        omitted_path = staging / "omitted_artifacts.json"
        atomic_json(omitted_path, omitted_payload)

        package_manifest: dict[str, Any] = {
            "schema": "sjaracne-brca100-pr67-threshold-sweep-compact-package-v1",
            "source_work_root": str(work_root),
            "sweep_design_sha256": analysis["sweep_design_sha256"],
            "analysis_manifest_sha256": sha256_file(
                analysis_root / "analysis_manifest.json"
            ),
            "copied_files": copied,
            "generated_files": {
                ".gitattributes": {
                    "bytes": attributes.stat().st_size,
                    "sha256": sha256_file(attributes),
                },
                "omitted_artifacts.json": {
                    "bytes": omitted_path.stat().st_size,
                    "sha256": sha256_file(omitted_path),
                },
            },
            "omitted_file_count": len(omitted_payload["files"]),
            "omitted_file_set_count": len(omitted_payload["file_sets"]),
            "optional_html_provenance_included": (
                work_root / "results" / "netbid2_qc_html_manifest.json"
            ).is_file(),
            "packager_sha256": sha256_file(Path(__file__).resolve()),
            "large_artifacts_copied": False,
            "checksum_policy": "SHA256SUMS covers every package file except itself",
        }
        package_manifest["fingerprint"] = json_fingerprint(package_manifest)
        atomic_json(staging / "package_manifest.json", package_manifest)
        assert_no_forbidden_payloads(staging)
        write_sha256s(staging)
        expected_checksum_entries = (
            sum(1 for path in staging.rglob("*") if path.is_file()) - 1
        )
        actual_checksum_entries = len(
            (staging / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
        )
        if actual_checksum_entries != expected_checksum_entries:
            raise RuntimeError("SHA256SUMS does not cover every non-self package file")
        os.replace(staging, output_root)
    except BaseException:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    return output_root


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output = package_results(args.work_root, args.output_root)
    print(f"Wrote compact threshold-sweep package to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
