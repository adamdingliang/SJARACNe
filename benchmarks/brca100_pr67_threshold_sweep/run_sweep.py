#!/usr/bin/env python3
"""Run a matched BRCA100 sweep of only PR67's per-subsample p-value.

The fixed design is PR67 commit 7633ebb, m=80 without replacement, Npar=40,
DPI epsilon 0, seeds 1..100, and consensus p=1e-5.  Only the seed-level
AP-MI p-value changes.  Every seed is inferred directly; adjacency replay is
not used because its rounded serialized MI values are not an exact substitute
for the in-memory threshold-plus-DPI path.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed


THIS_DIR = Path(__file__).resolve().parent
NETBID_BENCHMARK_DIR = THIS_DIR.parent / "brca100_netbid_qc"
sys.path.insert(0, str(NETBID_BENCHMARK_DIR))
import run_workflows as core  # noqa: E402


PR67_COMMIT = "7633ebb4a0d966dbda15a4e32d0efa492fb71aeb"
PR67_BUILD_STAGE = core.Stage(
    key="pr67_7633ebb",
    commit=PR67_COMMIT,
    sampling_args=("-u", "80%"),
    required_headers=(),
    null_model_relative=(
        "SJARACNe/config/apmi_null/apmi_null_m00080_npar040.model"
    ),
)


@dataclasses.dataclass(frozen=True)
class SweepPoint:
    key: str
    p_token: str
    label: str
    role: str

    @property
    def p_value(self) -> float:
        return float(self.p_token)


# Original coarse-first grid.  Its final point is the exact p at which the
# PR67 GPD cutoff equals PR66's legacy m=80 affine cutoff
# (0.17280321515749669).  Keep this tuple immutable: it is the exact v1 prefix
# accepted by the append-only design migration.
LEGACY_SWEEP_POINTS = (
    SweepPoint("p1e-07", "1e-7", "1e-7", "original-default"),
    SweepPoint("p1e-06", "1e-6", "1e-6", "extrapolated-grid"),
    SweepPoint("p1e-05", "1e-5", "1e-5", "extrapolated-grid"),
    SweepPoint("p2e-05", "2e-5", "2e-5", "direct-validation-boundary"),
    SweepPoint("p5e-05", "5e-5", "5e-5", "validated-grid"),
    SweepPoint("p1e-04", "1e-4", "1e-4", "validated-grid"),
    SweepPoint("p2e-04", "2e-4", "2e-4", "validated-grid"),
    SweepPoint("p3e-04", "3e-4", "3e-4", "validated-grid"),
    SweepPoint(
        "p_pr66_cutoff_match",
        "0.000352804562601613",
        "3.528045626e-4",
        "pr66-cutoff-match",
    ),
)
EXTENDED_SWEEP_POINTS = (
    SweepPoint("p4e-04", "4e-4", "4e-4", "extended-validated-grid"),
    SweepPoint("p5e-04", "5e-4", "5e-4", "extended-validated-grid"),
    SweepPoint("p7p5e-04", "7.5e-4", "7.5e-4", "extended-validated-grid"),
    SweepPoint("p1e-03", "1e-3", "1e-3", "extended-validated-grid"),
)
SWEEP_POINTS = LEGACY_SWEEP_POINTS + EXTENDED_SWEEP_POINTS


DRIVERS = core.DRIVER_CLASSES
MODEL_EXPECTED_SHA256 = (
    "e3a8522682a8ea239821aaa10b12db72d00e07bfdcad43599d8e76a06be80944"
)
PR66_MATCHED_CUTOFF = 0.17280321515749669
LEGACY_SWEEP_DESIGN_SCHEMA = "sjaracne-brca100-pr67-p-sweep-v1"
SWEEP_DESIGN_SCHEMA = "sjaracne-brca100-pr67-p-sweep-v2"
SWEEP_DESIGN_MIGRATION_SCHEMA = (
    "sjaracne-brca100-pr67-p-sweep-design-migration-v1"
)
SWEEP_DESIGN_HISTORY_DIRECTORY = "sweep_design_history"


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def serialized_json(payload: object) -> bytes:
    return (
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def atomic_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("wb") as handle:
        handle.write(payload)
    os.replace(temporary, path)


def atomic_json(path: Path, payload: object) -> None:
    atomic_bytes(path, serialized_json(payload))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def json_fingerprint(payload: object) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def parse_model(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise ValueError(f"Malformed model line {path}:{line_number}")
        key, value = line.split("=", 1)
        if not key or key in result:
            raise ValueError(f"Invalid/duplicate model key {key!r} in {path}")
        result[key] = value
    return result


def gpd_cutoff(model: dict[str, str], probability: float) -> float:
    threshold = float(model["tail_threshold"])
    tail_probability = float(model["tail_probability"])
    shape = float(model["tail_shape"])
    scale = float(model["tail_scale"])
    if not 0.0 < probability <= tail_probability:
        raise ValueError(
            f"p={probability:.17g} is outside the fitted tail (0,{tail_probability}]"
        )
    if abs(shape) < 1e-12:
        cutoff = threshold - scale * math.log(probability / tail_probability)
    else:
        cutoff = threshold + (scale / shape) * (
            (probability / tail_probability) ** (-shape) - 1.0
        )
    endpoint_text = model.get("tail_endpoint", "none")
    if endpoint_text != "none" and cutoff > float(endpoint_text) + 1e-12:
        raise ValueError(f"Cutoff {cutoff} exceeds fitted GPD endpoint")
    return cutoff


def point_design(point: SweepPoint, model: dict[str, str]) -> dict[str, object]:
    p_value = point.p_value
    cutoff = gpd_cutoff(model, p_value)
    validated_min = float(model["validated_p_min"])
    validated_max = float(model["validated_p_max"])
    if validated_min <= p_value <= validated_max:
        # This records the model's accepted continuous runtime range.  Reports
        # separately distinguish exact Gate-2 grid points from interpolation.
        validation_class = "independent-held-out-validated"
    elif p_value < validated_min:
        validation_class = "gpd-tail-extrapolation-below-held-out-range"
    else:
        validation_class = "modeled-but-above-held-out-validation-range"
    return {
        "key": point.key,
        "label": point.label,
        "role": point.role,
        "p_token": point.p_token,
        "p_value": p_value,
        "mi_cutoff": cutoff,
        "mi_cutoff_header": format(cutoff, ".6g"),
        "p_header": format(p_value, ".6g"),
        "tail_extrapolated": not (validated_min <= p_value <= validated_max),
        "validation_class": validation_class,
        "validated_p_min": validated_min,
        "validated_p_max": validated_max,
    }


def validation_stage(design: dict[str, object]) -> core.Stage:
    extrapolated = "yes" if design["tail_extrapolated"] else "no"
    return core.Stage(
        key=str(design["key"]),
        commit=PR67_COMMIT,
        sampling_args=("-u", "80%"),
        required_headers=(
            ">  MI threshold method estimator-matched AP-MI permutation-null GPD tail",
            ">  AP-MI null model m 80",
            ">  AP-MI null model Npar 40",
            f">  AP-MI cutoff tail extrapolated {extrapolated}",
            ">  Sampling method fixed-size without replacement",
            ">  Sampling request 80%",
            ">  Eligible observations 100",
            ">  Sampled observations 80",
            f">  MI threshold    {design['mi_cutoff_header']}",
            f">  MI P-value      {design['p_header']}",
        ),
        null_model_relative=PR67_BUILD_STAGE.null_model_relative,
    )


def command_for_job(
    *,
    design: dict[str, object],
    build: dict,
    driver: core.DriverClass,
    input_root: Path,
    seed: int,
    output: Path,
) -> list[str]:
    return [
        build["binary"],
        "-i",
        str(input_root / "BRCA100.exp"),
        "-l",
        str(input_root / driver.filename),
        "-s",
        str(input_root / driver.filename),
        "-p",
        str(design["p_token"]),
        "-e",
        "0",
        "-a",
        "adaptive_partitioning",
        "-H",
        build["config_directory"].rstrip("/") + "/",
        "-N",
        "40",
        "-S",
        str(seed),
        "-o",
        str(output),
        "-u",
        "80%",
        "-M",
        build["null_model"],
    ]


def point_manifest_payload(
    *, design: dict[str, object], build: dict, input_metadata: dict
) -> dict[str, object]:
    return {
        "schema": "sjaracne-brca100-pr67-p-sweep-point-v1",
        **design,
        "commit": PR67_COMMIT,
        "binary_sha256": build["binary_sha256"],
        "config_sha256": build["config_sha256"],
        "null_model_sha256": build["null_model_sha256"],
        "sampling": "fixed 80% without replacement",
        "m": 80,
        "npar": 40,
        "dpi_epsilon": 0,
        "consensus_p": 1e-5,
        "seeds": list(range(1, 101)),
        "inputs": input_metadata,
    }


def write_point_manifest(
    *, work_root: Path, design: dict[str, object], build: dict, input_metadata: dict
) -> None:
    manifest = point_manifest_payload(
        design=design, build=build, input_metadata=input_metadata
    )
    path = work_root / "results" / str(design["key"]) / "point_manifest.json"
    if path.is_file():
        existing = load_json(path)
        if existing != manifest:
            raise RuntimeError(f"Incompatible existing point manifest: {path}")
        return
    atomic_json(path, manifest)


def sweep_design_payload(
    *,
    schema: str,
    points: tuple[SweepPoint, ...],
    model: dict[str, str],
    build: dict,
    input_metadata: dict,
) -> dict[str, object]:
    return {
        "schema": schema,
        "commit": PR67_COMMIT,
        "binary_sha256": build["binary_sha256"],
        "config_sha256": build["config_sha256"],
        "null_model_sha256": build["null_model_sha256"],
        "all_points": [point_design(point, model) for point in points],
        "fixed_parameters": {
            "sampling": "fixed 80% without replacement",
            "m": 80,
            "npar": 40,
            "dpi_epsilon": 0,
            "consensus_p": 1e-5,
            "seeds": list(range(1, 101)),
        },
        "inputs": input_metadata,
    }


def validate_append_only_design_extension(
    legacy_design: dict[str, object], extended_design: dict[str, object]
) -> None:
    if legacy_design.get("schema") != LEGACY_SWEEP_DESIGN_SCHEMA:
        raise RuntimeError("Legacy sweep design does not have the exact v1 schema")
    if extended_design.get("schema") != SWEEP_DESIGN_SCHEMA:
        raise RuntimeError("Extended sweep design does not have the exact v2 schema")

    legacy_invariants = {
        key: value
        for key, value in legacy_design.items()
        if key not in {"schema", "all_points"}
    }
    extended_invariants = {
        key: value
        for key, value in extended_design.items()
        if key not in {"schema", "all_points"}
    }
    if legacy_invariants != extended_invariants:
        raise RuntimeError("Sweep extension changes a fixed v1 design invariant")

    legacy_points = legacy_design.get("all_points")
    extended_points = extended_design.get("all_points")
    if not isinstance(legacy_points, list) or not isinstance(extended_points, list):
        raise RuntimeError("Sweep design all_points must be a list")
    if extended_points[: len(legacy_points)] != legacy_points:
        raise RuntimeError("Sweep extension changes or reorders a legacy point")
    appended = extended_points[len(legacy_points) :]
    expected_appended_keys = [point.key for point in EXTENDED_SWEEP_POINTS]
    if [point.get("key") for point in appended] != expected_appended_keys:
        raise RuntimeError("Sweep extension is not the exact four-point append")
    if len(legacy_points) != len(LEGACY_SWEEP_POINTS):
        raise RuntimeError("Legacy sweep design is not the exact nine-point design")


def validate_legacy_point_manifests(
    *,
    work_root: Path,
    legacy_design: dict[str, object],
    build: dict,
    input_metadata: dict,
) -> None:
    expected = {
        str(point["key"]): point_manifest_payload(
            design=point, build=build, input_metadata=input_metadata
        )
        for point in legacy_design["all_points"]
    }
    results_root = work_root / "results"
    if not results_root.is_dir():
        return

    extension_keys = {point.key for point in EXTENDED_SWEEP_POINTS}
    for extension_key in sorted(extension_keys):
        extension_root = results_root / extension_key
        if extension_root.exists():
            raise RuntimeError(
                "A v1 sweep root already contains an extended-point path: "
                f"{extension_root}"
            )

    for manifest_path in sorted(results_root.glob("*/point_manifest.json")):
        point_key = manifest_path.parent.name
        if point_key not in expected:
            raise RuntimeError(
                f"Unexpected point manifest under the v1 sweep: {manifest_path}"
            )
        if load_json(manifest_path) != expected[point_key]:
            raise RuntimeError(
                f"Incompatible legacy point manifest: {manifest_path}"
            )
    for point_key in sorted(expected):
        point_root = results_root / point_key
        if point_root.is_dir() and not (point_root / "point_manifest.json").is_file():
            raise RuntimeError(
                f"Legacy point directory lacks its manifest: {point_root}"
            )


def ensure_exact_bytes(path: Path, payload: bytes, description: str) -> None:
    if path.exists():
        if not path.is_file() or path.read_bytes() != payload:
            raise RuntimeError(f"Incompatible existing {description}: {path}")
        return
    atomic_bytes(path, payload)


def migration_manifest_payload(
    *,
    work_root: Path,
    archive_path: Path,
    migration_path: Path,
    legacy_design: dict[str, object],
    extended_design: dict[str, object],
    legacy_sha256: str,
    extended_sha256: str,
) -> dict[str, object]:
    legacy_points = legacy_design["all_points"]
    extended_points = extended_design["all_points"]
    payload: dict[str, object] = {
        "schema": SWEEP_DESIGN_MIGRATION_SCHEMA,
        "operation": "append-only-point-extension",
        "from": {
            "schema": LEGACY_SWEEP_DESIGN_SCHEMA,
            "sha256": legacy_sha256,
            "archived_path": archive_path.relative_to(work_root).as_posix(),
            "point_keys": [point["key"] for point in legacy_points],
        },
        "to": {
            "schema": SWEEP_DESIGN_SCHEMA,
            "sha256": extended_sha256,
            "active_path": "sweep_design.json",
            "point_keys": [point["key"] for point in extended_points],
        },
        "appended_points": extended_points[len(legacy_points) :],
        "manifest_path": migration_path.relative_to(work_root).as_posix(),
    }
    payload["fingerprint"] = json_fingerprint(payload)
    return payload


def ensure_sweep_design(
    *,
    work_root: Path,
    legacy_design: dict[str, object],
    extended_design: dict[str, object],
    build: dict,
    input_metadata: dict,
) -> str:
    """Create v2 or atomically migrate an exact v1 design before point writes."""
    validate_append_only_design_extension(legacy_design, extended_design)
    design_path = work_root / "sweep_design.json"
    partial_path = design_path.with_name(design_path.name + ".partial")
    if partial_path.exists():
        raise RuntimeError(f"Unresolved partial sweep design: {partial_path}")

    if not design_path.exists():
        stray_manifests = sorted((work_root / "results").glob("*/point_manifest.json"))
        if stray_manifests:
            raise RuntimeError(
                "Point manifests exist without a sweep design; refusing to create v2: "
                f"{stray_manifests[0]}"
            )
        atomic_json(design_path, extended_design)
        return "created-v2"
    if not design_path.is_file():
        raise RuntimeError(f"Sweep design path is not a regular file: {design_path}")

    existing_bytes = design_path.read_bytes()
    existing = load_json(design_path)
    if existing == extended_design:
        if existing_bytes != serialized_json(extended_design):
            raise RuntimeError(
                f"Existing v2 sweep design is not canonical: {design_path}"
            )
        return "existing-v2"
    if existing != legacy_design:
        raise RuntimeError(f"Incompatible existing sweep design: {design_path}")
    if existing_bytes != serialized_json(legacy_design):
        raise RuntimeError(
            f"Existing v1 sweep design is not canonical: {design_path}"
        )

    validate_legacy_point_manifests(
        work_root=work_root,
        legacy_design=legacy_design,
        build=build,
        input_metadata=input_metadata,
    )
    legacy_sha256 = sha256_bytes(existing_bytes)
    extended_bytes = serialized_json(extended_design)
    extended_sha256 = sha256_bytes(extended_bytes)
    history_root = work_root / SWEEP_DESIGN_HISTORY_DIRECTORY
    archive_path = history_root / f"{legacy_sha256}.sweep_design.json"
    migration_path = (
        history_root / f"{legacy_sha256}_to_{extended_sha256}.migration.json"
    )
    manifest = migration_manifest_payload(
        work_root=work_root,
        archive_path=archive_path,
        migration_path=migration_path,
        legacy_design=legacy_design,
        extended_design=extended_design,
        legacy_sha256=legacy_sha256,
        extended_sha256=extended_sha256,
    )

    # Record immutable recovery/provenance evidence before promoting v2.  A
    # crash at either atomic write leaves v1 active and creates no new point
    # manifest; a retry validates/reuses these exact history records.
    ensure_exact_bytes(archive_path, existing_bytes, "archived v1 sweep design")
    ensure_exact_bytes(
        migration_path,
        serialized_json(manifest),
        "sweep-design migration manifest",
    )
    atomic_bytes(design_path, extended_bytes)
    if design_path.read_bytes() != extended_bytes:
        raise RuntimeError(f"Atomic sweep-design migration did not persist: {design_path}")
    return "migrated-v1-to-v2"


def run_seed_job(
    *,
    design: dict[str, object],
    build: dict,
    driver: core.DriverClass,
    input_root: Path,
    results_root: Path,
    seed: int,
    expression_ids: set[str],
    driver_ids: set[str],
) -> tuple[str, bool, dict]:
    point_key = str(design["key"])
    arm_root = results_root / point_key / driver.key
    adjacency_root = arm_root / "adjacency"
    log_root = arm_root / "logs"
    metadata_root = arm_root / "seed_metadata"
    adjacency_root.mkdir(parents=True, exist_ok=True)
    log_root.mkdir(parents=True, exist_ok=True)
    metadata_root.mkdir(parents=True, exist_ok=True)

    stem = f"TF_run_{seed:03d}"
    output = adjacency_root / f"{stem}.adj"
    marker = metadata_root / f"{stem}.json"
    preview = command_for_job(
        design=design,
        build=build,
        driver=driver,
        input_root=input_root,
        seed=seed,
        output=output,
    )
    fingerprint_payload = {
        "schema": "sjaracne-brca100-pr67-p-sweep-seed-v1",
        "point": design,
        "commit": PR67_COMMIT,
        "binary_sha256": build["binary_sha256"],
        "config_sha256": build["config_sha256"],
        "null_model_sha256": build["null_model_sha256"],
        "driver": driver.key,
        "driver_sha256": core.sha256_file(input_root / driver.filename),
        "expression_sha256": core.sha256_file(input_root / "BRCA100.exp"),
        "seed": seed,
        "command_without_output": [
            value if value != str(output) else "<OUTPUT>" for value in preview
        ],
    }
    fingerprint = json_fingerprint(fingerprint_payload)
    stage = validation_stage(design)

    if marker.is_file() and output.is_file():
        existing = load_json(marker)
        if existing.get("fingerprint") == fingerprint:
            stats = core.validate_adjacency(
                output,
                stage=stage,
                driver_ids=driver_ids,
                expression_ids=expression_ids,
            )
            if existing.get("adjacency", {}).get("full_sha256") == stats["full_sha256"]:
                return f"{point_key}/{driver.key}/{seed:03d}", True, existing
        raise RuntimeError(f"Stale or inconsistent completed seed {marker}")

    if marker.is_file() and not output.is_file():
        raise RuntimeError(f"Seed marker exists without adjacency: {marker}")

    # Recover the narrow crash window after atomic adjacency promotion but
    # before marker creation.  Every primary artifact is revalidated and the
    # recovered marker records that its wall-clock timestamps are unavailable.
    if output.is_file() and not marker.is_file():
        stdout_path = log_root / f"{stem}.stdout.log"
        stderr_path = log_root / f"{stem}.stderr.log"
        time_path = log_root / f"{stem}.time.txt"
        for required in (stdout_path, stderr_path, time_path):
            if not required.is_file():
                raise RuntimeError(
                    f"Unverifiable orphan adjacency {output}; missing {required}"
                )
        stats = core.validate_adjacency(
            output,
            stage=stage,
            driver_ids=driver_ids,
            expression_ids=expression_ids,
        )
        recovered = {
            "fingerprint": fingerprint,
            "schema": "sjaracne-brca100-pr67-p-sweep-seed-v1",
            "point": design,
            "commit": PR67_COMMIT,
            "driver": driver.key,
            "seed": seed,
            "command": preview,
            "binary_sha256": build["binary_sha256"],
            "config_sha256": build["config_sha256"],
            "null_model_sha256": build["null_model_sha256"],
            "recovered_orphan": True,
            "started_at_utc": None,
            "finished_at_utc": utc_now(),
            "high_resolution_wall_s": None,
            "gnu_time": core.parse_gnu_time(time_path),
            "stdout_sha256": core.sha256_file(stdout_path),
            "stderr_sha256": core.sha256_file(stderr_path),
            "stderr_bytes": stderr_path.stat().st_size,
            "adjacency": stats,
        }
        atomic_json(marker, recovered)
        return f"{point_key}/{driver.key}/{seed:03d}", True, recovered

    partial = arm_root / "work" / f"{stem}.adj.partial"
    partial.parent.mkdir(parents=True, exist_ok=True)
    if partial.exists():
        partial.unlink()
    stdout_path = log_root / f"{stem}.stdout.log"
    stderr_path = log_root / f"{stem}.stderr.log"
    time_path = log_root / f"{stem}.time.txt"
    command = command_for_job(
        design=design,
        build=build,
        driver=driver,
        input_root=input_root,
        seed=seed,
        output=partial,
    )
    timed_command = [
        "/usr/bin/time",
        "-f",
        "elapsed_s=%e\nuser_s=%U\nsystem_s=%S\nmax_rss_kib=%M",
        "-o",
        str(time_path),
        *command,
    ]
    started = utc_now()
    start_clock = time.perf_counter()
    try:
        core.checked_run(
            timed_command,
            cwd=arm_root,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
        )
    except Exception as error:
        raise RuntimeError(
            f"Seed failed ({point_key}/{driver.key}/{seed}); see {stderr_path}"
        ) from error
    high_resolution_wall = time.perf_counter() - start_clock
    stats = core.validate_adjacency(
        partial,
        stage=stage,
        driver_ids=driver_ids,
        expression_ids=expression_ids,
    )
    os.replace(partial, output)
    stats["full_sha256"] = core.sha256_file(output)
    record = {
        "fingerprint": fingerprint,
        "schema": "sjaracne-brca100-pr67-p-sweep-seed-v1",
        "point": design,
        "commit": PR67_COMMIT,
        "driver": driver.key,
        "seed": seed,
        "command": command_for_job(
            design=design,
            build=build,
            driver=driver,
            input_root=input_root,
            seed=seed,
            output=output,
        ),
        "binary_sha256": build["binary_sha256"],
        "config_sha256": build["config_sha256"],
        "null_model_sha256": build["null_model_sha256"],
        "started_at_utc": started,
        "finished_at_utc": utc_now(),
        "high_resolution_wall_s": high_resolution_wall,
        "gnu_time": core.parse_gnu_time(time_path),
        "stdout_sha256": core.sha256_file(stdout_path),
        "stderr_sha256": core.sha256_file(stderr_path),
        "stderr_bytes": stderr_path.stat().st_size,
        "adjacency": stats,
    }
    atomic_json(marker, record)
    return f"{point_key}/{driver.key}/{seed:03d}", False, record


def aggregate_seed_manifest(
    results_root: Path,
    designs: list[dict[str, object]],
    drivers: list[core.DriverClass],
) -> None:
    rows: list[dict[str, object]] = []
    for design in sorted(designs, key=lambda item: float(item["p_value"])):
        for driver in drivers:
            root = results_root / str(design["key"]) / driver.key / "seed_metadata"
            markers = sorted(root.glob("TF_run_*.json"))
            for marker in markers:
                record = load_json(marker)
                adjacency = record["adjacency"]
                timing = record["gnu_time"]
                rows.append(
                    {
                        "point": design["key"],
                        "p_value": design["p_value"],
                        "mi_cutoff": design["mi_cutoff"],
                        "validation_class": design["validation_class"],
                        "commit": PR67_COMMIT,
                        "driver": driver.key,
                        "seed": record["seed"],
                        "binary_sha256": record["binary_sha256"],
                        "elapsed_s": timing["elapsed_s"],
                        "user_s": timing["user_s"],
                        "system_s": timing["system_s"],
                        "max_rss_kib": timing["max_rss_kib"],
                        "edges": adjacency["edges"],
                        "source_rows": adjacency["source_rows"],
                        "adjacency_bytes": adjacency["bytes"],
                        "adjacency_sha256": adjacency["full_sha256"],
                        "data_sha256": adjacency["data_sha256"],
                        "stderr_bytes": record["stderr_bytes"],
                    }
                )
    output = results_root / "run_manifest.tsv"
    temporary = output.with_name(output.name + ".partial")
    fieldnames = list(rows[0]) if rows else []
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, output)


def select_points(specification: str) -> list[SweepPoint]:
    if specification == "all":
        return list(SWEEP_POINTS)
    requested = [item.strip() for item in specification.split(",") if item.strip()]
    known = {point.key: point for point in SWEEP_POINTS}
    unknown = sorted(set(requested) - set(known))
    if unknown:
        raise ValueError(f"Unknown sweep point(s): {', '.join(unknown)}")
    if not requested:
        raise ValueError("At least one sweep point must be requested")
    if len(set(requested)) != len(requested):
        raise ValueError("Duplicate sweep point requested")
    return [known[key] for key in requested]


def select_drivers(specification: str) -> list[core.DriverClass]:
    if specification == "all":
        return list(DRIVERS)
    requested = [item.strip() for item in specification.split(",") if item.strip()]
    known = {driver.key: driver for driver in DRIVERS}
    unknown = sorted(set(requested) - set(known))
    if unknown:
        raise ValueError(f"Unknown driver class(es): {', '.join(unknown)}")
    if not requested:
        raise ValueError("At least one driver class must be requested")
    if len(set(requested)) != len(requested):
        raise ValueError("Duplicate driver class requested")
    return [known[key] for key in requested]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("/mnt/d/GitHub/SJARACNe"))
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
    parser.add_argument(
        "--phase", choices=("prepare", "infer", "consensus", "all"), default="all"
    )
    parser.add_argument("--points", default="all")
    parser.add_argument("--drivers", default="all")
    parser.add_argument("--seed-start", type=int, default=1)
    parser.add_argument("--seed-end", type=int, default=100)
    parser.add_argument("--workers", type=int, default=12)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not (1 <= args.seed_start <= args.seed_end <= 100):
        raise ValueError("Seed range must be within 1..100")
    if args.workers < 1:
        raise ValueError("--workers must be positive")
    points = select_points(args.points)
    drivers = select_drivers(args.drivers)
    seeds = list(range(args.seed_start, args.seed_end + 1))
    args.work_root.mkdir(parents=True, exist_ok=True)

    input_metadata = core.stage_git_inputs(args.repo, args.work_root)
    expression_ids = core.parse_expression_ids(args.work_root / "inputs" / "BRCA100.exp")
    driver_ids = {
        driver.key: core.parse_driver_ids(args.work_root / "inputs" / driver.filename)
        for driver in drivers
    }
    build = core.extract_and_build_stage(args.repo, args.work_root, PR67_BUILD_STAGE)
    if build.get("commit") != PR67_COMMIT:
        raise RuntimeError("Unexpected PR67 build commit")
    if build.get("null_model_sha256") != MODEL_EXPECTED_SHA256:
        raise RuntimeError("Unexpected PR67 null-model SHA256")
    model = parse_model(Path(build["null_model"]))
    if model.get("m") != "80" or model.get("npar_limit") != "40":
        raise RuntimeError("The sweep requires the exact m=80, Npar=40 model")

    legacy_run_design = sweep_design_payload(
        schema=LEGACY_SWEEP_DESIGN_SCHEMA,
        points=LEGACY_SWEEP_POINTS,
        model=model,
        build=build,
        input_metadata=input_metadata,
    )
    run_design = sweep_design_payload(
        schema=SWEEP_DESIGN_SCHEMA,
        points=SWEEP_POINTS,
        model=model,
        build=build,
        input_metadata=input_metadata,
    )
    design_status = ensure_sweep_design(
        work_root=args.work_root,
        legacy_design=legacy_run_design,
        extended_design=run_design,
        build=build,
        input_metadata=input_metadata,
    )
    core.console(f"[DESIGN] {design_status}: {args.work_root / 'sweep_design.json'}")

    designs = [point_design(point, model) for point in points]
    matched = next(
        item for item in designs if item["key"] == "p_pr66_cutoff_match"
    ) if any(item["key"] == "p_pr66_cutoff_match" for item in designs) else None
    if matched is not None and not math.isclose(
        float(matched["mi_cutoff"]), PR66_MATCHED_CUTOFF, rel_tol=0.0, abs_tol=1e-14
    ):
        raise RuntimeError("The PR66 cutoff-match anchor no longer matches")
    for design in designs:
        write_point_manifest(
            work_root=args.work_root,
            design=design,
            build=build,
            input_metadata=input_metadata,
        )

    invocation_path = args.work_root / "invocations.json"
    invocations = load_json(invocation_path) if invocation_path.is_file() else {
        "schema": "sjaracne-brca100-pr67-p-sweep-invocations-v1",
        "invocations": [],
    }
    invocation = {
        "started_at_utc": utc_now(),
        "status": "running",
        "phase": args.phase,
        "points": [str(item["key"]) for item in designs],
        "drivers": [driver.key for driver in drivers],
        "seed_start": args.seed_start,
        "seed_end": args.seed_end,
        "workers": args.workers,
    }
    invocations["invocations"].append(invocation)
    atomic_json(invocation_path, invocations)

    if args.phase in ("infer", "all"):
        tasks = [
            (design, driver, seed)
            for seed in seeds
            for driver in drivers
            for design in designs
        ]
        completed = 0
        skipped = 0
        core.console(f"[INFER] {len(tasks)} jobs with {args.workers} workers")
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(
                    run_seed_job,
                    design=design,
                    build=build,
                    driver=driver,
                    input_root=args.work_root / "inputs",
                    results_root=args.work_root / "results",
                    seed=seed,
                    expression_ids=expression_ids,
                    driver_ids=driver_ids[driver.key],
                ): (design, driver, seed)
                for design, driver, seed in tasks
            }
            for future in as_completed(futures):
                label, was_skipped, record = future.result()
                completed += 1
                skipped += int(was_skipped)
                core.console(
                    f"[INFER {completed}/{len(tasks)}] {label} "
                    f"{'resume' if was_skipped else 'done'}; "
                    f"edges={record['adjacency']['edges']} "
                    f"wall={record['gnu_time']['elapsed_s']:.2f}s"
                )
        # Aggregate all completed point directories, not only this invocation.
        complete_designs = [
            point_design(point, model)
            for point in SWEEP_POINTS
            if (args.work_root / "results" / point.key).is_dir()
        ]
        aggregate_seed_manifest(args.work_root / "results", complete_designs, list(DRIVERS))
        invocation.update(
            {
                "inference_jobs": len(tasks),
                "inference_resumed_jobs": skipped,
                "inference_new_jobs": len(tasks) - skipped,
            }
        )

    if args.phase in ("consensus", "all"):
        if seeds != list(range(1, 101)):
            raise ValueError("Consensus requires the complete fixed seed range 1..100")
        # Process from strictest to loosest p, TF before SIG, while serialized.
        for design in sorted(designs, key=lambda item: float(item["p_value"])):
            stage = validation_stage(design)
            for driver in drivers:
                core.console(f"[CONSENSUS] {design['key']}/{driver.key}")
                record = core.run_consensus(
                    benchmark_repo=args.benchmark_repo,
                    work_root=args.work_root,
                    stage=stage,
                    driver=driver,
                    seeds=seeds,
                    expression_ids=expression_ids,
                    driver_ids=driver_ids[driver.key],
                )
                core.console(
                    f"[CONSENSUS] {design['key']}/{driver.key} "
                    f"edges={record['ncol']['edges']} "
                    f"rss={record['gnu_time']['max_rss_kib']} KiB"
                )

    invocation["status"] = "complete"
    invocation["finished_at_utc"] = utc_now()
    atomic_json(invocation_path, invocations)
    core.console(f"[DONE] {args.phase} at {args.work_root}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
