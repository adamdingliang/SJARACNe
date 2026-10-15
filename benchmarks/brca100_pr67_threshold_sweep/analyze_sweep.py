#!/usr/bin/env python3
"""Validate and summarize a matched PR67 per-subsample p-value sweep.

The topology screen may identify a provisional operating point, but it does
not establish biological optimality or empirical FDR.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
import math
import platform
from pathlib import Path
import sys
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "sjaracne-brca100-pr67-p-sweep-v1"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DRIVERS = ("tf", "sig")
CANDIDATE_FILES = {"tf": "BRCA100_TF.txt", "sig": "BRCA100_SIG.txt"}
EXPECTED_DRIVER_COUNTS = {"tf": 2608, "sig": 10680}
EXPECTED_EXPRESSION_COUNT = 28278
EXPECTED_CANDIDATE_PAIR_TESTS = {"tf": 73743244, "sig": 301984614}
PR67_COMMIT = "7633ebb4a0d966dbda15a4e32d0efa492fb71aeb"
NULL_MODEL_SHA256 = "e3a8522682a8ea239821aaa10b12db72d00e07bfdcad43599d8e76a06be80944"
EXPECTED_INPUT_SHA256 = {
    "BRCA100.exp": "ad8a334f5f8cdf46a1000d3ee259b35258a18b3da2e314bb3a0cf7a421d98bc8",
    "BRCA100_TF.txt": "9b1219a489b99432175e4c4ad46add7b06f25aae388ee8dd3261fa91e4c43ffd",
    "BRCA100_SIG.txt": "16ca27df655f16684f880a4ad719c4e2ae3f8dc0d7e6b9eccdd24cd97c40797c",
}
EXPECTED_POINTS = {
    "p1e-07": 1e-7,
    "p1e-06": 1e-6,
    "p1e-05": 1e-5,
    "p2e-05": 2e-5,
    "p5e-05": 5e-5,
    "p1e-04": 1e-4,
    "p2e-04": 2e-4,
    "p3e-04": 3e-4,
    "p_pr66_cutoff_match": 0.000352804562601613,
    "p4e-04": 4e-4,
    "p5e-04": 5e-4,
    "p7p5e-04": 7.5e-4,
    "p1e-03": 1e-3,
}
EXPECTED_SEEDS = tuple(range(1, 101))
HELD_OUT_P_MIN = 2e-5
HELD_OUT_P_MAX = 2e-3
EXACT_GATE2_GRID = frozenset((2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3))
ENGINEERING_FLOORS = {
    "active_driver_fraction": 0.90,
    "largest_weak_component_fraction_incident": 0.95,
    "incident_node_fraction_expression": 0.70,
}
NETWORK_FILENAME = "consensus_network_ncol_.txt"
SUPPORT_FILENAME = "consensus_support.tsv"
EXPECTED_NCOL_COLUMNS = (
    "source",
    "target",
    "source.symbol",
    "target.symbol",
    "MI",
    "pearson",
    "spearman",
    "slope",
    "p-value",
)
EXPECTED_SUPPORT_COLUMNS = (
    "source",
    "target",
    "consensus_MI",
    "support_count",
    "support_fraction",
    "mean_observed_MI",
    "consensus_MI_roundtrip_match",
)
TARGET_COUNT_THRESHOLDS = (10, 20, 30, 50, 100)
SUPPORT_FRACTION_THRESHOLDS = (0.20, 0.50, 0.80)
REPO_ROOT = Path(__file__).resolve().parents[2]
SVG_METADATA = {
    "Creator": "SJARACNe BRCA100 PR67 threshold sweep",
    "Date": None,
}


def compact_plot_p_label(p_value: float) -> str:
    """Format plot-only p labels compactly without changing stored values."""
    exponent = math.floor(math.log10(p_value))
    coefficient = p_value / (10.0 ** exponent)
    return f"{coefficient:.3g}e{exponent}"


def calibration_point_class(p_value: float) -> str:
    """Return the reporting class supported by the held-out null validation."""
    exact_grid = any(
        math.isclose(p_value, value, rel_tol=1e-14, abs_tol=0.0)
        for value in EXACT_GATE2_GRID
    )
    if exact_grid:
        return "exact-gate2-grid"
    if HELD_OUT_P_MIN <= p_value <= HELD_OUT_P_MAX:
        return "interpolation-within-accepted-range"
    return "outside-held-out-range"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_directory(path: Path) -> str:
    if not path.is_dir():
        raise ValueError(f"Missing directory to hash: {path}")
    digest = hashlib.sha256()
    for child in sorted(item for item in path.rglob("*") if item.is_file()):
        digest.update(child.relative_to(path).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(bytes.fromhex(sha256_file(child)))
    return digest.hexdigest()


def json_fingerprint(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def file_set_sha256(paths: list[Path], root: Path) -> str:
    digest = hashlib.sha256()
    for path in paths:
        digest.update(path.relative_to(root).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(bytes.fromhex(sha256_file(path)))
    return digest.hexdigest()


def require_sha256(value: Any, *, field: str, source: Path) -> str:
    text = str(value)
    if len(text) != 64 or any(character not in "0123456789abcdef" for character in text):
        raise ValueError(f"Invalid {field} in {source}: {value!r}")
    return text


def first_present(record: dict[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in record and record[key] is not None:
            return record[key]
    return None


def finite_float(value: Any, *, field: str, source: Path) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid {field} in {source}: {value!r}") from error
    if not math.isfinite(result):
        raise ValueError(f"Non-finite {field} in {source}: {value!r}")
    return result


def read_json_value(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"Cannot read JSON {path}: {error}") from error


def read_json(path: Path) -> dict[str, Any]:
    value = read_json_value(path)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def validate_input_metadata(
    input_root: Path,
    metadata: Any,
    *,
    source: Path,
) -> None:
    expected_keys = set(EXPECTED_INPUT_SHA256) | {"expression_id_count"}
    if not isinstance(metadata, dict) or set(metadata) != expected_keys:
        raise ValueError(f"Unexpected input metadata keys in {source}")
    if metadata.get("expression_id_count") != {"count": EXPECTED_EXPRESSION_COUNT}:
        raise ValueError(f"Unexpected expression ID count in {source}")
    for filename, expected_hash in EXPECTED_INPUT_SHA256.items():
        path = input_root / filename
        if not path.is_file():
            raise ValueError(f"Missing pinned input: {path}")
        record = metadata.get(filename)
        if not isinstance(record, dict):
            raise ValueError(f"Invalid {filename} metadata in {source}")
        actual_hash = sha256_file(path)
        if actual_hash != expected_hash or record.get("sha256") != expected_hash:
            raise ValueError(
                f"Pinned input SHA256 mismatch for {filename}: actual={actual_hash}, "
                f"metadata={record.get('sha256')}, expected={expected_hash}"
            )
        if int(record.get("bytes", -1)) != path.stat().st_size:
            raise ValueError(f"Input byte count mismatch for {filename} in {source}")


def validate_build(
    work_root: Path,
    *,
    binary_sha256: str,
    config_sha256: str,
) -> dict[str, Any]:
    build_root = work_root / "builds" / "pr67_7633ebb"
    manifest_path = build_root / "build_manifest.json"
    manifest = read_json(manifest_path)
    binary = build_root / "bin" / "sjaracne.exe"
    config_root = build_root / "source" / "SJARACNe" / "config"
    model = config_root / "apmi_null" / "apmi_null_m00080_npar040.model"
    checks = {
        "stage": "pr67_7633ebb",
        "commit": PR67_COMMIT,
        "binary_sha256": binary_sha256,
        "config_sha256": config_sha256,
        "null_model_sha256": NULL_MODEL_SHA256,
    }
    for field, expected in checks.items():
        if manifest.get(field) != expected:
            raise ValueError(
                f"Build manifest {field}={manifest.get(field)!r}, expected {expected!r}"
            )
    if not binary.is_file() or sha256_file(binary) != binary_sha256:
        raise ValueError(f"Built PR67 binary does not match its pinned hash: {binary}")
    if sha256_directory(config_root) != config_sha256:
        raise ValueError(f"PR67 config directory does not match its pinned hash: {config_root}")
    if not model.is_file() or sha256_file(model) != NULL_MODEL_SHA256:
        raise ValueError(f"PR67 AP-MI null model does not match its pinned hash: {model}")
    return {
        "manifest_path": manifest_path,
        "manifest_sha256": sha256_file(manifest_path),
        "binary_path": binary,
        "config_root": config_root,
        "model_path": model,
    }


def validate_sweep_design(
    work_root: Path,
) -> tuple[dict[str, Any], dict[str, dict[str, Any]], str, dict[str, Any]]:
    path = work_root / "sweep_design.json"
    design = read_json(path)
    if design.get("schema") != "sjaracne-brca100-pr67-p-sweep-v2":
        raise ValueError(f"Unexpected sweep design schema in {path}")
    if design.get("commit") != PR67_COMMIT:
        raise ValueError(f"Sweep design does not pin PR67 commit {PR67_COMMIT}")
    binary_hash = require_sha256(
        design.get("binary_sha256"), field="binary_sha256", source=path
    )
    config_hash = require_sha256(
        design.get("config_sha256"), field="config_sha256", source=path
    )
    if design.get("null_model_sha256") != NULL_MODEL_SHA256:
        raise ValueError(f"Sweep design does not pin null model {NULL_MODEL_SHA256}")
    fixed = design.get("fixed_parameters")
    expected_fixed = {
        "sampling": "fixed 80% without replacement",
        "m": 80,
        "npar": 40,
        "dpi_epsilon": 0,
        "consensus_p": 1e-5,
        "seeds": list(EXPECTED_SEEDS),
    }
    if fixed != expected_fixed:
        raise ValueError(
            f"Sweep design fixed parameters changed: observed={fixed!r}, "
            f"expected={expected_fixed!r}"
        )
    validate_input_metadata(work_root / "inputs", design.get("inputs"), source=path)
    raw_points = design.get("all_points")
    if not isinstance(raw_points, list) or len(raw_points) != len(EXPECTED_POINTS):
        raise ValueError(f"Unexpected all_points list in {path}")
    design_points: dict[str, dict[str, Any]] = {}
    for record in raw_points:
        if not isinstance(record, dict) or not isinstance(record.get("key"), str):
            raise ValueError(f"Invalid point entry in {path}")
        key = record["key"]
        if key in design_points:
            raise ValueError(f"Duplicate point {key} in {path}")
        design_points[key] = record
    if set(design_points) != set(EXPECTED_POINTS):
        raise ValueError(
            f"Sweep point keys changed: observed={sorted(design_points)}, "
            f"expected={sorted(EXPECTED_POINTS)}"
        )
    for key, expected_p in EXPECTED_POINTS.items():
        observed_p = finite_float(
            design_points[key].get("p_value"), field="p_value", source=path
        )
        if not math.isclose(observed_p, expected_p, rel_tol=1e-14, abs_tol=0.0):
            raise ValueError(f"Sweep point {key} changed from p={expected_p:.17g}")
    build = validate_build(
        work_root, binary_sha256=binary_hash, config_sha256=config_hash
    )
    return design, design_points, sha256_file(path), build


def point_metadata(
    results_root: Path,
    design: dict[str, Any],
    design_points: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    """Discover and verify every point against the immutable sweep design."""
    manifest_roots = {
        path.parent.name: path.parent
        for path in results_root.glob("*/point_manifest.json")
    }
    if set(manifest_roots) != set(design_points):
        raise ValueError(
            "Point manifests do not exactly match sweep_design.json: "
            f"observed={sorted(manifest_roots)}, expected={sorted(design_points)}"
        )
    points: list[dict[str, Any]] = []
    for point_key in sorted(manifest_roots):
        point_root = manifest_roots[point_key]
        manifest_path = point_root / "point_manifest.json"
        record = read_json(manifest_path)
        if record.get("schema") != "sjaracne-brca100-pr67-p-sweep-point-v1":
            raise ValueError(f"Unexpected point schema in {manifest_path}")
        if record.get("key") != point_root.name:
            raise ValueError(
                f"Point manifest key {record.get('key')!r} does not equal directory "
                f"key {point_root.name!r}"
            )
        for field, expected in design_points[point_key].items():
            if record.get(field) != expected:
                raise ValueError(
                    f"Point {point_key} field {field} disagrees with sweep_design.json"
                )
        fixed_checks = {
            "commit": PR67_COMMIT,
            "binary_sha256": design["binary_sha256"],
            "config_sha256": design["config_sha256"],
            "null_model_sha256": NULL_MODEL_SHA256,
            "sampling": "fixed 80% without replacement",
            "m": 80,
            "npar": 40,
            "dpi_epsilon": 0,
            "consensus_p": 1e-5,
            "seeds": list(EXPECTED_SEEDS),
            "inputs": design["inputs"],
        }
        for field, expected in fixed_checks.items():
            if record.get(field) != expected:
                raise ValueError(
                    f"Point {point_key} changes fixed field {field}: "
                    f"{record.get(field)!r} != {expected!r}"
                )
        p_value = finite_float(record.get("p_value"), field="p_value", source=manifest_path)
        if not 0.0 < p_value < 1.0:
            raise ValueError(f"p_value must be between zero and one in {manifest_path}")
        for driver in DRIVERS:
            if not (point_root / driver).is_dir():
                raise ValueError(f"Missing {driver} arm under {point_root}")
        mi_cutoff_value = first_present(record, "mi_cutoff", "MI_cutoff", "cutoff")
        mi_cutoff = (
            math.nan
            if mi_cutoff_value is None
            else finite_float(mi_cutoff_value, field="mi_cutoff", source=manifest_path)
        )
        points.append(
            {
                "p_key": point_root.name,
                "p_value": p_value,
                "log10_p": math.log10(p_value),
                "p_token": str(record.get("p_token", point_root.name)),
                "p_label": str(
                    record.get("label", record.get("p_token", point_root.name))
                ),
                "role": str(record.get("role", "")),
                "mi_cutoff": mi_cutoff,
                "validation_class": str(record.get("validation_class", "")),
                "calibration_point_class": calibration_point_class(p_value),
                "tail_extrapolated": bool(record.get("tail_extrapolated", False)),
                "commit": str(first_present(record, "commit", "commit_sha") or ""),
                "model_sha256": str(
                    first_present(record, "model_sha256", "null_model_sha256") or ""
                ),
                "manifest_sha256": sha256_file(manifest_path),
                "manifest": record,
                "root": point_root,
            }
        )
    points.sort(key=lambda point: (point["p_value"], point["p_key"]))
    values = [point["p_value"] for point in points]
    if len(values) != len(set(values)):
        raise ValueError(f"Duplicate p_value entries in {results_root}")
    return points


def candidate_lists(input_root: Path) -> dict[str, list[str]]:
    candidates: dict[str, list[str]] = {}
    for driver, filename in CANDIDATE_FILES.items():
        path = input_root / filename
        if not path.is_file():
            raise ValueError(f"Missing fixed candidate list: {path}")
        values = [
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        if (
            len(values) != EXPECTED_DRIVER_COUNTS[driver]
            or len(values) != len(set(values))
        ):
            raise ValueError(
                f"Expected {EXPECTED_DRIVER_COUNTS[driver]} unique {driver} candidates "
                f"in {path}, found {len(values)}"
            )
        candidates[driver] = values
    return candidates


def expression_annotations(path: Path) -> dict[str, str]:
    if not path.is_file():
        raise ValueError(f"Missing fixed expression matrix: {path}")
    annotations: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if header[:2] != ["isoformId", "geneSymbol"]:
            raise ValueError(f"Invalid expression-matrix header: {path}")
        for line_number, line in enumerate(handle, 2):
            fields = line.rstrip("\r\n").split("\t", 2)
            if len(fields) < 2:
                raise ValueError(f"Missing gene symbol at {path}:{line_number}")
            identifier, gene_symbol = fields[:2]
            if not identifier or identifier in annotations:
                raise ValueError(
                    f"Empty or duplicate expression ID at {path}:{line_number}"
                )
            annotations[identifier] = gene_symbol
    if len(annotations) != EXPECTED_EXPRESSION_COUNT:
        raise ValueError(
            f"Expected {EXPECTED_EXPRESSION_COUNT:,} expression IDs in {path}, "
            f"found {len(annotations)}"
        )
    return annotations


def eligible_candidate_pair_tests(
    candidates: list[str], expression_annotations: dict[str, str]
) -> int:
    """Count the directed MI tests SJARACNe can actually evaluate.

    SJARACNe suppresses a candidate paired with itself and with every other
    row carrying the same non-`---` geneSymbol. An empty geneSymbol therefore
    suppresses same-label rows, while `---` intentionally does not.
    """
    missing = set(candidates) - set(expression_annotations)
    if missing:
        raise ValueError(f"Candidate list contains missing expression IDs: {len(missing)}")
    label_counts = Counter(expression_annotations.values())
    expression_count = len(expression_annotations)
    result = 0
    for candidate in candidates:
        label = expression_annotations[candidate]
        same_gene_others = label_counts[label] - 1 if label != "---" else 0
        result += expression_count - 1 - same_gene_others
    return result


def arm_file(arm_root: Path, filename: str) -> Path:
    candidates = (arm_root / "consensus" / filename, arm_root / filename)
    matches = [path for path in candidates if path.is_file()]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one {filename} directly under {arm_root} or its "
            f"consensus directory; found {matches}"
        )
    return matches[0]


def load_network(arm_root: Path) -> pd.DataFrame:
    path = arm_file(arm_root, NETWORK_FILENAME)
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = tuple(handle.readline().rstrip("\r\n").split("\t"))
    if header != EXPECTED_NCOL_COLUMNS:
        raise ValueError(f"Unexpected consensus-network columns in {path}: {header}")
    frame = pd.read_csv(
        path,
        sep="\t",
        usecols=["source", "target", "MI"],
        dtype={"source": str, "target": str},
    )
    if frame[["source", "target"]].isna().any().any():
        raise ValueError(f"Missing edge endpoint in {path}")
    if frame.duplicated(["source", "target"]).any():
        raise ValueError(f"Duplicate directed edges in {path}")
    frame["MI"] = pd.to_numeric(frame["MI"], errors="raise")
    if len(frame) and (not np.isfinite(frame["MI"]).all() or (frame["MI"] <= 0).any()):
        raise ValueError(f"MI must be finite and positive in {path}")
    return frame.set_index(["source", "target"], drop=False).sort_index()


def load_support(arm_root: Path) -> pd.DataFrame:
    path = arm_file(arm_root, SUPPORT_FILENAME)
    frame = pd.read_csv(path, sep="\t", dtype={"source": str, "target": str})
    if tuple(frame.columns) != EXPECTED_SUPPORT_COLUMNS:
        raise ValueError(f"Unexpected support columns in {path}: {tuple(frame.columns)}")
    if frame[["source", "target"]].isna().any().any():
        raise ValueError(f"Missing support edge endpoint in {path}")
    if frame.duplicated(["source", "target"]).any():
        raise ValueError(f"Duplicate directed edges in {path}")
    for column in ("support_count", "support_fraction"):
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    if len(frame) and (
        not np.isfinite(frame[["support_count", "support_fraction"]]).all().all()
        or (frame["support_count"] < 1).any()
        or not frame["support_fraction"].between(0.0, 1.0).all()
    ):
        raise ValueError(f"Invalid support values in {path}")
    if "consensus_MI_roundtrip_match" in frame.columns and len(frame):
        if not (pd.to_numeric(frame["consensus_MI_roundtrip_match"]) == 1).all():
            raise ValueError(f"Consensus MI roundtrip mismatch in {path}")
    return frame.set_index(["source", "target"], drop=False).sort_index()


def merge_network_support(arm_root: Path) -> pd.DataFrame:
    network = load_network(arm_root)
    support = load_support(arm_root)
    if not network.index.equals(support.index):
        network_edges = set(network.index)
        support_edges = set(support.index)
        raise ValueError(
            f"Network/support edge mismatch in {arm_root}: "
            f"network-only={len(network_edges - support_edges)}, "
            f"support-only={len(support_edges - network_edges)}"
        )
    result = network.copy()
    result["support_count"] = support["support_count"]
    result["support_fraction"] = support["support_fraction"]
    if "consensus_MI" in support.columns and len(result):
        support_mi = pd.to_numeric(support["consensus_MI"], errors="raise")
        if not np.allclose(result["MI"], support_mi, rtol=0.0, atol=5.00001e-5):
            raise ValueError(f"Network/support MI mismatch in {arm_root}")
    return result


def command_flag(command: list[str], flag: str, *, source: Path) -> str:
    if command.count(flag) != 1:
        raise ValueError(f"Expected one {flag} in seed command from {source}")
    position = command.index(flag)
    if position + 1 >= len(command):
        raise ValueError(f"Missing {flag} value in seed command from {source}")
    return command[position + 1]


def validate_seed_command(
    command: Any,
    *,
    source: Path,
    work_root: Path,
    build: dict[str, Any],
    point: dict[str, Any],
    driver: str,
    seed: int,
    adjacency_path: Path,
) -> list[str]:
    if not isinstance(command, list) or not all(isinstance(value, str) for value in command):
        raise ValueError(f"Invalid seed command in {source}")
    expected_flags = (
        "-i", "-l", "-s", "-p", "-e", "-a", "-H", "-N", "-S", "-o", "-u", "-M"
    )
    if tuple(command[1::2]) != expected_flags or len(command) != 25:
        raise ValueError(f"Unexpected seed command structure in {source}: {command}")
    path_expectations = {
        "executable": (Path(command[0]), build["binary_path"]),
        "-i": (Path(command_flag(command, "-i", source=source)), work_root / "inputs/BRCA100.exp"),
        "-l": (
            Path(command_flag(command, "-l", source=source)),
            work_root / "inputs" / CANDIDATE_FILES[driver],
        ),
        "-s": (
            Path(command_flag(command, "-s", source=source)),
            work_root / "inputs" / CANDIDATE_FILES[driver],
        ),
        "-H": (Path(command_flag(command, "-H", source=source)), build["config_root"]),
        "-o": (Path(command_flag(command, "-o", source=source)), adjacency_path),
        "-M": (Path(command_flag(command, "-M", source=source)), build["model_path"]),
    }
    for label, (observed, expected) in path_expectations.items():
        if observed.resolve() != expected.resolve():
            raise ValueError(
                f"Seed command {label} path mismatch in {source}: "
                f"{observed} != {expected}"
            )
    scalar_expectations = {
        "-e": "0",
        "-a": "adaptive_partitioning",
        "-N": "40",
        "-S": str(seed),
        "-u": "80%",
    }
    for flag, expected in scalar_expectations.items():
        if command_flag(command, flag, source=source) != expected:
            raise ValueError(f"Seed command {flag} changed in {source}")
    observed_p = finite_float(
        command_flag(command, "-p", source=source), field="command -p", source=source
    )
    if not math.isclose(observed_p, point["p_value"], rel_tol=1e-14, abs_tol=0.0):
        raise ValueError(f"Seed command p-value disagrees with point in {source}")
    return command


def validate_completed_seeds(
    work_root: Path,
    points: list[dict[str, Any]],
    design_points: dict[str, dict[str, Any]],
    design: dict[str, Any],
    build: dict[str, Any],
) -> tuple[
    dict[tuple[str, str], set[int]],
    pd.DataFrame,
    dict[tuple[str, str], dict[str, Any]],
]:
    manifest_path = work_root / "results" / "run_manifest.tsv"
    if not manifest_path.is_file():
        raise ValueError(
            f"Missing completed-run evidence {manifest_path}; planned point seeds are "
            "not accepted as execution evidence"
        )
    run_manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
    required_columns = {
        "point", "p_value", "mi_cutoff", "validation_class", "commit", "driver",
        "seed", "binary_sha256", "edges", "source_rows", "adjacency_bytes",
        "adjacency_sha256", "data_sha256", "stderr_bytes",
    }
    if not required_columns.issubset(run_manifest.columns):
        raise ValueError(
            "run_manifest.tsv lacks columns: "
            f"{sorted(required_columns - set(run_manifest.columns))}"
        )
    expected_rows = len(points) * len(DRIVERS) * len(EXPECTED_SEEDS)
    if len(run_manifest) != expected_rows:
        raise ValueError(
            f"Expected exactly {expected_rows} completed run-manifest rows, "
            f"found {len(run_manifest)}"
        )
    if run_manifest.duplicated(["point", "driver", "seed"]).any():
        raise ValueError(f"Duplicate arm/seed rows in {manifest_path}")

    seed_sets: dict[tuple[str, str], set[int]] = {}
    summary_rows: list[dict[str, Any]] = []
    arms: dict[tuple[str, str], dict[str, Any]] = {}
    for point in points:
        point_key = point["p_key"]
        for driver in DRIVERS:
            arm_root = point["root"] / driver
            metadata_root = arm_root / "seed_metadata"
            adjacency_root = arm_root / "adjacency"
            metadata_paths = sorted(metadata_root.glob("TF_run_*.json"))
            adjacency_paths = sorted(adjacency_root.glob("TF_run_*.adj"))
            expected_metadata_names = [f"TF_run_{seed:03d}.json" for seed in EXPECTED_SEEDS]
            expected_adjacency_names = [f"TF_run_{seed:03d}.adj" for seed in EXPECTED_SEEDS]
            if [path.name for path in metadata_paths] != expected_metadata_names:
                raise ValueError(f"Seed metadata is not exactly seeds 1..100 in {metadata_root}")
            if [path.name for path in adjacency_paths] != expected_adjacency_names:
                raise ValueError(
                    f"Adjacency inputs are not exactly seeds 1..100 in "
                    f"{adjacency_root}"
                )
            subset = run_manifest[
                (run_manifest["point"] == point_key)
                & (run_manifest["driver"] == driver)
            ].copy()
            subset["seed_int"] = pd.to_numeric(subset["seed"], errors="raise").astype(int)
            subset = subset.sort_values("seed_int")
            if subset["seed_int"].tolist() != list(EXPECTED_SEEDS):
                raise ValueError(
                    f"run_manifest.tsv is not exactly seeds 1..100 for "
                    f"{point_key}/{driver}"
                )
            row_by_seed = {int(row.seed_int): row for row in subset.itertuples(index=False)}
            adjacency_hashes: list[str] = []
            for seed, metadata_path, adjacency_path in zip(
                EXPECTED_SEEDS, metadata_paths, adjacency_paths
            ):
                record = read_json(metadata_path)
                fixed_checks = {
                    "schema": "sjaracne-brca100-pr67-p-sweep-seed-v1",
                    "point": design_points[point_key],
                    "commit": PR67_COMMIT,
                    "driver": driver,
                    "seed": seed,
                    "binary_sha256": design["binary_sha256"],
                    "config_sha256": design["config_sha256"],
                    "null_model_sha256": NULL_MODEL_SHA256,
                }
                for field, expected in fixed_checks.items():
                    if record.get(field) != expected:
                        raise ValueError(
                            f"Seed metadata changes {field} in {metadata_path}: "
                            f"{record.get(field)!r} != {expected!r}"
                        )
                command = validate_seed_command(
                    record.get("command"),
                    source=metadata_path,
                    work_root=work_root,
                    build=build,
                    point=point,
                    driver=driver,
                    seed=seed,
                    adjacency_path=adjacency_path,
                )
                adjacency_hash = sha256_file(adjacency_path)
                adjacency = record.get("adjacency")
                if not isinstance(adjacency, dict):
                    raise ValueError(f"Missing adjacency metadata in {metadata_path}")
                if (
                    adjacency.get("full_sha256") != adjacency_hash
                    or int(adjacency.get("bytes", -1)) != adjacency_path.stat().st_size
                ):
                    raise ValueError(f"Adjacency file disagrees with {metadata_path}")
                fingerprint_payload = {
                    "schema": "sjaracne-brca100-pr67-p-sweep-seed-v1",
                    "point": design_points[point_key],
                    "commit": PR67_COMMIT,
                    "binary_sha256": design["binary_sha256"],
                    "config_sha256": design["config_sha256"],
                    "null_model_sha256": NULL_MODEL_SHA256,
                    "driver": driver,
                    "driver_sha256": EXPECTED_INPUT_SHA256[CANDIDATE_FILES[driver]],
                    "expression_sha256": EXPECTED_INPUT_SHA256["BRCA100.exp"],
                    "seed": seed,
                    "command_without_output": [
                        "<OUTPUT>" if value == str(adjacency_path) else value
                        for value in command
                    ],
                }
                if record.get("fingerprint") != json_fingerprint(fingerprint_payload):
                    raise ValueError(f"Seed fingerprint mismatch in {metadata_path}")
                stdout_path = arm_root / "logs" / f"TF_run_{seed:03d}.stdout.log"
                stderr_path = arm_root / "logs" / f"TF_run_{seed:03d}.stderr.log"
                if (
                    not stdout_path.is_file()
                    or not stderr_path.is_file()
                    or record.get("stdout_sha256") != sha256_file(stdout_path)
                    or record.get("stderr_sha256") != sha256_file(stderr_path)
                    or int(record.get("stderr_bytes", -1)) != stderr_path.stat().st_size
                ):
                    raise ValueError(f"Seed logs disagree with {metadata_path}")
                row = row_by_seed[seed]
                row_checks = {
                    "commit": PR67_COMMIT,
                    "validation_class": point["validation_class"],
                    "binary_sha256": design["binary_sha256"],
                    "adjacency_sha256": adjacency_hash,
                    "data_sha256": adjacency.get("data_sha256"),
                    "edges": str(adjacency.get("edges")),
                    "source_rows": str(adjacency.get("source_rows")),
                    "adjacency_bytes": str(adjacency.get("bytes")),
                    "stderr_bytes": str(record.get("stderr_bytes")),
                }
                for field, expected in row_checks.items():
                    if str(getattr(row, field)) != str(expected):
                        raise ValueError(
                            f"run_manifest.tsv {field} disagrees with {metadata_path}"
                        )
                if not math.isclose(float(row.p_value), point["p_value"], rel_tol=1e-14):
                    raise ValueError(
                        f"run_manifest.tsv p_value mismatch for "
                        f"{point_key}/{driver}/{seed}"
                    )
                if not math.isclose(float(row.mi_cutoff), point["mi_cutoff"], rel_tol=1e-14):
                    raise ValueError(
                        f"run_manifest.tsv mi_cutoff mismatch for "
                        f"{point_key}/{driver}/{seed}"
                    )
                adjacency_hashes.append(adjacency_hash)
            seed_set = set(EXPECTED_SEEDS)
            seed_sets[(point_key, driver)] = seed_set
            encoded = ",".join(str(seed) for seed in EXPECTED_SEEDS).encode("ascii")
            summary_rows.append(
                {
                    "p_key": point_key,
                    "p_value": point["p_value"],
                    "driver_class": driver,
                    "runs": len(EXPECTED_SEEDS),
                    "seed_min": EXPECTED_SEEDS[0],
                    "seed_max": EXPECTED_SEEDS[-1],
                    "seed_set_sha256": hashlib.sha256(encoded).hexdigest(),
                    "seed_source": str(manifest_path),
                    "command_p_fully_checked": True,
                }
            )
            arms[(point_key, driver)] = {
                "arm_root": arm_root,
                "adjacency_paths": adjacency_paths,
                "adjacency_hashes": adjacency_hashes,
                "adjacency_set_sha256": file_set_sha256(adjacency_paths, arm_root),
                "seed_metadata_set_sha256": file_set_sha256(metadata_paths, arm_root),
                "run_manifest_rows_sha256": json_fingerprint(
                    subset.drop(columns="seed_int").to_dict(orient="records")
                ),
            }
    return seed_sets, pd.DataFrame(summary_rows), arms


def validate_arm_artifacts(
    work_root: Path,
    points: list[dict[str, Any]],
    arms: dict[tuple[str, str], dict[str, Any]],
    *,
    sweep_design_sha256: str,
) -> dict[tuple[str, str], dict[str, Any]]:
    consensus_script = REPO_ROOT / "SJARACNe" / "bin" / "create_consensus_network.py"
    support_source = (
        REPO_ROOT
        / "benchmarks"
        / "brca100_netbid_qc"
        / "summarize_consensus_support.cpp"
    )
    support_binary = work_root / "tools" / "summarize_consensus_support"
    for path in (consensus_script, support_source, support_binary):
        if not path.is_file():
            raise ValueError(f"Missing provenance input: {path}")
    consensus_script_hash = sha256_file(consensus_script)
    support_source_hash = sha256_file(support_source)
    support_binary_hash = sha256_file(support_binary)

    aggregate_path = work_root / "results" / "support_summary_manifest.json"
    aggregate = read_json(aggregate_path)
    expected_point_order = [point["p_key"] for point in points]
    aggregate_checks = {
        "schema": "sjaracne-brca100-pr67-p-sweep-support-aggregate-v1",
        "sweep_design_sha256": sweep_design_sha256,
        "points": expected_point_order,
        "drivers": list(DRIVERS),
    }
    for field, expected in aggregate_checks.items():
        if aggregate.get(field) != expected:
            raise ValueError(
                f"Support aggregate {field} disagrees with the complete sweep: "
                f"{aggregate.get(field)!r} != {expected!r}"
            )
    aggregate_records = aggregate.get("records")
    if not isinstance(aggregate_records, list) or len(aggregate_records) != len(arms):
        raise ValueError(f"Support aggregate is incomplete: {aggregate_path}")
    aggregate_by_arm: dict[tuple[str, str], dict[str, Any]] = {}
    for record in aggregate_records:
        if not isinstance(record, dict):
            raise ValueError(f"Non-object support record in {aggregate_path}")
        key = (str(record.get("point")), str(record.get("driver")))
        if key in aggregate_by_arm:
            raise ValueError(f"Duplicate support aggregate arm {key}")
        aggregate_by_arm[key] = record
    if set(aggregate_by_arm) != set(arms):
        raise ValueError("Support aggregate arms do not exactly match the sweep arms")

    for point in points:
        point_key = point["p_key"]
        point_manifest_path = point["root"] / "point_manifest.json"
        point_manifest_hash = sha256_file(point_manifest_path)
        for driver in DRIVERS:
            key = (point_key, driver)
            arm = arms[key]
            arm_root = arm["arm_root"]
            pending_paths = (
                arm_root / "consensus_manifest.pending.json",
                arm_root / "support_summary_manifest.pending.json",
            )
            stale_pending = [path for path in pending_paths if path.exists()]
            if stale_pending:
                raise ValueError(f"Unresolved pending arm manifests: {stale_pending}")

            adjacency_hashes = arm["adjacency_hashes"]
            consensus_path = arm_root / "consensus" / NETWORK_FILENAME
            consensus_manifest_path = arm_root / "consensus_manifest.json"
            consensus_manifest = read_json(consensus_manifest_path)
            consensus_hash = sha256_file(consensus_path)
            expected_consensus_fingerprint = json_fingerprint(
                {
                    "stage": point_key,
                    "driver": driver,
                    "adjacency_hashes": adjacency_hashes,
                    "consensus_p": 1e-5,
                    "consensus_script_sha256": consensus_script_hash,
                }
            )
            if (
                consensus_manifest.get("stage") != point_key
                or consensus_manifest.get("driver") != driver
                or consensus_manifest.get("fingerprint")
                != expected_consensus_fingerprint
            ):
                raise ValueError(f"Consensus manifest provenance mismatch: {arm_root}")
            ncol = consensus_manifest.get("ncol")
            if not isinstance(ncol, dict) or (
                ncol.get("sha256") != consensus_hash
                or int(ncol.get("bytes", -1)) != consensus_path.stat().st_size
            ):
                raise ValueError(f"Consensus ncol hash/size mismatch: {arm_root}")
            consensus_root = arm_root / "consensus"
            auxiliary_hashes = {
                "consensus_3col_sha256": consensus_root / "consensus_network_3col_.txt",
                "parameter_info_sha256": consensus_root / "parameter_info_.txt",
                "bootstrap_info_sha256": consensus_root / "bootstrap_info_.txt",
            }
            for field, path in auxiliary_hashes.items():
                if not path.is_file() or consensus_manifest.get(field) != sha256_file(path):
                    raise ValueError(f"Consensus auxiliary hash mismatch for {path}")
            parameter_text = auxiliary_hashes["parameter_info_sha256"].read_text(
                encoding="utf-8"
            )
            if ">  Bootstrap No: 100" not in parameter_text:
                raise ValueError(f"Consensus does not record 100 inputs: {arm_root}")
            consensus_manifest_hash = sha256_file(consensus_manifest_path)

            support_path = arm_root / "consensus" / SUPPORT_FILENAME
            support_manifest_path = arm_root / "support_summary_manifest.json"
            support_manifest = read_json(support_manifest_path)
            support_hash = sha256_file(support_path)
            expected_support_fingerprint = json_fingerprint(
                {
                    "schema": "sjaracne-brca100-pr67-p-sweep-support-v1",
                    "point": point_key,
                    "p_value": point["p_value"],
                    "driver": driver,
                    "consensus_sha256": consensus_hash,
                    "consensus_fingerprint": expected_consensus_fingerprint,
                    "consensus_manifest_sha256": consensus_manifest_hash,
                    "point_manifest_sha256": point_manifest_hash,
                    "source_sha256": support_source_hash,
                    "binary_sha256": support_binary_hash,
                    "adjacency_sha256": adjacency_hashes,
                }
            )
            support_checks = {
                "schema": "sjaracne-brca100-pr67-p-sweep-support-v1",
                "fingerprint": expected_support_fingerprint,
                "point": point_key,
                "driver": driver,
                "source_sha256": support_source_hash,
                "binary_sha256": support_binary_hash,
                "point_manifest_sha256": point_manifest_hash,
                "consensus_sha256": consensus_hash,
                "consensus_fingerprint": expected_consensus_fingerprint,
                "consensus_manifest_sha256": consensus_manifest_hash,
                "adjacency_sha256": adjacency_hashes,
                "output_sha256": support_hash,
                "retained_edges": int(ncol.get("edges", -1)),
            }
            for field, expected in support_checks.items():
                if support_manifest.get(field) != expected:
                    raise ValueError(
                        f"Support manifest {field} mismatch for {point_key}/{driver}"
                    )
            if not math.isclose(
                float(support_manifest.get("p_value")),
                point["p_value"],
                rel_tol=1e-14,
                abs_tol=0.0,
            ):
                raise ValueError(f"Support manifest p-value mismatch: {arm_root}")
            if aggregate_by_arm[key] != support_manifest:
                raise ValueError(
                    f"Root and per-arm support manifests disagree for {point_key}/{driver}"
                )
            arm.update(
                {
                    "point_manifest_sha256": point_manifest_hash,
                    "consensus_sha256": consensus_hash,
                    "consensus_manifest_sha256": consensus_manifest_hash,
                    "consensus_fingerprint": expected_consensus_fingerprint,
                    "consensus_edges_manifest": int(ncol.get("edges", -1)),
                    "consensus_active_drivers_manifest": int(
                        ncol.get("active_drivers", -1)
                    ),
                    "consensus_incident_nodes_manifest": int(
                        ncol.get("incident_nodes", -1)
                    ),
                    "support_sha256": support_hash,
                    "support_manifest_sha256": sha256_file(support_manifest_path),
                    "support_fingerprint": expected_support_fingerprint,
                    "support_source_sha256": support_source_hash,
                    "support_binary_sha256": support_binary_hash,
                }
            )
    return arms


def weak_connectivity(edges: Iterable[tuple[str, str]]) -> tuple[int, int, int]:
    parent: dict[str, str] = {}
    size: dict[str, int] = {}

    def find(node: str) -> str:
        root = node
        while parent[root] != root:
            root = parent[root]
        while parent[node] != node:
            following = parent[node]
            parent[node] = root
            node = following
        return root

    def add(node: str) -> None:
        if node not in parent:
            parent[node] = node
            size[node] = 1

    for source, target in edges:
        add(source)
        add(target)
        left = find(source)
        right = find(target)
        if left == right:
            continue
        if size[left] < size[right]:
            left, right = right, left
        parent[right] = left
        size[left] += size[right]
    if not parent:
        return 0, 0, 0
    roots = [node for node in parent if parent[node] == node]
    return len(parent), len(roots), max(size[root] for root in roots)


def distribution_metrics(values: pd.Series, prefix: str) -> dict[str, float]:
    if values.empty:
        return {
            f"{prefix}_mean": math.nan,
            f"{prefix}_median": math.nan,
            f"{prefix}_q25": math.nan,
            f"{prefix}_q75": math.nan,
            f"{prefix}_min": math.nan,
            f"{prefix}_max": math.nan,
        }
    return {
        f"{prefix}_mean": float(values.mean()),
        f"{prefix}_median": float(values.median()),
        f"{prefix}_q25": float(values.quantile(0.25)),
        f"{prefix}_q75": float(values.quantile(0.75)),
        f"{prefix}_min": float(values.min()),
        f"{prefix}_max": float(values.max()),
    }


def zero_filled_target_sizes(
    frame: pd.DataFrame, candidates: list[str]
) -> pd.Series:
    """Count directed targets per candidate, retaining candidates with zero targets."""
    return (
        frame.groupby(level="source", sort=False)
        .size()
        .reindex(candidates, fill_value=0)
    )


def target_size_coverage_metrics(
    target_sizes: pd.Series, prefix: str
) -> dict[str, int | float]:
    """Summarize candidate-driver coverage at fixed target-count thresholds."""
    if target_sizes.empty:
        raise ValueError("Target-size coverage requires at least one candidate driver")
    result: dict[str, int | float] = {}
    for minimum_targets in TARGET_COUNT_THRESHOLDS:
        driver_count = int((target_sizes >= minimum_targets).sum())
        field = f"{prefix}_ge_{minimum_targets}"
        result[f"{field}_driver_count"] = driver_count
        result[f"{field}_driver_fraction"] = driver_count / len(target_sizes)
    return result


def support_threshold_token(threshold: float) -> str:
    """Return a stable field-name token for a configured support threshold."""
    percentage = threshold * 100.0
    rounded = round(percentage)
    if not math.isclose(percentage, rounded, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError(f"Support threshold is not an integer percentage: {threshold}")
    return f"{rounded:02d}pct"


def summarize_network(
    frame: pd.DataFrame,
    candidates: list[str],
    expression_ids: set[str],
    *,
    prefix: dict[str, Any],
    seed_runs: int,
    candidate_pair_tests: int,
) -> tuple[dict[str, Any], pd.Series]:
    source_set = set(frame["source"])
    unexpected = source_set - set(candidates)
    if unexpected:
        preview = sorted(unexpected)[:5]
        raise ValueError(f"Non-candidate sources for {prefix}: {preview}")
    incident_set = source_set | set(frame["target"])
    unexpected_incident = incident_set - expression_ids
    if unexpected_incident:
        preview = sorted(unexpected_incident)[:5]
        raise ValueError(
            f"Edges contain IDs absent from the expression matrix for {prefix}: "
            f"{preview}"
        )
    if not set(candidates).issubset(expression_ids):
        raise ValueError(
            f"Candidate list contains IDs absent from the expression matrix for {prefix}"
        )
    if seed_runs > 0 and len(frame):
        if (frame["support_count"] > seed_runs).any() or not np.allclose(
            frame["support_fraction"],
            frame["support_count"] / seed_runs,
            rtol=0.0,
            atol=1e-12,
        ):
            raise ValueError(
                f"Support counts/fractions disagree with {seed_runs} matched seeds for {prefix}"
            )
    target_sizes = zero_filled_target_sizes(frame, candidates)
    target_coverage = target_size_coverage_metrics(target_sizes, "target_size")
    recurrence_metrics: dict[str, Any] = {
        "recurrence_qualified_target_interpretation": (
            "zero-filled candidate-driver target counts restricted to consensus "
            "edges observed in at least the stated fraction of matched seed-level "
            "subnetworks; recurrence measures stability, not biological truth, an "
            "edge probability, or FDR"
        )
    }
    for support_threshold in SUPPORT_FRACTION_THRESHOLDS:
        threshold_token = support_threshold_token(support_threshold)
        qualified = frame.loc[frame["support_fraction"] >= support_threshold]
        qualified_target_sizes = zero_filled_target_sizes(qualified, candidates)
        metric_prefix = f"support_ge_{threshold_token}_target_size"
        recurrence_metrics[f"support_ge_{threshold_token}_edges"] = len(qualified)
        recurrence_metrics[f"support_ge_{threshold_token}_edge_fraction"] = (
            len(qualified) / len(frame) if len(frame) else 0.0
        )
        recurrence_metrics.update(
            {
                f"{metric_prefix}_mean_zero_filled": float(
                    qualified_target_sizes.mean()
                ),
                f"{metric_prefix}_median_zero_filled": float(
                    qualified_target_sizes.median()
                ),
                f"{metric_prefix}_q25_zero_filled": float(
                    qualified_target_sizes.quantile(0.25)
                ),
                f"{metric_prefix}_q75_zero_filled": float(
                    qualified_target_sizes.quantile(0.75)
                ),
                f"{metric_prefix}_max": int(qualified_target_sizes.max()),
            }
        )
        recurrence_metrics.update(
            target_size_coverage_metrics(qualified_target_sizes, metric_prefix)
        )
    incident_nodes, components, largest = weak_connectivity(frame.index)
    p_value = float(prefix.get("p_value", math.nan))
    row: dict[str, Any] = {
        **prefix,
        "seed_runs": seed_runs,
        "candidate_drivers": len(candidates),
        "candidate_pair_tests": candidate_pair_tests,
        "candidate_pair_tests_interpretation": (
            "exact directed candidate-target MI tests after SJARACNe same-accession "
            "and same-non-placeholder-geneSymbol suppression"
        ),
        "nominal_null_exceedances_before_DPI": candidate_pair_tests * p_value,
        "nominal_null_exceedances_interpretation": (
            "model-based independence-null expectation proxy per seed-level "
            "subnetwork before DPI; not a consensus-network expectation, post-DPI "
            "quantity, or FDR"
        ),
        "active_drivers": int((target_sizes > 0).sum()),
        "active_driver_fraction": float((target_sizes > 0).mean()),
        "consensus_edges": len(frame),
        "incident_nodes": incident_nodes,
        "expression_nodes": len(expression_ids),
        "incident_node_fraction_expression": incident_nodes / len(expression_ids),
        "weak_components": components,
        "largest_weak_component_nodes": largest,
        "largest_weak_component_fraction_incident": (
            largest / incident_nodes if incident_nodes else 0.0
        ),
        "target_size_mean_zero_filled": float(target_sizes.mean()),
        "target_size_median_zero_filled": float(target_sizes.median()),
        "target_size_q25_zero_filled": float(target_sizes.quantile(0.25)),
        "target_size_q75_zero_filled": float(target_sizes.quantile(0.75)),
        "target_size_max": int(target_sizes.max()),
        **target_coverage,
        **recurrence_metrics,
        **distribution_metrics(frame["support_fraction"], "support_fraction"),
        **distribution_metrics(frame["support_count"], "support_count"),
        **distribution_metrics(frame["MI"], "consensus_MI"),
    }
    return row, target_sizes


def safe_ratio(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else math.nan


def safe_correlation(left: pd.Series, right: pd.Series, method: str) -> float:
    if len(left) < 3 or left.nunique() < 2 or right.nunique() < 2:
        return math.nan
    return float(left.corr(right, method=method))


def overlap_row(
    reference: pd.DataFrame,
    comparison: pd.DataFrame,
    reference_sizes: pd.Series,
    comparison_sizes: pd.Series,
    *,
    prefix: dict[str, Any],
) -> dict[str, Any]:
    reference_edges = set(reference.index)
    comparison_edges = set(comparison.index)
    common = reference_edges & comparison_edges
    union = reference_edges | comparison_edges
    reference_active = set(reference["source"])
    comparison_active = set(comparison["source"])
    active_union = reference_active | comparison_active
    common_index = pd.MultiIndex.from_tuples(sorted(common), names=["source", "target"])
    if common:
        reference_mi = reference.loc[common_index, "MI"].reset_index(drop=True)
        comparison_mi = comparison.loc[common_index, "MI"].reset_index(drop=True)
        reference_support = reference.loc[common_index, "support_fraction"].reset_index(drop=True)
        comparison_support = comparison.loc[common_index, "support_fraction"].reset_index(drop=True)
    else:
        reference_mi = pd.Series(dtype=float)
        comparison_mi = pd.Series(dtype=float)
        reference_support = pd.Series(dtype=float)
        comparison_support = pd.Series(dtype=float)
    return {
        **prefix,
        "reference_edges": len(reference_edges),
        "comparison_edges": len(comparison_edges),
        "intersection_edges": len(common),
        "union_edges": len(union),
        "edge_jaccard": safe_ratio(len(common), len(union)),
        "reference_edge_retention": safe_ratio(len(common), len(reference_edges)),
        "comparison_edge_shared_fraction": safe_ratio(len(common), len(comparison_edges)),
        "lost_from_reference": len(reference_edges - comparison_edges),
        "gained_in_comparison": len(comparison_edges - reference_edges),
        "reference_active_drivers": len(reference_active),
        "comparison_active_drivers": len(comparison_active),
        "active_driver_intersection": len(reference_active & comparison_active),
        "active_driver_union": len(active_union),
        "active_driver_jaccard": safe_ratio(
            len(reference_active & comparison_active), len(active_union)
        ),
        "target_size_pearson_all": safe_correlation(reference_sizes, comparison_sizes, "pearson"),
        "target_size_spearman_all": safe_correlation(reference_sizes, comparison_sizes, "spearman"),
        "common_edge_mi_n": len(common),
        "common_edge_mi_pearson": safe_correlation(reference_mi, comparison_mi, "pearson"),
        "common_edge_mi_spearman": safe_correlation(reference_mi, comparison_mi, "spearman"),
        "common_edge_support_pearson": safe_correlation(
            reference_support, comparison_support, "pearson"
        ),
        "common_edge_support_spearman": safe_correlation(
            reference_support, comparison_support, "spearman"
        ),
    }


def locate_anchor(points: list[dict[str, Any]], requested: float) -> dict[str, Any]:
    matches = [
        point
        for point in points
        if math.isclose(point["p_value"], requested, rel_tol=1e-12, abs_tol=0.0)
    ]
    if len(matches) != 1:
        available = ", ".join(f"{point['p_value']:.12g}" for point in points)
        raise ValueError(f"Anchor p={requested:g} is absent or ambiguous; available: {available}")
    return matches[0]


def locate_pr66_stage(work_root: Path) -> Path:
    results_root = work_root / "results"
    candidates: list[Path] = []
    if all((results_root / driver).is_dir() for driver in DRIVERS):
        candidates.append(results_root)
    if results_root.is_dir():
        candidates.extend(
            path
            for path in sorted(results_root.iterdir())
            if path.is_dir()
            and path.name.lower().startswith("pr66")
            and all((path / driver).is_dir() for driver in DRIVERS)
        )
    unique = list(dict.fromkeys(candidates))
    if len(unique) != 1:
        raise ValueError(
            f"Could not identify exactly one PR66 results stage under {work_root}; "
            f"candidates={unique}"
        )
    return unique[0]


def validate_anchor_equivalence_evidence(
    work_root: Path,
    prior_work_root: Path,
    points: list[dict[str, Any]],
    *,
    sweep_design_sha256: str,
) -> dict[str, Any]:
    validation_root = work_root / "results" / "validation"
    table_path = validation_root / "anchor_seed_equivalence.tsv"
    manifest_path = validation_root / "anchor_seed_equivalence_manifest.json"
    manifest = read_json(manifest_path)
    table = pd.read_csv(table_path, sep="\t", dtype=str)
    expected_columns = (
        "sweep_point",
        "prior_stage",
        "driver",
        "seed",
        "data_sha256",
        "edges",
        "source_rows",
        "sweep_metadata_sha256",
        "prior_metadata_sha256",
    )
    if tuple(table.columns) != expected_columns or len(table) != 400:
        raise ValueError(
            f"Anchor-equivalence table must have the exact 400-row schema: {table_path}"
        )
    point_by_key = {point["p_key"]: point for point in points}
    anchor_stages = {
        "p1e-07": "pr67_7633ebb",
        "p_pr66_cutoff_match": "pr66_5809183",
    }
    expected_combinations = {
        (point, stage, driver, seed)
        for point, stage in anchor_stages.items()
        for driver in DRIVERS
        for seed in EXPECTED_SEEDS
    }
    observed_combinations: set[tuple[str, str, str, int]] = set()
    for row in table.itertuples(index=False):
        try:
            seed = int(row.seed)
            edges = int(row.edges)
            source_rows = int(row.source_rows)
        except (TypeError, ValueError) as error:
            raise ValueError(f"Invalid numeric anchor evidence row in {table_path}") from error
        combination = (row.sweep_point, row.prior_stage, row.driver, seed)
        if combination in observed_combinations:
            raise ValueError(f"Duplicate anchor-equivalence row: {combination}")
        observed_combinations.add(combination)
        if combination not in expected_combinations or edges < 0 or source_rows < 0:
            raise ValueError(f"Unexpected anchor-equivalence row: {combination}")
        for field in (
            "data_sha256",
            "sweep_metadata_sha256",
            "prior_metadata_sha256",
        ):
            require_sha256(getattr(row, field), field=field, source=table_path)
        stem = f"TF_run_{seed:03d}.json"
        sweep_metadata_path = (
            work_root
            / "results"
            / row.sweep_point
            / row.driver
            / "seed_metadata"
            / stem
        )
        prior_metadata_path = (
            prior_work_root
            / "results"
            / row.prior_stage
            / row.driver
            / "seed_metadata"
            / stem
        )
        if (
            sha256_file(sweep_metadata_path) != row.sweep_metadata_sha256
            or sha256_file(prior_metadata_path) != row.prior_metadata_sha256
        ):
            raise ValueError(
                f"Anchor metadata hash changed for {row.sweep_point}/{row.driver}/{seed}"
            )
        sweep_metadata = read_json(sweep_metadata_path)
        prior_metadata = read_json(prior_metadata_path)
        sweep_point = sweep_metadata.get("point")
        if (
            not isinstance(sweep_point, dict)
            or sweep_point.get("key") != row.sweep_point
            or sweep_metadata.get("driver") != row.driver
            or int(sweep_metadata.get("seed", -1)) != seed
            or prior_metadata.get("stage") != row.prior_stage
            or prior_metadata.get("driver") != row.driver
            or int(prior_metadata.get("seed", -1)) != seed
            or int(sweep_metadata.get("adjacency", {}).get("edges", -1)) != edges
            or int(prior_metadata.get("adjacency", {}).get("edges", -1)) != edges
            or int(sweep_metadata.get("adjacency", {}).get("source_rows", -1))
            != source_rows
            or int(prior_metadata.get("adjacency", {}).get("source_rows", -1))
            != source_rows
            or sweep_metadata.get("adjacency", {}).get("data_sha256")
            != row.data_sha256
            or prior_metadata.get("adjacency", {}).get("data_sha256")
            != row.data_sha256
        ):
            raise ValueError(
                f"Anchor data-section hash changed for "
                f"{row.sweep_point}/{row.driver}/{seed}"
            )
    if observed_combinations != expected_combinations:
        raise ValueError(f"Anchor-equivalence rows do not cover both anchors: {table_path}")

    expected_anchors = [
        {
            "sweep_point": point_key,
            "prior_stage": prior_stage,
            "expected_p": point_by_key[point_key]["p_value"],
            "expected_mi_cutoff": point_by_key[point_key]["mi_cutoff"],
            "point_manifest_sha256": point_by_key[point_key]["manifest_sha256"],
        }
        for point_key, prior_stage in anchor_stages.items()
    ]
    validator = (
        REPO_ROOT
        / "benchmarks"
        / "brca100_pr67_threshold_sweep"
        / "validate_anchor_equivalence.py"
    )
    expected_manifest = {
        "schema": "sjaracne-brca100-pr67-p-sweep-anchor-equivalence-v1",
        "anchors": expected_anchors,
        "drivers": list(DRIVERS),
        "seeds": list(EXPECTED_SEEDS),
        "comparisons": 400,
        "all_data_sections_equal": True,
        "table": str(table_path),
        "table_sha256": sha256_file(table_path),
        "script_sha256": sha256_file(validator),
        "sweep_design_sha256": sweep_design_sha256,
        "sweep_work_root": str(work_root),
        "prior_work_root": str(prior_work_root),
    }
    if manifest != expected_manifest:
        raise ValueError(f"Anchor-equivalence manifest mismatch: {manifest_path}")
    return {
        "manifest": str(manifest_path),
        "manifest_sha256": sha256_file(manifest_path),
        "table": str(table_path),
        "table_sha256": sha256_file(table_path),
        "comparisons": 400,
        "all_data_sections_equal": True,
    }


def save_plot(figure: plt.Figure, plots_root: Path, stem: str) -> None:
    figure.savefig(plots_root / f"{stem}.png", dpi=180)
    figure.savefig(plots_root / f"{stem}.svg", metadata=SVG_METADATA)


def plot_core_metrics(
    summary: pd.DataFrame,
    points: list[dict[str, Any]],
    anchor: dict[str, Any],
    plots_root: Path,
    context_summary: pd.DataFrame | None,
) -> None:
    metric_list = [
        ("consensus_edges", "Consensus edges", "log"),
        ("active_driver_fraction", "Active candidate drivers", "fraction"),
        (
            "incident_node_fraction_expression",
            "Incident nodes / expression nodes",
            "fraction",
        ),
        (
            "largest_weak_component_fraction_incident",
            "Largest weak component / incident nodes",
            "fraction",
        ),
        (
            "target_size_median_zero_filled",
            "Median targets / candidate (zeros included)",
            "symlog",
        ),
        ("support_fraction_median", "Median consensus support", "fraction"),
        ("consensus_MI_median", "Median consensus MI", "linear"),
    ]
    if "netbid2_scale_free_adjusted_r2" in summary and summary[
        "netbid2_scale_free_adjusted_r2"
    ].notna().any():
        metric_list.append(
            (
                "netbid2_scale_free_adjusted_r2",
                "NetBID2 scale-free adjusted R2",
                "linear",
            )
        )
    metrics = tuple(metric_list)
    x_values = np.asarray([point["log10_p"] for point in points])
    x_labels = [compact_plot_p_label(point["p_value"]) for point in points]
    figure, axes = plt.subplots(
        len(DRIVERS),
        len(metrics),
        figsize=(3.85 * len(metrics), 8.4),
        constrained_layout=True,
        squeeze=False,
    )
    colors = {"tf": "#0072B2", "sig": "#D55E00"}
    for row_index, driver in enumerate(DRIVERS):
        subset = summary[summary["driver_class"] == driver].sort_values("p_value")
        if len(subset) != len(points):
            raise ValueError(f"Incomplete summary rows for {driver}")
        for column_index, (metric, label, scale) in enumerate(metrics):
            axis = axes[row_index, column_index]
            y_values = subset[metric].to_numpy(dtype=float)
            axis.plot(
                x_values,
                y_values,
                marker="o",
                markersize=4.5,
                linewidth=1.8,
                color=colors[driver],
                label="PR67 sweep",
            )
            if metric in ("support_fraction_median", "consensus_MI_median"):
                prefix = "support_fraction" if metric.startswith("support") else "consensus_MI"
                axis.fill_between(
                    x_values,
                    subset[f"{prefix}_q25"].to_numpy(dtype=float),
                    subset[f"{prefix}_q75"].to_numpy(dtype=float),
                    color=colors[driver],
                    alpha=0.14,
                    linewidth=0.0,
                    label="IQR",
                )
            if context_summary is not None and metric in context_summary:
                context_value = float(
                    context_summary.loc[context_summary["driver_class"] == driver, metric].iloc[0]
                )
                if math.isfinite(context_value):
                    axis.axhline(
                        context_value,
                        color="#555555",
                        linestyle="--",
                        linewidth=1.2,
                        label="PR66 context",
                    )
            axis.axvline(anchor["log10_p"], color="#777777", linestyle=":", linewidth=1.0)
            if scale == "symlog":
                axis.set_yscale("symlog", linthresh=1)
                axis.set_ylim(bottom=0)
            elif scale == "log":
                axis.set_yscale("log")
            elif scale == "fraction":
                axis.set_ylim(-0.025, 1.025)
            title = label.replace(" / ", "\n/ ")
            if title == "NetBID2 scale-free adjusted R2":
                title = "NetBID2 scale-free\nadjusted R2"
            axis.set_title(title, fontsize=10)
            axis.grid(alpha=0.22)
            if row_index == len(DRIVERS) - 1:
                axis.set_xticks(x_values)
                axis.set_xticklabels(
                    x_labels,
                    rotation=30,
                    ha="right",
                    fontsize=6.5,
                )
                # After 1e-5 the log-spaced positions are too close for one
                # label row. Four deterministic tiers retain every point,
                # including the near-adjacent 3e-4 and cutoff-match points.
                for label_index, tick_label in enumerate(
                    axis.get_xticklabels()[2:], start=2
                ):
                    tier = (label_index - 2) % 4
                    tick_label.set_y(-0.06 * tier)
                axis.set_xlabel("Per-subsample p (x position = log10(p))")
            else:
                axis.set_xticks(x_values, [])
            if column_index == 0:
                axis.set_ylabel(driver.upper())
            if row_index == 0 and column_index == 0:
                axis.legend(frameon=False, fontsize=8)
    figure.suptitle(
        "BRCA100 PR67 per-subsample threshold sweep; dotted line is the "
        f"p={anchor['p_value']:.3g} anchor"
    )
    save_plot(figure, plots_root, "core_metrics_vs_log10_p")
    plt.close(figure)


def plot_overlaps(
    adjacent: pd.DataFrame,
    anchor_overlap: pd.DataFrame,
    plots_root: Path,
    anchor_p: float,
) -> None:
    figure, axes = plt.subplots(1, len(DRIVERS), figsize=(12.5, 4.8), constrained_layout=True)
    for axis, driver in zip(axes, DRIVERS):
        anchored = anchor_overlap[anchor_overlap["driver_class"] == driver].sort_values(
            "comparison_p_value"
        )
        neighboring = adjacent[adjacent["driver_class"] == driver].sort_values(
            "comparison_p_value"
        )
        axis.plot(
            anchored["comparison_log10_p"],
            anchored["edge_jaccard"],
            marker="o",
            color="#0072B2",
            label=f"vs p={anchor_p:.3g} anchor",
        )
        axis.plot(
            neighboring["comparison_log10_p"],
            neighboring["edge_jaccard"],
            marker="s",
            color="#D55E00",
            label="vs adjacent stricter p",
        )
        axis.set_ylim(-0.025, 1.025)
        axis.set_xlabel("log10(per-subsample p)")
        axis.set_ylabel("Directed-edge Jaccard")
        axis.set_title(driver.upper())
        axis.grid(alpha=0.22)
    axes[0].legend(frameon=False)
    figure.suptitle("Consensus-network edge stability across the PR67 threshold sweep")
    save_plot(figure, plots_root, "edge_overlap_vs_log10_p")
    plt.close(figure)


def integrate_optional_netbid2(
    summary: pd.DataFrame,
    points: list[dict[str, Any]],
    arms: dict[tuple[str, str], dict[str, Any]],
    *,
    sweep_design_sha256: str,
) -> tuple[pd.DataFrame, dict[tuple[str, str], dict[str, Any]]]:
    expected_keys = [(point["p_key"], driver) for point in points for driver in DRIVERS]
    paths = {
        key: arms[key]["arm_root"] / "netbid2_qc" / "network_summary.tsv"
        for key in expected_keys
    }
    present = {key for key, path in paths.items() if path.is_file()}
    netbid_evidence = any(
        arms[key]["arm_root"].joinpath("netbid2_qc").exists()
        or arms[key]["arm_root"].joinpath("netbid2_qc_manifest.json").exists()
        for key in expected_keys
    ) or points[0]["root"].parent.joinpath("netbid2_qc_manifest.json").exists()
    result = summary.copy()
    result["netbid2_scale_free_adjusted_r2"] = math.nan
    if not present and not netbid_evidence:
        return result, arms
    if present != set(expected_keys):
        raise ValueError(
            "NetBID2 QC is partially present; either provide all 26 arm summaries or none. "
            f"Missing={sorted(set(expected_keys) - present)}"
        )
    aggregate_path = points[0]["root"].parent / "netbid2_qc_manifest.json"
    aggregate = read_json(aggregate_path)
    environment = aggregate.get("environment")
    if not isinstance(environment, dict) or not all(
        isinstance(key, str) and isinstance(value, str)
        for key, value in environment.items()
    ):
        raise ValueError(f"Invalid NetBID2 environment in {aggregate_path}")
    expected_point_order = [point["p_key"] for point in points]
    expected_aggregate = {
        "schema": (
            "sjaracne-brca100-pr67-p-sweep-netbid2-summary-aggregate-v1"
        ),
        "sweep_design_sha256": sweep_design_sha256,
        "all_sweep_points": expected_point_order,
    }
    for field, expected in expected_aggregate.items():
        if aggregate.get(field) != expected:
            raise ValueError(f"NetBID2 aggregate {field} mismatch in {aggregate_path}")
    selection = aggregate.get("selection")
    if (
        not isinstance(selection, dict)
        or set(selection) != {"points", "drivers"}
        or selection.get("points") != expected_point_order
        or selection.get("drivers") != list(DRIVERS)
    ):
        raise ValueError(f"NetBID2 aggregate selection is not a complete summary run")
    expected_aggregate_fields = {
        "schema",
        "environment",
        "sweep_design_sha256",
        "all_sweep_points",
        "selection",
        "summary_runs",
        "fingerprint",
    }
    if set(aggregate) != expected_aggregate_fields:
        raise ValueError(
            f"Unexpected fields in immutable NetBID2 summary aggregate: {aggregate_path}"
        )
    aggregate_without_fingerprint = dict(aggregate)
    observed_aggregate_fingerprint = aggregate_without_fingerprint.pop(
        "fingerprint", None
    )
    if observed_aggregate_fingerprint != json_fingerprint(
        aggregate_without_fingerprint
    ):
        raise ValueError(f"NetBID2 aggregate fingerprint mismatch: {aggregate_path}")
    summary_runs = aggregate.get("summary_runs")
    if not isinstance(summary_runs, list) or len(summary_runs) != len(expected_keys):
        raise ValueError("NetBID2 aggregate does not contain all 26 summary runs")
    aggregate_summary_by_arm: dict[tuple[str, str], dict[str, Any]] = {}
    for record in summary_runs:
        if not isinstance(record, dict):
            raise ValueError(f"Invalid NetBID2 summary record in {aggregate_path}")
        key = (str(record.get("point")), str(record.get("driver")))
        if key in aggregate_summary_by_arm:
            raise ValueError(f"Duplicate NetBID2 summary record for {key}")
        aggregate_summary_by_arm[key] = record
    if set(aggregate_summary_by_arm) != set(expected_keys):
        raise ValueError("NetBID2 aggregate summary arms do not match the sweep")
    wrapper = REPO_ROOT / "benchmarks" / "brca100_netbid_qc" / "netbid2-r"
    r_script = REPO_ROOT / "benchmarks" / "brca100_pr67_threshold_sweep" / "run_netbid_qc.R"
    for path in (wrapper, r_script):
        if not path.is_file():
            raise ValueError(f"Missing NetBID2 provenance input: {path}")
    wrapper_hash = sha256_file(wrapper)
    r_script_hash = sha256_file(r_script)
    point_by_key = {point["p_key"]: point for point in points}
    metric_map = {
        "candidate_drivers": "candidate_drivers",
        "active_drivers": "active_drivers",
        "active_driver_fraction": "active_driver_fraction",
        "edges": "consensus_edges",
        "incident_nodes": "incident_nodes",
        "weak_components": "weak_components",
        "largest_weak_component": "largest_weak_component_nodes",
        "largest_weak_component_fraction": (
            "largest_weak_component_fraction_incident"
        ),
        "target_size_zero_mean": "target_size_mean_zero_filled",
        "target_size_zero_median": "target_size_median_zero_filled",
        "target_size_zero_q25": "target_size_q25_zero_filled",
        "target_size_zero_q75": "target_size_q75_zero_filled",
        "target_size_zero_max": "target_size_max",
    }
    for point_key, driver in expected_keys:
        path = paths[(point_key, driver)]
        point = point_by_key[point_key]
        arm_root = arms[(point_key, driver)]["arm_root"]
        manifest_path = arm_root / "netbid2_qc_manifest.json"
        pending_path = arm_root / "netbid2_qc_manifest.pending.json"
        partial_path = arm_root / "netbid2_qc.partial"
        if pending_path.exists() or partial_path.exists():
            raise ValueError(f"Unresolved NetBID2 state under {arm_root}")
        manifest = read_json(manifest_path)
        fingerprint_payload = {
            "schema": "sjaracne-brca100-pr67-p-sweep-netbid2-v1",
            "mode": "summary",
            "point": point_key,
            "p_value": point["p_value"],
            "mi_cutoff": point["mi_cutoff"],
            "point_manifest_sha256": arms[(point_key, driver)][
                "point_manifest_sha256"
            ],
            "sweep_design_sha256": sweep_design_sha256,
            "driver": driver,
            "driver_sha256": EXPECTED_INPUT_SHA256[CANDIDATE_FILES[driver]],
            "consensus_sha256": arms[(point_key, driver)]["consensus_sha256"],
            "consensus_manifest_sha256": arms[(point_key, driver)][
                "consensus_manifest_sha256"
            ],
            "r_script_sha256": r_script_hash,
            "wrapper_sha256": wrapper_hash,
            "environment": environment,
            "prefix": "TF_" if driver == "tf" else "SIG_",
        }
        for field, expected in fingerprint_payload.items():
            if manifest.get(field) != expected:
                raise ValueError(
                    f"NetBID2 manifest {field} mismatch for {point_key}/{driver}"
                )
        if manifest.get("fingerprint") != json_fingerprint(fingerprint_payload):
            raise ValueError(f"NetBID2 fingerprint mismatch for {point_key}/{driver}")
        netbid_root = path.parent
        inventory = [
            {
                "path": artifact.relative_to(netbid_root).as_posix(),
                "bytes": artifact.stat().st_size,
                "sha256": sha256_file(artifact),
            }
            for artifact in sorted(
                item for item in netbid_root.rglob("*") if item.is_file()
            )
        ]
        stdout_path = arm_root / "logs" / "netbid2_qc.stdout.log"
        stderr_path = arm_root / "logs" / "netbid2_qc.stderr.log"
        if (
            manifest.get("output_inventory") != inventory
            or not stdout_path.is_file()
            or not stderr_path.is_file()
            or manifest.get("stdout_sha256") != sha256_file(stdout_path)
            or manifest.get("stderr_sha256") != sha256_file(stderr_path)
            or int(manifest.get("stderr_bytes", -1)) != stderr_path.stat().st_size
        ):
            raise ValueError(f"NetBID2 output/log inventory mismatch for {point_key}/{driver}")
        if aggregate_summary_by_arm[(point_key, driver)] != manifest:
            raise ValueError(
                f"Root and per-arm NetBID2 manifests disagree for {point_key}/{driver}"
            )
        table = pd.read_csv(path, sep="\t", dtype=str)
        if set(table.columns) != {"metric", "value"} or table["metric"].duplicated().any():
            raise ValueError(f"Unexpected NetBID2 network summary: {path}")
        values: dict[str, float] = {}
        for row in table.itertuples(index=False):
            values[str(row.metric)] = (
                math.nan if str(row.value) == "NA" else float(row.value)
            )
        missing = (set(metric_map) | {"scale_free_adjusted_r2"}) - set(values)
        if missing:
            raise ValueError(f"Missing NetBID2 metrics in {path}: {sorted(missing)}")
        selector = (result["p_key"] == point_key) & (result["driver_class"] == driver)
        if int(selector.sum()) != 1:
            raise ValueError(f"Cannot align NetBID2 summary for {point_key}/{driver}")
        row = result.loc[selector].iloc[0]
        for netbid_metric, python_metric in metric_map.items():
            observed = values[netbid_metric]
            expected = float(row[python_metric])
            if not math.isclose(observed, expected, rel_tol=1e-10, abs_tol=1e-10):
                raise ValueError(
                    f"NetBID2 {netbid_metric}={observed} disagrees with "
                    f"edge-derived {python_metric}={expected} in {path}"
                )
        result.loc[selector, "netbid2_scale_free_adjusted_r2"] = values[
            "scale_free_adjusted_r2"
        ]
        netbid_paths = sorted(item for item in netbid_root.rglob("*") if item.is_file())
        arms[(point_key, driver)]["netbid2_qc_set_sha256"] = file_set_sha256(
            netbid_paths, netbid_root
        )
        arms[(point_key, driver)]["netbid2_network_summary_sha256"] = sha256_file(path)
        arms[(point_key, driver)]["netbid2_manifest_sha256"] = sha256_file(
            manifest_path
        )
        arms[(point_key, driver)]["netbid2_fingerprint"] = manifest["fingerprint"]
    return result, arms


def operating_point_screen(
    summary: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in summary.sort_values(["p_value", "driver_class"]).to_dict(
        orient="records"
    ):
        p_value = float(record["p_value"])
        within_held_out = HELD_OUT_P_MIN <= p_value <= HELD_OUT_P_MAX
        calibration_class = calibration_point_class(p_value)
        exact_grid = calibration_class == "exact-gate2-grid"
        if calibration_class != record["calibration_point_class"]:
            raise ValueError(
                f"Calibration class changed for {record['p_key']}: "
                f"{record['calibration_point_class']} != {calibration_class}"
            )
        coverage_pass = (
            float(record["active_driver_fraction"])
            >= ENGINEERING_FLOORS["active_driver_fraction"]
        )
        largest_component_pass = (
            float(record["largest_weak_component_fraction_incident"])
            >= ENGINEERING_FLOORS["largest_weak_component_fraction_incident"]
        )
        incident_pass = (
            float(record["incident_node_fraction_expression"])
            >= ENGINEERING_FLOORS["incident_node_fraction_expression"]
        )
        rows.append(
            {
                "p_key": record["p_key"],
                "p_value": p_value,
                "p_label": record["p_label"],
                "driver_class": record["driver_class"],
                "calibration_point_class": calibration_class,
                "within_held_out_range": within_held_out,
                "exact_gate2_grid_point": exact_grid,
                "active_driver_fraction": record["active_driver_fraction"],
                "active_driver_fraction_floor": ENGINEERING_FLOORS[
                    "active_driver_fraction"
                ],
                "active_driver_fraction_pass": coverage_pass,
                "largest_weak_component_fraction_incident": record[
                    "largest_weak_component_fraction_incident"
                ],
                "largest_weak_component_fraction_incident_floor": (
                    ENGINEERING_FLOORS[
                        "largest_weak_component_fraction_incident"
                    ]
                ),
                "largest_weak_component_fraction_incident_pass": (
                    largest_component_pass
                ),
                "incident_node_fraction_expression": record[
                    "incident_node_fraction_expression"
                ],
                "incident_node_fraction_expression_floor": ENGINEERING_FLOORS[
                    "incident_node_fraction_expression"
                ],
                "incident_node_fraction_expression_pass": incident_pass,
                "driver_network_floor_pass": (
                    coverage_pass and largest_component_pass and incident_pass
                ),
                "candidate_pair_tests": record["candidate_pair_tests"],
                "candidate_pair_tests_interpretation": record[
                    "candidate_pair_tests_interpretation"
                ],
                "nominal_null_exceedances_before_DPI": record[
                    "nominal_null_exceedances_before_DPI"
                ],
                "nominal_null_exceedances_interpretation": record[
                    "nominal_null_exceedances_interpretation"
                ],
            }
        )
    screen = pd.DataFrame(rows)
    joint_pass_by_key: dict[str, bool] = {}
    for p_key, group in screen.groupby("p_key", sort=False):
        joint_pass_by_key[p_key] = (
            set(group["driver_class"]) == set(DRIVERS)
            and bool(group["driver_network_floor_pass"].all())
            and bool(group["within_held_out_range"].all())
        )
    screen["joint_tf_sig_held_out_pass"] = screen["p_key"].map(joint_pass_by_key)

    candidates = (
        screen.loc[screen["joint_tf_sig_held_out_pass"]]
        .drop_duplicates(["p_key", "p_value", "calibration_point_class"])
        .sort_values("p_value")
    )
    selected: pd.Series | None = None
    selection_status = "no-passing-point"
    if not candidates.empty:
        selected = candidates.iloc[0]
        selection_status = (
            "selected-exact-gate2-grid-point"
            if bool(selected["exact_gate2_grid_point"])
            else "provisional-interpolated-operating-point"
        )
    selection = {
        "schema": "sjaracne-brca100-pr67-threshold-selection-v1",
        "selection_status": selection_status,
        "selected_p_key": None if selected is None else str(selected["p_key"]),
        "selected_p_value": None if selected is None else float(selected["p_value"]),
        "selected_calibration_point_class": (
            None if selected is None else str(selected["calibration_point_class"])
        ),
        "selected_point_kind": (
            None if selected is None else "provisional-topology-operating-point"
        ),
        "rule": (
            "Among points passing all three engineering floors for both TF and SIG "
            "inside the held-out range, choose the smallest p. Exact Gate-2 grid "
            "status versus within-range interpolation is reported but does not permit "
            "skipping a stricter passing point."
        ),
        "engineering_floors": ENGINEERING_FLOORS,
        "held_out_p_range": [HELD_OUT_P_MIN, HELD_OUT_P_MAX],
        "exact_gate2_grid": sorted(EXACT_GATE2_GRID),
        "scope": (
            "Provisional topology operating-point screen only; not biological "
            "validation, not an empirical FDR estimate, and not a regulon-density "
            "selection rule."
        ),
        "declaration_timing": (
            "Engineering floors were declared after observing the prior PR66/PR67 "
            "endpoint comparison but before observing intermediate sweep-point results."
        ),
    }
    return screen, selection


def plot_coverage_vs_null_burden(summary: pd.DataFrame, plots_root: Path) -> None:
    figure, axes = plt.subplots(1, len(DRIVERS), figsize=(13.0, 4.9), constrained_layout=True)
    metrics = (
        ("active_driver_fraction", "Active-driver fraction", "#0072B2", "o"),
        (
            "incident_node_fraction_expression",
            "Incident-node fraction",
            "#D55E00",
            "s",
        ),
        (
            "largest_weak_component_fraction_incident",
            "Largest-component fraction",
            "#009E73",
            "^",
        ),
    )
    for axis, driver in zip(axes, DRIVERS):
        subset = summary[summary["driver_class"] == driver].sort_values("p_value")
        x_values = subset["nominal_null_exceedances_before_DPI"].to_numpy(dtype=float)
        for metric, label, color, marker in metrics:
            axis.plot(
                x_values,
                subset[metric].to_numpy(dtype=float),
                color=color,
                marker=marker,
                linewidth=1.6,
                markersize=4.5,
                label=label,
            )
        annotations = zip(
            x_values,
            subset["active_driver_fraction"],
            subset["p_value"],
        )
        loose_label_start = max(0, len(subset) - 7)
        for point_index, (x_value, y_value, p_value) in enumerate(annotations):
            if point_index >= loose_label_start:
                # The seven loosest thresholds converge near y=1 and have
                # closely spaced log-x positions. Stack their labels below the
                # markers in deterministic tiers.
                tier = point_index - loose_label_start
                offset = (-3, -(8 + 10 * tier))
                horizontal_alignment = "right"
                vertical_alignment = "top"
            else:
                offset = (3, 3)
                horizontal_alignment = "left"
                vertical_alignment = "bottom"
            axis.annotate(
                compact_plot_p_label(float(p_value)),
                (x_value, y_value),
                xytext=offset,
                textcoords="offset points",
                fontsize=7,
                ha=horizontal_alignment,
                va=vertical_alignment,
            )
        axis.set_xscale("log")
        axis.set_ylim(-0.025, 1.025)
        axis.set_xlabel(
            "Nominal null exceedances per seed before DPI (model-based proxy)"
        )
        axis.set_ylabel("Coverage/connectivity fraction")
        axis.set_title(driver.upper())
        axis.grid(alpha=0.22)
    axes[0].legend(frameon=False, fontsize=8)
    figure.suptitle(
        "Coverage versus per-seed modeled null burden; proxy is not empirical FDR"
    )
    save_plot(figure, plots_root, "coverage_vs_nominal_null_burden")
    plt.close(figure)


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, float_format="%.12g")


def write_json(path: Path, value: object) -> None:
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    temporary.replace(path)


def software_versions() -> dict[str, str]:
    packages: dict[str, str] = {}
    for package in ("matplotlib", "numpy", "pandas", "scipy"):
        try:
            packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            packages[package] = "not-installed"
    return {
        "python": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        **packages,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--work-root",
        type=Path,
        required=True,
        help="Sweep root containing inputs/ and results/<p_key>/point_manifest.json",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        help="Output directory (default: <work-root>/results/analysis)",
    )
    parser.add_argument(
        "--anchor-p",
        type=float,
        default=1e-7,
        help="Sweep p used for anchor-overlap comparisons (default: 1e-7)",
    )
    parser.add_argument(
        "--pr66-work-root",
        type=Path,
        help="Optional prior matched BRCA100 work root used only as plot/table context",
    )
    args = parser.parse_args()

    work_root = args.work_root.resolve()
    results_root = work_root / "results"
    if not results_root.is_dir():
        raise ValueError(f"Missing results directory: {results_root}")
    final_output_root = (
        args.output_root.resolve()
        if args.output_root is not None
        else results_root / "analysis"
    )
    staging_output_root = final_output_root.with_name(final_output_root.name + ".partial")
    optional_context_outputs = (
        final_output_root / "pr66_context_summary.tsv",
        final_output_root / "pr66_context_overlap.tsv",
    )
    if args.pr66_work_root is None:
        stale = [path for path in optional_context_outputs if path.exists()]
        if stale:
            raise ValueError(
                "Stale PR66 context outputs exist but --pr66-work-root was omitted: "
                f"{stale}. Remove those exact files or rerun with the context argument."
            )
    if staging_output_root.exists():
        raise ValueError(
            f"Stale partial analysis directory exists: {staging_output_root}. "
            "Inspect and remove that exact directory before retrying."
        )
    if final_output_root.exists():
        raise ValueError(
            f"Analysis output already exists: {final_output_root}. To avoid mixing "
            "runs, remove that exact directory or choose a new --output-root."
        )

    design, design_points, sweep_design_hash, build = validate_sweep_design(work_root)
    points = point_metadata(results_root, design, design_points)
    anchor = locate_anchor(points, args.anchor_p)
    candidates = candidate_lists(work_root / "inputs")
    expression_path = work_root / "inputs" / "BRCA100.exp"
    annotations = expression_annotations(expression_path)
    expression_ids = set(annotations)
    pair_tests: dict[str, int] = {}
    for driver in DRIVERS:
        pair_tests[driver] = eligible_candidate_pair_tests(
            candidates[driver], annotations
        )
        if pair_tests[driver] != EXPECTED_CANDIDATE_PAIR_TESTS[driver]:
            raise ValueError(
                f"Exact {driver} candidate-pair test count changed: "
                f"observed={pair_tests[driver]}, "
                f"expected={EXPECTED_CANDIDATE_PAIR_TESTS[driver]}"
            )
    seed_sets, seed_summary, arms = validate_completed_seeds(
        work_root, points, design_points, design, build
    )
    arms = validate_arm_artifacts(
        work_root,
        points,
        arms,
        sweep_design_sha256=sweep_design_hash,
    )

    point_table = pd.DataFrame(
        [
            {
                key: point[key]
                for key in (
                    "p_key",
                    "p_value",
                    "log10_p",
                    "p_token",
                    "p_label",
                    "role",
                    "mi_cutoff",
                    "validation_class",
                    "calibration_point_class",
                    "tail_extrapolated",
                    "commit",
                    "model_sha256",
                    "manifest_sha256",
                )
            }
            for point in points
        ]
    )

    networks: dict[tuple[str, str], pd.DataFrame] = {}
    sizes: dict[tuple[str, str], pd.Series] = {}
    summary_rows: list[dict[str, Any]] = []
    for point in points:
        for driver in DRIVERS:
            key = (point["p_key"], driver)
            frame = merge_network_support(point["root"] / driver)
            networks[key] = frame
            row, target_sizes = summarize_network(
                frame,
                candidates[driver],
                expression_ids,
                prefix={
                    "p_key": point["p_key"],
                    "p_value": point["p_value"],
                    "log10_p": point["log10_p"],
                    "p_token": point["p_token"],
                    "p_label": point["p_label"],
                    "role": point["role"],
                    "mi_cutoff": point["mi_cutoff"],
                    "validation_class": point["validation_class"],
                    "calibration_point_class": point[
                        "calibration_point_class"
                    ],
                    "tail_extrapolated": point["tail_extrapolated"],
                    "driver_class": driver,
                },
                seed_runs=len(seed_sets[key]),
                candidate_pair_tests=pair_tests[driver],
            )
            summary_rows.append(row)
            sizes[key] = target_sizes
            provenance = arms[key]
            manifest_expectations = {
                "consensus_edges_manifest": len(frame),
                "consensus_active_drivers_manifest": row["active_drivers"],
                "consensus_incident_nodes_manifest": row["incident_nodes"],
            }
            for field, observed in manifest_expectations.items():
                if provenance[field] != observed:
                    raise ValueError(
                        f"Consensus manifest {field}={provenance[field]} disagrees "
                        f"with parsed value {observed} for {key}"
                    )
    summary = pd.DataFrame(summary_rows).sort_values(["driver_class", "p_value"])
    summary, arms = integrate_optional_netbid2(
        summary,
        points,
        arms,
        sweep_design_sha256=sweep_design_hash,
    )
    screen, selection = operating_point_screen(summary)

    adjacent_rows: list[dict[str, Any]] = []
    for reference_point, comparison_point in zip(points[:-1], points[1:]):
        for driver in DRIVERS:
            reference_key = (reference_point["p_key"], driver)
            comparison_key = (comparison_point["p_key"], driver)
            adjacent_rows.append(
                overlap_row(
                    networks[reference_key],
                    networks[comparison_key],
                    sizes[reference_key],
                    sizes[comparison_key],
                    prefix={
                        "driver_class": driver,
                        "reference_p_key": reference_point["p_key"],
                        "reference_p_value": reference_point["p_value"],
                        "reference_log10_p": reference_point["log10_p"],
                        "comparison_p_key": comparison_point["p_key"],
                        "comparison_p_value": comparison_point["p_value"],
                        "comparison_log10_p": comparison_point["log10_p"],
                    },
                )
            )
    adjacent = pd.DataFrame(adjacent_rows)

    anchor_rows: list[dict[str, Any]] = []
    for comparison_point in points:
        for driver in DRIVERS:
            reference_key = (anchor["p_key"], driver)
            comparison_key = (comparison_point["p_key"], driver)
            anchor_rows.append(
                overlap_row(
                    networks[reference_key],
                    networks[comparison_key],
                    sizes[reference_key],
                    sizes[comparison_key],
                    prefix={
                        "driver_class": driver,
                        "reference_p_key": anchor["p_key"],
                        "reference_p_value": anchor["p_value"],
                        "reference_log10_p": anchor["log10_p"],
                        "comparison_p_key": comparison_point["p_key"],
                        "comparison_p_value": comparison_point["p_value"],
                        "comparison_log10_p": comparison_point["log10_p"],
                    },
                )
            )
    anchor_overlap = pd.DataFrame(anchor_rows)

    context_summary: pd.DataFrame | None = None
    context_provenance: dict[str, Any] | None = None
    context_overlap: pd.DataFrame | None = None
    if args.pr66_work_root is not None:
        pr66_work_root = args.pr66_work_root.resolve()
        pr66_stage = locate_pr66_stage(pr66_work_root)
        prior_input_hashes: dict[str, str] = {}
        for filename, expected_hash in EXPECTED_INPUT_SHA256.items():
            prior_path = pr66_work_root / "inputs" / filename
            if not prior_path.is_file():
                raise ValueError(f"Missing required PR66 context input: {prior_path}")
            actual_hash = sha256_file(prior_path)
            if actual_hash != expected_hash:
                raise ValueError(
                    f"PR66 context input {filename} hash {actual_hash} does not match "
                    f"the pinned sweep hash {expected_hash}"
                )
            prior_input_hashes[filename] = actual_hash
        anchor_equivalence = validate_anchor_equivalence_evidence(
            work_root,
            pr66_work_root,
            points,
            sweep_design_sha256=sweep_design_hash,
        )
        context_networks: dict[str, pd.DataFrame] = {}
        context_sizes: dict[str, pd.Series] = {}
        context_rows: list[dict[str, Any]] = []
        for driver in DRIVERS:
            frame = merge_network_support(pr66_stage / driver)
            cutoff_match = networks[("p_pr66_cutoff_match", driver)]
            if not frame.equals(cutoff_match):
                common = len(set(frame.index) & set(cutoff_match.index))
                union = len(set(frame.index) | set(cutoff_match.index))
                raise ValueError(
                    "PR67 cutoff-match consensus is not exactly equal to the PR66 "
                    f"consensus for {driver}: edge_jaccard={safe_ratio(common, union)}"
                )
            context_networks[driver] = frame
            row, target_sizes = summarize_network(
                frame,
                candidates[driver],
                expression_ids,
                prefix={"context": "pr66", "driver_class": driver},
                seed_runs=0,
                candidate_pair_tests=pair_tests[driver],
            )
            context_rows.append(row)
            context_sizes[driver] = target_sizes
        context_summary = pd.DataFrame(context_rows)
        context_overlap_rows: list[dict[str, Any]] = []
        for comparison_point in points:
            for driver in DRIVERS:
                comparison_key = (comparison_point["p_key"], driver)
                context_overlap_rows.append(
                    overlap_row(
                        context_networks[driver],
                        networks[comparison_key],
                        context_sizes[driver],
                        sizes[comparison_key],
                        prefix={
                            "driver_class": driver,
                            "reference_context": "pr66",
                            "comparison_p_key": comparison_point["p_key"],
                            "comparison_p_value": comparison_point["p_value"],
                            "comparison_log10_p": comparison_point["log10_p"],
                        },
                    )
                )
        context_overlap = pd.DataFrame(context_overlap_rows)
        context_provenance = {
            "work_root": str(pr66_work_root),
            "stage_root": str(pr66_stage),
            "input_sha256": prior_input_hashes,
            "anchor_seed_equivalence": anchor_equivalence,
            "cutoff_match_consensus_exact": True,
            "arm_files": {
                driver: {
                    "network_sha256": sha256_file(
                        arm_file(pr66_stage / driver, NETWORK_FILENAME)
                    ),
                    "support_sha256": sha256_file(
                        arm_file(pr66_stage / driver, SUPPORT_FILENAME)
                    ),
                }
                for driver in DRIVERS
            },
        }

    # Do not emit partial analysis products until every provenance and content
    # check above has passed.
    output_root = staging_output_root
    plots_root = output_root / "plots"
    optional_context_outputs = (
        output_root / "pr66_context_summary.tsv",
        output_root / "pr66_context_overlap.tsv",
    )
    output_root.mkdir(parents=True, exist_ok=False)
    plots_root.mkdir(parents=True, exist_ok=True)
    write_tsv(seed_summary, output_root / "seed_manifest_summary.tsv")
    write_tsv(point_table, output_root / "point_manifest_summary.tsv")
    write_tsv(summary, output_root / "network_summary.tsv")
    write_tsv(adjacent, output_root / "adjacent_overlap.tsv")
    write_tsv(anchor_overlap, output_root / "anchor_overlap.tsv")
    write_tsv(screen, output_root / "operating_point_screen.tsv")
    write_json(output_root / "selection.json", selection)
    if context_summary is not None and context_overlap is not None:
        write_tsv(context_summary, optional_context_outputs[0])
        write_tsv(context_overlap, optional_context_outputs[1])

    arm_rows: list[dict[str, Any]] = []
    arm_hash_fields = (
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
    )
    for point in points:
        for driver in DRIVERS:
            provenance = arms[(point["p_key"], driver)]
            arm_rows.append(
                {
                    "p_key": point["p_key"],
                    "p_value": point["p_value"],
                    "driver_class": driver,
                    **{field: provenance.get(field, "") for field in arm_hash_fields},
                }
            )
    arm_provenance = pd.DataFrame(arm_rows)
    write_tsv(arm_provenance, output_root / "arm_provenance.tsv")

    plot_core_metrics(summary, points, anchor, plots_root, context_summary)
    plot_overlaps(adjacent, anchor_overlap, plots_root, anchor["p_value"])
    plot_coverage_vs_null_burden(summary, plots_root)

    expected_outputs = [
        output_root / "seed_manifest_summary.tsv",
        output_root / "point_manifest_summary.tsv",
        output_root / "network_summary.tsv",
        output_root / "adjacent_overlap.tsv",
        output_root / "anchor_overlap.tsv",
        output_root / "operating_point_screen.tsv",
        output_root / "selection.json",
        output_root / "arm_provenance.tsv",
        plots_root / "core_metrics_vs_log10_p.png",
        plots_root / "core_metrics_vs_log10_p.svg",
        plots_root / "edge_overlap_vs_log10_p.png",
        plots_root / "edge_overlap_vs_log10_p.svg",
        plots_root / "coverage_vs_nominal_null_burden.png",
        plots_root / "coverage_vs_nominal_null_burden.svg",
    ]
    if context_summary is not None:
        expected_outputs.extend(optional_context_outputs)
    output_hashes = {
        path.relative_to(output_root).as_posix(): {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in expected_outputs
    }

    analysis_manifest = {
        "schema": "sjaracne-brca100-pr67-threshold-sweep-analysis-v1",
        "work_root": str(work_root),
        "output_root": str(final_output_root),
        "anchor_p": anchor["p_value"],
        "anchor_p_key": anchor["p_key"],
        "p_keys_in_increasing_p_order": [point["p_key"] for point in points],
        "sweep_design_sha256": sweep_design_hash,
        "run_manifest_sha256": sha256_file(results_root / "run_manifest.tsv"),
        "support_aggregate_manifest_sha256": sha256_file(
            results_root / "support_summary_manifest.json"
        ),
        "netbid2_aggregate_manifest_sha256": (
            sha256_file(results_root / "netbid2_qc_manifest.json")
            if (results_root / "netbid2_qc_manifest.json").is_file()
            else None
        ),
        "build": {
            "commit": PR67_COMMIT,
            "binary_sha256": design["binary_sha256"],
            "config_sha256": design["config_sha256"],
            "null_model_sha256": NULL_MODEL_SHA256,
            "build_manifest_sha256": build["manifest_sha256"],
        },
        "candidate_files": {
            driver: {
                "path": str(work_root / "inputs" / filename),
                "sha256": sha256_file(work_root / "inputs" / filename),
                "count": len(candidates[driver]),
            }
            for driver, filename in CANDIDATE_FILES.items()
        },
        "expression_file": {
            "path": str(expression_path),
            "sha256": sha256_file(expression_path),
            "row_count": len(expression_ids),
        },
        "matched_seed_count": len(next(iter(seed_sets.values()))),
        "arms": {
            f"{row['p_key']}/{row['driver_class']}": {
                field: row[field] for field in arm_hash_fields if row[field]
            }
            for row in arm_rows
        },
        "analysis_software": {
            "script_sha256": sha256_file(Path(__file__).resolve()),
            "versions": software_versions(),
        },
        "output_files": output_hashes,
        "operating_point_selection": selection,
        "pr66_context": context_provenance,
        "interpretation_scope": (
            "Network/statistical QC and a provisional topology operating-point "
            "screen; floors were declared before intermediate sweep results but "
            "after the PR66/PR67 endpoint comparison. No biological-reference "
            "criterion or empirical FDR estimate was applied."
        ),
    }
    write_json(output_root / "analysis_manifest.json", analysis_manifest)
    output_root.replace(final_output_root)
    print(f"Wrote threshold-sweep analysis to {final_output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
