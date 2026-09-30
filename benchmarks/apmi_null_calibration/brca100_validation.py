#!/usr/bin/env python3
"""Held-out BRCA100 validation of estimator-matched SJARACNe AP-MI nulls.

The BRCA distribution is generated from independently permuted gene pairs and
is never used to fit the tail model.  Both BRCA and canonical synthetic values
come from executables linked to SJARACNe's shared production C++ AP-MI kernel.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import beta, ks_2samp


EXPECTED_KERNEL_SCHEMA = "sjaracne-apmi-v1"
MODEL_SCHEMA = "sjaracne-apmi-gpd-tail-v1"
REPORT_SCHEMA = "sjaracne-apmi-null-calibration-report-v1"
BRCA_BINARY_SCHEMA = "sjaracne-expression-permutation-null-binary-v1"
SYNTHETIC_BINARY_SCHEMA = "sjaracne-apmi-null-binary-v1"
DEFAULT_PROBABILITIES = (2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 5e-5, 2e-5)
QUANTILES = (0.5, 0.9, 0.95, 0.975, 0.99, 0.995, 0.999)


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    root = here.parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--expression", type=Path, default=root / "tests" / "inputs" / "BRCA100.exp"
    )
    parser.add_argument("--brca-generator", type=Path, required=True)
    parser.add_argument("--synthetic-generator", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--calibration-report", type=Path)
    parser.add_argument("--calibration-raw-dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--m", type=int, nargs="+", default=[80])
    parser.add_argument("--npar", type=int, default=40)
    parser.add_argument("--draws", type=int, default=100_000)
    parser.add_argument("--brca-seed", type=int, default=20260820)
    parser.add_argument("--synthetic-seed", type=int, default=20260821)
    parser.add_argument(
        "--probabilities", type=float, nargs="*", default=list(DEFAULT_PROBABILITIES)
    )
    parser.add_argument("--min-expected-exceedances", type=float, default=20.0)
    parser.add_argument("--confidence", type=float, default=0.95)
    parser.add_argument("--generator-equivalence-draws", type=int, default=10_000)
    parser.add_argument(
        "--execution", choices=("auto", "native", "wsl"), default="auto"
    )
    parser.add_argument("--reuse-raw", action="store_true")
    parser.add_argument("--allow-unaccepted-model", action="store_true")
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    for name in ("expression", "brca_generator", "synthetic_generator"):
        path = getattr(args, name)
        if not path.is_file():
            raise ValueError(f"{name.replace('_', ' ')} does not exist: {path}")
    if not args.model_dir.is_dir():
        raise ValueError(f"model directory does not exist: {args.model_dir}")
    if not args.m or any(m < 2 for m in args.m) or len(set(args.m)) != len(args.m):
        raise ValueError("--m values must be unique integers >= 2")
    args.m.sort()
    if args.npar < 1 or args.draws < 1:
        raise ValueError("--npar and --draws must be positive")
    if args.brca_seed < 0 or args.synthetic_seed < 0:
        raise ValueError("seeds must be nonnegative")
    if args.brca_seed == args.synthetic_seed:
        raise ValueError("BRCA and synthetic seeds must differ")
    if not args.probabilities or any(not 0.0 < p < 1.0 for p in args.probabilities):
        raise ValueError("--probabilities must contain values in (0,1)")
    if args.min_expected_exceedances < 1:
        raise ValueError("--min-expected-exceedances must be at least 1")
    if args.generator_equivalence_draws < 1:
        raise ValueError("--generator-equivalence-draws must be positive")
    if not 0.0 < args.confidence < 1.0:
        raise ValueError("--confidence must be in (0,1)")
    if args.calibration_report is None:
        args.calibration_report = args.model_dir / "calibration_report.json"
    if not args.calibration_report.is_file():
        raise ValueError(f"calibration report does not exist: {args.calibration_report}")
    if args.calibration_raw_dir is None:
        candidate = args.model_dir / "raw"
        args.calibration_raw_dir = candidate if candidate.is_dir() else None
    elif not args.calibration_raw_dir.is_dir():
        raise ValueError(f"calibration raw directory does not exist: {args.calibration_raw_dir}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_prefix(path: Path, byte_count: int) -> str:
    digest = hashlib.sha256()
    remaining = byte_count
    with path.open("rb") as handle:
        while remaining:
            block = handle.read(min(1024 * 1024, remaining))
            if not block:
                raise ValueError(f"{path} ended before the requested hash prefix")
            digest.update(block)
            remaining -= len(block)
    return digest.hexdigest()


def parse_metadata(path: Path) -> Dict[str, str]:
    result: Dict[str, str] = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise ValueError(f"{path}:{line_number}: expected key=value")
        key, value = line.split("=", 1)
        if not key or key in result:
            raise ValueError(f"{path}:{line_number}: empty or duplicate key")
        result[key] = value
    return result


def windows_path_to_wsl(path: Path) -> str:
    resolved = str(path.resolve())
    drive, tail = os.path.splitdrive(resolved)
    if not drive or not drive.endswith(":"):
        raise ValueError(f"cannot translate non-drive Windows path to WSL: {resolved}")
    normalized_tail = tail.lstrip("\\/").replace("\\", "/")
    return f"/mnt/{drive[0].lower()}/{normalized_tail}"


def execution_mode(requested: str) -> str:
    if requested == "auto":
        return "wsl" if os.name == "nt" else "native"
    return requested


def executable_command(executable: Path, arguments: Sequence[Tuple[str, object]], mode: str) -> List[str]:
    command = [windows_path_to_wsl(executable)] if mode == "wsl" else [str(executable)]
    if mode == "wsl":
        command.insert(0, "wsl")
    for option, value in arguments:
        command.append(option)
        if isinstance(value, Path):
            command.append(windows_path_to_wsl(value) if mode == "wsl" else str(value))
        else:
            command.append(str(value))
    return command


def run_checked(command: Sequence[str]) -> None:
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )


def expected_binary(path: Path, records: int, metadata_path: Path) -> bool:
    return path.is_file() and metadata_path.is_file() and path.stat().st_size == records * 8


def provenance_path(path: Path) -> Path:
    return Path(str(path) + ".provenance.json")


def provenance_matches(path: Path, expected: Dict[str, str]) -> bool:
    provenance = provenance_path(path)
    if not provenance.is_file():
        return False
    try:
        actual = json.loads(provenance.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return all(actual.get(key) == value for key, value in expected.items())


def write_provenance(path: Path, values: Dict[str, str]) -> None:
    provenance_path(path).write_text(
        json.dumps(values, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def load_binary(path: Path, records: int) -> np.ndarray:
    if path.stat().st_size != records * 8:
        raise ValueError(
            f"wrong binary length for {path}: {path.stat().st_size} != {records * 8}"
        )
    values = np.fromfile(path, dtype="<f8")
    if values.size != records or not np.all(np.isfinite(values)):
        raise ValueError(f"invalid or non-finite binary AP-MI values in {path}")
    return values


def run_brca_generator(args: argparse.Namespace, mode: str) -> Tuple[np.ndarray, Dict[str, str], Path]:
    raw_dir = args.output_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    output = raw_dir / "brca100_expression_permutation.bin"
    metadata_path = Path(str(output) + ".meta")
    records = args.draws * len(args.m)
    provenance = {
        "format": "sjaracne-brca-validation-raw-provenance-v1",
        "generator_sha256": sha256(args.brca_generator),
        "expression_sha256": sha256(args.expression),
    }
    if not (
        args.reuse_raw
        and expected_binary(output, records, metadata_path)
        and provenance_matches(output, provenance)
    ):
        command = executable_command(
            args.brca_generator,
            (
                ("--expression", args.expression),
                ("--m", ",".join(str(m) for m in args.m)),
                ("--draws", args.draws),
                ("--npar", args.npar),
                ("--seed", args.brca_seed),
                ("--output", output),
            ),
            mode,
        )
        run_checked(command)
        write_provenance(output, provenance)

    metadata = parse_metadata(metadata_path)
    expected = {
        "format": BRCA_BINARY_SCHEMA,
        "kernel_schema": EXPECTED_KERNEL_SCHEMA,
        "m_values": ",".join(str(m) for m in args.m),
        "draws": str(args.draws),
        "npar_limit": str(args.npar),
        "seed": str(args.brca_seed),
        "dtype": "float64",
        "byte_order": "little",
        "record_layout": "draw-major-m-values-order",
    }
    mismatches = {key: (metadata.get(key), value) for key, value in expected.items()
                  if metadata.get(key) != value}
    if mismatches:
        raise ValueError(f"BRCA generator metadata mismatch: {mismatches}")
    values = load_binary(output, records).reshape(args.draws, len(args.m))
    return values, metadata, output


def run_synthetic_generator(
    args: argparse.Namespace, m: int, mode: str
) -> Tuple[np.ndarray, Dict[str, str], Path]:
    raw_dir = args.output_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    output = raw_dir / f"synthetic_m{m:05d}.bin"
    metadata_path = Path(str(output) + ".meta")
    seed = args.synthetic_seed + m
    if seed > 2**32 - 1:
        raise ValueError("synthetic seed plus m exceeds uint32")
    provenance = {
        "format": "sjaracne-brca-validation-raw-provenance-v1",
        "generator_sha256": sha256(args.synthetic_generator),
    }
    if not (
        args.reuse_raw
        and expected_binary(output, args.draws, metadata_path)
        and provenance_matches(output, provenance)
    ):
        command = executable_command(
            args.synthetic_generator,
            (
                ("--m", m),
                ("--draws", args.draws),
                ("--npar", args.npar),
                ("--seed", seed),
                ("--format", "binary"),
                ("--output", output),
            ),
            mode,
        )
        run_checked(command)
        write_provenance(output, provenance)

    metadata = parse_metadata(metadata_path)
    expected = {
        "format": SYNTHETIC_BINARY_SCHEMA,
        "kernel_schema": EXPECTED_KERNEL_SCHEMA,
        "m": str(m),
        "draws": str(args.draws),
        "npar_limit": str(args.npar),
        "seed": str(seed),
        "dtype": "float64",
        "byte_order": "little",
    }
    mismatches = {key: (metadata.get(key), value) for key, value in expected.items()
                  if metadata.get(key) != value}
    if mismatches:
        raise ValueError(f"synthetic generator metadata mismatch for m={m}: {mismatches}")
    return load_binary(output, args.draws), metadata, output


def model_path(model_dir: Path, m: int, npar: int) -> Path:
    return model_dir / f"apmi_null_m{m:05d}_npar{npar:03d}.model"


def load_models(args: argparse.Namespace) -> Tuple[Dict[int, Dict[str, str]], Dict[int, bool], dict]:
    report = json.loads(args.calibration_report.read_text(encoding="utf-8"))
    if report.get("format") != REPORT_SCHEMA:
        raise ValueError(f"unsupported calibration report: {report.get('format')}")
    report_results = {int(result["m"]): result for result in report.get("results", [])}
    models: Dict[int, Dict[str, str]] = {}
    accepted: Dict[int, bool] = {}
    for m in args.m:
        if m not in report_results:
            raise ValueError(f"calibration report contains no result for m={m}")
        accepted[m] = bool(report_results[m].get("accepted"))
        if not accepted[m] and not args.allow_unaccepted_model:
            raise ValueError(
                f"calibration report rejects m={m}; refusing it without "
                "--allow-unaccepted-model"
            )
        path = model_path(args.model_dir, m, args.npar)
        if not path.is_file():
            raise ValueError(f"model does not exist: {path}")
        model = parse_metadata(path)
        expected = {
            "format": MODEL_SCHEMA,
            "kernel_schema": EXPECTED_KERNEL_SCHEMA,
            "sampling_null": "independent-uniform-rank-permutation",
            "rank_policy": "unique-ordinal-ranks",
            "m": str(m),
            "npar_limit": str(args.npar),
        }
        mismatches = {key: (model.get(key), value) for key, value in expected.items()
                      if model.get(key) != value}
        if mismatches:
            raise ValueError(f"model metadata mismatch for m={m}: {mismatches}")
        report_model = report_results[m].get("model")
        if not isinstance(report_model, dict):
            raise ValueError(f"calibration report has no model metadata for m={m}")
        report_mismatches = {}
        for key, report_value in report_model.items():
            file_value = model.get(key)
            if report_value is None:
                matches = file_value is not None and file_value.lower() in (
                    "none", "null", "na"
                )
            elif isinstance(report_value, (int, float)):
                try:
                    matches = float(file_value) == float(report_value)
                except (TypeError, ValueError):
                    matches = False
            else:
                matches = file_value == str(report_value)
            if not matches:
                report_mismatches[key] = (file_value, report_value)
        if report_mismatches:
            raise ValueError(
                f"model file disagrees with calibration report for m={m}: "
                f"{report_mismatches}"
            )
        models[m] = model
    return models, accepted, report


def verify_generator_equivalence(
    args: argparse.Namespace,
    mode: str,
    models: Dict[int, Dict[str, str]],
) -> Dict[int, dict]:
    """Verify a rebuilt generator against the exact calibration stream prefix."""

    current_hash = sha256(args.synthetic_generator)
    results: Dict[int, dict] = {}
    for m in args.m:
        model = models[m]
        calibrated_hash = model["generator_sha256"]
        if current_hash == calibrated_hash:
            results[m] = {
                "status": "executable-sha256-match",
                "current_generator_sha256": current_hash,
                "calibration_generator_sha256": calibrated_hash,
            }
            continue

        if args.calibration_raw_dir is None:
            results[m] = {
                "status": "executable-sha256-mismatch-not-stream-checked",
                "current_generator_sha256": current_hash,
                "calibration_generator_sha256": calibrated_hash,
            }
            continue

        calibration_values = (
            args.calibration_raw_dir / f"m{m:05d}_npar{args.npar:03d}_fit.f64"
        )
        calibration_metadata_path = Path(str(calibration_values) + ".meta")
        if not calibration_values.is_file() or not calibration_metadata_path.is_file():
            raise ValueError(
                f"generator hash changed, but calibration stream is unavailable for m={m}: "
                f"{calibration_values}"
            )
        calibration_metadata = parse_metadata(calibration_metadata_path)
        check_draws = min(args.generator_equivalence_draws, int(model["fit_draws"]))
        fit_seed = int(model["fit_seed"])
        expected_calibration = {
            "format": SYNTHETIC_BINARY_SCHEMA,
            "kernel_schema": EXPECTED_KERNEL_SCHEMA,
            "m": str(m),
            "npar_limit": str(args.npar),
            "seed": str(fit_seed),
        }
        mismatches = {
            key: (calibration_metadata.get(key), value)
            for key, value in expected_calibration.items()
            if calibration_metadata.get(key) != value
        }
        if mismatches or int(calibration_metadata.get("draws", "0")) < check_draws:
            raise ValueError(
                f"calibration-stream metadata mismatch for m={m}: {mismatches}"
            )

        output = args.output_dir / "raw" / f"generator_equivalence_m{m:05d}.bin"
        output.parent.mkdir(parents=True, exist_ok=True)
        output_metadata_path = Path(str(output) + ".meta")
        output_provenance = {
            "format": "sjaracne-brca-validation-raw-provenance-v1",
            "generator_sha256": current_hash,
        }
        if not (
            args.reuse_raw
            and expected_binary(output, check_draws, output_metadata_path)
            and parse_metadata(output_metadata_path).get("seed") == str(fit_seed)
            and provenance_matches(output, output_provenance)
        ):
            command = executable_command(
                args.synthetic_generator,
                (
                    ("--m", m),
                    ("--draws", check_draws),
                    ("--npar", args.npar),
                    ("--seed", fit_seed),
                    ("--format", "binary"),
                    ("--output", output),
                ),
                mode,
            )
            run_checked(command)
            write_provenance(output, output_provenance)

        byte_count = check_draws * 8
        with calibration_values.open("rb") as old_handle, output.open("rb") as new_handle:
            identical = old_handle.read(byte_count) == new_handle.read(byte_count)
        if not identical:
            raise ValueError(
                f"current generator differs from the calibration stream for m={m} "
                f"within the first {check_draws} draws"
            )
        results[m] = {
            "status": "executable-sha256-differs-but-calibration-prefix-byte-identical",
            "current_generator_sha256": current_hash,
            "calibration_generator_sha256": calibrated_hash,
            "draws_compared": check_draws,
            "fit_seed": fit_seed,
            "calibration_prefix_sha256": sha256_prefix(calibration_values, byte_count),
            "current_prefix_sha256": sha256_prefix(output, byte_count),
        }
    return results


def model_float(model: Dict[str, str], key: str) -> float:
    try:
        value = float(model[key])
    except (KeyError, ValueError) as error:
        raise ValueError(f"invalid model field {key!r}") from error
    if not math.isfinite(value):
        raise ValueError(f"non-finite model field {key!r}")
    return value


def model_cutoff(model: Dict[str, str], probability: float) -> float:
    p_min = model_float(model, "supported_p_min")
    p_max = model_float(model, "supported_p_max")
    tail_probability = model_float(model, "tail_probability")
    if not p_min <= probability <= min(p_max, tail_probability):
        raise ValueError(f"probability {probability:g} is outside model support")
    threshold = model_float(model, "tail_threshold")
    shape = model_float(model, "tail_shape")
    scale = model_float(model, "tail_scale")
    if scale <= 0.0:
        raise ValueError("tail_scale must be positive")
    ratio = probability / tail_probability
    if abs(shape) < 1e-10:
        cutoff = threshold + scale * math.log(1.0 / ratio)
    else:
        cutoff = threshold + scale / shape * (ratio ** (-shape) - 1.0)
    endpoint_text = model.get("tail_endpoint", "none").strip().lower()
    if endpoint_text not in ("", "none", "null", "na"):
        cutoff = min(cutoff, float(endpoint_text))
    if not math.isfinite(cutoff):
        raise ValueError("model produced a non-finite cutoff")
    return cutoff


def model_survival(model: Dict[str, str], values: np.ndarray) -> np.ndarray:
    threshold = model_float(model, "tail_threshold")
    p_tail = model_float(model, "tail_probability")
    shape = model_float(model, "tail_shape")
    scale = model_float(model, "tail_scale")
    scaled = 1.0 + shape * (values - threshold) / scale
    survival = np.full(values.shape, np.nan, dtype=float)
    valid = (values >= threshold) & (scaled > 0.0)
    if abs(shape) < 1e-10:
        survival[values >= threshold] = p_tail * np.exp(
            -(values[values >= threshold] - threshold) / scale
        )
    else:
        survival[valid] = p_tail * scaled[valid] ** (-1.0 / shape)
    return survival


def quantile(values: np.ndarray, q: Sequence[float]) -> np.ndarray:
    try:
        return np.quantile(values, q, method="linear")
    except TypeError:
        return np.quantile(values, q, interpolation="linear")


def clopper_pearson(count: int, total: int, confidence: float) -> Tuple[float, float]:
    alpha = 1.0 - confidence
    lower = 0.0 if count == 0 else float(beta.ppf(alpha / 2.0, count, total - count + 1))
    upper = 1.0 if count == total else float(
        beta.ppf(1.0 - alpha / 2.0, count + 1, total - count)
    )
    return lower, upper


def distribution_rows(m: int, source: str, values: np.ndarray) -> dict:
    q_values = quantile(values, QUANTILES)
    row = {
        "m": m,
        "source": source,
        "draws": int(values.size),
        "mean": float(np.mean(values)),
        "sd": float(np.std(values, ddof=1)),
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
    }
    for q, value in zip(QUANTILES, q_values):
        row[f"q{q:g}"] = float(value)
    return row


def csv_write(path: Path, rows: Sequence[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def survival_points(values: np.ndarray, maximum_points: int = 5000) -> Tuple[np.ndarray, np.ndarray]:
    ordered = np.sort(values)
    n = ordered.size
    if n > maximum_points:
        indices = np.unique(np.linspace(0, n - 1, maximum_points).astype(int))
        ordered = ordered[indices]
    else:
        indices = np.arange(n)
    survival = (n - indices) / float(n)
    return ordered, survival


def plot_ecdf(
    path_base: Path,
    m_values: Sequence[int],
    brca: Dict[int, np.ndarray],
    synthetic: Dict[int, np.ndarray],
    models: Dict[int, Dict[str, str]],
) -> None:
    figure, axes = plt.subplots(1, len(m_values), figsize=(6.2 * len(m_values), 4.8),
                                squeeze=False)
    for axis, m in zip(axes[0], m_values):
        for values, label, color in (
            (brca[m], "Held-out BRCA permutation", "#167D9A"),
            (synthetic[m], "Canonical rank permutation", "#D55E00"),
        ):
            x, y = survival_points(values)
            axis.plot(x, y, label=label, color=color, linewidth=1.8)
        model = models[m]
        threshold = model_float(model, "tail_threshold")
        x_max = max(float(np.max(brca[m])), float(np.max(synthetic[m])),
                    model_cutoff(model, max(model_float(model, "supported_p_min"), 1e-7)))
        x_model = np.linspace(threshold, x_max, 400)
        y_model = model_survival(model, x_model)
        valid = np.isfinite(y_model) & (y_model > 0)
        axis.plot(x_model[valid], y_model[valid], color="#4D4D4D", linestyle="--",
                  linewidth=1.7, label="Fitted GPD tail")
        axis.set_yscale("log")
        axis.set_ylim(max(0.5 / brca[m].size, 1e-7), 1.0)
        axis.set_xlabel("AP-MI")
        axis.set_ylabel("Empirical survival probability")
        axis.set_title(f"m={m}, Npar={models[m]['npar_limit']}")
        axis.grid(True, which="both", alpha=0.22)
        axis.legend(frameon=False, fontsize=8)
    figure.suptitle("BRCA100 held-out permutation null vs estimator-matched null")
    figure.tight_layout()
    for suffix in ("png", "svg"):
        figure.savefig(path_base.with_suffix(f".{suffix}"), dpi=220, bbox_inches="tight",
                       metadata={"Date": None})
    plt.close(figure)


def plot_qq(path_base: Path, m_values: Sequence[int], brca: Dict[int, np.ndarray],
            synthetic: Dict[int, np.ndarray]) -> None:
    figure, axes = plt.subplots(1, len(m_values), figsize=(5.5 * len(m_values), 4.8),
                                squeeze=False)
    tail_probabilities = np.unique(np.geomspace(1e-4, 0.5, 240))
    q = 1.0 - tail_probabilities
    for axis, m in zip(axes[0], m_values):
        x = quantile(synthetic[m], q)
        y = quantile(brca[m], q)
        limit = max(float(np.max(x)), float(np.max(y)))
        axis.scatter(x, y, c=-np.log10(tail_probabilities), cmap="viridis",
                     s=14, alpha=0.85, linewidths=0)
        axis.plot([0, limit], [0, limit], color="#555555", linestyle="--", linewidth=1.2)
        axis.set_xlabel("Canonical-null AP-MI quantile")
        axis.set_ylabel("BRCA-null AP-MI quantile")
        axis.set_title(f"m={m}")
        axis.grid(True, alpha=0.22)
    figure.suptitle("AP-MI null quantile agreement (yellow = deeper tail)")
    figure.tight_layout()
    for suffix in ("png", "svg"):
        figure.savefig(path_base.with_suffix(f".{suffix}"), dpi=220, bbox_inches="tight",
                       metadata={"Date": None})
    plt.close(figure)


def plot_exceedance(path_base: Path, rows: Sequence[dict], m_values: Sequence[int]) -> None:
    figure, axes = plt.subplots(1, len(m_values), figsize=(5.5 * len(m_values), 4.8),
                                squeeze=False)
    for axis, m in zip(axes[0], m_values):
        subset = [row for row in rows if row["m"] == m]
        p = np.asarray([row["nominal_p"] for row in subset])
        brca_ratio = np.asarray([row["brca_ratio_to_nominal"] for row in subset])
        synthetic_ratio = np.asarray([row["synthetic_ratio_to_nominal"] for row in subset])
        brca_lower = np.asarray([row["brca_ci_lower"] for row in subset]) / p
        brca_upper = np.asarray([row["brca_ci_upper"] for row in subset]) / p
        synthetic_lower = np.asarray([row["synthetic_ci_lower"] for row in subset]) / p
        synthetic_upper = np.asarray([row["synthetic_ci_upper"] for row in subset]) / p
        axis.errorbar(
            p,
            brca_ratio,
            yerr=np.vstack((brca_ratio - brca_lower, brca_upper - brca_ratio)),
            marker="o",
            color="#167D9A",
            capsize=3,
            label="Held-out BRCA",
        )
        axis.errorbar(
            p,
            synthetic_ratio,
            yerr=np.vstack(
                (synthetic_ratio - synthetic_lower, synthetic_upper - synthetic_ratio)
            ),
            marker="s",
            color="#D55E00",
            capsize=3,
            label="Canonical null",
        )
        axis.axhline(1.0, color="#555555", linestyle="--", linewidth=1.2)
        axis.set_xscale("log")
        axis.invert_xaxis()
        axis.set_xticks(p)
        axis.set_xticklabels([f"{value:g}" for value in p])
        axis.set_xlabel("Nominal fitted-tail probability")
        axis.set_ylabel("Observed / nominal exceedance ratio")
        axis.set_title(f"m={m}")
        axis.grid(True, which="both", alpha=0.22)
        axis.legend(frameon=False, fontsize=8)
    figure.suptitle("Held-out exceedance calibration at directly resolvable probabilities")
    figure.tight_layout()
    for suffix in ("png", "svg"):
        figure.savefig(path_base.with_suffix(f".{suffix}"), dpi=220, bbox_inches="tight",
                       metadata={"Date": None})
    plt.close(figure)


def main() -> int:
    args = parse_args()
    validate_args(args)
    script_path = Path(__file__).resolve()
    benchmark_dir = script_path.parent
    generator_source = benchmark_dir / "brca100_validation_generator.cpp"
    apmi_source = benchmark_dir.parent.parent / "SJARACNe" / "src" / "apmi.cpp"
    if not generator_source.is_file() or not apmi_source.is_file():
        raise ValueError("unable to locate tracked BRCA generator or shared AP-MI source")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    mode = execution_mode(args.execution)
    models, model_accepted, calibration_report = load_models(args)
    generator_equivalence = verify_generator_equivalence(args, mode, models)
    brca_matrix, brca_metadata, brca_raw = run_brca_generator(args, mode)

    brca = {m: brca_matrix[:, index].copy() for index, m in enumerate(args.m)}
    synthetic: Dict[int, np.ndarray] = {}
    synthetic_metadata: Dict[int, dict] = {}
    synthetic_raw: Dict[int, Path] = {}
    for m in args.m:
        values, metadata, path = run_synthetic_generator(args, m, mode)
        synthetic[m] = values
        synthetic_metadata[m] = metadata
        synthetic_raw[m] = path

    distribution_summary: List[dict] = []
    quantile_comparison: List[dict] = []
    exceedance_rows: List[dict] = []
    comparison_summary: List[dict] = []
    directly_resolvable_min = args.min_expected_exceedances / args.draws

    for m in args.m:
        for source, values in (("brca", brca[m]), ("synthetic", synthetic[m])):
            tolerance = 64 * np.finfo(np.float64).eps
            if np.any(values < -tolerance) or np.any(values > math.log(m) + tolerance):
                raise ValueError(f"{source} values for m={m} fall outside [0, log(m)]")
            distribution_summary.append(distribution_rows(m, source, values))

        brca_q = quantile(brca[m], QUANTILES)
        synthetic_q = quantile(synthetic[m], QUANTILES)
        for q, b_value, s_value in zip(QUANTILES, brca_q, synthetic_q):
            quantile_comparison.append(
                {
                    "m": m,
                    "quantile": q,
                    "brca_mi": float(b_value),
                    "synthetic_mi": float(s_value),
                    "brca_minus_synthetic": float(b_value - s_value),
                    "brca_to_synthetic": (
                        float(b_value / s_value) if s_value != 0 else None
                    ),
                }
            )

        model = models[m]
        probabilities = sorted(
            {
                float(p)
                for p in args.probabilities
                if p >= directly_resolvable_min
                and model_float(model, "supported_p_min") <= p
                <= min(model_float(model, "supported_p_max"),
                       model_float(model, "tail_probability"))
            },
            reverse=True,
        )
        if not probabilities:
            raise ValueError(
                f"no requested probability is directly resolvable for m={m}; "
                f"need p >= {directly_resolvable_min:g}"
            )
        for p in probabilities:
            cutoff = model_cutoff(model, p)
            b_count = int(np.count_nonzero(brca[m] >= cutoff))
            s_count = int(np.count_nonzero(synthetic[m] >= cutoff))
            b_rate = b_count / args.draws
            s_rate = s_count / args.draws
            b_lower, b_upper = clopper_pearson(b_count, args.draws, args.confidence)
            s_lower, s_upper = clopper_pearson(s_count, args.draws, args.confidence)
            exceedance_rows.append(
                {
                    "m": m,
                    "nominal_p": p,
                    "cutoff": cutoff,
                    "expected_exceedances": p * args.draws,
                    "brca_exceedances": b_count,
                    "brca_observed_p": b_rate,
                    "brca_ci_lower": b_lower,
                    "brca_ci_upper": b_upper,
                    "brca_nominal_in_pointwise_ci": b_lower <= p <= b_upper,
                    "brca_ratio_to_nominal": b_rate / p,
                    "synthetic_exceedances": s_count,
                    "synthetic_observed_p": s_rate,
                    "synthetic_ci_lower": s_lower,
                    "synthetic_ci_upper": s_upper,
                    "synthetic_nominal_in_pointwise_ci": s_lower <= p <= s_upper,
                    "synthetic_ratio_to_nominal": s_rate / p,
                    "brca_ratio_to_synthetic": (
                        b_rate / s_rate if s_rate != 0 else None
                    ),
                }
            )

        ks = ks_2samp(brca[m], synthetic[m], alternative="two-sided", method="asymp")
        m_exceedances = [row for row in exceedance_rows if row["m"] == m]
        comparison_summary.append(
            {
                "m": m,
                "model_accepted_by_calibration": model_accepted[m],
                "ks_statistic": float(ks.statistic),
                "ks_pvalue": float(ks.pvalue),
                "max_abs_reported_quantile_delta_q50_to_q999": float(
                    np.max(np.abs(brca_q - synthetic_q))
                ),
                "mean_difference": float(np.mean(brca[m]) - np.mean(synthetic[m])),
                "mean_relative_difference": (
                    float(np.mean(brca[m]) / np.mean(synthetic[m]) - 1.0)
                    if np.mean(synthetic[m]) != 0.0
                    else None
                ),
                "all_brca_nominal_in_pointwise_ci": all(
                    row["brca_nominal_in_pointwise_ci"] for row in m_exceedances
                ),
                "all_synthetic_nominal_in_pointwise_ci": all(
                    row["synthetic_nominal_in_pointwise_ci"] for row in m_exceedances
                ),
                "direct_validation_p_min": min(probabilities),
                "default_model_p": model_float(model, "default_p"),
                "default_p_directly_tested": False,
            }
        )

    csv_write(args.output_dir / "distribution_summary.csv", distribution_summary)
    csv_write(args.output_dir / "quantile_comparison.csv", quantile_comparison)
    csv_write(args.output_dir / "model_exceedance.csv", exceedance_rows)
    csv_write(args.output_dir / "comparison_summary.csv", comparison_summary)

    matplotlib.rcParams["svg.hashsalt"] = "sjaracne-brca-null-v1"
    plot_ecdf(args.output_dir / "brca100_null_survival", args.m, brca, synthetic, models)
    plot_qq(args.output_dir / "brca100_null_qq", args.m, brca, synthetic)
    plot_exceedance(args.output_dir / "brca100_model_exceedance", exceedance_rows, args.m)

    report = {
        "format": "sjaracne-brca-heldout-null-validation-v1",
        "command": sys.argv,
        "cwd": str(Path.cwd()),
        "design": (
            "BRCA100 is held out from fitting. Each draw selects two distinct, "
            "nonconstant, tie-free BRCA rows, selects m observations without "
            "replacement, permutes gene Y within that selected subset, constructs "
            "SJARACNe ordinal ranks, and calls the shared production C++ AP-MI kernel."
        ),
        "interpretation": (
            "For tie-free continuous data, ranks make the expression-permutation "
            "null distribution-free, so agreement with the canonical rank null is "
            "expected. This validates implementation and fitted-tail behavior at "
            "Monte Carlo-resolvable probabilities; it does not refit the model."
        ),
        "explicit_limitation": (
            f"With {args.draws} draws and a minimum of "
            f"{args.min_expected_exceedances:g} expected exceedances, direct checks "
            f"stop at p={directly_resolvable_min:g}. The default p=1e-7 is a tail "
            "extrapolation and is not directly validated here."
        ),
        "expression": str(args.expression.resolve()),
        "expression_sha256": sha256(args.expression),
        "m_values": args.m,
        "npar_limit": args.npar,
        "draws_per_distribution_per_m": args.draws,
        "brca_seed": args.brca_seed,
        "synthetic_seed_base": args.synthetic_seed,
        "confidence": args.confidence,
        "minimum_expected_exceedances": args.min_expected_exceedances,
        "execution_mode": mode,
        "brca_generator": str(args.brca_generator.resolve()),
        "brca_generator_sha256": sha256(args.brca_generator),
        "validation_script_sha256": sha256(script_path),
        "brca_generator_source_sha256": sha256(generator_source),
        "shared_apmi_source_sha256": sha256(apmi_source),
        "synthetic_generator": str(args.synthetic_generator.resolve()),
        "synthetic_generator_sha256": sha256(args.synthetic_generator),
        "generator_equivalence_to_calibration": generator_equivalence,
        "calibration_report": str(args.calibration_report.resolve()),
        "calibration_report_sha256": sha256(args.calibration_report),
        "calibration_report_format": calibration_report["format"],
        "calibration_raw_dir": (
            str(args.calibration_raw_dir.resolve())
            if args.calibration_raw_dir is not None
            else None
        ),
        "brca_raw": str(brca_raw.resolve()),
        "brca_raw_sha256": sha256(brca_raw),
        "brca_metadata": brca_metadata,
        "synthetic_raw": {str(m): str(synthetic_raw[m].resolve()) for m in args.m},
        "synthetic_raw_sha256": {str(m): sha256(synthetic_raw[m]) for m in args.m},
        "synthetic_metadata": synthetic_metadata,
        "models": {
            str(m): {
                "path": str(model_path(args.model_dir, m, args.npar).resolve()),
                "sha256": sha256(model_path(args.model_dir, m, args.npar)),
                "accepted_by_calibration": model_accepted[m],
                "metadata": models[m],
            }
            for m in args.m
        },
        "comparison_summary": comparison_summary,
        "python": sys.version,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
    }
    (args.output_dir / "validation_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(comparison_summary, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
