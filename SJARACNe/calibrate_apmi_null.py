#!/usr/bin/env python3
"""Fit and validate exact-m SJARACNe AP-MI permutation-null tails.

Every null MI is produced by ``apmi_null_generator.exe``, which links the same
C++ adaptive-partitioning kernel as the network executable.  This script only
orchestrates independent fit/validation simulations and fits a peaks-over-
threshold generalized Pareto tail.  Biological expression data are deliberately
not used for fitting; they belong in separate held-out validation.
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
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import scipy
from scipy.stats import beta, genpareto


MODEL_SCHEMA = "sjaracne-apmi-gpd-tail-v1"
CALIBRATOR_SCHEMA = "sjaracne-apmi-gpd-calibrator-v1"
EXPECTED_KERNEL_SCHEMA = "sjaracne-apmi-v1"
DEFAULT_P = 1e-7
DEFAULT_VALIDATION_P = (1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 5e-5, 2e-5, 1e-5)
STABILITY_QUANTILES = (0.9925, 0.995, 0.9975)


@dataclass(frozen=True)
class TailFit:
    threshold_quantile: float
    threshold: float
    tail_count: int
    tail_probability: float
    shape: float
    scale: float
    endpoint: Optional[float]
    cutoff_at_stability_probability: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generator", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--m", type=int, nargs="+", required=True)
    parser.add_argument("--npar", type=int, default=40)
    parser.add_argument("--fit-draws", type=int, default=5_000_000)
    parser.add_argument("--validation-draws", type=int, default=5_000_000)
    parser.add_argument("--fit-seed", type=int, default=20260813)
    parser.add_argument("--validation-seed", type=int, default=20260814)
    parser.add_argument("--tail-quantile", type=float, default=0.995)
    parser.add_argument("--min-tail-count", type=int, default=20_000)
    parser.add_argument("--p-min", type=float, default=1e-7)
    parser.add_argument("--p-max", type=float, default=1e-2)
    parser.add_argument("--default-p", type=float, default=DEFAULT_P)
    parser.add_argument(
        "--validation-p",
        type=float,
        nargs="*",
        default=list(DEFAULT_VALIDATION_P),
    )
    parser.add_argument("--stability-relative-tolerance", type=float, default=0.10)
    parser.add_argument("--min-stability-fits", type=int, default=3)
    parser.add_argument("--validation-family-confidence", type=float, default=0.99)
    parser.add_argument(
        "--validation-family-comparisons",
        type=int,
        default=len(DEFAULT_VALIDATION_P),
    )
    parser.add_argument("--reuse-raw", action="store_true")
    parser.add_argument("--allow-rejected-calibration", action="store_true")
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if not args.generator.is_file():
        raise ValueError(f"null generator does not exist: {args.generator}")
    if not args.m or any(m < 2 for m in args.m) or len(set(args.m)) != len(args.m):
        raise ValueError("--m values must be unique integers >= 2")
    if args.npar < 1:
        raise ValueError("--npar must be positive")
    if args.fit_draws < 1 or args.validation_draws < 1:
        raise ValueError("fit and validation draw counts must be positive")
    if not 0.5 < args.tail_quantile < 1.0:
        raise ValueError("--tail-quantile must be in (0.5, 1)")
    if args.min_tail_count < 100:
        raise ValueError("--min-tail-count must be at least 100")
    if not 0.0 < args.p_min < args.p_max <= 0.05:
        raise ValueError("require 0 < --p-min < --p-max <= 0.05")
    if not args.p_min <= args.default_p <= args.p_max:
        raise ValueError("--default-p must be within [p-min, p-max]")
    if args.p_min != args.default_p:
        raise ValueError("v1 requires --p-min to equal --default-p")
    if any(not args.p_min <= p <= args.p_max for p in args.validation_p):
        raise ValueError("every --validation-p must be within [p-min, p-max]")
    if args.fit_seed < 0 or args.validation_seed < 0:
        raise ValueError("seeds must be nonnegative")
    if args.fit_seed == args.validation_seed:
        raise ValueError("fit and validation seeds must differ")
    fit_seeds = {args.fit_seed + m for m in args.m}
    validation_seeds = {args.validation_seed + 1_000_000_000 + m for m in args.m}
    if fit_seeds.intersection(validation_seeds):
        raise ValueError("derived fit and validation seed streams overlap across m")
    if not 0.0 < args.stability_relative_tolerance < 1.0:
        raise ValueError("--stability-relative-tolerance must be in (0,1)")
    if args.min_stability_fits < 2:
        raise ValueError("--min-stability-fits must be at least 2")
    if not 0.0 < args.validation_family_confidence < 1.0:
        raise ValueError("--validation-family-confidence must be in (0,1)")
    if args.validation_family_comparisons < len(args.validation_p):
        raise ValueError(
            "--validation-family-comparisons must cover every requested validation p"
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def calibrator_sha256() -> str:
    return sha256(Path(__file__).resolve())


def parse_metadata(path: Path) -> Dict[str, str]:
    result = {}  # type: Dict[str, str]
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw or raw.startswith("#"):
            continue
        if "=" not in raw:
            raise ValueError(f"{path}:{line_number}: expected key=value metadata")
        key, value = raw.split("=", 1)
        if not key or key in result:
            raise ValueError(f"{path}:{line_number}: duplicate or empty metadata key")
        result[key] = value
    return result


def run_generator(
    generator: Path,
    output: Path,
    *,
    m: int,
    draws: int,
    npar: int,
    seed: int,
    reuse: bool,
) -> Tuple[np.ndarray, Dict[str, str]]:
    metadata_path = Path(str(output) + ".meta")
    provenance_path = Path(str(output) + ".provenance.json")
    current_generator_sha256 = sha256(generator)
    expected_bytes = draws * np.dtype("<f8").itemsize
    reusable = (
        reuse
        and output.is_file()
        and metadata_path.is_file()
        and provenance_path.is_file()
        and output.stat().st_size == expected_bytes
    )
    if reusable:
        try:
            provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
            reusable = (
                provenance.get("format") == "sjaracne-apmi-null-raw-provenance-v1"
                and provenance.get("generator_sha256") == current_generator_sha256
            )
        except (OSError, ValueError, TypeError):
            reusable = False
    if not reusable:
        output.parent.mkdir(parents=True, exist_ok=True)
        command = [
            str(generator),
            "--m",
            str(m),
            "--draws",
            str(draws),
            "--npar",
            str(npar),
            "--seed",
            str(seed),
            "--format",
            "binary",
            "--output",
            str(output),
        ]
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        if completed.returncode != 0:
            raise RuntimeError(
                f"null generator failed ({completed.returncode}): {' '.join(command)}\n"
                f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
            )
        write_text_lf(
            provenance_path,
            json.dumps(
                {
                    "format": "sjaracne-apmi-null-raw-provenance-v1",
                    "generator_sha256": current_generator_sha256,
                },
                indent=2,
                allow_nan=False,
            )
            + "\n",
        )

    metadata = parse_metadata(metadata_path)
    expected = {
        "format": "sjaracne-apmi-null-binary-v1",
        "kernel_schema": EXPECTED_KERNEL_SCHEMA,
        "m": str(m),
        "draws": str(draws),
        "npar_limit": str(npar),
        "seed": str(seed),
        "dtype": "float64",
        "byte_order": "little",
        "record_bytes": "8",
    }
    mismatches = {
        key: (metadata.get(key), value)
        for key, value in expected.items()
        if metadata.get(key) != value
    }
    if mismatches:
        raise ValueError(f"null-generator metadata mismatch for {output}: {mismatches}")
    if output.stat().st_size != expected_bytes:
        raise ValueError(
            f"wrong binary length for {output}: expected {expected_bytes}, "
            f"found {output.stat().st_size}"
        )

    values = np.fromfile(output, dtype="<f8")
    if values.size != draws:
        raise ValueError(f"wrong draw count in {output}: {values.size} != {draws}")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"non-finite AP-MI in {output}")
    tolerance = 64 * np.finfo(np.float64).eps
    if np.any(values < -tolerance) or np.any(values > math.log(m) + tolerance):
        raise ValueError(f"AP-MI outside [0, log(m)] in {output}")
    return values, metadata


def higher_quantile(values: np.ndarray, quantile: float) -> float:
    # ``method`` replaced ``interpolation`` in NumPy 1.22.  Retain compatibility
    # with the project's NumPy >=1.20 requirement.
    try:
        return float(np.quantile(values, quantile, method="higher"))
    except TypeError:
        return float(np.quantile(values, quantile, interpolation="higher"))


def gpd_cutoff(fit: TailFit, probability: float) -> float:
    if not 0.0 < probability <= fit.tail_probability:
        raise ValueError(
            f"probability {probability:g} is outside fitted tail "
            f"(0,{fit.tail_probability:g}]"
        )
    ratio = probability / fit.tail_probability
    if abs(fit.shape) < 1e-10:
        cutoff = fit.threshold + fit.scale * math.log(1.0 / ratio)
    else:
        cutoff = fit.threshold + fit.scale / fit.shape * (
            math.pow(ratio, -fit.shape) - 1.0
        )
    if fit.endpoint is not None:
        cutoff = min(cutoff, fit.endpoint)
    return cutoff


def fit_gpd(values: np.ndarray, quantile: float, p_min: float, m: int) -> TailFit:
    threshold = higher_quantile(values, quantile)
    excesses = values[values > threshold] - threshold
    if excesses.size < 100:
        raise ValueError(
            f"only {excesses.size} strict exceedances above q={quantile:g}; "
            "cannot fit a stable tail"
        )
    shape, location, scale = genpareto.fit(excesses, floc=0.0)
    shape = float(shape)
    location = float(location)
    scale = float(scale)
    if location != 0.0 or not math.isfinite(shape) or not math.isfinite(scale) or scale <= 0:
        raise ValueError("invalid generalized-Pareto fit")
    tail_probability = (int(excesses.size) + 1.0) / (values.size + 1.0)
    endpoint = threshold - scale / shape if shape < 0.0 else None
    theoretical_maximum = math.log(m)
    support_tolerance = 64 * np.finfo(np.float64).eps * max(1.0, theoretical_maximum)
    if endpoint is not None and endpoint + support_tolerance < float(np.max(values)):
        raise ValueError("fitted GPD endpoint is below an observed AP-MI value")
    if endpoint is not None and endpoint > theoretical_maximum + support_tolerance:
        raise ValueError(
            "fitted GPD endpoint exceeds the theoretical AP-MI maximum log(m)"
        )
    placeholder = TailFit(
        threshold_quantile=quantile,
        threshold=threshold,
        tail_count=int(excesses.size),
        tail_probability=tail_probability,
        shape=shape,
        scale=scale,
        endpoint=endpoint,
        cutoff_at_stability_probability=math.nan,
    )
    cutoff = gpd_cutoff(placeholder, p_min)
    if (
        not math.isfinite(cutoff)
        or cutoff < threshold
        or cutoff > theoretical_maximum + support_tolerance
    ):
        raise ValueError("tail model produced an invalid extreme cutoff")
    return TailFit(
        **{
            **asdict(placeholder),
            "cutoff_at_stability_probability": cutoff,
        }
    )


def clopper_pearson(count: int, total: int, confidence: float = 0.99) -> Tuple[float, float]:
    alpha = 1.0 - confidence
    lower = 0.0 if count == 0 else float(beta.ppf(alpha / 2.0, count, total - count + 1))
    upper = 1.0 if count == total else float(
        beta.ppf(1.0 - alpha / 2.0, count + 1, total - count)
    )
    return lower, upper


def validate_tail(
    fit: TailFit,
    values: np.ndarray,
    probabilities: Iterable[float],
    confidence: float,
) -> List[Dict[str, object]]:
    rows = []  # type: List[Dict[str, object]]
    for probability in sorted(set(probabilities), reverse=True):
        if probability > fit.tail_probability:
            continue
        cutoff = gpd_cutoff(fit, probability)
        count = int(np.count_nonzero(values >= cutoff))
        observed = count / values.size
        lower, upper = clopper_pearson(count, values.size, confidence)
        rows.append(
            {
                "nominal_p": probability,
                "cutoff": cutoff,
                "exceedances": count,
                "draws": int(values.size),
                "observed_p": observed,
                "ratio_observed_to_nominal": observed / probability,
                "interval_lower": lower,
                "interval_upper": upper,
                "interval_confidence": confidence,
                "nominal_in_interval": lower <= probability <= upper,
            }
        )
    return rows


def write_text_lf(path: Path, text: str) -> None:
    # pathlib.Path.write_text(newline=...) is unavailable on Python 3.7.
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def write_model(path: Path, model: Dict[str, object]) -> None:
    ordered_keys = (
        "format",
        "kernel_schema",
        "estimator",
        "sampling_null",
        "rank_policy",
        "m",
        "npar_limit",
        "tail_model",
        "tail_threshold_quantile",
        "tail_threshold",
        "tail_probability",
        "tail_shape",
        "tail_scale",
        "tail_endpoint",
        "calibration_status",
        "calibrator_schema",
        "calibrator_sha256",
        "validation_method",
        "stability_probability",
        "stability_relative_range",
        "stability_relative_tolerance",
        "validation_family_confidence",
        "validation_point_confidence",
        "supported_p_min",
        "supported_p_max",
        "validated_p_min",
        "validated_p_max",
        "default_p",
        "default_p_cutoff",
        "fit_draws",
        "validation_draws",
        "fit_seed",
        "validation_seed",
        "rng",
        "generator_sha256",
        "fit_values_sha256",
        "validation_values_sha256",
        "scipy_version",
    )
    lines = []
    for key in ordered_keys:
        value = model[key]
        if value is None:
            value = "none"
        elif isinstance(value, float):
            value = format(value, ".17g")
        lines.append(f"{key}={value}")
    temporary = path.with_name(path.name + ".tmp")
    write_text_lf(temporary, "\n".join(lines) + "\n")
    os.replace(str(temporary), str(path))


def contiguous_validated_minimum(validation: List[Dict[str, object]]) -> Optional[float]:
    """Return the lowest p in the uninterrupted, descending passing range."""
    validated = None  # type: Optional[float]
    for row in sorted(validation, key=lambda item: float(item["nominal_p"]), reverse=True):
        if not bool(row["nominal_in_interval"]):
            break
        validated = float(row["nominal_p"])
    return validated


def calibrate_one(args: argparse.Namespace, m: int) -> Dict[str, object]:
    raw_dir = args.output_dir / "raw"
    fit_path = raw_dir / f"m{m:05d}_npar{args.npar:03d}_fit.f64"
    validation_path = raw_dir / f"m{m:05d}_npar{args.npar:03d}_validation.f64"
    model_path = args.output_dir / f"apmi_null_m{m:05d}_npar{args.npar:03d}.model"
    if model_path.exists():
        model_path.unlink()

    # Use disjoint numeric namespaces so a fit stream at one m can never become
    # a validation stream at another m in a dense sweep.
    fit_seed = args.fit_seed + m
    validation_seed = args.validation_seed + 1_000_000_000 + m
    if fit_seed > 2**32 - 1 or validation_seed > 2**32 - 1:
        raise ValueError("derived generator seed exceeds uint32")
    if fit_seed == validation_seed:
        raise ValueError("derived fit and validation seeds must differ")

    fit_values, fit_metadata = run_generator(
        args.generator,
        fit_path,
        m=m,
        draws=args.fit_draws,
        npar=args.npar,
        seed=fit_seed,
        reuse=args.reuse_raw,
    )
    validation_values, validation_metadata = run_generator(
        args.generator,
        validation_path,
        m=m,
        draws=args.validation_draws,
        npar=args.npar,
        seed=validation_seed,
        reuse=args.reuse_raw,
    )

    # The model is intended to serve the default p, so acceptance must test
    # threshold sensitivity at that extrapolated probability itself.  Testing at
    # an easier, empirically resolved p can hide an unstable default cutoff.
    stability_probability = args.default_p
    primary_fit = fit_gpd(fit_values, args.tail_quantile, stability_probability, m)
    if primary_fit.tail_count < args.min_tail_count:
        raise ValueError(
            f"m={m}: tail has {primary_fit.tail_count} points, below "
            f"--min-tail-count={args.min_tail_count}"
        )

    stability_fits = []  # type: List[TailFit]
    stability_quantiles = sorted(set(STABILITY_QUANTILES + (args.tail_quantile,)))
    stability_tail_minimum = max(100, args.min_tail_count // 10)
    skipped_stability_quantiles = []  # type: List[Dict[str, object]]
    for quantile in stability_quantiles:
        try:
            candidate = fit_gpd(fit_values, quantile, stability_probability, m)
        except ValueError as error:
            skipped_stability_quantiles.append(
                {"quantile": quantile, "reason": str(error)}
            )
            continue
        if candidate.tail_count >= stability_tail_minimum:
            stability_fits.append(candidate)
        else:
            skipped_stability_quantiles.append(
                {
                    "quantile": quantile,
                    "reason": "tail count {} is below stability minimum {}".format(
                        candidate.tail_count, stability_tail_minimum
                    ),
                }
            )
    relative_stability_range = None  # type: Optional[float]
    if len(stability_fits) >= 2:
        cutoffs = np.asarray(
            [fit.cutoff_at_stability_probability for fit in stability_fits]
        )
        cutoff_median = float(np.median(cutoffs))
        if cutoff_median != 0.0:
            relative_stability_range = float(np.ptp(cutoffs) / cutoff_median)
    stability_pass = (
        len(stability_fits) >= args.min_stability_fits
        and relative_stability_range is not None
        and math.isfinite(relative_stability_range)
        and relative_stability_range <= args.stability_relative_tolerance
    )

    directly_testable_floor = 100.0 / args.validation_draws
    modeled_validation_probabilities = sorted(
        set(
            probability
            for probability in args.validation_p
            if probability <= primary_fit.tail_probability
        ),
        reverse=True,
    )
    directly_testable_count = sum(
        probability >= directly_testable_floor
        for probability in modeled_validation_probabilities
    )
    # Keep one model's acceptance independent of which other m values happen to be
    # calibrated in the same invocation. The fixed family size covers the full
    # predeclared validation grid, including levels too extreme to test directly.
    total_family_comparisons = args.validation_family_comparisons
    point_confidence = 1.0 - (
        (1.0 - args.validation_family_confidence) / total_family_comparisons
    )
    validation = validate_tail(
        primary_fit,
        validation_values,
        modeled_validation_probabilities,
        point_confidence,
    )
    directly_testable = [
        row
        for row in validation
        if float(row["nominal_p"]) >= directly_testable_floor
    ]
    validation_pass = bool(directly_testable) and all(
        bool(row["nominal_in_interval"]) for row in directly_testable
    )
    validated_p_min = contiguous_validated_minimum(directly_testable)
    validated_p_max = (
        max(float(row["nominal_p"]) for row in directly_testable)
        if validated_p_min is not None
        else None
    )
    accepted = stability_pass and validation_pass

    model = {
        "format": MODEL_SCHEMA,
        "kernel_schema": fit_metadata["kernel_schema"],
        "estimator": fit_metadata["estimator"],
        "sampling_null": "independent-uniform-rank-permutation",
        "rank_policy": "unique-ordinal-ranks",
        "m": m,
        "npar_limit": args.npar,
        "tail_model": "generalized-pareto-mle-floc0",
        "tail_threshold_quantile": primary_fit.threshold_quantile,
        "tail_threshold": primary_fit.threshold,
        "tail_probability": primary_fit.tail_probability,
        "tail_shape": primary_fit.shape,
        "tail_scale": primary_fit.scale,
        "tail_endpoint": primary_fit.endpoint,
        "calibration_status": "accepted" if accepted else "rejected",
        "calibrator_schema": CALIBRATOR_SCHEMA,
        "calibrator_sha256": calibrator_sha256(),
        "validation_method": "independent-rank-permutation-stream",
        "stability_probability": stability_probability,
        "stability_relative_range": relative_stability_range,
        "stability_relative_tolerance": args.stability_relative_tolerance,
        "validation_family_confidence": args.validation_family_confidence,
        "validation_point_confidence": point_confidence,
        "supported_p_min": args.p_min,
        "supported_p_max": min(args.p_max, primary_fit.tail_probability),
        "validated_p_min": validated_p_min,
        "validated_p_max": validated_p_max,
        "default_p": args.default_p,
        "default_p_cutoff": gpd_cutoff(primary_fit, args.default_p),
        "fit_draws": args.fit_draws,
        "validation_draws": args.validation_draws,
        "fit_seed": fit_seed,
        "validation_seed": validation_seed,
        "rng": fit_metadata["rng"],
        "generator_sha256": sha256(args.generator),
        "fit_values_sha256": sha256(fit_path),
        "validation_values_sha256": sha256(validation_path),
        "scipy_version": scipy.__version__,
    }
    # A rejected fit is diagnostic evidence, not a usable runtime model.  The
    # override controls the process exit status only; it never emits a model that
    # could be mistaken for accepted calibration.
    if accepted:
        write_model(model_path, model)
    else:
        model_path = None

    return {
        "m": m,
        "accepted": accepted,
        "model_path": None if model_path is None else str(model_path),
        "model": model,
        "primary_fit": asdict(primary_fit),
        "stability_fits": [asdict(fit) for fit in stability_fits],
        "skipped_stability_quantiles": skipped_stability_quantiles,
        "relative_stability_range": relative_stability_range,
        "stability_probability": stability_probability,
        "stability_pass": stability_pass,
        "direct_validation_floor": directly_testable_floor,
        "validation_family_confidence": args.validation_family_confidence,
        "validation_point_confidence": point_confidence,
        "validation_family_comparisons": total_family_comparisons,
        "validation_pass": validation_pass,
        "validation": validation,
        "fit_summary": {
            "minimum": float(np.min(fit_values)),
            "maximum": float(np.max(fit_values)),
            "mean": float(np.mean(fit_values)),
            "zero_fraction": float(np.mean(fit_values == 0.0)),
            "unique_values": int(np.unique(fit_values).size),
        },
        "fit_metadata": fit_metadata,
        "validation_metadata": validation_metadata,
    }


def main() -> int:
    try:
        args = parse_args()
        validate_args(args)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        results = [calibrate_one(args, m) for m in sorted(args.m)]
    except (OSError, RuntimeError, ValueError) as error:
        print("AP-MI null calibration failed: {}".format(error), file=sys.stderr)
        return 1

    report = {
        "format": "sjaracne-apmi-null-calibration-report-v1",
        "command": sys.argv,
        "cwd": os.getcwd(),
        "generator": str(args.generator.resolve()),
        "generator_sha256": sha256(args.generator),
        "tail_design": (
            "Predeclared peaks-over-threshold generalized Pareto fit at the "
            "requested quantile. Fit and validation use independent deterministic "
            "permutation streams. The default p=1e-7 is extrapolated and is not "
            "directly resolved by the Monte Carlo sample."
        ),
        "results": results,
    }
    report_path = args.output_dir / "calibration_report.json"
    write_text_lf(report_path, json.dumps(report, indent=2, allow_nan=False) + "\n")

    summary_path = args.output_dir / "calibration_summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "m",
            "npar_limit",
            "accepted",
            "fit_draws",
            "validation_draws",
            "shape",
            "scale",
            "tail_threshold",
            "tail_probability",
            "cutoff_p1e-7",
            "relative_stability_range",
            "validation_pass",
            "validated_p_min",
            "validated_p_max",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            if "model" not in result:
                writer.writerow(
                    {
                        "m": result["m"],
                        "accepted": False,
                    }
                )
                continue
            model = result["model"]
            writer.writerow(
                {
                    "m": result["m"],
                    "npar_limit": model["npar_limit"],
                    "accepted": result["accepted"],
                    "fit_draws": model["fit_draws"],
                    "validation_draws": model["validation_draws"],
                    "shape": model["tail_shape"],
                    "scale": model["tail_scale"],
                    "tail_threshold": model["tail_threshold"],
                    "tail_probability": model["tail_probability"],
                    "cutoff_p1e-7": model["default_p_cutoff"],
                    "relative_stability_range": result["relative_stability_range"],
                    "validation_pass": result["validation_pass"],
                    "validated_p_min": model["validated_p_min"],
                    "validated_p_max": model["validated_p_max"],
                }
            )

    failures = [result["m"] for result in results if not bool(result["accepted"])]
    print(f"wrote {report_path}")
    print(f"wrote {summary_path}")
    if failures and not args.allow_rejected_calibration:
        print(f"calibration rejected for m={failures}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
