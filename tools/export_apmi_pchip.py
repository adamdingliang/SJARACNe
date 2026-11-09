#!/usr/bin/env python3
"""Mechanically export the frozen PCHIP surface; never recalibrate or select.

SciPy/NumPy are offline export dependencies only. The input is an already
selected, immutable model and its immutable summary, not null draws. Output
coefficients describe exactly the existing first-stage interpolation curves.
The Python callable and native executable consume only the generated data.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import scipy
from scipy.interpolate import PchipInterpolator


FORMAT = "sjaracne-apmi-pchip-data-v1"
STATUS = "experimental-not-uniformly-validated"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_locked(path):
    path = Path(path).resolve(strict=True)
    digest = sha256(path)
    if path.with_suffix(".sha256").read_text(encoding="ascii").strip() != digest:
        raise ValueError("frozen input hash mismatch: " + str(path))
    value = json.loads(path.read_text(encoding="utf-8"),
                       parse_constant=lambda text: (_ for _ in ()).throw(ValueError("nonfinite JSON: " + text)))
    if not isinstance(value, dict):
        raise ValueError("frozen model inputs must be JSON objects")
    return value, digest


def export_data(primary_path, summary_path):
    primary, primary_hash = read_locked(primary_path)
    summary, summary_hash = read_locked(summary_path)
    if (primary.get("format") != "apmi-two-surface-frozen-primary-v1"
            or primary.get("selected_model") != "pchip" or summary.get("primary") != primary
            or summary.get("selection", {}).get("selected_model") != "pchip"):
        raise ValueError("inputs do not identify the same frozen PCHIP primary")
    model = primary.get("surface", {}).get("fit_models", {})
    if model.get("schema") != "sjaracne-apmi-empirical-surface-v1" or model.get("method") != "pchip":
        raise ValueError("unsupported frozen PCHIP surface schema")
    expected_m = [m for m in range(50, 10001, 50) if m % 500 != 350]
    ms = model.get("train_m")
    ps = model.get("p_knots")
    columns = model.get("columns")
    if (ms != expected_m or primary.get("training_m") != expected_m or not isinstance(ps, list)
            or len(ps) != 25 or len(columns or []) != 25
            or any(not isinstance(p, (int, float)) or isinstance(p, bool) or not math.isfinite(p) for p in ps)
            or any(a <= b for a, b in zip(ps, ps[1:]))
            or ps[0] != 0.05 or ps[-1] != 1e-7 or sum(p >= 0.0001 for p in ps) != 22):
        raise ValueError("frozen dimensions or construction knots differ from the 180x25 primary")
    log_min = math.log(50.0)
    log_width = math.log(10000.0) - log_min
    xs = (np.log(np.asarray(ms, dtype=np.float64)) - log_min) / log_width
    observed, log_observed, coefficients = [], [], []
    for p, column in zip(ps, columns):
        obs, logs = column.get("observed"), column.get("log_values")
        hexadecimal = column.get("observed_hex")
        if (column.get("p") != p or column.get("available") is not True
                or len(obs or []) != 180 or len(logs or []) != 180 or len(hexadecimal or []) != 180):
            raise ValueError("frozen PCHIP column is unavailable or malformed")
        if any(not math.isfinite(q) or not 0 < q <= math.log(m) or q.hex() != text
               for q, text, m in zip(obs, hexadecimal, ms)):
            raise ValueError("literal frozen quantile or its hexadecimal representation is invalid")
        if not np.all(np.isfinite(logs)):
            raise ValueError("nonfinite frozen log-cutoff")
        observed.append(list(obs))
        log_observed.append(list(logs))
        coefficients.append(PchipInterpolator(xs, np.asarray(logs, dtype=np.float64), extrapolate=False).c.tolist())
    result = {
        "format": FORMAT, "status": STATUS,
        "kernel_schema": "sjaracne-apmi-v1", "estimator": "sjaracne-adaptive-partitioning",
        "sampling_null": "independent-uniform-rank-permutation", "rank_policy": "unique-ordinal-ranks",
        "npar": 40, "m_min": 50, "m_max": 10000,
        "moderate_p_min": 0.0005, "p_max": 0.05, "diagnostic_p_min": 1e-7,
        "scope_boundary": 0.0001, "moderate_knot_count": 22,
        "numeric_policy": "binary64-m-first-p-second-literal-empirical-knots-v1",
        "LOG_M_MIN": log_min, "LOG_M_WIDTH": log_width,
        "NM": 180, "NP": 25, "M": expected_m, "X": xs.tolist(), "P": list(ps),
        "LOG_P": np.log(np.asarray(ps, dtype=np.float64)).tolist(),
        "OBSERVED": observed, "LOG_OBSERVED": log_observed, "COEFF": coefficients,
        "PRIMARY_SHA256": primary_hash, "SUMMARY_SHA256": summary_hash,
        "generation_plan_sha256": primary["generation_plan_sha256"],
        "source_repository": "https://github.com/adamdingliang/sjaracne_benchmark",
        "source_artifact_directory": "benchmarks/apmi_null_m_sweep/broad_grid_10m/results_2026-10-17",
        "export_numpy_version": np.__version__, "export_scipy_version": scipy.__version__,
        "exporter_sha256": sha256(__file__),
        "no_new_draws_fitting_selection_or_recalibration": True,
    }
    return result


def _cpp_number(value):
    text = format(float(value), ".17g")
    return text if "." in text or "e" in text.lower() else text + ".0"


def _initializer(values, integers=False, level=0):
    if not isinstance(values[0], list):
        pieces = [str(value) if integers else _cpp_number(value) for value in values]
        lines = [", ".join(pieces[start:start + 6]) for start in range(0, len(pieces), 6)]
        return "{\n" + ",\n".join("   " * (level + 1) + line for line in lines) + "\n" + "   " * level + "}"
    return "{\n" + ",\n".join("   " * (level + 1) + _initializer(row, integers, level + 1) for row in values) + "\n" + "   " * level + "}"


def header(data, data_hash):
    lines = ["// Generated by tools/export_apmi_pchip.py from frozen model parameters.",
             "// Offline representation export only; no raw data, fitting, or recalibration.",
             "#ifndef SJARACNE_APMI_PCHIP_DATA_H", "#define SJARACNE_APMI_PCHIP_DATA_H", "",
             "namespace ApmiPchipData {", "static const int NM = 180;", "static const int NP = 25;"]
    for key in ("PRIMARY_SHA256", "SUMMARY_SHA256"):
        lines.append('static const char ' + key + '[] = "' + data[key] + '";')
    lines.append('static const char DATA_SHA256[] = "' + data_hash + '";')
    for name in ("LOG_M_MIN", "LOG_M_WIDTH"):
        lines.append("static const double " + name + " = " + _cpp_number(data[name]) + ";")
    dimensions = {"M": "[NM]", "X": "[NM]", "P": "[NP]", "LOG_P": "[NP]",
                  "OBSERVED": "[NP][NM]", "LOG_OBSERVED": "[NP][NM]", "COEFF": "[NP][4][NM-1]"}
    for name, size in dimensions.items():
        lines.append("static const " + ("int" if name == "M" else "double") + " " + name + size + " = "
                     + _initializer(data[name], integers=name == "M") + ";")
    lines.extend(("} // namespace ApmiPchipData", "", "#endif", ""))
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--out-root", required=True, type=Path)
    args = parser.parse_args()
    root = args.out_root.resolve(strict=True)
    if not (root / "SJARACNe" / "src").is_dir() or not (root / "SJARACNe" / "config").is_dir():
        raise ValueError("out-root must be the intended SJARACNe source checkout")
    json_path = root / "SJARACNe" / "config" / "apmi_pchip.json"
    header_path = root / "SJARACNe" / "src" / "apmi_pchip_data.h"
    if json_path.exists() or header_path.exists():
        raise ValueError("refuse to overwrite exported model data")
    data = export_data(args.primary, args.summary)
    payload = (json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("ascii")
    data_hash = hashlib.sha256(payload).hexdigest()
    with json_path.open("xb") as handle:
        handle.write(payload)
    with header_path.open("x", encoding="ascii", newline="\n") as handle:
        handle.write(header(data, data_hash))
    print(json.dumps({"json": str(json_path), "header": str(header_path), "DATA_SHA256": data_hash,
                      "PRIMARY_SHA256": data["PRIMARY_SHA256"], "SUMMARY_SHA256": data["SUMMARY_SHA256"]}, indent=2))


if __name__ == "__main__":
    main()
