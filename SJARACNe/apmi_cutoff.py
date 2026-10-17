"""Dependency-free callable for the frozen, experimental AP-MI PCHIP surface.

The runtime reads only bundled, hash-pinned model parameters. It never fits,
selects, simulates, calls SciPy or imports the private benchmark repository.
"""
from __future__ import annotations

from bisect import bisect_right
from functools import lru_cache
import hashlib
import json
import math
from numbers import Integral, Real
from pathlib import Path
import warnings


__all__ = ["get_mi_cutoff"]
_DATA_SHA256 = "5ae44f3ec251282b242e1e25e5695fe8c6d80a0e7a3f156170dbb72c29e09f59"
_DATA_PATH = Path(__file__).resolve().parent / "config" / "apmi_pchip.json"


@lru_cache(maxsize=1)
def _load_data():
    payload = _DATA_PATH.read_bytes()
    if hashlib.sha256(payload).hexdigest() != _DATA_SHA256:
        raise ValueError("Bundled AP-MI PCHIP data hash mismatch; reinstall the matching package.")
    result = json.loads(payload.decode("ascii"))
    if (result.get("format") != "sjaracne-apmi-pchip-data-v1"
            or result.get("status") != "experimental-not-uniformly-validated"
            or result.get("NM") != 180 or result.get("NP") != 25 or result.get("npar") != 40
            or len(result.get("M", [])) != 180 or len(result.get("X", [])) != 180
            or len(result.get("P", [])) != 25 or len(result.get("LOG_P", [])) != 25
            or any(len(result.get(name, [])) != 25 for name in ("OBSERVED", "LOG_OBSERVED", "COEFF"))):
        raise ValueError("Unsupported bundled AP-MI PCHIP data schema.")
    return result


def _sign(value):
    return 1 if value > 0 else (-1 if value < 0 else 0)


def _endpoint(h0, h1, slope0, slope1):
    value = ((2 * h0 + h1) * slope0 - h0 * slope1) / (h0 + h1)
    if _sign(value) != _sign(slope0):
        return 0.0
    if _sign(slope0) != _sign(slope1) and abs(value) > 3 * abs(slope0):
        return 3 * slope0
    return value


def _pchip_value(xs, ys, query):
    """Fritsch-Butland PCHIP with SciPy-compatible endpoint and cubic rules."""
    if len(xs) != len(ys) or len(xs) < 2 or not xs[0] <= query <= xs[-1]:
        raise ValueError("PCHIP interpolation cannot extrapolate.")
    h = [b - a for a, b in zip(xs, xs[1:])]
    if any(not math.isfinite(value) or value <= 0 for value in h):
        raise ValueError("PCHIP coordinates must increase strictly.")
    slopes = [(b - a) / width for a, b, width in zip(ys, ys[1:], h)]
    if len(xs) == 2:
        derivatives = [slopes[0], slopes[0]]
    else:
        derivatives = [_endpoint(h[0], h[1], slopes[0], slopes[1])]
        for i in range(1, len(xs) - 1):
            left, right = slopes[i - 1], slopes[i]
            if left == 0 or right == 0 or _sign(left) != _sign(right):
                derivatives.append(0.0)
            else:
                w1, w2 = 2 * h[i] + h[i - 1], h[i] + 2 * h[i - 1]
                weighted = (w1 / left + w2 / right) / (w1 + w2)
                derivatives.append(1.0 / weighted)
        derivatives.append(_endpoint(h[-1], h[-2], slopes[-1], slopes[-2]))
    index = min(bisect_right(xs, query) - 1, len(xs) - 2)
    width = h[index]
    slope = slopes[index]
    t = (derivatives[index] + derivatives[index + 1] - 2 * slope) / width
    coefficients = (t / width, (slope - derivatives[index]) / width - t,
                    derivatives[index], ys[index])
    dx = query - xs[index]
    return ((coefficients[0] * dx + coefficients[1]) * dx + coefficients[2]) * dx + coefficients[3]


def _predict(data, m, p):
    training = data["M"]
    exact = bisect_right(training, m) - 1
    at_training = exact >= 0 and training[exact] == m
    if at_training:
        knots = [column[exact] for column in data["OBSERVED"]]
        logs = [column[exact] for column in data["LOG_OBSERVED"]]
    else:
        x = (math.log(m) - data["LOG_M_MIN"]) / data["LOG_M_WIDTH"]
        index = bisect_right(data["X"], x) - 1
        if index < 0 or index >= data["NM"] - 1:
            raise ValueError("PCHIP sample-count interpolation cannot extrapolate.")
        dx = x - data["X"][index]
        knots, logs = [], []
        for coefficients in data["COEFF"]:
            log_value = ((coefficients[0][index] * dx + coefficients[1][index]) * dx + coefficients[2][index]) * dx + coefficients[3][index]
            value = math.exp(log_value)
            knots.append(value)
            # Preserve the original exp -> log operation between the two axes.
            logs.append(math.log(value))
    indices = list(range(22)) if p >= data["scope_boundary"] else list(range(21, 25))
    selected = [knots[i] for i in indices]
    if any(not math.isfinite(value) or not 0 < value <= math.log(m) for value in selected):
        raise ValueError("PCHIP produced an invalid cutoff within the requested probability scope.")
    if any(a > b for a, b in zip(selected, selected[1:])):
        raise ValueError("PCHIP cutoff knots cross; no repair or sorting is permitted.")
    if p in data["P"]:
        # Exact-p return avoids a second log/exp round trip at discrete MI atoms.
        value = knots[data["P"].index(p)]
    else:
        xs = [data["LOG_P"][i] for i in reversed(indices)]
        ys = [logs[i] for i in reversed(indices)]
        value = math.exp(_pchip_value(xs, ys, math.log(p)))
    if not math.isfinite(value) or not 0 < value <= math.log(m):
        raise ValueError("PCHIP produced a nonfinite, nonpositive or out-of-support cutoff.")
    return value


def get_mi_cutoff(m, p, npar=40):
    """Return the experimental frozen PCHIP MI cutoff for actual sample count m.

    m must be an integer in [50,10000], npar must be 40, and the required p
    must lie in [1e-7,0.05]. This full computational range needs no additional
    acknowledgment flag. The principal assessment range [0.0005,0.05] is an
    evidence/reporting distinction, not a computational failure boundary.
    Off-grid m, extreme-tail Monte Carlo precision and incomplete uniform
    calibration are warned explicitly.
    No legacy/GPD fallback, extrapolation or clipping is performed.
    """
    if isinstance(m, bool) or not isinstance(m, Integral) or not 50 <= int(m) <= 10000:
        raise ValueError("m must be the actual integer sample count within [50,10000].")
    if isinstance(npar, bool) or not isinstance(npar, Integral) or int(npar) != 40:
        raise ValueError("The frozen AP-MI PCHIP surface requires npar=40.")
    if isinstance(p, bool) or not isinstance(p, Real) or not math.isfinite(float(p)) or not 1e-7 <= float(p) <= 0.05:
        raise ValueError("p must be finite and within [1e-7,0.05].")
    m, p = int(m), float(p)
    data = _load_data()
    warnings.warn("The frozen AP-MI PCHIP surface is experimental and not uniformly validated; it is not a biological false-positive guarantee.",
                  RuntimeWarning, stacklevel=2)
    if m % 50:
        warnings.warn("This off-grid m was not directly confirmed in the broad-grid study.", RuntimeWarning, stacklevel=2)
    if p < 0.0005:
        warnings.warn("This p is outside the principal assessment range [0.0005,0.05]; this is an evidence/reporting distinction, not a computational failure boundary.",
                      RuntimeWarning, stacklevel=2)
    if p <= 1e-5:
        warnings.warn("Extreme-tail Monte Carlo precision is limited: 10M validation draws expect 100, 10 and 1 exceedances at p=1e-5, 1e-6 and 1e-7. A returned cutoff does not establish 10%-accuracy tail calibration.",
                      RuntimeWarning, stacklevel=2)
    return _predict(data, m, p)
