"""Shared scientific profiles for the public CLI and WDL task preflight.

These are workflow defaults, not changes to the native executable's low-level
defaults. Engine choice never selects the scientific profile implicitly.
"""
import math
import re


OPTION_NAMES = {
    "profile", "network_type", "sampling_method", "bootstrap_num", "subsample_spec", "depth",
    "p_value_bootstrap", "apmi_null_model", "min_recurrence", "p_value_consensus", "tolerance",
    "dynamic_dpi", "dpi_href", "dpi_slope", "dpi_max",
}
NETWORK_PROBABILITIES = {"TF": 0.001, "SIG": 0.0005}


def integer(value, name, minimum=1, maximum=2147483647):
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ValueError("{} must be an integer within [{},{}]".format(name, minimum, maximum))
    return value


def finite_number(value, name, lower, upper, lower_inclusive=True, upper_inclusive=True):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(name + " must be numeric")
    value = float(value)
    if (not math.isfinite(value)
            or (value < lower if lower_inclusive else value <= lower)
            or (value > upper if upper_inclusive else value >= upper)):
        raise ValueError(name + " is outside its supported range")
    return value


def selected_count(spec, observation_count):
    if not isinstance(spec, str):
        raise ValueError("subsample_spec must be a count or explicit percentage string")
    if re.fullmatch(r"[0-9]+", spec):
        m = int(spec)
    elif re.fullmatch(r"(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)%", spec):
        percentage = float(spec[:-1])
        if not 0 < percentage <= 100:
            raise ValueError("subsample percentage must be within (0,100]")
        # Match the exact operation order in Parameter::resolveSubsampleSize.
        m = math.ceil(percentage * observation_count / 100.0)
    else:
        raise ValueError("subsample_spec must be a count or explicit percentage string")
    if not 2 <= m <= observation_count:
        raise ValueError("selected sample count must be within [2,N]")
    return m


def resolve_options(parameters, observation_count):
    """Resolve omissions once and reject conflicts before launching any runs.

The legacy profile restores pre-sampling workflow settings, not an obsolete
binary. Explicit historical pb/pc/depth/fixed-tolerance tuning remains possible.
"""
    if not isinstance(parameters, dict) or set(parameters) - OPTION_NAMES:
        raise ValueError("Unknown workflow parameters")
    values = {key: value for key, value in parameters.items() if value is not None}
    profile = values.get("profile", "current")
    if profile not in ("current", "legacy"):
        raise ValueError("profile must be current or legacy")
    network = values.get("network_type")
    if isinstance(network, str):
        network = network.upper()
    if network not in (None, "TF", "SIG"):
        raise ValueError("network_type must be TF or SIG")
    if profile == "current" and network is None:
        raise ValueError("The current profile requires network_type TF or SIG, even with explicit pb")
    n = integer(values.get("bootstrap_num", 100), "bootstrap_num")
    depth = integer(values.get("depth", 40), "depth")
    method = values.get("sampling_method", "with-replacement" if profile == "legacy" else "without-replacement")
    if method not in ("with-replacement", "without-replacement"):
        raise ValueError("sampling_method must be with-replacement or without-replacement")
    model = values.get("apmi_null_model", "affine" if profile == "legacy" else "pchip")
    if model not in ("affine", "pchip"):
        raise ValueError("apmi_null_model must be affine or pchip")
    dynamic = values.get("dynamic_dpi", profile == "current" and "tolerance" not in values)
    if not isinstance(dynamic, bool):
        raise ValueError("dynamic_dpi must be boolean")
    if dynamic and "tolerance" in values:
        raise ValueError("dynamic_dpi conflicts with an explicit fixed tolerance")
    if profile == "legacy":
        if method != "with-replacement" or "subsample_spec" in values:
            raise ValueError("The legacy profile requires full-size sampling with replacement; omit subsampling")
        if model != "affine" or dynamic:
            raise ValueError("The legacy profile requires affine MI cutoffs and fixed DPI")
        if "min_recurrence" in values:
            raise ValueError("The legacy profile uses p_value_consensus, not minimum recurrence")
    if method == "with-replacement":
        if "subsample_spec" in values:
            raise ValueError("Replacement sampling conflicts with subsample_spec")
        spec, m = None, observation_count
        if observation_count < 2:
            raise ValueError("At least two observations are required")
    else:
        spec = values.get("subsample_spec", "80%")
        m = selected_count(spec, observation_count)
    pb_default = 1e-7 if profile == "legacy" else NETWORK_PROBABILITIES[network]
    pb = finite_number(values.get("p_value_bootstrap", pb_default), "p_value_bootstrap", 0, 1,
                       lower_inclusive=False)
    if model == "pchip" and pb != 1:
        if method != "without-replacement":
            raise ValueError("PCHIP requires sampling without replacement; choose affine for replacement sampling")
        if depth != 40:
            raise ValueError("PCHIP requires depth 40")
        from .apmi_cutoff import get_mi_cutoff
        get_mi_cutoff(m, pb)
    k, pc = values.get("min_recurrence"), values.get("p_value_consensus")
    if k is not None and pc is not None:
        raise ValueError("min_recurrence and p_value_consensus are mutually exclusive")
    if profile == "legacy" and pc is None:
        pc = 1e-5
    if pc is None:
        k = integer(6 if k is None else k, "min_recurrence")
        if k > n:
            raise ValueError("minimum recurrence exceeds bootstrap_num")
    else:
        pc = finite_number(pc, "p_value_consensus", 0, 1, lower_inclusive=False)
    tuning = {key: values[key] for key in ("dpi_href", "dpi_slope", "dpi_max") if key in values}
    if tuning and not dynamic:
        raise ValueError("DPI tuning parameters require dynamic_dpi")
    if "dpi_href" in tuning:
        integer(tuning["dpi_href"], "dpi_href")
    if "dpi_slope" in tuning:
        finite_number(tuning["dpi_slope"], "dpi_slope", 0, math.inf, lower_inclusive=False)
    if "dpi_max" in tuning:
        finite_number(tuning["dpi_max"], "dpi_max", 0, 1, upper_inclusive=False)
    if not dynamic:
        values["tolerance"] = finite_number(values.get("tolerance", 0), "tolerance", 0, 1)
    values.update(profile=profile, network_type=network, sampling_method=method,
                  bootstrap_num=n, depth=depth, subsample_spec=spec,
                  p_value_bootstrap=pb, dynamic_dpi=dynamic,
                  min_recurrence=k, p_value_consensus=pc, apmi_null_model=model)
    return values
