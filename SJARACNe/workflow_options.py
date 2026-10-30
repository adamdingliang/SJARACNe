"""Resolve scientific profiles shared by the CLI and WDL task preflight.

These profiles supply workflow settings, not new defaults for direct native
execution. Engine selection and scheduler resources are configured elsewhere;
neither selects a scientific profile implicitly.

Attributes:
    OPTION_NAMES (set[str]): Scientific keys accepted in profile requests.
    NETWORK_PROBABILITIES (dict[str, float]): Current-profile default null-tail
        probabilities for TF and SIG networks, not biological FDR estimates.
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
    """Validate an integer without accepting booleans or coercing other types.

    Args:
        value (int): Candidate value.
        name (str): Parameter name used in diagnostics.
        minimum (int): Inclusive lower bound.
        maximum (int): Inclusive upper bound.

    Returns:
        int: The unchanged value.

    Raises:
        ValueError: If the value is not a non-boolean integer within the bounds.
    """
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ValueError("{} must be an integer within [{},{}]".format(name, minimum, maximum))
    return value


def finite_number(value, name, lower, upper, lower_inclusive=True, upper_inclusive=True):
    """Convert an integer or float to a finite float within explicit bounds.

    Args:
        value (int | float): Candidate value; booleans are rejected.
        name (str): Parameter name used in diagnostics.
        lower (float): Lower bound.
        upper (float): Upper bound.
        lower_inclusive (bool): Whether equality with the lower bound is valid.
        upper_inclusive (bool): Whether equality with the upper bound is valid.

    Returns:
        float: The validated numeric value.

    Raises:
        ValueError: If the value has an unsupported type, is non-finite, or lies
            outside the requested interval.
    """
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(name + " must be numeric")
    value = float(value)
    if (not math.isfinite(value)
            or (value < lower if lower_inclusive else value <= lower)
            or (value > upper if upper_inclusive else value >= upper)):
        raise ValueError(name + " is outside its supported range")
    return value


def selected_count(spec, observation_count):
    """Resolve a count or percentage for sampling without replacement.

    Percentage requests use the native estimator's ceiling rule and floating-
    point operation order so that preflight checks use the same sample count.

    Args:
        spec (str): Integer count such as ``"80"`` or explicit percentage such
            as ``"80%"``. A bare decimal fraction is not a percentage request.
        observation_count (int): Number of observations in the input matrix.

    Returns:
        int: Selected observation count within ``[2, observation_count]``.

    Raises:
        ValueError: If the specification is malformed, its percentage is outside
            ``(0, 100]``, or the resulting count is outside the allowed range.
    """
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
    """Resolve profile defaults and reject incompatible scientific settings.

    The current profile requires an explicit TF/SIG network type and defaults
    to 80% sampling without replacement, PCHIP cutoffs, dynamic DPI and K=6.
    The legacy profile defaults to full-size sampling with replacement, affine
    cutoffs with ``pb=1e-7``, strict DPI and consensus ``pc=1e-5``. It selects
    retained compatibility behavior in the current executable, not an old
    executable or a different workflow engine.

    Explicit compatible overrides are retained. In particular, historical
    probability, depth and fixed-tolerance tuning is allowed in the legacy
    profile; mixing that profile with modern sampling, PCHIP, dynamic DPI or
    recurrence-based consensus is rejected.

    Args:
        parameters (dict[str, object]): Scientific request using only
            ``OPTION_NAMES``. Omitted keys and values of ``None`` select profile
            defaults. Engine and resource settings are not accepted here.
        observation_count (int): Full input observation count, used to resolve
            the selected count and validate estimator support.

    Returns:
        dict[str, object]: A new mapping with resolved scientific settings.
            Exactly one of ``min_recurrence`` and ``p_value_consensus`` is active;
            replacement sampling has a null ``subsample_spec``. Optional
            dynamic-DPI tuning keys remain absent unless explicitly supplied.

    Raises:
        ValueError: If keys, values, profile combinations, sampling counts or
            estimator requests are unsupported or inconsistent.

    Note:
        PCHIP support is checked at the actual selected sample count; the cutoff
        itself is not stored in the returned mapping. Its experimental and
        off-grid warnings may be emitted during this check.
    """
    if not isinstance(parameters, dict) or set(parameters) - OPTION_NAMES:
        raise ValueError("Unknown workflow parameters")
    # WDL optional inputs arrive as null; treat them like omitted CLI options.
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
    # A legacy label must not silently combine incompatible scientific modes.
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
    # pb=1 bypasses the null-cutoff lookup in the native executable as well.
    if model == "pchip" and pb != 1:
        if method != "without-replacement":
            raise ValueError("PCHIP requires sampling without replacement; choose affine for replacement sampling")
        if depth != 40:
            raise ValueError("PCHIP requires depth 40")
        from .apmi_cutoff import get_mi_cutoff
        get_mi_cutoff(m, pb)
    # Resolve consensus as either a recurrence count or a legacy probability.
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
