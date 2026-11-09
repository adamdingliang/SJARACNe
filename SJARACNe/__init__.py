"""SJARACNe public API; the cutoff helper needs only the Python standard library."""

from .apmi_cutoff import get_mi_cutoff

__all__ = [
    'get_mi_cutoff', 'ch_line_ending', 'os', 'sys', 'argparse', 'main',
    'create_consensus_network', 'math', 'logging', 'np', 'pathlib', 're',
    'stats', 'pd', 'create_enhanced_consensus_network', 'uprob',
    'numbers', 'DEFAULT_MIN_RECURRENCE', 'minimum_recurrence',
    'consensus_probability',
]


def __dir__():
    return sorted(set(globals()) | set(__all__))


def __getattr__(name):
    # Preserve access to historical helper exports without importing SciPy/pandas
    # when callers only need the bundled MI-cutoff function.
    from importlib import import_module
    helpers = import_module('.bin', __name__)
    try:
        return getattr(helpers, name)
    except AttributeError:
        raise AttributeError("module {!r} has no attribute {!r}".format(__name__, name)) from None
