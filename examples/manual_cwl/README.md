# Manual local CWL examples — not automated regression tests

These retained launchers run commands; they are **not asserting tests** and are
never launched by `tests/run_regression.py`. They were moved out of `tests/` to
make that distinction explicit. Their BRCA inputs remain in `tests/inputs/`
because the main README also uses those bundled demonstration inputs.

Run a local example explicitly with `bash examples/manual_cwl/cwltool/<name>.sh`.
The scripts locate the repository root themselves. Install the package first.
The job templates use explicit historical affine/consensus probabilities; they
are compatibility examples, not recommended current PCHIP/K=6 settings. The
workflow template uses two resamples and explicitly selects legacy p-value
consensus; it does not test the default K=6 behavior.

The obsolete CWLEXEC cluster launchers and their configuration templates have
been removed. These remaining examples use local `cwltool` only.
**Do not include these launchers in routine CI.** Current parameter examples are in the main
[README](../../README.md); small asserting tests are described in
[tests/README.md](../../tests/README.md).
