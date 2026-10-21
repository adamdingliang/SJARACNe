# Manual CWL and LSF examples — not automated regression tests

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

The `cwlexec/` examples additionally require a functioning LSF cluster and
site-specific queue configuration. The separate `../lsf/run_real_data.sh` can
submit 100 resampling jobs and requires explicit input environment variables.
**Do not include these launchers in routine CI.** LSF execution was not validated
by the local regression cleanup. Current parameter examples are in the main
[README](../../README.md); small asserting tests are described in
[tests/README.md](../../tests/README.md).
