# Small workflow regression tests

The suite checks functional correctness with small fixtures. It does not run
BRCA100/DREAM benchmark batches or require a cluster.

## Complete check (Linux or WSL)

Use Python **3.11 or 3.14**, uv **0.11.16**, Make, a C++ compiler, Node.js,
`cwltool`, Bats and **Sprocket 0.32.0**. Install this checkout in the active
environment so CWL tasks import the revision being tested:

```bash
# For the compatibility lane, change both Python selectors to 3.11.
uv sync --locked --python 3.14 --group regression --extra cwl
uv run --locked --python 3.14 --no-sync python tests/run_regression.py --jobs 4
```

Keep that `--python` selection on later `uv run` commands, or set
`UV_PYTHON=3.11` for a compatibility-lane session.

On Ubuntu, install system tools with
`sudo apt-get install build-essential nodejs bats`. The regression group and
CWL extra supply Python tools; install Sprocket separately using the
[workflow guide](../SJARACNe/wdl/README.md#local-use).
`--locked` preserves the committed dependency resolution. The runner installs
nothing and downloads neither data nor dependencies.

Prepare a dedicated environment before submitting jobs. Do not synchronize it
while other jobs use it. Environment placement and deployment limits are in
[the environment guide](../docs/uv.md#hpc-and-cloud-precautions).

The runner prints a fresh evidence directory under the system temp directory.
To choose a new directory:

```bash
uv run --locked --no-sync python tests/run_regression.py --jobs 4 --output-dir /tmp/sjaracne-regression-check
```

It checks dependencies and installed-runtime bytes, force-rebuilds both native
programs, runs all `test_*.py` tests, validates production CWL, runs Bats utilities
and executes real default workflows plus WDL/CWL parity. Missing dependencies,
empty discovery, failures, skips, expected failures or changed source/fixtures
prevent a successful result. Logs, versions and source/native hashes are retained.

## Coverage and fixtures

| Coverage | What is checked |
| --- | --- |
| Scientific compatibility | Tiny expression/hub fixtures, sampling, MI conversion, DPI, explicit K6 and historical answerkeys. |
| Current defaults | Real Sprocket TF/SIG runs with twenty resamples and omitted MI/DPI/K options; checks native consensus-first order, default K20, role-specific probabilities, outputs and provenance. |
| Legacy profile | Replacement sampling, affine MI, strict DPI and probability consensus through WDL and CWL. |
| Native inference | Seeds, rank caches, AP-MI workspace/null generation, PCHIP bounds/warnings, expression/hub validation, adjacency replay and alias protection. |
| Dynamic DPI | Formula, offset, hub count including empty rows, cap, fixed-epsilon equivalence, conflicting options and witness diagnostics. |
| Native consensus | Inclusive support, supporting-run mean MI, direction/zero MI, deterministic manifest order, one-file-at-a-time reading, malformed graphs, allocation failures and atomic output publication. |
| Consensus and enhancement | K20 boundaries, duplicate handling, explicit smaller K and probability selection, missing/empty inputs, symbols and valid statistics. |
| Workflow plumbing | Option forwarding, literal argv, line-ending copies, runtime/input/settings provenance, default and opt-in raw-graph hash verification. |
| Utilities | CWL array types/order and isolated Bats fixtures with byte-exact output checks. |
| Packaging | Package-data contracts, fresh archives/wheels/editable installs, correct build failures and installed runtime assets. |
| uv deployment | Disposable editable/non-editable installations, warmed-cache offline reuse, native/configuration refresh and missing-tool errors. |

The T-cell `cnn_5.txt` network, legacy answerkeys and rank-cache inputs are
compatibility fixtures. Bundled BRCA100 inputs support manual examples; routine
tests do not run the full resampling examples.

For a focused check:

```bash
python -m unittest discover -s tests -p 'test_create_consensus_network.py' -v
```

Direct discovery can skip unavailable native programs. Set
`SJARACNE_REQUIRE_NATIVE_TESTS=1` to fail instead. Invalid explicit
`SJARACNE_TEST_EXE` or `SJARACNE_APMI_NULL_EXE` paths always fail. The complete
runner pins both executables and the installed runtime to this checkout;
`PYTHONPATH` alone does not control imports in CWL subprocesses.

After building the native binary and installing this checkout, dynamic-DPI
checks can also be run directly:

```bash
python -m unittest discover -s tests -p 'test_sjaracne_dynamic_dpi.py' -v
python -m unittest discover -s tests -p 'test_dynamic_dpi_workflow.py' -v
```

## Removed or separated during maintenance

Unused cluster logs/templates, non-asserting launchers, stale acceptance outputs
and obsolete CI files remain in Git history. Active answerkeys are retained;
scientific benchmarks and large artifacts belong in
[sjaracne_benchmark](https://github.com/adamdingliang/sjaracne_benchmark).
The supported local CWL workflows remain under [`SJARACNe/cwl`](../SJARACNe/cwl/).

## WDL 1.3 engine parity

Run the tiny real CWL/Sprocket comparison directly with:

```bash
python tests/wdl/run_parity.py --sprocket /absolute/path/sprocket \
  --output-dir /absolute/path/new-parity-evidence
```

The complete runner invokes this after unit/CWL/Bats checks. Sprocket 0.32.0
is required; a missing engine fails rather than skipping. Test fixtures use two
worker slots, a 2-GiB backend budget and 1-GiB task requests. The generated
configuration and hash are recorded. These resource overrides do not change
scientific defaults.

## GitHub Actions

[regression.yml](../.github/workflows/regression.yml) runs on pull requests to
`develop` or `master`, pushes to either branch, and manual dispatch. Its Ubuntu
24.04 matrix tests Python **3.11 and 3.14** independently with uv **0.11.16**,
the committed lock and pinned Sprocket. Python 3.14 adds the optional coverage
group; Python 3.11 is uninstrumented.

Each lane prepares one environment under `RUNNER_TEMP`, then uses `--no-sync`.
It checks fresh pip source/wheel/editable installation and uv offline/cache
behavior before the strict regression. These probes use disposable source
copies, not the real checkout. Dependency caches and mutable environments are
not shared between lanes.

The job limit is 20 minutes, with 15 minutes for regression. Superseded runs are
cancelled. Actions are pinned by commit SHA; requirements exports, versions,
logs and source/native hashes are retained as seven-day per-version artifacts.
PR events test the merge revision, recorded in run metadata. No full biological
benchmark batch or cluster job is submitted.

## Execution coverage

Python 3.14 collects separate Python and C++ core reports. The optional coverage
group pins coverage.py **7.17.0** and gcovr **8.8**; neither is a runtime
requirement. Use a dedicated checkout and a new environment outside it:

```bash
export UV_PROJECT_ENVIRONMENT=/absolute/path/sjaracne-coverage-env
uv sync --locked --python 3.14 --group regression --group coverage --extra cwl
uv run --locked --no-sync python tests/run_regression.py --coverage --jobs 4 \
  --output-dir /absolute/path/new-coverage-evidence
```

Coverage mode builds native programs with GCC `-O0`, `--coverage` and atomic
counter updates. Rebuild normally before using that checkout for analyses;
the runner does not restore `-O3` binaries. A checked temporary startup hook in
the external virtual environment captures engine-launched Python, including
CWL tasks. The hook is removed on exit and is never installed in system Python.

Successful runs write `coverage/python.xml`, `coverage/cpp.xml` and
`coverage/summary.json` under the evidence directory. Python measures package
runtime modules. C++ measures `SJARACNe/src` and executable header lines,
excluding the generated coefficient table and separate null-generator command
source. Tests, dependencies, WDL/CWL data and fresh packaging environments are
outside these denominators.

Coverage measures execution, not biological accuracy or production `-O3`
numerical equivalence. GCC branch counts include compiler-generated and
exception paths, so Python and C++ branch percentages are not directly comparable.

CI uploads explicit XML files using Codecov flags `python` and `cpp` after
successful regression. OIDC needs `id-token: write`, not a repository-write
permission or upload secret. External fork PRs preserve reports but skip
uploads. Upload failures and coverage statuses are informational; comments and
stale flag carry-forward are disabled.

For a fork without a Codecov project, enable the
[Codecov GitHub app](https://github.com/apps/codecov) as needed. The dashboard
and README badges use `develop`, with separate flags; a missing upload can show
as unknown. See Codecov's [setup guide](https://docs.codecov.com/docs/quick-start)
and [flag documentation](https://docs.codecov.com/docs/flags).

## Ordinary pip and fresh-source installation

To prepare a separate Python 3.11 or 3.14 environment using the full lock export:

```bash
uv export --locked --group regression --extra cwl --no-emit-project --no-hashes \
  --output-file /tmp/sjaracne-locked-requirements.txt
python -m pip install -r /tmp/sjaracne-locked-requirements.txt
python -m pip install --no-deps --no-build-isolation -e .
```

The export excludes the checkout and includes Python-version markers. Hashes
are omitted so it can also be used as pip constraints; artifact hashes remain
in `uv.lock`. Historical `requirements-*.txt` snapshots do not lock all
transitive dependencies and are not CI's dependency source.

The clean-source check performs installation:

```bash
python tests/check_source_install.py \
  --requirements /tmp/sjaracne-locked-requirements.txt \
  --build-requirements /tmp/sjaracne-locked-requirements.txt \
  --output-dir /tmp/sjaracne-new-source-install-check
```

The uv check needs a new cache and evidence directory. Initial preparation may
download dependencies; subsequent installations use the warmed cache offline:

```bash
python tests/check_uv_install.py \
  --uv "$(command -v uv)" --python "$(command -v python)" \
  --cache-dir /tmp/sjaracne-new-uv-install-cache \
  --output-dir /tmp/sjaracne-new-uv-install-check
```

It requires compiler/Make but not CWL, Sprocket or biological data. It checks the
current host, not arbitrary offline clusters or architectures.

## Cross-Python scientific output comparison

Install this checkout in two environments, then run the tiny compatibility
comparison. Do not run it while a regression rebuilds the same native executable:

```bash
python tests/python_compatibility.py \
  --baseline-python /absolute/path/py311/bin/python \
  --candidate-python /absolute/path/py314/bin/python \
  --sprocket /absolute/path/sprocket \
  --output-dir /absolute/path/new-python-comparison
```

The TF/SIG cases use current-profile sampling/MI settings with explicit
**DPI-first** order; the third case uses tuned, nonempty legacy probability
consensus. Native consensus-first defaults are checked separately by the real
default-workflow tests.

Expression, hubs, native executable, seeds and Sprocket are fixed. Native edges,
MI and selected consensus identities must agree exactly; correlations and
regression statistics allow declared float64 tolerance. Selection is checked
independently. Versions, commands, logs and differences are retained. This is
small-fixture compatibility evidence, not a biological benchmark or universal
bitwise-identity claim.
