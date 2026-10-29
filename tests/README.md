# Small workflow regression tests

The goal is functional correctness of the complete workflow after code changes,
not biological accuracy or performance benchmarking. Routine tests do **not**
launch BRCA100/DREAM resampling batches or require an LSF cluster.

## Complete check (Linux or WSL)

Use Python **3.11 or 3.14** with the matching pinned dependency set, a C++
compiler/Make, Node.js, `cwltool`, Bats and Sprocket 0.32.0. Install this checkout
in the active environment so Python-module CWL commands cannot silently import
another revision:

```bash
python -m pip install -r tests/requirements-build.txt
# Python 3.14; use tests/requirements-ci.txt instead for Python 3.11.
python -m pip install -r tests/requirements-ci-py314.txt
python -m pip install --no-deps --no-build-isolation -e .
python tests/run_regression.py --jobs 4
```

On Ubuntu, the system prerequisites can be installed with
`sudo apt-get install build-essential nodejs bats`. The CI manifests include
the optional CWL dependencies. Install the separately distributed
Sprocket executable as described in the [WDL guide](../SJARACNe/wdl/README.md).
The runner itself installs nothing and downloads neither data nor dependencies.

The command can be invoked from another working directory using its absolute
path. It creates a fresh evidence directory under the system temporary
directory and prints the location. To choose a **new** directory explicitly:

```bash
python tests/run_regression.py --jobs 4 --output-dir /tmp/sjaracne-regression-check
```

It performs the following checks:

1. Fail early for missing dependencies or an installed Python/CWL/configuration
   runtime whose bytes differ from this checkout. CWL jobs normally do not
   inherit the parent's `PYTHONPATH`; setting it alone is insufficient.
2. Force-rebuild both `sjaracne.exe` and `apmi_null_generator.exe`. The force
   rebuild matters because the native Makefile does not track every header.
3. Run every `test_*.py` test with the freshly built native executables and
   required-native mode. Empty/duplicate discovery, failures, errors, skips and
   expected failures do not count as a complete successful regression.
4. Validate all production CWL documents and run the isolated Bats CWL utility
   tests, requiring a nonempty, all-passing plan with no skips.
5. Record commands, versions, source/fixture hashes, native binary hashes,
   counts and logs. Reject source/fixture changes during the check.
6. Run real current TF/SIG defaults and the legacy profile through the CLI,
   plus the paired WDL/CWL engine matrix. Missing Sprocket is a failure, not a
   skipped test or a reason to switch engines silently.

## Coverage and fixtures

| Coverage | What is checked |
| --- | --- |
| Scientific core compatibility | Generated six-gene/100-observation expression, three hubs, six resamples: input QC, 80% sampling, PCHIP, real DPI removals, K=6 and enhanced consensus. Strict-DPI baseline settings are explicit in the historical CWL/CLI cases; answerkeys are not regenerated to match new defaults. |
| Current default workflow | Actual Sprocket execution of the CLI for TF and SIG, without explicit MI/DPI settings; verifies role-specific pb, PCHIP, dynamic DPI, without-replacement sampling, K=6, published outputs and provenance. |
| Pre-sampling compatibility | Explicit legacy profile with replacement, affine conversion, fixed DPI and legacy consensus through WDL and CWL; the native low-level defaults remain unchanged. |
| Selection compatibility | Tiny affine/PCHIP × K=6/explicit legacy `-pc` matrix; explicit probability behavior remains supported and must not be deleted as obsolete. |
| Native correctness | Sampling/seed behavior, rank caches, AP-MI workspace/null generator, PCHIP bounds/provenance/warnings, expression/hub validation, adjacency replay and DPI diagnostics/file-alias protection. |
| Experimental dynamic DPI | Native formula/anchor/doubling/cap, exact fixed-epsilon edge equivalence, unique constructed source counting (including empty rows), invalid inputs and fixed/dynamic conflicts, opt-out compatibility, provenance and witness diagnostics. Wrapper/CWL tests cover forwarding, real dynamic-versus-fixed inference, and two scattered resamples through consensus. All use tiny synthetic data. |
| Consensus and enrichment | Recurrence boundaries, directed duplicate handling, legacy probability behavior, missing/empty inputs, output schemas, symbols and finite/bounded enhanced statistics. |
| Utilities | Byte-exact line endings; checked CWL array types/order, empty arrays, output file contents and copies. Bats fixtures/results live in per-test temporary directories, not the checkout. |
| Standalone WDL task plumbing | Parameters/defaults/conflicts, literal argv, safe line-ending copies, runtime/input/settings hashes, seed provenance and exact consensus staging. These ordinary tests require no WDL engine. |
| Packaging and Python compatibility | Build-backend/package-data contracts and strict cross-version output-comparator tests. CI separately builds a source archive, compiles and installs its platform wheel in a fresh environment, and checks installed assets and native execution. |

The small T-cell `cnn_5.txt` golden network and the active legacy consensus
answerkeys are intentionally retained as compatibility fixtures. The BRCA100
input files remain bundled demonstration data used by the main README; keeping
them does not make full resampling part of this suite. Rank-cache fixtures are
native correctness inputs, not benchmark datasets.

For a focused developer check, use, for example:

```bash
python -m unittest discover -s tests -p 'test_create_consensus_network.py' -v
```

Direct unit-test discovery may skip unavailable native programs for convenience.
That is **not** equivalent to the strict complete runner. Set
`SJARACNE_REQUIRE_NATIVE_TESTS=1` to make native discovery fail rather than skip.
An invalid explicitly configured `SJARACNE_TEST_EXE` or
`SJARACNE_APMI_NULL_EXE` always fails; it never silently selects another binary.
The complete runner pins both variables and PATH to this checkout.

The dynamic-DPI modules are discovered by the same `test_*.py` pattern and
therefore run in GitHub Actions without a separate benchmark job. For a focused
check after building the native binary and installing this checkout:

```bash
python -m unittest discover -s tests -p 'test_sjaracne_dynamic_dpi.py' -v
python -m unittest discover -s tests -p 'test_dynamic_dpi_workflow.py' -v
```

## Removed or separated during maintenance

Unused 2019 cluster logs, absolute-path job files, unreferenced 100-run
acceptance snapshots, the disabled enhanced-output test/golden, and obsolete
Travis configuration were removed. They remain recoverable from Git history;
active answerkeys were not regenerated to make tests pass.

Non-asserting local CWL launchers and their jobs remain under
[`examples/manual_cwl`](../examples/manual_cwl/README.md), outside automated test
discovery. Obsolete CWLEXEC/LSF launchers and configuration templates were
removed; their history remains in Git. No cluster jobs are submitted by tests.

Scientific benchmarks and large run artifacts remain in
[`sjaracne_benchmark`](https://github.com/adamdingliang/sjaracne_benchmark).

## WDL 1.3 engine parity

The [WDL guide](../SJARACNe/wdl/README.md) documents the tiny real CWL/Sprocket
comparison, which can also be run directly:

```bash
python tests/wdl/run_parity.py --sprocket /absolute/path/sprocket \
  --output-dir /absolute/path/new-parity-evidence
```

It requires Sprocket 0.32.0 and fails for missing dependencies, rather than
skipping. Its standalone runner is outside `test_*.py` discovery, but the full
regression invokes it explicitly after the ordinary tests and CWL utilities.
GitHub Actions installs a checksum-verified Sprocket and therefore executes
both real default CLI runs and the paired engine matrix. No production dataset
benchmarks are launched by either check.

## GitHub Actions

`.github/workflows/regression.yml` runs this same strict small regression on
pull requests targeting `develop` or `master`, and on pushes to either branch.
The workflow is introduced on `develop`; `master` gains it and the refreshed
runner together through the next release merge, not by silently promoting
scientific changes now. This does not change branch-protection requirements.

The hosted matrix uses standard Ubuntu 24.04 with **Python 3.11 and 3.14**.
The 3.11 lane retains the numerical/CWL pins in `requirements-ci.txt`; the 3.14
lane uses `requirements-ci-py314.txt`. Both use `requirements-build.txt` for
build tools, and both must pass independently (`fail-fast: false`, no allowed
failures). These are tested environment pins, not new runtime dependency
minimums. Sprocket's external binary is pinned to 0.32.0 and its release-archive
SHA256. Top-level Python pins do not lock every transitive dependency: the resolved
`pip freeze`, tool versions, test logs, counts and source/binary hashes are saved
as a version-specific seven-day Actions artifact. First-party actions are pinned by commit SHA,
permissions are read-only, and no repository secret is required.

The job is bounded to 20 minutes (15 minutes for regression), cancels superseded
runs, and does not invoke manual examples or large biological benchmark batches.
Manual dispatch is also available once the workflow is present on the default
branch. GitHub's regular PR event tests its merge revision; the run metadata
records the tested SHA and event.

## Cross-Python scientific output comparison

After installing this checkout in two separate environments, run the additional
tiny compatibility comparison. Do not run it while a regression rebuilds the
same native executable:

```bash
python tests/python_compatibility.py \
  --baseline-python /absolute/path/py311/bin/python \
  --candidate-python /absolute/path/py314/bin/python \
  --sprocket /absolute/path/sprocket \
  --output-dir /absolute/path/new-python-comparison
```

It holds expression values, hubs, native executable, seeds and Sprocket fixed.
Current TF/SIG settings and an explicitly tuned, nonempty legacy-consensus case
are compared. Native edges/MI and selected consensus identities must agree
exactly; correlations and regression statistics use a small declared float64
tolerance. Selection is independently checked rather than treating agreement
between environments as sufficient evidence. Versions, commands, logs and
differences are retained. This is separate from a biological benchmark or a
claim of bitwise identity for arbitrary datasets.
