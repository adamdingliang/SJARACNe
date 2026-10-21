# Small workflow regression tests

The goal is functional correctness of the complete workflow after code changes,
not biological accuracy or performance benchmarking. Routine tests do **not**
launch BRCA100/DREAM resampling batches or require an LSF cluster.

## Complete check (Linux or WSL)

Use Python 3.8+ with compatible package dependencies, a C++ compiler/Make,
Node.js, `cwltool`, and Bats. Install this checkout in the active environment so
Python-module CWL commands cannot silently import another revision:

```bash
python -m pip install -e .
python tests/run_regression.py --jobs 4
```

On Ubuntu, the system prerequisites can be installed with
`sudo apt-get install build-essential nodejs bats`. Python dependencies come
from the package installation. The runner itself installs nothing and downloads
neither data nor dependencies.

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

## Coverage and fixtures

| Coverage | What is checked |
| --- | --- |
| Current end-to-end workflow | Generated six-gene/100-observation expression, three hubs, six resamples: real input QC, default 80% sampling, PCHIP, actual DPI removals, default K=6, and nonempty enhanced consensus. Real serial CWL, parallel CWL and CLI outputs must agree. Expected recurrence/MI is derived from the actual native networks, not a frozen scientific network. |
| Selection compatibility | Tiny affine/PCHIP × K=6/explicit legacy `-pc` matrix; explicit probability behavior remains supported and must not be deleted as obsolete. |
| Native correctness | Sampling/seed behavior, rank caches, AP-MI workspace/null generator, PCHIP bounds/provenance/warnings, expression/hub validation, adjacency replay and DPI diagnostics/file-alias protection. |
| Consensus and enrichment | Recurrence boundaries, directed duplicate handling, legacy probability behavior, missing/empty inputs, output schemas, symbols and finite/bounded enhanced statistics. |
| Utilities | Byte-exact line endings; checked CWL array types/order, empty arrays, output file contents and copies. Bats fixtures/results live in per-test temporary directories, not the checkout. |

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

## Removed or separated during maintenance

Unused 2019 cluster logs, absolute-path job files, unreferenced 100-run
acceptance snapshots, the disabled enhanced-output test/golden, and obsolete
Travis configuration were removed. They remain recoverable from Git history;
active answerkeys were not regenerated to make tests pass.

Non-asserting CWL/LSF launchers and their jobs are retained under
[`examples/manual_cwl`](../examples/manual_cwl/README.md) and `examples/lsf/`,
outside automated test discovery. Some launchers can submit large batches;
they require explicit manual invocation and are not local CI validation of LSF.

Scientific benchmarks and large run artifacts remain in
[`sjaracne_benchmark`](https://github.com/adamdingliang/sjaracne_benchmark).
No GitHub Actions workflow is enabled by this cleanup.
