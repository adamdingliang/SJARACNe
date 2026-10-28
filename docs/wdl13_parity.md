# Standalone WDL 1.3: local scientific parity

Validation date: 2026-10-28. Implementation branch:
`feat/wdl13-scientific-parity`, based on `develop` commit
`4ef87a01698394dd020aa6e2c0b66c6cab230ec2`.

The local parity matrix passed: **259 assertions, zero failed assertions**.
The complete small regression also passed: **259 Python tests, ten CWL document
validations and four Bats tests**, with no skipped or expected-failure tests.

## Design

The [workflow](../SJARACNe/wdl/sjaracne.wdl) uses WDL **1.3**, tested with
Sprocket **0.32.0** on the WSL/Linux local host backend. The local example has
two task slots, no container isolation, no retries and no call caching.
The existing Python CLI and production CWL documents are unchanged.

Both engines use the same freshly built native executable, configuration,
Python/numerical environment and synthetic inputs. The fixture is the existing
`write_multihub_dpi_inputs()` generator: six genes, 100 observations and three
annotated sources with genuine DPI pruning opportunities. Seeds are paired
exactly, starting at one, with 80% sampling without replacement and Npar=40.
The same full source list supplies both native `-s` and `-l`.

The paired matrix covers:

| Configuration | Resamples per engine | Consensus |
| --- | ---: | --- |
| PCHIP, explicit pb=0.001, strict DPI | 6 | Omitted selection resolves to K=6 |
| Same PCHIP, dynamic DPI with Href=2, slope=0.1, cap=0.3 | 6 | Default K=6 |
| Legacy affine, omitted pb resolves to 1e-7, strict DPI | 6 | Default K=6 |
| PCHIP, pb=0.001, strict DPI | 3 | Explicit legacy pc=0.01 |

All four configurations matched exactly on the compared scientific outputs.
The three K=6 configurations each retained six consensus edges in both engines.
The three-resample legacy pc=0.01 case retained zero edges in both engines;
this expected selection result was checked independently, not counted as a
nonempty-network demonstration.

The dynamic fixture deliberately has H=3 > Href=2, so epsilon is positive
(approximately 0.05849625), rather than merely exercising a zero-epsilon path.
It is a functional parity configuration, not a recommended biological anchor.

Every seed's adjacency body is compared byte-for-byte. All native scientific
headers must agree; only localized `Input file` and `Subnetwork file` values
may differ. Three-column consensus, enhanced network and bootstrap metadata
must agree exactly. Parameter metadata additionally permits the localized
`Output network` path to differ. An independent calculation counts distinct
ordered edges per run and checks recurrence selection and mean printed MI.
Enhanced correlations/slopes use the full expression matrix in both engines.

Additional WDL checks compare CRLF inputs to the LF baseline without modifying
original bytes, and require missing PCHIP probability or zero hub overlap to
fail before native inference. Existing CWL's CRLF staging defect is not used as
a correctness oracle. Native binary, runtime source, configuration, normalized
input and prepared-settings hashes protect the task boundaries.

## Evidence and reproduction

| Check | Result |
| --- | --- |
| WDL 1.3 analysis and standard formatting | Passed with Sprocket 0.32.0 |
| Four real paired CWL/WDL configurations | Exact per-seed and consensus agreement |
| Positive dynamic epsilon | Approximately 0.05849625, with actual DPI removals in every tested run |
| CRLF WDL versus LF WDL | Same normalized inputs and scientific outputs; originals unchanged |
| Missing PCHIP probability and zero hub overlap | Both rejected before native inference |
| Independent recurrence/mean-MI checks and metadata checks | Passed |
| Entire small regression | 259 Python tests, ten CWL validations, four Bats tests passed |

Environment: Python 3.9.25, NumPy 1.20.1, pandas 1.2.3, SciPy 1.7.3,
cwltool 3.1.20211107152837, Node.js 18.19.1, Sprocket 0.32.0
(`70e29d31d`). The tested native executable SHA256 was
`dbdbe074f1e6263cddf296dd0f9ebbacc437edbb47e0904ffc79f34a119a4fde`.

Full local evidence is not versioned:

- Engine comparison: `D:\GitHub\local_runs\wdl13_parity_20261028_r2\summary.json`
  and its `commands/` and `cases/` directories.
- Complete regression: `D:\GitHub\local_runs\wdl13_regression_20261028_r4\manifest.json`,
  `python_tests.json` and logs.

The regression finished before a version-string parser correction to the separate
engine harness. No scientific/runtime files changed afterward; the corrected
harness was then executed in full and froze its own tested source hashes.

The [WDL guide](../SJARACNe/wdl/README.md) gives inputs, outputs and exact local
commands. The engine comparison is separate from ordinary regression discovery:

```bash
python tests/wdl/run_parity.py --sprocket /absolute/path/sprocket \
  --output-dir /absolute/path/new-parity-evidence
python tests/run_regression.py --jobs 2 \
  --output-dir /absolute/path/new-regression-evidence
```

Both runners require dependencies already installed and preserve logs and hashes.
Neither runs BRCA100, LUAD, DREAM5 or any full benchmark batch. Ordinary regression
includes 25 new task-plumbing tests but does not install or invoke Sprocket.

## Conditional conclusion

With the same native binary, inputs, configuration and numerical environment,
the standalone WDL 1.3 workflow preserves the tested CWL scientific operations
and defaults on this tiny local fixture. This supports continuing the opt-in
workflow migration, not replacing CWL universally yet. It does not independently validate
AP-MI, DPI, PCHIP calibration or biological edge accuracy. It does not establish
production-scale memory behavior, cluster/cloud portability, container execution,
caching or failure recovery. CWL remains the supported fallback while those later
migration steps are evaluated.
