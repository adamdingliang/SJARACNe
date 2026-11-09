# Native consensus-first workflow validation

This document retains dated validation records and their original receipts.
The [sequential-input revision](#sequential-input-revision) describes the current
C++ read design and its focused local checks. The earlier
[workflow refresh](#workflow-validation-refresh) and the initial snapshot below
record the implementations tested at their stated dates.

## Sequential-input revision

Revision date: **2026-11-09**. Native consensus now reads each complete raw-MI
file in manifest order before opening the next. It keeps one raw input and one
1 MiB read buffer open, replacing the initial design's source-by-source reads
across all resample files. No SSD is required by this access pattern. This is
an I/O design change; no cold-HDD or end-to-end speedup result is claimed here.

Support and MI sums are accumulated in RAM using dense per-source arrays:
`H * G * (sizeof(double) + sizeof(uint32_t))` bytes, normally 12 bytes per
source-target pair. The current large SIG panels need approximately **3.16 GiB**
(LUAD) and **3.14 GiB** (BLCA) for these arrays alone. Expression data, the retained
graph, parser buffers and allocator overhead require additional memory. Each
source's arrays are released while the K-supported graph is built. No aggregation
cache, adjacency copy or repeated-input-pass fallback is introduced. The default
local consensus scheduling request remains **8 GiB**; larger panels may require
higher `--consensus-memory` and overall `--memory` budgets.

Focused local checks for this revision passed on Linux/WSL:

- A fresh native build and all **22 native consensus tests** passed without
  skips. New checks cover 80 resamples under a 32-file descriptor limit, complete
  validation of one file before the next, nonmonotonic source subsets, finite-MI
  sum overflow and an explicit accumulator-allocation failure.
- The descriptor-limit and complete-file-order tests both failed on the previous
  interleaved-reader binary, confirming that they detect the access-pattern
  regression rather than only checking final numerical results.
- Three long-row comparisons produced byte-identical complete adjacency outputs
  against the original native implementation: dynamic DPI offsets 0 and 0.025,
  and fixed epsilon 1. The fixture has 50,002 genes, rows larger than 1 MiB,
  40,000 retained targets, 10,000 emitted-zero edges and an empty source.

The pull-request CI runs the complete small-workflow regression under Python
3.11 and 3.14, including real WDL/CWL execution. The dated test counts, hashes
and JSON receipts below remain evidence for earlier revisions. No production
bulk resampling or cold-HDD speed benchmark was run for this revision.

## Initial validation snapshot

Validation date: 2026-11-08. Historical feature branch: `feat/native-consensus-before-dpi`,
based on fork develop `78096288ab35f629452ca3b9339c65175200dec0`.

## Implemented behavior

The current profile runs 100 seeded 80% resamples, keeps above-cutoff raw MI,
selects edges supported in at least 20 runs, computes supporting-run mean MI,
and applies the existing native DPI implementation once to that graph. Dynamic
DPI uses Href=479, slope=0.05, cap=0.30 and default additive offset=0.025.
`--dpi-offset` changes the offset; `--consensus-order dpi-first` selects the
previous current order with an omitted offset resolving to zero. The explicit
legacy profile retains its original scientific settings.

The new native modes retain 17-digit MI and every matched source row, including
zero-target rows. The final enhanced native-order output preserves the exact MI
tokens. Historical enhancement formatting remains four decimal places.

## External dataset benchmarks

Dataset-level method comparisons, result tables, saved-data equivalence evidence
and the replay driver are maintained in the dedicated benchmark repository:
[study and results](https://github.com/adamdingliang/sjaracne_benchmark/blob/b195f8331496c645121fc614c930ce1f58049300/benchmarks/consensus_before_dpi/README.md) and
[equivalence receipt](https://github.com/adamdingliang/sjaracne_benchmark/blob/b195f8331496c645121fc614c930ce1f58049300/benchmarks/consensus_before_dpi/native_equivalence.json).
This document retains software correctness, build and workflow regression checks.

## Functional checks

The strict regression's forced native rebuild produced SHA-256
`e9fc297e3025b8aeb6d236a57290cb05edfbdd3f23e2f429b4f3ecf2d4cd9475`.

- 17 focused new native tests passed: K boundaries, 300-run support counts,
  supporting-run means, deterministic summation order, empty rows, invalid
  inputs, duplicates/hardlinks, output aliases and failed writes.
- Real new-default TF and SIG workflows, explicit legacy workflows, and
  WDL/CWL parity at offsets 0 and 0.025 passed. The precision/retention rerun
  verified exact MI tokens in both final formats and retained CWL resample
  graphs, receipts and logs after engine cleanup.
- Historical compatibility tests explicitly select DPI-first. Their golden
  outputs and independent selection oracles were not changed.
- A file-size-limited native failure test isolates its gcov output in a private
  directory. A disposable instrumented check confirmed that shared coverage
  counters remain unchanged.

The complete strict regression passed under Python 3.11.15 on Linux/WSL:

| Check | Result |
| --- | --- |
| Fresh native build | Passed; binary hash recorded above |
| Python/native/CLI and real workflow tests | 388 passed, zero failures/errors/skips |
| Production CWL documents | 14 validated |
| Bats utility tests | 4 passed |
| Independent WDL/CWL parity matrix | 563 assertions passed, zero failed |

The focused tests above are included in the complete test count; they are not
additional tests to add to 388. The independent matrix retains the old-order
oracles; the ordinary test suite separately exercises the new consensus-first
default and both offsets through WDL and CWL. Original acceptance goldens were
preserved by making their intended DPI-first order explicit.

[Validation receipt](native_consensus_validation.json) records source hashes,
compiler outputs, dependency versions and exact local evidence-file hashes.
It is a dated regression snapshot from before the CWL naming cleanup; its older
CWL paths and hashes intentionally describe the files tested at that time,
including wrappers subsequently removed.
The full retained run is `native_consensus_regression_20261108_py311_r2` in the
local benchmark-results directory.

## Engineering review and limits

Two independent design reviews, followed by implementation cross-review,
checked support semantics, direction, empty-source H, precision, deterministic
order, memory scaling and publication evidence. At this **2026-11-08 snapshot**,
the implementation used O(G) per-row accumulators plus the retained graph,
expression data and input buffers; total memory was not O(G). That design kept
all resample streams open and is superseded by the sequential-input revision
above. The snapshot removed Python edge aggregation, adjacency staging copies
and the mean-candidate write/read before DPI.
Distributed raw scatter files and their transfers remain necessary. No
end-to-end speedup measurement is claimed.

Completion checks bind the runtime, prepared settings, normalized and original
input hashes, seed set and output hashes. Each CWL run receives a new owned
scratch child. Failed publication is recorded as failure, with evidence kept.

These local checks validate Linux/WSL execution. They do not establish native
Windows, production multi-node HPC or hosted CI execution. The complete local
regression uses Python 3.11; no new Python 3.14 run is claimed here. Dataset-level
scientific interpretation and method comparisons are documented externally.

## Workflow validation refresh

Refresh date: **2026-11-08 (America/New_York)**. Tested runtime revision:
[`8acd9c6`](https://github.com/adamdingliang/SJARACNe/commit/8acd9c6a1a0399eb369705df62711f72b39fb23b).
The complete strict regression rebuilt both native executables and ran against
an unchanged source tree on Linux/WSL with Python **3.11.15**,
Sprocket **0.32.0** and cwltool **3.1.20240508115724**.
NumPy **1.26.4**, pandas **2.2.3** and
SciPy **1.13.1** matched the Python 3.11 environment pins.

| Check | Result |
| --- | --- |
| Fresh native build | Passed |
| Python/native/CLI/workflow tests | 388 passed; zero failures, errors or skips |
| Remaining CWL documents | 12 validated |
| Bats utility tests | 4 passed |
| Independent WDL/CWL parity matrix | 563 assertions passed; zero failures |
| Additional SIG consensus-first WDL/CWL comparison | Passed; identical final network/consensus files and exact MI tokens |

The ordinary tests exercise real new-order TF and SIG workflows, WDL/CWL
consensus-first equality at offsets 0 and 0.025, legacy-profile equality and
invalid-input handling. The independent matrix checks the explicit DPI-first
and legacy paths against its support/mean-MI oracles, plus CRLF normalization.
The extra SIG comparison used three synthetic resamples, K=2, pb=0.0005 and
offset=0.025; it also verified retained raw graphs, provenance and native logs.
These are small software checks, separate from dataset benchmarking.

Fresh native SHA-256:
`8c6ac630d698ccd6e74a29c69da76c2772fd83a5cc444d29219be2551305a894`.

The tested runtime was frozen throughout the regression. Subsequent changes
update only documentation and its badge-link assertions; the related checks
were rerun after applying the README corrections. The README's Bash examples
were syntax-checked and its `sjaracne` options parsed without executing install
commands. This refresh does not claim a new Python 3.14 or hosted-CI result.
The badges remain explicitly scoped to develop.

Retained local evidence: `workflow_readme_validation_20261108_8acd9c6/`, including
`regression/manifest.json`, `regression/python_tests.json`,
`regression/wdl_parity/summary.json` and `sig_engine_parity/receipt.json`.
The full suite can be reproduced in the prepared Python 3.11 environment with:

```bash
python tests/run_regression.py --jobs 2 --output-dir /absolute/path/new-evidence
```
