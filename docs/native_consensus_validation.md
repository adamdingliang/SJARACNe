# Native consensus-first workflow validation

Validation date: 2026-11-08. Feature branch: `feat/native-consensus-before-dpi`,
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
The full retained run is `native_consensus_regression_20261108_py311_r2` in the
local benchmark-results directory.

## Engineering review and limits

Two independent design reviews, followed by implementation cross-review,
checked support semantics, direction, empty-source H, precision, deterministic
order, memory scaling and publication evidence. The final implementation uses
O(G) per-row accumulators plus the retained graph, expression data and input
buffers; total memory is not O(G). It removes Python edge aggregation,
adjacency staging copies and the mean-candidate write/read before DPI.
Distributed raw scatter files and their transfers remain necessary. No
end-to-end speedup measurement is claimed.

Completion checks bind the runtime, prepared settings, normalized and original
input hashes, seed set and output hashes. Each CWL run receives a new owned
scratch child. Failed publication is recorded as failure, with evidence kept.

These local checks validate Linux/WSL execution. They do not establish native
Windows, production multi-node HPC or hosted CI execution. The complete local
regression uses Python 3.11; no new Python 3.14 run is claimed here. Dataset-level
scientific interpretation and method comparisons are documented externally.
