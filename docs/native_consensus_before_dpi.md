# Native consensus before DPI

## Process and interface

The current workflow resolves to 100 distinct seeded 80% resamples, PCHIP MI
thresholding, inclusive K20 recurrence, arithmetic mean MI over supporting runs,
and a single final DPI pass. Href remains 479, slope 0.05 and cap 0.30. The
adjustable additive offset defaults to 0.025 in this workflow. No additional MI
threshold is applied to the consensus graph.

```text
epsilon(H) = min(dpi_max, dpi_offset + dpi_slope * max(0, log2(H/dpi_href)))
```

The native inference flag `--raw-mi` writes all matched source rows in hub-file
order, including empty rows, and MI with 17 significant digits. It disables DPI
for that resample. Ordinary native inference output retains its historical
serialization unless this flag is selected.

The final native invocation consumes a text manifest containing one raw graph
path per line, in deterministic seed order. Relative paths are resolved against
the manifest directory. Entries are literal filesystem paths, not shell expressions. For example:

```bash
sjaracne.exe -i expression.exp -s hubs.txt -l hubs.txt \
  --consensus-inputs raw-inputs.txt --min-recurrence 20 \
  --dynamic-dpi --dpi-href 479 --dpi-slope 0.05 \
  --dpi-offset 0.025 --dpi-max 0.30 -o consensus.adj
```

This mode requires the canonical complete row layout emitted by `--raw-mi`.
Older saved adjacency files may omit empty sources or use a different row
order. They need explicit, audited canonicalization before reuse; do not
silently interpret a truncated graph as an empty source. Canonicalization of
six-digit historical MI cannot recover the original unrounded values.

## Engineering choices

Each input is read sequentially. Corresponding source rows are accumulated in
gene-indexed counts and double-precision sums. An ordered source-target edge is
counted once per run; support K is inclusive. The mean denominator is that
edge's support count, so absent runs are not zero-filled. MI zero is still a
present edge. Surviving edges are inserted directly into the existing native
Matrix and pruned with the existing DPI implementation. Source rows with zero
survivors remain part of H.

The accumulator uses O(G) storage for G expression genes. Total memory also
includes the retained graph and DPI indices, expression data and input row
buffers; it is not O(G) overall. A loose K can produce a large retained graph.
The process opens the listed run files concurrently. Extremely large resample
counts can exceed the host's file descriptor limit; ordinary B=100 is the
intended operating scale.

This removes Python per-edge aggregation caches, adjacency staging copies and
the candidate graph write/read between consensus and DPI. Distributed scatter
still materializes and transfers raw graphs, and provenance verification still
reads their bytes. Full-precision raw MI increases file size compared with the
historical six-digit output. This change does not reduce MI estimation work,
and no end-to-end speedup is claimed without a measurement.

WDL remains responsible for the independent resample scatter. One post-scatter
task invokes native consensus/DPI, then enhances only the final pruned edge
list using the full expression matrix. Native-order three-column and enhanced
outputs preserve the native MI token precision; historical enhancement defaults
retain their four-decimal format. The scientific calculations remain
outside WDL. CWL consensus-first execution uses the same Python adapter.

## Validation and compatibility

Inputs must have consistent canonical source order and valid expression IDs.
Missing, duplicated, reordered or extra source rows, duplicate targets,
self-pairs and invalid MI values are errors. The workflow additionally checks
all successful seed receipts, input/runtime/settings fingerprints and graph
hashes before aggregation. Seed order determines floating-point addition order.

The native offset defaults to zero to preserve explicit historical native
calls. The current consensus-first workflow supplies 0.025, whereas explicit
DPI-first defaults to zero. Fixed DPI and dynamic tuning options conflict.
Legacy remains full-size replacement sampling, affine MI, strict per-run DPI
and its original probability consensus. New modes publish their native output
only after successful processing; the wrapper publishes final workflow files
only after a successful workflow.

See [the benchmark report](../benchmarks/consensus_before_dpi/README.md) for TP,
FP, precision, recall, F1, size-weighted summaries and median-target retention.
Those saved scientific comparisons are separate from implementation equivalence
checks and do not establish that the selected default wins on every dataset.
