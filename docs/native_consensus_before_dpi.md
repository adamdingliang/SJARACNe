# Native consensus before DPI

## Process and interface

The current workflow uses 100 seeded 80% resamples, PCHIP cutoffs, inclusive K20
recurrence, supporting-run arithmetic mean MI and one final DPI pass. Dynamic
DPI uses Href=479, slope=0.05, cap=0.30 and workflow offset=0.025:

```text
epsilon(H) = min(dpi_max, dpi_offset + dpi_slope * max(0, log2(H/dpi_href)))
```

`--raw-mi` disables per-resample DPI and writes every matched source row in
hub-file order, including empty rows, with 17-significant-digit MI.
The consensus invocation reads a manifest of one raw-graph path per line, in
seed order. Relative paths resolve against the manifest directory; entries are
literal paths, not shell expressions:

```bash
sjaracne.exe -i expression.exp -s hubs.txt -l hubs.txt \
  --consensus-inputs raw-inputs.txt --min-recurrence 20 \
  --dynamic-dpi --dpi-href 479 --dpi-slope 0.05 \
  --dpi-offset 0.025 --dpi-max 0.30 -o consensus.adj
```

Inputs must use the complete row layout emitted by `--raw-mi`. Older graphs
with omitted/reordered sources need audited canonicalization before reuse;
rounding lost in historical six-digit MI cannot be recovered.

## Optional raw graph hash verification

`--verify-raw-hashes` is an operational workflow CLI option, not a C++ flag.
Its Boolean input `verify_raw_hashes` defaults to `false`, outside scientific
settings. WDL/Sprocket supports either order; CWL supports consensus-first and
rejects the option for DPI-first, including legacy.

Producers always compute graph SHA-256 and save it in seed receipts. Consensus
checks ordered successful seeds and matching runtime/input/settings fingerprints.
By default it records producer hashes without recomputing graph hashes; enable
the option for a separate sequential pass comparing current bytes with receipts.
The receipt records whether re-verification ran.

Input/runtime/native/configuration, settings and final-output hash checks remain
mandatory. Native parsing still rejects invalid row layouts, IDs, duplicate or
self edges, and invalid MI. Verification changes no scientific calculation for
unchanged inputs.

Skipping re-verification removes one extra read of all graphs, not producer
hashing, scatter writes/transfers or the aggregation read. For 206 GB of graphs,
that avoids another 206 GB read. Hashing is not fused into C++ aggregation;
opt-in verification still uses a separate pass.

## Engineering choices

C++ reads each raw graph completely before opening the next, using one 1 MiB
input buffer. Sequential access avoids alternating small reads across resamples
on mechanical disks and does not require SSD storage.

Dense arrays hold a double-precision MI sum and uint32 support count for each
requested source-target pair. Additions follow manifest order. An ordered edge
counts once per run; zero MI still counts as present. Support K is inclusive,
and absent runs do not enter the mean denominator. After structural validation,
retained edges enter the existing native Matrix directly. Each source's arrays
are released as its edges enter Matrix, before DPI. Empty source rows count in H.

For H matched sources and G expression genes, arrays need
`H * G * (sizeof(double) + sizeof(uint32_t))` bytes, normally **12 bytes per pair**.
This cost is independent of resample count. For example, 9,464 sources and
29,841 genes require about **3.16 GiB** for arrays alone. Expression data, the
retained graph, DPI indices, parser buffers and allocator overhead are additional;
a smaller K can retain a larger graph.

Allocation sizes are checked for overflow; allocation errors report dimensions
and bytes. Aggregation does not spill to disk or reread inputs in source blocks.
The default **8 GiB consensus request** is a scheduling reservation, not a peak
estimate or hard limit. Large panels may need higher `--consensus-memory` and
`--memory` budgets; see [workflow resources](../SJARACNe/wdl/README.md#resources-outputs-and-provenance).

The design removes Python aggregation caches, adjacency staging copies and the
candidate-graph write/read. Scatter still writes/localizes raw graphs, and
full-precision MI uses more disk space than historical six-digit output.
End-to-end time depends on the workload and storage.

WDL/CWL schedule resamples independently. One final task runs native
consensus/DPI, then Python adds gene symbols and correlation/regression statistics
to the pruned graph using the full matrix. Consensus-first outputs preserve MI
tokens; the older order keeps four-decimal formatting. See the
[WDL](../SJARACNe/wdl/README.md) and [CWL](../SJARACNe/cwl/README.md) guides for
stage names and inputs.

## Validation and compatibility

Each graph must contain exactly the canonical source rows and valid unique
non-self targets with finite, nonnegative MI. Missing, extra or reordered rows
are errors. Seed order fixes floating-point addition order; workflow receipt and
hash policy checks are described above.

Native offset defaults to zero; the current consensus-first workflow supplies
0.025, and explicit DPI-first defaults to zero. Dynamic tuning conflicts with
fixed DPI. Legacy retains full-size replacement sampling, affine MI, strict
per-run DPI and probability consensus; see the
[compatibility guide](BACKWARD_COMPATIBILITY.md).

Native output is published only after successful processing; the wrapper also
requires successful workflow provenance and final hashes before publication.
The [benchmark report](https://github.com/adamdingliang/sjaracne_benchmark/blob/b195f8331496c645121fc614c930ce1f58049300/benchmarks/consensus_before_dpi/README.md)
records TP, FP, precision, recall, F1, size-weighted summaries and median-target
retention. Those scientific comparisons are separate from implementation
equivalence and do not establish a universal default optimum.
