# BRCA100 SIG `K_DPI` witness screen

This is a cheap, matched first test of whether requiring more than one qualifying
DPI intermediary reduces the dependence of DPI pruning on the SIG hub-list size.
It does **not** select a production default and it does not test biological edge
correctness.

## Fixed design

- Network: BRCA100 SIG only.
- Hub/annotation panels (`-s` and `-l` are the same nested panel):
  `H=1,335`, `H=5,340`, and `H=10,680`.
- Seeds: `1..10`, with the frozen pilot's exact fixed 80-of-100 samples.
- AP-MI per-subsample cutoff: `p_b=5e-4`, corresponding to
  `MI=0.1644671599536221` for `m=80`, `Npar=40`.
- DPI tolerance: `epsilon=0`.
- Candidate witness thresholds: `K_DPI=1,2,3,5,10`.
- Direct inference: one candidate run per panel and seed (30 runs), not one run
  per `K_DPI` value.
- Consensus recurrence (`K_edge`) is not run or analyzed.

`K_DPI=1` is the current DPI rule: one qualifying intermediary is sufficient to
remove an edge. The instrumented candidate binary accepts `-W <sidecar.tsv>`,
writes `# key<TAB>value` provenance records, and then writes one row per source
index under this header:

```text
source_index  pre_edges  witnesses_ge_1  witnesses_ge_2  witnesses_ge_3  witnesses_ge_5  witnesses_ge_10
```

Here `witnesses_ge_k` is the number of that source row's pre-DPI edges having at
least `k` qualifying DPI witnesses. It is therefore the exact inferred number of
edges that would be removed from that source row at `K_DPI=k`. The candidate
adjacency itself remains the ordinary `K_DPI=1` output.

## Why two source populations are reported

The all-source statistic uses every source row in each panel. That is useful but
the compared source population changes with `H`.

The primary common-source statistic always uses the 1,335 source indices from
the smallest panel, including when analyzing the 5,340- and 10,680-hub runs.
This holds the evaluated source genes fixed while changing the number of
available reconstructed/annotated hubs. It is the cleaner readout of the
technical hub-list-size effect.

For each seed, source population, panel, and `K_DPI`, the pruning rate is

```text
sum(edges with >= K_DPI witnesses) / sum(pre-DPI edges).
```

The full-minus-small gap is the 10,680-panel pruning rate minus the 1,335-panel
rate. Both its signed and absolute values are retained. The slope divides the
gap by 9.345 thousand added hubs. The full/small ratio is also reported. The
normalized gap is computed within each seed as

```text
gap(K_DPI) / gap(K_DPI=1).
```

A normalized gap of `0` would eliminate the measured K=1 size gap; `1` leaves
it unchanged. This metric can be unstable when a seed's K=1 gap is near zero,
so the absolute paired gap, ratio, and seed-level direction counts remain the
primary evidence. `pruning_retained_vs_k1` is the fraction of K=1 removals that
would still be removed at a higher witness threshold.

## Required candidate interface

The candidate native CLI must accept the optional argument below without
changing ordinary `K_DPI=1` network output:

```text
-W <sidecar.tsv>
```

The sidecar must have the expected provenance schema and exact header shown
above, integer counts, one unique row for every source in `-s`, and monotonically
nonincreasing witness columns.

## Run

Run from WSL. The default paths point to the existing isolated worktree and the
frozen 2026-08-25 pilot:

```bash
python3 benchmarks/brca100_kdpi_witness_screen/run_screen.py \
  --phase prepare

python3 benchmarks/brca100_kdpi_witness_screen/run_screen.py \
  --phase infer --workers 10

python3 benchmarks/brca100_kdpi_witness_screen/analyze_screen.py \
  --require-complete
```

The `prepare` phase fingerprints a candidate worktree snapshot, builds its
binary in the external work root, verifies that its config/null-model hashes
match the frozen baseline, and freezes `screen_design.json`. Both inference and
analysis revalidate the frozen harness hashes. The `infer` phase is resumable: a
seed is reused only after its output, sidecar, log, fingerprint, and all
validation checks pass again.

Useful narrow invocations:

```bash
python3 benchmarks/brca100_kdpi_witness_screen/run_screen.py \
  --phase infer --hub-counts 1335 --seed-start 1 --seed-end 2 --workers 2

python3 -m unittest \
  benchmarks/brca100_kdpi_witness_screen/test_screen.py -v
```

Default external roots:

- Baseline: `/home/adam/sjaracne-benchmarks/brca100-hub-size-dpi-pilot-20260825`
- New screen: `/home/adam/sjaracne-benchmarks/brca100-kdpi-witness-screen-20260826`

## Validation gates

Every candidate run must reproduce the frozen `K_DPI=1` baseline exactly in:

- sampled observation indices;
- global pre/pruned/post DPI counts;
- adjacency data rows and MI strings (header-only differences are tolerated);
- sidecar `sum(pre_edges)` and `sum(witnesses_ge_1)` accounting.

The harness also checks exact source-index coverage, nested common-source
membership, nonnegative integer counts, witness-count monotonicity, and the
dynamic input/panel/annotation/network paths recorded by the sidecar. The
recorded network path is the deterministic temporary path actually executed;
it remains unchanged when the validated files are atomically published. A
report is complete only when all 30 matched candidate runs pass.

## Decision rule and limitation

`K_DPI>1` helps only if it materially reduces the common-source paired
full-minus-small gap across most seeds while retaining a nontrivial share of
K=1 pruning. Making DPI nearly inert is not a solution.

Ten seeds and one deterministic nested panel trajectory are a screening design,
not final calibration. If the signal is promising, repeat the selected narrow
grid with all 100 seeds and additional panel compositions before considering a
code default. This screen says nothing about whether any individual removed or
retained edge is biologically correct.

## Completed result (2026-08-26)

The screen is complete. All 30 matched runs and all nine validation gates
passed. `K_DPI=2`, `3`, and `5` increased the fixed-source full-minus-small
pruning gap. `K_DPI=10` reduced the absolute gap only by retaining 1.66% of
small-panel `K_DPI=1` pruning, versus 40.9% in the full panel. A fixed
`K_DPI` is therefore rejected as a hub-size-bias mitigation in this design.

See the [compact evidence package](results_2026-08-26/README.md). Raw
adjacencies, sidecars, logs, and seed metadata remain under
`/home/adam/sjaracne-benchmarks/brca100-kdpi-witness-screen-20260826`.
