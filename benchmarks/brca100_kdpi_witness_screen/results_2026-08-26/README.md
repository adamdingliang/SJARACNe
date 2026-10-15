# Compact evidence package

## Outcome

A fixed minimum witness count does not correct the SIG hub-list-size effect in
this screen. `K_DPI=2`, `3`, and `5` increase the matched full-minus-small
pruning gap. `K_DPI=10` reduces that absolute gap only by making small-panel DPI
nearly inert: it retains 1.66% of small-panel `K_DPI=1` pruning, versus 40.9%
in the full panel, and raises the full/small pruning ratio to 43.35.

The defensible decision is to reject fixed `K_DPI` as a hub-size-bias
mitigation. It could still be studied as a generic conservatism control for a
fixed panel, but that is a different goal.

## Matched design

- BRCA100 SIG networks only.
- Nested hub/annotation panels: `H=1,335`, `5,340`, and `10,680`.
- Ten matched seeds using the frozen pilot's exact 80-of-100 samples.
- Fixed `p_b=5e-4`, AP-MI cutoff `0.1644671599536221`, `Npar=40`, and
  `epsilon_DPI=0`.
- Tested `K_DPI=1,2,3,5,10`; consensus recurrence `K_edge` was not run.
- Primary analysis uses the same 1,335 source rows at every `H`.

Because current DPI marks edges without deleting them or excluding them as
later evidence, an edge is pruned at threshold `k` exactly when its number of
distinct eligible witnesses is at least `k`. One instrumented inference run per
panel/seed therefore gives exact pruning counts for every recorded `K_DPI`.

See [common_source_decision_table.tsv](common_source_decision_table.tsv) for the
primary medians, [REPORT.md](REPORT.md) for the full readable result,
[paired_hub_size_effects.tsv](paired_hub_size_effects.tsv) and
[source_group_summary.tsv](source_group_summary.tsv) for detailed summaries,
[seed_level_metrics.tsv](seed_level_metrics.tsv) for all 300 analysis rows, and
[validation_gates.json](validation_gates.json) for the hard gates. The frozen
[screen design](screen_design.json) and [build manifest](build_manifest.json)
make the compact package independently auditable.

## Validation

All nine required gates passed:

- 30/30 matched runs and 300/300 analysis rows;
- exact `K_DPI=1` sample, DPI-count, and adjacency-data reproduction;
- exact sidecar accounting and frozen candidate provenance;
- 26,700/26,700 common-source pre-DPI edge-count comparisons, with zero mismatches;
- monotone pruning with increasing `K_DPI` and complete paired effects.

The diagnostic passed 32 focused native tests, including exact 3-, 5-, 9-, and
10-witness boundaries. The harness passed nine synthetic tests.

## Provenance and retained raw data

Final external work root:

```text
/home/adam/sjaracne-benchmarks/brca100-kdpi-witness-screen-20260826
```

It retains the frozen build, all adjacency files, witness sidecars, logs, seed
metadata, seed-level metrics, and the exact design. The preliminary three-arm
smoke run is preserved separately at
`/home/adam/sjaracne-benchmarks/brca100-kdpi-witness-screen-20260826-smoke-prehardening`.
The final 10-worker inference span was 501.1 seconds (8 minutes 21 seconds), and
the retained external tree is approximately 120 MiB.

Key identifiers:

- base commit: `32fe12c168ef80291e487dbf4045f430b9c5d90a`
- source-provenance fingerprint: `0af71bd617b5c1d7025a648c6c9b7afcc2b50a281cd8a9028e319d773427ba6b`
- built source-directory SHA256: `d662e2100f8e8f608a00677d4fb2363c471095f195181912627573db9340cf6c`
- candidate binary: `32bdcc55bc834bcea37dd44c94068445279963dde4331e7f34ec4bcbb1d4b621`
- frozen design: `34af423605f6ece9319dc833da20be5ba48b633e1218dd73554110297868323e`

The TSV files in this compact package are LF-normalized text copies; the exact
external originals remain under the work root. `MANIFEST.sha256` fingerprints
the local package.

## Limits

This is a 10-seed screen using one deterministic nested panel trajectory. It
does not emit separate `K_DPI>1` adjacency networks, test consensus topology,
evaluate downstream activity, or establish biological correctness for any
individual edge. Expanding to 100 seeds is not justified for parameter
selection because the tested direction is already consistently unfavorable.
