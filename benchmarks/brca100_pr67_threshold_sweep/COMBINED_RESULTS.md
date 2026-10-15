# BRCA100 PR67 per-subsample threshold sweep: combined results

## Bottom line

This matched sweep tested whether PR67's per-subsample AP-MI p-value explains
the severe sparsity observed at the current `p=1e-7` default, and then mapped
the tradeoff among network topology, regulon density, and recurrence as that
single parameter was relaxed. It was not designed to optimize the consensus
threshold, DPI, or any other network setting.

The answer to the first question is clear for BRCA100: `p=1e-7` is too
stringent if the engineering goal is broad driver coverage, connectivity, and
nontrivial target sets. The original nine-point sweep found that `p=3e-4` was
the first tested point to pass the declared topology floors for both TF and
SIG, but its median target count was only 11 in both classes. The four-point
extension increased the zero-filled median target count to 21 / 21 at `p=5e-4`
and 51 / 49 at `p=1e-3` for TF / SIG.

The stability result is less favorable. At `p=1e-3`, restricting to edges
seen in at least 20 of 100 matched seed-level subnetworks reduces those medians
to 9 / 9. At 50% and 80% support, the zero-filled median is zero for both
classes at every extended point. Loosening the threshold therefore adds many
edges, but most of the added density is low recurrence.

That does not prove that the added edges are false. Recurrence is a subsample
stability measure, not biological truth, an edge probability, or FDR. It also
does not establish that the added edges improve activity estimation. No
biological-reference, independent-cohort, or downstream activity-robustness
criterion was evaluated. Consequently, this sweep supports rejecting
`p=1e-7` as a useful BRCA100 topology/density operating point, but it does
**not** support a new global default.

## Goal of the sweep

The work had two linked phases:

1. The original nine-point sweep asked whether changing only the
   per-subsample PR67 p-value could recover driver coverage and graph
   connectivity, and where that topology transition occurred.
2. The four-point extension asked how much further relaxation was needed to
   increase target-set sizes that may be more useful for downstream activity
   estimation, while reporting how much of that density was recurrent across
   subsamples.

The extension was descriptive, not a retrospective optimization. No target
density, recurrence, or downstream-utility rule was declared that would
select one of `4e-4`, `5e-4`, `7.5e-4`, or `1e-3`. The topology selector still
returns the smallest point that passes its previously declared engineering
floors, `p=3e-4`; by construction, adding looser points cannot change that
answer.

## Matched design and provenance

Only the seed-level PR67 p-value changed. Every arm used:

- PR67 commit `7633ebb4a0d966dbda15a4e32d0efa492fb71aeb`;
- the estimator-matched `m=80`, `Npar=40` null model;
- BRCA100 with 2,608 candidate TFs, 10,680 candidate SIGs, and 28,278
  expression nodes;
- fixed 80% sampling without replacement and seeds 1 through 100;
- adaptive-partitioning MI with `Npar=40`;
- DPI pruning with epsilon 0; and
- consensus-network `p=1e-5`.

The threshold grid contained 13 points. Each point was run independently for
TF and SIG at all 100 seeds, for 2,600 seed-level inference jobs. Seed identity
and every non-threshold setting were matched across arms. All threshold-point
networks were inferred directly; none was produced by filtering a looser
post-DPI adjacency.

NetBID2 QC was run on the final consensus networks after seed-level MI
filtering and DPI pruning. Thus, the reported topology and target-size metrics
are post-DPI, post-consensus quantities. The NetBID2 environment was pinned to
R 4.4.3, NetBID2 2.2.0 at commit
`5defa454d600b94f5dd6d1f9f4428f99759a6821`, and igraph 2.3.3.

The four new probabilities added 800 direct inference jobs without changing
the original 1,800-job result set. The design migration is append-only: the
original design SHA-256
`9fd4db540575d03d0fb50aeb7b4860f3c3afd24a5190186d04ae2523eb4ceefb`
is archived with the migration record, and the active 13-point design SHA-256
is
`2eca0f9ef388fdf48632f8292c226bac59e0e325d9641574b314689651c6a8e5`.
The compact artifacts retain job ledgers, point and arm manifests, build and
environment provenance, omitted-artifact hashes, and package-wide SHA-256
checksums.

The two result packages are preserved separately:

- [Original nine-point package](results_2026-08-19/) and its
  [results report](results_2026-08-19/RESULTS.md)
- [Extended 13-point package](results_extended_2026-08-19/) and its
  [results report](results_extended_2026-08-19/RESULTS.md)

All 400 anchor comparisons passed. The `p=1e-7` data sections reproduce the
prior PR67 run, and the `p=0.000352804562601613` data sections reproduce the
PR66 cutoff-match run for both TF and SIG across all 100 seeds.

## Calibration provenance

Gate 2 is the independent held-out calibration check of the fitted AP-MI
independence-null tail. It is not an independent BRCA cohort or biological
validation.

- `p=2e-5`, `5e-5`, `1e-4`, `2e-4`, `5e-4`, and `1e-3` are exact Gate-2
  grid probabilities.
- `p=3e-4`, the PR66 cutoff match, `4e-4`, and `7.5e-4` are interpolations
  inside the accepted held-out range of `2e-5` through `2e-3`.
- `p=1e-7`, `1e-6`, and `1e-5` extrapolate below the directly testable
  held-out range.

Calibration class records how the AP-MI cutoff was obtained. It does not rank
the biological quality or downstream usefulness of a resulting network.

## Full 13-point density and topology trajectory

Values are shown as TF / SIG. Active is the percentage of candidate drivers
with at least one consensus edge. Incident is the percentage of all 28,278
expression nodes incident to an edge. LCC is the percentage of incident nodes
in the largest weak component. Median targets are zero-filled over all
candidate drivers. Adjusted R-squared is NetBID2's scale-free-fit diagnostic;
it is reported as QC and was not used to select a threshold.

| Per-seed p | MI cutoff | Calibration | Consensus edges | Median targets | Active | Incident | LCC | NetBID2 adjusted R-squared |
|---:|---:|---|---:|---:|---:|---:|---:|---:|
| `1e-7` | 0.322465 | GPD extrapolation | 4,672 / 13,183 | 0 / 0 | 10.31% / 13.93% | 6.08% / 11.56% | 77.25% / 82.32% | 0.737 / 0.932 |
| `1e-6` | 0.288001 | GPD extrapolation | 7,555 / 22,287 | 0 / 0 | 15.91% / 21.22% | 8.92% / 17.72% | 80.89% / 85.32% | 0.756 / 0.940 |
| `1e-5` | 0.247951 | GPD extrapolation | 15,258 / 44,165 | 0 / 0 | 31.37% / 38.36% | 16.45% / 32.94% | 89.72% / 88.53% | 0.778 / 0.932 |
| `2e-5` | 0.234669 | exact Gate 2; boundary | 19,910 / 62,077 | 0 / 1 | 41.60% / 50.89% | 21.25% / 43.71% | 88.79% / 88.49% | 0.813 / 0.935 |
| `5e-5` | 0.216162 | exact Gate 2 | 30,216 / 90,860 | 1 / 2 | 62.54% / 68.69% | 31.59% / 59.75% | 90.02% / 91.11% | 0.815 / 0.945 |
| `1e-4` | 0.201409 | exact Gate 2 | 39,731 / 125,572 | 2 / 3 | 79.56% / 84.67% | 41.54% / 75.02% | 91.37% / 95.21% | 0.819 / 0.941 |
| `2e-4` | 0.185973 | exact Gate 2 | 67,194 / 196,306 | 7 / 7 | 97.55% / 97.16% | 61.05% / 90.68% | 98.02% / 99.31% | 0.830 / 0.941 |
| `3e-4` | 0.176615 | within-range interpolation | 84,437 / 252,471 | 11 / 11 | 99.54% / 99.61% | 71.00% / 96.35% | 99.51% / 99.95% | 0.842 / 0.938 |
| `3.528045626e-4` | 0.172803 | within-range interpolation | 92,944 / 281,490 | 13 / 13 | 99.77% / 99.85% | 75.22% / 97.79% | 99.77% / 99.98% | 0.850 / 0.932 |
| `4e-4` | 0.169824 | within-range interpolation | 101,434 / 309,153 | 15 / 16 | 99.96% / 99.95% | 78.53% / 98.55% | 99.91% / 100.00% | 0.848 / 0.939 |
| `5e-4` | 0.164467 | exact Gate 2 | 122,895 / 379,053 | 21 / 21 | 100.00% / 99.99% | 85.33% / 99.56% | 99.98% / 100.00% | 0.866 / 0.925 |
| `7.5e-4` | 0.154532 | within-range interpolation | 175,795 / 546,441 | 36 / 35 | 100.00% / 100.00% | 94.19% / 99.96% | 100.00% / 100.00% | 0.871 / 0.880 |
| `1e-3` | 0.147322 | exact Gate 2 | 224,608 / 699,783 | 51 / 49 | 100.00% / 100.00% | 97.81% / 99.99% | 100.00% / 100.00% | 0.884 / 0.804 |

The full trajectory is visualized in the
[13-point core-metrics plot](results_extended_2026-08-19/analysis/plots/core_metrics_vs_log10_p.png)
([SVG](results_extended_2026-08-19/analysis/plots/core_metrics_vs_log10_p.svg)).
The original result remains available as the
[nine-point core-metrics plot](results_2026-08-19/analysis/plots/core_metrics_vs_log10_p.png)
([SVG](results_2026-08-19/analysis/plots/core_metrics_vs_log10_p.svg)).
Exact source values are in the extended
[network summary](results_extended_2026-08-19/analysis/network_summary.tsv).

The trajectory separates two effects. Driver activation and SIG graph
connectivity are nearly saturated by `p=3e-4` to `4e-4`, whereas TF
incident-node coverage and target-set density continue increasing. Topology
alone therefore cannot choose among the extended points.

### Extension-point interpretation

`p=5e-4` is the first extended point at which both classes have a median of at
least 20 targets. `p=1e-3` is the only tested point at which both are near 50.
Those are density descriptions, not evidence that either point is optimal.
The SIG scale-free adjusted R-squared decreases from 0.939 at `p=4e-4` to
0.804 at `p=1e-3`, illustrating that adding density does not improve every QC
metric monotonically.

Across all 12 adjacent threshold comparisons and both driver classes, every
edge in the stricter consensus network was also present at the next looser
point. The exact comparisons are in
[adjacent_overlap.tsv](results_extended_2026-08-19/analysis/adjacent_overlap.tsv)
and the
[edge-overlap plot](results_extended_2026-08-19/analysis/plots/edge_overlap_vs_log10_p.png)
([SVG](results_extended_2026-08-19/analysis/plots/edge_overlap_vs_log10_p.svg)).
This is an observed BRCA100 result, not a general guarantee of the consensus
procedure.

## Representative NetBID2 network-QC reports

Two single-file HTML reports are included so the summarized metrics can be
inspected in NetBID2's full network-QC output:

- [`p=1e-3` TF network QC](representative_netbid2_qc/p1e-03_tf_netbid2_qc.html)
- [`p=5e-4` SIG network QC](representative_netbid2_qc/p5e-04_sig_netbid2_qc.html)

Their hashes, source consensus-network hashes, pinned environment, and source
arm paths are recorded in the
[representative-report README](representative_netbid2_qc/README.md) and
[`SHA256SUMS`](representative_netbid2_qc/SHA256SUMS). GitHub may display an
HTML file as source rather than execute it; in that case, download the raw
file and open it locally. The reports preserve the generated bytes and do not
need the original WSL result directory, but the standard R Markdown template
still attempts to load MathJax from `mathjax.rstudio.com`. These are post-DPI,
post-consensus QC reports. They are representative views, not biological
validation.

## What recurrence-qualified targets mean

For each consensus edge, support is the fraction of the same 100 seed-level
subnetworks in which that edge appears. A support-qualified target count first
discards consensus edges below a stated support threshold, then recounts
targets and zero-fills all candidate drivers. It asks how much of the apparent
regulon density is recurrent under the fixed BRCA100 subsampling design.

The table reports, for TF / SIG, the percentage of all consensus edges meeting
each support threshold and the corresponding zero-filled median targets.

| Per-seed p | Raw median | Edge share at >=20% support | Median at >=20% | Edge share at >=50% support | Median at >=50% | Edge share at >=80% support | Median at >=80% |
|---:|---:|---:|---:|---:|---:|---:|---:|
| `4e-4` | 15 / 16 | 31.24% / 24.75% | 2 / 3 | 5.96% / 3.12% | 0 / 0 | 1.36% / 0.56% | 0 / 0 |
| `5e-4` | 21 / 21 | 29.73% / 23.79% | 3 / 4 | 5.32% / 2.80% | 0 / 0 | 1.17% / 0.47% | 0 / 0 |
| `7.5e-4` | 36 / 35 | 28.31% / 22.88% | 6 / 6 | 4.62% / 2.44% | 0 / 0 | 0.95% / 0.37% | 0 / 0 |
| `1e-3` | 51 / 49 | 26.94% / 21.96% | 9 / 9 | 4.05% / 2.16% | 0 / 0 | 0.78% / 0.31% | 0 / 0 |

Absolute recurrence-qualified edge counts do increase. From `p=4e-4` to
`p=1e-3`, TF edges at >=20% support increase from 31,683 to 60,510 and SIG
edges increase from 76,524 to 153,642. Total edge counts grow faster, however:
TF consensus edges increase 121.4%, while TF edges at >=80% support increase
26.8%; SIG consensus edges increase 126.4%, while SIG edges at >=80% support
increase 24.1%.

At `p=1e-3`, 48.35% of candidate TFs and 46.45% of candidate SIGs have at
least 10 targets at >=20% support. At >=50% support, those fractions are 9.51%
and 2.21%; at >=80%, they are 1.50% and 0.06%. The raw medians near 50
therefore cannot be described as medians near 50 stable targets.

Support is not a posterior confidence score and does not estimate biological
truth or FDR. Low-recurrence edges may include sample-specific biology as well
as noise, and high recurrence does not by itself prove a direct regulatory
relationship.

## Null-burden context

After same-accession and same-symbol exclusions, each seed-level TF arm tests
73,743,244 directed candidate pairs and each SIG arm tests 301,984,614.
Multiplying those counts by the nominal p-value yields a global-independence
null exceedance proxy before DPI. At `p=1e-3`, that proxy is 73,743.2 for TF
and 301,984.6 for SIG per seed-level subnetwork.

These values are not expected post-DPI false edges, expected false consensus
edges, or empirical FDR. Expression correlation, DPI, and recurrence-based
consensus prevent that translation. They only make explicit that loosening the
per-pair threshold increases the pre-DPI null burden. The complete comparison
is in the
[coverage-versus-null-burden plot](results_extended_2026-08-19/analysis/plots/coverage_vs_nominal_null_burden.png)
([SVG](results_extended_2026-08-19/analysis/plots/coverage_vs_nominal_null_burden.svg)).

## Limitations

- This is one dataset, BRCA100, under one `m=80`, 80%-subsampling design. It
  does not establish a universal threshold across cohorts or sample sizes.
- The sweep changes only the per-subsample PR67 threshold. Consensus
  `p=1e-5`, DPI epsilon 0, estimator settings, and all other choices remain
  fixed, so their interactions were not swept.
- Gate 2 validates observable independence-null exceedance rates. It is not a
  held-out patient cohort, biological validation, or an FDR calculation.
- The topology floors were endpoint-informed. Although declared before the
  intermediate sweep was inspected, the `p=3e-4` TF incident-node result
  passed the 70% floor by only 1.00 percentage point. It is an engineering
  screen, not a biological optimum.
- No biological-reference edges, independent-cohort reproducibility, or
  downstream NetBID activity robustness were evaluated. The sweep cannot show
  that larger regulons yield more accurate or more stable activity estimates.
- No density or downstream-utility criterion was declared before the extended
  results were examined. Choosing a default now solely because it produces a
  preferred median target count would be post hoc.
- Scale-free adjusted R-squared and recurrence are QC summaries, not ground
  truth. Neither can select a biologically correct network by itself.

## Defensible conclusion

The matched evidence supports three conclusions:

1. PR67's `p=1e-7` setting causes severe BRCA100 sparsity under this design.
   Relaxing that single parameter restores driver coverage, connectivity, and
   target-set size.
2. Topology and density saturate at different rates. `p=3e-4` is the first
   tested point passing the declared topology floors, `p=5e-4` is the first
   extended point with median target count >=20 in both classes, and `p=1e-3`
   produces medians near 50. These are three different criteria.
3. Most of the density added at the loosest points has low recurrence. That is
   a warning, not proof of false edges, and it prevents a density-only result
   from being presented as a reliability improvement.

No new default is justified. A production choice requires a declared utility
criterion and the deferred downstream or biological validation. Until that is
available, the four extended points should be presented as a measured
density-stability tradeoff, not ranked as a winner.
