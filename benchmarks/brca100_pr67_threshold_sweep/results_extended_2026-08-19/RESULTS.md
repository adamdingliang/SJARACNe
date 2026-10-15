# BRCA100 PR67 extended per-subsample threshold sweep

## Bottom line

The extended sweep did what it was intended to do mechanically: looser PR67
per-subsample thresholds increased regulon size substantially. Median
zero-filled target counts rose from 15 TF / 16 SIG at `p=4e-4`, to 21 / 21 at
`p=5e-4`, 36 / 35 at `p=7.5e-4`, and 51 / 49 at `p=1e-3`. Of these four
points, `p=5e-4` is the first at which both classes reached a median of at
least 20 targets.

The recurrence-qualified result is much less favorable. When targets are
restricted to edges present in at least 20 of 100 seed-level subnetworks, the
TF / SIG medians rise only from 2 / 3 at `p=4e-4` to 9 / 9 at `p=1e-3`.
At support thresholds of 50% and 80%, the zero-filled median remains zero for
both driver classes at every new point. Thus, most of the added density is
carried by low-recurrence edges.

That is not proof that the added edges are false: recurrence measures
subsample stability, not biological truth, an edge probability, or FDR. But
there is also no evidence here that the low-recurrence additions improve
activity estimation. Biological-reference and downstream activity-robustness
tests were deliberately deferred.

No new default is supported by this sweep. No density or downstream-utility
criterion was declared, so there is no defensible rule for choosing among
`4e-4`, `5e-4`, `7.5e-4`, and `1e-3`.

## Matched design and provenance

Only the seed-level PR67 p-value varied. Every arm used:

- PR67 commit `7633ebb4a0d966dbda15a4e32d0efa492fb71aeb`;
- the estimator-matched `m=80`, `Npar=40` null model;
- fixed 80% sampling without replacement;
- BRCA100, with 2,608 candidate TFs, 10,680 candidate SIGs, and 28,278
  expression nodes;
- seeds 1 through 100;
- adaptive-partitioning MI, `Npar=40`, and DPI epsilon 0; and
- consensus-network `p=1e-5`.

All four new probabilities were inferred directly for TF and SIG at every
seed: 4 points x 2 driver classes x 100 seeds = 800 new seed-level networks.
They were not manufactured by filtering a looser post-DPI adjacency file.
Together with the nine earlier points, the extended artifact contains 2,600
validated inference jobs and 26 matched consensus/support/NetBID2 arms.

The sweep design was extended append-only. The prior design SHA-256 is
archived with an explicit migration record under
[`provenance/sweep_design_history`](provenance/sweep_design_history), and the
active design SHA-256 is
`2eca0f9ef388fdf48632f8292c226bac59e0e325d9641574b314689651c6a8e5`.
The completed job ledger is
[`provenance/run_manifest.tsv`](provenance/run_manifest.tsv).

## Calibration status

All four new probabilities lie inside the independently accepted held-out
range, `2e-5` through `2e-3`, and none uses GPD-tail extrapolation.

- `p=5e-4` and `p=1e-3` are exact Gate-2 grid points: those probabilities
  were evaluated directly in the independent held-out null stream.
- `p=4e-4` and `p=7.5e-4` are interpolations inside the accepted range. They
  inherit within-range calibration status but were not evaluated as direct
  Gate-2 grid probabilities.

Calibration class describes null-calibration provenance. It does not rank the
biological quality of the resulting networks.

## Four-point density and topology

Active is the percentage of candidate drivers with at least one consensus
edge. Incident is the percentage of all 28,278 expression nodes incident to
an edge. LCC is the percentage of incident nodes in the largest weak
component. Median targets are zero-filled over all candidate drivers.
Adjusted R-squared is the NetBID2 scale-free-fit diagnostic, not a selection
criterion.

| Per-seed p | Class | MI cutoff | Calibration | Edges | Active | Incident | LCC | Median targets | NetBID2 adjusted R-squared |
|---:|---|---:|---|---:|---:|---:|---:|---:|---:|
| `4e-4` | TF | 0.169824 | interpolation | 101,434 | 99.96% | 78.53% | 99.91% | 15 | 0.848 |
| `4e-4` | SIG | 0.169824 | interpolation | 309,153 | 99.95% | 98.55% | 100.00% | 16 | 0.939 |
| `5e-4` | TF | 0.164467 | exact Gate 2 | 122,895 | 100.00% | 85.33% | 99.98% | 21 | 0.866 |
| `5e-4` | SIG | 0.164467 | exact Gate 2 | 379,053 | 99.99% | 99.56% | 100.00% | 21 | 0.925 |
| `7.5e-4` | TF | 0.154532 | interpolation | 175,795 | 100.00% | 94.19% | 100.00% | 36 | 0.871 |
| `7.5e-4` | SIG | 0.154532 | interpolation | 546,441 | 100.00% | 99.96% | 100.00% | 35 | 0.880 |
| `1e-3` | TF | 0.147322 | exact Gate 2 | 224,608 | 100.00% | 97.81% | 100.00% | 51 | 0.884 |
| `1e-3` | SIG | 0.147322 | exact Gate 2 | 699,783 | 100.00% | 99.99% | 100.00% | 49 | 0.804 |

Driver activation and SIG topology are already close to saturated at `4e-4`.
The main topology gain at looser thresholds is TF incident-node coverage,
which rises from 78.53% to 97.81%. The SIG scale-free adjusted R-squared moves
in the opposite direction, from 0.939 at `4e-4` to 0.804 at `1e-3`; density is
not a monotonic improvement on every QC metric.

Every stricter network was an edge subset of its next looser neighbor in these
BRCA100 results. That observed nesting is documented in
[`adjacent_overlap.tsv`](analysis/adjacent_overlap.tsv); it is not a general
guarantee of the consensus procedure. The full metric panel is
[`core_metrics_vs_log10_p.png`](analysis/plots/core_metrics_vs_log10_p.png),
with an [SVG version](analysis/plots/core_metrics_vs_log10_p.svg).

## Candidate-driver target coverage

The table reports the percentage of all candidate drivers with at least the
stated number of consensus targets. Drivers with no edges are included as
zeros.

| Per-seed p | Class | Median | >=10 | >=20 | >=30 | >=50 | >=100 |
|---:|---|---:|---:|---:|---:|---:|---:|
| `4e-4` | TF | 15 | 68.6% | 42.0% | 30.6% | 20.8% | 10.7% |
| `4e-4` | SIG | 16 | 69.2% | 42.9% | 31.5% | 17.6% | 4.6% |
| `5e-4` | TF | 21 | 84.5% | 52.7% | 37.7% | 24.4% | 13.4% |
| `5e-4` | SIG | 21 | 84.6% | 53.4% | 39.0% | 22.9% | 6.5% |
| `7.5e-4` | TF | 36 | 98.8% | 80.7% | 59.2% | 37.1% | 19.5% |
| `7.5e-4` | SIG | 35 | 99.0% | 80.6% | 58.3% | 36.3% | 12.6% |
| `1e-3` | TF | 51 | 100.0% | 96.7% | 81.6% | 52.1% | 25.3% |
| `1e-3` | SIG | 49 | 100.0% | 96.7% | 80.2% | 49.6% | 18.8% |

The source values, including exact counts and unrounded fractions, are in
[`network_summary.tsv`](analysis/network_summary.tsv).

## Recurrence-qualified target coverage

Support is the fraction of the 100 matched seed-level subnetworks containing
a consensus edge. For each support threshold below, `edge share` is the
percentage of consensus edges meeting that threshold; `median` and the target
coverage columns are then recomputed after discarding lower-support edges and
zero-filling all candidate drivers.

### TF

| Per-seed p | Edge support | Edge share | Median | >=10 | >=20 | >=30 | >=50 | >=100 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `4e-4` | >=20% | 31.24% | 2 | 23.70% | 16.18% | 11.73% | 6.52% | 2.30% |
| `4e-4` | >=50% | 5.96% | 0 | 6.37% | 3.03% | 1.76% | 0.69% | 0.19% |
| `4e-4` | >=80% | 1.36% | 0 | 1.27% | 0.35% | 0.23% | 0.08% | 0.00% |
| `5e-4` | >=20% | 29.73% | 3 | 27.80% | 18.56% | 13.57% | 7.75% | 2.72% |
| `5e-4` | >=50% | 5.32% | 0 | 6.79% | 3.14% | 1.76% | 0.73% | 0.19% |
| `5e-4` | >=80% | 1.17% | 0 | 1.27% | 0.35% | 0.23% | 0.08% | 0.00% |
| `7.5e-4` | >=20% | 28.31% | 6 | 38.34% | 23.81% | 18.17% | 11.00% | 3.87% |
| `7.5e-4` | >=50% | 4.62% | 0 | 8.78% | 3.99% | 2.07% | 0.81% | 0.19% |
| `7.5e-4` | >=80% | 0.95% | 0 | 1.42% | 0.46% | 0.23% | 0.08% | 0.00% |
| `1e-3` | >=20% | 26.94% | 9 | 48.35% | 29.10% | 21.82% | 13.50% | 4.52% |
| `1e-3` | >=50% | 4.05% | 0 | 9.51% | 4.52% | 2.22% | 0.84% | 0.19% |
| `1e-3` | >=80% | 0.78% | 0 | 1.50% | 0.54% | 0.23% | 0.08% | 0.00% |

### SIG

| Per-seed p | Edge support | Edge share | Median | >=10 | >=20 | >=30 | >=50 | >=100 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `4e-4` | >=20% | 24.75% | 3 | 22.95% | 10.65% | 5.24% | 1.42% | 0.05% |
| `4e-4` | >=50% | 3.12% | 0 | 1.38% | 0.19% | 0.06% | 0.00% | 0.00% |
| `4e-4` | >=80% | 0.56% | 0 | 0.06% | 0.00% | 0.00% | 0.00% | 0.00% |
| `5e-4` | >=20% | 23.79% | 4 | 27.06% | 12.84% | 6.46% | 1.72% | 0.07% |
| `5e-4` | >=50% | 2.80% | 0 | 1.51% | 0.19% | 0.06% | 0.00% | 0.00% |
| `5e-4` | >=80% | 0.47% | 0 | 0.06% | 0.00% | 0.00% | 0.00% | 0.00% |
| `7.5e-4` | >=20% | 22.88% | 6 | 37.62% | 18.97% | 10.28% | 2.88% | 0.10% |
| `7.5e-4` | >=50% | 2.44% | 0 | 2.06% | 0.26% | 0.06% | 0.00% | 0.00% |
| `7.5e-4` | >=80% | 0.37% | 0 | 0.06% | 0.00% | 0.00% | 0.00% | 0.00% |
| `1e-3` | >=20% | 21.96% | 9 | 46.45% | 24.08% | 13.55% | 4.08% | 0.12% |
| `1e-3` | >=50% | 2.16% | 0 | 2.21% | 0.28% | 0.06% | 0.00% | 0.00% |
| `1e-3` | >=80% | 0.31% | 0 | 0.06% | 0.00% | 0.00% | 0.00% | 0.00% |

The absolute number of recurrence-qualified edges does rise as `p` is
loosened, but total edge count rises faster. For example, TF edges with at
least 20% support increase from 31,683 at `4e-4` to 60,510 at `1e-3`, while
their share falls from 31.24% to 26.94%. The SIG count increases from 76,524
to 153,642 while its share falls from 24.75% to 21.96%.

At the loosest point, only 48.35% of TFs and 46.45% of SIGs have at least 10
targets at 20% support. At 50% support those values are 9.51% and 2.21%; at
80% support they are 1.50% and 0.06%. The raw increase to medians of roughly
50 targets therefore should not be described as a similarly large increase
in stable targets.

## Null-burden proxy

After SJARACNe's same-accession and same-symbol exclusions, each seed-level TF
network tests 73,743,244 directed candidate pairs and each SIG network tests
301,984,614. Multiplying those exact test counts by `p` gives the following
global-independence-null expectation before DPI:

| Per-seed p | TF nominal exceedances | SIG nominal exceedances |
|---:|---:|---:|
| `4e-4` | 29,497.3 | 120,793.8 |
| `5e-4` | 36,871.6 | 150,992.3 |
| `7.5e-4` | 55,307.4 | 226,488.5 |
| `1e-3` | 73,743.2 | 301,984.6 |

These are not expected consensus false edges, post-DPI quantities, or an
empirical FDR. They expose the increasing pre-DPI null burden; correlation,
DPI, and recurrence-based consensus prevent translating them directly into a
false-edge count. The corresponding plot is
[`coverage_vs_nominal_null_burden.png`](analysis/plots/coverage_vs_nominal_null_burden.png),
with an [SVG version](analysis/plots/coverage_vs_nominal_null_burden.svg).

## Why the topology selector still returns `p=3e-4`

The existing automated selector chooses the smallest probability inside the
held-out range that passes the predeclared active-driver, incident-node, and
largest-component floors for both TF and SIG. `p=3e-4` already passes those
topology floors, so adding looser points cannot displace it under that rule.

That selector is explicitly not a regulon-density rule. Its retained
`p=3e-4` result does not answer the present target-size question and should
not be treated as evidence against, or in favor of, one of the four extended
points. The rule and its endpoint-informed declaration timing are recorded in
[`selection.json`](analysis/selection.json).

## Integrity checks

The compact package records:

- 2,600 seed-level jobs: 13 thresholds x 2 driver classes x 100 seeds;
- 26 consensus, support-summary, and NetBID2 summary arms;
- 400 successful seed-level anchor comparisons reproducing the prior PR67
  default and PR66 cutoff-match data sections;
- the archived v1 design, v1-to-v2 design-migration ledger, and migrated
  NetBID2 manifest history; and
- hashes and sizes for 13,464 omitted raw artifacts in
  [`omitted_artifacts.json`](omitted_artifacts.json).

Large seed-level and consensus artifacts are intentionally omitted from this
repository package. [`SHA256SUMS`](SHA256SUMS),
[`package_manifest.json`](package_manifest.json), and the per-arm provenance
records authenticate the compact evidence set.

## Limitations and decision status

- This is one BRCA100 benchmark, not a multi-cohort replication.
- The swept `p` is a per-pair seed-level null threshold. It is not an edge FDR
  or FWER guarantee.
- Support is a stability diagnostic, not a truth label. Low-support edges may
  still help activity estimation, but that has not been tested here.
- Biological-reference accuracy, independent-cohort reproducibility, and
  downstream TF/SIG activity robustness were not evaluated.
- The consensus `p=1e-5` and all non-threshold inference settings were held
  fixed; this sweep does not establish how conclusions change under another
  consensus rule.
- Scale-free adjusted R-squared is only one NetBID2 QC statistic. In
  particular, the SIG value declines as these networks become denser.
- No minimum useful target count, maximum acceptable null burden, or
  reliability-weighted density objective was declared before this sweep.

The defensible conclusion is narrow: `5e-4` is the first of the four tested
extensions to reach median target size 21 in both driver classes, while
`7.5e-4` and `1e-3` provide much larger raw regulons at the cost of a growing
low-recurrence edge majority and higher nominal null burden. Choosing a new
default requires the deferred downstream robustness or biological validation,
or an explicitly declared engineering tradeoff. This sweep alone does not
choose one.
