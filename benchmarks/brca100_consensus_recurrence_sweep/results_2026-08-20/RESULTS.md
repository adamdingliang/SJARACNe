# BRCA100 consensus-recurrence sweep

## Outcome

This benchmark does **not** statistically or biologically identify a unique
consensus cutoff. It quantifies the engineering trade-off after fixing the
per-subsample networks at the two denser PR67 operating points selected for
follow-up. The later software decision proposes `K=6` as a density-favoring
engineering rule; this report records the supporting evidence and its limits:

- TF: per-subsample `p_b=1e-3`, AP-MI cutoff `0.14732247558240297`;
- SIG: per-subsample `p_b=5e-4`, AP-MI cutoff `0.1644671599536221`;
- 100 fixed BRCA100 seeds, `m=80`, `Npar=40`, fixed-size sampling without
  replacement, and DPI epsilon 0.

Only the post-DPI recurrence requirement changes. We denote this consensus-edge
threshold by `K_edge` to distinguish it from proposed DPI-specific recurrence
parameters; tables and filenames abbreviate it as `K`. An ordered edge is
retained when it occurs in at least `K_edge` of the 100 fixed subnetworks.
The sweep covers every integer `K_edge=6..20`; the historical consensus setting
`p_c=1e-5` maps to `K=9` in both arms.

Loosening from `K=9` to `K=8` is the smallest tested change. It adds 44,686 TF
and 83,046 SIG edges, all with support exactly 8/100. Median zero-filled target
counts increase from 51 to 65 for TF and from 21 to 27 for SIG. Loosening to
`K=6` nearly doubles the `K=9` edge count and raises the medians to 110.5 TF
and 47 SIG, but it admits edges with recurrence as low as 6/100. That is a
real loss of resampling stability, not free density.

## Full trajectory

`Drivers >=20` is zero-filled over all 2,608 TF or 10,680 SIG candidates.
The plug-in Poisson-binomial (PB) tail is an untruncated model reference, not
an FDR or probability that an edge is true. NetBID2 adjusted R2 is descriptive
and was not optimized. TF and SIG have different candidate-driver lists and
different per-subsample thresholds, so their raw edge totals should not be
compared as though they shared one opportunity set; the valid comparison is
within each arm as `K` changes.

| K | TF: edges / median targets / drivers >=20 | SIG: edges / median targets / drivers >=20 | Plug-in PB tail TF / SIG | NetBID2 adj. R2 TF / SIG |
|---:|---:|---:|---:|---:|
| 6 | 416,408 / 110.5 / 100.0% | 739,958 / 47 / 93.8% | 0.0327 / 0.0177 | 0.879 / 0.864 |
| 7 | 330,184 / 83 / 100.0% | 575,657 / 35 / 79.2% | 0.0101 / 0.00478 | 0.880 / 0.895 |
| 8 | 269,294 / 65 / 99.7% | 462,099 / 27 / 64.1% | 0.00276 / 0.00113 | 0.887 / 0.914 |
| 9 | 224,608 / 51 / 96.7% | 379,053 / 21 / 53.4% | 0.000666 / 0.000237 | 0.884 / 0.925 |
| 10 | 190,551 / 42 / 89.2% | 316,865 / 17 / 45.3% | 0.000144 / 4.43e-5 | 0.893 / 0.931 |
| 11 | 163,819 / 35 / 79.7% | 268,943 / 14 / 39.7% | 2.81e-5 / 7.50e-6 | 0.893 / 0.936 |
| 12 | 142,500 / 29 / 68.1% | 231,273 / 11 / 34.8% | 5.00e-6 / 1.15e-6 | 0.892 / 0.935 |
| 13 | 125,055 / 25 / 59.4% | 200,728 / 10 / 30.4% | 8.13e-7 / 1.62e-7 | 0.887 / 0.942 |
| 14 | 110,774 / 21 / 52.3% | 175,654 / 8 / 27.0% | 1.22e-7 / 2.10e-8 | 0.891 / 0.939 |
| 15 | 98,786 / 18 / 46.4% | 154,946 / 7 / 23.7% | 1.68e-8 / 2.51e-9 | 0.896 / 0.944 |
| 16 | 88,583 / 15 / 42.1% | 137,739 / 6 / 21.0% | 2.16e-9 / 2.78e-10 | 0.903 / 0.942 |
| 17 | 80,006 / 13.5 / 38.2% | 122,976 / 5 / 18.6% | 2.58e-10 / 2.87e-11 | 0.894 / 0.945 |
| 18 | 72,699 / 12 / 34.4% | 110,389 / 5 / 16.6% | 2.87e-11 / 2.77e-12 | 0.896 / 0.947 |
| 19 | 66,169 / 10 / 32.1% | 99,526 / 4 / 14.7% | 3.00e-12 / 2.50e-13 | 0.906 / 0.948 |
| 20 | 60,510 / 9 / 29.1% | 90,179 / 4 / 12.8% | 2.95e-13 / 2.12e-14 | 0.904 / 0.947 |

The complete machine-readable results include target-count quartiles,
fractions with at least 1/10/20/50/100 targets, incident nodes, weak
components, giant-component coverage, MI summaries, support summaries, and
the two tail calculations:

- [network summary](analysis/network_summary.tsv)
- [driver target coverage](analysis/driver_target_coverage.tsv)
- [density and coverage plot](analysis/plots/recurrence_density_coverage.png)

## Why a nominal consensus p-value is misleading here

For the fixed observed union `U`, the historical code sets
`q_i=E_i/U`, calculates a normal approximation to the recurrence count, and
filters its upper tail. Under the *same* plug-in occupancy assumptions, the
exact Poisson-binomial tail is much larger at `K=9`:

| Quantity | TF | SIG |
|---|---:|---:|
| Observed union `U` | 5,666,377 | 13,243,781 |
| Mean recurrence `mu` | 2.3812 | 2.0618 |
| SD recurrence `sigma` | 1.5246 | 1.4210 |
| Legacy normal tail at `K=9` | 7.08e-6 | 5.23e-7 |
| Exact plug-in PB tail at `K=9` | 6.66e-4 | 2.37e-4 |
| First `K` with plug-in PB tail below `1e-5` | 12 | 11 |
| First plug-in Bonferroni `K` (`0.05/U`) | 16 | 15 |

Thus `p_c=1e-5` does not represent a `1e-5` recurrence tail under the plug-in
model used by the historical code. At `K=9`, the normal approximation
understates the corresponding exact Poisson-binomial tail by about 94-fold for
TF and 452-fold for SIG. The first exact plug-in tail below `1e-5` occurs at
`K=12` for TF and `K=11` for SIG, where the networks are even sparser. These
are model-reference boundaries, not recommended cutoffs.

Replacing the normal approximation with an exact Poisson-binomial calculation
would fix that numerical approximation, but it would not produce a calibrated
edge p-value or FDR. The `q_i` values are plug-in occupancies estimated from the
same observed networks and union, and the inference runs share samples and
data, so the exchangeability and independence assumptions are not established.
Filtering directly on the observable rule `S_e >= K_edge` is therefore more
transparent, deterministic, and auditable than relabeling `K` through this
probability model. It remains an engineering stability rule, so reports should
state both `K_edge` and the number of subsamples `B` (here, `K_edge=6` of
`B=100`).

## High-recurrence cores

Loosening `K` only adds lower-recurrence edges; it does not remove the stable
core. The absolute core sizes are invariant for every displayed network:

| Driver class | Support >=20/100 | Support >=50/100 | Support >=80/100 |
|---|---:|---:|---:|
| TF | 60,510 | 9,100 | 1,752 |
| SIG | 90,179 | 10,597 | 1,800 |

This is why absolute core counts are more informative than their fraction of
the total network, which mechanically falls when lower-support edges are
added.

## Defensible conclusion

The sweep establishes a deterministic, nested density-versus-recurrence curve
for these two BRCA100 arms. It supports the following density-favoring operating
points for the BRCA100 development workflow; it does not prove biological
correctness or a globally optimal cutoff.

| Network | Per-subsample `p_b` | AP-MI cutoff (`m=80`, `Npar=40`) | Consensus |
|---|---:|---:|---:|
| TF | `1e-3` | `0.147322` | `K=6` |
| SIG | `5e-4` | `0.164467` | `K=6` |

The `p_b` threshold is applied to each subsample network before DPI; `K` is
applied afterward to the post-DPI recurrence count across the 100 fixed
subsamples. At `K=6`, median target counts are 110.5 for TF and 47 for SIG,
all TF drivers and 93.8% of SIG drivers have at least 20 targets, and the
NetBID2 adjusted R2 values are 0.879 and 0.864. Those density and topology
metrics are compatible with the intended activity-calculation use case.

The cost is explicit: `K=6` admits edges observed in only 6% of the subsamples,
and its plug-in Poisson-binomial tails are 0.0327 for TF and 0.0177 for SIG.
The choice is therefore an engineering density-versus-stability trade-off, not
a significance threshold, FDR guarantee, biological validation, or universal
optimum. Although the AP-MI cutoff mappings came from a held-out null range,
selection of these operating points used BRCA100 network behavior, making
BRCA100 a development dataset rather than independent validation.

The generic software default for `p_b` remains `1e-7` for retrospective
compatibility. TF and SIG use the proposed `p_b` values through separate,
explicit invocations; the software does not infer a network type or silently
switch defaults. Downstream biological-reference and activity-robustness
validation remain deferred.

## Validation and cost

The one-pass aggregator read each of the 200 source adjacency files once. All
30 materialized edge sets were nested, and the `K=9` TF/SIG networks and
NetBID2 summaries reproduced the previous outputs byte-for-byte. Independent
audits recomputed all support counts, exact tails, target coverage, topology,
and 90 NetBID2 output hashes. All 30 NetBID2 stderr logs were empty; the two
aggregator stderr logs contained progress messages only.

The complete run took about 54 minutes on the 16-CPU, 62-GB WSL host without
rerunning SJARACNe inference: about 73 seconds for aggregation, 22 minutes for
network materialization/annotation, 26 minutes for 30 serial NetBID2 summaries,
and roughly 1.5 minutes for the joined analysis.

The package preserves the compact tables, plots, manifests, scripts, pinned
R/NetBID2 environment records, hashes for omitted large networks, frozen
aggregate/exact-tail records, and source `K=9` anchor manifests documenting the
byte-for-byte reproduction.
