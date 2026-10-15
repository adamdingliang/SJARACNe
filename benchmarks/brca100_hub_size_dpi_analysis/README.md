# BRCA100 hub-list size and DPI pruning

## Question

Does the number of supplied hub genes change DPI pruning under the current
SJARACNe workflow, even when the expression matrix, target universe,
per-subsample AP-MI cutoff, samples, and other inference parameters are fixed?

The answer in this BRCA100 analysis is **yes**. Hub-list size is a material
technical covariate of DPI pruning. This does not establish that the additional
edges removed with a larger list are biologically wrong, nor does it identify
an optimal hub count.

## Why hub-list size can create technical bias

For a candidate source-target edge $A-B$, SJARACNe removes the edge when at
least one eligible intermediary $C$ completes a stronger DPI triangle. With
DPI tolerance $\epsilon_{\mathrm{DPI}}$, define

$$
T_{AB}=\frac{I(A;B)}{1-\epsilon_{\mathrm{DPI}}}.
$$

An intermediary qualifies when both

$$
I(A;C)\gt T_{AB}
\qquad\text{and}\qquad
I(B;C)\gt T_{AB},
$$

and the annotation-protection rule permits $C$. The current rule is therefore
a one-witness gate:

$$
\text{prune }A-B
\quad\Longleftrightarrow\quad
\left|\mathcal{W}_{AB}\right|\ge 1,
$$

where $\mathcal{W}_{AB}$ is the set of eligible qualifying intermediaries
available to the run.

The hub list affects $\mathcal{W}_{AB}$ through two coupled implementation
channels:

1. **Source-row availability.** The $A-C$ value comes from source row $A$.
   The $B-C$ lookup prefers row $B$, if it exists, and otherwise can fall
   back to the reverse $C-B$ value from row $C$. Only reconstructed source
   rows can supply those stored relationships. Adding hubs therefore makes
   more triangles testable and can also change direct-versus-reverse lookup
   availability.
2. **Annotation eligibility.** If either endpoint of $A-B$ is annotated, an
   unannotated intermediary cannot remove the edge. The standard workflow uses
   the same list for reconstruction (**-s**) and annotation (**-l**). Increasing
   that list therefore also converts more potential intermediaries into
   eligible annotated witnesses.

The resulting opportunity effect is directional: when one qualifying witness
is sufficient, a run with more available and eligible intermediaries has more
ways to remove a fixed candidate edge.

~~~text
Candidate edge:        A -------- B       weaker
                       \        /
                        \      /
                         C                both A-C and B-C stronger

Smaller panel: C lacks the required source-row/annotation status
               -> the triangle is unavailable or protected -> retain A-B

Larger panel:  C becomes available and eligible
               -> the same A-B edge has a qualifying witness -> prune A-B
~~~

This is a technical context effect because the decision for $A-B$ can change
when other hubs are added, without changing the samples or the pre-DPI
$A-B$ edge. The additional witness can still be biologically meaningful.
The analysis establishes sensitivity to the supplied opportunity set, not
edge-level biological error.

## Matched 100-seed pilot

The workflow-faithful pilot held the 28,278-gene BRCA100 expression and target
universe fixed. It used deterministic nested TF and SIG panels, the same
80-of-100 sample indices at each matched seed, $N_{\mathrm{par}}=40$,
$\epsilon_{\mathrm{DPI}}=0$, and the operating-point AP-MI filters from PR
#72. The same panel was supplied to **-s** and **-l**.

| Network | Small hubs | Middle hubs | Full hubs | Median pruning fraction over 100 seeds, small / middle / full | Difference of arm medians, full minus small |
|---|---:|---:|---:|---:|---:|
| TF | 326 | 1,304 | 2,608 | 0.1079 / 0.1917 / 0.2420 | +0.1341 |
| SIG | 1,335 | 5,340 | 10,680 | 0.2249 / 0.3432 / 0.4166 | +0.1917 |

Every adjacent-size paired difference was positive in all 100 matched seeds.
The median increases were:

- TF, 326 to 1,304: +8.243 percentage points;
- TF, 1,304 to 2,608: +5.115 points;
- SIG, 1,335 to 5,340: +11.749 points;
- SIG, 5,340 to 10,680: +7.555 points.

These are not small fluctuations around zero. Under this one nested panel
trajectory, the list size materially changes how much of each pre-DPI network
is removed.

## Fixed-source mechanism audit

The pilot changes both the evaluated source population and the surrounding
hub/annotation panel. A separate ten-seed SIG audit therefore held the
evaluated source population fixed at the same 1,335 genes while expanding the
available reconstructed/annotated panel to 1,335, 5,340, and 10,680 hubs.

A **source row** contains the candidate edges from one source gene to its
targets. At every panel size, the audit counts only edges from the original
1,335 source genes; their targets still span the same expression-gene universe.
This is not a subnetwork restricted to 1,335 genes at both endpoints. Added
source rows remain available to the DPI algorithm, but their own edge counts
are excluded from both the numerator and denominator of the reported fraction.

For each matched seed $s$ and available panel size $H$, calculate

$$
r_{s,H}=\frac{\text{edges pruned from the fixed 1,335 source rows}}
{\text{pre-DPI edges in those same 1,335 source rows}}.
$$

This is a ratio of summed edge counts, not a median of per-gene pruning
fractions. The table then takes the median of these per-seed ratios over seeds
1 through 10, using the native $K_{\mathrm{DPI}}=1$ rule:

| Available SIG panel | Fixed evaluated sources | Median pruning fraction over 10 seeds |
|---:|---:|---:|
| 1,335 | 1,335 | 0.2264 |
| 5,340 | 1,335 | 0.3372 |
| 10,680 | 1,335 | 0.3959 |

The two sections use different seed windows. The primary pilot value 0.2249
for 1,335 SIG hubs is the median over seeds 1 through 100
(0.2249030821). The fixed-source value 0.2264 is the median over seeds 1
through 10 (0.2264073464). Recomputing the primary pilot from only those first
ten seeds gives exactly 0.2264073464. At this smallest panel, the fixed 1,335
sources are also the complete source set, and the audit exactly reproduced the
pilot's per-seed pre-DPI counts, pruned counts, and fractions. The numerical
difference is therefore caused only by summarizing 10 rather than 100 seeds,
not by a different denominator or DPI calculation.

### What the 26,700 comparisons check

For each source gene and seed, the audit compares that row's pre-DPI edge count
in each larger panel with its count in the smallest panel:

$$
1{,}335\ \text{sources}\times 10\ \text{seeds}\times
2\ \text{larger-versus-small comparisons}=26{,}700.
$$

All comparisons had identical counts, with zero mismatches. Therefore the
total evaluated pre-DPI edge count also stayed fixed across panel sizes within
each seed. Counts can differ between seeds; these are matched count checks,
not 26,700 independent experiments.

For example, the recorded **seed 1** results are:

| Available SIG panel | Evaluated source rows | Pre-DPI edges in those rows | Pruned edges in those rows | Pruning fraction |
|---:|---:|---:|---:|---:|
| 1,335 | 1,335 | 55,894 | 12,302 | 22.01% |
| 5,340 | 1,335 | 55,894 | 18,222 | 32.60% |
| 10,680 | 1,335 | 55,894 | 21,583 | 38.61% |

The denominator is unchanged, but more edges from the original source rows
are pruned as the surrounding panel expands.

### What the 0.1692 increase means

First calculate the paired full-minus-small difference within each seed:

$$
\Delta_s=r_{s,10680}-r_{s,1335}.
$$

Then take the median of those ten differences: 0.169159, rounded to 0.1692,
or **16.92 percentage points**. This is not a 16.92% relative increase and is
not calculated by subtracting the two separately reported arm medians. The
full-panel pruning fraction exceeded the small-panel fraction in all ten seeds.

With an all-source summary, a larger panel could change the aggregate pruning
fraction simply by adding source genes with different pruning rates, even if
pruning in the original rows were unchanged. This fixed-source audit excludes
that change in the measured source population as the sole explanation:
the original rows themselves have higher aggregate pruning, with the same
pre-DPI edge count within each matched seed. This supports **panel-context
sensitivity in the fixed source population** under this nested-panel design.

The audit checked pre-DPI edge **counts**, not target identities or MI values
across panels. Equal counts do not establish identical candidate-edge sets.
Source-row availability and annotation eligibility also change together, so
the audit does not apportion their effects or establish whether the additional
removed edges are biologically correct or incorrect.

The example and paired summary are recorded in
[seed-level metrics](../brca100_kdpi_witness_screen/results_2026-08-26/seed_level_metrics.tsv)
and
[paired hub-size effects](../brca100_kdpi_witness_screen/results_2026-08-26/paired_hub_size_effects.tsv).

## Relationship to consensus recurrence

Consensus recurrence $K_{\mathrm{edge}}$ is applied after DPI and is not the
cause of this effect. The pilot independently aggregated support
$K_{\mathrm{edge}}\ge 6$ across the 100 post-DPI networks as a supplementary
topology check:

| Network | Hubs, small / middle / full | $K_{\mathrm{edge}}\ge 6$ edges | Zero-filled median targets |
|---|---:|---:|---:|
| TF | 326 / 1,304 / 2,608 | 61,621 / 227,784 / 416,408 | 115 / 113 / 110.5 |
| SIG | 1,335 / 5,340 / 10,680 | 127,581 / 424,006 / 739,958 | 48 / 49 / 47 |

The full-size counts exactly reproduced the independent anchors used in PR
#73: 416,408 TF edges and 739,958 SIG edges. This package used an independent
benchmark aggregator rather than the implementation proposed in PR #73.
Consensus filtering cannot recover edges that DPI already removed, so
$K_{\mathrm{edge}}$ does not correct the upstream panel-size dependence.

## What the witness-threshold screen showed

The diagnostic also counted how many eligible DPI witnesses supported each
removal. Requiring a fixed $K_{\mathrm{DPI}}\gt 1$ did not mitigate the
hub-size effect:

- $K_{\mathrm{DPI}}=2,3,5$ increased the fixed-source full-minus-small gap;
- $K_{\mathrm{DPI}}=10$ reduced the absolute gap only by retaining 1.66% of
  small-panel pruning while retaining 40.9% in the full panel.

A fixed witness-count threshold is therefore rejected as a hub-size-bias
correction in this screen. It can make DPI more conservative for a fixed panel,
but that is a different objective.

## Why call this a technical-bias risk

The observed size trajectory is systematic rather than seed noise. In the
fixed-source audit, expanding the surrounding panel changed pruning while the
samples, evaluated source population, and every matched source-row pre-DPI edge
count remained fixed. Therefore hub-list size is a technical exposure of the
current method. Because edge identities and MI values were not independently
hashed across panels, this audit should not be read as proving that the input
candidate-edge sets were byte-identical.

Calling every extra removal a biological error would overstate the evidence.
Some added intermediaries may reveal genuine indirect relationships. The bias
risk is that larger panels receive systematically more opportunities to find
at least one qualifying witness, so pruning fractions and retained networks
are not directly comparable across panels as if hub size were biologically
neutral.

In particular, raw TF-versus-SIG pruning differences cannot be attributed only
to biology: the two networks use different hub universes and per-subsample
filters as well as different list sizes.

## Defensible conclusion

This analysis supports the following limited conclusions:

- record hub-list size and composition as part of network provenance;
- treat DPI pruning fraction as panel-dependent, not an intrinsic property of
  the dataset;
- use matched hub-size sensitivity checks before comparing networks built from
  substantially different regulator universes;
- do not infer that extra pruning from a larger panel is biologically correct
  or incorrect without external evidence;
- do not claim an optimal hub count or extrapolate these whole-transcriptome
  BRCA100 results to Xenium 5k.

The analysis identifies a technical dependency. It does not yet provide a
validated correction. Dynamic $\epsilon_{\mathrm{DPI}}$ testing and
edge-level biological adjudication are intentionally outside this PR.

## Audit instrumentation

The native changes in this branch expose measurements; they do not change the
default DPI decision rule:

- every run emits one machine-readable **DPI_STATS** record containing exact
  pre-DPI, pruned, and retained edge counts;
- optional **-W** writes per-source counts of edges with at least 1, 2, 3, 5,
  or 10 qualifying witnesses;
- the witness path fails closed on invalid DPI state, exact CLI-path collisions
  with known input/output files, malformed accounting, or any disagreement
  between the reconstructed one-witness counts and native pruning marks.

The path-collision guard compares the path strings supplied to the CLI; it does
not resolve lexical aliases or symbolic links. The diagnostic should therefore
be written to a dedicated output path.

All 30 witness-audit runs exactly reproduced the frozen native
$K_{\mathrm{DPI}}=1$ sample selection, DPI totals, and adjacency data. The
sidecar is diagnostic output, not a production $K_{\mathrm{DPI}}\gt 1$ pruning
option.

## Evidence and provenance

- [100-seed pilot design and runner](../brca100_hub_size_dpi_pilot/)
- [Compact pilot evidence](../brca100_hub_size_dpi_pilot/results_2026-08-25/)
- [Witness diagnostic and ten-seed fixed-source screen](../brca100_kdpi_witness_screen/)
- [Compact witness-screen evidence](../brca100_kdpi_witness_screen/results_2026-08-26/)

All required validation gates passed in both analyses. The pilot contains 600
matched inference runs; the witness screen contains 30 matched candidate runs
and 300 source-group/threshold analysis rows. The frozen manifests record the
exact source snapshots, binaries, panels, samples, inputs, commands, and
artifact hashes. Large adjacency files and logs remain in the recorded external
WSL work roots and are deliberately excluded from Git.
