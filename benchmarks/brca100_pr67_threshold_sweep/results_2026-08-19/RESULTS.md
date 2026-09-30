# BRCA100 PR67 per-subsample threshold sweep

## Bottom line

The original PR67 setting, `p=1e-7`, is too stringent for BRCA100 if the goal
is to preserve broad driver coverage and a connected network for NetBID2. The
first tested operating point that passed the endpoint-informed topology floors
for both TF and SIG was `p=3e-4`, corresponding to an AP-MI cutoff of
`0.176615102822`.

That result is deliberately labeled **provisional**:

- `p=3e-4` was interpolated inside the accepted calibration range, rather than
  being one of the exact probabilities used in Gate-2 held-out testing;
- its limiting TF incident-node coverage was 71.00%, only 1.00 percentage point
  above the declared 70% floor;
- the floors were chosen after seeing the earlier PR66/PR67 endpoints, although
  before seeing the intermediate sweep results; and
- this is BRCA100 topology and NetBID2 QC, not biological validation or an
  empirical false-discovery-rate analysis.

The result therefore supports refining the interval between `2e-4` and `3e-4`.
It does **not** justify changing PR67's global default directly to `3e-4`.

## Matched design

Only the seed-level PR67 p-value varied. Every arm used:

- PR67 commit `7633ebb4a0d966dbda15a4e32d0efa492fb71aeb`;
- the estimator-matched `m=80`, `Npar=40` model;
- fixed 80% sampling without replacement;
- BRCA100, with separate TF and SIG driver lists;
- seeds 1 through 100;
- adaptive-partitioning MI, `Npar=40`, and DPI epsilon 0; and
- consensus-network `p=1e-5`.

All nine probabilities were inferred directly. No network was manufactured by
filtering a looser adjacency file after DPI.

Gate 2 refers to an independent held-out stream of canonical rank permutations
under the AP-MI independence null. It is not a held-out BRCA cohort. The exact
Gate-2 probabilities in this sweep were `2e-5`, `5e-5`, `1e-4`, and `2e-4`.
The two looser points were interpolations inside the accepted probability
range.

## Topology results

Active is the percentage of candidate drivers with at least one retained edge.
Incident is the percentage of all 28,278 expression nodes incident to an edge.
LCC is the percentage of incident nodes in the largest weak component.

| Per-seed p | MI cutoff | Calibration class | TF edges | TF active | TF incident | TF LCC | SIG edges | SIG active | SIG incident | SIG LCC |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `1e-7` | 0.322465 | extrapolated below held-out range | 4,672 | 10.3% | 6.1% | 77.3% | 13,183 | 13.9% | 11.6% | 82.3% |
| `1e-6` | 0.288001 | extrapolated below held-out range | 7,555 | 15.9% | 8.9% | 80.9% | 22,287 | 21.2% | 17.7% | 85.3% |
| `1e-5` | 0.247951 | extrapolated below held-out range | 15,258 | 31.4% | 16.5% | 89.7% | 44,165 | 38.4% | 32.9% | 88.5% |
| `2e-5` | 0.234669 | exact Gate-2 point | 19,910 | 41.6% | 21.3% | 88.8% | 62,077 | 50.9% | 43.7% | 88.5% |
| `5e-5` | 0.216162 | exact Gate-2 point | 30,216 | 62.5% | 31.6% | 90.0% | 90,860 | 68.7% | 59.8% | 91.1% |
| `1e-4` | 0.201409 | exact Gate-2 point | 39,731 | 79.6% | 41.5% | 91.4% | 125,572 | 84.7% | 75.0% | 95.2% |
| `2e-4` | 0.185973 | exact Gate-2 point | 67,194 | 97.5% | 61.1% | 98.0% | 196,306 | 97.2% | 90.7% | 99.3% |
| **`3e-4`** | **0.176615** | **within-range interpolation** | **84,437** | **99.5%** | **71.0%** | **99.5%** | **252,471** | **99.6%** | **96.3%** | **99.9%** |
| `3.528045626e-4` | 0.172803 | PR66-cutoff interpolation | 92,944 | 99.8% | 75.2% | 99.8% | 281,490 | 99.9% | 97.8% | 100.0% |

The full metric panel is [here](analysis/plots/core_metrics_vs_log10_p.png),
with its [SVG version](analysis/plots/core_metrics_vs_log10_p.svg). The source
tables are [network_summary.tsv](analysis/network_summary.tsv) and
[operating_point_screen.tsv](analysis/operating_point_screen.tsv).

## Why `p=3e-4` was nominated

The endpoint-informed engineering screen required both driver classes to have:

- at least 90% active candidate drivers;
- at least 70% of expression nodes incident to an edge; and
- at least 95% of incident nodes in the largest weak component.

The loosest exact Gate-2 point, `p=2e-4`, passed active-driver and LCC floors in
both networks. It failed only TF incident coverage: 61.05% versus the 70% floor.
The next point, `p=3e-4`, passed all three floors in both networks, so the
selection rule returned it as a provisional interpolated fallback.

This boundary is not robust enough to call final. TF incident coverage at
`p=3e-4` is 70.995%, a margin of only 0.995 percentage points. A denser grid
between `2e-4` and `3e-4` is warranted before proposing a production default.

## Relationship to PR66

The PR67 probability `0.000352804562601613` gives the same MI cutoff as PR66's
legacy affine threshold at `m=80`: `0.1728032151574967`. At that point, the TF
and SIG 3-column networks, enhanced networks, and support tables are
byte-identical to PR66.

At the nominated `p=3e-4` point:

- all retained edges are PR66 edges;
- 90.85% of PR66 TF edges and 89.69% of PR66 SIG edges are retained;
- median zero-filled targets per candidate driver are 11 for both TF and SIG,
  versus 13 at the exact PR66-cutoff match;
- median support fractions are 0.14 for TF and 0.13 for SIG; and
- median consensus MI values are 0.2175 for TF and 0.2168 for SIG.

Across adjacent thresholds, every stricter consensus edge was retained at the
next looser threshold in this experiment. From `2e-4` to `3e-4`, directed-edge
Jaccard was 0.796 for TF and 0.778 for SIG. From `3e-4` to the PR66 cutoff match,
it was 0.908 and 0.897, respectively. These are observed BRCA100 nesting
results, not a general guarantee of the consensus procedure.

The edge-overlap curve is [here](analysis/plots/edge_overlap_vs_log10_p.png),
with its [SVG version](analysis/plots/edge_overlap_vs_log10_p.svg). Exact
comparisons are in [adjacent_overlap.tsv](analysis/adjacent_overlap.tsv) and
[pr66_context_overlap.tsv](analysis/pr66_context_overlap.tsv).

## Null-burden proxy

After SJARACNe's same-accession and same-symbol exclusions, each seed-level TF
network tests 73,743,244 directed candidate pairs and each SIG network tests
301,984,614. Multiplying these counts by `p` gives a model-based expectation
under a global independence null before DPI:

| Per-seed p | TF nominal exceedances | SIG nominal exceedances |
|---:|---:|---:|
| `1e-7` | 7.37 | 30.20 |
| `2e-4` | 14,748.65 | 60,396.92 |
| `3e-4` | 22,122.97 | 90,595.38 |
| `3.528045626e-4` | 26,016.95 | 106,541.55 |

These numbers are **not** expected consensus false edges, post-DPI counts, or
an empirical FDR. Correlated expression, DPI, and the recurrence-based
consensus step prevent that interpretation. They make the cost of loosening a
per-pair threshold explicit and strengthen the case for later biological or
independent-cohort validation.

The coverage-versus-burden plot is
[here](analysis/plots/coverage_vs_nominal_null_burden.png), with its
[SVG version](analysis/plots/coverage_vs_nominal_null_burden.svg).

## Integrity checks

Two independent read-only audits found zero discrepancies.

- All 1,800 inference jobs were present: 9 thresholds × 2 driver classes × 100
  seeds. All inference stderr files were empty, and no partial artifacts
  remained.
- All 18 consensus, support, and NetBID2 summary runs validated. Six optional
  NetBID2 HTML reports were also validated.
- All 400 anchor comparisons passed: the `p=1e-7` data sections reproduce the
  prior PR67 run, and the cutoff-match data sections reproduce PR66.
- Raw topology, support, MI, target-size, NetBID2, overlap, and null-burden
  values were independently recomputed with zero numerical discrepancies.
- `SHA256SUMS` authenticates the compact package. Large raw artifacts are
  omitted, while [omitted_artifacts.json](omitted_artifacts.json) records their
  independently revalidated sizes and SHA-256 values.

## Interpretation and next step

The sweep answers the immediate question: PR67's current `p=1e-7` operating
point is responsible for severe BRCA100 sparsity, and relaxing only that
per-subsample threshold restores driver coverage and connectivity.

It does not yet identify a universal threshold. The defensible next step is a
focused, otherwise identical sweep between `2e-4` and `3e-4`, followed by the
deferred biological-reference or independent-cohort NetBID reproducibility
test. Until then, `p=3e-4` should be described only as a provisional BRCA100
topology operating point.
