# BRCA100 consensus-recurrence sweep: compact evidence

## Scope

This directory is the compact review package for the proposed
minimum-recurrence consensus rule. We denote consensus-edge recurrence by
\(K_{\mathrm{edge}}\) to distinguish it from proposed DPI-specific recurrence
parameters; tables and filenames abbreviate it as \(K\). The package preserves
the complete \(K_{\mathrm{edge}}=6\ldots20\)
density/topology trajectory, frozen provenance, exact sweep scripts, and
representative NetBID2 QC for the proposed \(K=6\) TF and SIG networks. It
excludes the 30 large materialized networks and redundant per-\(K\) NetBID2
tables.

The sweep did not rerun inference. It reused the same 100 frozen, post-DPI
adjacency networks from each selected arm of the BRCA100 per-subsample
threshold sweep and changed only the consensus recurrence requirement.

## Matched design

The two fixed source arms are:

- TF: \(p_b=10^{-3}\), AP-MI cutoff 0.14732247558240297;
- SIG: \(p_b=5\times10^{-4}\), AP-MI cutoff 0.1644671599536221;
- BRCA100, \(m=80\) observations sampled without replacement, 100 fixed seeds,
  \(N_{\mathrm{par}}=40\), and DPI epsilon 0;
- the same expression matrix and TF/SIG driver lists used in the 13-point
  per-subsample threshold sweep; and
- SJARACNe source commit
  `7633ebb4a0d966dbda15a4e32d0efa492fb71aeb`.

The source commit is retained because it identifies the binary and files that
actually generated the evidence. Its source tree is identical to the
`develop` reapplication in PR #71, but the replacement commit hash must not
be substituted into historical provenance.

For each ordered edge \(e\), the sweep records

\[
S_e=\sum_{i=1}^{B}X_{ei},
\]

where \(X_{ei}=1\) when \(e\) appears in post-DPI network \(i\) and \(B=100\).
The materialized networks use the inclusive rule
\(S_e\ge K_{\mathrm{edge}}\) for every integer
\(K_{\mathrm{edge}}=6,\ldots,20\). Repeated copies of an ordered edge within one adjacency file
count once; reverse directions remain distinct.

## Why use \(K_{\mathrm{edge}}\) directly

The historical consensus gate sets \(q_i=E_i/U\), where \(E_i\) is the number
of ordered edges in run \(i\) and \(U\) is the observed union. It treats a
candidate edge as an independent Bernoulli event with probability \(q_i\) in
run \(i\). Because \(q_i\) differs among runs, the exact recurrence law under
those assumptions is Poisson-binomial:

\[
S_e\sim\operatorname{PB}(q_1,\ldots,q_B).
\]

The legacy implementation instead calculates

\[
\mu=\sum_i q_i,\qquad
\sigma^2=\sum_i q_i(1-q_i)
\]

and substitutes a continuous normal upper tail for this discrete distribution,
without a continuity correction. That approximation is inaccurate in the
far-right tail used by the consensus gate. At the historical
\(p_c=10^{-5}\) boundary, which maps to \(K=9\) in both selected arms:

| Network | Legacy normal tail | Exact plug-in Poisson-binomial tail |
|---|---:|---:|
| TF | \(7.08\times10^{-6}\) | \(6.66\times10^{-4}\) |
| SIG | \(5.23\times10^{-7}\) | \(2.37\times10^{-4}\) |

The normal approximation understates the tail under the same plug-in model by
about 94-fold for TF and 452-fold for SIG. Thus, the historical
\(p_c=10^{-5}\) label does not represent a \(10^{-5}\) recurrence tail even
under its own occupancy assumptions.

Replacing the normal tail with an exact Poisson-binomial calculation would fix
only that numerical approximation. It would not create a calibrated edge
p-value: \(U\) and \(E_i\) are estimated from the same observed networks, the
tested union is selected by observed recurrence, candidate ordered edges are
treated as exchangeable within a run, and overlapping resamples do not provide
independent biological evidence. The plug-in tail is not an FDR, FWER,
posterior probability, or probability that an edge is biologically correct.

The direct rule \(S_e\ge K_{\mathrm{edge}}\) records exactly what the software enforces and
avoids attaching an unsupported probability interpretation. It is
deterministic and auditable, but it is still an engineering stability rule.
Downstream analysis consumes the retained network, not either tail probability;
direct-\(K\) mode does not calculate or emit that probability.
Because \(K_{\mathrm{edge}}\) is an absolute count, both
\(K_{\mathrm{edge}}\) and \(K_{\mathrm{edge}}/B\) must be recorded and
the choice must be reconsidered when \(B\) changes.

## Density and topology evidence

Values below are TF / SIG. Median targets are zero-filled over all candidate
drivers. NetBID2 adjusted \(R^2\) is descriptive scale-free-fit QC, not an
optimizer.

| \(K\) | Consensus edges | Median targets | Drivers with at least 20 targets | NetBID2 adjusted \(R^2\) |
|---:|---:|---:|---:|---:|
| 6 | 416,408 / 739,958 | 110.5 / 47 | 100.0% / 93.8% | 0.879 / 0.864 |
| 8 | 269,294 / 462,099 | 65 / 27 | 99.7% / 64.1% | 0.887 / 0.914 |
| 9 | 224,608 / 379,053 | 51 / 21 | 96.7% / 53.4% | 0.884 / 0.925 |

The [full trajectory](results_2026-08-20/analysis/network_summary.tsv),
[target-coverage table](results_2026-08-20/analysis/driver_target_coverage.tsv),
and [density/coverage plot](results_2026-08-20/analysis/plots/recurrence_density_coverage.png)
contain all \(K=6\ldots20\) results. Loosening \(K\) adds only lower-recurrence
edges; the absolute TF/SIG cores with support at least 20, 50, and 80 remain
unchanged.

## Compact package contents

- [Results and interpretation](results_2026-08-20/RESULTS.md)
- [Complete machine-readable trajectory](results_2026-08-20/analysis/network_summary.tsv)
- [Driver target coverage](results_2026-08-20/analysis/driver_target_coverage.tsv)
- [Frozen sweep design](results_2026-08-20/provenance/design.json)
- [Original 13-point source-sweep design and run manifest](results_2026-08-20/provenance/source_sweep/)
- [Build manifest](results_2026-08-20/provenance/build_manifest.json)
- [Frozen aggregate and exact-tail records](results_2026-08-20/provenance/aggregate/)
- [Source `K=9` reproduction anchors](results_2026-08-20/provenance/k9_anchor/)
- [Recurrence-side `K=9` network and NetBID2 manifests](results_2026-08-20/provenance/arms/)
- [Hashes and sizes for omitted large artifacts](results_2026-08-20/omitted_artifacts.json)
- Exact aggregation, packaging, QC-generation scripts, and focused tests
- [Compact inventory](PACKAGE_MANIFEST.json) and [SHA-256 checksums](SHA256SUMS)

Regenerate and immediately verify the compact inventory with
`python build_compact_manifest.py`; verify an existing package without changing
it with `python build_compact_manifest.py --verify`.

The compact manifest deliberately excludes redundant per-\(K\) NetBID2 tables
and logs other than the small \(K=9\) hash-chain records, the materialized
networks, the \(K=8\) HTML reports, and generated caches. The full source
directory remains available locally but is not proposed for the review diff.

## Representative NetBID2 QC

The two included reports correspond exactly to the proposed rows below:

- [TF: \(p_b=10^{-3}\), AP-MI 0.147322, \(K=6\)](representative_netbid2_qc/k006_tf_netbid2_qc.html)
- [SIG: \(p_b=5\times10^{-4}\), AP-MI 0.164467, \(K=6\)](representative_netbid2_qc/k006_sig_netbid2_qc.html)

They were generated with R 4.4.3, NetBID2 2.2.0 at commit
`5defa454d600b94f5dd6d1f9f4428f99759a6821`, and igraph 2.3.3. The HTML
files embed their plots and have no local asset dependency. GitHub displays
HTML as source, so download the raw file and open it locally. MathJax is loaded
from `mathjax.rstudio.com`; equations may not typeset without network access.

## Reproduction and validation

The evidence utilities and their tests require Python 3.10 or newer; the
one-pass aggregator additionally requires a C++11 compiler. HTML regeneration
also requires the pinned R/NetBID2 environment recorded below.

The one-pass aggregator read each of the 200 source adjacency files once. All
30 materialized edge sets were nested, and the \(K=9\) TF/SIG networks and
NetBID2 summaries reproduced the previous outputs byte-for-byte. Independent
audits recomputed support counts, plug-in tails, target coverage, topology, and
90 NetBID2 output hashes.

The included scripts can reproduce the analysis when the omitted BRCA100 input
and source adjacency files are restored at paths supplied to
`run_recurrence_sweep.py`. Absolute paths in frozen manifests are historical
records, not portable dependencies. This compact repository package alone is
not a data-self-contained rerun bundle.

The frozen analysis manifest records that the sweep itself selected no
operating \(K\). That record is intentionally unchanged. The proposed
\(K=6\) setting below is a later engineering decision based on the observed
density/topology trade-off.

## Defensible conclusion

For the BRCA100 development workflow, the proposed density-favoring operating
points are:

| Network | Per-subsample \(p_b\) | AP-MI cutoff (\(m=80,\ N_{\mathrm{par}}=40\)) | Consensus |
|---|---:|---:|---:|
| TF | \(10^{-3}\) | 0.147322 | \(K=6\) |
| SIG | \(5\times10^{-4}\) | 0.164467 | \(K=6\) |

The per-subsample \(p_b\) filter is applied before DPI; \(K\) is applied
afterward across the 100 post-DPI networks. \(K=6\) raises the zero-filled
median target count to 110.5 for TF and 47 for SIG, but it admits edges present
in only 6% of resamples. At \(K=6\), the plug-in Poisson-binomial tails are
0.0327 for TF and 0.0177 for SIG, which reinforces that \(K=6\) is not a
significance threshold.

These are proposed BRCA100 development settings, not statistically or
biologically validated optima, FDR guarantees, or universal recommendations.
The two \(p_b\)-to-AP-MI mappings lie inside the independently held-out
calibration range, but selection of the operating points and \(K=6\) was
informed by BRCA100 topology and target coverage. The generic software
\(p_b\) default remains \(10^{-7}\) for compatibility; TF and SIG require
separate explicit invocations for the proposed settings. Biological-reference
recovery and downstream activity robustness remain deferred.
