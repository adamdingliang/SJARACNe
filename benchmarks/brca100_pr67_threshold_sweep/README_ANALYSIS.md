# PR67 threshold-sweep analysis contract

`analyze_sweep.py` validates and summarizes only the PR67 per-subsample p-value
sweep. It applies a provisional topology operating-point screen, but it does
not establish biological optimality or estimate empirical FDR.

## Required evidence

The analyzer is fail-closed. Planned seeds in `point_manifest.json` are not
execution evidence. It requires all 13 sweep points, both TF and SIG arms,
and exactly seeds 1 through 100 in both `run_manifest.tsv` and the per-seed
metadata directories.

```text
<work-root>/
  sweep_design.json
  builds/pr67_7633ebb/{build_manifest.json,bin/,source/}
  inputs/{BRCA100.exp,BRCA100_TF.txt,BRCA100_SIG.txt}
  results/
    run_manifest.tsv
    support_summary_manifest.json
    netbid2_qc_manifest.json                 # optional immutable 26-arm summary aggregate
    netbid2_qc_html_manifest.json            # optional HTML aggregate; not analysis input
    <p_key>/
      point_manifest.json
      <tf|sig>/
        adjacency/TF_run_001.adj ... TF_run_100.adj
        seed_metadata/TF_run_001.json ... TF_run_100.json
        consensus_manifest.json
        support_summary_manifest.json
        consensus/
          consensus_network_ncol_.txt
          consensus_support.tsv
        netbid2_qc/                          # optional; required for every arm if any
        netbid2_qc_manifest.json             # optional
```

Before reading network metrics, the analyzer verifies:

- the exact PR67 commit, m=80/Npar=40 null-model hash, 80% sampling without
  replacement, DPI epsilon 0, consensus p=1e-5, binary/config hashes, and
  pinned BRCA100 input hashes;
- point-directory identity and agreement with the immutable sweep design;
- all 2,600 completed seed records, commands, fingerprints, adjacency hashes,
  and run-manifest rows;
- consensus and support fingerprints, manifests, aggregate manifests, and
  output hashes; and
- if NetBID2 QC exists, all 26 per-arm manifests, output inventories, the
  immutable root summary-aggregate fingerprint, and agreement between NetBID2
  and edge-derived QC. The separate HTML aggregate is presentation provenance
  and is not an input to this analysis.

## Run in the pinned environment

```bash
ROOT=/mnt/d/GitHub/SJARACNe-brca100-netbid-qc
RUN="$ROOT/benchmarks/brca100_netbid_qc/netbid2-r"
WORK="$HOME/sjaracne-benchmarks/brca100-pr67-threshold-sweep-20260818"

"$RUN" python \
  "$ROOT/benchmarks/brca100_pr67_threshold_sweep/analyze_sweep.py" \
  --work-root "$WORK"
```

The default output is `<work-root>/results/analysis`. Use
`--pr66-work-root /path/to/prior-brca100-work-root` to add PR66 horizontal
references and PR66-to-sweep overlap tables. The prior expression and driver
files must have the pinned BRCA100 hashes. PR66 is context, not ground truth.
Supplying this argument also requires the complete 400-row
`results/validation/anchor_seed_equivalence.tsv` and its manifest. The analyzer
requires the cutoff-match consensus edge set, MI, and support values to equal
PR66 exactly.
If stale optional PR66 tables exist without that argument, analysis stops
instead of silently mixing contexts.

Outputs are first rendered under `analysis.partial` and renamed only after the
analysis manifest is complete. The analyzer refuses to overwrite an existing
final or partial directory; remove the exact stale directory or choose a new
`--output-root` before rerunning.

## Outputs and screen

Machine-readable outputs include network summaries, adjacent and p=1e-7
anchor overlaps, seed/point/arm provenance, `operating_point_screen.tsv`, and
`selection.json`. The analysis manifest records per-arm hashes, script and
package versions, and hashes for every emitted output other than the manifest
itself.

`network_summary.tsv` reports target-size coverage against the complete,
zero-filled TF or SIG candidate list. For each target-count threshold
`N` in `{10, 20, 30, 50, 100}`, the fields
`target_size_ge_N_driver_count` and `target_size_ge_N_driver_fraction` give the
number and fraction of candidate drivers with at least `N` consensus targets.
The fraction denominator is always `candidate_drivers`, including candidates
with no retained targets.

Recurrence-qualified diagnostics repeat these summaries after restricting
edges to inclusive support-fraction thresholds `S` in `{20pct, 50pct, 80pct}`:

- `support_ge_S_edges` and `support_ge_S_edge_fraction` count qualifying edges
  and divide by all consensus edges in that arm;
- `support_ge_S_target_size_{mean,median,q25,q75}_zero_filled` and
  `support_ge_S_target_size_max` summarize qualifying targets per candidate;
  and
- `support_ge_S_target_size_ge_N_driver_{count,fraction}` gives candidate-driver
  coverage for every `S` and `N` combination above.

Support is the fraction of the 100 matched seed-level subnetworks containing an
edge. These recurrence-qualified fields measure stability conditional on the
edge surviving consensus; they are not biological-truth probabilities, edge
posterior probabilities, or FDR estimates. They are diagnostics rather than
selection gates.

The topology engineering floors apply to both TF and SIG:

- active-driver fraction at least 0.90;
- largest weak-component fraction among incident nodes at least 0.95; and
- incident-node fraction of the expression universe at least 0.70.

A point is eligible only inside the held-out range `[2e-5, 2e-3]`. The topology
screen chooses the smallest passing p; exact Gate-2-grid status (`2e-5`,
`5e-5`, `1e-4`, `2e-4`, `5e-4`, or `1e-3`) versus within-range interpolation
is reported but does not permit skipping a stricter passing point. If no point
passes, the selected point is null. Any selection is a provisional topology
operating point, not a biologically optimal threshold or a regulon-density
selection. These floors were declared after seeing the prior PR66 and PR67
endpoints, but before seeing results for the intermediate sweep points.
Gate 2 is a second, independent canonical rank-permutation stream under the
AP-MI independence null.  It is not a held-out BRCA cohort or biological
validation.

`candidate_pair_tests` is the exact directed test count after SJARACNe's
same-accession and same-non-placeholder-gene-symbol suppression. The reported
`candidate_pair_tests * p` quantity is a model-based independence-null
exceedance proxy per seed-level subnetwork before DPI. It is not a consensus
network expectation, a post-DPI quantity, or an empirical FDR.
Plots are deterministic PNG/SVG files versus `log10(p)` and include coverage
versus this nominal null-burden proxy.
Support distributions are conditional on edges retained by consensus.  Their
IQR bands summarize retained edges; they are not confidence intervals,
bootstrap uncertainty, or variation across independent benchmark cohorts.
