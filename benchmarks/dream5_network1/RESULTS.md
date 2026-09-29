# Initial DREAM5 Network 1 run (PR #70 code)

The opt-in runner was exercised on the PR #70 `develop` source at
`8c5f26c4fb5e05c7b280bcaa7d5a6461125a6aa3`. The Linux binary used for
this run has SHA256
`8264118ea83b8c37b16000ce0c36087e3fe03a28cb1687eba5329370abb8c33e`;
its source files and threshold configuration match this commit. The input
archive and three file checksums are pinned in `benchmark.py`; the generated
`.exp` SHA256 was
`b329115cc9355d1d949748eaa31270e94651e0f42df6974c0a380b70011ba692`.
The three source inputs, prepared `.exp`, TF hub list, and input provenance
are bundled under `data/Network1`. The 60 native `.adj` files, full challenge
archive, and local run logs are not in this PR; the runner regenerates the
network files.

All arms use 1,643 genes, 195 supplied TF hubs, 805 available expression
profiles, `Npar=40`, seeds 1–3, and nominal `p=1e-7`. Legacy uses 805 draws
with replacement; PR #70 uses 516/644/725 distinct profiles at 64/80/90%.
The final network uses DPI tolerance zero. Each value below is the median of
three **per-seed** statistics, not a pooled network. Counts and rates are
rounded independently.

| Measure | Legacy `-r 1` | 64% `-u` | 80% `-u` | 90% `-u` |
|---|---:|---:|---:|---:|
| Nominal-p pre-DPI edges | 145,762 | 22,362 | 33,501 | 40,896 |
| Fraction of nominal-p pre-DPI edges pruned | 91.67% | 68.08% | 77.65% | 81.66% |
| Final edges (all) | 11,475 | 7,012 | 7,487 | 7,500 |
| Final edges in labeled universe | 10,247 | 6,491 | 6,826 | 6,804 |
| True-positive final edges | 1,080 | 1,183 | 1,249 | 1,262 |
| Precision, labeled final edges | 10.71% | 18.35% | 18.29% | 18.55% |
| Sensitivity / recall | 26.92% | 29.49% | 31.13% | 31.46% |
| Specificity | 96.67% | 98.07% | 97.97% | 97.98% |
| AUROC, all unthresholded pre-DPI MI | 0.7424 | 0.7494 | 0.7586 | 0.7640 |
| Average precision, all unthresholded pre-DPI MI | 0.1552 | 0.1792 | 0.1862 | 0.1872 |

The 278,392 explicitly labeled pairs contain 4,012 positives and 274,380
negatives (1.44% prevalence). Unlabeled output pairs are excluded, not called
false positives. Raw-MI ranking assigns zero to omitted zero-MI labeled
pairs. AUROC gives half credit to positive-negative ties; average precision
groups tied scores. These are our defined metrics, not official DREAM5
challenge scores.

At a **common numeric MI threshold** of `0.03759353615996272`, selected
before this PR from the original 80% affine-cutoff pilot, the median final
edge counts are 11,475 / 7,850 / 7,487 / 7,267; precision is
10.71% / 16.79% / 18.29% / 18.95%; and sensitivity is
26.92% / 29.99% / 31.13% / 31.31% (legacy / 64 / 80 / 90%). Thus the
nominal-p edge-count differences cannot be read as a pure sampling effect.
The common cutoff is an exploratory control, **not** a gold-optimized cutoff
or an accepted PR #71 calibration.

In this simulated dataset, PR #70's without-replacement arms outperform the
legacy path on the reported median precision, recall, raw-MI AUROC, and raw-MI
average precision. The gains are modest in rank AUROC; the 80→90% gain is
especially small and not monotonic for every seed/metric. This benchmark does
not isolate sample count from sample identity, validate a new default
fraction, establish a biologically correct human regulatory network, or
measure consensus-network or downstream activity quality.
