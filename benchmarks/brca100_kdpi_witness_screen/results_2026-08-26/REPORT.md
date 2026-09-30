# BRCA100 SIG `K_DPI` witness screen

This 10-seed matched screen changes only the number of qualifying DPI witnesses required to remove an edge. It does not retune `p_b`, AP-MI, epsilon, panels, or samples.

## Primary full-minus-small result

`common_1335_sources` evaluates the same 1,335 source rows at every H and is the primary technical-bias readout. `all_sources` changes the evaluated source population with H.

| source group | K_DPI | paired seeds | median signed full-small gap | median absolute gap | absolute slope / 1,000 hubs | full/small ratio | normalized absolute gap vs K=1 | small-panel pruning retained | full-panel pruning retained | seeds gap reduced |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| all_sources | 1 | 10 | 0.19564 | 0.19564 | 0.020935 | 1.868 | 1 | 1 | 1 | 0/10 |
| all_sources | 2 | 10 | 0.21123 | 0.21123 | 0.022603 | 2.6452 | 1.093 | 0.57409 | 0.82018 | 0/10 |
| all_sources | 3 | 10 | 0.2152 | 0.2152 | 0.023028 | 3.5473 | 1.1136 | 0.382 | 0.72441 | 0/10 |
| all_sources | 5 | 10 | 0.21101 | 0.21101 | 0.02258 | 5.9722 | 1.087 | 0.18509 | 0.60595 | 1/10 |
| all_sources | 10 | 10 | 0.18077 | 0.18077 | 0.019345 | 49.294 | 0.90929 | 0.016579 | 0.43695 | 10/10 |
| common_1335_sources | 1 | 10 | 0.16916 | 0.16916 | 0.018102 | 1.7497 | 1 | 1 | 1 | 0/10 |
| common_1335_sources | 2 | 10 | 0.18487 | 0.18487 | 0.019783 | 2.4555 | 1.1143 | 0.57409 | 0.80611 | 0/10 |
| common_1335_sources | 3 | 10 | 0.18876 | 0.18876 | 0.020199 | 3.2346 | 1.137 | 0.382 | 0.70618 | 0/10 |
| common_1335_sources | 5 | 10 | 0.18302 | 0.18302 | 0.019585 | 5.4331 | 1.1027 | 0.18509 | 0.5806 | 1/10 |
| common_1335_sources | 10 | 10 | 0.15915 | 0.15915 | 0.01703 | 43.347 | 0.92592 | 0.016579 | 0.40895 | 9/10 |

## Pruning trajectories

The pruning fraction is the number of source-row pre-DPI edges with at least K qualifying witnesses divided by all pre-DPI edges in the stated source population.

| source group | H | K_DPI | seeds | median pruning fraction [IQR] | median pruning retained vs K=1 |
|---|---:|---:|---:|---:|---:|
| all_sources | 1,335 | 1 | 10 | 0.2264 [0.2163, 0.2297] | 1 |
| all_sources | 1,335 | 2 | 10 | 0.1309 [0.1243, 0.1364] | 0.57409 |
| all_sources | 1,335 | 3 | 10 | 0.0869 [0.0807, 0.0909] | 0.382 |
| all_sources | 1,335 | 5 | 10 | 0.0418 [0.0388, 0.0466] | 0.18509 |
| all_sources | 1,335 | 10 | 10 | 0.0037 [0.0029, 0.0055] | 0.016579 |
| all_sources | 5,340 | 1 | 10 | 0.3449 [0.3261, 0.3536] | 1 |
| all_sources | 5,340 | 2 | 10 | 0.2656 [0.2463, 0.2700] | 0.76679 |
| all_sources | 5,340 | 3 | 10 | 0.2218 [0.2035, 0.2248] | 0.63839 |
| all_sources | 5,340 | 5 | 10 | 0.1669 [0.1511, 0.1699] | 0.47933 |
| all_sources | 5,340 | 10 | 10 | 0.0931 [0.0830, 0.0991] | 0.26334 |
| all_sources | 10,680 | 1 | 10 | 0.4213 [0.4028, 0.4321] | 1 |
| all_sources | 10,680 | 2 | 10 | 0.3464 [0.3280, 0.3542] | 0.82018 |
| all_sources | 10,680 | 3 | 10 | 0.3068 [0.2890, 0.3125] | 0.72441 |
| all_sources | 10,680 | 5 | 10 | 0.2572 [0.2402, 0.2612] | 0.60595 |
| all_sources | 10,680 | 10 | 10 | 0.1863 [0.1714, 0.1910] | 0.43695 |
| common_1335_sources | 1,335 | 1 | 10 | 0.2264 [0.2163, 0.2297] | 1 |
| common_1335_sources | 1,335 | 2 | 10 | 0.1309 [0.1243, 0.1364] | 0.57409 |
| common_1335_sources | 1,335 | 3 | 10 | 0.0869 [0.0807, 0.0909] | 0.382 |
| common_1335_sources | 1,335 | 5 | 10 | 0.0418 [0.0388, 0.0466] | 0.18509 |
| common_1335_sources | 1,335 | 10 | 10 | 0.0037 [0.0029, 0.0055] | 0.016579 |
| common_1335_sources | 5,340 | 1 | 10 | 0.3372 [0.3213, 0.3410] | 1 |
| common_1335_sources | 5,340 | 2 | 10 | 0.2561 [0.2413, 0.2593] | 0.75611 |
| common_1335_sources | 5,340 | 3 | 10 | 0.2124 [0.1984, 0.2160] | 0.62385 |
| common_1335_sources | 5,340 | 5 | 10 | 0.1586 [0.1474, 0.1630] | 0.46167 |
| common_1335_sources | 5,340 | 10 | 10 | 0.0895 [0.0812, 0.0943] | 0.26163 |
| common_1335_sources | 10,680 | 1 | 10 | 0.3959 [0.3785, 0.4028] | 1 |
| common_1335_sources | 10,680 | 2 | 10 | 0.3198 [0.3045, 0.3239] | 0.80611 |
| common_1335_sources | 10,680 | 3 | 10 | 0.2802 [0.2661, 0.2850] | 0.70618 |
| common_1335_sources | 10,680 | 5 | 10 | 0.2309 [0.2179, 0.2350] | 0.5806 |
| common_1335_sources | 10,680 | 10 | 10 | 0.1652 [0.1539, 0.1665] | 0.40895 |

## Validation

- **PASS -- complete_matched_runs**: runs=30/30; metrics=300/300; errors=0
- **PASS -- exact_k1_reproduction**: sample/DPI/adjacency-data failures=0
- **PASS -- exact_sidecar_accounting**: sum(pre) or sum(witnesses_ge_1) failures=0
- **PASS -- frozen_candidate_provenance**: marker provenance failures=0
- **PASS -- nested_common_source_coverage**: nested-panel failures=0
- **PASS -- identical_common_source_pre_dpi_edge_counts**: source-row comparisons=26700/26700; mismatches=0; missing panel/seed pairs=0
- **PASS -- per_seed_pruning_monotone_with_kdpi**: source-group trajectories with violations=0
- **PASS -- small_panel_all_equals_common**: mismatched seed/K rows=0
- **PASS -- complete_paired_full_minus_small_effects**: paired rows=10/10

## Defensible interpretation

A higher K_DPI is promising only if it materially lowers the common-source paired gap and does so consistently across seeds while retaining a nontrivial fraction of K=1 pruning. A flatter trajectory caused by eliminating nearly all DPI pruning is not a useful correction. This screen cannot choose a biological default.
