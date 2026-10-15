# BRCA100 hub-size DPI pilot

Primary question: with the whole-transcriptome BRCA100 target universe and all other inference settings fixed, how does hub-list size change DPI pruning?

| Network | Hubs | Fraction | Median pre-DPI | Median pruned | Median post-DPI | Median pruned fraction (IQR) | K>=6 edges | Median targets, zero-filled | Active hubs |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SIG | 1,335 | 12.5% | 55,888 | 12,576 | 43,410 | 0.2249 (0.2124-0.2324) | 127,581 | 48.0 | 100.0% |
| SIG | 5,340 | 50.0% | 227,010 | 77,790 | 149,195 | 0.3432 (0.3281-0.3539) | 424,006 | 49.0 | 100.0% |
| SIG | 10,680 | 100.0% | 469,846 | 196,302 | 272,918 | 0.4166 (0.4027-0.4309) | 739,958 | 47.0 | 100.0% |
| TF | 326 | 12.5% | 22,337 | 2,450 | 19,856 | 0.1079 (0.1011-0.1157) | 61,621 | 115.0 | 100.0% |
| TF | 1,304 | 50.0% | 89,574 | 17,289 | 72,653 | 0.1917 (0.1798-0.2003) | 227,784 | 113.0 | 100.0% |
| TF | 2,608 | 100.0% | 177,684 | 42,996 | 134,840 | 0.2420 (0.2300-0.2529) | 416,408 | 110.5 | 100.0% |

## Validation gates

- **PASS -- complete_matched_runs**: rows=600 expected=600; each arm must contain seeds 1..100 exactly once
- **PASS -- empirically_matched_sampling**: invalid_records=0; mismatched_seeds=0; each seed must report the same 80 original indices in all six arms
- **PASS -- frozen_source_and_binary**: all 600 runs must share the design commit, source snapshot, binary, and panel manifest
- **PASS -- exact_dpi_accounting**: accounting_errors=0; fraction_errors=0
- **PASS -- full_size_anchor**: evaluated=200/200; full-size seed adjacency data must match the prior operating-point sweep
- **PASS -- independent_k6_aggregation**: six benchmark-only direct K>=6 outputs verified
- **PASS -- full_size_k6_anchor**: expected={'tf': 416408, 'sig': 739958}; observed={'sig': 739958, 'tf': 416408}; must reproduce the independent prior direct-recurrence sweep

## Interpretation boundary

This is a one-panel, workflow-faithful screen. A rising pruning fraction would justify the fixed-annotation control and replicated panels; it would not by itself identify a biologically optimal hub count or validate Xenium 5k.

The K>=6 results are independently counted benchmark outputs and must not be represented as the production minimum-recurrence implementation.
