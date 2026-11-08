# Historical consensus-order benchmark results

These saved results compare the original SJARACNe workflow, K20 consensus controls, and consensus-before-DPI experiments. They were produced before the native C++ consensus feature in this branch. Implementation equivalence was subsequently checked separately, as documented below. No MI inference was rerun for that check; no speedup measurement is implied.

## Saved tables

- [per_dataset.tsv](per_dataset.tsv): 90 directed results: 18 DREAM settings multiplied by five methods. Includes TP, FP, FN, precision, recall, F1, target medians and graph hashes. Precision, recall and F1 are fractions, not percentages.
- [weighted_summary.tsv](weighted_summary.tsv): 20 rows covering five methods and four weighting schemes. Values are fractions. The median-retention summary is a weighted mean of within-dataset ratios.
- [bulk_metrics.tsv](bulk_metrics.tsv): 18 existing bulk QC rows: BRCA100, LUAD517 and BLCA406, with TF and SIG networks and three methods. Bulk data have no edge-level reference F1 in these experiments.
- [provenance.json](provenance.json): exact source paths and SHA-256 hashes, historical engine metadata, graph receipts and export checks. Raw graphs and expression data are not included.

Numerical strings in the tables are preserved from the saved source tables. FN is copied from the corresponding original comparison row. This export does not recalculate or round the saved metrics.

## Methods and controlled comparisons

| Table identifier | Inference and pruning | Consensus |
| --- | --- | --- |
| `legacy` | Original engine; 100 full-size bootstrap samples with replacement; MI p=1e-7, original affine calibration; fixed strict DPI epsilon=0 within each run | Original probability gate, pc=1e-5 |
| `legacy_k20` | Same saved legacy per-run DPI networks | At least 20 supporting runs instead of the original probability gate |
| `develop` | 100 unique 80% subsamples; MI p=0.001, PCHIP calibration; per-run dynamic DPI | At least 20 supporting runs |
| `new` | Same saved modern raw MI networks; K20 candidate edges, supporting-run mean MI, then DPI | K20 is applied before DPI |
| `new_relaxed` | Same candidates and weights as `new`; resolved DPI epsilon increased by 0.025, capped at 0.30 | Unchanged K20 candidate gate |

All runs use Npar=40 and seeds 1-100. Modern dynamic tolerance is `min(0.30, 0.05 * max(0, log2(H / 479)))`, where H includes every requested source, even sources with no retained targets. It resolves to zero for N1, N3 and N4 and approximately 0.0530951 for SysGenA. The relaxed trial uses epsilon=0.025 and approximately 0.0780951, respectively.

Supporting-run mean MI divides the sum by the number of supporting runs; absent runs are not assigned zero MI. The paired modern order comparison shares resamples, MI values, cutoff, K20 and the baseline tolerance schedule. It changes the joint consensus/weighting/DPI stage, so it does not isolate operation order from averaging or witness availability. Comparing with original legacy also changes sampling and MI calibration. The relaxed trial changes tolerance as well as the order relative to `develop`.

## F1 and target coverage

Median targets counts unique outgoing nonself targets over the **fixed requested source list, including zero-target sources**. Every emitted target counts, including benchmark-unlabelled pairs. The active-source median is supplementary. Thus this statistic describes regulon size, not reference-edge recall. Adding false edges can increase the median.

For a method and dataset, median-target retention is its median divided by the original legacy median. Every original legacy median is positive. The aggregate is a weighted mean of those ratios, not a pooled median or a fraction of true edges recovered.

The primary larger-network summary first averages the five distinct graphs within each SysGen sample-size group. It then weights N1, N3, N4, A100, A300 and A999 by input gene count. Equivalently, each graph has weight G divided by its group's number of graphs. This grouping is an analyst choice made after the experiments, not a pre-registered rule.

| Method | Six-group gene-weighted F1 (%) | Median-target retention vs. original legacy (%) |
| --- | ---: | ---: |
| Original legacy | 8.1718 | 100.00 |
| Legacy + K20, strict DPI | 10.0228 | 55.52 |
| Develop: DPI then K20 | 9.3851 | 56.99 |
| K20 then mean MI then DPI | 9.7085 | 32.72 |
| New order, epsilon +0.025 | 9.7353 | 40.35 |

Under this summary, the modern methods trade F1 against target coverage. The relaxed order has slightly higher F1 than develop but lower median-target retention. Legacy + K20 has the highest summary F1, but uses the legacy bootstrap and strict DPI. No scalar score combining F1 and regulon size is imposed.

### Sensitivity to weighting

| Method | Equal-network F1 (%) | Gene-weighted 18-network F1 (%) | Hub-weighted 18-network F1 (%) |
| --- | ---: | ---: | ---: |
| Original legacy | 8.7933 | 8.2813 | 8.4740 |
| Legacy + K20, strict DPI | 14.2179 | 12.0630 | 14.4106 |
| Develop | 14.0467 | 11.6942 | 14.3369 |
| New order | 13.7672 | 11.5851 | 13.7770 |
| Relaxed new order | 13.8607 | 11.6710 | 13.9216 |

Without grouping, gene weighting leaves develop slightly ahead of both new-order variants. These results do not establish a weighting-independent F1 winner among modern methods. All four modified methods improve F1 over original legacy in 17 of 18 settings and lose on N4.

### Classic DREAM target medians

| Dataset | Genes | Candidate sources | Legacy | Legacy + K20 | Develop | New order | Relaxed new order |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| DREAM5 N1 | 1643 | 195 | 63 | 44 | 53 | 30 | 35 |
| DREAM5 N3 | 4511 | 334 | 113.5 | 73 | 65 | 37.5 | 46 |
| DREAM5 N4 | 5950 | 333 | 126 | 80 | 89 | 49 | 63 |

## Reference and scope limits

- N3 and N4 have incomplete biological references. Unlabelled predictions are excluded from reference-defined FP and precision, but included in target medians. A reference-defined FP is not proof of a biologically false interaction.
- SysGenA truth is complete within the simulation. Each group contains five different 1000-gene generating networks. All 15 share a simulation framework; they are not independent biological cohorts. Sample-size groups use different graphs, so they do not isolate sample size. Expression-only inference omits the supplied genotypes.
- SysGenA uses all 1000 genes as candidate sources; these are not 1000 validated TFs. Therefore hub-count weights do not have a common biological interpretation across these datasets.
- The six-group weighting assigns N4 approximately 39.39% of the total weight. Gene count does not resolve reference incompleteness or simulation differences. Win counts and weighting summaries are descriptive, without a significance claim.
- Bulk topology QC, including adjusted scale-free R-squared, is not an accuracy measurement. Bulk relaxed-order results do not exist in this snapshot. BLCA406 is the earlier 406-patient cohort, not the separate MIBC372 run.
- Reused N1/N3 MI output retains its historical serialization precision. The newly run settings used an existing output-precision-only build with direct-DPI/replay controls. Historical validation did not independently reaggregate every raw resample.

## Native C++ implementation validation

The native implementation reproduced **all 18 saved DREAM datasets at both offset 0 and offset 0.025**: 36 successful comparisons with exact directed edge sets and binary64 numerical MI values. The checks also verified B=100, K=20, the resolved tolerance and every requested source row, including empty rows. This reused the existing 100 raw MI graphs per dataset; **zero MI inference jobs** were launched.

[Native equivalence receipt](native_equivalence.json) records the compiled native SHA-256, native source hashes, frozen benchmark manifest hashes, graph hashes, inherited metrics and source-raw receipt digests. [validate_saved_raw.py](validate_saved_raw.py) provides the reusable validation driver; its arguments select the original completed runs and a fresh output directory. Full raw graphs and local execution logs remain outside the repository.

Historical raw files were explicitly canonicalized into private temporary files by ordering source rows and adding omitted empty rows. Their numeric MI tokens were unchanged. N1/N3 original final files rounded MI to six digits; their full-precision saved baseline replays are the weight references. The initial rounded-output diagnostic and its resolution are recorded in the receipt. Input precision lost in old raw files cannot be recovered, so this is equivalence on those saved numerical inputs, not a claim about newly recomputed unrounded N1/N3 MI.

TP, FP, precision, recall and F1 are inherited from the frozen scored results **after exact graph equivalence**; this check did not independently rescore the gold standards. Native tests and real workflow integration checks are reported separately in the [feature validation document](../../docs/native_consensus_validation.md). The native binary for these checks has SHA-256 `e9fc297e3025b8aeb6d236a57290cb05edfbdd3f23e2f429b4f3ecf2d4cd9475`.

