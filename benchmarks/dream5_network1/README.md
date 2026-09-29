# DREAM5 Network 1: opt-in PR #70 benchmark

This measures direct, directed TF-to-target edge recovery on DREAM5 Network 1,
an **in-silico** network with a known simulation graph. It compares the PR #70
binary's retained legacy sampling path (`-r 1`, full-size with replacement)
with fixed-size sampling without replacement (`-u 64%`, `80%`, `90%`). It does
not test PR #71's estimator-matched null, consensus recurrence, biological
validity, or downstream NetBID activity. It is an external benchmark, **not**
a unit-test fixture or a new default-parameter recommendation.

The [DREAM5 archive](https://zenodo.org/records/17854236) is not committed.
Download `1_Challenge_Data_Supplement.zip` from that record (or use an already
extracted `Network1` directory). The script validates the published archive
MD5 and three input-file SHA256 digests and reads exact ZIP members without
extracting arbitrary archive paths. Keep data and generated networks out of
Git; `_work/` here is ignored. The archive's license is listed as CC BY-ND 4.0
at Zenodo; consult its terms before redistribution.

From the repository root, on a Linux environment with Python 3 and a built
PR #70 `sjaracne.exe` (for example `make -C SJARACNe`):

```sh
python3 benchmarks/dream5_network1/benchmark.py prepare \
  --archive /path/to/1_Challenge_Data_Supplement.zip \
  --out benchmarks/dream5_network1/_work
python3 benchmarks/dream5_network1/benchmark.py run \
  --out benchmarks/dream5_network1/_work \
  --binary SJARACNe/bin/sjaracne.exe --config SJARACNe/config
python3 benchmarks/dream5_network1/benchmark.py score \
  --archive /path/to/1_Challenge_Data_Supplement.zip \
  --out benchmarks/dream5_network1/_work
python3 -m unittest discover -s benchmarks/dream5_network1 -p 'test_*.py'
```

`--data-root /path/to/Network1` can replace `--archive` for `prepare` and
`score`. Defaults are three seeds (`1 2 3`), fractions (`64 80 90`),
`Npar=40`, nominal `p=1e-7`, DPI tolerance zero for final networks, and all
195 supplied TF hubs. `--seeds` and `--fractions` can make a smaller smoke run.
The runner also makes a nominal-p pre-DPI network for a DPI-pruning count.
An optional `--fixed-mi 0.03759353615996272` adds pre-/post-DPI runs at a
common explicit MI threshold, as a cutoff-sensitivity control; this number
was the 80% arm's affine cutoff in the initial pilot and is **not** a
gold-optimized or PR #71-calibrated threshold.
Choose a new output directory for a rerun; existing adjacencies are never
overwritten. `provenance.json`, `run_manifest.json`, logs, native adjacencies,
and `metrics.json` are written below `--out`.

The expression table is transposed to SJARACNe's gene-by-sample `.exp`
format, preserving all 1,643 genes and 805 observations. The script does not
select hubs using the gold labels. The official Network 1 label file has
278,392 explicitly labeled directed pairs: 4,012 positives and 274,380
negatives. Only those pairs enter scoring. A native output edge outside this
universe is reported as unlabeled, **not** counted as a false positive.

The no-cutoff pre-DPI run (`-p 1 -t 0 -e 1`) supplies the complete MI ranking;
its absent pairs receive MI zero. Final calls use nominal `p=1e-7` followed
by DPI (`-e 0`). The report gives precision, sensitivity, specificity, and
tie-aware AUROC and average precision (AP) for both MI rankings. AP is the
stepwise precision-recall area, grouping tied scores; this custom scorer is
**not** the DREAM5 challenge's official scoring script. Its binary final-edge
AUROC is merely `(sensitivity + specificity) / 2`.

Interpret nominal-p network density cautiously: PR #70 changes sampling
size, which changes its legacy affine p-to-MI cutoff. The legacy draw has 805
draws with replacement while the 64/80/90% arms have 516/644/725 distinct
profiles. Sharing a seed does **not** match their actual observations. Raw-MI
AUROC and AP avoid the nominal threshold but still conflate which observations
were sampled and how many were sampled. Do not tune a cutoff against this gold
file and then describe performance on the same file as independent validation.

The initial three-seed, four-arm exploratory results and these limitations are
documented in the PR. They are evidence about this simulated network, not a
claim that DREAM5 Network 1 is biological ground truth for human cancers.
