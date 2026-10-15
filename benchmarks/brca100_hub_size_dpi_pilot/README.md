# BRCA100 hub-size DPI screening pilot

This benchmark asks one narrow question:

> With the complete 28,278-gene BRCA100 expression matrix and all other inference settings fixed, how does reducing the TF or SIG hub list affect DPI pruning?

It does **not** test Xenium 5k, choose a biologically optimal hub set, or isolate source-row availability from DPI annotation protection. It is the deliberately small, workflow-faithful screen requested before any deeper study.

## Fixed design

| Network | Per-seed p-value | AP-MI cutoff | Hub counts |
|---|---:|---:|---:|
| TF | `1e-3` | `0.14732247558240297` | 326, 1,304, 2,608 |
| SIG | `5e-4` | `0.1644671599536221` | 1,335, 5,340, 10,680 |

All six arms use:

- one deterministic, nested hub panel;
- rank-based expression-variance quintiles, SHA-256 ordering within each quintile, and round-robin quintile interleaving;
- the original input order when a selected panel is serialized;
- the same subset for `-s` and `-l`, matching the current workflow;
- seeds 1 through 100 and fixed 80-of-100 sampling without replacement;
- verbose sampling traces (`-v on`), requiring the same 80 original observation indices in all six arms for each seed;
- adaptive partitioning, `Npar=40`, and DPI epsilon 0;
- exactly one machine-readable stdout record per run:

  ```text
  [DPI_STATS] pre_edges=<n> pruned_edges=<n> post_edges=<n> dpi_applied=1
  ```

The principal outcome is `pruned_edges / pre_edges`. The supplementary consensus quantities use a direct benchmark-only `support_count >= 6` calculation. That aggregator is independent validation scaffolding, **not** the future production minimum-recurrence implementation.

## Mapping to plan steps 1-6

1. **Freeze implementation:** `run_pilot.py` snapshots and hashes the exact dirty `SJARACNe/` source tree when instrumentation is uncommitted. It records the base commit, Git status, tracked source patch hash, per-file source hashes, source-tree fingerprint, compiler, configuration, null model, and binary hash. A changed source tree cannot resume the same work root.
2. **Exact DPI accounting:** every completed seed must satisfy `pre = pruned + post`, and `post` must equal the validated adjacency edge count.
3. **Three-point matched pilot:** six arms, 100 identical seed values, one deterministic nested panel.
4. **Primary outputs:** pre/pruned/post edge quartiles, pruning-fraction median and IQR, provisional K>=6 edge count, zero-filled median target count, and active-hub fraction.
5. **Minimum recurrence:** seed inference is independent of consensus; the included direct K>=6 aggregator is provisional and labeled accordingly.
6. **Decision gates:** complete matched seeds, identical recorded sample indices across arms, frozen provenance, exact DPI accounting, prior full-size adjacency equivalence, six independently verified K>=6 outputs, and reproduction of the prior full-size K>=6 edge counts (416,408 TF and 739,958 SIG).

## Run from WSL

Use an external work root. Generated seed networks and logs do not belong in Git.
For this exact anchor replication, every phase must run through the existing
`netbid2-r` environment wrapper (GCC 15.3 and its pinned libraries). A generic
user may substitute another fully recorded environment, but that run cannot
claim exact equivalence to the historical anchors.

```bash
cd /mnt/d/GitHub/SJARACNe-hub-dpi
BENCH_ENV=/mnt/d/GitHub/SJARACNe-brca100-netbid-qc/benchmarks/brca100_netbid_qc/netbid2-r

"$BENCH_ENV" python3 benchmarks/brca100_hub_size_dpi_pilot/run_pilot.py \
  --repo /mnt/d/GitHub/SJARACNe-hub-dpi \
  --work-root ~/sjaracne-benchmarks/brca100-hub-size-dpi-pilot \
  --phase prepare
```

The prepare phase builds the exact source snapshot and generates all six panel and arm manifests. A command-only smoke check can then be run without inference:

```bash
"$BENCH_ENV" python3 benchmarks/brca100_hub_size_dpi_pilot/run_pilot.py \
  --repo /mnt/d/GitHub/SJARACNe-hub-dpi \
  --work-root ~/sjaracne-benchmarks/brca100-hub-size-dpi-pilot \
  --phase infer --drivers tf --hub-counts 326 \
  --seed-start 1 --seed-end 3 --dry-run
```

Run the complete matched inference with the existing operating-point sweep as its mandatory full-size anchor:

```bash
"$BENCH_ENV" python3 benchmarks/brca100_hub_size_dpi_pilot/run_pilot.py \
  --repo /mnt/d/GitHub/SJARACNe-hub-dpi \
  --work-root ~/sjaracne-benchmarks/brca100-hub-size-dpi-pilot \
  --phase infer --workers 12 \
  --anchor-root ~/sjaracne-benchmarks/brca100-pr67-threshold-sweep-20260818
```

The inference phase is resumable at the seed level. A completed marker is reused only after the command fingerprint, source/binary hashes, output hash, adjacency validation, and DPI accounting all agree.

After all 600 seed runs finish:

```bash
"$BENCH_ENV" python3 benchmarks/brca100_hub_size_dpi_pilot/run_pilot.py \
  --repo /mnt/d/GitHub/SJARACNe-hub-dpi \
  --work-root ~/sjaracne-benchmarks/brca100-hub-size-dpi-pilot \
  --phase aggregate

"$BENCH_ENV" python3 benchmarks/brca100_hub_size_dpi_pilot/analyze_pilot.py \
  --work-root ~/sjaracne-benchmarks/brca100-hub-size-dpi-pilot \
  --require-anchor
```

## Outputs

The external work root contains:

- `pilot_design.json`: frozen complete design and source/build/input provenance;
- `panels/panel_manifest.json`: deterministic membership, nesting, quintile counts, and hashes;
- `results/<arm>/seed_metadata/*.json`: command, timing, hashes, the exact sampled indices, adjacency statistics, DPI statistics, and optional anchor comparison;
- `results/run_manifest.tsv`: one compact row per completed seed;
- `results/<arm>/provisional_k6/`: direct support-count output and manifest;
- `results/dpi_summary.tsv`: six-arm numerical summary;
- `results/validation_gates.json`: machine-readable gate results;
- `results/REPORT.md`: compact human-readable result table and interpretation boundary.

## Unit tests

```bash
python3 -m unittest discover \
  -s benchmarks/brca100_hub_size_dpi_pilot \
  -p 'test_pilot.py' -v
```

The tests cover deterministic balanced/nested panels using the actual BRCA100 inputs, strict DPI-record parsing and accounting, disk-backed recurrence counting at the K=6 boundary, aggregation-output resume validation, and the complete synthetic six-arm analysis gates. Seed-level inference resume is enforced by the run harness but is not simulated by this unit-test module.

## Expansion rule

This one-panel pilot is descriptive. Expand only if its pruning-fraction trajectory warrants it:

- add the fixed-annotation control (`-s=subset`, `-l=full`) to explain mechanism;
- add two panel replicates if membership sensitivity matters;
- add smaller or 75% sizes if the curve needs resolution;
- add broader topology/NetBID2 QC only if retained target coverage becomes problematic.

BRCA100 retains the whole-transcriptome target universe. These results cannot be presented as a Xenium 5k validation.
