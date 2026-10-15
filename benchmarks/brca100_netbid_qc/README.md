# BRCA100 NetBID2 network QC benchmark

This benchmark compares cumulative SJARACNe methods on the same BRCA100 input:

1. `12113fb`: pre-PR66 baseline, using the legacy N-out-of-N bootstrap with replacement and legacy affine MI cutoff.
2. `5809183`: PR66, using fixed 80% subsampling without replacement and the legacy affine cutoff evaluated at m=80.
3. `7633ebb`: PR67 stacked on PR66, using fixed 80% subsampling without replacement and the estimator-matched m=80, Npar=40 AP-MI null.

The comparison must keep expression data, TF/SIG lists, 100 resampling seeds,
`p=1e-7`, `Npar=40`, DPI, consensus, compiler, and workflow settings fixed.
NetBID2 QC should be run separately for TF and SIG consensus networks.

NetBID2's plots are descriptive QC, not proof that one inferred network is
biologically more correct. The benchmark should therefore also report edge,
node, and driver counts; edge Jaccard; per-driver target-size correlations;
retained-edge MI/rank correlations; target-size distributions; and consensus
support distributions.

## Isolated environment

The worktree-local environment does not modify the existing ArchR environment.
It pins R 4.4.3 through Micromamba and installs NetBID2 2.2.0 from exact commit
`5defa454d600b94f5dd6d1f9f4428f99759a6821`.

Create the environment from WSL:

```bash
ROOT=/mnt/d/GitHub/SJARACNe-brca100-netbid-qc
mkdir -p "$ROOT/.micromamba/bin" "$ROOT/.micromamba/home" "$ROOT/.micromamba/cache"
cp /mnt/d/GitHub/Xenium_benchmark/.micromamba/bin/micromamba \
  "$ROOT/.micromamba/bin/micromamba"
"$ROOT/.micromamba/bin/micromamba" create -y \
  -r "$ROOT/.micromamba" \
  -p "$ROOT/.micromamba/envs/netbid2-r44" \
  -f "$ROOT/benchmarks/brca100_netbid_qc/environment-netbid2-r44.yml"
"$ROOT/benchmarks/brca100_netbid_qc/netbid2-r" Rscript \
  "$ROOT/benchmarks/brca100_netbid_qc/install_netbid2.R"
```

For exact reproduction of the solved Linux packages, use the checked-in
explicit lock instead of re-solving the YAML:

```bash
ROOT=/mnt/d/GitHub/SJARACNe-brca100-netbid-qc
"$ROOT/.micromamba/bin/micromamba" create -y \
  -r "$ROOT/.micromamba" \
  -p "$ROOT/.micromamba/envs/netbid2-r44" \
  -f "$ROOT/benchmarks/brca100_netbid_qc/environment-netbid2-r44-linux-64.lock"
"$ROOT/benchmarks/brca100_netbid_qc/netbid2-r" Rscript \
  "$ROOT/benchmarks/brca100_netbid_qc/install_netbid2.R"
```

Verify the pinned package and required network-QC exports:

```bash
ROOT=/mnt/d/GitHub/SJARACNe-brca100-netbid-qc
"$ROOT/benchmarks/brca100_netbid_qc/netbid2-r" Rscript \
  "$ROOT/benchmarks/brca100_netbid_qc/verify_netbid2.R"
```

Exercise the actual NetBID2 importer and HTML QC renderer on its bundled TF
consensus network:

```bash
ROOT=/mnt/d/GitHub/SJARACNe-brca100-netbid-qc
"$ROOT/benchmarks/brca100_netbid_qc/netbid2-r" Rscript \
  "$ROOT/benchmarks/brca100_netbid_qc/smoke_test_netbid2_qc.R"
```

After any environment change, regenerate both pinned records:

```bash
ROOT=/mnt/d/GitHub/SJARACNe-brca100-netbid-qc
"$ROOT/benchmarks/brca100_netbid_qc/netbid2-r" python \
  "$ROOT/benchmarks/brca100_netbid_qc/refresh_environment_records.py" \
  --benchmark-repo "$ROOT"
```

The generated Linux package lock and R session information are recorded after
environment creation. Large environment, work, and QC output directories are
ignored by Git.

## Matched workflow

The resumable runner builds clean snapshots at the three exact commits, stages
LF-normalized Git inputs on WSL storage, and invokes the native executable for
seeds 1 through 100. Each completed adjacency is structurally validated and
checksummed before its `.partial` file is atomically renamed. The adjacency
directories contain no logs or manifests because the consensus implementation
opens every entry in its input directory. The Python entry points default to
`Path.home()/sjaracne-benchmarks/brca100-netbid-qc-20260817-rerun`, which is
persistent across WSL restarts; `/tmp` must not be used for this multi-hour
workflow. The commands below also pass the invoking WSL user's `$HOME`
explicitly because the isolated environment wrapper supplies its own `HOME`.

```bash
ROOT=/mnt/d/GitHub/SJARACNe-brca100-netbid-qc
WORK="$HOME/sjaracne-benchmarks/brca100-netbid-qc-20260817-rerun"
RUN="$ROOT/benchmarks/brca100_netbid_qc/netbid2-r"
mkdir -p "$WORK"

"$RUN" python "$ROOT/benchmarks/brca100_netbid_qc/run_workflows.py" \
  --phase build --work-root "$WORK"
"$RUN" python "$ROOT/benchmarks/brca100_netbid_qc/run_workflows.py" \
  --phase infer --seed-start 1 --seed-end 100 --workers 12 \
  --work-root "$WORK"
```

Consensus is intentionally serialized because the stock implementation keeps
the union of all seed-level edges in memory. The runner processes PR67 TF,
PR67 SIG, PR66 TF, PR66 SIG, baseline TF, and baseline SIG in that order, so
the largest expected union is last. Run the six arms, reconstruct the support
of every retained edge, render NetBID2 QC, and create cross-stage tables and
plots with:

```bash
"$RUN" python "$ROOT/benchmarks/brca100_netbid_qc/run_workflows.py" \
  --phase consensus --workers 1 --work-root "$WORK"
"$RUN" python "$ROOT/benchmarks/brca100_netbid_qc/run_support_summaries.py" \
  --work-root "$WORK"
"$RUN" python "$ROOT/benchmarks/brca100_netbid_qc/run_netbid_qc_all.py" \
  --work-root "$WORK"
"$RUN" python "$ROOT/benchmarks/brca100_netbid_qc/compare_networks.py" \
  --work-root "$WORK"
```

The comparison step writes machine-readable network, inference, pairwise,
driver-size, support, and directed-edge-membership tables. Its plots cover
network size, target-size and consensus-support ECDFs, paired target sizes,
common-edge MI, and all seven exact three-stage directed-edge membership
regions. PNG and SVG versions are emitted with deterministic metadata; the
compressed TSVs use a fixed gzip timestamp.

The repository-local `benchmarks/brca100_netbid_qc/outputs/checkpoints/`
directory is ignored by Git and is reserved for manual durable checkpoint
exports. The workflow does not copy its active WSL-home work tree there
automatically.

PR67 must receive its packaged `m=80`, `Npar=40` model explicitly. At the kept
default `p=1e-7`, the cutoff is a GPD-tail extrapolation below the independently
validated probability range (`p>=2e-5`); the expected warning and model
provenance are retained in every seed log/header. This benchmark therefore
describes the resulting topology and NetBID2 QC. It does not by itself prove
that the PR67 network is biologically superior or that `p=1e-7` is calibrated.

Per-seed elapsed times come from a mixed 12-process workload and are reported
only as execution diagnostics. They are not a controlled runtime benchmark.

The completed matched-workflow report, comparison tables, plots, and compact
provenance package are in
[`results_2026-08-17/`](results_2026-08-17/RESULTS.md).
