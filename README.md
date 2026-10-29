# SJARACNe

[![CI](https://img.shields.io/github/actions/workflow/status/adamdingliang/SJARACNe/regression.yml?branch=develop&event=push&label=CI)](https://github.com/adamdingliang/SJARACNe/actions/workflows/regression.yml?query=branch%3Adevelop "Regression status for develop")

SJARACNe reconstructs candidate regulatory associations using
adaptive-partitioning mutual information (AP-MI), DPI pruning and repeated
observation resampling. It is also integrated into
[scMINER](https://jyyulab.github.io/scMINER/).

## Current workflow defaults

This feature branch makes **WDL 1.3 with Sprocket** the default workflow,
with **PCHIP MI cutoffs**, **dynamic DPI tolerance**, and **K=6 consensus**.
PCHIP and dynamic DPI remain experimental scientific choices despite becoming
defaults. Promoting them does not establish universal calibration or improved
biological accuracy.

| Setting | Current profile |
| --- | --- |
| Workflow / engine | WDL 1.3 / Sprocket 0.32.0 |
| Network type | Explicit `--network-type TF` or `SIG`; never inferred from a filename |
| Per-resample probability `pb` | TF: `0.001`; SIG: `0.0005`; override with `-pb` |
| MI conversion | Bundled empirical PCHIP surface |
| Sampling | 80% without replacement, native ceiling rule |
| Resampled networks | 100, seeds 1 through 100 |
| AP-MI Npar limit | 40 |
| DPI | Dynamic epsilon; Href=479, slope=0.05, maximum=0.30 |
| Consensus | Edge appears in at least six distinct networks, inclusive |

The same hub list supplies reconstruction sources (`-s`) and regulator
annotations (`-l`) to the native program. Network type selects a probability
default; it does not select, infer or modify the supplied hub list.

Backward-compatible settings and the CWL route are retained. See
[BACKWARD_COMPATIBILITY.md](docs/BACKWARD_COMPATIBILITY.md) for the explicit
pre-sampling `legacy` profile, older develop settings and native API behavior.
Changing the engine alone does **not** select a historical scientific profile.

## Installation

Use Linux or WSL with Python, Make, a C++ compiler, Bash and the separately
installed [Sprocket 0.32.0 executable](https://github.com/stjude-rust-labs/sprocket/releases/tag/v0.32.0).
The Python runtime minimum remains 3.7.6; local validation uses Python 3.9 and
GitHub Actions uses Python 3.11. NumPy, pandas and SciPy are installed by pip.

These instructions describe the feature branch being reviewed, not an older
PyPI release or the current release branches:

```bash
git clone --branch feat/wdl13-scientific-parity https://github.com/adamdingliang/SJARACNe.git
cd SJARACNe
make -C SJARACNe
python -m pip install -e .
sprocket --version
```

Sprocket is not installed automatically by pip or downloaded when a workflow
starts. A missing or unsupported engine fails clearly; there is no silent CWL
fallback. The [WDL guide](SJARACNe/wdl/README.md) contains the verified binary
checksum and standalone execution instructions. The default local backend
runs trusted task commands directly on the host, without container isolation.

For CWL compatibility or the complete regression, also install the optional
CWL dependencies and Node.js:

```bash
python -m pip install -e '.[cwl]'
```

CWLEXEC is no longer supported or required. Its legacy CLI route, configuration
templates and launchers have been removed; local CWL compatibility remains.

## Inputs

The expression matrix is tab-separated, with header columns `isoformId`,
`geneSymbol`, followed by observations. Hub IDs must match the matrix's **first
column**, not necessarily its symbol column. Bundled
[TF](SJARACNe/config/TF_list.txt) and [SIG](SJARACNe/config/SIG_list.txt) lists
are reference annotations to match to your matrix, not a guarantee of the
latest annotation release.

The workflow validates hub overlap and stages LF-normalized copies. It does
not log-transform, filter or otherwise preprocess expression values for you.
Original inputs are not modified. Use a new or empty output directory.

### Outputs

The familiar result files are published at the output root:

- `consensus_network_ncol_.txt`: source/target IDs and symbols, MI, Pearson and
  Spearman correlations, regression slope and regression p-value.
- `consensus_network_3col_.txt`: source, target and consensus MI.
- `bootstrap_info_.txt`: selection rule, number of networks and recurrence or
  legacy-consensus metadata; the filename is retained for compatibility.
- `parameter_info_.txt`: native scientific parameters and consensus settings.

Sprocket attempts, per-seed adjacency files, logs and provenance remain in the
engine directory recorded by `sjaracne_run.json`. Correlations and regressions
use the full expression matrix, not just a resampled subset. The enhanced
network's regression p-value is **not** an MI or consensus-edge probability.

## Examples

Run from the source checkout after installation. `local` remains an accepted
subcommand but now uses WDL/Sprocket unless `--engine cwl` is supplied.

### Full TF and SIG networks

These commands use the defaults in the table above. The role must still be
supplied when `-pb` is explicitly overridden.

```bash
sjaracne local --network-type TF \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_TF.txt \
  -o results/brca100_tf

sjaracne local --network-type SIG \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_SIG.txt \
  -o results/brca100_sig
```

The provisional TF/SIG probabilities are operating choices, not universal
biological thresholds. For a functional smoke test only, add `-n 2 -k 1`;
default K=6 requires at least six input networks.

### Explicit overrides

```bash
# Current profile, but strict fixed DPI instead of the dynamic default.
sjaracne local --network-type TF --fixed-dpi \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_TF.txt \
  -o results/tf_strict

# Current scientific defaults through the retained CWL engine.
sjaracne local --engine cwl --network-type SIG \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_SIG.txt \
  -o results/sig_cwl

# Pre-sampling scientific profile through the default WDL engine.
sjaracne local --profile legacy \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_TF.txt \
  -o results/tf_legacy
```

| Option | Meaning |
| --- | --- |
| `--network-type TF` / `SIG` | Required for the current profile; selects its default `pb`. |
| `-pb` | Override the per-resample MI null-tail probability. |
| `-M pchip` / `affine` | Explicit MI conversion; affine does not automatically restore other historical settings. |
| `--fixed-dpi` | Use fixed epsilon, zero unless `--tolerance` is supplied. |
| `--tolerance 0.1` | Select fixed epsilon 0.1; conflicts with an explicit `--dynamic-dpi`. |
| `--dpi-href`, `--dpi-slope`, `--dpi-max` | Tune the current dynamic schedule; not allowed with fixed DPI. |
| `--subsample-fraction 0.8` / `--subsample-size 80` | Sampling without replacement; fraction and exact size are mutually exclusive. |
| `-n`, `-k` | Number of networks and inclusive minimum recurrence, respectively. |
| `-pc` | Explicit legacy normal-approximation consensus filter; mutually exclusive with `-k`. |
| `--jobs`, `--serial` | Sprocket task budget; serial means one task slot. CWL also accepts `--serial`, but not `--jobs`. |
| `--sprocket-config` | Explicit trusted engine configuration; overrides generated local scheduling configuration. |
| `-tmp` | Optional parent for a fresh engine working directory; existing unrelated files are not removed. |

Run `sjaracne local --help` for resource options and the full interface.
Scheduler backends require explicit Sprocket configuration. LSF, Slurm and cloud
execution are not established by local parity tests; the retired `lsf` command
does not submit jobs or silently fall back to local execution.

## Scientific details and limitations

PCHIP interpolates the bundled empirical null surface at selected sample count
`m` and probability `pb`. It requires Npar=40, `m=50..10000` and computable
`pb=1e-7..0.05`, except the native preserve-all request `p=1`. Unsupported inputs
fail rather than falling back to affine. Requests below 0.0005 warn that they
are outside the main assessment range; extreme tails have limited Monte Carlo
precision. No biological FDR guarantee is claimed. See the null-model
[process](https://github.com/adamdingliang/sjaracne_benchmark/tree/main/benchmarks/apmi_null_m_sweep/broad_grid_10m)
and [results](https://github.com/adamdingliang/sjaracne_benchmark/blob/main/benchmarks/apmi_null_m_sweep/broad_grid_10m/RESULTS.md).

Dynamic DPI sets one scalar per reconstruction:

```text
epsilon(H) = min(0.30, 0.05 * max(0, log2(H / 479)))
```

H counts unique reconstructed source rows before DPI, including empty rows;
it is not the annotation-list length or emitted-edge count. Larger epsilon
relaxes pruning. At H<=479 it is zero. Dynamic DPI is an engineering choice
motivated by hub-context sensitivity, not proof of more accurate edges.
DREAM5 did not establish a dynamic-specific accuracy advantage. See the
[benchmark summaries](https://github.com/adamdingliang/sjaracne_benchmark/blob/main/docs/benchmark_summaries/README.md).
Benchmark reports are in a separate private repository and require access;
SJARACNe has no runtime dependency on that repository.

The direct native executable retains its existing low-level defaults and
special modes. The public workflow passes its current profile explicitly.
See the compatibility guide before interpreting a bare native invocation as
equivalent to a full workflow run.

## Regression and validation

```bash
python tests/run_regression.py --jobs 4
```

The complete small regression requires the compiler, Sprocket, optional CWL
dependencies, Node.js and Bats. It tests native correctness, real default CLI
execution, CWL compatibility and paired WDL/CWL scientific parity. It does not
launch BRCA100, LUAD or DREAM5 benchmark batches. GitHub Actions installs a
checksum-pinned Sprocket and runs this suite for PRs into and pushes to
`develop` and `master` when this workflow is present on the branch.

See [tests/README.md](tests/README.md), the
[current-default validation report](docs/default_profiles_validation.md), and
the [WDL guide](SJARACNe/wdl/README.md). Local validation does not establish
production-scale resource behavior, remote scheduling, caching or recovery.

## Reference

Alireza Khatamian, Evan O. Paull, Andrea Califano* & Jiyang Yu*.
SJARACNe: a scalable software tool for gene network reverse engineering from
big data. Bioinformatics (2018). *Corresponding authors.
