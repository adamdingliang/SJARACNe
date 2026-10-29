# SJARACNe

SJARACNe is a command-line tool for reconstructing gene regulatory networks
from gene expression data. It builds transcription factor (TF) and signaling
(SIG) networks for exploring candidate regulators and their target genes.

[![CI](https://img.shields.io/github/actions/workflow/status/adamdingliang/SJARACNe/regression.yml?branch=develop&event=push&label=CI)](https://github.com/adamdingliang/SJARACNe/actions/workflows/regression.yml?query=branch%3Adevelop "Regression status for develop")

[Installation](#installation) · [Quick start](#quick-start) ·
[Inputs](#inputs) · [Configuration](#configuration) ·
[Documentation](#documentation) · [Citation](#citation)

> This README describes `feat/wdl13-scientific-parity`, not the published PyPI
> release. Its default workflow uses WDL 1.3 and Sprocket. PCHIP MI cutoffs and
> dynamic DPI tolerance are experimental, even though they are defaults here.

## Features

- Reconstruct TF and SIG networks from a prepared expression matrix and a
  user-supplied regulator list.
- Run the complete reconstruction and consensus workflow from one command.
- Export edge lists with gene annotations, correlation statistics and recorded
  run settings.
- Select a backward-compatible scientific profile or the retained local CWL
  workflow when needed.

## Installation

### Prerequisites

- Linux, or Linux through WSL on Windows. These are the validated workflow
  environments. On Windows, open **Ubuntu in WSL**, not Git Bash or an
  MSYS2/UCRT64 shell; native Windows execution is not covered here.
- Python **3.14 recommended**; Python 3.11 remains a tested compatibility
  option. uv can download the selected interpreter if it is missing.
- Git, Bash, Make and a C++ compiler for installation from source. On Ubuntu,
  the one-time setup below installs the build tools together.
- **uv 0.11.16** for the locked installation below; see the
  [environment guide](docs/uv.md) for installation and the ordinary pip alternative.
- **Sprocket 0.32.0**, installed separately using the
  [official installation guide](https://sprocket.bio/getting-started/installation)
  and available on `PATH`. The CLI currently requires this exact version.

The default local workflow does not require Docker. It runs trusted task
commands directly on the host, without container isolation.

### One-time setup on Ubuntu or Ubuntu in WSL

Run these commands **inside Ubuntu**:

```bash
sudo apt update
sudo apt install --yes build-essential git curl
curl -LsSf https://astral.sh/uv/0.11.16/install.sh | sh
```

Open a new Ubuntu terminal so `uv` is available on `PATH`. Other Linux
distributions should install Make and a C++ compiler through their own package
manager. Install Sprocket separately as described above.

### Install from source

Run the following in a Linux or WSL terminal. The committed `uv.lock` selects
the validated dependency versions; `.python-version` selects Python 3.14:

```bash
git clone --branch feat/wdl13-scientific-parity https://github.com/adamdingliang/SJARACNe.git
cd SJARACNe

uv sync --locked
source .venv/bin/activate

sprocket --version
sjaracne local --help
```

`uv sync --locked` installs the Python dependencies and automatically builds
the native C++ executables. A separate `make` command is not needed for normal
installation. Missing build tools produce an installation error with setup
instructions; the installer does not run `sudo` or install system packages for
you. See the [native installation validation](docs/native_build_installation.md).

For Python 3.11, use `uv sync --locked --python 3.11` instead. Use separate
environments when working with both versions; do not replace an environment
while jobs use it. Python 3.11 and 3.14 are tested; the lock's 3.11–3.14 range
does not establish validation of every intermediate version or platform.
See the [environment guide](docs/uv.md) and
[Python compatibility report](docs/python_compatibility.md).

uv installs NumPy, pandas and SciPy, but **does not install Sprocket**. A missing
or incompatible engine is an error; there is no automatic download or silent
fallback. See the [WDL guide](SJARACNe/wdl/README.md#local-use) for the verified
Sprocket archive checksum and standalone workflow instructions.

## Quick start

From the installed source checkout, run a reduced workflow check with the
bundled BRCA100 expression matrix and TF list:

```bash
sjaracne local --network-type TF \
  -e tests/inputs/BRCA100.exp \
  -g tests/inputs/BRCA100_TF.txt \
  -o results/brca100_tf_check \
  -n 2 -k 1 \
  --jobs 2 --memory "8 GiB"
```

- `-e`: prepared expression matrix.
- `-g`: regulator IDs to reconstruct as sources, also called hubs.
- `-o`: a **new or empty** output directory; existing results are not overwritten.
- `-n 2 -k 1`: two resampled networks, retaining edges present in at least one.
  This is a functional check, not a recommended scientific configuration. It
  still processes the full bundled matrix and TF list.
- `--jobs 2 --memory "8 GiB"`: up to two concurrent tasks within an 8-GiB
  scheduling budget. This fits two 4-GiB inference requests or the subsequent
  8-GiB consensus request; it is not a measured RAM requirement or a hard limit.

A successful run publishes network files under `results/brca100_tf_check/`
and records its settings and completion state in `sjaracne_run.json`.

### Full TF and SIG networks

Omit `-n` and `-k` to use **100 resampled networks and minimum recurrence K=6**.
These examples retain the same resource budgets:

```bash
sjaracne local --network-type TF \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_TF.txt \
  -o results/brca100_tf --jobs 2 --memory "8 GiB"

sjaracne local --network-type SIG \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_SIG.txt \
  -o results/brca100_sig --jobs 2 --memory "8 GiB"
```

For your own dataset, replace the example matrix and hub-list paths with your
prepared inputs. `TF` selects the default per-resample probability
`pb=0.001`; `SIG` selects `pb=0.0005`. It does **not** choose or modify your hub
list. The network type is required for the current profile, even with an
explicit `-pb` override.

## Inputs

### Expression matrix

Use a tab-separated text file: one gene or feature per row, with the first two
header fields exactly `isoformId` and `geneSymbol`, followed by observations.
The following shows the column layout, not a runnable input file. The actual
file must use tabs between fields:

| isoformId | geneSymbol | sample_1 | sample_2 | sample_3 |
| --- | --- | --- | --- | --- |
| GENE_001 | STAT1 | 2.1 | 3.4 | 1.8 |
| GENE_002 | MYD88 | 1.5 | 2.2 | 2.9 |

Values must already be prepared for network inference. SJARACNe does **not**
normalize, log-transform or filter genes for you. With the default PCHIP model,
the selected observation count `m` must be within **50–10,000**; the default is
`m = ceil(0.8 * N)`, where `N` is the number of input observations.

### Regulator / hub list

Supply one source ID per line, without a header:

```text
GENE_001
GENE_002
```

IDs must match the expression matrix's **first column**, not necessarily
`geneSymbol`. The same list supplies native reconstruction sources (`-s`) and
regulator annotations (`-l`). Bundled [TF](SJARACNe/config/TF_list.txt) and
[SIG](SJARACNe/config/SIG_list.txt) lists are reference annotations to match to
your matrix, not automatically selected or guaranteed-current annotations.

The workflow reports missing hubs and refuses to launch if none match. It
normalizes line endings on staged copies without modifying the original inputs.

## Configuration

### Current workflow defaults

| Setting | Default |
| --- | --- |
| Workflow / engine | WDL 1.3 / Sprocket 0.32.0 |
| Resampling | 80% of observations without replacement; 100 networks, seeds 1–100 |
| Per-resample probability `pb` | TF: `0.001`; SIG: `0.0005` |
| Mutual information (MI) cutoff | Experimental empirical PCHIP conversion; Npar=40 |
| Data processing inequality (DPI) pruning | Experimental dynamic tolerance; Href=479, slope=0.05, cap=0.30 |
| Consensus | Retain edges present in at least six distinct resampled networks (K=6) |
| Local scheduling | Four task slots; 16 GiB memory budget |

These are **workflow defaults**, not the defaults of every bare native
executable invocation. Memory settings are scheduling requests, not a guarantee
that a dataset fits. Larger datasets may need higher budgets and task memory
requests; see [resources and outputs](SJARACNe/wdl/README.md#resources-outputs-and-provenance).

### Common overrides

```bash
# Example: change pb and use strict DPI, keeping the other current defaults.
sjaracne local --network-type TF -pb 0.002 --fixed-dpi \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_TF.txt \
  -o results/brca100_tf_strict
```

| Option | Meaning |
| --- | --- |
| `-pb VALUE` | Override the per-resample MI null-tail probability, not a biological false-discovery rate. |
| `-n COUNT`, `-k COUNT` | Number of resampled networks and minimum edge recurrence; `1 <= k <= n`. |
| `--subsample-fraction 0.8` or `--subsample-size 80` | Select a fraction or exact number of observations without replacement; use only one. |
| `--fixed-dpi` or `--tolerance 0.1` | Use strict DPI (epsilon=0), or fixed positive epsilon. Larger epsilon relaxes pruning. |
| `--dpi-href`, `--dpi-slope`, `--dpi-max` | Adjust the dynamic tolerance schedule; incompatible with fixed DPI. |
| `--jobs 8`, `--memory "32 GiB"` | Change the local Sprocket task and memory budgets; `--serial` uses one task slot. |
| `-tmp PATH` | Choose a scratch parent; only a fresh child is used, without deleting existing contents. |

Run `sjaracne local --help` for the full interface. The
[WDL guide](SJARACNe/wdl/README.md) explains the dynamic formula, standalone WDL
inputs, task memory requests and explicit `--sprocket-config` settings. Local
validation does not establish support for production LSF, Slurm or cloud runs.

### Backward compatibility

Use the explicit legacy profile for the pre-sampling scientific settings,
still running through Sprocket by default:

```bash
sjaracne local --profile legacy \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_TF.txt \
  -o results/brca100_tf_legacy
```

Local CWL is optional. Install Node.js and add the extra to the prepared
environment:

```bash
uv sync --locked --extra cwl
```

Then add `--engine cwl` to a workflow command. Engine choice and scientific
profile are independent; choosing CWL alone does not restore legacy settings.
CWLEXEC and the old `lsf` command have been removed.

See [Backward compatibility](docs/BACKWARD_COMPATIBILITY.md) for the profile
comparison, historical overrides and native API behavior.

## Outputs

Final files are published directly in the chosen output directory:

| File | Contents |
| --- | --- |
| `consensus_network_ncol_.txt` | Source/target IDs and symbols, MI, Pearson and Spearman correlations, regression slope and regression p-value. |
| `consensus_network_3col_.txt` | Source ID, target ID and consensus MI. |
| `bootstrap_info_.txt` | Number of networks and edge-selection metadata; historical filename retained for compatibility. |
| `parameter_info_.txt` | Native reconstruction and consensus parameters. |
| `sjaracne_run.json` | Resolved options, runtime fingerprints, engine output location and run state. |

Per-resample adjacency files, task logs and provenance are retained in the
engine directory recorded in `sjaracne_run.json`. Correlations and regressions
use the full input matrix. The output regression p-value is **not** an MI or
consensus-edge probability.

## Scientific scope

Inferred edges are candidate associations, not verified direct or causal
regulation. The provisional TF/SIG probabilities are operating choices, not
universal optima or biological false-discovery-rate guarantees.

PCHIP converts sample count and probability to an MI cutoff using a bundled
empirical null surface. It requires Npar=40, `m=50..10000` and probabilities
`1e-7..0.05` (with `pb=1` reserved for preserving all MI candidates). Unsupported
inputs fail rather than silently falling back. Probabilities below `0.0005`
warn that they are outside the main assessment range; extreme-tail assessment
has limited Monte Carlo precision.

Dynamic DPI relaxes pruning as the reconstructed hub count increases above its
reference. It addresses hub-context sensitivity as an engineering choice;
DREAM5 tests did not establish a dynamic-specific accuracy advantage. See the
linked benchmark reports for the evidence and its limitations.

## Documentation

| Resource | What it covers |
| --- | --- |
| [Workflow guide](SJARACNe/wdl/README.md) | WDL inputs, dynamic DPI formula, resources, outputs and provenance. |
| [Backward compatibility](docs/BACKWARD_COMPATIBILITY.md) | Legacy profile, local CWL and low-level native behavior. |
| [Regression guide](tests/README.md) | Test prerequisites, execution and retained evidence. |
| [Environment and deployment guide](docs/uv.md) | uv lockfile, pip alternative, native builds and offline/HPC/cloud precautions. |
| [Native installation validation](docs/native_build_installation.md) | Automatic C++ builds during source installation and missing-tool diagnostics. |
| [Current validation report](docs/default_profiles_validation.md) | Default-profile and engine-parity checks, with validation limits. |
| [Python compatibility report](docs/python_compatibility.md) | Python 3.11/3.14 dependencies, clean source installations and scientific-output comparisons. |
| [Benchmark index](https://github.com/adamdingliang/sjaracne_benchmark) | DREAM5, THP-1, BRCA100, LUAD and calibration analyses. |
| [Null calibration process](https://github.com/adamdingliang/sjaracne_benchmark/tree/main/benchmarks/apmi_null_m_sweep/broad_grid_10m) / [results](https://github.com/adamdingliang/sjaracne_benchmark/blob/main/benchmarks/apmi_null_m_sweep/broad_grid_10m/RESULTS.md) | PCHIP calibration design and measured accuracy. |
| [DPI benchmark summaries](https://github.com/adamdingliang/sjaracne_benchmark/blob/main/docs/benchmark_summaries/README.md) | Fixed and dynamic tolerance comparisons. |

The benchmark repository is private and requires access; it is **not a runtime
dependency** of SJARACNe. Related software: [scMINER](https://jyyulab.github.io/scMINER/).

## Development and testing

Install the optional CWL dependencies, Node.js and Bats in addition to the
normal prerequisites, then run the small regression suite:

```bash
uv sync --locked --group regression --extra cwl
uv run --locked --no-sync python tests/run_regression.py --jobs 4
```

Add `--python 3.11` to the sync command for the compatibility lane. The
`regression` group adds the pinned build/test tools; the `cwl` extra adds CWL
support. Sprocket, Node.js, Bats and the compiler remain system prerequisites.

It tests native correctness, real CLI execution and WDL/CWL scientific parity,
without launching full dataset benchmark batches. GitHub Actions runs the
suite on Python 3.11 and 3.14 for PRs into and pushes to `develop` and `master` when the workflow is
present on the branch. The badge above reports **develop's** push status, not
this feature branch's validation.

## Support and contributions

Report bugs or suggest improvements through
[Issues](https://github.com/adamdingliang/SJARACNe/issues); proposed changes can
be submitted as [pull requests](https://github.com/adamdingliang/SJARACNe/pulls)
against `develop`. Include the revision, command, engine version, input
dimensions and relevant logs. Remove sensitive data before sharing reports.

## Citation

Khatamian A, Paull EO, Califano A, Yu J.
[SJARACNe: a scalable software tool for gene network reverse engineering from big data](https://doi.org/10.1093/bioinformatics/bty907).
*Bioinformatics*, 35(12):2165–2166 (2019; published online 2018).

## License

See [LICENSE.md](LICENSE.md) for the software's license terms.
