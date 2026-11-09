# SJARACNe2

SJARACNe2 is a command-line tool for reconstructing gene regulatory networks
from gene expression data. It builds transcription factor (TF) and signaling
(SIG) networks for exploring candidate regulators and their target genes.

[![Status](https://img.shields.io/badge/status-active-success.svg)](https://github.com/adamdingliang/SJARACNe2 "Active development")
[![Issues](https://img.shields.io/github/issues/adamdingliang/SJARACNe2)](https://github.com/adamdingliang/SJARACNe2/issues "Open issues")
[![Pull requests](https://img.shields.io/github/issues-pr/adamdingliang/SJARACNe2)](https://github.com/adamdingliang/SJARACNe2/pulls "Open pull requests")
[![CI](https://img.shields.io/github/actions/workflow/status/adamdingliang/SJARACNe2/regression.yml?branch=develop&event=push&label=CI)](https://github.com/adamdingliang/SJARACNe2/actions/workflows/regression.yml?query=branch%3Adevelop "Regression status for develop")
[![Python coverage](https://img.shields.io/codecov/c/github/adamdingliang/SJARACNe2/develop?flag=python&label=Python)](https://codecov.io/gh/adamdingliang/SJARACNe2 "Python runtime coverage for develop")
[![C++ core coverage](https://img.shields.io/codecov/c/github/adamdingliang/SJARACNe2/develop?flag=cpp&label=C%2B%2B)](https://codecov.io/gh/adamdingliang/SJARACNe2 "C++ core coverage for develop")

[Installation](#installation) · [Quick start](#quick-start) ·
[Inputs](#inputs) · [Configuration](#configuration) ·
[Outputs](#outputs) · [Documentation](#documentation) · [Citation](#citation)

> This guide covers the retained **DPI-before-consensus** `master` snapshot.
> For the current native consensus-before-DPI workflow, use the
> [`develop` installation and guide](https://github.com/adamdingliang/SJARACNe2/blob/develop/README.md).

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

- Linux or **Ubuntu in WSL** on Windows. Native Windows/Git Bash/MSYS2
  workflow execution is not validated.
- Python **3.14 recommended**; Python 3.11 remains a tested compatibility
  option. uv can download the selected interpreter if it is missing.
- Git, Bash, Make and a C++ compiler for installation from source. On Ubuntu,
  the one-time setup below installs the build tools together.
- **uv 0.11.16** for the locked installation below; see the
  [environment guide](docs/uv.md) for installation and the ordinary pip alternative.
- **Sprocket 0.32.0**, installed separately using the
  [official installation guide](https://sprocket.bio/getting-started/installation)
  and available on `PATH`. The CLI currently requires this exact version.

The local workflow runs on the host without Docker or container isolation.

### One-time setup on Ubuntu or Ubuntu in WSL

Run these commands **inside Ubuntu**:

```bash
sudo apt update
sudo apt install --yes build-essential git curl
curl -LsSf https://astral.sh/uv/0.11.16/install.sh | sh
```

Open a new terminal if `uv` is not on `PATH`. Other Linux distributions need
Make and a C++ compiler from their package manager. Install Sprocket separately.

### Install from source

Run the following in a Linux or WSL terminal. The committed `uv.lock` selects
the validated dependency versions; `.python-version` selects Python 3.14:

```bash
git clone --branch master https://github.com/adamdingliang/SJARACNe2.git
cd SJARACNe2

uv sync --locked
source .venv/bin/activate

sprocket --version
sjaracne local --help
```

`uv sync --locked` installs dependencies and builds the native executables.
No separate Make step is needed. See [native installation](docs/native_build_installation.md)
for build-tool requirements and installation checks.

For Python 3.11, use `uv sync --locked --python 3.11` in a separate environment.
Both 3.11 and 3.14 are tested; other resolvable versions are not validated.
See the [environment guide](docs/uv.md) for version selection.

**Install Sprocket separately.** The [workflow guide](SJARACNe/wdl/README.md#local-use)
provides the pinned archive checksum and standalone commands.

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

Omit `-n` and `-k` to use **100 resampled networks and minimum recurrence K=20**.
These examples request three task slots and a 12-GiB scheduling budget.
Use the smaller quick-start settings on memory-limited hosts.

```bash
sjaracne local --network-type TF \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_TF.txt \
  -o results/brca100_tf --jobs 3 --memory "12 GiB"

sjaracne local --network-type SIG \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_SIG.txt \
  -o results/brca100_sig --jobs 3 --memory "12 GiB"
```

Run resource-heavy examples one at a time on memory-limited machines. Each
Sprocket invocation has its own budget, not a shared machine-wide memory limit.

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

Values must already be prepared for network inference. SJARACNe2 does **not**
normalize, log-transform or filter genes for you. PCHIP cutoff lookup (`pb != 1`)
requires a selected observation count `m` within **50–10,000**; the default is
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
| Consensus | Retain edges present in at least 20 distinct resampled networks (K=20 of 100) |
| Local scheduling | Four task slots; 16 GiB memory budget |

See [default parameter values](SJARACNe/wdl/README.md#default-parameter-values)
for interpretation and [workflow resources](SJARACNe/wdl/README.md#resources-outputs-and-provenance)
for scheduling and memory requirements.

### Common overrides

```bash
# Example: change pb and use strict DPI; keep other current scientific defaults.
sjaracne local --network-type TF -pb 0.002 --fixed-dpi \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_TF.txt \
  -o results/brca100_tf_strict --jobs 3 --memory "12 GiB"
```

| Option | Meaning |
| --- | --- |
| `-pb VALUE` | Override the per-resample MI null-tail probability, not a biological false-discovery rate. |
| `-n COUNT` | Number of resampled networks (default: 100). |
| `-k COUNT` or `--min-recurrence COUNT` | Retain edges present in at least `COUNT` distinct resampled networks. Current-profile default: **K=20**; requires `1 <= K <= n`. If `n < 20`, supply a smaller K explicitly. |
| `--subsample-fraction 0.8` or `--subsample-size 80` | Select a fraction or exact number of observations without replacement; use only one. |
| `--fixed-dpi` or `--tolerance 0.1` | Use strict DPI (epsilon=0), or fixed positive epsilon. Larger epsilon relaxes pruning. |
| `--dpi-href`, `--dpi-slope`, `--dpi-max` | Adjust the dynamic tolerance schedule; incompatible with fixed DPI. |
| `--jobs 8`, `--memory "32 GiB"` | Change the local Sprocket task and memory budgets; `--serial` uses one task slot. |
| `-tmp PATH` | Choose a scratch parent; only a fresh child is used, without deleting existing contents. |

Run `sjaracne local --help` for all options. The [workflow guide](SJARACNe/wdl/README.md)
covers standalone inputs, resources and custom Sprocket configurations.

### Backward compatibility

Use the legacy profile for historical sampling, MI, DPI and consensus settings.
It uses Sprocket by default; resources remain explicit:

```bash
sjaracne local --profile legacy \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_TF.txt \
  -o results/brca100_tf_legacy --jobs 3 --memory "12 GiB"
```

If retrying a failed run, choose a new output directory; failure logs are retained.

For the retained local CWL engine, install Node.js and the optional CWL dependency:

```bash
uv sync --locked --extra cwl
```

Select `--engine cwl`; omit Sprocket resource options and use `--serial` for
one task at a time. Engine choice does not change the scientific profile.
CWLEXEC and the old `lsf` command are removed.

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

## Documentation

| Guide | What it covers |
| --- | --- |
| [Default parameters](SJARACNe/wdl/README.md#default-parameter-values) | Default values, recurrence cutoff guidance and dynamic DPI tolerance. |
| [Workflow guide](SJARACNe/wdl/README.md) | WDL inputs, dynamic DPI formula, resources, outputs and provenance. |
| [Backward compatibility](docs/BACKWARD_COMPATIBILITY.md) | Legacy profile, local CWL and low-level native behavior. |
| [Regression guide](tests/README.md) | Test prerequisites, execution and retained evidence. |
| [Execution coverage](tests/README.md#execution-coverage) | Separate Python/C++ coverage, Codecov setup and measurement limits. |
| [Environment and deployment guide](docs/uv.md) | uv lockfile, pip alternative, native builds and offline/HPC/cloud precautions. |
| [Native installation validation](docs/native_build_installation.md) | Automatic C++ builds during source installation and missing-tool diagnostics. |
| [Historical profile validation](docs/default_profiles_validation.md) | Dated default-profile and engine-parity evidence. |
| [Python compatibility report](docs/python_compatibility.md) | Dated installation and cross-version output checks. |
| [Benchmark index](https://github.com/adamdingliang/sjaracne_benchmark) | DREAM5, THP-1, BRCA100, LUAD and calibration analyses. |
| [Null calibration process](https://github.com/adamdingliang/sjaracne_benchmark/tree/main/benchmarks/apmi_null_m_sweep/broad_grid_10m) / [results](https://github.com/adamdingliang/sjaracne_benchmark/blob/main/benchmarks/apmi_null_m_sweep/broad_grid_10m/RESULTS.md) | PCHIP calibration design and measured accuracy. |
| [DPI benchmark summaries](https://github.com/adamdingliang/sjaracne_benchmark/blob/main/docs/benchmark_summaries/README.md) | Fixed and dynamic tolerance comparisons. |

Dataset benchmarks are maintained in
[sjaracne_benchmark](https://github.com/adamdingliang/sjaracne_benchmark), separately
from runtime code and correctness fixtures. Related software: [scMINER](https://jyyulab.github.io/scMINER/).

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

The suite checks native correctness, CLI execution and WDL/CWL parity without
full dataset benchmark runs. CI tests Python 3.11 and 3.14 on both maintained
branches. Badges report **develop** CI and execution coverage; see the
[regression guide](tests/README.md) for coverage scope and evidence.

## Support and contributions

Submit changes as [pull requests](https://github.com/adamdingliang/SJARACNe2/pulls)
against `develop`. Include the revision, command, engine version, input
dimensions and relevant logs; remove sensitive data before sharing.

## Citation

Khatamian A, Paull EO, Califano A, Yu J.
[SJARACNe: a scalable software tool for gene network reverse engineering from big data](https://doi.org/10.1093/bioinformatics/bty907).
*Bioinformatics*, 35(12):2165–2166 (2019; published online 2018).

## License

See [LICENSE.md](LICENSE.md) for the software's license terms.
