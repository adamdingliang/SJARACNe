# SJARACNe2

SJARACNe2 is a command-line tool for reconstructing gene regulatory networks
from gene expression data. It builds transcription factor (TF) and signaling
(SIG) networks for exploring candidate regulators and their target genes.

[![Status](https://img.shields.io/badge/status-active-success.svg)](https://github.com/adamdingliang/SJARACNe2 "Active development")
[![Issues](https://img.shields.io/github/issues/adamdingliang/SJARACNe2)](https://github.com/adamdingliang/SJARACNe2/issues "Open issues")
[![Pull requests](https://img.shields.io/github/issues-pr/adamdingliang/SJARACNe2)](https://github.com/adamdingliang/SJARACNe2/pulls "Open pull requests")
[![CI](https://img.shields.io/github/actions/workflow/status/adamdingliang/SJARACNe2/regression.yml?branch=develop&event=push&label=develop%20CI)](https://github.com/adamdingliang/SJARACNe2/actions/workflows/regression.yml?query=branch%3Adevelop "Regression status for develop")
[![Python coverage](https://img.shields.io/codecov/c/github/adamdingliang/SJARACNe2/develop?flag=python&label=develop%20Python)](https://codecov.io/gh/adamdingliang/SJARACNe2 "Python runtime coverage for develop")
[![C++ core coverage](https://img.shields.io/codecov/c/github/adamdingliang/SJARACNe2/develop?flag=cpp&label=develop%20C%2B%2B)](https://codecov.io/gh/adamdingliang/SJARACNe2 "C++ core coverage for develop")

[Installation](#installation) · [Quick start](#quick-start) ·
[Inputs](#inputs) · [Configuration](#configuration) ·
[Outputs](#outputs) · [Documentation](#documentation) · [Citation](#citation)

> Install [SJARACNe2 v0.4.0](https://github.com/adamdingliang/SJARACNe2/releases/tag/v0.4.0) from the tagged source below.

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
- Python **3.11 or 3.14**, both tested in CI. The checkout prefers 3.14;
  the examples below explicitly select 3.11.
- Git, Bash, Make and a C++ compiler for installation from source. On Ubuntu,
  the one-time setup below installs the build tools together.
- **uv 0.11.16** for the locked installation below; see the
  [environment guide](docs/uv.md) for installation and the ordinary pip alternative.
- **Sprocket 0.32.0** for the default WDL engine, installed separately using the
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

The committed `uv.lock` pins dependencies. In a Linux or WSL terminal:

```bash
git clone --branch v0.4.0 --depth 1 https://github.com/adamdingliang/SJARACNe2.git
cd SJARACNe2

uv sync --locked --python 3.11
source .venv/bin/activate

sprocket --version
sjaracne local --help
```

`uv sync` installs Python dependencies and builds the native executables.
Install Sprocket separately and make it available on `PATH`. See the
[environment guide](docs/uv.md) for other Python versions and pip installation,
or [native installation](docs/native_build_installation.md) for build checks.

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

- `-e` and `-g`: expression matrix and regulator IDs.
- `-o`: a new or empty output directory.
- `-n 2 -k 1`: a short workflow check using the full bundled matrix and hub list.
- `--jobs 2 --memory "8 GiB"`: scheduling budget, not a measured RAM requirement.

A successful run publishes network files under `results/brca100_tf_check/`
and records its settings and completion state in `sjaracne_run.json`.

### Full TF and SIG networks

Omit `-n` and `-k` to use **100 resamples and K=20**. These examples request
three task slots and a 12-GiB scheduling budget; adjust them to fit the
resources available to Linux/WSL.

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
normalize, log-transform or filter genes for you. For PCHIP cutoff lookup (`pb != 1`),
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

The workflow reports missing hubs and refuses to launch if none match. The
current dynamic DPI mode requires at least two matched source IDs, including
sources whose final regulons are empty. Use fixed DPI (`--fixed-dpi`) for a
single-source workflow. Input preparation normalizes line endings on staged
copies without modifying the original files.

## Configuration

### Current workflow defaults

| Setting | Default |
| --- | --- |
| Workflow / engine | WDL 1.3 / Sprocket 0.32.0 |
| Resampling | 80% of observations without replacement; 100 networks, seeds 1–100 |
| Per-resample probability `pb` | TF: `0.001`; SIG: `0.0005` |
| Mutual information (MI) cutoff | PCHIP calibration; Npar=40 |
| Operation order | MI resampling → consensus → DPI |
| DPI pruning | Dynamic; Href=479, slope=0.05, offset=0.025, cap=0.30 |
| Consensus | K=20 distinct resamples; mean MI over supporting runs |
| Local scheduling | Four task slots; 16 GiB memory budget |

See [default parameter values](SJARACNe/wdl/README.md#default-parameter-values)
for consensus and DPI settings. The [native consensus guide](docs/native_consensus_before_dpi.md)
explains aggregation, memory requirements and hash verification.

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
| `--dpi-offset VALUE` | Add a finite offset in `[0,1)` before the dynamic cap; default 0.025 for consensus-first, zero for DPI-first. Requires dynamic DPI. |
| `--dpi-href`, `--dpi-slope`, `--dpi-max` | Adjust the dynamic tolerance schedule; incompatible with fixed DPI. Href remains 479 by default. |
| `--consensus-order consensus-first` or `dpi-first` | Select operation order. Current defaults to consensus-first; legacy requires DPI-first. |
| `--verify-raw-hashes` | Recompute every resample graph's SHA-256 before consensus and compare it with its producer receipt; default **off**. This is an operational check, not a scientific parameter. |
| `--jobs 8`, `--memory "32 GiB"` | Change the local Sprocket task and memory budgets; `--serial` uses one task slot. |
| `-tmp PATH` | Choose a scratch parent; only a fresh child is used, without deleting existing contents. |

Run `sjaracne local --help` for all options. The [workflow guide](SJARACNe/wdl/README.md)
covers tolerance formulas, verification options and task resources.

### Backward compatibility

Use the explicit legacy profile for the pre-sampling scientific settings,
still running through Sprocket by default. The profile does not change
scheduling budgets; specify resources separately:

```bash
sjaracne local --profile legacy \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_TF.txt \
  -o results/brca100_tf_legacy --jobs 3 --memory "12 GiB"
```

If retrying a failed run, choose a new output directory; failure logs are retained.

CWL is a retained alternative for running the workflow locally with `cwltool`
instead of WDL/Sprocket. Skip this setup if you use the default Sprocket engine.
To use CWL, install Node.js separately and add the optional `cwltool` dependency
to the prepared environment:

```bash
uv sync --locked --python 3.11 --extra cwl
```

Select `--engine cwl` and omit the Sprocket-only `--jobs` and memory options
from the examples above. Use `--serial` if you want CWL tasks to run one at a
time. For example:

```bash
sjaracne local --engine cwl --network-type TF \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_TF.txt \
  -o results/brca100_tf_cwl_check -n 2 -k 1
```

Engine choice and scientific profile are independent: choosing CWL alone
does not restore legacy settings. CWLEXEC and the old `lsf` command have been
removed.

See [Backward compatibility](docs/BACKWARD_COMPATIBILITY.md) for profile
settings and native options, or the [CWL guide](SJARACNe/cwl/README.md) for
workflow entry points and steps.

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
| [Default parameters](SJARACNe/wdl/README.md#default-parameter-values) | Values, target coverage and calibration limits. |
| [Workflow](SJARACNe/wdl/README.md) | WDL setup, inputs, resources and provenance. |
| [Native consensus](docs/native_consensus_before_dpi.md) | Aggregation, memory and graph-hash verification. |
| [Backward compatibility](docs/BACKWARD_COMPATIBILITY.md) | Legacy profile, CWL and native options. |
| [Regression](tests/README.md) | Dependencies and test commands. |
| [Execution coverage](tests/README.md#execution-coverage) | Python/C++ coverage and Codecov. |
| [Environments](docs/uv.md) | uv, pip and deployment options. |
| [Native installation validation](docs/native_build_installation.md) | Build checks and missing-tool handling. |
| [Native consensus validation](docs/native_consensus_validation.md) | Dated correctness and workflow checks. |
| [Historical profile validation](docs/default_profiles_validation.md) | K6 profile and engine-parity results. |
| [Python compatibility](docs/python_compatibility.md) | Dated cross-version measurements. |
| [Benchmark index](https://github.com/adamdingliang/sjaracne_benchmark) | Datasets and method comparisons. |
| [Null calibration](https://github.com/adamdingliang/sjaracne_benchmark/tree/main/benchmarks/apmi_null_m_sweep/broad_grid_10m) | Calibration design and results. |
| [DPI benchmarks](https://github.com/adamdingliang/sjaracne_benchmark/blob/main/docs/benchmark_summaries/README.md) | Fixed and dynamic DPI comparisons. |

Dataset benchmark reports, result tables, provenance and dataset replay drivers
belong in [sjaracne_benchmark](https://github.com/adamdingliang/sjaracne_benchmark).
Keep the algorithm implementation, usage documentation and small correctness
fixtures in this software repository. The benchmark repository is **not a runtime
dependency** of SJARACNe2. Related software: [scMINER](https://jyyulab.github.io/scMINER/).

## Development and testing

Install the optional CWL dependencies, Node.js and Bats in addition to the
normal prerequisites, then run the small regression suite:

```bash
uv sync --locked --python 3.11 --group regression --extra cwl
uv run --locked --no-sync --python 3.11 python tests/run_regression.py --jobs 4
```

Use the same explicit Python selection for both `uv sync` and `uv run`; this
keeps the checkout's 3.14 preference from selecting a different environment.
The `regression` group adds the pinned build/test tools; the `cwl` extra adds CWL
support. Sprocket, Node.js, Bats and the compiler remain system prerequisites.

The suite covers native correctness, real CLI workflows and WDL/CWL parity.
GitHub Actions runs it on Python 3.11 and 3.14 for PRs into and pushes to
`develop` and `master`. The badges show develop CI and execution coverage;
see the [coverage guide](tests/README.md#execution-coverage) for details.
[Validation records](docs/native_consensus_validation.md) preserve dated results.

## Support and contributions

Propose changes through [pull requests](https://github.com/adamdingliang/SJARACNe2/pulls)
against `develop`. Include the revision, command, engine version, input
dimensions and relevant logs. Remove sensitive data before sharing reports.

## Citation

Khatamian A, Paull EO, Califano A, Yu J.
[SJARACNe: a scalable software tool for gene network reverse engineering from big data](https://doi.org/10.1093/bioinformatics/bty907).
*Bioinformatics*, 35(12):2165–2166 (2019; published online 2018).

## License

See [LICENSE.md](LICENSE.md) for the software's license terms.
