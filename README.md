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
[Documentation](#documentation) · [Citation](#citation)

> This README describes the `develop` branch of SJARACNe2.
> Clone `develop` using the command below to get this implementation.
> The default engine is WDL 1.3 / Sprocket. The current profile applies native
> recurrence consensus before DPI, with a configurable tolerance offset of
> 0.025. PCHIP cutoffs and dynamic DPI remain experimental.

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
- Python **3.11** for the current workflow's local validation and the commands
  below. The project's `.python-version` still selects 3.14; the earlier
  [compatibility report](docs/python_compatibility.md) covers historical 3.11/3.14
  checks. uv can download the explicitly selected interpreter if it is missing.
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

Run the following in a Linux or WSL terminal. The committed `uv.lock` selects
the dependency versions. Select Python 3.11 explicitly to reproduce the
recorded local workflow validation, overriding `.python-version`:

```bash
git clone --branch develop https://github.com/adamdingliang/SJARACNe2.git
cd SJARACNe2

uv sync --locked --python 3.11
source .venv/bin/activate

sprocket --version
sjaracne local --help
```

`uv sync --locked` installs the Python dependencies and automatically builds
the native C++ executables. A separate `make` command is not needed for normal
installation. Missing build tools produce an installation error with setup
instructions; the installer does not run `sudo` or install system packages for
you. See the [native installation validation](docs/native_build_installation.md).

The native consensus-first workflow has local Python 3.11 validation. The
[Python compatibility report](docs/python_compatibility.md) records earlier
3.11/3.14 checks and does not establish a new 3.14 result for this implementation.
For a separate Python 3.14 environment, select `--python 3.14` explicitly;
use separate environments and do not replace one while jobs use it. The lock's
version range is not validation of every interpreter or platform. See the
[environment guide](docs/uv.md).

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

Omit `-n` and `-k` to use **100 resampled networks and minimum recurrence K=20**.
These examples allow up to three concurrent inference tasks with a 12-GiB
scheduling budget. Use them on a host with at least three available CPUs and
more than 12 GiB of RAM available to Linux/WSL; otherwise keep the smaller
`--jobs 2 --memory "8 GiB"` settings above.

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
| Mutual information (MI) cutoff | Experimental empirical PCHIP conversion; Npar=40 |
| Operation order | Raw MI resamples → native C++ K20 consensus → supporting-run mean MI → DPI |
| Data processing inequality (DPI) pruning | Dynamic tolerance after consensus; Href=479, slope=0.05, offset=0.025, cap=0.30 |
| Consensus | Retain directed edges present in at least 20 distinct resampled networks; average MI over supporting runs only |
| Local scheduling | Four task slots; 16 GiB memory budget |

These are **workflow defaults**, not the defaults of every bare native
executable invocation. Memory settings are scheduling requests, not a guarantee
that a dataset fits. Larger datasets may need higher budgets and task memory
requests; see [resources and outputs](SJARACNe/wdl/README.md#resources-outputs-and-provenance).

K=20 is a practical consensus setting, not a biological false-discovery-rate
guarantee or an optimum for every dataset. Inspect target coverage; the bundled
BRCA100 example becomes sparse at K20. For **scRNA-seq**, a smaller user-selected
K may be appropriate, but neither a single-cell default nor a dynamic K rule
has been validated yet. Override with `-k COUNT` as needed; `-k 6` restores the
previous recurrence cutoff, not the full legacy workflow. See the
[benchmark summary](https://github.com/adamdingliang/sjaracne_benchmark/pull/4).

The current process is **100 × 80% resamples → K20 candidates → supporting-run
mean MI → DPI**. The final native call finishes each raw-MI file before opening
the next, using one 1 MiB input buffer. It accumulates support and MI sums in RAM
and passes the retained graph directly to the existing DPI code, without an
intermediate candidate graph or aggregation cache. Dense accumulator memory
scales with the number of hubs multiplied by the number of expression genes;
larger panels can need higher `--consensus-memory` and `--memory` budgets. See
[the algorithm and memory requirements](docs/native_consensus_before_dpi.md#engineering-choices).
Python adds symbols and correlation/regression statistics to the final pruned
edge list.

Raw-graph hash re-verification at consensus is **off by default**. Producers
still compute and record each resample graph's SHA-256; consensus records those
hashes and checks seed, runtime, input and settings provenance without rereading
every graph solely to hash it. Add `--verify-raw-hashes` when a fresh byte-level
comparison with the producer hashes is required. Native format checks and
input/runtime/configuration, settings and final-output hash checks remain in
place. This operational option does not change scientific settings or results
for unchanged inputs. See [verification behavior](docs/native_consensus_before_dpi.md#optional-raw-graph-hash-verification).

```text
epsilon(H) = min(0.30, dpi_offset + 0.05 * max(0, log2(H / 479)))
```

`H` counts the requested matched source rows, including zero-target hubs. The
workflow default is `--dpi-offset 0.025`. Set it to zero to test the original
dynamic schedule while keeping consensus first. The offset is additive and the
cap applies after addition. Fixed DPI does not use an offset.

The previous current order remains available as `--consensus-order dpi-first`;
its omitted offset resolves to zero. The legacy profile retains its historical
sampling, strict DPI and probability consensus. See the
[native design and interface](docs/native_consensus_before_dpi.md) and
[benchmark results](https://github.com/adamdingliang/sjaracne_benchmark/blob/b195f8331496c645121fc614c930ce1f58049300/benchmarks/consensus_before_dpi/README.md).

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

For an explicit offset experiment:

```bash
sjaracne local --network-type TF --dpi-offset 0.025 \
  -e /absolute/path/input.exp -g /absolute/path/hubs.txt \
  -o results/tf_consensus_first

# Reproduce the previous current order and its original tolerance schedule:
sjaracne local --network-type TF --consensus-order dpi-first \
  -e /absolute/path/input.exp -g /absolute/path/hubs.txt \
  -o results/tf_dpi_first
```

WDL/Sprocket supports `--verify-raw-hashes` with either processing order. CWL
supports the option for consensus-first; the CLI rejects it with CWL DPI-first,
including CWL legacy runs, rather than silently ignoring it. Use WDL/Sprocket
when fresh graph-hash verification is needed with DPI-first.

Run `sjaracne local --help` for the full interface. The
[WDL guide](SJARACNe/wdl/README.md) explains the dynamic formula, standalone WDL
inputs, task memory requests and explicit `--sprocket-config` settings. Local
validation does not establish support for production LSF, Slurm or cloud runs.

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

Workflow and step names describe the scientific operations. The default CWL
entry point is [consensus_before_dpi.cwl](SJARACNe/cwl/consensus_before_dpi.cwl);
explicit DPI-first runs use
[dpi_before_consensus.cwl](SJARACNe/cwl/dpi_before_consensus.cwl).
Both engines name the per-seed stage `resample_and_estimate_mi`. In the default
order it estimates and significance-filters raw MI without DPI; the final stage
selects recurrent edges, averages supporting-run MI and applies DPI.
The [CWL workflow index](SJARACNe/cwl/README.md) explains the steps and the
replacement paths for removed import-only aliases.

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

Dynamic DPI applies the configured baseline offset and further relaxes pruning
as the reconstructed hub count increases above its reference. It addresses hub-context sensitivity as an engineering choice;
DREAM5 tests did not establish a dynamic-specific accuracy advantage. See the
linked benchmark reports for the evidence and its limitations.

## Documentation

| Resource | What it covers |
| --- | --- |
| [Workflow guide](SJARACNe/wdl/README.md) | WDL inputs, dynamic DPI formula, resources, outputs and provenance. |
| [Backward compatibility](docs/BACKWARD_COMPATIBILITY.md) | Legacy profile, local CWL and low-level native behavior. |
| [Regression guide](tests/README.md) | Test prerequisites, execution and retained evidence. |
| [Execution coverage](tests/README.md#execution-coverage) | Separate Python/C++ coverage, Codecov setup and measurement limits. |
| [Environment and deployment guide](docs/uv.md) | uv lockfile, pip alternative, native builds and offline/HPC/cloud precautions. |
| [Native installation validation](docs/native_build_installation.md) | Automatic C++ builds during source installation and missing-tool diagnostics. |
| [Native consensus validation](docs/native_consensus_validation.md) | Native correctness, build and default-workflow regression checks. |
| [Historical profile validation](docs/default_profiles_validation.md) | Earlier default-profile and engine-parity snapshots. |
| [Python compatibility report](docs/python_compatibility.md) | Python 3.11/3.14 dependencies, clean source installations and scientific-output comparisons. |
| [Benchmark index](https://github.com/adamdingliang/sjaracne_benchmark) | DREAM5, THP-1, BRCA100, LUAD and calibration analyses. |
| [Null calibration process](https://github.com/adamdingliang/sjaracne_benchmark/tree/main/benchmarks/apmi_null_m_sweep/broad_grid_10m) / [results](https://github.com/adamdingliang/sjaracne_benchmark/blob/main/benchmarks/apmi_null_m_sweep/broad_grid_10m/RESULTS.md) | PCHIP calibration design and measured accuracy. |
| [DPI benchmark summaries](https://github.com/adamdingliang/sjaracne_benchmark/blob/main/docs/benchmark_summaries/README.md) | Fixed and dynamic tolerance comparisons. |

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

The [current workflow validation](docs/native_consensus_validation.md#workflow-validation-refresh)
records the tested runtime revision, local Python/tool versions and results.
The checks include real consensus-first and DPI-first workflows, TF/SIG roles,
legacy settings, and WDL/CWL output comparisons.

The suite tests native correctness, real CLI execution and WDL/CWL parity,
without launching full dataset benchmark batches. GitHub Actions runs the
suite on Python 3.11 and 3.14 for PRs into and pushes to `develop` and `master`
when the workflow is present on the branch. The CI badge reports **develop's**
push status; the separate Python and C++ badges report its uploaded execution
coverage. These badges do not measure biological accuracy. See the
[coverage guide](tests/README.md#execution-coverage) for the execution-coverage
scope and reporting setup.

## Support and contributions

Report bugs or suggest improvements through
[Issues](https://github.com/adamdingliang/SJARACNe2/issues); proposed changes can
be submitted as [pull requests](https://github.com/adamdingliang/SJARACNe2/pulls)
against `develop`. Include the revision, command, engine version, input
dimensions and relevant logs. Remove sensitive data before sharing reports.

## Citation

Khatamian A, Paull EO, Califano A, Yu J.
[SJARACNe: a scalable software tool for gene network reverse engineering from big data](https://doi.org/10.1093/bioinformatics/bty907).
*Bioinformatics*, 35(12):2165–2166 (2019; published online 2018).

## License

See [LICENSE.md](LICENSE.md) for the software's license terms.
