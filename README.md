# SJARACNe

[![CI](https://img.shields.io/github/actions/workflow/status/adamdingliang/SJARACNe/regression.yml?branch=develop&event=push&label=CI)](https://github.com/adamdingliang/SJARACNe/actions/workflows/regression.yml?query=branch%3Adevelop "Regression status for develop")

SJARACNe is a scalable solution of ARACNe that dramatically improves the computational 
performance, especially on the memory usage to allow even researchers with modest 
computational power to generate networks from thousands of samples. The algorithm uses 
adaptive-partitioning mutual information to estimate statistical dependence
between genes and infer candidate regulatory associations.
SJARACNe is integrated into [scMINER](https://jyyulab.github.io/scMINER/),
for single-cell RNA-seq data analysis.
Check out scMINER to learn how to use SJARACNe for analyzing single-cell RNA-seq data.

## Download

```bash
git clone --branch develop https://github.com/adamdingliang/SJARACNe.git
```


## Prerequisites

* [Python>=3.7.6](https://www.python.org/downloads/)
    * [numpy>=1.20.1](https://www.scipy.org/scipylib/download.html)
    * [scipy>=1.7.3](https://www.scipy.org/install.html)
    * [pandas>=1.2.3](https://pandas.pydata.org/)
    * [cwltool>=3.0.20201117141248](https://github.com/common-workflow-language/cwltool/releases)
* [Node.js](https://nodejs.org/) for local CWL JavaScript expressions
* Make and a C++ compiler (GCC or Clang) when building from source
* [cwlexec==0.2.2](https://github.com/IBMSpectrumComputing/cwlexec/releases) only for the legacy IBM LSF integration

These Python package minimums follow `setup.py`. The automated regression
runs on Linux; on Windows, run the source-build and workflow commands in WSL.


## Create a virtual environment (recommended)

### Using conda to create a virtual environment 

The recommended method of setting up the required Python environment and dependencies is to use the
[conda](https://conda.io/en/latest/) dependency manager:

The example uses Python 3.11, matching the GitHub Actions test environment:

```bash
conda create -n sjaracne python=3.11
conda activate sjaracne
```

## Installation

Depends on the runtime environment, [node.js](https://nodejs.org/en/download/) may be installed manually to run 
cwltool locally; [cwlexec](https://github.com/yuch7/cwlexec) may be installed manually to run on IBM LSF platform. 

There are two options to install SJARACNe and its dependencies:

### (Option 1) Install via pip

```bash
python -m pip install SJARACNe
```

### (Option 2) Install from source

```bash
git clone --branch develop https://github.com/adamdingliang/SJARACNe.git
cd SJARACNe
python setup.py build     # build SJARACNe binary
python setup.py install
```

The examples below use this source installation of `develop`; the PyPI release
may not include the latest subsampling, PCHIP, minimum-recurrence or dynamic-DPI
options.

### Install optional packages depends on runtime platform

SJARACNe workflow is implemented in [Common Workflow Language](https://www.commonwl.org/). 
Install node.js for running locally using cwltool; install cwlexec to run on IBM LSF platform. 
Users may check [Common Workflow Language](https://www.commonwl.org/) site for available workflow engines 
to run on other platforms, e.g., [Toil](https://toil.readthedocs.io/en/latest/).


## Usage

```text
$ sjaracne
usage: sjaracne [-h] {local,lsf} ...

SJARACNe is a scalable tool for gene network reverse engineering.

optional arguments:
  -h, --help   show this help message and exit

Subcommands:
  {local,lsf}  platforms
    local      run cwltool in a local workstation
    lsf        run cwlexec as in a IBM LSf cluster
```
```sjaracne``` workflow is implemented with [CWL](https://www.commonwl.org/). It supports multiple
 computing platforms. The current small regression tests execute real local
 workflows using [cwltool](https://github.com/common-workflow-language/cwltool).
 The legacy IBM LSF integration uses [cwlexec](https://github.com/IBMSpectrumComputing/cwlexec);
 current LSF option forwarding is tested with mocks, not by submitting jobs
 to a cluster.
 For the convenience, a python wrapper is developed for you to choose computing platform using ```subcommand```.
 
The local mode (sjaracne local) runs in parallel by default using cwltool's --parallel option. To run it in serial, 
use --serial option.

To use LSF mode, editing the LSF-specific configuration file SJARACNe/config/config_cwlexec.json to change the default 
queue and adjust memory reservation for each step is necessary. Consider increasing memory reservation for the resampling
step and consensus step if the dimension of your expression matrix file is large.


### Inputs
The main input is a tab-separated gene/protein-by-sample/cell expression matrix
whose first two columns contain ID and symbol. The second required input is
a list of hub/source IDs matching the expression matrix's **first column**.
Gene symbols work as hub IDs only when they are used in that first column.
Bundled reference lists are available in [TF_list.txt](SJARACNe/config/TF_list.txt)
and [SIG_list.txt](SJARACNe/config/SIG_list.txt); map them to the matrix's ID type
before use. They are not a guarantee of the latest annotation release.
An output directory is required for storing output files. Additional parameters (e.g., LSF queue) for running on different platforms are required. 
Those are available in the helping information of the corresponding subcommands, e.g., ```sjaracne lsf -h```.


### Outputs
The main output of SJARACNe is a network file, which is a tab delimited text file with the following columns: source,
target, mutual information, Pearson and Spearman correlations coefficients, regression line slope and p-value. SJARACNe
also outputs two meta information files: parameter_info_.txt and bootstrap_info_.txt, which store SJARACNe
input parameters and resampling metadata, respectively. `bootstrap_info_.txt` is retained as a legacy filename.
The network's regression p-value is not a consensus-edge probability;
recurrence selection is recorded in the metadata.

### Experimental PCHIP MI cutoffs (opt-in)

Use `-M pchip` / `--apmi-null-model pchip` to replace the legacy affine
p-to-MI conversion **for that run** with the bundled empirical PCHIP surface.
No `-M` means legacy affine, with its unchanged workflow `-pb` default of
`1e-7`. **With `-M pchip`, the workflow requires an explicit `-pb`; it does
not inherit that default or choose a TF/SIG value.** Direct CWL jobs must
likewise supply `p_value_bootstrap`. Native `-p 1` still preserves all MI.
PCHIP requires `Npar=40`, an actual selected sample count `m=50..10000`,
and a probability in the computable range `1e-7..0.05`.

See the [TF and SIG workflow examples](#full-local-runs-with-pchip) below.

The same callable is available without running inference or importing SciPy.
Here `m` is the selected sample count and `p` is the requested MI null-tail probability:

```python
from SJARACNe import get_mi_cutoff

tf_cutoff = get_mi_cutoff(m=80, p=0.001)       # approximately 0.147831
sig_cutoff = get_mi_cutoff(m=80, p=0.0005)     # approximately 0.164452
```

For one native inference run, without the workflow's repeated resampling or
consensus steps:

```bash
./SJARACNe/bin/sjaracne.exe -i ./tests/inputs/BRCA100.exp \
  -s ./tests/inputs/BRCA100_TF.txt -l ./tests/inputs/BRCA100_TF.txt \
  -N 40 -e 0 -u 80 -M pchip -p 0.001 -o brca100_tf_run.adj
```

Native flags differ from workflow flags: `-i` is the expression file, `-s`
selects source genes, and `-l` supplies regulator annotations for DPI.
`-N 40` sets the AP-MI Npar limit, `-e 0` enables zero-tolerance DPI, `-u 80` selects
80 distinct samples, `-p 0.001` requests the MI null-tail probability, and `-o`
names the adjacency output file.

Native-only `-W dpi_witnesses.tsv` adds a per-source DPI witness-count diagnostic
without changing the pruning rule; it requires active DPI and at least two
available source rows. See the [historical hub-list/DPI analysis](https://github.com/adamdingliang/sjaracne_benchmark/blob/main/benchmarks/brca100_hub_size_dpi_analysis/README.md)
for its scope and provenance; those measured pruning rates predate PCHIP.

Both TF and SIG use the same function; network identity does not choose a
different default. These probabilities are provisional BRCA100 operating
points, not universal biological recommendations; they have not been
revalidated as optimal settings under PCHIP.

PCHIP remains **experimental**; null calibration does not establish uniform
accuracy or biological FDR. See the
[process](https://github.com/adamdingliang/sjaracne_benchmark/tree/main/benchmarks/apmi_null_m_sweep/broad_grid_10m)
and [results](https://github.com/adamdingliang/sjaracne_benchmark/blob/main/benchmarks/apmi_null_m_sweep/broad_grid_10m/RESULTS.md)
for calibration, runtime policy, provenance and limitations.

### Experimental dynamic DPI tolerance (opt-in)

`--dynamic-dpi` enables a hub-count-dependent tolerance in the native executable
and the Python/CWL workflow. It changes one scalar epsilon per reconstructed
network, not the MI estimator, cutoff or DPI pruning algorithm:

```text
epsilon(H) = min(dpi_max, dpi_slope * max(0, log2(H / dpi_href)))
```

| Option | Default when enabled | Meaning |
| --- | ---: | --- |
| `--dpi-href` | 479 | Positive integer anchor. At or below it, epsilon is zero (strict DPI). |
| `--dpi-slope` | 0.05 | Finite positive epsilon increase per doubling above the anchor. |
| `--dpi-max` | 0.30 | Finite cap in `[0,1)`; zero is allowed and keeps strict DPI. |

The three tuning options require `--dynamic-dpi`. Larger epsilon relaxes DPI.
`H` is counted automatically from unique reconstructed source rows before DPI,
including rows with zero edges after MI thresholding. Missing IDs and
unconstructed rows do not count; it is not the annotation-list size or the
number of nonempty rows. This initial option requires a native `-s` source
list and at least two reconstructed sources. Adjacency replay (`-j`), native
conditional reconstruction (`-c`), and single-hub or all-gene modes are not
supported in this first experimental version. These are **native** flag names:
workflow `-c` still selects a configuration directory, `-j` selects an LSF
scheduler configuration, and `-s` selects serial local execution.

Existing behavior is unchanged without opt-in: the workflow uses strict DPI,
while the bare native executable defaults to no DPI (`-e 1`). Native
`--dynamic-dpi` cannot be combined with an explicit fixed `-e`, including
`-e 0`. Direct CWL jobs use `dynamic_dpi: true` and optional `dpi_href`,
`dpi_slope`, `dpi_max`; do not also supply `tolerance`. Native adjacency headers
record the mode, actual H, parameters and resolved decimal/hex epsilon; the
workflow also retains this metadata in `parameter_info_.txt`.

PCHIP and dynamic DPI are independent opt-ins. Dynamic DPI also works with
legacy affine MI cutoffs; enabling it does not enable PCHIP or change `-pb`.

For one native BRCA100 inference, using the experimental preset:

```bash
./SJARACNe/bin/sjaracne.exe -i ./tests/inputs/BRCA100.exp \
  -s ./tests/inputs/BRCA100_TF.txt -l ./tests/inputs/BRCA100_TF.txt \
  -N 40 -u 80 -M pchip -p 0.001 --dynamic-dpi \
  --dpi-href 479 --dpi-slope 0.05 --dpi-max 0.30 \
  -o brca100_tf_dynamic.adj
```

For the full workflow, add the same four long options to a local or LSF
example below. Workflow `-e` still means the expression file, not fixed DPI
tolerance. For example, this two-resample smoke run uses dynamic DPI with
the preset and an explicit PCHIP probability:

```bash
sjaracne local \
  -e ./tests/inputs/BRCA100.exp -g ./tests/inputs/BRCA100_TF.txt \
  -n 2 -k 1 --subsample-size 80 -M pchip -pb 0.001 \
  --dynamic-dpi --dpi-href 479 --dpi-slope 0.05 --dpi-max 0.30 \
  -o ./results/tf_dynamic_smoke -tmp ./results/tmp_tf_dynamic/
```

This is an **experimental engineering option**, not a new default or a claim
of better biological edges. BRCA100/LUAD support a possible pruning/consistency
balance; strong fixed tolerance remains competitive. DREAM5 did not show a
dynamic-specific accuracy advantage and used different anchors. See the
[benchmark summaries and conditional conclusions](https://github.com/adamdingliang/sjaracne_benchmark/blob/bd329ebf1ebdfb27307a81c583fe8301b62961be/docs/benchmark_summaries/README.md).
Benchmark reports are in a separate private repository and require access;
SJARACNe has no runtime dependency on it.

## Examples

Run these commands from the repository root after installing `develop`.
Use an output directory without a trailing slash and a temporary directory
ending in `/`. Local runs are parallel by default; add `--serial` if needed.

Parameters used below:

| Parameter | Meaning |
| --- | --- |
| `-e` | Input expression matrix. |
| `-g` | Hub/source IDs matching the expression matrix's first column; use the TF or SIG list for the desired network. |
| `-n 100` | Generate 100 resampled networks (the default). |
| `--subsample-size 80` | Use 80 distinct samples per network, sampled without replacement. Alternatively, use `--subsample-fraction 0.8` for 80% (the default); do not combine both options. |
| `-M pchip` | Opt in to PCHIP MI cutoffs; omit it to retain legacy affine conversion. |
| `-d 40` / `--depth 40` | AP-MI Npar limit, passed to native `-N`; the workflow default is 40. |
| `-pb` | Requested per-resample MI null-tail probability; required explicitly with PCHIP. Without PCHIP, the default is `1e-7`. |
| `--dynamic-dpi` | Opt in to experimental hub-count-dependent DPI tolerance; omitted means strict DPI in the workflow. |
| `--dpi-href`, `--dpi-slope`, `--dpi-max` | Optional dynamic-DPI tuning; require `--dynamic-dpi`. See [defaults and definitions](#experimental-dynamic-dpi-tolerance-opt-in). |
| `-k 6` / `--min-recurrence 6` | Retain ordered edges present in at least six input networks (inclusive; the default). Must be between 1 and the number of networks. |
| `-pc` | Explicitly select the deprecated legacy normal-approximation consensus filter instead of recurrence; cannot be combined with `-k`. Separate from `-pb`. |
| `-o`, `-tmp` | Final output directory and temporary working directory, respectively; both are required. |
| `-j` | LSF mode only: path to the scheduler configuration JSON. |

### Quick local smoke test

This uses only two resampled networks and explicitly retains edges seen in at
least one: `-n 2 -k 1`, not scientific-analysis settings. Default `-k 6` requires
at least six networks. The smoke test uses the default
80% subsampling and legacy affine MI cutoffs (`-pb 1e-7`).

```bash
sjaracne local \
  -e ./tests/inputs/BRCA100.exp -g ./tests/inputs/BRCA100_TF.txt \
  -n 2 -k 1 \
  -o ./results/tf_smoke -tmp ./results/tmp_smoke/
```

### Full local runs with PCHIP

These BRCA100 examples use 100 resampled networks of 80 samples each. The TF
`-pb 0.001` and SIG `-pb 0.0005` values are exploratory choices, **not automatic
defaults or universally optimal settings**. Both use the default consensus rule,
`-k 6`: an edge must appear in at least 6 of the 100 input networks. This is an
engineering policy, not a p-value or FDR guarantee; report both K and the number
of networks. See the current PCHIP K6
[process](https://github.com/adamdingliang/sjaracne_benchmark/blob/d39103d44b03e6c7244417bea37aa9034db7ecc2/benchmarks/brca100_pchip_recurrence/README.md)
and [results with NetBID2 QC](https://github.com/adamdingliang/sjaracne_benchmark/blob/d39103d44b03e6c7244417bea37aa9034db7ecc2/benchmarks/brca100_pchip_recurrence/results_2026-10-18/RESULTS.md).

Both commands below use **strict DPI**, the unchanged workflow default. Add
`--dynamic-dpi` to opt in to the experimental preset. Optional tuning defaults
are `--dpi-href 479`, `--dpi-slope 0.05` and `--dpi-max 0.30`.

```bash
# TF network
sjaracne local \
  -e ./tests/inputs/BRCA100.exp -g ./tests/inputs/BRCA100_TF.txt \
  -n 100 --subsample-size 80 -M pchip -pb 0.001 -k 6 \
  -o ./results/tf_pchip -tmp ./results/tmp_tf_pchip/

# SIG network
sjaracne local \
  -e ./tests/inputs/BRCA100.exp -g ./tests/inputs/BRCA100_SIG.txt \
  -n 100 --subsample-size 80 -M pchip -pb 0.0005 -k 6 \
  -o ./results/sig_pchip -tmp ./results/tmp_sig_pchip/
```

The [historical 13-point MI-threshold sweep](https://github.com/adamdingliang/sjaracne_benchmark/blob/c763a2617216db71f28ef7b99c72492f35b1bd87/benchmarks/brca100_pr67_threshold_sweep/COMBINED_RESULTS.md)
and [fresh 200-network PCHIP inference report](https://github.com/adamdingliang/sjaracne_benchmark/blob/d39103d44b03e6c7244417bea37aa9034db7ecc2/benchmarks/brca100_pchip_refresh/results_2026-10-17/RESULTS.md)
document the starting probabilities and original inference provenance. Their
older consensus results are comparison context, not PCHIP + K6 outcomes.

### Running on an IBM LSF cluster

After configuring the LSF queue and resource settings, this runs the TF workflow
with the same PCHIP and strict-DPI settings as the local example. Dynamic mode
can be requested with `--dynamic-dpi`; current regression does not exercise an
actual LSF cluster:

```bash
sjaracne lsf -j ./SJARACNe/config/config_cwlexec.json \
  -e ./tests/inputs/BRCA100.exp -g ./tests/inputs/BRCA100_TF.txt \
  -n 100 --subsample-size 80 -M pchip -pb 0.001 -k 6 \
  -o ./results/tf_lsf -tmp ./results/tmp_tf_lsf/
```

### Explicit legacy consensus mode

For retrospective comparisons, replace `-k 6` with an explicit `-pc`, such as
`-pc 1e-5`. Do not pass both. This example retains PCHIP for MI cutoffs while
using the old consensus-probability filter:

```bash
sjaracne local \
  -e ./tests/inputs/BRCA100.exp -g ./tests/inputs/BRCA100_TF.txt \
  -n 100 --subsample-size 80 -M pchip -pb 0.001 -pc 1e-5 \
  -o ./results/tf_pchip_legacy_consensus -tmp ./results/tmp_tf_pchip_legacy/
```


## Regression tests

An experimental [standalone WDL 1.3 workflow](SJARACNe/wdl/README.md) is being
developed for local Sprocket execution with the same scientific steps. It does
not replace the current Python/CWL CLI; see its documentation for scope and
validation status.

Run `python tests/run_regression.py` on Linux/WSL after installing this checkout.
It builds both native programs and checks the small end-to-end workflow,
Python tests and CWL utilities—not full BRCA100/DREAM resampling benchmarks.
See [tests/README.md](tests/README.md) for prerequisites, coverage and logs.
GitHub Actions uses the same small suite for PRs into and pushes to `develop`
and `master`, once the workflow and runner are present on the branch.

## Reference
Alireza Khatamian, Evan O. Paull, Andrea Califano* & Jiyang Yu*. SJARACNe: a scalable 
software tool for gene network reverse engineering from big data. Bioinformatics (2018). *Corresponding authors.
