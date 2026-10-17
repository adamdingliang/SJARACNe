# SJARACNe

SJARACNe is a scalable solution of ARACNe that dramatically improves the computational 
performance, especially on the memory usage to allow even researchers with modest 
computational power to generate networks from thousands of samples. The algorithm uses 
adaptive partitioning mutual information to calculate the correlation between all pairs 
of genes to reconstruct the regulatory network. 
SJARACNe is now integrated into our latest framework, [scMINER](https://jyyulab.github.io/scMINER/), 
for single-cell RNA-seq data analysis.
Check out scMINER to learn how to use SJARACNe for analyzing single-cell RNA-seq data.

## Download
```git clone https://github.com/jyyulab/SJARACNe  # Clone the repo```


## Prerequisites
* [Python>=3.7.6](https://www.python.org/downloads/)
    * [numpy>=1.20.1](https://www.scipy.org/scipylib/download.html)
    * [scipy>=1.6.1](https://www.scipy.org/install.html)
    * [pandas>=1.2.3](https://pandas.pydata.org/)
    * [cwltool>=3.0.20201117141248](https://github.com/common-workflow-language/cwltool/releases)
    * [node.js>=4.4.4](https://nodejs.org/fa/blog/release/v4.4.4/) (required by cwltool to run locally)
* [cwlexec==0.2.2](https://github.com/IBMSpectrumComputing/cwlexec/releases) (CWL engine to run on IBM LSF)


## Create a virtual environment (recommended)
### Using conda to create a virtual environment 
The recommended method of setting up the required Python environment and dependencies is to use the
[conda](https://conda.io/en/latest/) dependency manager:

```bash
$ conda create -n py392 python=3.9.2
$ source activate py392
```

## Installation
Depends on the runtime environment, [node.js](https://nodejs.org/en/download/) may be installed manually to run 
cwltool locally; [cwlexec](https://github.com/yuch7/cwlexec) may be installed manually to run on IBM LSF platform. 

There are two options to install SJARACNe and its dependencies:
### (Option 1) Install via pip
```$ pip install SJARACNe```

### (Option 2) Install from source
```bash
$ git clone https://github.com/jyyulab/SJARACNe
$ cd SJARACNe
$ python setup.py build     # build SJARACNe binary
$ python setup.py install
```

### Install optional packages depends on runtime platform
SJARACNe workflow is implemented in [Common Workflow Language](https://www.commonwl.org/). 
Install node.js for running locally using cwltool; install cwlexec to run on IBM LSF platform. 
Users may check [Common Workflow Language](https://www.commonwl.org/) site for available workflow engines 
to run on other platforms, e.g., [Toil](https://toil.readthedocs.io/en/latest/).


## Usage
```$ sjaracne 
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
 computing platforms. We have tested it locally using [cwltool](https://github.com/common-workflow-language/cwltool) 
 and on an IBM LSF cluster using [cwlexec](https://github.com/IBMSpectrumComputing/cwlexec). 
 For the convenience, a python wrapper is developed for you to choose computing platform using ```subcommand```.
 
The local mode (sjaracne local) runs in parallel by default using cwltool's --parallel option. To run it in serial, 
use --serial option.

To use LSF mode, editing the LSF-specific configuration file SJARACNe/config/config_cwlexec.json to change the default 
queue and adjust memory reservation for each step is necessary. Consider increasing memory reservation for the resampling
step and consensus step if the dimension of your expression matrix file is large.


### Inputs
The main input for SJARACNe is a tab-separated genes/protein by cells/samples expression matrix
with the first two columns being ID and symbol. The second required input file is the list of
significant genes/proteins IDs to be considered as hubs in the reconstructed network (**the most recent version of curated 
transcription factors and signaling proteins can be found in ./SJARACNe/config/TF_list.txt and ./SJARACNe/config/SIG_list.txt, respectively**). 
An output directory is required for storing output files. Additional parameters (e.g., LSF queue) for running on different platforms are required. 
Those are available in the helping information of the corresponding subcommands, e.g., ```sjaracne lsf -h```.


### Outputs
The main output of SJARACNe is a network file, which is a tab delimited text file with the following columns: source,
target, mutual information, Pearson and Spearman correlations coefficients, regression line slope and p-value. SJARACNe
also outputs two meta information files: parameter_info_.txt and bootstrap_info_.txt, which store SJARACNe
input parameters and resampling metadata, respectively. `bootstrap_info_.txt` is retained as a legacy filename.

### Repeated subsampling of observations

The standard CWL/Python workflow now builds each input network from a fixed-size
subset sampled **without replacement**. By default, each network uses
`m = ceil(0.8 * N)` distinct observations, where `N` is the number of eligible
observations after any conditional selection. Thus, BRCA100 uses 80 of its 100
samples. The MI significance threshold and optional noise-correction variance are
calculated using the actual sampled `m`.

For every seed, SJARACNe draws a uniform fixed-size subset with a partial
Fisher-Yates shuffle, keeps the selected columns in their original order, and
recomputes ranks on that subset before adaptive-partitioning MI. Because an
observation can occur at most once, resampling no longer creates duplicated
joint expression points that can spuriously drive adaptive partitions.

Use `--subsample-fraction` to choose another fraction or `--subsample-size` to
set an exact count. For example:

```bash
sjaracne local -e expression.exp -g hubs.txt -n 100 \
  --subsample-fraction 0.8 -o results -tmp tmp

sjaracne local -e expression.exp -g hubs.txt -n 100 \
  --subsample-size 80 -o results -tmp tmp
```

At the native executable level, `-u 80%` and `-u 80` request the corresponding
without-replacement samples. Native runs use all observations when `-u` is
omitted. The old `-r` full-size bootstrap-with-replacement path remains available
only to reproduce legacy analyses; `-r` and `-u` cannot be combined.

In a manual CWL job, `subsample_spec` is a string, so quote both percentages and
exact counts (for example, `subsample_spec: "80%"` or `subsample_spec: "80"`).
The low-level `sjaracne.cwl` no longer supplies its historical implicit `-r 1`;
callers must choose `subsample_spec`, explicitly request legacy `sample_number`,
or omit both to analyze all observations. Legacy workflow names such as
`bootstrap_num` and the `TF_run_*.adj` filenames are retained for compatibility.

The 80% default is a pragmatic starting point, not a universal optimum. For an
important dataset, compare a small sensitivity range (for example, 64%, 80%, and
90%) and recalibrate any sample-size-dependent MI threshold for each choice.
Dataset-level benchmark processes, results and data provenance, including the
simulated DREAM5 Network 1 comparison, are maintained separately in
[sjaracne_benchmark](https://github.com/adamdingliang/sjaracne_benchmark)
(currently private). Historical implementation-performance evidence and its
reusable benchmarking tools are maintained there as well. SJARACNe does not
require that repository to build, run or test. Small engine correctness fixtures
remain here under `tests/fixtures/`; the root `benchmarks/` directory has been
retired. Benchmark results are not a new default-fraction recommendation.

### Experimental PCHIP MI cutoffs (opt-in)

Use `-M pchip` / `--apmi-null-model pchip` to replace the legacy affine
p-to-MI conversion **for that run** with the bundled empirical PCHIP surface.
Existing defaults remain unchanged: no `-M` means legacy affine, workflow
`-pb` remains `1e-7`, and native `-p 1` preserves all MI. PCHIP requires
`Npar=40`, an actual selected sample count `m=50..10000`, and ordinarily
`p=0.0005..0.05`. Set `-pb` explicitly when opting in:

```bash
sjaracne local -e expression.exp -g hubs.txt -n 100 \
  --subsample-size 80 -M pchip -pb 0.001 \
  -o results -tmp tmp
```

The same callable is available without running inference or importing SciPy:

```python
from SJARACNe import get_mi_cutoff

tf_cutoff = get_mi_cutoff(m=80, p=0.001)       # approximately 0.147831
sig_cutoff = get_mi_cutoff(m=80, p=0.0005)     # approximately 0.164452
```

Native equivalent: `sjaracne.exe -i expression.exp -s hubs.txt -N 40 -u 80 -M pchip -p 0.001`.
Both TF and SIG use the same function; network identity does not choose a
different default. These probabilities are provisional BRCA100 operating
points, not universal biological recommendations.

PCHIP is **experimental**, not a uniformly validated probability guarantee.
The independent-rank study had 6599/6600 moderate observed-rate checks within
10%, but simultaneous intervals confirmed only 784/800 core cases, with 16
unresolved. Off-grid sample counts such as `m=80` are interpolated and were
not tested in the 50-step broad-grid study. Rare probabilities down to `1e-7`
require `--allow-rare-tail` (native `-R`; Python `allow_rare_tail=True`) and
emit a warning: ten million validation draws cannot establish 10% accuracy
at `1e-7`. No biological FDR control is claimed.
The null assumes exchangeable observations and unique ordinal ranks; the
calibration is not an independent validation of zero-inflated or heavily tied
expression data.

Unsupported inputs, crossed cutoff knots and estimator mismatches fail rather
than extrapolating or silently falling back. Legacy with-replacement `-r` and
adjacency replay `-j` are not supported when PCHIP is active. An explicit
native `-t`, including zero, overrides the p-value conversion.

The model and native coefficients are bundled with SJARACNe; there is no
runtime dependency on the benchmark repository. Network headers retain source
hashes, actual `m`, rare-tail/off-grid status and decimal/hex binary64 cutoffs.
Literal training-size/probability knots retain their exact stored binary64
values. Between knots, library/compiler arithmetic can differ from the study's
NumPy/SciPy implementation by a few ULPs (15 in the recorded port check).
Numerical closeness does not prove identical exceedance counts near MI ties;
the published null-rate summaries describe the original benchmark callable.
See the [process](https://github.com/adamdingliang/sjaracne_benchmark/tree/main/benchmarks/apmi_null_m_sweep/broad_grid_10m)
and [results](https://github.com/adamdingliang/sjaracne_benchmark/blob/main/benchmarks/apmi_null_m_sweep/broad_grid_10m/RESULTS.md).
The earlier GPD calibrator and exact-size models are retained only as
[exploratory benchmark history](https://github.com/adamdingliang/sjaracne_benchmark/tree/main/benchmarks/apmi_null_gpd_exploratory),
not installed as a runtime option.


## Examples to create a transcription factor network
**Note:** for testing purpose, the number of resampled networks (legacy option ```-n```) is set to 2, the consensus p-value threshold
```-pc``` is set to 1.0 in the following examples. ```-n 100``` and ```-pc 1e-5``` are recommended for real 
applications. Note that there is no / at the end of the -o option but there is a / at the end of the -tmp option.
The default ```P-value``` for sjaracne is ```1e-7```. The minimum P-value accepted with the ```-pb argument is 1e-10```.

### Running on a single machine (Linux/OSX) 
```sjaracne local -e ./tests/inputs/BRCA100.exp -g ./tests/inputs/BRCA100_TF.txt -n 2 -o ./results/SJARACNE_out.final -pc 1.0 -tmp ./results/tmp/```

### Running on an IBM LSF cluster
```sjaracne lsf -j ./SJARACNe/config/config_cwlexec.json -e ./tests/inputs/BRCA100.exp -g ./tests/inputs/BRCA100_TF.txt -n 2 -o ./results/SJARACNE_out.final -pc 1.0```


## Reference
Alireza Khatamian, Evan O. Paull, Andrea Califano* & Jiyang Yu*. SJARACNe: a scalable 
software tool for gene network reverse engineering from big data. Bioinformatics (2018). *Corresponding authors.
