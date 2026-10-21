# SJARACNe

develop: [![Regression](https://github.com/adamdingliang/SJARACNe/actions/workflows/regression.yml/badge.svg?branch=develop&event=push)](https://github.com/adamdingliang/SJARACNe/actions/workflows/regression.yml?query=branch%3Adevelop)

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
queue and adjust memory reservation for each step is necessary. Consider increasing memory reservation for bootstrap 
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
also outputs two meta information files: parameter_info_.txt and bootstrap_info_.txt, which stores SJARACNe 
input parameters and bootstrap parameters respectively.


## Examples to create a transcription factor network
**Note:** for testing purpose, the number of bootstraps (```-n```) is set to 2, the consensus p-value threshold 
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
