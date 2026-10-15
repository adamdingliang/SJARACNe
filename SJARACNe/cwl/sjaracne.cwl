#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: Scalable solution of ARACNe that dramatically improves the computational performance.

baseCommand: sjaracne.exe

inputs:
  preflight_report:
    type: File?
    doc: Successful hub-overlap validation report; establishes a workflow dependency before resampling
  exp_file:
    type: File
    inputBinding:
      position: 1
      prefix: -i
    doc: Input gene expression profile dataset (required)
  probe_file_tf:
    type: File
    inputBinding:
      position: 2
      prefix: -l
    doc: File containing a list of probes annotated as transcription factors in the input dataset
  probe_file_subnetwork:
    type: File
    inputBinding:
      position: 3
      prefix: -s
    doc: File containing a list of probes for which a subnetwork will be constructed
  apmi_null_model:
    type: File?
    inputBinding:
      position: 4
      prefix: -M
    doc: Optional exact-m, exact-Npar estimator-matched AP-MI GPD-tail model
  p_value:
    type: float
    inputBinding:
      position: 5
      valueFrom: |
        ${
          return inputs.p_value.toFixed(10); // Format to 10 decimal places (1e-10)
        }
      prefix: -p
  tolerance:
    type: int
    default: 0
    inputBinding:
      position: 6
      prefix: -e
    doc: DPI tolerance
  algorithm:
    type: string
    default: adaptive_partitioning
    inputBinding:
      position: 7
      prefix: -a
    doc: algorithm
  sample_number:
    type: int?
    inputBinding:
      position: 8
      prefix: -r
    doc: Deprecated legacy full-size bootstrap with replacement; use subsample_spec for new analyses
  subsample_spec:
    type: string?
    inputBinding:
      position: 8
      prefix: -u
    doc: Fixed-size sampling without replacement as an exact count or explicit percentage, for example 80 or 80%
  aracne_config_dir:
    type: Directory
    inputBinding:
      position: 9
      prefix: -H
    doc: Directory containing ARACNe configuration files, default is current working directory
  npar_limit:
    type: int
    default: 20
    inputBinding:
      position: 10
      prefix: -N
    doc: Maximum allowed value of npar
  output_file_name:
    type: string
    inputBinding:
      position: 11
      prefix: -o
    doc: Output file name (optional)
  seed:
    type: int
    default: 1
    inputBinding:
      position: 12
      prefix: -S
    doc: Initial seed for random number generator

outputs:
  out_adj:
    type: File
    outputBinding:
      glob: $(inputs.output_file_name)
