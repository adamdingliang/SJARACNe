#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: Scalable solution of ARACNe that dramatically improves the computational performance.

baseCommand: sjaracne.exe

requirements:
  - class: InlineJavascriptRequirement

arguments:
  - position: 6
    valueFrom: |
      ${
        function hasValue(value) {
          return value !== null && typeof value !== 'undefined';
        }
        var fixed = inputs.tolerance;
        var tuning = hasValue(inputs.dpi_href) || hasValue(inputs.dpi_slope) || hasValue(inputs.dpi_max);
        if (inputs.dynamic_dpi) {
          if (hasValue(fixed)) {
            throw new Error('dynamic_dpi conflicts with an explicit fixed tolerance');
          }
          return [];
        }
        if (tuning) {
          throw new Error('dpi_href, dpi_slope and dpi_max require dynamic_dpi');
        }
        return ['-e', String(hasValue(fixed) ? fixed : 0)];
      }

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
    type: string?
    inputBinding:
      position: 4
      prefix: -M
    doc: Optional bundled model name pchip; experimental, requires Npar 40
  p_value:
    type: double
    inputBinding:
      position: 5
      prefix: -p
  tolerance:
    type: double?
    doc: Optional fixed DPI tolerance; omitted resolves to strict zero unless dynamic_dpi is enabled
  dynamic_dpi:
    type: boolean
    default: false
    inputBinding:
      position: 6
      prefix: --dynamic-dpi
    doc: Explicit experimental dynamic epsilon opt-in; conflicts with tolerance
  dpi_href:
    type: int?
    inputBinding:
      position: 6
      prefix: --dpi-href
    doc: Optional dynamic reference source count; requires dynamic_dpi; native default 479
  dpi_slope:
    type: double?
    inputBinding:
      position: 6
      prefix: --dpi-slope
    doc: Optional dynamic epsilon increase per source-count doubling; requires dynamic_dpi; native default 0.05
  dpi_max:
    type: double?
    inputBinding:
      position: 6
      prefix: --dpi-max
    doc: Optional dynamic epsilon cap in [0,1); requires dynamic_dpi; native default 0.30
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
