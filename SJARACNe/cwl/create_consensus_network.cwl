#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: [python3, -m, SJARACNe.bin.create_consensus_network]

inputs:
  adjmat_dir:
    type: Directory
    inputBinding:
      position: 1
      prefix: -a
    doc: directory with adjacent matrix

  min_recurrence:
    type: int?
    inputBinding:
      position: 2
      prefix: -k
    doc: Minimum number of distinct bootstrap networks containing an edge; defaults to 6 when both cutoff inputs are null

  p_thresh_arg:
    type: float?
    inputBinding:
      position: 3
      prefix: -p
    doc: Deprecated legacy normal-approximation P value threshold

  exp_mat:
    type: File
    inputBinding:
      position: 4
      prefix: -e
    doc: expression matrix file

  output_dir:
    type: string
    inputBinding:
      position: 5
      prefix: -o
    doc: output directory name

  subnet:
    type: File?
    inputBinding:
      position: 6
      prefix: -s
    doc: file with gene symbols of interest to build a subnet

outputs:
  out_dir:
    type: File
    outputBinding:
      glob: $(inputs.output_dir)/consensus_network_ncol_.txt

  bootstrap_info:
    type: File
    outputBinding:
      glob: $(inputs.output_dir)/bootstrap_info_.txt

  parameter_info:
    type: File
    outputBinding:
      glob: $(inputs.output_dir)/parameter_info_.txt
