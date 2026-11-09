#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: Build a historical probability consensus network
doc: |
  Retained low-level compatibility tool for the historical probability-based
  consensus script. It merges supplied adjacency graphs and does not estimate
  MI or apply DPI. New workflow entry points use the documented pipeline tools.

baseCommand: get_consensus_network.py

inputs:
  adjmat_dir:
    type: Directory
    inputBinding:
      position: 1
    doc: directory with adjacent matrix

  p_thresh_arg:
    type: float
    inputBinding:
      position: 2
    doc: P value threshold

  output_dir:
    type: string
    inputBinding:
      position: 3
    doc: output directory name

outputs:
  out_dir:
    type: Directory
    outputBinding:
      glob: $(inputs.output_dir)
