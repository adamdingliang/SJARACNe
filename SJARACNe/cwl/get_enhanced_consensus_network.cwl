#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: Annotate an existing consensus network
doc: |
  Retained low-level compatibility tool that annotates an existing
  three-column consensus network using the supplied expression matrix.
  This stage neither reconstructs MI edges nor applies DPI.

requirements:
  - class: InlineJavascriptRequirement

baseCommand: get_enhanced_consensus_network.py

inputs:
  exp_file:
    type: File
    inputBinding:
      position: 1
    doc: Input gene expression profile dataset (required)

  consensus_network_3col:
    type: Directory
    inputBinding:
      valueFrom: $(self.path)/consensus_network_3col_.txt
      position: 2
    doc: Consensus network file with 3 columns

  output_dir:
    type: string
    inputBinding:
      position: 3
    doc: Output directory name

outputs:
  out_dir:
    type: Directory
    outputBinding:
      glob: $(inputs.output_dir)
