#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: Validate expression and hub inputs
doc: |
  Check expression matrix structure, hub identifiers and hub-expression
  overlap. The validation report is a hard dependency of the DPI-first
  resample jobs, so invalid inputs fail before native inference starts.

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.exp_file)
      - $(inputs.probe_file)

baseCommand: [python3, -m, SJARACNe.bin.QC_input]

inputs:
  exp_file:
    type: File
    inputBinding:
      position: 1
      prefix: -e
      valueFrom: $(self.basename)
  probe_file:
    type: File
    inputBinding:
      position: 2
      prefix: -g
      valueFrom: $(self.basename)
  output_file:
    type: string
    default: hub_overlap_validation.txt
    inputBinding:
      position: 3
      prefix: -o

outputs:
  validation_report:
    type: File
    outputBinding:
      glob: $(inputs.output_file)
