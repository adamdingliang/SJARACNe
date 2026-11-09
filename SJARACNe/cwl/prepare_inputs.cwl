#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
id: main
label: Validate and prepare network inputs
doc: |
  Validate expression and hub IDs, resolve scientific settings and stage
  LF-normalized input copies without modifying the original files. Record
  input, runtime and configuration fingerprints for all downstream tasks.
  The JSON request is also used by WDL; shell interpolation is not required.
baseCommand: [python3, -m, SJARACNe.bin.wdl_workflow]
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: request.json
        entry: |
          ${ return JSON.stringify({exp_file: inputs.exp_file.path, probe_file: inputs.probe_file.path, aracne_config_dir: inputs.aracne_config_dir.path, parameters: inputs.scientific_options, bootstrap_num: inputs.bootstrap_num}); }
# Adapter subcommands are stable interfaces shared with WDL.
arguments: [prepare, request.json]

inputs:
  exp_file: File
  probe_file: File
  aracne_config_dir: Directory
  scientific_options: Any
  bootstrap_num: int
outputs:
  expression:
    type: File
    outputBinding: {glob: expression.exp}
  hubs:
    type: File
    outputBinding: {glob: hubs.txt}
  parameters:
    type: File
    outputBinding: {glob: parameters.json}
  validation_report:
    type: File
    outputBinding: {glob: hub_overlap_validation.txt}
