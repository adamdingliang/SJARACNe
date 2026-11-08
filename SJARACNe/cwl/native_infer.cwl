#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: [python3, -m, SJARACNe.bin.wdl_workflow]
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: request.json
        entry: |
          ${ return JSON.stringify({exp_file: inputs.exp_file.path, probe_file: inputs.probe_file.path, aracne_config_dir: inputs.aracne_config_dir.path, settings: inputs.settings.path, seed: inputs.seed}); }
arguments: [infer, request.json]

inputs:
  exp_file: File
  probe_file: File
  aracne_config_dir: Directory
  settings: File
  seed: int
outputs:
  adjacency:
    type: File
    outputBinding: {glob: 'TF_run_*.adj'}
  provenance:
    type: File
    outputBinding: {glob: run_provenance.json}
  native_stdout:
    type: File
    outputBinding: {glob: native.stdout.txt}
  native_stderr:
    type: File
    outputBinding: {glob: native.stderr.txt}
