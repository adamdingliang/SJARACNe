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
          ${ return JSON.stringify({exp_file: inputs.exp_file.path, probe_file: inputs.probe_file.path, aracne_config_dir: inputs.aracne_config_dir.path, settings: inputs.settings.path, adjacencies: inputs.adjacencies.map(function(f) { return f.path; }), run_provenance: inputs.run_provenance.map(function(f) { return f.path; })}); }
arguments: [consensus, request.json]

inputs:
  exp_file: File
  probe_file: File
  aracne_config_dir: Directory
  settings: File
  adjacencies: File[]
  run_provenance: File[]
outputs:
  enhanced_network:
    type: File
    outputBinding: {glob: results/consensus_network_ncol_.txt}
  network:
    type: File
    outputBinding: {glob: results/consensus_network_3col_.txt}
  bootstrap_info:
    type: File
    outputBinding: {glob: results/bootstrap_info_.txt}
  parameter_info:
    type: File
    outputBinding: {glob: results/parameter_info_.txt}
  provenance:
    type: File
    outputBinding: {glob: consensus_provenance.json}
  native_stdout:
    type: File
    outputBinding: {glob: consensus.native.stdout.txt}
  native_stderr:
    type: File
    outputBinding: {glob: consensus.native.stderr.txt}
