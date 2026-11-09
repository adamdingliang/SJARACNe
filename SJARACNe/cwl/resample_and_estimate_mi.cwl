#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
id: main
label: Resample observations and estimate significant MI edges
doc: |
  Run one seed: verify prepared inputs, resample observations, estimate MI
  and retain edges above the configured significance cutoff. With this
  workflow's consensus-first settings, the adapter passes --raw-mi to the
  C++ executable so DPI is disabled in every resample. Retain the raw graph,
  its provenance and native logs for the post-scatter consensus task.
  The internal adapter command remains infer for interface compatibility.
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
  # Historical TF_run_* filenames are used for TF and SIG roles alike.
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
