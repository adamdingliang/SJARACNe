#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
id: main
label: Build a consensus graph and apply DPI
doc: |
  After every resample succeeds, verify the complete seed set and matching
  input/runtime fingerprints. With consensus-first settings, one C++ call
  retains edges with at least K supporting runs, averages MI over supporting
  runs only, and applies DPI directly to that consensus graph in memory.
  Annotate only the final pruned network using the full expression matrix.
  The internal adapter command remains consensus for interface compatibility.
baseCommand: [python3, -m, SJARACNe.bin.wdl_workflow]
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: request.json
        entry: |
          ${ return JSON.stringify({exp_file: inputs.exp_file.path, probe_file: inputs.probe_file.path, aracne_config_dir: inputs.aracne_config_dir.path, settings: inputs.settings.path, adjacencies: inputs.adjacencies.map(function(f) { return f.path; }), run_provenance: inputs.run_provenance.map(function(f) { return f.path; }), verify_raw_hashes: inputs.verify_raw_hashes}); }
arguments: [consensus, request.json]

inputs:
  exp_file: File
  probe_file: File
  aracne_config_dir: Directory
  settings: File
  adjacencies: File[]
  run_provenance: File[]
  verify_raw_hashes:
    type: boolean
    default: false
    doc: Verify raw graph SHA-256 hashes; metadata and seed checks always run
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
