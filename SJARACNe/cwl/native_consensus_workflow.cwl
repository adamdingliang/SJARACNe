#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow
# The CLI supplies resolved scientific settings. The same adapter as WDL
# independently revalidates them; CWL schedules the scatter and localizes files.
requirements:
  ScatterFeatureRequirement: {}
inputs:
  exp_file: File
  probe_file: File
  aracne_config_dir: Directory
  scientific_options: Any
  bootstrap_num: int
outputs:
  out_dir: {type: File, outputSource: consensus/enhanced_network}
  consensus_network: {type: File, outputSource: consensus/network}
  bootstrap_info: {type: File, outputSource: consensus/bootstrap_info}
  parameter_info: {type: File, outputSource: consensus/parameter_info}
  parameters: {type: File, outputSource: prepare/parameters}
  # Export scatter artifacts into the owned engine output directory so default
  # cwltool cleanup does not remove the auditable resample evidence.
  adjacencies: {type: 'File[]', outputSource: infer/adjacency}
  run_provenance: {type: 'File[]', outputSource: infer/provenance}
  native_stdout: {type: 'File[]', outputSource: infer/native_stdout}
  native_stderr: {type: 'File[]', outputSource: infer/native_stderr}
  consensus_provenance: {type: File, outputSource: consensus/provenance}
  consensus_native_stdout: {type: File, outputSource: consensus/native_stdout}
  consensus_native_stderr: {type: File, outputSource: consensus/native_stderr}
steps:
  prepare:
    run: native_prepare.cwl
    in:
      exp_file: exp_file
      probe_file: probe_file
      aracne_config_dir: aracne_config_dir
      scientific_options: scientific_options
      bootstrap_num: bootstrap_num
    out: [expression, hubs, parameters, validation_report]
  seeds:
    run: int_to_int_array.cwl
    in: {number: bootstrap_num}
    out: [int_array]
  infer:
    run: native_infer.cwl
    in:
      exp_file: prepare/expression
      probe_file: prepare/hubs
      aracne_config_dir: aracne_config_dir
      settings: prepare/parameters
      seed: seeds/int_array
    scatter: seed
    out: [adjacency, provenance, native_stdout, native_stderr]
  consensus:
    run: native_consensus.cwl
    in:
      exp_file: prepare/expression
      probe_file: prepare/hubs
      aracne_config_dir: aracne_config_dir
      settings: prepare/parameters
      adjacencies: infer/adjacency
      run_provenance: infer/provenance
    out: [enhanced_network, network, bootstrap_info, parameter_info, provenance, native_stdout, native_stderr]
