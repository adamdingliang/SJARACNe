#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow
id: main
label: Reconstruct a regulatory network with consensus before DPI
doc: |
  Prepare inputs, scatter seeded resamples to estimate significant MI edges,
  then retain recurrent edges, average MI over supporting runs and apply DPI
  once to the consensus graph. Finally annotate the pruned network using the
  full expression matrix.

  This is the public CLI's consensus-first CWL entry point. The CLI supplies
  resolved scientific_options with consensus_order set to consensus-first;
  direct callers must supply the same order. The shared adapter revalidates
  settings and provenance. CWL schedules jobs and localizes their input files.
  Use dpi_before_consensus.cwl for the historical DPI-first processing order.
requirements:
  ScatterFeatureRequirement: {}
inputs:
  exp_file: File
  probe_file: File
  aracne_config_dir: Directory
  scientific_options: Any
  bootstrap_num: int
  verify_raw_hashes:
    type: boolean
    default: false
    doc: Re-read raw graph contents to verify producer SHA-256 hashes before consensus
# Public output keys and historical filenames are preserved. out_dir is an
# enhanced network file, despite its historical directory-like name.
outputs:
  out_dir: {type: File, outputSource: consensus_and_dpi/enhanced_network}
  consensus_network: {type: File, outputSource: consensus_and_dpi/network}
  bootstrap_info: {type: File, outputSource: consensus_and_dpi/bootstrap_info}
  parameter_info: {type: File, outputSource: consensus_and_dpi/parameter_info}
  parameters: {type: File, outputSource: prepare_inputs/parameters}
  # Export scatter artifacts into the owned engine output directory so default
  # cwltool cleanup does not remove the auditable resample evidence.
  adjacencies: {type: 'File[]', outputSource: resample_and_estimate_mi/adjacency}
  run_provenance: {type: 'File[]', outputSource: resample_and_estimate_mi/provenance}
  native_stdout: {type: 'File[]', outputSource: resample_and_estimate_mi/native_stdout}
  native_stderr: {type: 'File[]', outputSource: resample_and_estimate_mi/native_stderr}
  consensus_provenance: {type: File, outputSource: consensus_and_dpi/provenance}
  consensus_native_stdout: {type: File, outputSource: consensus_and_dpi/native_stdout}
  consensus_native_stderr: {type: File, outputSource: consensus_and_dpi/native_stderr}
steps:
  # Runs once: validate IDs/settings and prepare LF-normalized copies.
  prepare_inputs:
    label: Validate and prepare expression and hub inputs
    run: prepare_inputs.cwl
    in:
      exp_file: exp_file
      probe_file: probe_file
      aracne_config_dir: aracne_config_dir
      scientific_options: scientific_options
      bootstrap_num: bootstrap_num
    out: [expression, hubs, parameters, validation_report]
  # Stable seeds 1..bootstrap_num identify independent resample jobs.
  create_resample_seeds:
    label: Create deterministic resample seeds
    run: int_to_int_array.cwl
    in: {number: bootstrap_num}
    out: [int_array]
  # One job per seed: sample observations, estimate MI and apply its cutoff.
  # consensus-first uses --raw-mi, so these graphs have not undergone DPI.
  resample_and_estimate_mi:
    label: Resample observations and estimate significant MI edges
    run: resample_and_estimate_mi.cwl
    in:
      exp_file: prepare_inputs/expression
      probe_file: prepare_inputs/hubs
      aracne_config_dir: aracne_config_dir
      settings: prepare_inputs/parameters
      seed: create_resample_seeds/int_array
    scatter: seed
    out: [adjacency, provenance, native_stdout, native_stderr]
  # Wait for every seed and verify its provenance before aggregation.
  # One C++ invocation performs K support -> supporting-run mean MI -> DPI.
  consensus_and_dpi:
    label: Build the consensus graph and apply DPI
    run: consensus_and_dpi.cwl
    in:
      exp_file: prepare_inputs/expression
      probe_file: prepare_inputs/hubs
      aracne_config_dir: aracne_config_dir
      settings: prepare_inputs/parameters
      adjacencies: resample_and_estimate_mi/adjacency
      run_provenance: resample_and_estimate_mi/provenance
      verify_raw_hashes: verify_raw_hashes
    out: [enhanced_network, network, bootstrap_info, parameter_info, provenance, native_stdout, native_stderr]
