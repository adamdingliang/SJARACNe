#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
  - class: ScatterFeatureRequirement
  - class: InlineJavascriptRequirement

inputs:
  exp_file:
    type: File
    label: expression matrix file, row indexes are used as the nodes in the network
  probe_file:
    type: File
    label: file with a list of symbols annotated as transcription factors (hub genes) for constructing subnetworks
  min_recurrence:
    type: int?
    label: Minimum number of distinct resampled networks containing an edge; defaults to 6 when both cutoff inputs are null
  p_value_consensus:
    type: float?
    label: Deprecated legacy normal-approximation P-value threshold in building consensus network
  p_value_bootstrap:
    type: float
    default: 1e-7
    label: P-value threshold in building individual resampled networks
  apmi_null_model:
    type: File?
    label: Optional exact-m and exact-depth estimator-matched AP-MI null-tail model
  depth:
    type: int
    default: 40
    label: maximum partitioning depth
  aracne_config_dir:
    type: Directory
    label: Directory containing ARACNe configuration files, default is current working directory
  bootstrap_num:
    type: int
    default: 100
    label: Number of resampled networks to generate (legacy input name)
  subsample_spec:
    type: string
    default: "80%"
    label: Fixed observation count or percentage sampled without replacement in every network
  final_out_dir_name:
    type: string
    label: final output directory name

outputs:
  out_dir:
    type: File
    outputSource: consensus/out_dir
  bootstrap_info:
    type: File
    outputSource: consensus/bootstrap_info
  parameter_info:
    type: File
    outputSource: consensus/parameter_info

steps:
  # Step 0: validate input file
  validate_files:
    run: QC_input.cwl
    in:
      exp_file: exp_file
      probe_file: probe_file
    out: [validation_report]

  # Step 1: create seeds from the resampled-network count (legacy input name)
  create_seeds:
    run: int_to_int_array.cwl
    in:
      number: bootstrap_num
    out: [int_array]

  # Step 2: create adjacency file names from the resampled-network count
  create_adjmat_names:
    run: int_to_str_array.cwl
    in:
      number: bootstrap_num
    out: [str_array]

  # Step 3: change expression file line ending
  ch_ending_exp:
    run: ch_line_ending.cwl
    in:
      input_file: exp_file
    out: [out_file]

  # Step 4: change probe file line ending
  ch_ending_probe:
    run: ch_line_ending.cwl
    in:
      input_file: probe_file
    out: [out_file]

  # Step 5: fixed-size subsampling without replacement using different seeds
  bootstrap:
    run: sjaracne.cwl
    in:
      preflight_report: validate_files/validation_report
      exp_file: ch_ending_exp/out_file
      probe_file_tf: ch_ending_probe/out_file
      probe_file_subnetwork: ch_ending_probe/out_file
      p_value: p_value_bootstrap
      apmi_null_model: apmi_null_model
      aracne_config_dir: aracne_config_dir
      npar_limit: depth
      subsample_spec: subsample_spec
      output_file_name: create_adjmat_names/str_array
      seed: create_seeds/int_array
    scatter: [output_file_name, seed]
    scatterMethod: dotproduct
    out: [out_adj]

  # Step 6: copy output adjacent matrix files to final output directory
  copy_to_dir:
    run: copy_files_to_dir.cwl
    in:
      input_files: bootstrap/out_adj
      dirname: final_out_dir_name
    out: [out_dir]

  # Step 7: generate a consensus network
  consensus:
    run: create_consensus_network.cwl
    in:
      adjmat_dir: copy_to_dir/out_dir
      min_recurrence: min_recurrence
      p_thresh_arg: p_value_consensus
      exp_mat: ch_ending_exp/out_file
      output_dir: final_out_dir_name
    out: [out_dir, bootstrap_info, parameter_info]
