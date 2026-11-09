#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
id: main
label: Reconstruct a regulatory network with DPI before consensus
doc: |
  Compatibility workflow for the historical processing order: each scattered
  job resamples observations, estimates MI, filters edges by significance and
  applies DPI. Python then merges these already-pruned graphs by recurrence
  or the deprecated consensus probability and annotates the final network.

  There is no additional DPI pass after consensus in this workflow. Sampling
  and MI settings remain explicit inputs; the public CLI resolves the chosen
  current or legacy profile before invoking this entry point. For the new
  consensus-first order, use consensus_before_dpi.cwl instead.

requirements:
  - class: ScatterFeatureRequirement
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

inputs:
  exp_file:
    type: File
    label: Expression matrix with row IDs used as network node IDs
  probe_file:
    type: File
    label: Hub IDs matching the expression matrix, for either TF or SIG networks
  min_recurrence:
    type: int?
    label: Minimum number of distinct resampled networks containing an edge; defaults to 20 when both cutoff inputs are null
  p_value_consensus:
    type: float?
    label: Deprecated legacy normal-approximation P-value threshold in building consensus network
  p_value_bootstrap:
    type: double?
    label: Explicitly required with PCHIP; omitted legacy probability resolves to 1e-7
  apmi_null_model:
    type: string?
    label: Optional bundled experimental AP-MI PCHIP model name pchip
  tolerance:
    type: double?
    label: Optional fixed DPI tolerance; omitted means strict zero unless dynamic_dpi is enabled
  dynamic_dpi:
    type: boolean
    default: false
    label: Experimental dynamic DPI opt-in; cannot be combined with tolerance
  dpi_href:
    type: int?
    label: Optional dynamic reference source count; native default 479
  dpi_slope:
    type: double?
    label: Optional dynamic epsilon increase per source-count doubling; native default 0.05
  dpi_max:
    type: double?
    label: Optional dynamic epsilon cap; native default 0.30
  dpi_offset:
    type: double?
    label: Additive dynamic-DPI offset before the cap; requires dynamic_dpi
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
    type: string?
    label: Fixed observation count or percentage sampled without replacement in every network
  sample_number:
    type: int?
    label: Explicit legacy full-size replacement sampling; set subsample_spec to null when used
  final_out_dir_name:
    type: string
    label: final output directory name

# Public output keys and historical filenames are preserved. out_dir is an
# enhanced network file, despite its historical directory-like name.
outputs:
  consensus_network:
    type: File
    outputSource: build_consensus_network/consensus_network
  out_dir:
    type: File
    outputSource: build_consensus_network/out_dir
  bootstrap_info:
    type: File
    outputSource: build_consensus_network/bootstrap_info
  parameter_info:
    type: File
    outputSource: build_consensus_network/parameter_info

steps:
  # Validate expression structure, hub IDs and overlap before any resample job.
  validate_inputs:
    label: Validate expression and hub inputs
    run: QC_input.cwl
    in:
      exp_file: exp_file
      probe_file: probe_file
    out: [validation_report]

  # Stable seeds 1..bootstrap_num are paired with the output filenames below.
  create_resample_seeds:
    label: Create deterministic resample seeds
    run: int_to_int_array.cwl
    in:
      number: bootstrap_num
    out: [int_array]

  # Keep historical TF_run_*.adj filenames for both TF and SIG networks.
  create_resample_graph_names:
    label: Create resample adjacency filenames
    run: int_to_str_array.cwl
    in:
      number: bootstrap_num
    out: [str_array]

  # Preserve the historical line-ending check; its omitted output_file exposes
  # the staged input, rather than promising a normalized returned file.
  check_expression_line_endings:
    label: Run the historical expression line-ending check
    run: ch_line_ending.cwl
    in:
      input_file: exp_file
    out: [out_file]

  # Use the same historical check for hub IDs; an explicit output filename
  # is needed when callers require the converted LF copy as the returned file.
  check_hub_line_endings:
    label: Run the historical hub-list line-ending check
    run: ch_line_ending.cwl
    in:
      input_file: probe_file
    out: [out_file]

  # One job per seed: resample -> estimate MI -> significance filter -> DPI.
  # The sampling inputs select unique observations or legacy replacement draws.
  resample_estimate_mi_and_apply_dpi:
    label: Resample, estimate MI and apply DPI
    run: sjaracne.cwl
    in:
      preflight_report: validate_inputs/validation_report
      exp_file: check_expression_line_endings/out_file
      probe_file_tf: check_hub_line_endings/out_file
      probe_file_subnetwork: check_hub_line_endings/out_file
      p_value:
        source: p_value_bootstrap
        valueFrom: |
          ${
            if (self === null || typeof self === 'undefined') {
              if (inputs.apmi_null_model === 'pchip') {
                throw new Error('PCHIP requires explicit p_value_bootstrap; no PCHIP probability default is chosen');
              }
              return 1e-7;
            }
            return self;
          }
      apmi_null_model: apmi_null_model
      tolerance: tolerance
      dynamic_dpi: dynamic_dpi
      dpi_href: dpi_href
      dpi_slope: dpi_slope
      dpi_max: dpi_max
      dpi_offset: dpi_offset
      aracne_config_dir: aracne_config_dir
      npar_limit: depth
      subsample_spec:
        source: subsample_spec
        valueFrom: |
          ${
            var hasSpec = self !== null && typeof self !== 'undefined';
            var replacement = inputs.sample_number !== null &&
                              typeof inputs.sample_number !== 'undefined' &&
                              inputs.sample_number > 0;
            if (replacement) {
              if (hasSpec) {
                throw new Error('Replacement sampling conflicts with subsample_spec');
              }
              return null;
            }
            return hasSpec ? self : '80%';
          }
      sample_number: sample_number
      output_file_name: create_resample_graph_names/str_array
      seed: create_resample_seeds/int_array
    scatter: [output_file_name, seed]
    scatterMethod: dotproduct
    out: [out_adj]

  # Collect all already-pruned resample graphs into the consensus input directory.
  collect_resample_networks:
    label: Collect the DPI-pruned resample networks
    run: copy_files_to_dir.cwl
    in:
      input_files: resample_estimate_mi_and_apply_dpi/out_adj
      dirname: final_out_dir_name
    out: [out_dir]

  # Merge the DPI-pruned graphs in Python and annotate from full expression.
  # DPI has already run in each resample; this step does not apply it again.
  build_consensus_network:
    label: Build and annotate consensus from pruned networks
    run: create_consensus_network.cwl
    in:
      adjmat_dir: collect_resample_networks/out_dir
      min_recurrence: min_recurrence
      p_thresh_arg: p_value_consensus
      exp_mat: check_expression_line_endings/out_file
      output_dir: final_out_dir_name
    out: [out_dir, consensus_network, bootstrap_info, parameter_info]
