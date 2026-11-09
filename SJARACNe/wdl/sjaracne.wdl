version 1.3

# Reconstruct a regulatory network from independently resampled MI estimates.
# Current default: 100 seeded 80% resamples -> significant MI edges -> inclusive
# K20 consensus -> supporting-run mean MI -> one final DPI pass. The dynamic
# tolerance adds 0.025 before the cap. Explicit dpi-first and legacy settings
# instead prune each resample before building consensus; see the task comments.
#
# WDL schedules tasks and their dependencies. The compiled executable performs
# MI estimation and default consensus/DPI; the shared Python adapter validates
# task boundaries and enhances the final graph from the full expression matrix.
# Each task receives a JSON request, so user values are not shell-interpolated.

# Scientific settings travel together, independently of the execution engine.
# Optional fields remain unset until workflow_options.resolve_options applies
# the selected profile and rejects incompatible combinations.
struct WorkflowOptions {
    String profile
    String? network_type
    String? sampling_method
    Int bootstrap_num
    String? subsample_spec
    Int depth
    Float? p_value_bootstrap
    String? apmi_null_model
    Int? min_recurrence
    Float? p_value_consensus
    Float? tolerance
    Boolean? dynamic_dpi
    Int? dpi_href
    Float? dpi_slope
    Float? dpi_max
    Float? dpi_offset
    String? consensus_order
}

# Keep the SJARACNe workflow identifier stable for existing input/output JSON
# keys (SJARACNe.*). The workflow and task display names describe their actions.
workflow SJARACNe {
    input {
        # Prepared expression matrix and source IDs; this workflow does not
        # normalize expression values, filter genes or choose the hub list.
        File exp_file
        File probe_file
        Directory? aracne_config_dir

        # "legacy" selects historical settings in the current executable,
        # not an old executable or a different execution engine.
        String profile = "current"
        # Current profile requires TF/SIG to choose its default null probability;
        # the supplied probe_file, not this label, determines the actual sources.
        String? network_type

        # Profile defaults: current=80% without replacement; legacy=N draws with
        # replacement. An explicit subsample_spec is a count or percentage string.
        String? sampling_method
        # Total networks, not concurrent workers. Seeds are 1..bootstrap_num.
        Int bootstrap_num = 100
        String? subsample_spec
        # The frozen PCHIP model requires the AP-MI Npar limit to remain 40.
        Int depth = 40

        # Null-tail probability, not biological FDR: current TF=0.001,
        # SIG=0.0005; legacy=1e-7. Explicit supported values override defaults.
        Float? p_value_bootstrap
        String? apmi_null_model

        # Mutually exclusive: current defaults to K=20 distinct-network support
        # (inclusive), whereas legacy uses the historical pc=1e-5 calculation.
        Int? min_recurrence
        Float? p_value_consensus

        # Larger epsilon relaxes DPI pruning. Fixed tolerance and dynamic=true
        # conflict; omitted current settings select dynamic, legacy selects zero.
        Float? tolerance
        Boolean? dynamic_dpi
        # Native dynamic schedule: reference hub count, log2 slope and upper cap.
        Int? dpi_href
        Float? dpi_slope
        Float? dpi_max
        # Added before the cap; consensus-first defaults to 0.025, dpi-first to zero.
        Float? dpi_offset
        # Omitted current order is consensus-first; legacy requires dpi-first.
        String? consensus_order

        # Per-task scheduling requests, not the backend's total memory budget.
        # The CLI's --jobs/--memory configure Sprocket separately, not WDL inputs.
        # On the trusted local host backend these are not hard process RAM limits.
        String prepare_memory = "2 GiB"
        String inference_memory = "4 GiB"
        String consensus_memory = "8 GiB"
        # Operational integrity check, separate from scientific parameters.
        # False preserves metadata/seed checks without a full raw-file reread.
        Boolean verify_raw_hashes = false
    }

    # Preserve omissions so direct WDL requests and the CLI share one resolver.
    WorkflowOptions options = WorkflowOptions {
        profile: profile,
        network_type: network_type,
        sampling_method: sampling_method,
        bootstrap_num: bootstrap_num,
        subsample_spec: subsample_spec,
        depth: depth,
        p_value_bootstrap: p_value_bootstrap,
        apmi_null_model: apmi_null_model,
        min_recurrence: min_recurrence,
        p_value_consensus: p_value_consensus,
        tolerance: tolerance,
        dynamic_dpi: dynamic_dpi,
        dpi_href: dpi_href,
        dpi_slope: dpi_slope,
        dpi_max: dpi_max,
        dpi_offset: dpi_offset,
        consensus_order: consensus_order,
    }

    # Prepare once: validate expression/hub IDs and resolve the scientific
    # options, then stage LF-normalized copies without changing numeric values.
    # parameters.json binds those options to input and runtime/config hashes.
    call prepare_inputs { input:
        exp_file = exp_file,
        probe_file = probe_file,
        aracne_config_dir = aracne_config_dir,
        options = options,
        task_memory = prepare_memory,
    }

    # One job per seed resamples observations, estimates AP-MI and retains
    # edges passing the configured null-tail cutoff (PCHIP by current default).
    # Consensus-first uses --raw-mi: these outputs have NO per-resample DPI.
    # Explicit dpi-first/legacy applies DPI here before writing each graph.
    # Backend CPU, memory and scatter limits determine actual concurrency.
    scatter (index in range(bootstrap_num)) {
        call resample_and_estimate_mi { input:
            exp_file = prepare_inputs.expression,
            probe_file = prepare_inputs.hubs,
            aracne_config_dir = aracne_config_dir,
            settings = prepare_inputs.parameters,
            seed = index + 1,
            task_memory = inference_memory,
        }
    }

    # Wait for every seed and verify matching graph/provenance receipts.
    # Consensus-first: retain edges with support >= K, average MI only over
    # supporting runs, then apply DPI inside one C++ invocation. No candidate
    # consensus graph is written and reread between aggregation and DPI.
    # DPI-first/legacy: Python merges graphs already pruned in the scatter;
    # there is no additional DPI pass here. Both orders enhance the final graph.
    call build_consensus_network { input:
        exp_file = prepare_inputs.expression,
        probe_file = prepare_inputs.hubs,
        aracne_config_dir = aracne_config_dir,
        settings = prepare_inputs.parameters,
        adjacencies = resample_and_estimate_mi.adjacency,
        run_provenance = resample_and_estimate_mi.provenance,
        verify_raw_hashes = verify_raw_hashes,
        task_memory = consensus_memory,
    }

    output {
        # Historical name retained for compatibility: out_dir is an enhanced
        # network FILE, not a directory. The CLI publishes these four final files.
        File out_dir = build_consensus_network.enhanced_network
        File consensus_network = build_consensus_network.network
        File bootstrap_info = build_consensus_network.bootstrap_info
        File parameter_info = build_consensus_network.parameter_info
        # Per-seed artifacts and prepared settings remain available for auditing.
        Array[File] adjacencies = resample_and_estimate_mi.adjacency
        Array[File] run_provenance = resample_and_estimate_mi.provenance
        Array[File] native_stdout = resample_and_estimate_mi.native_stdout
        Array[File] native_stderr = resample_and_estimate_mi.native_stderr
        File validation_report = prepare_inputs.validation_report
        File parameters = prepare_inputs.parameters
        File consensus_provenance = build_consensus_network.provenance
        File consensus_native_stdout = build_consensus_network.native_stdout
        File consensus_native_stderr = build_consensus_network.native_stderr
    }

    meta {
        name: "Reconstruct a regulatory network"
        description: "Prepare inputs, resample and estimate MI, then build the consensus network. Current defaults perform K20/supporting-run mean MI before DPI; explicit dpi-first and legacy profiles prune each resample before consensus."
    }
}

# Validate the input bundle once before starting any resample job. Resolve
# profile/order settings, stage safe LF copies and fingerprint the Python
# adapter, compiled executable, input files and model/config. No expression
# normalization, gene filtering or hub selection is performed by this task.
task prepare_inputs {
    input {
        File exp_file
        File probe_file
        Directory? aracne_config_dir
        WorkflowOptions options
        String task_memory = "2 GiB"
    }

    # Pass data via a localized JSON file, not shell-expanded user arguments.
    env File request = write_json(object {
        exp_file: exp_file,
        probe_file: probe_file,
        aracne_config_dir: aracne_config_dir,
        parameters: options,
    })

    command <<<
        set -euo pipefail
        python3 -m SJARACNe.bin.wdl_workflow prepare "$request"
    >>>

    output {
        File expression = "expression.exp"
        File hubs = "hubs.txt"
        File parameters = "parameters.json"
        File validation_report = "hub_overlap_validation.txt"
    }

    requirements {
        cpu: 1
        memory: task_memory
    }

    meta {
        name: "Prepare and validate network inputs"
        description: "Validate IDs and scientific settings, stage expression and hub files, and record input/runtime fingerprints before resampling."
    }
}

# Resample observations and estimate significant MI edges for one seed.
# Current defaults select 80% without replacement; legacy samples N observations
# with replacement. The shared adapter rechecks prepared hashes before invoking
# the executable and saves its exact command, output hashes and stdout/stderr.
# In consensus-first mode --raw-mi skips DPI and preserves all source rows,
# including empty rows, plus full-precision MI. With explicit dpi-first/legacy,
# this same task also applies the configured per-resample DPI before output.
task resample_and_estimate_mi {
    input {
        File exp_file
        File probe_file
        Directory? aracne_config_dir
        File settings
        Int seed
        String task_memory = "4 GiB"
    }

    env File request = write_json(object {
        exp_file: exp_file,
        probe_file: probe_file,
        aracne_config_dir: aracne_config_dir,
        settings: settings,
        seed: seed,
    })

    command <<<
        set -euo pipefail
        python3 -m SJARACNe.bin.wdl_workflow infer "$request"
    >>>

    output {
        # TF_run_* is the retained native filename convention for both TF and SIG.
        File adjacency = glob("TF_run_*.adj")[0]
        File provenance = "run_provenance.json"
        File native_stdout = "native.stdout.txt"
        File native_stderr = "native.stderr.txt"
    }

    requirements {
        cpu: 1
        memory: task_memory
    }

    meta {
        name: "Resample observations and estimate MI"
        description: "Estimate significant MI edges for one seeded resample. Consensus-first emits raw MI without DPI; explicit dpi-first and legacy profiles also prune this resample."
    }
}

# Build the final network after all resample jobs have succeeded. Verify the
# complete seed set and matching input/runtime/graph receipts before aggregation.
# Consensus-first runs K-support selection -> supporting-run arithmetic mean MI
# -> DPI in C++, using localized graph paths without adapter staging copies or
# candidate-weight files. Missing-run edges do not contribute zeros to the mean.
# DPI-first/legacy uses Python consensus on already-pruned graphs instead.
# Finally annotate the retained edges using the full input expression matrix;
# final filenames and the public workflow output contract remain unchanged.
task build_consensus_network {
    input {
        File exp_file
        File probe_file
        Directory? aracne_config_dir
        File settings
        Array[File] adjacencies
        Array[File] run_provenance
        Boolean verify_raw_hashes = false
        String task_memory = "8 GiB"
    }

    env File request = write_json(object {
        exp_file: exp_file,
        probe_file: probe_file,
        aracne_config_dir: aracne_config_dir,
        settings: settings,
        adjacencies: adjacencies,
        run_provenance: run_provenance,
        verify_raw_hashes: verify_raw_hashes,
    })

    command <<<
        set -euo pipefail
        python3 -m SJARACNe.bin.wdl_workflow consensus "$request"
    >>>

    output {
        File enhanced_network = "results/consensus_network_ncol_.txt"
        File network = "results/consensus_network_3col_.txt"
        File bootstrap_info = "results/bootstrap_info_.txt"
        File parameter_info = "results/parameter_info_.txt"
        File provenance = "consensus_provenance.json"
        File native_stdout = "consensus.native.stdout.txt"
        File native_stderr = "consensus.native.stderr.txt"
    }

    requirements {
        cpu: 1
        memory: task_memory
    }

    meta {
        name: "Build the consensus network"
        description: "Verify all resample receipts and build the final annotated network. Consensus-first performs K-support and supporting-run mean MI followed by DPI; dpi-first merges previously pruned graphs."
    }
}
