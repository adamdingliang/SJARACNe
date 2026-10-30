version 1.3

# Orchestration only: the native executable computes AP-MI and DPI, while the
# existing Python consensus module combines and enhances the inferred networks.
# Each task sends a JSON request to the shared Python adapter; user-provided
# values are not interpolated into shell commands.

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
}

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

        # Mutually exclusive: current defaults to K=6 distinct-network support
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

        # Per-task scheduling requests, not the backend's total memory budget.
        # The CLI's --jobs/--memory configure Sprocket separately, not WDL inputs.
        # On the trusted local host backend these are not hard process RAM limits.
        String prepare_memory = "2 GiB"
        String inference_memory = "4 GiB"
        String consensus_memory = "8 GiB"
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
    }

    # Validate and stage LF-normalized copies without modifying original inputs.
    # parameters.json binds resolved options to input and runtime/config hashes.
    call Prepare { input:
        exp_file = exp_file,
        probe_file = probe_file,
        aracne_config_dir = aracne_config_dir,
        options = options,
        task_memory = prepare_memory,
    }

    # Independent seeded reconstructions consume the same prepared inputs.
    # Sprocket's CPU, memory and scatter limits determine actual concurrency.
    scatter (index in range(bootstrap_num)) {
        call Infer { input:
            exp_file = Prepare.expression,
            probe_file = Prepare.hubs,
            aracne_config_dir = aracne_config_dir,
            settings = Prepare.parameters,
            seed = index + 1,
            task_memory = inference_memory,
        }
    }

    # Wait for every reconstruction. The adapter checks the complete ordered
    # seed set and provenance before selecting consensus edges, then enhances
    # them using the full staged expression matrix rather than one resample.
    call Consensus { input:
        exp_file = Prepare.expression,
        aracne_config_dir = aracne_config_dir,
        settings = Prepare.parameters,
        adjacencies = Infer.adjacency,
        run_provenance = Infer.provenance,
        task_memory = consensus_memory,
    }

    output {
        # Historical name retained for compatibility: out_dir is an enhanced
        # network FILE, not a directory. The CLI publishes these four final files.
        File out_dir = Consensus.enhanced_network
        File consensus_network = Consensus.network
        File bootstrap_info = Consensus.bootstrap_info
        File parameter_info = Consensus.parameter_info
        # Per-seed artifacts and prepared settings remain available for auditing.
        Array[File] adjacencies = Infer.adjacency
        Array[File] run_provenance = Infer.provenance
        Array[File] native_stdout = Infer.native_stdout
        Array[File] native_stderr = Infer.native_stderr
        File validation_report = Prepare.validation_report
        File parameters = Prepare.parameters
    }

    meta {
        description: "Default SJARACNe WDL workflow, with explicit current and legacy scientific profiles."
    }
}

# Runs once before the scatter: validate IDs/options, stage safe copies and
# fingerprint the installed Python sources, native executable and model/config.
task Prepare {
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
}

# Runs once per seed: recheck prepared hashes, invoke the native executable and
# retain stdout/stderr plus provenance linking the adjacency to its inputs.
task Infer {
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
}

# Runs after all seeds: accept only matching successful adjacencies/provenance,
# stage adjacency files alone and call the shared consensus/enhancement functions.
task Consensus {
    input {
        File exp_file
        Directory? aracne_config_dir
        File settings
        Array[File] adjacencies
        Array[File] run_provenance
        String task_memory = "8 GiB"
    }

    env File request = write_json(object {
        exp_file: exp_file,
        aracne_config_dir: aracne_config_dir,
        settings: settings,
        adjacencies: adjacencies,
        run_provenance: run_provenance,
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
    }

    requirements {
        cpu: 1
        memory: task_memory
    }
}
