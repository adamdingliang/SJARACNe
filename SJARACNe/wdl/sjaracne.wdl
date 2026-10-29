version 1.3

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
        File exp_file
        File probe_file
        Directory? aracne_config_dir
        String profile = "current"
        String? network_type
        String? sampling_method
        Int bootstrap_num = 100
        String? subsample_spec
        Int depth = 40
        Float? p_value_bootstrap
        String? apmi_null_model
        Int? min_recurrence
        Float? p_value_consensus
        Float? tolerance
        Boolean? dynamic_dpi
        Int? dpi_href
        Float? dpi_slope
        Float? dpi_max
        String prepare_memory = "2 GiB"
        String inference_memory = "4 GiB"
        String consensus_memory = "8 GiB"
    }

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

    call Prepare { input:
        exp_file = exp_file,
        probe_file = probe_file,
        aracne_config_dir = aracne_config_dir,
        options = options,
        task_memory = prepare_memory,
    }

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

    call Consensus { input:
        exp_file = Prepare.expression,
        aracne_config_dir = aracne_config_dir,
        settings = Prepare.parameters,
        adjacencies = Infer.adjacency,
        run_provenance = Infer.provenance,
        task_memory = consensus_memory,
    }

    output {
        File out_dir = Consensus.enhanced_network
        File consensus_network = Consensus.network
        File bootstrap_info = Consensus.bootstrap_info
        File parameter_info = Consensus.parameter_info
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

task Prepare {
    input {
        File exp_file
        File probe_file
        Directory? aracne_config_dir
        WorkflowOptions options
        String task_memory = "2 GiB"
    }

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
