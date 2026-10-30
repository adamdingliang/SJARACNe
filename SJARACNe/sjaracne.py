#!/usr/bin/env python3
"""SJARACNe workflow CLI: WDL/Sprocket by default, explicit CWL compatibility."""
import argparse
from contextvars import ContextVar
import json
import logging
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import threading

from SJARACNe.bin.QC_input import validate_inputs
from SJARACNe.bin.wdl_workflow import runtime_fingerprint, sha256
from SJARACNe.workflow_options import resolve_options

PACKAGE = Path(__file__).resolve().parent
DEFAULT_MIN_RECURRENCE = 6
SPROCKET_VERSION = "0.32.0"
_EXECUTION_CONTEXT = ContextVar("sjaracne_cli_execution", default=None)
SCIENTIFIC_FIELDS = (
    "profile", "network_type", "bootstrap_num", "subsample_spec", "depth",
    "p_value_bootstrap", "apmi_null_model", "min_recurrence", "p_value_consensus",
    "tolerance", "dynamic_dpi", "dpi_href", "dpi_slope", "dpi_max", "sampling_method",
)
FINAL_OUTPUTS = {
    "out_dir": "consensus_network_ncol_.txt", "consensus_network": "consensus_network_3col_.txt",
    "bootstrap_info": "bootstrap_info_.txt", "parameter_info": "parameter_info_.txt",
}


def dpi_reference_hubs(value):
    """Parse the positive, signed-32-bit reference hub count.

    Args:
        value: Integer text supplied by argparse.

    Returns:
        int: Reference hub count used by dynamic DPI.

    Raises:
        ValueError: The value cannot be converted to an integer.
        argparse.ArgumentTypeError: The count is outside the supported range.
    """
    parsed = int(value)
    if not 1 <= parsed <= 2147483647:
        raise argparse.ArgumentTypeError("DPI reference hub count must be within [1,2147483647]")
    return parsed


def dpi_positive_slope(value):
    """Parse a finite, positive dynamic-DPI slope.

    Args:
        value: Numeric text supplied by argparse.

    Returns:
        float: Validated slope.

    Raises:
        ValueError: The value cannot be converted to a float.
        argparse.ArgumentTypeError: The slope is nonfinite or nonpositive.
    """
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0:
        raise argparse.ArgumentTypeError("DPI slope must be finite and positive")
    return parsed


def dpi_maximum(value):
    """Parse the dynamic-DPI cap in the half-open interval [0, 1).

    Args:
        value: Numeric text supplied by argparse.

    Returns:
        float: Validated DPI cap.

    Raises:
        ValueError: The value cannot be converted to a float.
        argparse.ArgumentTypeError: The cap is nonfinite or outside [0, 1).
    """
    parsed = float(value)
    if not math.isfinite(parsed) or not 0 <= parsed < 1:
        raise argparse.ArgumentTypeError("DPI maximum must be finite and within [0,1)")
    return parsed


def sampling_fraction(value):
    """Parse a subsampling fraction in the interval (0, 1].

    Args:
        value: Numeric text supplied by argparse.

    Returns:
        float: Validated fraction of observations to sample.

    Raises:
        ValueError: The value cannot be converted to a float.
        argparse.ArgumentTypeError: The fraction is nonfinite or out of range.
    """
    parsed = float(value)
    if not math.isfinite(parsed) or not 0 < parsed <= 1:
        raise argparse.ArgumentTypeError("subsample fraction must be within (0,1]")
    return parsed


def sampling_size(value):
    """Parse an exact subsample size of at least two observations.

    Args:
        value: Integer text supplied by argparse.

    Returns:
        int: Requested size; dataset-dependent bounds are checked later.

    Raises:
        ValueError: The value cannot be converted to an integer.
        argparse.ArgumentTypeError: The size is less than two.
    """
    parsed = int(value)
    if parsed < 2:
        raise argparse.ArgumentTypeError("subsample size must be at least 2")
    return parsed


def minimum_recurrence(value):
    """Parse a positive minimum number of distinct supporting networks.

    Args:
        value: Integer text supplied by argparse.

    Returns:
        int: Recurrence threshold; its upper bound is resolved with the profile.

    Raises:
        ValueError: The value cannot be converted to an integer.
        argparse.ArgumentTypeError: The count is less than one.
    """
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("minimum recurrence must be at least 1")
    return parsed


def bootstrap_count(value):
    """Parse a positive integer for count-valued CLI options.

    Args:
        value: Integer text supplied by argparse for networks, depth, or jobs.

    Returns:
        int: Validated positive count.

    Raises:
        ValueError: The value cannot be converted to an integer.
        argparse.ArgumentTypeError: The count is less than one.
    """
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("bootstrap count must be at least 1")
    return parsed


def consensus_probability(value):
    """Parse a legacy consensus probability in the interval (0, 1].

    Args:
        value: Numeric text supplied by argparse.

    Returns:
        float: Validated probability, not a recurrence threshold.

    Raises:
        ValueError: The value cannot be converted to a float.
        argparse.ArgumentTypeError: The probability is nonfinite or out of range.
    """
    parsed = float(value)
    if not math.isfinite(parsed) or not 0 < parsed <= 1:
        raise argparse.ArgumentTypeError("consensus p-value must be within (0,1]")
    return parsed


def build_parser():
    """Build the CLI without resolving dataset-dependent scientific options.

    Scientific profiles, workflow engines, and scheduling resources are separate
    choices. Additional cross-option and engine constraints are checked in main.

    Returns:
        argparse.ArgumentParser: Parser for the local command and its run alias.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-e", "--exp-file", required=True, help="Expression matrix (ID, symbol, observations)")
    common.add_argument("-g", "--hub-genes", required=True, help="Source IDs matching the matrix's first column")
    common.add_argument("-o", "--output-dir", required=True, help="New or empty final output directory")
    common.add_argument("-tmp", "--tmpdir-prefix", help="Scratch parent; only a new child is used, existing contents are not deleted")
    common.add_argument("--profile", choices=("current", "legacy"), default="current")
    common.add_argument("--network-type", type=str.upper, choices=("TF", "SIG"), help="Required for current profile")
    common.add_argument("--engine", choices=("sprocket", "cwl"), help="Local default sprocket; cwl retains compatibility")
    common.add_argument("-M", "--apmi-null-model", choices=("pchip", "affine"), help="Current default PCHIP")
    common.add_argument("-pb", "--p-value-bootstrap", type=float, help="Current TF=.001, SIG=.0005; legacy=1e-7")
    common.add_argument("-n", "--bootstrap-num", type=bootstrap_count, help="Number of networks (default100)")
    common.add_argument("-d", "--depth", type=bootstrap_count, help="AP-MI Npar limit (default40)")
    common.add_argument("-c", "--config-dir", help="ARACNe configuration directory; default bundled")
    consensus = common.add_mutually_exclusive_group()
    consensus.add_argument("-k", "--min-recurrence", type=minimum_recurrence, help="Current default K=6")
    consensus.add_argument("-pc", "--p-value-consensus", type=consensus_probability, help="Legacy normal-approximation consensus probability")
    sampling = common.add_mutually_exclusive_group()
    sampling.add_argument("-sf", "--subsample-fraction", type=sampling_fraction, help="Fraction without replacement; current default.8")
    sampling.add_argument("-sm", "--subsample-size", type=sampling_size, help="Exact size without replacement")
    dpi = common.add_mutually_exclusive_group()
    dpi.add_argument("--dynamic-dpi", dest="dynamic_dpi", action="store_true", help="Current default: hub-count-dependent DPI")
    dpi.add_argument("--fixed-dpi", dest="dynamic_dpi", action="store_false", help="Fixed DPI; omitted tolerance means strict zero")
    common.set_defaults(dynamic_dpi=None)
    common.add_argument("--tolerance", type=float, help="Fixed DPI epsilon; implies --fixed-dpi")
    common.add_argument("--dpi-href", type=dpi_reference_hubs)
    common.add_argument("--dpi-slope", type=dpi_positive_slope)
    common.add_argument("--dpi-max", type=dpi_maximum)
    common.add_argument("--sprocket", default=os.environ.get("SJARACNE_SPROCKET", "sprocket"), help="Sprocket0.32.0 executable, or SJARACNE_SPROCKET")
    common.add_argument("--sprocket-config", help="Explicit trusted backend configuration")
    common.add_argument("--jobs", type=bootstrap_count, help="Sprocket-only local task budget (default4); custom config owns concurrency")
    common.add_argument("--memory", default="16 GiB", help="Sprocket-only local scheduling memory budget (default16 GiB)")
    common.add_argument("--prepare-memory", default="2 GiB", help="Sprocket-only input-preparation memory request")
    common.add_argument("--inference-memory", default="4 GiB", help="Sprocket-only per-reconstruction memory request")
    common.add_argument("--consensus-memory", default="8 GiB", help="Sprocket-only consensus memory request")
    commands = parser.add_subparsers(dest="subcommand")
    local = commands.add_parser("local", aliases=["run"], parents=[common], help="WDL/Sprocket local workflow")
    local.add_argument("-s", "--serial", action="store_true", help="One local task at a time")
    return parser


def require_tool(name):
    """Find the requested executable without selecting a fallback engine.

    Args:
        name: Executable name or path accepted by shutil.which.

    Returns:
        str: Located executable path.

    Raises:
        ValueError: The executable cannot be found.
    """
    found = shutil.which(str(name))
    if not found:
        raise ValueError("Required executable not found: {}. Install it explicitly; no engine fallback is performed.".format(name))
    return found


def check_sprocket_version(executable):
    """Require the tested Sprocket version from its combined version output.

    Args:
        executable: Path to the Sprocket executable.

    Returns:
        str: Stripped stdout and stderr from the version command.

    Raises:
        ValueError: The reported version does not match the pinned version.
        OSError: The executable cannot be started.
        subprocess.CalledProcessError: The version command exits unsuccessfully.
    """
    result = subprocess.run([str(executable), "--version"], check=True, capture_output=True, text=True)
    version = (result.stdout + result.stderr).strip()
    if not re.search(r"^sprocket v?" + re.escape(SPROCKET_VERSION) + r"(?:\s|$)", version):
        raise ValueError("This workflow is tested with Sprocket {}; found {}".format(SPROCKET_VERSION, version))
    return version


def memory_bytes(value):
    """Parse a positive memory quantity for scheduling-budget comparisons.

    Args:
        value: Quantity with a supported, case-sensitive binary or decimal unit.

    Returns:
        float: Quantity converted to bytes, without checking available host RAM.

    Raises:
        ValueError: The quantity has invalid syntax or is not positive.
    """
    match = re.fullmatch(r"\s*([0-9]+(?:\.[0-9]+)?)\s*(B|KiB|MiB|GiB|TiB|KB|MB|GB|TB)\s*", value)
    if not match or float(match.group(1)) <= 0:
        raise ValueError("Memory must be a positive quantity such as '4 GiB'")
    units = {"B": 1, "KiB": 1024, "MiB": 1024**2, "GiB": 1024**3, "TiB": 1024**4,
             "KB": 1000, "MB": 1000**2, "GB": 1000**3, "TB": 1000**4}
    return float(match.group(1)) * units[match.group(2)]


def local_config(jobs, memory):
    """Render a trusted-host Sprocket backend and scatter configuration.

    Jobs and total memory belong to engine configuration, not WDL inputs.
    Per-task memory requests are supplied separately; neither these reservations
    nor the backend budget impose a hard process-memory limit.

    Args:
        jobs: Positive local CPU and scatter-slot budget, already validated.
        memory: Validated total scheduling-memory quantity, including its unit.

    Returns:
        str: TOML configuration with retries and task caching disabled.
    """
    return ('# Generated trusted-host configuration; no container isolation.\n'
            '[run.backends.default]\ntype = "local"\ncpu = {}\nmemory = {}\n\n'
            '[run.workflow.scatter]\nconcurrency = {}\n\n'
            '[run.task]\nretries = 0\ncache = "off"\n').format(jobs, json.dumps(memory), jobs)


def execution_env():
    """Prepare engine subprocesses to use this package and numerical thread caps.

    The parent environment is copied, not modified. Package binaries and the
    active Python directory lead PATH, and the package parent leads PYTHONPATH.

    Returns:
        dict: Environment requesting one thread from OpenBLAS, OpenMP, and MKL.
    """
    environment = os.environ.copy()
    environment["PATH"] = os.pathsep.join((str(PACKAGE / "bin"), str(Path(sys.executable).parent), environment.get("PATH", "")))
    environment["PYTHONPATH"] = os.pathsep.join(filter(None, (str(PACKAGE.parent), environment.get("PYTHONPATH", ""))))
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        environment[name] = "1"
    return environment


def run_shell_command_call(cmd):
    """Run a command without a shell, retaining configured engine stream logs.

    An execution context supplies the environment and log paths. Without log
    paths, the child inherits the normal streams. With logs, both pipes are
    drained concurrently; console failures do not discard retained output.

    Args:
        cmd: Argument sequence, or command text split into arguments by shlex.

    Returns:
        int: Zero after a successful command and successful stream retention.

    Raises:
        OSError: The child cannot start, or a configured stream log cannot be kept.
        subprocess.CalledProcessError: The command exits with a nonzero status.
    """
    argv = shlex.split(cmd) if isinstance(cmd, str) else [str(value) for value in cmd]
    context = _EXECUTION_CONTEXT.get() or {}
    if context.get("stdout") is None:
        return subprocess.check_call(argv, env=context.get("env"))
    with Path(context["stdout"]).open("w", encoding="utf-8") as stdout_log, Path(context["stderr"]).open("w", encoding="utf-8") as stderr_log:
        with subprocess.Popen(argv, env=context.get("env"), stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, text=True, bufsize=1) as process:
            errors = []
            def tee(source, log, destination):
                """Drain one child pipe, keeping its log even if the console fails.

                Args:
                    source: Readable child-process text stream, closed on exit.
                    log: Retained UTF-8 log stream.
                    destination: Console text stream, or None to disable echoing.
                """
                try:
                    for line in source:
                        log.write(line)
                        log.flush()
                        if destination is not None:
                            try:
                                destination.write(line)
                                destination.flush()
                            except (OSError, UnicodeError):
                                # A closed console pipe must not stop draining
                                # the child or discard the retained UTF-8 log.
                                destination = None
                except (OSError, UnicodeError) as error:
                    errors.append(error)
                    process.terminate()
                finally:
                    source.close()
            # Drain stdout and stderr together so a full pipe cannot stall the child.
            threads = [threading.Thread(target=tee, args=values) for values in (
                (process.stdout, stdout_log, sys.stdout), (process.stderr, stderr_log, sys.stderr))]
            for thread in threads:
                thread.start()
            returncode = process.wait()
            for thread in threads:
                thread.join()
    if errors:
        raise OSError("Could not retain workflow stream log: " + str(errors[0]))
    if returncode:
        raise subprocess.CalledProcessError(returncode, argv)
    return returncode


def ensure_empty_output(path):
    """Reject output paths that would reuse existing content or a symlink.

    Args:
        path: pathlib.Path naming a new or empty output directory.

    Raises:
        ValueError: The path is a symlink, a non-directory, or a nonempty directory.
        OSError: Existing directory contents cannot be inspected.
    """
    if path.is_symlink() or (path.exists() and (not path.is_dir() or any(path.iterdir()))):
        raise ValueError("Output must be a new or empty directory; refusing to overwrite: " + str(path))


def normalized_input(path, name):
    """Resolve an input path and require an existing regular file.

    Args:
        path: Input path, optionally containing a user-home abbreviation.
        name: Human-readable input label used in validation errors.

    Returns:
        pathlib.Path: Resolved absolute input path.

    Raises:
        ValueError: The resolved input is not an existing file.
    """
    candidate = Path(path).expanduser().resolve()
    if not candidate.is_file():
        raise ValueError(name + " does not exist: " + str(candidate))
    return candidate


def science_request(args):
    """Extract explicit scientific options independently of engine resources.

    None-valued options remain absent for resolve_options to fill from the profile;
    parser-supplied defaults such as the profile itself are retained.
    An explicit tolerance selects fixed DPI unless dynamic DPI was also requested.

    Args:
        args: Parsed CLI namespace, including sampling and DPI options.

    Returns:
        dict: Scientific request with a normalized subsampling specification.

    Raises:
        ValueError: Explicit fixed tolerance conflicts with explicit dynamic DPI.
    """
    values = {name: getattr(args, name, None) for name in SCIENTIFIC_FIELDS}
    if args.subsample_size is not None:
        values["subsample_spec"] = str(args.subsample_size)
    elif args.subsample_fraction is not None:
        values["subsample_spec"] = "{:.12g}%".format(args.subsample_fraction * 100)
    if args.tolerance is not None:
        if args.dynamic_dpi is True:
            raise ValueError("--dynamic-dpi conflicts with an explicit fixed tolerance")
        values["dynamic_dpi"] = False
    return {key: value for key, value in values.items() if value is not None}


def make_cwl_job(args, parameters, expression, hubs, config):
    """Map resolved scientific settings into the retained CWL input contract.

    Engine selection does not select the scientific profile. CWL's sampling and
    null-model encodings are translated here; Sprocket budgets are not CWL inputs.

    Args:
        args: Parsed CLI namespace providing the final output directory name.
        parameters: Validated scientific options returned by resolve_options.
        expression: Resolved expression-matrix path.
        hubs: Resolved source-ID file path.
        config: Resolved ARACNe configuration directory.

    Returns:
        dict: CWL inputs with File and Directory objects and scientific settings.
    """
    job = {"exp_file": {"class": "File", "path": str(expression)},
           "probe_file": {"class": "File", "path": str(hubs)},
           "aracne_config_dir": {"class": "Directory", "path": str(config)},
           "final_out_dir_name": Path(args.output_dir).expanduser().absolute().name,
           "bootstrap_num": parameters["bootstrap_num"], "depth": parameters["depth"],
           "p_value_bootstrap": parameters["p_value_bootstrap"],
           "apmi_null_model": "pchip" if parameters["apmi_null_model"] == "pchip" else None,
           "dynamic_dpi": parameters["dynamic_dpi"], "subsample_spec": parameters.get("subsample_spec"),
           "sample_number": 1 if parameters["sampling_method"] == "with-replacement" else None}
    for name in ("min_recurrence", "p_value_consensus", "tolerance", "dpi_href", "dpi_slope", "dpi_max"):
        if parameters.get(name) is not None:
            job[name] = parameters[name]
    return job


def write_json(path, data):
    """Write deterministic, UTF-8 JSON without nonstandard numeric values.

    Args:
        path: Destination file; parent directories must already exist.
        data: JSON-serializable value to write with sorted keys and a final newline.

    Raises:
        OSError: The destination cannot be written.
        TypeError: A value cannot be serialized as JSON.
        ValueError: A value is circular or contains NaN or infinity.
    """
    with Path(path).open("w", encoding="utf-8", newline="\n") as output:
        json.dump(data, output, sort_keys=True, indent=2, allow_nan=False)
        output.write("\n")


def path_within(path, root):
    """Check resolved path containment, including symlink resolution.

    Args:
        path: Candidate file or directory path.
        root: Expected containing directory.

    Returns:
        bool: Whether the candidate resolves inside or exactly to root; false for
            paths whose common root cannot be compared, such as different drives.
    """
    try:
        return os.path.commonpath((str(Path(path).resolve()), str(Path(root).resolve()))) == str(Path(root).resolve())
    except ValueError:
        return False


def publish_sprocket_outputs(engine_root, output, expected_runtime=None, expected_parameters=None):
    """Validate one completed execution and copy its final files without overwrite.

    Every required source must resolve within this execution. When supplied,
    expected runtime and scientific settings are compared with task-prepared
    provenance before publication. All sources are checked before copying, but
    publication is not transactional: a copy failure may leave partial evidence.

    Args:
        engine_root: Fresh Sprocket execution directory as a pathlib.Path.
        output: Final output directory as a pathlib.Path.
        expected_runtime: Optional CLI-preflight runtime fingerprint to verify.
        expected_parameters: Optional resolved scientific settings to verify.

    Returns:
        tuple: Engine outputs.json path and published-file paths with SHA256 hashes.

    Raises:
        ValueError: Completion, provenance, containment, or overwrite checks fail.
        OSError: Engine evidence cannot be read or final files cannot be copied.
    """
    runs = engine_root / "runs/SJARACNe"
    candidates = [path / "outputs.json" for path in runs.iterdir()
                  if path.is_dir() and not path.name.startswith("_") and (path / "outputs.json").is_file()] if runs.is_dir() else []
    if len(candidates) != 1:
        raise ValueError("Expected exactly one completed Sprocket workflow outputs.json")
    outputs_json = candidates[0]
    outputs = json.loads(outputs_json.read_text(encoding="utf-8"))
    if expected_runtime is not None or expected_parameters is not None:
        prepared_path = outputs.get("SJARACNe.parameters")
        if not isinstance(prepared_path, str) or not Path(prepared_path).is_file() or not path_within(prepared_path, engine_root):
            raise ValueError("Missing prepared-parameter provenance inside this execution")
        prepared = json.loads(Path(prepared_path).read_text(encoding="utf-8"))
        if expected_runtime is not None and prepared.get("runtime") != expected_runtime:
            raise ValueError("Executed SJARACNe runtime differs from the CLI preflight")
        if expected_parameters is not None and prepared.get("parameters") != expected_parameters:
            raise ValueError("Executed scientific parameters differ from the CLI request")
    sources = {}
    for key, filename in FINAL_OUTPUTS.items():
        value = outputs.get("SJARACNe." + key)
        if not isinstance(value, str):
            raise ValueError("Missing or invalid workflow file output: " + key)
        source = Path(value)
        if not source.is_file() or not path_within(source, engine_root):
            raise ValueError("Workflow output is absent or outside this fresh execution: " + value)
        destination = output / filename
        if destination.exists() or destination.is_symlink():
            raise ValueError("Refusing to overwrite published output: " + str(destination))
        sources[key] = (source, destination)
    published = {}
    # Exclusive creation protects each destination even after the earlier check.
    for key, (source, destination) in sources.items():
        with source.open("rb") as reader, destination.open("xb") as writer:
            shutil.copyfileobj(reader, writer)
        published[key] = {"path": str(destination), "sha256": sha256(destination)}
    write_json(output / "sprocket_outputs.json", outputs)
    return outputs_json, published


def validate_cwl_outputs(output):
    """Require and hash CWL's final files rather than trusting engine exit alone.

    Args:
        output: Directory where CWL was instructed to publish final files.

    Returns:
        dict: Required output keys mapped to their file paths and SHA256 hashes.

    Raises:
        ValueError: A required final file is absent.
        OSError: A final file cannot be read for hashing.
    """
    published = {}
    for key, filename in FINAL_OUTPUTS.items():
        path = Path(output) / filename
        if not path.is_file():
            raise ValueError("CWL did not publish the required final output: " + str(path))
        published[key] = {"path": str(path), "sha256": sha256(path)}
    return published


def main(argv=None):
    """Validate, execute, and publish one workflow while retaining run evidence.

    Preflight resolves scientific settings separately from engine configuration
    and checks inputs before creating final output artifacts. Once execution is
    prepared, failures retain the manifest and any engine logs or work files
    already created.
    Existing output content is never deliberately cleared for a retry.

    Args:
        argv: Optional argument sequence; None uses sys.argv without the program.

    Returns:
        int: Zero for a completed workflow or an invocation with no arguments.

    Raises:
        SystemExit: Argparse handles help, invalid input, or a recorded run failure.
        OSError: Output setup or failure-evidence writing cannot be completed.
    """
    parser = build_parser()
    arguments = list(sys.argv[1:] if argv is None else argv)
    if not arguments:
        parser.print_help()
        return 0
    if arguments[0] == "lsf":
        parser.error("The old 'lsf' command was removed. Use WDL/Sprocket with an explicit "
                     "--sprocket-config; scheduler execution remains unvalidated.")
    if arguments[0] not in ("local", "run", "-h", "--help"):
        arguments.insert(0, "local")
    args = parser.parse_args(arguments)
    logging.basicConfig(level=logging.INFO)
    try:
        engine = args.engine or "sprocket"
        # Reject incompatible engine/concurrency options instead of ignoring them.
        if engine == "cwl" and args.sprocket_config is not None:
            raise ValueError("--sprocket-config is only valid with --engine sprocket")
        if engine == "cwl" and args.jobs is not None:
            raise ValueError("--jobs is only supported by Sprocket; use --serial for serial CWL execution")
        if args.sprocket_config is not None and (args.jobs is not None or args.serial):
            raise ValueError("A supplied --sprocket-config owns concurrency; do not combine it with --jobs or --serial")
        expression = normalized_input(args.exp_file, "Expression matrix")
        hubs = normalized_input(args.hub_genes, "Hub file")
        config = Path(args.config_dir).expanduser().resolve() if args.config_dir else PACKAGE / "config"
        output = Path(args.output_dir).expanduser().absolute()
        ensure_empty_output(output)
        if args.tmpdir_prefix:
            scratch_parent = Path(args.tmpdir_prefix).expanduser().absolute()
            if scratch_parent.exists() and not scratch_parent.is_dir():
                raise ValueError("Scratch parent is not a directory: " + str(scratch_parent))
        with expression.open(encoding="utf-8-sig") as handle:
            observations = len(handle.readline().rstrip("\r\n").split("\t")) - 2
        parameters = resolve_options(science_request(args), observations)
        # A custom engine config owns its backend budget; task requests stay inputs.
        for value in (args.memory, args.prepare_memory, args.inference_memory, args.consensus_memory):
            memory_bytes(value)
        if engine == "sprocket" and args.sprocket_config is None:
            if any(memory_bytes(value) > memory_bytes(args.memory) for value in (
                    args.prepare_memory, args.inference_memory, args.consensus_memory)):
                raise ValueError("A task memory request exceeds the local memory scheduling budget")
        tool = require_tool(args.sprocket if engine == "sprocket" else "cwltool")
        version = check_sprocket_version(tool) if engine == "sprocket" else None
        if args.sprocket_config is not None:
            normalized_input(args.sprocket_config, "Sprocket configuration")
        # Capture runtime/config identity for evidence and Sprocket publication checks.
        runtime = runtime_fingerprint(str(config))
        with tempfile.TemporaryDirectory(prefix="sjaracne-preflight-") as temporary:
            validate_inputs(str(expression), str(hubs), str(Path(temporary) / "validation.txt"))
    except (ValueError, OSError, subprocess.SubprocessError) as error:
        parser.error(str(error))

    # Recheck after preflight without deleting anything created in the meantime.
    ensure_empty_output(output)
    output.mkdir(parents=True, exist_ok=True)
    manifest_path = output / "sjaracne_run.json"
    record = {"schema": "sjaracne-cli-run-v1", "state": "prepared", "engine": engine,
              "scientific_options": parameters, "runtime": runtime, "sprocket_version": version,
              "inputs": {"expression": str(expression), "hubs": str(hubs)}, "published_outputs": {}}
    environment = execution_env()
    try:
        job_file = output / "sjaracne_workflow.yml"
        if engine == "sprocket":
            if args.tmpdir_prefix:
                scratch = Path(args.tmpdir_prefix).expanduser().absolute()
                scratch.mkdir(parents=True, exist_ok=True)
                engine_root = Path(tempfile.mkdtemp(prefix="sjaracne-sprocket-", dir=str(scratch)))
            else:
                engine_root = output / ".sprocket"
            record["engine_output_dir"] = str(engine_root)
            job = {"SJARACNe.exp_file": str(expression), "SJARACNe.probe_file": str(hubs), "SJARACNe.aracne_config_dir": str(config)}
            job.update({"SJARACNe." + name: value for name, value in parameters.items() if name in SCIENTIFIC_FIELDS})
            # Only per-task requests enter WDL; --jobs/--memory configure the engine.
            for name in ("prepare_memory", "inference_memory", "consensus_memory"):
                job["SJARACNe." + name] = getattr(args, name)
            if args.sprocket_config:
                backend_config = Path(args.sprocket_config).expanduser().resolve()
            else:
                backend_config = output / "sprocket.toml"
                backend_config.write_text(local_config(1 if args.serial else (args.jobs or 4), args.memory), encoding="utf-8")
            command = [tool, "--skip-config-search", "--config", str(backend_config), "run",
                       str(PACKAGE / "wdl/sjaracne.wdl"), "@" + str(job_file), "--output-dir", str(engine_root),
                       "--no-call-cache", "--disable-retries", "--show-task-stderr"]
            record["engine_configuration_sha256"] = sha256(backend_config)
        else:
            job = make_cwl_job(args, parameters, expression, hubs, config)
            scratch = Path(args.tmpdir_prefix).expanduser().absolute() if args.tmpdir_prefix else output / ".cwl_tmp"
            scratch.mkdir(parents=True, exist_ok=True)
            command = [tool, "--preserve-environment", "PYTHONPATH", "--tmpdir-prefix", str(scratch) + os.sep, "--outdir", str(output)]
            if not args.serial:
                command += ["--parallel"]
            command += [str(PACKAGE / "cwl/sjaracne_workflow.cwl"), str(job_file)]
        write_json(job_file, job)
        record.update(state="running", argv=command, job_file=str(job_file))
        write_json(manifest_path, record)
        logging.info("Profile=%s engine=%s; evidence: %s", parameters["profile"], record["engine"], output)
        token = _EXECUTION_CONTEXT.set({"env": environment, "stdout": output / "engine.stdout.log", "stderr": output / "engine.stderr.log"})
        try:
            run_shell_command_call(command)
        finally:
            _EXECUTION_CONTEXT.reset(token)
        if engine == "sprocket":
            outputs_json, published = publish_sprocket_outputs(engine_root, output, runtime, parameters)
            record.update(engine_outputs_json=str(outputs_json), published_outputs=published)
        else:
            record["published_outputs"] = validate_cwl_outputs(output)
        record["state"] = "complete"
        write_json(manifest_path, record)
    except (ValueError, OSError, subprocess.SubprocessError) as error:
        # Keep failed-run artifacts for diagnosis; a retry needs a new/empty output.
        record.update(state="failed", error=str(error))
        write_json(manifest_path, record)
        parser.exit(1, "Workflow failed: {}\nEvidence retained: {}\n".format(error, output))
    logging.info("All done: %s", output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
