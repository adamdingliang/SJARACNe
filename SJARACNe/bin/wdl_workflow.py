"""Task-boundary plumbing for the standalone WDL workflow, not an inference engine.

Scientific work stays in the existing native program, input QC and consensus
module. Requests arrive as JSON, never as shell-interpolated user arguments.
"""
import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from .QC_input import validate_inputs
from ..workflow_options import OPTION_NAMES, integer, resolve_options


PACKAGE = Path(__file__).resolve().parents[1]
SCHEMA = "sjaracne-wdl-parameters-v2"


def sha256(path):
    """Hash file bytes without loading the whole task artifact into memory.

    Args:
        path (str or Path): File to fingerprint.

    Returns:
        str: Hexadecimal SHA-256 digest of the file contents.
    """
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    """Read a task-boundary JSON object, rejecting non-finite constants.

    Args:
        path (str or Path): UTF-8 JSON request, settings or provenance file.

    Returns:
        dict: Decoded JSON object.

    Raises:
        ValueError: The JSON is malformed, contains NaN/Infinity constants or
            has a non-object root.
    """
    def reject_constant(value):
        raise ValueError("Non-finite JSON number: " + value)
    with Path(path).open(encoding="utf-8") as handle:
        value = json.load(handle, parse_constant=reject_constant)
    if not isinstance(value, dict):
        raise ValueError("A JSON object is required")
    return value


def write_json(path, value):
    """Create a JSON artifact without overwriting an existing task output.

    Args:
        path (str or Path): Destination file in the task's working directory.
        value (object): JSON-serializable data, excluding non-finite numbers.

    Raises:
        FileExistsError: The destination already exists.
        ValueError: The data contains non-finite numbers.
    """
    with Path(path).open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def validate_parameters(parameters, observation_count):
    """Resolve scientific options with the same profile rules as the CLI.

    Args:
        parameters (dict): Requested profile and scientific overrides.
        observation_count (int): Number of expression-matrix observations.

    Returns:
        dict: Validated options with shared defaults resolved.

    Raises:
        ValueError: An option is unsupported, conflicting or out of range.
    """
    return resolve_options(parameters, observation_count)


def normalize_lf(source, destination):
    """Stage an LF-normalized copy without changing source bytes or row order.

    Args:
        source (str or Path): Input file whose CRLF and lone CR become LF.
        destination (str or Path): New output file; no decoding is performed.

    Raises:
        FileExistsError: The destination already exists.
    """
    with Path(source).open("rb") as source_handle, Path(destination).open("xb") as output:
        pending = b""
        for block in iter(lambda: source_handle.read(1024 * 1024), b""):
            block = pending + block
            # A trailing CR may belong to a CRLF pair split across chunks.
            pending = b"\r" if block.endswith(b"\r") else b""
            if pending:
                block = block[:-1]
            output.write(block.replace(b"\r\n", b"\n").replace(b"\r", b"\n"))
        output.write(pending.replace(b"\r", b"\n"))


def runtime_fingerprint(config_dir=None):
    """Fingerprint package code, native inference and selected calibration files.

    Args:
        config_dir (str or Path, optional): Native configuration directory.
            Defaults to the installed package's configuration directory.

    Returns:
        dict: Hashes of package-root/bin Python sources, the native executable,
            two native configuration files and the bundled PCHIP model data.
            Interpreter and dependency versions are not included.

    Raises:
        ValueError: The installed native executable is missing or lacks execute
            permission.
    """
    binary = PACKAGE / "bin/sjaracne.exe"
    if not binary.is_file() or not os.access(str(binary), os.X_OK):
        raise ValueError("Build the native executable in this installed SJARACNe package first")
    config = Path(config_dir) if config_dir is not None else PACKAGE / "config"
    paths = list(PACKAGE.glob("*.py")) + list((PACKAGE / "bin").glob("*.py"))
    return {
        "python_runtime": {str(path.relative_to(PACKAGE)): sha256(path) for path in sorted(paths)},
        "native_sha256": sha256(binary),
        "config_sha256": {name: sha256(config / name)
                          for name in ("config_threshold.txt", "config_kernel.txt")},
        "pchip_data_sha256": sha256(PACKAGE / "config/apmi_pchip.json"),
    }


def prepare(request):
    """Validate inputs and write normalized files plus resolved task settings.

    Writes the QC report, ``expression.exp``, ``hubs.txt`` and
    ``parameters.json`` in the task's working directory. Settings bind the
    resolved scientific profile to runtime hashes and both original and staged
    input hashes; source inputs are not modified.

    Args:
        request (dict): Original ``exp_file`` and ``probe_file`` paths,
            ``parameters`` overrides and optional ``aracne_config_dir``.

    Raises:
        ValueError: Scientific options or the installed runtime are invalid.
        SystemExit: Shared input QC rejects the expression matrix or hub list.
        FileExistsError: A normalized input or settings output already exists.
    """
    expression, hubs = Path(request["exp_file"]), Path(request["probe_file"])
    with expression.open(encoding="utf-8-sig") as handle:
        observations = len(handle.readline().rstrip("\r\n").split("\t")) - 2
    parameters = validate_parameters(request["parameters"], observations)
    fingerprint = runtime_fingerprint(request.get("aracne_config_dir"))
    validate_inputs(str(expression), str(hubs), "hub_overlap_validation.txt")
    normalize_lf(expression, "expression.exp")
    normalize_lf(hubs, "hubs.txt")
    write_json("parameters.json", {
        "schema": SCHEMA, "parameters": parameters, "observation_count": observations,
        "runtime": fingerprint,
        "input_sha256": {"exp_file": sha256("expression.exp"), "probe_file": sha256("hubs.txt")},
        "original_input_sha256": {"exp_file": sha256(expression), "probe_file": sha256(hubs)},
    })


def checked_settings(request):
    """Reject settings or localized inputs that differ from task preflight.

    Args:
        request (dict): ``settings`` JSON path, optional ``aracne_config_dir``
            and any staged ``exp_file``/``probe_file`` paths to verify.

    Returns:
        dict: Prepared settings after schema, runtime and supplied-input checks.

    Raises:
        ValueError: The schema is unsupported, the runtime/configuration has
            changed or a supplied input differs from its staged hash.
    """
    settings = read_json(request["settings"])
    if settings.get("schema") != SCHEMA:
        raise ValueError("Unsupported prepared parameter schema")
    if settings["runtime"] != runtime_fingerprint(request.get("aracne_config_dir")):
        raise ValueError("SJARACNe runtime/configuration differs from the preflight task")
    for key in ("exp_file", "probe_file"):
        if key in request and sha256(request[key]) != settings["input_sha256"][key]:
            raise ValueError("Prepared input changed: " + key)
    return settings


def adjacency_name(seed):
    """Return the shared adjacency filename for a positive integer seed.

    Args:
        seed (int): Native random seed, also used as the reconstruction index.

    Returns:
        str: ``TF_run_NNN.adj`` filename, padded to at least three digits.

    Raises:
        ValueError: The seed is outside the shared integer range.
    """
    return "TF_run_{:03d}.adj".format(integer(seed, "seed"))


def native_command(request, parameters):
    """Build literal native argv from resolved options and a bounded seed.

    The same hub file supplies reconstruction sources (``-s``) and regulator
    annotations (``-l``). Explicit seed and output naming bind each scattered
    task to its reconstruction index; argv uses resolved profile settings,
    without shell interpolation.

    Args:
        request (dict): Staged input paths, ``seed`` and optional native config.
        parameters (dict): Scientific options resolved during preflight.

    Returns:
        list: Native executable and arguments for one reconstruction.

    Raises:
        ValueError: The seed is not within ``1..bootstrap_num``.
    """
    seed = integer(request["seed"], "seed", maximum=parameters["bootstrap_num"])
    config = request.get("aracne_config_dir") or str(PACKAGE / "config")
    command = [str(PACKAGE / "bin/sjaracne.exe"),
               "-i", request["exp_file"], "-l", request["probe_file"], "-s", request["probe_file"],
               "-p", str(parameters["p_value_bootstrap"]), "-a", "adaptive_partitioning", "-H", config,
               "-N", str(parameters["depth"]), "-o", adjacency_name(seed), "-S", str(seed)]
    if parameters["sampling_method"] == "with-replacement":
        command += ["-r", "1"]
    else:
        command += ["-u", parameters["subsample_spec"]]
    if parameters["apmi_null_model"] == "pchip":
        command += ["-M", "pchip"]
    if parameters["dynamic_dpi"]:
        command += ["--dynamic-dpi"]
        for name, flag in (("dpi_href", "--dpi-href"), ("dpi_slope", "--dpi-slope"), ("dpi_max", "--dpi-max")):
            if name in parameters:
                command += [flag, str(parameters[name])]
    else:
        command += ["-e", str(parameters["tolerance"])]
    return command


def infer(request):
    """Run one seeded reconstruction and retain native logs and provenance.

    Creates the adjacency, ``native.stdout.txt``, ``native.stderr.txt`` and
    ``run_provenance.json`` in the task directory. After the native process
    returns, provenance records argv, exit status and settings/runtime/input
    hashes; successful output also receives an adjacency hash.

    Args:
        request (dict): Prepared settings and input paths, ``seed`` and optional
            ``aracne_config_dir``.

    Raises:
        ValueError: Preflight checks fail, the seed is invalid or the adjacency
            destination already exists.
        FileExistsError: A native log or provenance file already exists.
        RuntimeError: Native inference returns a nonzero exit status.
    """
    settings = checked_settings(request)
    command = native_command(request, settings["parameters"])
    output = Path(adjacency_name(request["seed"]))
    if output.exists():
        raise ValueError("Refusing to overwrite an adjacency file: " + str(output))
    with Path("native.stdout.txt").open("x") as stdout, Path("native.stderr.txt").open("x") as stderr:
        result = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
    provenance = {"seed": request["seed"], "argv": command, "returncode": result.returncode,
                  "runtime": settings["runtime"], "input_sha256": settings["input_sha256"],
                  "prepared_settings_sha256": sha256(request["settings"])}
    if result.returncode == 0:
        provenance["adjacency_sha256"] = sha256(output)
    write_json("run_provenance.json", provenance)
    if result.returncode:
        print(Path("native.stderr.txt").read_text(), file=sys.stderr)
        raise RuntimeError("Native inference failed; see native.stdout.txt/native.stderr.txt")


def consensus(request):
    """Verify reconstruction provenance before staging shared consensus inputs.

    Requires one adjacency/provenance pair per seed, ordered from 1 through
    ``bootstrap_num``. Copies only validated adjacencies into a new staging
    directory, then writes shared consensus and enhancement outputs to
    ``results`` using the prepared profile's recurrence or legacy cutoff.

    Args:
        request (dict): Prepared ``settings``, staged ``exp_file``, ordered
            ``adjacencies`` and ``run_provenance`` lists, and optional config.

    Raises:
        ValueError: Settings/input checks fail, the reconstruction count is
            wrong or an adjacency does not match its successful provenance.
        FileExistsError: The adjacency staging directory already exists.
    """
    settings = checked_settings(request)
    parameters = settings["parameters"]
    adjacencies, provenance_files = request["adjacencies"], request["run_provenance"]
    n = parameters["bootstrap_num"]
    settings_sha256 = sha256(request["settings"])
    if len(adjacencies) != n or len(provenance_files) != n:
        raise ValueError("Expected exactly bootstrap_num adjacencies and provenance records")
    # Reject incomplete or foreign runs before staging any consensus inputs.
    for seed, (adjacency, record) in enumerate(zip(adjacencies, provenance_files), 1):
        provenance = read_json(record)
        if (Path(adjacency).name != adjacency_name(seed) or provenance.get("seed") != seed
                or provenance.get("returncode") != 0 or provenance.get("runtime") != settings["runtime"]
                or provenance.get("input_sha256") != settings["input_sha256"]
                or provenance.get("prepared_settings_sha256") != settings_sha256
                or provenance.get("adjacency_sha256") != sha256(adjacency)):
            raise ValueError("Mismatched adjacency/provenance for seed {}".format(seed))
    staging = Path("adjacencies")
    staging.mkdir(exist_ok=False)
    for seed, adjacency in enumerate(adjacencies, 1):
        shutil.copyfile(adjacency, staging / adjacency_name(seed))
    # No logs, provenance or parameter files enter this directory. The shared
    # scientific function reads every entry and counts header-only runs too.
    module = importlib.import_module("SJARACNe.bin.create_consensus_network")
    selection = ({"min_recurrence": parameters["min_recurrence"]}
                 if parameters["p_value_consensus"] is None
                 else {"p_value": parameters["p_value_consensus"]})
    network = module.create_consensus_network(str(staging), out_dir="results", **selection)
    module.create_enhanced_consensus_network(request["exp_file"], network, "results")


def main():
    """Dispatch an engine-localized JSON request to one WDL task stage."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "infer", "consensus"))
    parser.add_argument("request", help="Engine-localized JSON request file")
    args = parser.parse_args()
    {"prepare": prepare, "infer": infer, "consensus": consensus}[args.stage](read_json(args.request))


if __name__ == "__main__":
    main()
