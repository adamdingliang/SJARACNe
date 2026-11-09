"""Task-boundary plumbing for the standalone WDL workflow, not an inference engine.

Scientific work stays in the existing native program, input QC and consensus
module. Requests arrive as JSON, never as shell-interpolated user arguments.
"""
import argparse
import csv
import math
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
SCHEMA = "sjaracne-wdl-parameters-v3"


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
    if "bootstrap_num" in request and integer(request["bootstrap_num"], "bootstrap_num") != parameters["bootstrap_num"]:
        raise ValueError("Scatter bootstrap_num differs from the scientific settings")
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
    if parameters["consensus_order"] == "consensus-first":
        command += ["--raw-mi"]
        return command
    if parameters["dynamic_dpi"]:
        command += ["--dynamic-dpi"]
        for name, flag in (("dpi_href", "--dpi-href"), ("dpi_slope", "--dpi-slope"), ("dpi_max", "--dpi-max"),
                           ("dpi_offset", "--dpi-offset")):
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
    ``bootstrap_num``. Consensus-first passes the existing localized files to
    the native accumulator and DPI in one invocation. DPI-first retains the
    historical Python selector. Only final networks are enhanced from the full
    expression matrix. Successful completion hashes all four retained outputs.

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
    validated_adjacencies = []
    for seed, (adjacency, record) in enumerate(zip(adjacencies, provenance_files), 1):
        provenance = read_json(record)
        if (Path(adjacency).name != adjacency_name(seed) or provenance.get("seed") != seed
                or provenance.get("returncode") != 0 or provenance.get("runtime") != settings["runtime"]
                or provenance.get("input_sha256") != settings["input_sha256"]
                or provenance.get("prepared_settings_sha256") != settings_sha256
                or provenance.get("adjacency_sha256") != sha256(adjacency)):
            raise ValueError("Mismatched adjacency/provenance for seed {}".format(seed))
        validated_adjacencies.append({"seed": seed, "sha256": provenance["adjacency_sha256"]})
    module = importlib.import_module("SJARACNe.bin.create_consensus_network")
    native_order = parameters["consensus_order"] == "consensus-first"
    provenance = {"schema": "sjaracne-consensus-v1", "runtime": settings["runtime"],
                  "parameters": parameters, "prepared_settings_sha256": settings_sha256,
                  "input_sha256": settings["input_sha256"],
                  "adjacencies": validated_adjacencies, "returncode": 0}
    if native_order:
        command = native_consensus_command(request, parameters)
        provenance["argv"] = command
        with Path("consensus.native.stdout.txt").open("x", encoding="utf-8") as stdout, \
                Path("consensus.native.stderr.txt").open("x", encoding="utf-8") as stderr:
            result = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
        provenance["returncode"] = result.returncode
        if result.returncode:
            write_json("consensus_provenance.json", provenance)
            raise RuntimeError("Native consensus/DPI failed; see consensus.native.stdout.txt/consensus.native.stderr.txt")
        network = export_native_consensus("consensus.adj", "results", parameters)
        provenance["native_adjacency_sha256"] = sha256("consensus.adj")
    else:
        staging = Path("adjacencies")
        staging.mkdir(exist_ok=False)
        for seed, adjacency in enumerate(adjacencies, 1):
            shutil.copyfile(adjacency, staging / adjacency_name(seed))
        selection = ({"min_recurrence": parameters["min_recurrence"]}
                     if parameters["p_value_consensus"] is None
                     else {"p_value": parameters["p_value_consensus"]})
        network = module.create_consensus_network(str(staging), out_dir="results", **selection)
        for name in ("consensus.native.stdout.txt", "consensus.native.stderr.txt"):
            Path(name).touch(exist_ok=False)
    if native_order:
        module.create_enhanced_consensus_network(request["exp_file"], network, "results", preserve_mi=True)
    else:
        module.create_enhanced_consensus_network(request["exp_file"], network, "results")
    provenance["output_sha256"] = {name: sha256(Path("results") / name) for name in (
        "consensus_network_3col_.txt", "consensus_network_ncol_.txt", "bootstrap_info_.txt", "parameter_info_.txt")}
    write_json("consensus_provenance.json", provenance)


def native_consensus_command(request, parameters):
    """Write a path-only manifest and build one native consensus/DPI invocation.

    The manifest preserves seed order and refers to already localized files;
    neither raw adjacencies nor candidate mean weights are copied or serialized.

    Args:
        request (dict): Localized expression, hubs, adjacency paths and config.
        parameters (dict): Resolved recurrence and DPI settings.

    Returns:
        list[str]: Literal native argv for consensus followed by DPI.

    Raises:
        ValueError: An adjacency path contains a line break.
        FileExistsError: The manifest already exists.
    """
    manifest = Path("consensus_inputs.txt")
    with manifest.open("x", encoding="utf-8", newline="\n") as output:
        for path in request["adjacencies"]:
            path = str(Path(path).resolve())
            if "\n" in path or "\r" in path:
                raise ValueError("Adjacency paths must not contain line breaks")
            output.write(path + "\n")
    config = request.get("aracne_config_dir") or str(PACKAGE / "config")
    command = [str(PACKAGE / "bin/sjaracne.exe"), "--consensus-inputs", str(manifest.resolve()),
               "--min-recurrence", str(parameters["min_recurrence"]),
               "-i", request["exp_file"], "-s", request["probe_file"], "-l", request["probe_file"],
               "-H", config, "-N", str(parameters["depth"]), "-o", "consensus.adj"]
    if parameters["dynamic_dpi"]:
        command += ["--dynamic-dpi"]
        for name, flag in (("dpi_href", "--dpi-href"), ("dpi_slope", "--dpi-slope"),
                           ("dpi_max", "--dpi-max"), ("dpi_offset", "--dpi-offset")):
            if name in parameters:
                command += [flag, str(parameters[name])]
    else:
        command += ["-e", str(parameters["tolerance"])]
    return command


def export_native_consensus(adjacency, output_dir, parameters):
    """Stream the final pruned adjacency into the retained three-column format.

    Only the final graph is converted. Native headers remain intact in
    parameter_info_.txt, including recurrence and resolved tolerance metadata.

    Args:
        adjacency (str or Path): Completed native consensus/DPI output.
        output_dir (str or Path): New directory for retained output filenames.
        parameters (dict): Requested run count and minimum recurrence.

    Returns:
        Path: Final three-column network for expression-based enhancement.

    Raises:
        ValueError: Native rows or completion counters are inconsistent.
        FileExistsError: The output directory already exists.
    """
    directory = Path(output_dir)
    directory.mkdir(exist_ok=False)
    network = directory / "consensus_network_3col_.txt"
    headers = []
    edges = 0
    with Path(adjacency).open(encoding="utf-8") as source, network.open("x", encoding="utf-8", newline="\n") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(("source", "target", "MI"))
        seen_sources = set()
        for line in source:
            if line.startswith(">"):
                if seen_sources:
                    raise ValueError("Native consensus header follows a data row")
                headers.append(line)
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) % 2 != 1 or fields[0] in seen_sources:
                raise ValueError("Malformed native consensus row")
            seen_sources.add(fields[0])
            targets = set()
            for index in range(1, len(fields), 2):
                target, weight = fields[index:index + 2]
                value = float(weight)
                if target == fields[0] or target in targets or not math.isfinite(value) or value < 0:
                    raise ValueError("Malformed native consensus edge")
                targets.add(target)
                writer.writerow((fields[0], target, weight))
                edges += 1
    def counter(name):
        prefix = ">  " + name + " "
        values = [line[len(prefix):].strip() for line in headers if line.startswith(prefix)]
        if len(values) != 1:
            raise ValueError("Missing or duplicate native consensus metadata: " + name)
        return int(values[0])
    if (counter("Consensus input runs B") != parameters["bootstrap_num"]
            or counter("Consensus minimum recurrence K") != parameters["min_recurrence"]
            or counter("Consensus final edges") != edges
            or counter("Consensus source rows H") != len(seen_sources)):
        raise ValueError("Native consensus metadata differs from the final graph or requested settings")
    with (directory / "parameter_info_.txt").open("x", encoding="utf-8", newline="\n") as output:
        output.writelines(headers)
        output.write(">  Workflow consensus order consensus-first\n")
        output.write(">  Final stage input already-thresholded mean MI; no MI rethresholding\n")
        output.write(">  Raw inference workflow options JSON " + json.dumps(parameters, sort_keys=True, allow_nan=False) + "\n")
        output.write(">  Source: sjaracne2\n")
    with (directory / "bootstrap_info_.txt").open("x", encoding="utf-8", newline="\n") as output:
        output.write("Consensus selection: minimum recurrence before DPI\n")
        output.write("Bootstrap networks: {}\nMinimum recurrence: {}\n".format(
            parameters["bootstrap_num"], parameters["min_recurrence"]))
        output.write("MI aggregation: supporting-run arithmetic mean\n")
        for line in headers:
            if line.startswith(">  Consensus ") and not line.startswith(">  Consensus input manifest "):
                output.write(line)
    return network



def main():
    """Dispatch an engine-localized JSON request to one WDL task stage."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "infer", "consensus"))
    parser.add_argument("request", help="Engine-localized JSON request file")
    args = parser.parse_args()
    {"prepare": prepare, "infer": infer, "consensus": consensus}[args.stage](read_json(args.request))


if __name__ == "__main__":
    main()
