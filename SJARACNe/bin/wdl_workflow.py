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
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    def reject_constant(value):
        raise ValueError("Non-finite JSON number: " + value)
    with Path(path).open(encoding="utf-8") as handle:
        value = json.load(handle, parse_constant=reject_constant)
    if not isinstance(value, dict):
        raise ValueError("A JSON object is required")
    return value


def write_json(path, value):
    with Path(path).open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def validate_parameters(parameters, observation_count):
    """Use the same profile resolution as the public CLI."""
    return resolve_options(parameters, observation_count)


def normalize_lf(source, destination):
    """Copy bytes, converting CRLF/lone CR only; never mutate source or reorder."""
    with Path(source).open("rb") as source_handle, Path(destination).open("xb") as output:
        pending = b""
        for block in iter(lambda: source_handle.read(1024 * 1024), b""):
            block = pending + block
            pending = b"\r" if block.endswith(b"\r") else b""
            if pending:
                block = block[:-1]
            output.write(block.replace(b"\r\n", b"\n").replace(b"\r", b"\n"))
        output.write(pending.replace(b"\r", b"\n"))


def runtime_fingerprint(config_dir=None):
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
    return "TF_run_{:03d}.adj".format(integer(seed, "seed"))


def native_command(request, parameters):
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
    settings = checked_settings(request)
    parameters = settings["parameters"]
    adjacencies, provenance_files = request["adjacencies"], request["run_provenance"]
    n = parameters["bootstrap_num"]
    settings_sha256 = sha256(request["settings"])
    if len(adjacencies) != n or len(provenance_files) != n:
        raise ValueError("Expected exactly bootstrap_num adjacencies and provenance records")
    # Validate the full expected seed set before creating the consensus directory.
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "infer", "consensus"))
    parser.add_argument("request", help="Engine-localized JSON request file")
    args = parser.parse_args()
    {"prepare": prepare, "infer": infer, "consensus": consensus}[args.stage](read_json(args.request))


if __name__ == "__main__":
    main()
