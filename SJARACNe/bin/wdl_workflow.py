"""Task-boundary plumbing for the standalone WDL workflow, not an inference engine.

Scientific work stays in the existing native program, input QC and consensus
module. Requests arrive as JSON, never as shell-interpolated user arguments.
"""
import argparse
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from .QC_input import validate_inputs


PACKAGE = Path(__file__).resolve().parents[1]
SCHEMA = "sjaracne-wdl-parameters-v1"
OPTION_NAMES = {
    "bootstrap_num", "subsample_spec", "depth", "p_value_bootstrap",
    "apmi_null_model", "min_recurrence", "p_value_consensus", "tolerance",
    "dynamic_dpi", "dpi_href", "dpi_slope", "dpi_max",
}


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


def integer(value, name, minimum=1, maximum=2147483647):
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ValueError("{} must be an integer within [{},{}]".format(name, minimum, maximum))
    return value


def finite_number(value, name, lower, upper, lower_inclusive=True, upper_inclusive=True):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(name + " must be numeric")
    value = float(value)
    if (not math.isfinite(value)
            or (value < lower if lower_inclusive else value <= lower)
            or (value > upper if upper_inclusive else value >= upper)):
        raise ValueError(name + " is outside its supported range")
    return value


def validate_parameters(parameters, observation_count):
    """Resolve workflow defaults and reject conflicts before scattering tasks."""
    if not isinstance(parameters, dict) or set(parameters) - OPTION_NAMES:
        raise ValueError("Unknown workflow parameters")
    values = {key: value for key, value in parameters.items() if value is not None}
    n = integer(values.get("bootstrap_num", 100), "bootstrap_num")
    depth = integer(values.get("depth", 40), "depth")
    spec = values.get("subsample_spec", "80%")
    if not isinstance(spec, str):
        raise ValueError("subsample_spec must be a count or explicit percentage string")
    if re.fullmatch(r"[0-9]+", spec):
        m = int(spec)
    elif re.fullmatch(r"(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)%", spec):
        percentage = float(spec[:-1])
        if not 0 < percentage <= 100:
            raise ValueError("subsample percentage must be within (0,100]")
        # Same operation order as Parameter::resolveSubsampleSize in C++.
        m = math.ceil(percentage * observation_count / 100.0)
    else:
        raise ValueError("subsample_spec must be a count or explicit percentage string")
    if not 2 <= m <= observation_count:
        raise ValueError("selected sample count must be within [2,N]")
    model = values.get("apmi_null_model")
    if model not in (None, "pchip"):
        raise ValueError("apmi_null_model must be omitted or pchip")
    if model == "pchip" and "p_value_bootstrap" not in values:
        raise ValueError("PCHIP requires explicit p_value_bootstrap; no TF/SIG default is chosen")
    pb = finite_number(values.get("p_value_bootstrap", 1e-7), "p_value_bootstrap", 0, 1,
                       lower_inclusive=False)
    if model == "pchip" and pb != 1:
        if depth != 40:
            raise ValueError("PCHIP requires depth 40")
        # The shared callable validates sample-count/probability bounds. It does
        # not set the native cutoff: the executable still resolves its own MI.
        from ..apmi_cutoff import get_mi_cutoff
        get_mi_cutoff(m, pb)
    k, pc = values.get("min_recurrence"), values.get("p_value_consensus")
    if k is not None and pc is not None:
        raise ValueError("min_recurrence and p_value_consensus are mutually exclusive")
    if pc is None:
        k = integer(6 if k is None else k, "min_recurrence")
        if k > n:
            raise ValueError("minimum recurrence exceeds bootstrap_num")
    else:
        pc = finite_number(pc, "p_value_consensus", 0, 1, lower_inclusive=False)
    dynamic = values.get("dynamic_dpi", False)
    if not isinstance(dynamic, bool):
        raise ValueError("dynamic_dpi must be boolean")
    if dynamic and "tolerance" in values:
        raise ValueError("dynamic_dpi conflicts with an explicit fixed tolerance")
    tuning = {key: values[key] for key in ("dpi_href", "dpi_slope", "dpi_max") if key in values}
    if tuning and not dynamic:
        raise ValueError("DPI tuning parameters require dynamic_dpi")
    if "dpi_href" in tuning:
        integer(tuning["dpi_href"], "dpi_href")
    if "dpi_slope" in tuning:
        finite_number(tuning["dpi_slope"], "dpi_slope", 0, math.inf, lower_inclusive=False)
    if "dpi_max" in tuning:
        finite_number(tuning["dpi_max"], "dpi_max", 0, 1, upper_inclusive=False)
    if not dynamic:
        values["tolerance"] = finite_number(values.get("tolerance", 0), "tolerance", 0, 1)
    values.update(bootstrap_num=n, depth=depth, subsample_spec=spec,
                  p_value_bootstrap=pb, dynamic_dpi=dynamic,
                  min_recurrence=k, p_value_consensus=pc, apmi_null_model=model)
    return values


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
               "-p", str(parameters["p_value_bootstrap"]), "-a", "adaptive_partitioning",
               "-u", parameters["subsample_spec"], "-H", config,
               "-N", str(parameters["depth"]), "-o", adjacency_name(seed), "-S", str(seed)]
    if parameters["apmi_null_model"] is not None:
        command += ["-M", parameters["apmi_null_model"]]
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
