#!/usr/bin/env python3
"""Run tiny, real CWL/WDL 1.3 parity checks; never run dataset benchmarks.

Requires the already-built native executable, Sprocket 0.32.0, cwltool, Node.js,
and this checkout's Python dependencies. No tools are installed by this script.
Execution artifacts are retained in a new output directory, including failures.
The local Sprocket backend runs trusted task commands directly on the host.
"""

import argparse
import csv
import datetime
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time
import traceback
from urllib.parse import unquote, urlparse


ROOT = Path(__file__).resolve().parents[2]
WDL = ROOT / "SJARACNe/wdl/sjaracne.wdl"
CWL = ROOT / "SJARACNe/cwl/sjaracne_workflow.cwl"
CONFIG = ROOT / "SJARACNe/wdl/sprocket.local.toml"
NATIVE = ROOT / "SJARACNe/bin/sjaracne.exe"
SPROCKET_VERSION = "0.32.0"
PATH_HEADER = re.compile(r"^(>\s+(?:Input file|Subnetwork file)\s+).*$")
PARAMETER_PATH_HEADER = re.compile(
    r"^(>\s+(?:Input file\s+|Subnetwork file\s+|Output network:\s+)).*$"
)
DPI_STATS = re.compile(
    r"\[DPI_STATS\] pre_edges=(\d+) pruned_edges=(\d+) "
    r"post_edges=(\d+) dpi_applied=(\d+)"
)


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    with Path(path).open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def source_hashes():
    paths = [ROOT / "setup.py", ROOT / "SJARACNe/Makefile"]
    for folder in ("SJARACNe/src", "SJARACNe/config", "SJARACNe/cwl", "SJARACNe/wdl"):
        paths.extend(path for path in (ROOT / folder).rglob("*") if path.is_file())
    paths.extend((ROOT / "SJARACNe").glob("*.py"))
    paths.extend((ROOT / "SJARACNe/bin").glob("*.py"))
    paths.extend([Path(__file__), ROOT / "tests/test_sjaracne_acceptance.py"])
    return {str(path.relative_to(ROOT)): sha256(path) for path in sorted(set(paths))}


def file_path(value):
    """Accept standard output strings and CWL/extended-file output objects."""
    if isinstance(value, dict):
        value = value.get("path") or value.get("location") or value.get("uri")
    if not isinstance(value, str):
        raise AssertionError("Expected a file output, got {!r}".format(value))
    parsed = urlparse(value)
    if parsed.scheme:
        if parsed.scheme != "file" or parsed.netloc not in ("", "localhost"):
            raise AssertionError("Expected a local file output, got " + value)
        value = unquote(parsed.path)
    path = Path(value)
    if not path.is_absolute() or not path.is_file():
        raise AssertionError("Output is not an existing absolute file: " + str(path))
    return path


def header_value(text, key):
    match = re.search(r"^>\s+" + re.escape(key) + r"\s+([^\r\n]+)$", text, re.M)
    if match is None:
        raise AssertionError("Missing scientific header {!r}".format(key))
    return match.group(1).strip()


def native_sections(path):
    """Retain exact body bytes and all headers except two path values."""
    data = Path(path).read_bytes()
    headers, body = [], []
    for line in data.splitlines(keepends=True):
        if line.startswith(b">"):
            decoded = line.decode("utf-8")
            headers.append(PATH_HEADER.sub(r"\1<localized-path>", decoded.rstrip("\n")) + "\n")
        else:
            body.append(line)
    if not headers:
        raise AssertionError("Missing native headers in " + str(path))
    return "".join(headers).encode("utf-8"), b"".join(body)


def normalized_parameters(path):
    return "".join(
        PARAMETER_PATH_HEADER.sub(r"\1<localized-path>", line.rstrip("\n")) + "\n"
        for line in Path(path).read_text(encoding="utf-8").splitlines(keepends=True)
    ).encode("utf-8")


def collect_cwl_adjacencies(folder, count):
    collected = {}
    for path in sorted(folder.rglob("TF_run_*.adj")):
        data = path.read_bytes()
        if path.name in collected and data != collected[path.name].read_bytes():
            raise AssertionError("CWL retained unequal copies of " + path.name)
        collected[path.name] = path
    expected = {"TF_run_{:03d}.adj".format(seed) for seed in range(1, count + 1)}
    if set(collected) != expected:
        raise AssertionError("CWL did not retain exactly the expected seed set: " + repr(sorted(collected)))
    return [collected[name] for name in sorted(collected)]


def unique_identical_file(folder, name):
    paths = sorted(folder.rglob(name))
    if not paths:
        raise AssertionError("Missing retained file " + name)
    reference = paths[0].read_bytes()
    if any(path.read_bytes() != reference for path in paths[1:]):
        raise AssertionError("Retained copies differ: " + name)
    return paths[0]


def sprocket_outputs(folder):
    # Every invocation receives a new -o directory. Find only a direct workflow
    # run's outputs, never nested task outputs. Do not depend on _latest being
    # a usable symlink on a Windows-mounted filesystem.
    candidates = {path.resolve() for path in (folder / "runs/SJARACNe").glob("*/outputs.json")
                  if path.parent.name != "_latest" and path.is_file()}
    if len(candidates) != 1:
        raise AssertionError("Expected exactly one Sprocket workflow outputs.json: " + repr(sorted(map(str, candidates))))
    path = next(iter(candidates))
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise AssertionError("Invalid Sprocket output object")
    if "outputs" in value and not any(key.startswith("SJARACNe.") for key in value):
        value = value["outputs"]
    return path, value


def adjacency_edges(path):
    edges = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line or line.startswith(">"):
            continue
        fields = line.split("\t")
        if len(fields) % 2 != 1 or not fields[0]:
            raise AssertionError("Malformed adjacency row in " + str(path))
        for index in range(1, len(fields), 2):
            edge = (fields[0], fields[index])
            mi = float(fields[index + 1])
            if edge in edges or edge[0] == edge[1] or not math.isfinite(mi):
                raise AssertionError("Invalid or duplicate adjacency edge in " + str(path))
            edges[edge] = mi
    return edges


def consensus_edges(path):
    edges = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"source", "target", "MI"}.issubset(reader.fieldnames):
            raise AssertionError("Missing consensus column headers in " + str(path))
        for row in reader:
            edge = (row["source"], row["target"])
            if edge in edges or not math.isfinite(float(row["MI"])):
                raise AssertionError("Invalid or duplicate consensus edge in " + str(path))
            edges[edge] = row["MI"]
    if list(edges) != sorted(edges):
        raise AssertionError("Consensus edge ordering is not deterministic in " + str(path))
    return edges


def failed_wdl_diagnostic(result):
    # Sprocket may point to task stderr rather than repeat it on its own
    # stderr. Check the retained task diagnostics too, not only the wrapper.
    paths = [Path(result["command"]["stdout"]), Path(result["command"]["stderr"])]
    paths.extend(path for path in Path(result["execution"]).rglob("*")
                 if path.is_file() and "stderr" in path.name.lower())
    return "\n".join(path.read_text(encoding="utf-8", errors="replace") for path in sorted(set(paths)))


class ParityRun:
    def __init__(self, sprocket, output):
        self.sprocket = Path(sprocket).expanduser().resolve()
        self.output = Path(output).expanduser().absolute()
        self.output.mkdir(parents=True, exist_ok=False)
        self.env = os.environ.copy()
        self.env.update({
            "PYTHONPATH": str(ROOT), "PYTHONDONTWRITEBYTECODE": "1", "PYTHONNOUSERSITE": "1",
            "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        })
        self.env["PATH"] = os.pathsep.join((
            str(NATIVE.parent), str(Path(sys.executable).absolute().parent), self.env.get("PATH", ""),
        ))
        self.summary = {
            "schema": "sjaracne-wdl13-parity-v2", "started_utc": now(), "complete": False,
            "passed": False, "failed_checks": 0, "checks": [], "commands": [], "cases": [],
            "source_root": str(ROOT), "python": sys.version,
            "scope": "tiny functional scientific parity only; no dataset benchmark",
            "path_normalization": {
                "native_headers": ["Input file", "Subnetwork file"],
                "parameter_info": ["Input file", "Subnetwork file", "Output network"],
                "other_bytes": "Compared exactly, including adjacency bodies and scientific headers",
            },
        }
        self.persist()

    def persist(self):
        write_json(self.output / "summary.json", self.summary)

    def check(self, name, condition, detail=None):
        entry = {"name": name, "passed": bool(condition)}
        if detail is not None:
            entry["detail"] = detail
        self.summary["checks"].append(entry)
        if not condition:
            self.summary["failed_checks"] += 1
        self.persist()
        if not condition:
            raise AssertionError(name + (": " + str(detail) if detail is not None else ""))

    def command(self, name, argv, expect_failure=False):
        folder = self.output / "commands" / name
        folder.mkdir(parents=True, exist_ok=False)
        entry = {"name": name, "argv": [str(item) for item in argv], "cwd": str(ROOT),
                 "started_utc": now(), "expected_failure": expect_failure,
                 "stdout": str(folder / "stdout.txt"), "stderr": str(folder / "stderr.txt")}
        self.summary["commands"].append(entry)
        self.persist()
        print("START " + name, flush=True)
        started = time.monotonic()
        with (folder / "stdout.txt").open("wb") as stdout, (folder / "stderr.txt").open("wb") as stderr:
            process = subprocess.run(entry["argv"], cwd=ROOT, env=self.env, stdout=stdout, stderr=stderr)
        entry.update(returncode=process.returncode, elapsed_seconds=time.monotonic() - started,
                     stdout_sha256=sha256(folder / "stdout.txt"), stderr_sha256=sha256(folder / "stderr.txt"),
                     completed_utc=now())
        self.persist()
        self.check(name + ":exit", (process.returncode != 0) if expect_failure else (process.returncode == 0),
                   "returncode={}; inspect {}".format(process.returncode, folder))
        print("PASS " + name, flush=True)
        return entry

    def sprocket_command(self, action):
        return [self.sprocket, "--skip-config-search", "--config", CONFIG, action]

    def preflight(self):
        if platform.system() != "Linux":
            raise RuntimeError("Run this parity harness in Linux/WSL with Bash and the Linux native binary")
        for path in (self.sprocket, NATIVE):
            if not path.is_file() or not os.access(str(path), os.X_OK):
                raise RuntimeError("Required executable is absent or not executable: " + str(path))
        tools = {}
        for name in ("python3", "cwltool", "bash"):
            tools[name] = shutil.which(name, path=self.env["PATH"])
            if not tools[name]:
                raise RuntimeError("Missing required command: " + name)
        tools["node"] = next((shutil.which(name, path=self.env["PATH"])
                              for name in ("node", "nodejs") if shutil.which(name, path=self.env["PATH"])), None)
        if not tools["node"]:
            raise RuntimeError("Node.js is required for the real CWL workflow")
        for path in (WDL, CWL, CONFIG):
            if not path.is_file():
                raise RuntimeError("Missing workflow/configuration file: " + str(path))
        self.summary["versions"] = {
            name: importlib.metadata.version(name) for name in ("numpy", "pandas", "scipy", "cwltool")
        }
        self.summary.update(tools=tools, source_sha256=source_hashes(), native_sha256=sha256(NATIVE),
                            sprocket_sha256=sha256(self.sprocket))
        entry = self.command("sprocket_version", [self.sprocket, "--version"])
        version = Path(entry["stdout"]).read_text(encoding="utf-8").strip()
        self.check("pinned_sprocket_version",
                   bool(re.search(r"^sprocket v?" + re.escape(SPROCKET_VERSION) + r"(?:\s|$)", version)), version)
        self.summary["versions"]["sprocket"] = version
        entry = self.command("node_version", [tools["node"], "--version"])
        self.summary["versions"]["node"] = Path(entry["stdout"]).read_text(encoding="utf-8").strip()
        entry = self.command("cwltool_version", [tools["cwltool"], "--version"])
        self.summary["versions"]["cwltool_command"] = (
            Path(entry["stdout"]).read_text(encoding="utf-8")
            + Path(entry["stderr"]).read_text(encoding="utf-8")
        ).strip()
        self.command("sprocket_check", self.sprocket_command("check") + [WDL])
        self.command("cwl_validate", [tools["cwltool"], "--validate", CWL])
        self.persist()

    def fixtures(self):
        sys.path.insert(0, str(ROOT / "tests"))
        from test_sjaracne_acceptance import write_multihub_dpi_inputs
        fixture = self.output / "fixtures" / "paths with spaces"
        fixture.mkdir(parents=True)
        expression, hubs = write_multihub_dpi_inputs(fixture)
        crlf = self.output / "fixtures" / "crlf"
        crlf.mkdir()
        crlf_expression, crlf_hubs = crlf / expression.name, crlf / hubs.name
        crlf_expression.write_bytes(expression.read_bytes().replace(b"\n", b"\r\n"))
        crlf_hubs.write_bytes(hubs.read_bytes().replace(b"\n", b"\r\n"))
        bad_hubs = fixture / "zero_overlap.txt"
        bad_hubs.write_text("NOT_AN_EXPRESSION_GENE\n", encoding="ascii")
        self.fixture_hashes = {
            str(path): sha256(path) for path in (expression, hubs, crlf_expression, crlf_hubs, bad_hubs)
        }
        self.summary["fixture_sha256"] = self.fixture_hashes
        self.persist()
        return expression, hubs, crlf_expression, crlf_hubs, bad_hubs

    def run_wdl(self, name, expression, hubs, options, expect_failure=False):
        folder = self.output / "cases" / name / "wdl"
        folder.mkdir(parents=True)
        job = {"SJARACNe.exp_file": str(expression), "SJARACNe.probe_file": str(hubs),
               "SJARACNe.aracne_config_dir": str(ROOT / "SJARACNe/config"),
               "SJARACNe.prepare_memory": "1 GiB", "SJARACNe.inference_memory": "1 GiB",
               "SJARACNe.consensus_memory": "1 GiB"}
        job.update({"SJARACNe." + key: value for key, value in options.items()})
        inputs = folder / "inputs.json"
        write_json(inputs, job)
        if not expect_failure:
            self.command(name + "_wdl_validate", self.sprocket_command("validate") + [WDL, "@" + str(inputs)])
        output = folder / "execution"
        entry = self.command(name + "_wdl", self.sprocket_command("run") + [WDL, "@" + str(inputs), "-o", output],
                             expect_failure=expect_failure)
        if expect_failure:
            artifacts = [str(path.relative_to(output)) for path in output.rglob("*")
                         if path.name.startswith("TF_run_") or path.name in ("run_provenance.json", "native.stdout.txt")]
            infer_calls = [str(path.relative_to(output)) for path in output.rglob("*")
                           if path.is_dir() and re.fullmatch(r"(?:call[-_])?Infer(?:[-_].*)?", path.name)]
            self.check(name + ":no_native_artifacts", not artifacts, artifacts)
            self.check(name + ":no_infer_calls", not infer_calls, infer_calls)
            return {"command": entry, "execution": str(output), "failed_as_expected": True}
        output_json, values = sprocket_outputs(output)
        required = ("out_dir", "consensus_network", "bootstrap_info", "parameter_info", "parameters",
                    "validation_report", "adjacencies", "run_provenance", "native_stdout", "native_stderr")
        self.check(name + ":wdl_output_contract", all("SJARACNe." + key in values for key in required), sorted(values))
        result = {key: values["SJARACNe." + key] for key in required}
        result.update(command=entry, outputs_json=str(output_json), engine="wdl")
        for key in ("out_dir", "consensus_network", "bootstrap_info", "parameter_info", "parameters", "validation_report"):
            result[key] = file_path(result[key])
        for key in ("adjacencies", "run_provenance", "native_stdout", "native_stderr"):
            result[key] = [file_path(value) for value in result[key]]
        count = options["bootstrap_num"]
        expected_names = ["TF_run_{:03d}.adj".format(seed) for seed in range(1, count + 1)]
        self.check(name + ":wdl_seed_order", [path.name for path in result["adjacencies"]] == expected_names)
        for key in ("run_provenance", "native_stdout", "native_stderr"):
            self.check(name + ":wdl_" + key + "_count", len(result[key]) == count)
        for seed, (adjacency, provenance, stdout) in enumerate(zip(
                result["adjacencies"], result["run_provenance"], result["native_stdout"]), 1):
            record = json.loads(provenance.read_text(encoding="utf-8"))
            self.check(name + ":seed{}:provenance".format(seed),
                       record.get("seed") == seed and record.get("returncode") == 0
                       and record.get("adjacency_sha256") == sha256(adjacency)
                       and record.get("runtime", {}).get("native_sha256") == self.summary["native_sha256"])
            statistics = DPI_STATS.findall(stdout.read_text(encoding="utf-8"))
            self.check(name + ":seed{}:real_dpi".format(seed), len(statistics) == 1
                       and int(statistics[0][3]) == 1
                       and int(statistics[0][2]) > 0
                       and (options.get("profile") == "legacy" or int(statistics[0][1]) > 0)
                       and int(statistics[0][0]) == int(statistics[0][1]) + int(statistics[0][2]))
        return result

    def run_cwl(self, name, expression, hubs, options):
        folder = self.output / "cases" / name / "cwl"
        folder.mkdir(parents=True)
        (folder / "work").mkdir()
        (folder / "retained").mkdir()
        job = {"exp_file": {"class": "File", "path": str(expression)},
               "probe_file": {"class": "File", "path": str(hubs)},
               "aracne_config_dir": {"class": "Directory", "path": str(ROOT / "SJARACNe/config")},
               "final_out_dir_name": "results"}
        # Resolve the selected scientific profile explicitly for the retained
        # compatibility DAG; choosing CWL must not choose a different profile.
        from SJARACNe.workflow_options import resolve_options
        with expression.open(encoding="utf-8-sig") as handle:
            observations = len(handle.readline().rstrip("\r\n").split("\t")) - 2
        parameters = resolve_options(options, observations)
        for key in ("bootstrap_num", "depth", "p_value_bootstrap", "min_recurrence",
                    "p_value_consensus", "dynamic_dpi", "dpi_href", "dpi_slope", "dpi_max", "tolerance"):
            if parameters.get(key) is not None:
                job[key] = parameters[key]
        job["apmi_null_model"] = (None if parameters["apmi_null_model"] == "affine"
                                  else parameters["apmi_null_model"])
        if parameters["sampling_method"] == "with-replacement":
            job.update(sample_number=1, subsample_spec=None)
        else:
            job["subsample_spec"] = parameters["subsample_spec"]
        inputs = folder / "inputs.json"
        write_json(inputs, job)
        output = folder / "published"
        entry = self.command(name + "_cwl", [
            self.summary["tools"]["cwltool"], "--leave-tmpdir", "--preserve-environment", "PYTHONPATH",
            "--outdir", output, "--tmpdir-prefix", str(folder / "work") + "/",
            "--tmp-outdir-prefix", str(folder / "retained") + "/", CWL, inputs,
        ])
        values = json.loads(Path(entry["stdout"]).read_text(encoding="utf-8"))
        required = ("out_dir", "bootstrap_info", "parameter_info")
        self.check(name + ":cwl_output_contract", all(key in values for key in required), sorted(values))
        result = {key: file_path(values[key]) for key in required}
        result.update(command=entry, consensus_network=unique_identical_file(folder, "consensus_network_3col_.txt"),
                      adjacencies=collect_cwl_adjacencies(folder, options["bootstrap_num"]), engine="cwl")
        return result

    def independent_consensus(self, name, result, options):
        # Neither engine is the oracle: independently count distinct per-seed
        # edges and average their printed native MI values in seed order.
        networks = [adjacency_edges(path) for path in result["adjacencies"]]
        support, total_mi = {}, {}
        for network in networks:
            for edge, mi in network.items():
                support[edge] = support.get(edge, 0) + 1
                total_mi[edge] = total_mi.get(edge, 0.0) + mi
        metadata = dict(line.split(": ", 1) for line in result["bootstrap_info"].read_text(encoding="utf-8").splitlines())
        prefix = name + ":" + result["engine"] + ":independent:"
        self.check(prefix + "union_size", metadata.get("Total edge tested") == str(len(support)))
        count = options["bootstrap_num"]
        pc = options.get("p_value_consensus", 1e-5 if options.get("profile") == "legacy" else None)
        if pc is None:
            k = options.get("min_recurrence", 6)
            expected_metadata = {
                "Total edge tested": str(len(support)), "Consensus selection": "minimum recurrence",
                "Bootstrap networks": str(count), "Minimum recurrence": str(k),
                "Minimum recurrence fraction": "{:.12g}".format(k / count),
            }
            self.check(prefix + "recurrence_metadata", metadata == expected_metadata)
            selected = {edge for edge, frequency in support.items() if frequency >= k}
        else:
            edge_count = len(support)
            probabilities = [len(network) / edge_count for network in networks] if edge_count else []
            mu = sum(probabilities)
            sigma = math.sqrt(sum(value * (1 - value) for value in probabilities))
            self.check(prefix + "legacy_mu", math.isclose(float(metadata["mu"]), mu, rel_tol=0, abs_tol=1e-14))
            self.check(prefix + "legacy_sigma", math.isclose(float(metadata["sigma"]), sigma, rel_tol=0, abs_tol=1e-14))
            # The historical implementation uses an approximation to the
            # Gaussian upper tail. An independent erfc check is conclusive for
            # this fixture only when no edge is close to the decision boundary.
            tails = {edge: 0.5 * math.erfc(((frequency - mu) / sigma if sigma else 100) / math.sqrt(2))
                     for edge, frequency in support.items()}
            self.check(prefix + "legacy_gate_not_near_approximation_boundary",
                       all(abs(value - pc) > max(1e-12, 0.01 * pc) for value in tails.values()))
            selected = {edge for edge, value in tails.items() if value < pc}
        expected = {edge: "{:.4f}".format(total_mi[edge] / support[edge]) for edge in sorted(selected)}
        for key in ("consensus_network", "out_dir"):
            self.check(prefix + key + ":recurrence_and_mean_mi", consensus_edges(result[key]) == expected)
        self.summary.setdefault("independent_consensus", {})[name + ":" + result["engine"]] = {
            "resamples": count, "union_edges": len(support), "retained_edges": len(expected),
            "edge_support": {source + "\t" + target: support[(source, target)] for source, target in sorted(support)},
        }
        self.persist()

    def compare(self, name, left, right, options):
        count = options["bootstrap_num"]
        self.check(name + ":same_seed_count", len(left["adjacencies"]) == len(right["adjacencies"]) == count)
        files = {}
        for seed, (a, b) in enumerate(zip(left["adjacencies"], right["adjacencies"]), 1):
            ah, ab = native_sections(a)
            bh, bb = native_sections(b)
            self.check(name + ":seed{}:raw_body_bytes".format(seed), ab == bb)
            self.check(name + ":seed{}:scientific_headers".format(seed), ah == bh)
            files["seed{:03d}".format(seed)] = {
                "left": str(a), "right": str(b), "left_sha256": sha256(a), "right_sha256": sha256(b),
                "body_sha256": hashlib.sha256(ab).hexdigest(),
                "scientific_headers_sha256": hashlib.sha256(ah).hexdigest(),
            }
        for key in ("consensus_network", "out_dir", "bootstrap_info"):
            self.check(name + ":" + key + ":bytes", left[key].read_bytes() == right[key].read_bytes())
            files[key] = {"left": str(left[key]), "right": str(right[key]), "sha256": sha256(left[key])}
        self.check(name + ":parameter_info:scientific_bytes",
                   normalized_parameters(left["parameter_info"]) == normalized_parameters(right["parameter_info"]))
        self.independent_consensus(name, left, options)
        self.independent_consensus(name, right, options)
        self.summary["cases"].append({"name": name, "passed": True, "resamples": count, "files": files})
        self.persist()

    def run(self):
        try:
            self.preflight()
            expression, hubs, crlf_expression, crlf_hubs, bad_hubs = self.fixtures()
            cases = (
                ("pchip_strict_default_k6", {"bootstrap_num": 6, "network_type": "TF", "apmi_null_model": "pchip",
                                             "p_value_bootstrap": 0.001, "dynamic_dpi": False}),
                ("pchip_positive_dynamic", {"bootstrap_num": 6, "network_type": "TF", "apmi_null_model": "pchip", "p_value_bootstrap": 0.001,
                                            "dynamic_dpi": True, "dpi_href": 2, "dpi_slope": 0.1, "dpi_max": 0.3}),
                ("affine_explicit_compatibility_k6", {"bootstrap_num": 6, "network_type": "TF", "apmi_null_model": "affine",
                                                       "p_value_bootstrap": 1e-7, "dynamic_dpi": False}),
                ("explicit_legacy_consensus", {"bootstrap_num": 3, "network_type": "TF", "apmi_null_model": "pchip", "p_value_bootstrap": 0.001,
                                               "dynamic_dpi": False, "p_value_consensus": 0.01}),
                ("current_tf_defaults", {"bootstrap_num": 6, "network_type": "TF"}),
                ("current_sig_defaults", {"bootstrap_num": 6, "network_type": "SIG"}),
                ("presampling_legacy_profile", {"bootstrap_num": 6, "profile": "legacy"}),
            )
            baseline = None
            for name, options in cases:
                cwl = self.run_cwl(name, expression, hubs, options)
                wdl = self.run_wdl(name, expression, hubs, options)
                self.compare(name, cwl, wdl, options)
                if name == "pchip_strict_default_k6":
                    baseline = wdl
                    self.check(name + ":nonempty_default_consensus", len(wdl["consensus_network"].read_bytes().splitlines()) > 1)
                    for path in wdl["adjacencies"]:
                        text = path.read_text(encoding="utf-8")
                        self.check(name + ":" + path.name + ":pchip",
                                   header_value(text, "MI threshold method") == "experimental PCHIP empirical cutoff"
                                   and header_value(text, "Sampling request") == "80%"
                                   and header_value(text, "AP-MI actual m") == "80")
                elif name == "pchip_positive_dynamic":
                    expected = min(0.3, 0.1 * max(0, math.log2(3 / 2)))
                    for path in wdl["adjacencies"]:
                        value = float(header_value(path.read_text(encoding="utf-8"), "DPI epsilon binary64"))
                        self.check(name + ":" + path.name + ":positive_epsilon", value > 0
                                   and math.isclose(value, expected, rel_tol=0, abs_tol=1e-15))
                elif name == "affine_explicit_compatibility_k6":
                    for path in wdl["adjacencies"]:
                        text = path.read_text(encoding="utf-8")
                        self.check(name + ":" + path.name + ":compatibility_pb",
                                   float(header_value(text, "MI P-value")) == 1e-7 and "PCHIP" not in text)
                elif name in ("current_tf_defaults", "current_sig_defaults"):
                    expected_pb = 0.001 if name == "current_tf_defaults" else 0.0005
                    parameters = json.loads(wdl["parameters"].read_text(encoding="utf-8"))["parameters"]
                    self.check(name + ":resolved_defaults", parameters["apmi_null_model"] == "pchip"
                               and parameters["p_value_bootstrap"] == expected_pb
                               and parameters["dynamic_dpi"] is True
                               and parameters["sampling_method"] == "without-replacement"
                               and parameters["subsample_spec"] == "80%"
                               and parameters["min_recurrence"] == 6)
                    for adjacency, provenance in zip(wdl["adjacencies"], wdl["run_provenance"]):
                        text = adjacency.read_text(encoding="utf-8")
                        argv = json.loads(provenance.read_text(encoding="utf-8"))["argv"]
                        self.check(name + ":" + adjacency.name + ":native_defaults",
                                   float(header_value(text, "MI P-value")) == expected_pb
                                   and header_value(text, "MI threshold method") == "experimental PCHIP empirical cutoff"
                                   and header_value(text, "DPI mode") == "experimental-dynamic"
                                   and header_value(text, "AP-MI actual m") == "80"
                                   and header_value(text, "DPI H_ref") == "479"
                                   and "--dynamic-dpi" in argv and "-e" not in argv and "-r" not in argv)
                elif name == "presampling_legacy_profile":
                    parameters = json.loads(wdl["parameters"].read_text(encoding="utf-8"))["parameters"]
                    self.check(name + ":resolved_defaults", parameters["apmi_null_model"] == "affine"
                               and parameters["p_value_bootstrap"] == 1e-7
                               and parameters["dynamic_dpi"] is False
                               and parameters["sampling_method"] == "with-replacement"
                               and parameters["subsample_spec"] is None
                               and parameters["p_value_consensus"] == 1e-5
                               and parameters["min_recurrence"] is None)
                    expected_cutoff = (1.062 - math.log(1e-7)) / (48.7 + 0.634 * 100)
                    for adjacency, provenance in zip(wdl["adjacencies"], wdl["run_provenance"]):
                        text = adjacency.read_text(encoding="utf-8")
                        argv = json.loads(provenance.read_text(encoding="utf-8"))["argv"]
                        self.check(name + ":" + adjacency.name + ":native_legacy_sampling",
                                   header_value(text, "Sampling method") == "legacy bootstrap with replacement"
                                   and header_value(text, "Resample draws") == "100"
                                   and header_value(text, "DPI tolerance") == "0"
                                   and math.isclose(float(header_value(text, "MI threshold")), expected_cutoff,
                                                    rel_tol=0, abs_tol=5e-7)
                                   and "-r" in argv and argv[argv.index("-r") + 1] == "1"
                                   and "-u" not in argv and "-M" not in argv and "--dynamic-dpi" not in argv
                                   and argv[argv.index("-N") + 1] == "40")

            crlf_options = dict(cases[0][1])
            crlf = self.run_wdl("wdl_crlf_normalization", crlf_expression, crlf_hubs, crlf_options)
            self.compare("wdl_crlf_vs_lf", baseline, crlf, crlf_options)
            lf_settings = json.loads(baseline["parameters"].read_text(encoding="utf-8"))
            crlf_settings = json.loads(crlf["parameters"].read_text(encoding="utf-8"))
            self.check("crlf_normalized_input_hashes", lf_settings["input_sha256"] == crlf_settings["input_sha256"])
            self.check("crlf_original_bytes_differ", lf_settings["original_input_sha256"] != crlf_settings["original_input_sha256"])
            self.summary["crlf_scope"] = (
                "WDL CRLF results compared to LF baseline; CWL's historical in-place line-ending helper is not a correctness oracle"
            )

            missing_type = self.run_wdl("invalid_current_missing_network_type", expression, hubs,
                                        {"bootstrap_num": 6}, expect_failure=True)
            message = failed_wdl_diagnostic(missing_type)
            self.check("invalid_current_missing_network_type:diagnostic", "network_type" in message)
            zero_overlap = self.run_wdl("invalid_zero_overlap", expression, bad_hubs, dict(cases[0][1]), expect_failure=True)
            message = failed_wdl_diagnostic(zero_overlap)
            self.check("invalid_zero_overlap:diagnostic", "zero hub genes" in message)
            self.check("fixture_inputs_unchanged", all(sha256(path) == digest for path, digest in self.fixture_hashes.items()))
            self.check("tested_sources_unchanged", self.summary["source_sha256"] == source_hashes())
            self.check("tested_native_unchanged", self.summary["native_sha256"] == sha256(NATIVE))
            self.check("tested_sprocket_unchanged", self.summary["sprocket_sha256"] == sha256(self.sprocket))
            self.summary.update(passed=True, complete=True, completed_utc=now())
            self.persist()
            print("PASS: {} checks, 0 failed; {}".format(len(self.summary["checks"]), self.output / "summary.json"), flush=True)
            return 0
        except Exception as error:
            self.summary.update(passed=False, complete=True, completed_utc=now(),
                                error={"type": type(error).__name__, "message": str(error), "traceback": traceback.format_exc()})
            if not self.summary["failed_checks"]:
                self.summary["failed_checks"] = 1
            self.persist()
            print("FAIL: {}; inspect {}".format(error, self.output / "summary.json"), file=sys.stderr, flush=True)
            return 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sprocket", required=True, help="Path to the pinned Sprocket 0.32.0 executable")
    parser.add_argument("--output-dir", required=True, help="New directory for retained execution artifacts and summary.json")
    args = parser.parse_args()
    try:
        run = ParityRun(args.sprocket, args.output_dir)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return run.run()


if __name__ == "__main__":
    sys.exit(main())
