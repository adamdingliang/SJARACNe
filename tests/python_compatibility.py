#!/usr/bin/env python3
"""Compare tiny scientific outputs across two installed Python environments.

This is a compatibility check, not a dataset benchmark. Both environments run
the same checkout, native executable, fixture, seeds and Sprocket executable.
Nothing is installed and no existing output directory is reused.
"""
import argparse
import csv
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

from test_sjaracne_acceptance import write_multihub_dpi_inputs
from wdl.run_parity import (adjacency_edges, consensus_edges, file_path,
                            native_sections, normalized_parameters, source_hashes)


ROOT = Path(__file__).resolve().parents[1]
NATIVE = ROOT / "SJARACNe/bin/sjaracne.exe"
NUMERIC_COLUMNS = ("pearson", "spearman", "slope", "p-value")
ENHANCED_COLUMNS = ("source", "target", "source.symbol", "target.symbol", "MI") + NUMERIC_COLUMNS
RTOL = 1e-12
ATOL = 1e-14
CASES = (
    ("current_tf", ("--network-type", "TF")),
    ("current_sig", ("--network-type", "SIG")),
    # Explicit legacy tuning keeps this small compatibility fixture nonempty;
    # this is not the legacy profile's default pb/pc operating point.
    ("legacy_probability", ("--profile", "legacy", "-pb", "1", "-pc", "0.1")),
)


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
                          encoding="utf-8")


def close_number(left, right, rtol=RTOL, atol=ATOL):
    """Undefined statistics must remain undefined; finite values stay close."""
    left, right = float(left), float(right)
    if math.isnan(left) or math.isnan(right):
        return math.isnan(left) and math.isnan(right)
    if not math.isfinite(left) or not math.isfinite(right):
        return left == right
    return math.isclose(left, right, rel_tol=rtol, abs_tol=atol)


def enhanced_rows(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ENHANCED_COLUMNS:
            raise AssertionError("Unexpected enhanced-network columns: " + str(path))
        rows = list(reader)
    keys = [(row["source"], row["target"]) for row in rows]
    if keys != sorted(set(keys)):
        raise AssertionError("Enhanced edges must be unique and sorted: " + str(path))
    return rows


def compare_enhanced(left, right):
    """Exact identities/MI; small float64 tolerance only for enhancement stats."""
    left_rows, right_rows = enhanced_rows(left), enhanced_rows(right)
    if not left_rows or len(left_rows) != len(right_rows):
        raise AssertionError("Enhanced outputs must contain the same nonzero edge count")
    differences = {name: {"max_abs_difference": 0.0, "matched_nonfinite": 0}
                   for name in NUMERIC_COLUMNS}
    for left_row, right_row in zip(left_rows, right_rows):
        for name in ENHANCED_COLUMNS[:5]:
            if left_row[name] != right_row[name]:
                raise AssertionError("Enhanced identity/MI differs in column " + name)
        for name in NUMERIC_COLUMNS:
            a, b = float(left_row[name]), float(right_row[name])
            if not close_number(a, b):
                raise AssertionError("Enhanced {} differs for {} -> {}: {} versus {}".format(
                    name, left_row["source"], left_row["target"], a, b))
            if math.isfinite(a) and math.isfinite(b):
                differences[name]["max_abs_difference"] = max(
                    differences[name]["max_abs_difference"], abs(a - b))
            else:
                differences[name]["matched_nonfinite"] += 1
    return {"edges": len(left_rows), "statistics": differences,
            "rtol": RTOL, "atol": ATOL}


def compare_metadata(left, right):
    def fields(path):
        return dict(line.split(": ", 1) for line in Path(path).read_text(encoding="utf-8").splitlines())
    a, b = fields(left), fields(right)
    if set(a) != set(b):
        raise AssertionError("Consensus metadata fields differ")
    for name in a:
        if name in ("mu", "sigma", "Bonferroni corrected (0.05) alpha"):
            if not close_number(a[name], b[name]):
                raise AssertionError("Consensus numeric metadata differs: " + name)
        elif a[name] != b[name]:
            raise AssertionError("Consensus metadata differs: " + name)


def independent_selection(adjacencies, options, consensus):
    """Check selected identities and mean MI without using either run as oracle."""
    networks = [adjacency_edges(path) for path in adjacencies]
    counts, total = {}, {}
    for network in networks:
        for edge, mi in network.items():
            counts[edge] = counts.get(edge, 0) + 1
            total[edge] = total.get(edge, 0.0) + mi
    if not counts:
        raise AssertionError("Compatibility fixture emitted no native edges")
    probability = options["p_value_consensus"]
    if probability is None:
        selected = {edge for edge, count in counts.items() if count >= options["min_recurrence"]}
        selection = {"mode": "minimum recurrence", "K": options["min_recurrence"]}
    else:
        rates = [len(network) / len(counts) for network in networks]
        mu = sum(rates)
        sigma = math.sqrt(sum(rate * (1 - rate) for rate in rates))
        tails = {edge: 0.5 * math.erfc(((count - mu) / sigma if sigma else 100) / math.sqrt(2))
                 for edge, count in counts.items()}
        # The production legacy gate uses its historical Gaussian approximation.
        # erfc is an independent oracle only away from that approximation's
        # decision boundary, so require a generous 1% probability margin.
        if any(abs(tail - probability) <= max(1e-12, probability * 0.01) for tail in tails.values()):
            raise AssertionError("Legacy fixture edge is too close to the approximation decision boundary")
        selected = {edge for edge, tail in tails.items() if tail < probability}
        selection = {"mode": "legacy normal probability", "pc": probability,
                     "mu": mu, "sigma": sigma,
                     "minimum_probability_margin": min(abs(tail - probability) for tail in tails.values())}
    expected = {edge: "{:.4f}".format(total[edge] / counts[edge]) for edge in sorted(selected)}
    if not expected or consensus_edges(consensus) != expected:
        raise AssertionError("Consensus does not match independent nonempty selection and mean MI")
    selection.update(union_edges=len(counts), retained_edges=len(expected),
                     source_target_support={"\t".join(edge): count for edge, count in sorted(counts.items())})
    return selection


class CompatibilityRun:
    def __init__(self, baseline_python, candidate_python, sprocket, output):
        # Do not resolve Python symlinks: that would bypass the requested venv.
        self.pythons = {"baseline": Path(baseline_python).expanduser().absolute(),
                        "candidate": Path(candidate_python).expanduser().absolute()}
        self.sprocket = Path(sprocket).expanduser().resolve()
        self.output = Path(output).expanduser().absolute()
        self.output.mkdir(parents=True, exist_ok=False)
        self.summary = {"schema": "sjaracne-python-compatibility-v1", "started_utc": now(),
                        "complete": False, "passed": False, "cases": [], "commands": [],
                        "scope": "tiny cross-Python compatibility; no dataset benchmark",
                        "numeric_policy": {"exact": "native body/scientific headers, selected identities and MI",
                                           "tolerant": list(NUMERIC_COLUMNS) + ["legacy mu/sigma"],
                                           "rtol": RTOL, "atol": ATOL,
                                           "reason": "small float64 library/version roundoff only; no relaxed edge selection"}}
        self.persist()

    def persist(self):
        write_json(self.output / "summary.json", self.summary)

    def environment(self, python):
        env = os.environ.copy()
        env.update(PYTHONPATH=str(ROOT), PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1",
                   OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        env["PATH"] = os.pathsep.join((str(python.parent), str(NATIVE.parent),
                                       str(self.sprocket.parent), env.get("PATH", "")))
        # The active interpreter is selected by absolute path/PATH, never an
        # inherited virtualenv activation from the coordinator.
        env.pop("VIRTUAL_ENV", None)
        env.pop("PYTHONHOME", None)
        return env

    def command(self, name, argv, env):
        folder = self.output / "commands"
        folder.mkdir(exist_ok=True)
        entry = {"name": name, "argv": list(map(str, argv)), "started_utc": now(),
                 "stdout": str(folder / (name + ".stdout.txt")),
                 "stderr": str(folder / (name + ".stderr.txt"))}
        self.summary["commands"].append(entry)
        self.persist()
        started = time.monotonic()
        with Path(entry["stdout"]).open("w") as stdout, Path(entry["stderr"]).open("w") as stderr:
            result = subprocess.run(entry["argv"], cwd=str(ROOT), env=env, stdout=stdout,
                                    stderr=stderr, timeout=600, check=False)
        entry.update(returncode=result.returncode, elapsed_seconds=time.monotonic() - started)
        self.persist()
        if result.returncode:
            raise RuntimeError("{} failed; retained diagnostics: {}".format(name, entry["stderr"]))
        return entry

    def run_case(self, label, name, options, expression, hubs):
        python = self.pythons[label]
        output = self.output / "cases" / name / label
        argv = [python, "-m", "SJARACNe.sjaracne", "local", "--sprocket", self.sprocket,
                "--jobs", "2", "-e", expression, "-g", hubs, "-n", "6", "-o", output]
        self.command(name + "_" + label, argv + list(options), self.environment(python))
        manifest = json.loads((output / "sjaracne_run.json").read_text(encoding="utf-8"))
        if manifest["state"] != "complete" or manifest["engine"] != "sprocket":
            raise AssertionError("Incomplete or incorrect workflow engine")
        outputs = json.loads(Path(manifest["engine_outputs_json"]).read_text(encoding="utf-8"))
        adjacencies = [file_path(value) for value in outputs["SJARACNe.adjacencies"]]
        expected_names = ["TF_run_{:03d}.adj".format(seed) for seed in range(1, 7)]
        if [path.name for path in adjacencies] != expected_names:
            raise AssertionError("Unexpected native seed order")
        records = [file_path(value) for value in outputs["SJARACNe.run_provenance"]]
        if len(records) != len(adjacencies):
            raise AssertionError("Missing seed provenance")
        for seed, (path, record) in enumerate(zip(adjacencies, records), 1):
            provenance = json.loads(record.read_text(encoding="utf-8"))
            if (provenance["seed"] != seed or provenance["returncode"] != 0
                    or provenance["adjacency_sha256"] != sha256(path)
                    or provenance["runtime"] != manifest["runtime"]):
                raise AssertionError("Invalid native seed provenance")
        paths = {key: output / filename for key, filename in (
            ("consensus", "consensus_network_3col_.txt"), ("enhanced", "consensus_network_ncol_.txt"),
            ("metadata", "bootstrap_info_.txt"), ("parameters", "parameter_info_.txt"))}
        selection = independent_selection(adjacencies, manifest["scientific_options"], paths["consensus"])
        return {"manifest": manifest, "adjacencies": adjacencies, "paths": paths,
                "independent_selection": selection}

    def compare(self, name, baseline, candidate):
        if baseline["manifest"]["scientific_options"] != candidate["manifest"]["scientific_options"]:
            raise AssertionError("Resolved scientific options differ")
        if baseline["manifest"]["runtime"] != candidate["manifest"]["runtime"]:
            raise AssertionError("Source/native/configuration fingerprints differ")
        native = []
        for seed, (a, b) in enumerate(zip(baseline["adjacencies"], candidate["adjacencies"]), 1):
            ah, ab = native_sections(a)
            bh, bb = native_sections(b)
            if ah != bh or ab != bb:
                raise AssertionError("Native scientific headers/body differ for seed " + str(seed))
            native.append({"seed": seed, "body_sha256": hashlib.sha256(ab).hexdigest(),
                           "scientific_headers_sha256": hashlib.sha256(ah).hexdigest()})
        a, b = baseline["paths"], candidate["paths"]
        if a["consensus"].read_bytes() != b["consensus"].read_bytes():
            raise AssertionError("Consensus identities/MI differ")
        if normalized_parameters(a["parameters"]) != normalized_parameters(b["parameters"]):
            raise AssertionError("Normalized scientific parameter metadata differs")
        compare_metadata(a["metadata"], b["metadata"])
        enhanced = compare_enhanced(a["enhanced"], b["enhanced"])
        self.summary["cases"].append({"name": name, "passed": True, "resamples_per_environment": 6,
                                      "native": native, "enhanced": enhanced,
                                      "independent_selection": baseline["independent_selection"],
                                      "files": {label: {key: {"path": str(path), "sha256": sha256(path)}
                                                        for key, path in result["paths"].items()}
                                                for label, result in (("baseline", baseline), ("candidate", candidate))}})
        self.persist()

    def run(self):
        try:
            self.summary.update(source_sha256=source_hashes(), comparator_sha256=sha256(__file__),
                                native_sha256=sha256(NATIVE),
                                sprocket_sha256=sha256(self.sprocket))
            version_code = ("import importlib.metadata as m,json,platform,sys; "
                            "print(json.dumps({'python':platform.python_version(),'executable':sys.executable,"
                            "'packages':{n:m.version(n) for n in ('numpy','pandas','scipy')}}))")
            versions = {}
            for label, python in self.pythons.items():
                if not python.is_file() or not os.access(str(python), os.X_OK):
                    raise ValueError("Missing Python executable: " + str(python))
                entry = self.command(label + "_versions", [python, "-c", version_code], self.environment(python))
                versions[label] = json.loads(Path(entry["stdout"]).read_text(encoding="utf-8"))
            self.summary["environments"] = versions
            if versions["baseline"]["python"] == versions["candidate"]["python"]:
                raise ValueError("Cross-Python comparison requires two different interpreter versions")
            fixture = self.output / "fixture"
            fixture.mkdir()
            expression, hubs = write_multihub_dpi_inputs(fixture)
            hashes = {str(path): sha256(path) for path in (expression, hubs)}
            self.summary["fixture_sha256"] = hashes
            for name, options in CASES:
                results = {label: self.run_case(label, name, options, expression, hubs)
                           for label in self.pythons}
                self.compare(name, results["baseline"], results["candidate"])
            if hashes != {path: sha256(path) for path in hashes}:
                raise AssertionError("Fixture inputs changed")
            if (self.summary["source_sha256"] != source_hashes()
                    or self.summary["comparator_sha256"] != sha256(__file__)
                    or self.summary["native_sha256"] != sha256(NATIVE)
                    or self.summary["sprocket_sha256"] != sha256(self.sprocket)):
                raise AssertionError("Tested sources, native program or Sprocket changed during comparison")
            self.summary.update(complete=True, passed=True, finished_utc=now())
            self.persist()
            return 0
        except Exception:
            self.summary.update(complete=True, passed=False, finished_utc=now(), error=traceback.format_exc())
            self.persist()
            print(self.summary["error"], file=sys.stderr)
            return 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-python", required=True)
    parser.add_argument("--candidate-python", required=True)
    parser.add_argument("--sprocket", required=True)
    parser.add_argument("--output-dir", required=True, help="New evidence directory; must not already exist")
    args = parser.parse_args()
    return CompatibilityRun(args.baseline_python, args.candidate_python,
                            args.sprocket, args.output_dir).run()


if __name__ == "__main__":
    sys.exit(main())
