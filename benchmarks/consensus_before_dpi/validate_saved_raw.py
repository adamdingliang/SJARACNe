#!/usr/bin/env python3
"""Compare native fused consensus/DPI with frozen DREAM outputs, without MI jobs.

Run under Linux/WSL with the retained raw graphs and the two completed benchmark
roots. Historical raw files are canonicalized only in private temporary files:
source-row order and omitted empty rows change, numeric MI tokens do not.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import time


def sha256(path):
    """Return a file's SHA-256 without loading it into memory."""
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def local_path(value):
    """Resolve recorded Windows D-drive paths when executing in WSL."""
    value = str(value)
    if len(value) > 2 and value[1] == ":":
        value = "/mnt/" + value[0].lower() + "/" + value[3:].replace("\\", "/")
    return Path(value)


def read_tsv(path):
    """Read a small metadata table, preserving exact numerical strings."""
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_json(path, value):
    """Create a deterministic LF JSON receipt, refusing overwrites."""
    with Path(path).open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def canonicalize(source, destination, hubs):
    """Reorder legacy raw source rows and add explicit missing empty rows.

    Args:
        source (Path): Frozen thresholded, unpruned raw native adjacency.
        destination (Path): New temporary canonical adjacency path.
        hubs (list[str]): Authoritative source IDs in native input order.

    Returns:
        dict: Source/canonical hashes and an exact description of the conversion.

    Raises:
        ValueError: A source is duplicated/unknown or a row is malformed.
    """
    headers, rows = [], {}
    hub_set = set(hubs)
    with source.open(encoding="utf-8-sig") as stream:
        for line in stream:
            if line.startswith(">"):
                if rows:
                    raise ValueError("Header after data in frozen raw graph")
                headers.append(line.rstrip("\r\n"))
                continue
            line = line.rstrip("\r\n")
            fields = line.split("\t")
            if len(fields) % 2 != 1 or fields[0] not in hub_set or fields[0] in rows:
                raise ValueError("Malformed/duplicate/unknown raw source row")
            rows[fields[0]] = line
    with destination.open("x", encoding="utf-8", newline="\n") as stream:
        for line in headers:
            stream.write(line + "\n")
        for hub in hubs:
            stream.write(rows.get(hub, hub) + "\n")
    return {"source": str(source), "source_sha256": sha256(source),
            "canonical_sha256": sha256(destination), "source_rows": len(rows),
            "canonical_rows": len(hubs), "MI_tokens_changed": False,
            "added_empty_rows": len(hubs) - len(rows)}


def load_graph(path):
    """Read a native graph and reject duplicate directed edges."""
    edges, headers = {}, []
    with path.open(encoding="utf-8-sig") as stream:
        first = stream.readline().rstrip("\r\n").split("\t")
        stream.seek(0)
        if {"source", "target", "MI"}.issubset(first):
            for row in csv.DictReader(stream, delimiter="\t"):
                key = row["source"], row["target"]
                if key in edges:
                    raise ValueError("Duplicate directed reference edge")
                edges[key] = float(row["MI"])
            return edges, headers
        for line in stream:
            if line.startswith(">"):
                headers.append(line.rstrip("\r\n"))
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) % 2 != 1:
                raise ValueError("Malformed output adjacency")
            for index in range(1, len(fields), 2):
                key = fields[0], fields[index]
                if key in edges:
                    raise ValueError("Duplicate directed output edge")
                edges[key] = float(fields[index + 1])
    return edges, headers


def raw_inputs(dataset, comparison, classic_raw):
    """Locate the exact 100 frozen raw files and their completion hashes."""
    inputs = []
    for seed in range(1, 101):
        if dataset in ("DREAM5_N1", "DREAM5_N3"):
            receipt = classic_raw / "jobs" / dataset / ("pre_dpi_%03d" % seed) / "complete.json"
            record = json.loads(receipt.read_text(encoding="utf-8"))["graph"]
        else:
            receipt = comparison / "work" / dataset / "current" / ("seed_%03d" % seed) / "complete.json"
            item = json.loads(receipt.read_text(encoding="utf-8"))
            if item["seed"] != seed or item["kind"] != "current":
                raise ValueError("Wrong frozen seed/kind")
            record = next(row for row in item["outputs"] if row["role"] == "current_raw")
        path = local_path(record["path"])
        if sha256(path) != record["sha256"]:
            raise ValueError("Frozen raw graph hash mismatch: " + str(path))
        inputs.append((seed, path, receipt, sha256(receipt)))
    return inputs


def run_dataset(descriptor, args, rows, binary_sha):
    """Validate both offsets for one dataset using the same canonical raw files."""
    dataset = descriptor["id"]
    case = args.output / dataset
    case.mkdir()
    hubs_path = local_path(descriptor["hubs"])
    expression = local_path(descriptor["expression"])
    hubs = hubs_path.read_text(encoding="utf-8-sig").splitlines()
    if len(hubs) != descriptor["H"] or len(hubs) != len(set(hubs)):
        raise ValueError("Wrong authoritative source list")
    for name, path in (("expression.exp", expression), ("hubs.txt", hubs_path)):
        if sha256(path) != descriptor["input_records"][name]["sha256"]:
            raise ValueError("Frozen prepared input hash mismatch: " + str(path))
    frozen_receipt_path = args.relaxed / "results" / dataset / "complete.json"
    frozen_receipt = json.loads(frozen_receipt_path.read_text(encoding="utf-8"))
    result = {"dataset": dataset, "native_sha256": binary_sha,
              "expression_sha256": sha256(expression), "hubs_sha256": sha256(hubs_path),
              "H": len(hubs), "K": 20, "B": 100, "MI_jobs": 0,
              "source_raw_receipts": [], "comparisons": []}
    started = time.monotonic()
    # Only this owned directory is retired. The original raw files never change.
    with tempfile.TemporaryDirectory(prefix="canonical_", dir=str(case)) as temporary:
        canonical_dir = Path(temporary)
        canonical = []
        for seed, source, receipt, receipt_sha in raw_inputs(dataset, args.comparison, args.classic_raw):
            destination = canonical_dir / ("seed_%03d.adj" % seed)
            info = canonicalize(source, destination, hubs)
            info.update(seed=seed, receipt=str(receipt), receipt_sha256=receipt_sha)
            result["source_raw_receipts"].append(info)
            canonical.append(destination)
        manifest = canonical_dir / "inputs.txt"
        manifest.write_text("".join(str(path) + "\n" for path in canonical), encoding="utf-8")
        result["canonical_manifest_sha256"] = sha256(manifest)
        for method, offset in (("new", 0.0), ("new_relaxed", 0.025)):
            expected = next(row for row in rows if row["dataset"] == dataset and row["method"] == method)
            # The original N1/N3 final graphs used six-digit output. The saved
            # full-precision baseline replay preserves the same edge IDs and
            # candidate mean weights, so use that pinned graph for both checks.
            stage = "baseline_replay" if method == "new" else "relaxed"
            expected_path = args.relaxed / "results" / dataset / stage / "network.adj"
            expected_hash = frozen_receipt["baseline_replay_sha256" if method == "new" else "trial_sha256"]
            if sha256(expected_path) != expected_hash:
                raise ValueError("Frozen reference graph hash mismatch")
            output = case / (method + ".adj")
            command = [str(args.native), "-i", str(expression), "-s", str(hubs_path),
                       "-l", str(hubs_path), "--consensus-inputs", str(manifest),
                       "--min-recurrence", "20", "--dynamic-dpi", "--dpi-href", "479",
                       "--dpi-slope", "0.05", "--dpi-max", "0.30", "--dpi-offset", str(offset),
                       "-o", str(output)]
            tick = time.monotonic()
            with (case / (method + ".stdout.txt")).open("x") as stdout, (case / (method + ".stderr.txt")).open("x") as stderr:
                process = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
            if process.returncode:
                raise RuntimeError("Native consensus failed: " + str(case / (method + ".stderr.txt")))
            actual_edges, headers = load_graph(output)
            expected_edges, _ = load_graph(expected_path)
            header_values = {}
            for line in headers:
                for name in ("Consensus source rows H", "DPI source rows H", "DPI epsilon binary64"):
                    prefix = ">  " + name + " "
                    if line.startswith(prefix):
                        if name in header_values:
                            raise ValueError("Duplicate native metadata")
                        header_values[name] = line[len(prefix):].strip()
            if (int(header_values["Consensus source rows H"]) != len(hubs)
                    or int(header_values["DPI source rows H"]) != len(hubs)
                    or float(header_values["DPI epsilon binary64"]) != descriptor["epsilon_base" if method == "new" else "epsilon_trial"]):
                raise AssertionError("Resolved H/tolerance differs from the frozen schedule")
            if actual_edges != expected_edges:
                missing = len(expected_edges.keys() - actual_edges.keys())
                extra = len(actual_edges.keys() - expected_edges.keys())
                weights = sum(actual_edges[key] != expected_edges[key] for key in actual_edges.keys() & expected_edges.keys())
                raise AssertionError("%s %s mismatch: missing=%s extra=%s weights=%s" % (dataset, method, missing, extra, weights))
            result["comparisons"].append({"method": method, "offset": offset, "argv": command,
                "returncode": process.returncode, "elapsed_seconds": time.monotonic()-tick,
                "expected_path": str(expected_path), "expected_sha256": sha256(expected_path),
                "actual_path": str(output), "actual_sha256": sha256(output),
                "directed_edges": len(actual_edges), "edge_IDs_exact": True, "MI_binary64_numeric_exact": True,
                "metric_origin": "inherited frozen metrics after exact graph equivalence; not independently rescored",
                "resolved_H_and_epsilon_exact": True,
                "headers": headers, "TP": int(expected["TP"]), "FP": int(expected["FP"]),
                "precision": float(expected["precision"]), "recall": float(expected["recall"]),
                "f1": float(expected["f1"])})
    result["elapsed_seconds"] = time.monotonic()-started
    result["canonical_temporary_files_retired"] = True
    write_json(case / "validation.json", result)
    return result


def main():
    """Validate saved datasets in a fresh output directory and write receipts."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--relaxed", type=Path, required=True)
    parser.add_argument("--classic-raw", type=Path, required=True)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--datasets", nargs="*")
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    args.output.mkdir(parents=True, exist_ok=False)
    descriptors = json.loads((args.relaxed / "manifest.json").read_text(encoding="utf-8"))["datasets"]
    if args.datasets:
        requested = set(args.datasets)
        descriptors = [item for item in descriptors if item["id"] in requested]
        if {item["id"] for item in descriptors} != requested:
            parser.error("Unknown dataset")
    rows = read_tsv(args.relaxed / "full_comparison.tsv")
    binary_sha = sha256(args.native)
    result = {"schema": "native-consensus-saved-raw-validation-v1", "complete": False,
              "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "driver_sha256": sha256(__file__), "native_sha256": binary_sha,
              "comparison_manifest_sha256": sha256(args.comparison / "manifest.json"),
              "relaxed_manifest_sha256": sha256(args.relaxed / "manifest.json"),
              "MI_jobs": 0, "datasets": [], "errors": []}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_dataset, descriptor, args, rows, binary_sha): descriptor["id"] for descriptor in descriptors}
        for future in as_completed(futures):
            dataset = futures[future]
            try:
                receipt = future.result()
                result["datasets"].append(receipt)
                print(dataset + ": both offsets exact", flush=True)
            except Exception as error:
                result["errors"].append({"dataset": dataset, "type": type(error).__name__, "message": str(error)})
                print(dataset + ": " + str(error), flush=True)
    if sha256(args.native) != binary_sha:
        result["errors"].append({"message": "Native binary changed during validation"})
    result["datasets"].sort(key=lambda item: item["dataset"])
    result["complete"] = not result["errors"] and len(result["datasets"]) == len(descriptors)
    result["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    write_json(args.output / "validation.json", result)
    if not result["complete"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
