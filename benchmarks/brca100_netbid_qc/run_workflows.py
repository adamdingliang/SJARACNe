#!/usr/bin/env python3
"""Run the matched BRCA100 baseline, PR66, and PR67 inference workflows.

The runner deliberately invokes frozen native binaries instead of CWL so that
the 600 seed jobs are bounded, resumable, and individually checksummed.  It
reproduces the workflow's seeds (1..100), native arguments, and consensus code.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed


@dataclasses.dataclass(frozen=True)
class Stage:
    key: str
    commit: str
    sampling_args: tuple[str, ...]
    required_headers: tuple[str, ...]
    null_model_relative: str | None = None


@dataclasses.dataclass(frozen=True)
class DriverClass:
    key: str
    filename: str
    expected_count: int


STAGES = (
    Stage(
        key="baseline_12113fb",
        commit="12113fbc80d753d945598ffc2c7d9e45787bc8e0",
        sampling_args=("-r", "1"),
        required_headers=(
            ">  MI P-value      1e-07",
            ">  MI threshold    0.153257",
        ),
    ),
    Stage(
        key="pr66_5809183",
        commit="58091832848b2eaf2ae08f6f69482357b6b9b72c",
        sampling_args=("-u", "80%"),
        required_headers=(
            ">  Sampling method fixed-size without replacement",
            ">  Sampling request 80%",
            ">  Eligible observations 100",
            ">  Sampled observations 80",
            ">  MI threshold    0.172803",
        ),
    ),
    Stage(
        key="pr67_7633ebb",
        commit="7633ebb4a0d966dbda15a4e32d0efa492fb71aeb",
        sampling_args=("-u", "80%"),
        null_model_relative=(
            "SJARACNe/config/apmi_null/apmi_null_m00080_npar040.model"
        ),
        required_headers=(
            ">  MI threshold method estimator-matched AP-MI permutation-null GPD tail",
            ">  AP-MI null model m 80",
            ">  AP-MI null model Npar 40",
            ">  AP-MI cutoff tail extrapolated yes",
            ">  Sampling method fixed-size without replacement",
            ">  Sampling request 80%",
            ">  Eligible observations 100",
            ">  Sampled observations 80",
            ">  MI threshold    0.322465",
        ),
    ),
)

DRIVER_CLASSES = (
    DriverClass("tf", "BRCA100_TF.txt", 2608),
    DriverClass("sig", "BRCA100_SIG.txt", 10680),
)

# The stock consensus implementation retains the union of all seed-level
# edges in Python dictionaries.  Process the sparsest arms first and leave the
# largest baseline/SIG union until last so a full run proves the low-memory
# path before reaching the peak-memory arm.
CONSENSUS_ARM_ORDER = (
    ("pr67_7633ebb", "tf"),
    ("pr67_7633ebb", "sig"),
    ("pr66_5809183", "tf"),
    ("pr66_5809183", "sig"),
    ("baseline_12113fb", "tf"),
    ("baseline_12113fb", "sig"),
)

_EXPECTED_CONSENSUS_ARMS = {
    (stage.key, driver.key) for stage in STAGES for driver in DRIVER_CLASSES
}
if (
    len(CONSENSUS_ARM_ORDER) != len(_EXPECTED_CONSENSUS_ARMS)
    or set(CONSENSUS_ARM_ORDER) != _EXPECTED_CONSENSUS_ARMS
):
    raise RuntimeError("CONSENSUS_ARM_ORDER must cover every stage/driver arm once")

INPUT_FILES = (
    "BRCA100.exp",
    "BRCA100_TF.txt",
    "BRCA100_SIG.txt",
)

EXPECTED_LF_SHA256 = {
    "BRCA100.exp": "ad8a334f5f8cdf46a1000d3ee259b35258a18b3da2e314bb3a0cf7a421d98bc8",
    "BRCA100_TF.txt": "9b1219a489b99432175e4c4ad46add7b06f25aae388ee8dd3261fa91e4c43ffd",
    "BRCA100_SIG.txt": "16ca27df655f16684f880a4ad719c4e2ae3f8dc0d7e6b9eccdd24cd97c40797c",
}

CONSOLE_LOCK = threading.Lock()


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def console(message: str) -> None:
    with CONSOLE_LOCK:
        print(message, flush=True)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_directory(path: Path) -> str:
    digest = hashlib.sha256()
    for child in sorted(item for item in path.rglob("*") if item.is_file()):
        digest.update(child.relative_to(path).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(bytes.fromhex(sha256_file(child)))
    return digest.hexdigest()


def json_fingerprint(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    os.replace(temporary, path)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object in {path}")
    return value


def checked_run(
    command: list[str],
    *,
    cwd: Path | None = None,
    stdout_path: Path | None = None,
    stderr_path: Path | None = None,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess:
    stdout_handle = (
        stdout_path.open("w", encoding="utf-8", newline="\n")
        if stdout_path is not None
        else None
    )
    stderr_handle = (
        stderr_path.open("w", encoding="utf-8", newline="\n")
        if stderr_path is not None
        else None
    )
    try:
        return subprocess.run(
            command,
            cwd=str(cwd) if cwd is not None else None,
            env=env,
            stdout=stdout_handle,
            stderr=stderr_handle,
            check=True,
            text=True,
        )
    finally:
        if stdout_handle is not None:
            stdout_handle.close()
        if stderr_handle is not None:
            stderr_handle.close()


def git_command(repo: Path, *arguments: str) -> list[str]:
    return [
        "git",
        "-c",
        f"safe.directory={repo}",
        "-C",
        str(repo),
        *arguments,
    ]


def stage_git_inputs(repo: Path, work_root: Path) -> dict[str, dict[str, object]]:
    input_root = work_root / "inputs"
    input_root.mkdir(parents=True, exist_ok=True)
    metadata: dict[str, dict[str, object]] = {}
    baseline = STAGES[0].commit

    for filename in INPUT_FILES:
        destination = input_root / filename
        expected = EXPECTED_LF_SHA256[filename]
        if not destination.is_file() or sha256_file(destination) != expected:
            temporary = destination.with_name(destination.name + ".partial")
            if temporary.exists():
                temporary.unlink()
            with temporary.open("wb") as handle:
                subprocess.run(
                    git_command(
                        repo,
                        "show",
                        f"{baseline}:tests/inputs/{filename}",
                    ),
                    stdout=handle,
                    check=True,
                )
            actual = sha256_file(temporary)
            if actual != expected:
                temporary.unlink()
                raise RuntimeError(
                    f"Unexpected LF input SHA256 for {filename}: {actual}"
                )
            os.replace(temporary, destination)

        metadata[filename] = {
            "path": str(destination),
            "sha256": expected,
            "bytes": destination.stat().st_size,
        }

    expression_ids: set[str] = set()
    expression_path = input_root / "BRCA100.exp"
    with expression_path.open("r", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header[:2] != ["isoformId", "geneSymbol"]:
            raise ValueError(f"Unexpected expression header in {expression_path}")
        for line_number, line in enumerate(handle, 2):
            accession = line.split("\t", 1)[0]
            if not accession or accession in expression_ids:
                raise ValueError(
                    f"Invalid/duplicate expression accession at line {line_number}"
                )
            expression_ids.add(accession)
    if len(expression_ids) != 28278:
        raise ValueError(f"Expected 28,278 expression IDs, got {len(expression_ids)}")

    for driver in DRIVER_CLASSES:
        path = input_root / driver.filename
        ids = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
        if len(ids) != driver.expected_count or len(set(ids)) != len(ids):
            raise ValueError(
                f"Unexpected {driver.key} driver count/duplicates: {len(ids)}"
            )
        missing = set(ids) - expression_ids
        if missing:
            raise ValueError(f"{driver.key} contains {len(missing)} missing IDs")

    metadata["expression_id_count"] = {"count": len(expression_ids)}
    return metadata


def extract_and_build_stage(repo: Path, work_root: Path, stage: Stage) -> dict:
    build_root = work_root / "builds" / stage.key
    source_root = build_root / "source"
    binary = build_root / "bin" / "sjaracne.exe"
    manifest_path = build_root / "build_manifest.json"

    if manifest_path.is_file() and binary.is_file():
        manifest = load_json(manifest_path)
        config_root = source_root / "SJARACNe" / "config"
        config_matches = (
            config_root.is_dir()
            and manifest.get("config_sha256") == sha256_directory(config_root)
        )
        model_matches = True
        if stage.null_model_relative:
            model = source_root / stage.null_model_relative
            model_matches = (
                model.is_file()
                and manifest.get("null_model_sha256") == sha256_file(model)
            )
        if (
            manifest.get("commit") == stage.commit
            and manifest.get("binary_sha256") == sha256_file(binary)
            and config_matches
            and model_matches
        ):
            return manifest
        raise RuntimeError(f"Stale or inconsistent build at {build_root}")

    if source_root.exists() or binary.exists():
        raise RuntimeError(
            f"Incomplete existing build at {build_root}; remove that exact stage directory"
        )

    build_root.mkdir(parents=True, exist_ok=True)
    archive_path = build_root / "source.tar"
    temporary_source = build_root / "source.partial"
    temporary_source.mkdir()
    try:
        checked_run(
            git_command(repo, "archive", "--format=tar", "-o", str(archive_path), stage.commit)
        )
        with tarfile.open(archive_path, "r") as archive:
            archive.extractall(temporary_source, filter="data")
        os.replace(temporary_source, source_root)
    finally:
        if archive_path.exists():
            archive_path.unlink()

    build_log = build_root / "build.stdout.log"
    build_error = build_root / "build.stderr.log"
    checked_run(
        ["make", "-C", str(source_root / "SJARACNe"), "-j8", "bin/sjaracne.exe"],
        stdout_path=build_log,
        stderr_path=build_error,
    )

    built_binary = source_root / "SJARACNe" / "bin" / "sjaracne.exe"
    if not built_binary.is_file():
        raise RuntimeError(f"Build did not create {built_binary}")
    binary.parent.mkdir(parents=True)
    shutil.copy2(built_binary, binary)
    binary.chmod(0o755)

    compiler = subprocess.run(
        ["g++", "--version"], capture_output=True, text=True, check=True
    ).stdout.splitlines()[0]
    config_root = source_root / "SJARACNe" / "config"
    manifest = {
        "stage": stage.key,
        "commit": stage.commit,
        "binary": str(binary),
        "binary_sha256": sha256_file(binary),
        "compiler": compiler,
        "config_directory": str(config_root),
        "config_sha256": sha256_directory(config_root),
        "built_at_utc": utc_now(),
    }
    if stage.null_model_relative:
        model = source_root / stage.null_model_relative
        if not model.is_file():
            raise RuntimeError(f"Missing PR67 null model {model}")
        manifest["null_model"] = str(model)
        manifest["null_model_sha256"] = sha256_file(model)
    atomic_json(manifest_path, manifest)
    return manifest


def parse_expression_ids(expression_path: Path) -> set[str]:
    ids: set[str] = set()
    with expression_path.open("r", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header[:2] != ["isoformId", "geneSymbol"]:
            raise ValueError(f"Unexpected expression header in {expression_path}")
        for line in handle:
            ids.add(line.split("\t", 1)[0])
    return ids


def parse_driver_ids(driver_path: Path) -> set[str]:
    return {
        line.strip()
        for line in driver_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }


def validate_adjacency(
    path: Path,
    *,
    stage: Stage,
    driver_ids: set[str],
    expression_ids: set[str],
) -> dict[str, object]:
    headers: list[str] = []
    sources: set[str] = set()
    edges: set[tuple[str, str]] = set()
    data_digest = hashlib.sha256()
    mi_min = math.inf
    mi_max = -math.inf

    with path.open("rb") as raw:
        for line_number, raw_line in enumerate(raw, 1):
            if raw_line.startswith(b">"):
                headers.append(raw_line.decode("utf-8").rstrip("\r\n"))
                continue
            data_digest.update(raw_line)
            line = raw_line.decode("utf-8").rstrip("\r\n")
            if not line:
                raise ValueError(f"Blank data row in {path}:{line_number}")
            tokens = line.split("\t")
            if len(tokens) < 3 or len(tokens) % 2 == 0:
                raise ValueError(f"Malformed adjacency row in {path}:{line_number}")
            source = tokens[0]
            if source in sources:
                raise ValueError(f"Duplicate source row {source} in {path}")
            if source not in driver_ids:
                raise ValueError(f"Non-driver source {source} in {path}")
            sources.add(source)
            for index in range(1, len(tokens), 2):
                target = tokens[index]
                if target not in expression_ids:
                    raise ValueError(f"Unknown target {target} in {path}")
                if source == target:
                    raise ValueError(f"Self edge {source} in {path}")
                edge = (source, target)
                if edge in edges:
                    raise ValueError(f"Duplicate edge {source}->{target} in {path}")
                edges.add(edge)
                mi = float(tokens[index + 1])
                if not math.isfinite(mi) or mi <= 0.0:
                    raise ValueError(f"Invalid MI {mi} for {source}->{target} in {path}")
                mi_min = min(mi_min, mi)
                mi_max = max(mi_max, mi)

    header_text = "\n".join(headers)
    for expected in stage.required_headers:
        if expected not in header_text:
            raise ValueError(f"Missing header {expected!r} in {path}")
    if not headers:
        raise ValueError(f"No adjacency headers in {path}")

    return {
        "full_sha256": sha256_file(path),
        "data_sha256": data_digest.hexdigest(),
        "bytes": path.stat().st_size,
        "header_lines": len(headers),
        "source_rows": len(sources),
        "edges": len(edges),
        "mi_min": None if not edges else mi_min,
        "mi_max": None if not edges else mi_max,
    }


def parse_gnu_time(path: Path) -> dict[str, float | int]:
    values: dict[str, float | int] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        key, value = line.split("=", 1)
        values[key] = int(value) if key == "max_rss_kib" else float(value)
    return values


def command_for_job(
    stage: Stage,
    build: dict,
    driver: DriverClass,
    input_root: Path,
    seed: int,
    output: Path,
) -> list[str]:
    command = [
        build["binary"],
        "-i",
        str(input_root / "BRCA100.exp"),
        "-l",
        str(input_root / driver.filename),
        "-s",
        str(input_root / driver.filename),
        "-p",
        "0.0000001",
        "-e",
        "0",
        "-a",
        "adaptive_partitioning",
        "-H",
        build["config_directory"].rstrip("/") + "/",
        "-N",
        "40",
        "-S",
        str(seed),
        "-o",
        str(output),
        *stage.sampling_args,
    ]
    if stage.null_model_relative:
        command.extend(["-M", build["null_model"]])
    return command


def run_seed_job(
    *,
    stage: Stage,
    build: dict,
    driver: DriverClass,
    input_root: Path,
    results_root: Path,
    seed: int,
    expression_ids: set[str],
    driver_ids: set[str],
) -> tuple[str, bool, dict]:
    arm_root = results_root / stage.key / driver.key
    adjacency_root = arm_root / "adjacency"
    log_root = arm_root / "logs"
    metadata_root = arm_root / "seed_metadata"
    adjacency_root.mkdir(parents=True, exist_ok=True)
    log_root.mkdir(parents=True, exist_ok=True)
    metadata_root.mkdir(parents=True, exist_ok=True)

    stem = f"TF_run_{seed:03d}"
    output = adjacency_root / f"{stem}.adj"
    marker = metadata_root / f"{stem}.json"
    command_preview = command_for_job(
        stage, build, driver, input_root, seed, output
    )
    fingerprint_payload = {
        "stage": stage.key,
        "commit": stage.commit,
        "binary_sha256": build["binary_sha256"],
        "driver": driver.key,
        "driver_sha256": sha256_file(input_root / driver.filename),
        "expression_sha256": sha256_file(input_root / "BRCA100.exp"),
        "seed": seed,
        "command_without_output": [
            value if value != str(output) else "<OUTPUT>" for value in command_preview
        ],
    }
    legacy_fingerprint = json_fingerprint(fingerprint_payload)
    fingerprint_payload["config_sha256"] = build["config_sha256"]
    fingerprint_payload["null_model_sha256"] = build.get("null_model_sha256")
    fingerprint = json_fingerprint(fingerprint_payload)

    if marker.is_file() and output.is_file():
        existing = load_json(marker)
        if existing.get("fingerprint") in (fingerprint, legacy_fingerprint):
            stats = validate_adjacency(
                output,
                stage=stage,
                driver_ids=driver_ids,
                expression_ids=expression_ids,
            )
            if existing.get("adjacency", {}).get("full_sha256") == stats["full_sha256"]:
                if existing.get("fingerprint") == legacy_fingerprint:
                    existing["fingerprint"] = fingerprint
                    existing["config_sha256"] = build["config_sha256"]
                    existing["null_model_sha256"] = build.get("null_model_sha256")
                    atomic_json(marker, existing)
                return f"{stage.key}/{driver.key}/{seed:03d}", True, existing
        raise RuntimeError(f"Stale or inconsistent completed seed {marker}")

    partial = arm_root / "work" / f"{stem}.adj.partial"
    partial.parent.mkdir(parents=True, exist_ok=True)
    if partial.exists():
        partial.unlink()
    stdout_path = log_root / f"{stem}.stdout.log"
    stderr_path = log_root / f"{stem}.stderr.log"
    time_path = log_root / f"{stem}.time.txt"
    command = command_for_job(stage, build, driver, input_root, seed, partial)
    timed_command = [
        "/usr/bin/time",
        "-f",
        "elapsed_s=%e\nuser_s=%U\nsystem_s=%S\nmax_rss_kib=%M",
        "-o",
        str(time_path),
        *command,
    ]

    started = utc_now()
    start_clock = time.perf_counter()
    try:
        checked_run(
            timed_command,
            cwd=arm_root,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
        )
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            f"Seed failed ({stage.key}/{driver.key}/{seed}): {error.returncode}; "
            f"see {stderr_path}"
        ) from error
    high_resolution_wall = time.perf_counter() - start_clock
    stats = validate_adjacency(
        partial,
        stage=stage,
        driver_ids=driver_ids,
        expression_ids=expression_ids,
    )
    os.replace(partial, output)
    stats["full_sha256"] = sha256_file(output)

    record = {
        "fingerprint": fingerprint,
        "stage": stage.key,
        "commit": stage.commit,
        "driver": driver.key,
        "seed": seed,
        "command": command_for_job(stage, build, driver, input_root, seed, output),
        "binary_sha256": build["binary_sha256"],
        "config_sha256": build["config_sha256"],
        "null_model_sha256": build.get("null_model_sha256"),
        "started_at_utc": started,
        "finished_at_utc": utc_now(),
        "high_resolution_wall_s": high_resolution_wall,
        "gnu_time": parse_gnu_time(time_path),
        "stdout_sha256": sha256_file(stdout_path),
        "stderr_sha256": sha256_file(stderr_path),
        "stderr_bytes": stderr_path.stat().st_size,
        "adjacency": stats,
    }
    atomic_json(marker, record)
    return f"{stage.key}/{driver.key}/{seed:03d}", False, record


def aggregate_seed_manifest(results_root: Path, selected_stages: list[Stage], selected_drivers: list[DriverClass]) -> None:
    rows: list[dict[str, object]] = []
    for stage in selected_stages:
        for driver in selected_drivers:
            metadata_root = results_root / stage.key / driver.key / "seed_metadata"
            for marker in sorted(metadata_root.glob("TF_run_*.json")):
                record = load_json(marker)
                adjacency = record["adjacency"]
                timing = record["gnu_time"]
                rows.append(
                    {
                        "stage": stage.key,
                        "commit": stage.commit,
                        "driver": driver.key,
                        "seed": record["seed"],
                        "binary_sha256": record["binary_sha256"],
                        "elapsed_s": timing["elapsed_s"],
                        "user_s": timing["user_s"],
                        "system_s": timing["system_s"],
                        "max_rss_kib": timing["max_rss_kib"],
                        "edges": adjacency["edges"],
                        "source_rows": adjacency["source_rows"],
                        "adjacency_bytes": adjacency["bytes"],
                        "adjacency_sha256": adjacency["full_sha256"],
                        "data_sha256": adjacency["data_sha256"],
                        "stderr_bytes": record["stderr_bytes"],
                    }
                )

    output = results_root / "run_manifest.tsv"
    temporary = output.with_name(output.name + ".partial")
    fieldnames = list(rows[0]) if rows else []
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, output)


def validate_consensus_ncol(
    path: Path,
    *,
    driver_ids: set[str],
    expression_ids: set[str],
) -> dict[str, object]:
    expected_header = [
        "source",
        "target",
        "source.symbol",
        "target.symbol",
        "MI",
        "pearson",
        "spearman",
        "slope",
        "p-value",
    ]
    edges: set[tuple[str, str]] = set()
    sources: set[str] = set()
    nodes: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != expected_header:
            raise ValueError(f"Unexpected consensus header in {path}: {reader.fieldnames}")
        for row_number, row in enumerate(reader, 2):
            source = row["source"]
            target = row["target"]
            edge = (source, target)
            if source not in driver_ids or source not in expression_ids:
                raise ValueError(f"Invalid source at {path}:{row_number}")
            if target not in expression_ids or source == target or edge in edges:
                raise ValueError(f"Invalid/duplicate target at {path}:{row_number}")
            mi = float(row["MI"])
            if not math.isfinite(mi) or mi <= 0.0:
                raise ValueError(f"Invalid consensus MI at {path}:{row_number}")
            edges.add(edge)
            sources.add(source)
            nodes.update(edge)
    if not edges:
        raise ValueError(f"Consensus network is empty: {path}")
    return {
        "sha256": sha256_file(path),
        "bytes": path.stat().st_size,
        "edges": len(edges),
        "active_drivers": len(sources),
        "incident_nodes": len(nodes),
    }


def run_consensus(
    *,
    benchmark_repo: Path,
    work_root: Path,
    stage: Stage,
    driver: DriverClass,
    seeds: list[int],
    expression_ids: set[str],
    driver_ids: set[str],
) -> dict:
    if seeds != list(range(1, 101)):
        raise ValueError("Consensus requires the complete ordered seed range 1..100")
    results_root = work_root / "results"
    arm_root = results_root / stage.key / driver.key
    adjacency_root = arm_root / "adjacency"
    files = sorted(adjacency_root.glob("TF_run_*.adj"))
    if len(files) != 100:
        raise ValueError(f"Expected 100 adjacency files in {adjacency_root}, got {len(files)}")
    if any(item.suffix != ".adj" for item in adjacency_root.iterdir()):
        raise ValueError(f"Non-adjacency entry found in {adjacency_root}")

    adjacency_hashes = [sha256_file(path) for path in files]
    fingerprint = json_fingerprint(
        {
            "stage": stage.key,
            "driver": driver.key,
            "adjacency_hashes": adjacency_hashes,
            "consensus_p": 1e-5,
            "consensus_script_sha256": sha256_file(
                benchmark_repo / "SJARACNe" / "bin" / "create_consensus_network.py"
            ),
        }
    )
    consensus_root = arm_root / "consensus"
    marker_path = arm_root / "consensus_manifest.json"
    pending_marker_path = arm_root / "consensus_manifest.pending.json"
    ncol_path = consensus_root / "consensus_network_ncol_.txt"
    stdout_path = arm_root / "logs" / "consensus.stdout.log"
    stderr_path = arm_root / "logs" / "consensus.stderr.log"
    time_path = arm_root / "logs" / "consensus.time.txt"
    expression = work_root / "inputs" / "BRCA100.exp"
    command = [
        sys.executable,
        "-m",
        "SJARACNe.bin.create_consensus_network",
        "-a",
        str(adjacency_root),
        "-p",
        "0.00001",
        "-e",
        str(expression),
        "-o",
        str(arm_root / "consensus.partial"),
    ]

    def validated_consensus_record(root: Path) -> dict:
        stats = validate_consensus_ncol(
            root / "consensus_network_ncol_.txt",
            driver_ids=driver_ids,
            expression_ids=expression_ids,
        )
        parameter_info = root / "parameter_info_.txt"
        if ">  Bootstrap No: 100" not in parameter_info.read_text(encoding="utf-8"):
            raise ValueError(f"Consensus did not record 100 inputs: {parameter_info}")
        return {
            "fingerprint": fingerprint,
            "stage": stage.key,
            "driver": driver.key,
            "command": command,
            "finished_at_utc": utc_now(),
            "gnu_time": parse_gnu_time(time_path),
            "ncol": stats,
            "consensus_3col_sha256": sha256_file(
                root / "consensus_network_3col_.txt"
            ),
            "parameter_info_sha256": sha256_file(parameter_info),
            "bootstrap_info_sha256": sha256_file(root / "bootstrap_info_.txt"),
        }

    if marker_path.is_file() and ncol_path.is_file():
        marker = load_json(marker_path)
        if marker.get("fingerprint") == fingerprint:
            stats = validate_consensus_ncol(
                ncol_path,
                driver_ids=driver_ids,
                expression_ids=expression_ids,
            )
            if stats["sha256"] == marker.get("ncol", {}).get("sha256"):
                if pending_marker_path.exists():
                    pending_marker_path.unlink()
                return marker
        raise RuntimeError(f"Stale or inconsistent consensus at {consensus_root}")
    if consensus_root.exists() and not marker_path.exists():
        if not pending_marker_path.is_file():
            raise RuntimeError(
                f"Unverifiable orphan consensus directory: {consensus_root}"
            )
        pending = load_json(pending_marker_path)
        actual = validated_consensus_record(consensus_root)
        comparable_keys = (
            "fingerprint",
            "ncol",
            "consensus_3col_sha256",
            "parameter_info_sha256",
            "bootstrap_info_sha256",
        )
        if any(pending.get(key) != actual.get(key) for key in comparable_keys):
            raise RuntimeError(f"Orphan consensus fails pending manifest: {consensus_root}")
        pending["recovered_after_interrupted_manifest"] = True
        atomic_json(marker_path, pending)
        pending_marker_path.unlink()
        return pending

    temporary = arm_root / "consensus.partial"
    if consensus_root.exists():
        raise RuntimeError(f"Unexpected existing consensus path under {arm_root}")
    if temporary.exists():
        shutil.rmtree(temporary)
    if pending_marker_path.exists():
        pending_marker_path.unlink()
    temporary.mkdir()
    timed_command = [
        "/usr/bin/time",
        "-f",
        "elapsed_s=%e\nuser_s=%U\nsystem_s=%S\nmax_rss_kib=%M",
        "-o",
        str(time_path),
        *command,
    ]
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(benchmark_repo)
    checked_run(
        timed_command,
        cwd=benchmark_repo,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        env=environment,
    )
    record = validated_consensus_record(temporary)
    atomic_json(pending_marker_path, record)
    os.replace(temporary, consensus_root)
    atomic_json(marker_path, record)
    pending_marker_path.unlink()
    return record


def select_by_key(values: tuple, requested: str, label: str) -> list:
    if requested == "all":
        return list(values)
    keys = [part.strip() for part in requested.split(",") if part.strip()]
    mapping = {value.key: value for value in values}
    unknown = set(keys) - set(mapping)
    if unknown:
        raise ValueError(f"Unknown {label}: {sorted(unknown)}")
    return [mapping[key] for key in keys]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo",
        type=Path,
        default=Path("/mnt/d/GitHub/SJARACNe"),
        help="Main Git checkout used for exact git archive/show operations",
    )
    parser.add_argument(
        "--benchmark-repo",
        type=Path,
        default=Path("/mnt/d/GitHub/SJARACNe-brca100-netbid-qc"),
        help="Benchmark worktree containing current consensus code",
    )
    parser.add_argument(
        "--work-root",
        type=Path,
        default=(
            Path.home()
            / "sjaracne-benchmarks"
            / "brca100-netbid-qc-20260817-rerun"
        ),
    )
    parser.add_argument(
        "--phase",
        choices=("build", "infer", "consensus", "all"),
        default="all",
    )
    parser.add_argument("--stages", default="all")
    parser.add_argument("--drivers", default="all")
    parser.add_argument("--seed-start", type=int, default=1)
    parser.add_argument("--seed-end", type=int, default=100)
    parser.add_argument("--workers", type=int, default=8)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not (1 <= args.seed_start <= args.seed_end <= 100):
        raise ValueError("Seed range must be within 1..100")
    if args.workers < 1:
        raise ValueError("--workers must be positive")
    selected_stages = select_by_key(STAGES, args.stages, "stage")
    selected_drivers = select_by_key(DRIVER_CLASSES, args.drivers, "driver")
    seeds = list(range(args.seed_start, args.seed_end + 1))

    args.work_root.mkdir(parents=True, exist_ok=True)
    input_metadata = stage_git_inputs(args.repo, args.work_root)
    expression_ids = parse_expression_ids(args.work_root / "inputs" / "BRCA100.exp")
    driver_id_map = {
        driver.key: parse_driver_ids(args.work_root / "inputs" / driver.filename)
        for driver in selected_drivers
    }
    builds: dict[str, dict] = {}

    if args.phase in ("build", "infer", "all"):
        for stage in selected_stages:
            console(f"[BUILD] {stage.key} {stage.commit}")
            builds[stage.key] = extract_and_build_stage(args.repo, args.work_root, stage)

    run_metadata = {
        "schema": "sjaracne-brca100-netbid-qc-v1",
        "generated_at_utc": utc_now(),
        "benchmark_repo": str(args.benchmark_repo),
        "work_root": str(args.work_root),
        "stages": [dataclasses.asdict(stage) for stage in selected_stages],
        "drivers": [dataclasses.asdict(driver) for driver in selected_drivers],
        "seeds": seeds,
        "workers": args.workers,
        "inputs": input_metadata,
        "builds": builds,
        "common_parameters": {
            "bootstrap_p": 1e-7,
            "consensus_p": 1e-5,
            "dpi_epsilon": 0,
            "npar": 40,
            "algorithm": "adaptive_partitioning",
        },
    }
    invocation = {
        "started_at_utc": utc_now(),
        "status": "running",
        "phase": args.phase,
        "stages": [stage.key for stage in selected_stages],
        "drivers": [driver.key for driver in selected_drivers],
        "seed_start": args.seed_start,
        "seed_end": args.seed_end,
        "workers": args.workers,
    }
    invocation_completion: dict[str, object] = {}
    metadata_path = args.work_root / "run_metadata.json"
    if metadata_path.is_file():
        existing_metadata = load_json(metadata_path)
        if (
            existing_metadata.get("schema") != run_metadata["schema"]
            or existing_metadata.get("inputs") != run_metadata["inputs"]
            or existing_metadata.get("common_parameters")
            != run_metadata["common_parameters"]
        ):
            raise RuntimeError(f"Incompatible existing run metadata: {metadata_path}")
        existing_metadata.setdefault("invocations", []).append(invocation)
        if args.phase in ("infer", "all"):
            existing_metadata["inference_workers"] = args.workers
        atomic_json(metadata_path, existing_metadata)
    else:
        run_metadata["invocations"] = [invocation]
        if args.phase in ("infer", "all"):
            run_metadata["inference_workers"] = args.workers
        atomic_json(metadata_path, run_metadata)

    if args.phase in ("infer", "all"):
        tasks = [
            (stage, driver, seed)
            for seed in seeds
            for driver in selected_drivers
            for stage in selected_stages
        ]
        completed = 0
        skipped = 0
        console(f"[INFER] {len(tasks)} jobs with {args.workers} workers")
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(
                    run_seed_job,
                    stage=stage,
                    build=builds[stage.key],
                    driver=driver,
                    input_root=args.work_root / "inputs",
                    results_root=args.work_root / "results",
                    seed=seed,
                    expression_ids=expression_ids,
                    driver_ids=driver_id_map[driver.key],
                ): (stage, driver, seed)
                for stage, driver, seed in tasks
            }
            for future in as_completed(futures):
                label, was_skipped, record = future.result()
                completed += 1
                skipped += int(was_skipped)
                console(
                    f"[INFER {completed}/{len(tasks)}] {label} "
                    f"{'resume' if was_skipped else 'done'}; "
                    f"edges={record['adjacency']['edges']} "
                    f"wall={record['gnu_time']['elapsed_s']:.2f}s"
                )
        aggregate_seed_manifest(
            args.work_root / "results", selected_stages, selected_drivers
        )
        invocation_completion.update(
            {
                "inference_jobs": len(tasks),
                "inference_resumed_jobs": skipped,
                "inference_new_jobs": len(tasks) - skipped,
            }
        )
        console(f"[INFER] complete; resumed={skipped}, new={len(tasks) - skipped}")

    if args.phase in ("consensus", "all"):
        if not builds:
            for stage in selected_stages:
                builds[stage.key] = extract_and_build_stage(args.repo, args.work_root, stage)
        selected_stage_map = {stage.key: stage for stage in selected_stages}
        selected_driver_map = {driver.key: driver for driver in selected_drivers}
        consensus_arms = [
            (selected_stage_map[stage_key], selected_driver_map[driver_key])
            for stage_key, driver_key in CONSENSUS_ARM_ORDER
            if stage_key in selected_stage_map and driver_key in selected_driver_map
        ]
        for stage, driver in consensus_arms:
            console(f"[CONSENSUS] {stage.key}/{driver.key}")
            record = run_consensus(
                benchmark_repo=args.benchmark_repo,
                work_root=args.work_root,
                stage=stage,
                driver=driver,
                seeds=seeds,
                expression_ids=expression_ids,
                driver_ids=driver_id_map[driver.key],
            )
            console(
                f"[CONSENSUS] {stage.key}/{driver.key} "
                f"edges={record['ncol']['edges']} "
                f"rss={record['gnu_time']['max_rss_kib']} KiB"
            )

    completed_metadata = load_json(metadata_path)
    for recorded_invocation in reversed(completed_metadata["invocations"]):
        if recorded_invocation.get("started_at_utc") == invocation["started_at_utc"]:
            recorded_invocation.update(invocation_completion)
            recorded_invocation["status"] = "complete"
            recorded_invocation["finished_at_utc"] = utc_now()
            break
    else:
        raise RuntimeError("Current invocation disappeared from run metadata")
    atomic_json(metadata_path, completed_metadata)

    console(f"[DONE] {args.phase} at {args.work_root}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        console(f"[ERROR] {type(error).__name__}: {error}")
        raise
