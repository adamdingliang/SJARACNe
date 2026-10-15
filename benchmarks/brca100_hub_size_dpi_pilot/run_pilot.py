#!/usr/bin/env python3
"""Run the matched, three-point BRCA100 hub-size DPI screening pilot.

The six inference arms change only the TF/SIG hub list used for both ``-s``
and ``-l``.  Each arm uses seeds 1..100, fixed 80-of-100 sampling, Npar=40,
DPI epsilon 0, and the class-specific AP-MI operating point.  A machine
``[DPI_STATS]`` line is mandatory for every completed seed.

The K>=6 aggregation in this benchmark is deliberately independent and
provisional.  It verifies direct recurrence counts; it is not the production
minimum-recurrence implementation planned for SJARACNe.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed


THIS_DIR = Path(__file__).resolve().parent
NETBID_BENCHMARK_DIR = THIS_DIR.parent / "brca100_netbid_qc"
sys.path.insert(0, str(NETBID_BENCHMARK_DIR))
import run_workflows as core  # noqa: E402

from pilot_common import (  # noqa: E402
    DRIVER_BY_KEY,
    DRIVERS,
    EXPRESSION_FILENAME,
    FRACTIONS,
    K_MINIMUM_RECURRENCE,
    SEEDS,
    arm_key,
    atomic_bytes,
    atomic_json,
    create_panel_files,
    fingerprint,
    iter_adjacency_edges,
    load_json,
    median,
    parse_dpi_stats,
    parse_sampling_indices,
    read_nonempty_unique_ids,
    sha256_file,
    utc_now,
)


MODEL_RELATIVE = "SJARACNe/config/apmi_null/apmi_null_m00080_npar040.model"
MODEL_EXPECTED_SHA256 = (
    "e3a8522682a8ea239821aaa10b12db72d00e07bfdcad43599d8e76a06be80944"
)
SCHEMA = "sjaracne-brca100-hub-size-dpi-pilot-v1"


def repo_git_command(
    repo: Path, *arguments: str, normalize_worktree: bool = True
) -> list[str]:
    """Use a linked worktree's own index even when .git has a Windows path."""

    dot_git = repo / ".git"
    if dot_git.is_file():
        declaration = dot_git.read_text(encoding="utf-8").strip()
        if not declaration.startswith("gitdir: "):
            raise ValueError(f"Malformed linked-worktree pointer: {dot_git}")
        git_dir_text = declaration[len("gitdir: ") :].replace("\\", "/")
        if os.name != "nt" and len(git_dir_text) >= 3 and git_dir_text[1:3] == ":/":
            git_dir_text = (
                f"/mnt/{git_dir_text[0].lower()}/{git_dir_text[3:]}"
            )
        git_dir = Path(git_dir_text)
        if not git_dir.is_absolute():
            git_dir = (repo / git_dir).resolve()
        command = [
            "git",
            "-c",
            f"safe.directory={repo}",
        ]
        if normalize_worktree:
            command.extend(
                ["-c", "core.autocrlf=true", "-c", "core.filemode=false"]
            )
        command.extend(
            [f"--git-dir={git_dir}", f"--work-tree={repo}", *arguments]
        )
        return command
    return core.git_command(repo, *arguments)


def resolve_commit(repo: Path, specification: str) -> str:
    completed = subprocess.run(
        repo_git_command(repo, "rev-parse", "--verify", f"{specification}^{{commit}}"),
        capture_output=True,
        text=True,
        check=True,
    )
    commit = completed.stdout.strip()
    if len(commit) != 40:
        raise RuntimeError(f"Unexpected commit resolution for {specification}: {commit}")
    return commit


def git_bytes(repo: Path, *arguments: str) -> bytes:
    return subprocess.run(
        repo_git_command(repo, *arguments), capture_output=True, check=True
    ).stdout


def stage_inputs(repo: Path, work_root: Path) -> dict[str, dict[str, object]]:
    """Stage canonical LF Git blobs while supporting Windows-created worktrees."""

    input_root = work_root / "inputs"
    input_root.mkdir(parents=True, exist_ok=True)
    metadata: dict[str, dict[str, object]] = {}
    baseline = core.STAGES[0].commit
    for filename in core.INPUT_FILES:
        destination = input_root / filename
        expected = core.EXPECTED_LF_SHA256[filename]
        if not destination.is_file() or sha256_file(destination) != expected:
            temporary = destination.with_name(destination.name + ".partial")
            if temporary.exists():
                temporary.unlink()
            with temporary.open("wb") as handle:
                subprocess.run(
                    repo_git_command(
                        repo,
                        "show",
                        f"{baseline}:tests/inputs/{filename}",
                        normalize_worktree=False,
                    ),
                    stdout=handle,
                    check=True,
                )
            actual = sha256_file(temporary)
            if actual != expected:
                temporary.unlink()
                raise RuntimeError(
                    f"Unexpected canonical input SHA256 for {filename}: {actual}"
                )
            os.replace(temporary, destination)
        metadata[filename] = {
            "path": str(destination),
            "sha256": expected,
            "bytes": destination.stat().st_size,
        }
    expression_ids = core.parse_expression_ids(input_root / EXPRESSION_FILENAME)
    if len(expression_ids) != 28278:
        raise ValueError(f"Expected 28,278 expression IDs, got {len(expression_ids)}")
    for driver in DRIVERS:
        ids = read_nonempty_unique_ids(
            input_root / driver.filename, expected_count=driver.full_count
        )
        missing = set(ids) - expression_ids
        if missing:
            raise ValueError(f"{driver.key} contains {len(missing)} missing expression IDs")
    metadata["expression_id_count"] = {"count": len(expression_ids)}
    return metadata


def working_source_provenance(repo: Path, base_commit: str) -> dict[str, object]:
    """Fingerprint tracked and untracked, non-ignored files under SJARACNe/."""

    status = git_bytes(
        repo,
        "status",
        "--porcelain=v1",
        "--untracked-files=all",
        "--",
        "SJARACNe",
    )
    diff = git_bytes(repo, "diff", "--binary", base_commit, "--", "SJARACNe")
    changed = git_bytes(
        repo, "diff", "--name-only", base_commit, "--", "SJARACNe"
    ).decode("utf-8").splitlines()
    untracked = git_bytes(
        repo,
        "ls-files",
        "--others",
        "--exclude-standard",
        "--",
        "SJARACNe",
    ).decode("utf-8").splitlines()
    listed = git_bytes(
        repo,
        "ls-files",
        "--cached",
        "--others",
        "--exclude-standard",
        "--",
        "SJARACNe",
    ).decode("utf-8").splitlines()
    entries: list[dict[str, object]] = []
    has_stats_marker = False
    for relative_text in sorted(set(listed)):
        relative = Path(relative_text)
        path = repo / relative
        if not path.is_file():
            entries.append({"path": relative.as_posix(), "exists": False})
            continue
        digest = sha256_file(path)
        entries.append(
            {
                "path": relative.as_posix(),
                "exists": True,
                "bytes": path.stat().st_size,
                "sha256": digest,
            }
        )
        if b"[DPI_STATS]" in path.read_bytes():
            has_stats_marker = True
    payload = {
        "base_commit": base_commit,
        "git_status_sha256": hashlib.sha256(status).hexdigest(),
        "git_status": status.decode("utf-8").splitlines(),
        "tracked_diff_sha256": hashlib.sha256(diff).hexdigest(),
        "tracked_diff_bytes": len(diff),
        "overlay_files": [
            next(entry for entry in entries if entry["path"] == Path(path).as_posix())
            for path in sorted(set(changed + untracked))
        ],
        "source_files": entries,
    }
    payload["source_tree_fingerprint"] = fingerprint(payload)
    payload["dirty"] = bool(status or diff)
    payload["has_dpi_stats_marker"] = has_stats_marker
    return payload


def extract_or_snapshot_and_build(
    repo: Path, work_root: Path, base_commit: str
) -> dict:
    """Build a commit archive or a frozen dirty-worktree source snapshot."""

    provenance = working_source_provenance(repo, base_commit)
    if not provenance["has_dpi_stats_marker"]:
        raise RuntimeError("Selected source tree does not contain [DPI_STATS] instrumentation")
    snapshot_id = str(provenance["source_tree_fingerprint"])
    source_mode = (
        "frozen-worktree-snapshot" if provenance["dirty"] else "commit-archive"
    )
    build_prefix = "worktree" if provenance["dirty"] else "commit"
    build_root = work_root / "builds" / f"{build_prefix}_{snapshot_id[:12]}"
    source_root = build_root / "source"
    binary = build_root / "bin" / "sjaracne.exe"
    manifest_path = build_root / "build_manifest.json"
    provenance_path = build_root / "source_provenance.json"
    diff_path = build_root / "tracked_source.patch"
    if manifest_path.is_file() and binary.is_file():
        manifest = load_json(manifest_path)
        if (
            manifest.get("source_provenance") == provenance
            and manifest.get("binary_sha256") == sha256_file(binary)
            and provenance_path.is_file()
            and load_json(provenance_path) == provenance
        ):
            return manifest
        raise RuntimeError(f"Stale or inconsistent worktree build: {build_root}")
    if source_root.exists() or binary.exists():
        raise RuntimeError(f"Incomplete worktree build: {build_root}")

    temporary_source = build_root / "source.partial"
    temporary_source.mkdir(parents=True)
    archive_path = build_root / "base_source.tar"
    try:
        core.checked_run(
            repo_git_command(
                repo,
                "archive",
                "--format=tar",
                "-o",
                str(archive_path),
                base_commit,
                "SJARACNe",
                normalize_worktree=False,
            )
        )
        with tarfile.open(archive_path, "r") as archive:
            archive.extractall(temporary_source, filter="data")
        for entry in provenance["overlay_files"]:
            relative = Path(str(entry["path"]))
            source = repo / relative
            destination = temporary_source / relative
            if entry.get("exists") and source.is_file():
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, destination)
            elif destination.exists():
                destination.unlink()
        snapshot_source_sha256 = core.sha256_directory(
            temporary_source / "SJARACNe"
        )
        os.replace(temporary_source, source_root)
    finally:
        if archive_path.exists():
            archive_path.unlink()
        if temporary_source.exists():
            shutil.rmtree(temporary_source)

    build_log = build_root / "build.stdout.log"
    build_error = build_root / "build.stderr.log"
    core.checked_run(
        ["make", "-C", str(source_root / "SJARACNe"), "-j8", "bin/sjaracne.exe"],
        stdout_path=build_log,
        stderr_path=build_error,
    )
    built_binary = source_root / "SJARACNe" / "bin" / "sjaracne.exe"
    if not built_binary.is_file():
        raise RuntimeError(f"Worktree build did not create {built_binary}")
    binary.parent.mkdir(parents=True)
    shutil.copy2(built_binary, binary)
    binary.chmod(0o755)
    compiler = subprocess.run(
        ["g++", "--version"], capture_output=True, text=True, check=True
    ).stdout.splitlines()[0]
    config_root = source_root / "SJARACNe" / "config"
    model = source_root / MODEL_RELATIVE
    manifest = {
        "stage": f"{build_prefix}_{snapshot_id[:12]}",
        "commit": base_commit,
        "source_mode": source_mode,
        "source_provenance": provenance,
        "snapshot_source_sha256": snapshot_source_sha256,
        "binary": str(binary),
        "binary_sha256": sha256_file(binary),
        "compiler": compiler,
        "config_directory": str(config_root),
        "config_sha256": core.sha256_directory(config_root),
        "null_model": str(model),
        "null_model_sha256": sha256_file(model),
        "built_at_utc": utc_now(),
    }
    atomic_json(provenance_path, provenance)
    atomic_bytes(
        diff_path,
        git_bytes(repo, "diff", "--binary", base_commit, "--", "SJARACNe"),
    )
    atomic_json(manifest_path, manifest)
    return manifest


def build_stage(commit: str) -> core.Stage:
    return core.Stage(
        key=f"hub_dpi_{commit[:12]}",
        commit=commit,
        sampling_args=("-u", "80%"),
        required_headers=(),
        null_model_relative=MODEL_RELATIVE,
    )


def validation_stage(commit: str, driver_key: str) -> core.Stage:
    driver = DRIVER_BY_KEY[driver_key]
    return core.Stage(
        key=arm_key(driver, driver.full_count),
        commit=commit,
        sampling_args=("-u", "80%"),
        required_headers=(
            ">  MI threshold method estimator-matched AP-MI permutation-null GPD tail",
            ">  AP-MI null model m 80",
            ">  AP-MI null model Npar 40",
            ">  AP-MI cutoff tail extrapolated no",
            ">  Sampling method fixed-size without replacement",
            ">  Sampling request 80%",
            ">  Eligible observations 100",
            ">  Sampled observations 80",
            f">  MI threshold    {format(driver.mi_cutoff, '.6g')}",
            f">  MI P-value      {format(driver.p_value, '.6g')}",
            ">  DPI tolerance   0",
        ),
        null_model_relative=MODEL_RELATIVE,
    )


def parse_model(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise ValueError(f"Malformed model line {path}:{line_number}")
        key, value = line.split("=", 1)
        if not key or key in result:
            raise ValueError(f"Invalid/duplicate model key {key!r} in {path}")
        result[key] = value
    return result


def gpd_cutoff(model: dict[str, str], probability: float) -> float:
    threshold = float(model["tail_threshold"])
    tail_probability = float(model["tail_probability"])
    shape = float(model["tail_shape"])
    scale = float(model["tail_scale"])
    if not 0.0 < probability <= tail_probability:
        raise ValueError(f"p={probability} is outside fitted model tail")
    if abs(shape) < 1e-12:
        return threshold - scale * math.log(probability / tail_probability)
    return threshold + (scale / shape) * (
        (probability / tail_probability) ** (-shape) - 1.0
    )


def verify_model(build: dict) -> None:
    if build.get("null_model_sha256") != MODEL_EXPECTED_SHA256:
        raise RuntimeError("Unexpected AP-MI null-model SHA256")
    model = parse_model(Path(build["null_model"]))
    if model.get("m") != "80" or model.get("npar_limit") != "40":
        raise RuntimeError("Pilot requires the exact m=80, Npar=40 null model")
    for driver in DRIVERS:
        observed = gpd_cutoff(model, driver.p_value)
        if not math.isclose(observed, driver.mi_cutoff, rel_tol=0.0, abs_tol=1e-14):
            raise RuntimeError(
                f"{driver.key} cutoff drift: expected {driver.mi_cutoff}, got {observed}"
            )


def selected_drivers(specification: str) -> list:
    if specification == "all":
        return list(DRIVERS)
    requested = [item.strip() for item in specification.split(",") if item.strip()]
    unknown = sorted(set(requested) - set(DRIVER_BY_KEY))
    if unknown or not requested or len(set(requested)) != len(requested):
        raise ValueError(f"Invalid/duplicate driver selection: {specification}")
    return [DRIVER_BY_KEY[key] for key in requested]


def selected_counts(driver, specification: str) -> list[int]:
    if specification == "all":
        return list(driver.counts)
    requested = [int(item.strip()) for item in specification.split(",") if item.strip()]
    unknown = sorted(set(requested) - set(driver.counts))
    if unknown or not requested or len(set(requested)) != len(requested):
        raise ValueError(
            f"Invalid/duplicate hub counts for {driver.key}: {specification}"
        )
    return requested


def command_for_job(
    *, build: dict, driver, panel: Path, seed: int, output: Path
) -> list[str]:
    return [
        build["binary"],
        "-i",
        str(panel.parents[3] / "inputs" / EXPRESSION_FILENAME),
        "-l",
        str(panel),
        "-s",
        str(panel),
        "-p",
        driver.p_token,
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
        "-v",
        "on",
        "-o",
        str(output),
        "-u",
        "80%",
        "-M",
        build["null_model"],
    ]


def arm_manifest(
    *, work_root: Path, source_commit: str, build: dict, panel_manifest: dict, driver, count: int
) -> dict[str, object]:
    panel_record = next(
        item
        for item in panel_manifest["drivers"][driver.key]["panels"]
        if int(item["hub_count"]) == count
    )
    return {
        "schema": "sjaracne-brca100-hub-size-dpi-arm-v1",
        "arm": arm_key(driver, count),
        "source_commit": source_commit,
        "source_mode": build["source_mode"],
        "source_tree_fingerprint": build["source_provenance"]["source_tree_fingerprint"],
        "binary_sha256": build["binary_sha256"],
        "config_sha256": build["config_sha256"],
        "null_model_sha256": build["null_model_sha256"],
        "driver": driver.key,
        "hub_count": count,
        "hub_fraction": count / driver.full_count,
        "panel": panel_record,
        "source_and_annotation_lists": "same subset (-s=subset, -l=subset)",
        "p_value": driver.p_value,
        "mi_cutoff": driver.mi_cutoff,
        "sampling": "fixed 80% without replacement",
        "m": 80,
        "npar": 40,
        "dpi_epsilon": 0,
        "minimum_recurrence": K_MINIMUM_RECURRENCE,
        "consensus_implementation": "benchmark-only direct support count",
        "seeds": list(SEEDS),
        "expression_sha256": sha256_file(work_root / "inputs" / EXPRESSION_FILENAME),
    }


def ensure_exact_json(path: Path, payload: dict) -> None:
    if path.is_file():
        if load_json(path) != payload:
            raise RuntimeError(f"Existing manifest differs from fixed design: {path}")
        return
    atomic_json(path, payload)


def prepare(
    *, repo: Path, work_root: Path, source_commit: str, drivers: list, counts_spec: str
) -> tuple[dict, dict, dict]:
    input_metadata = stage_inputs(repo, work_root)
    panel_manifest = create_panel_files(work_root / "inputs", work_root / "panels")
    build = extract_or_snapshot_and_build(repo, work_root, source_commit)
    if build.get("commit") != source_commit:
        raise RuntimeError("Build manifest commit mismatch")
    verify_model(build)

    arms = []
    # The design is always the complete fixed six-arm pilot.  CLI selectors
    # only control the current invocation and cannot mutate its provenance.
    for driver in DRIVERS:
        for count in driver.counts:
            manifest = arm_manifest(
                work_root=work_root,
                source_commit=source_commit,
                build=build,
                panel_manifest=panel_manifest,
                driver=driver,
                count=count,
            )
            root = work_root / "results" / arm_key(driver, count)
            ensure_exact_json(root / "arm_manifest.json", manifest)
            arms.append(manifest)

    harness_hashes = {
        path.name: sha256_file(path)
        for path in (THIS_DIR / "pilot_common.py", THIS_DIR / "run_pilot.py")
    }
    analysis_path = THIS_DIR / "analyze_pilot.py"
    if analysis_path.is_file():
        harness_hashes[analysis_path.name] = sha256_file(analysis_path)
    design = {
        "schema": SCHEMA,
        "source_commit": source_commit,
        "source_mode": build["source_mode"],
        "source_tree_fingerprint": build["source_provenance"]["source_tree_fingerprint"],
        "build": {
            key: build[key]
            for key in (
                "stage",
                "commit",
                "source_mode",
                "binary_sha256",
                "config_sha256",
                "null_model_sha256",
                "compiler",
            )
        },
        "source_provenance": build["source_provenance"],
        "inputs": input_metadata,
        "panel_manifest_sha256": sha256_file(work_root / "panels" / "panel_manifest.json"),
        "harness_sha256": harness_hashes,
        "arms": arms,
        "fixed_parameters": {
            "expression_gene_count": 28278,
            "sampling": "fixed 80% without replacement",
            "m": 80,
            "npar": 40,
            "dpi_epsilon": 0,
            "minimum_recurrence": K_MINIMUM_RECURRENCE,
            "seeds": list(SEEDS),
        },
    }
    ensure_exact_json(work_root / "pilot_design.json", design)
    return build, panel_manifest, design


def anchor_path(anchor_root: Path, driver_key: str, seed: int) -> Path:
    point = "p1e-03" if driver_key == "tf" else "p5e-04"
    return (
        anchor_root
        / "results"
        / point
        / driver_key
        / "adjacency"
        / f"TF_run_{seed:03d}.adj"
    )


def validate_output(
    *, path: Path, stdout_path: Path, source_commit: str, driver, panel_ids: set[str], expression_ids: set[str]
) -> tuple[dict, dict, dict]:
    adjacency = core.validate_adjacency(
        path,
        stage=validation_stage(source_commit, driver.key),
        driver_ids=panel_ids,
        expression_ids=expression_ids,
    )
    dpi = parse_dpi_stats(stdout_path, require_applied=True)
    sampling = parse_sampling_indices(stdout_path)
    if adjacency["edges"] != dpi["post_edges"]:
        raise ValueError(
            f"Adjacency/DPI mismatch for {path}: {adjacency['edges']} != {dpi['post_edges']}"
        )
    return adjacency, dpi, sampling


def run_seed_job(
    *,
    work_root: Path,
    build: dict,
    source_commit: str,
    driver,
    count: int,
    seed: int,
    expression_ids: set[str],
    panel_ids: set[str],
    anchor_root: Path | None,
) -> tuple[str, bool, dict]:
    arm = arm_key(driver, count)
    arm_root = work_root / "results" / arm
    panel = work_root / "panels" / driver.key / f"h{count:05d}" / driver.filename
    adjacency_root = arm_root / "adjacency"
    log_root = arm_root / "logs"
    marker_root = arm_root / "seed_metadata"
    for path in (adjacency_root, log_root, marker_root, arm_root / "work"):
        path.mkdir(parents=True, exist_ok=True)
    stem = f"TF_run_{seed:03d}"
    output = adjacency_root / f"{stem}.adj"
    marker = marker_root / f"{stem}.json"
    stdout_path = log_root / f"{stem}.stdout.log"
    stderr_path = log_root / f"{stem}.stderr.log"
    time_path = log_root / f"{stem}.time.txt"
    preview = command_for_job(
        build=build, driver=driver, panel=panel, seed=seed, output=output
    )
    fp_payload = {
        "schema": "sjaracne-brca100-hub-size-dpi-seed-v1",
        "source_commit": source_commit,
        "source_mode": build["source_mode"],
        "source_tree_fingerprint": build["source_provenance"]["source_tree_fingerprint"],
        "binary_sha256": build["binary_sha256"],
        "config_sha256": build["config_sha256"],
        "null_model_sha256": build["null_model_sha256"],
        "arm": arm,
        "driver": driver.key,
        "hub_count": count,
        "panel_sha256": sha256_file(panel),
        "expression_sha256": sha256_file(work_root / "inputs" / EXPRESSION_FILENAME),
        "seed": seed,
        "command_without_output": [
            "<OUTPUT>" if value == str(output) else value for value in preview
        ],
    }
    run_fingerprint = fingerprint(fp_payload)

    if marker.is_file() and output.is_file():
        existing = load_json(marker)
        adjacency, dpi, sampling = validate_output(
            path=output,
            stdout_path=stdout_path,
            source_commit=source_commit,
            driver=driver,
            panel_ids=panel_ids,
            expression_ids=expression_ids,
        )
        if (
            existing.get("fingerprint") == run_fingerprint
            and existing.get("adjacency", {}).get("full_sha256") == adjacency["full_sha256"]
            and existing.get("dpi") == dpi
            and existing.get("sampling") == sampling
        ):
            return f"{arm}/{seed:03d}", True, existing
        raise RuntimeError(f"Stale or inconsistent completed seed: {marker}")
    if marker.is_file() and not output.is_file():
        raise RuntimeError(f"Seed marker exists without adjacency: {marker}")
    if output.is_file() and not marker.is_file():
        raise RuntimeError(
            f"Orphan adjacency requires manual inspection (marker missing): {output}"
        )

    partial = arm_root / "work" / f"{stem}.adj.partial"
    if partial.exists():
        partial.unlink()
    command = command_for_job(
        build=build, driver=driver, panel=panel, seed=seed, output=partial
    )
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
        core.checked_run(
            timed_command,
            cwd=arm_root,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
        )
    except Exception as error:
        raise RuntimeError(f"Seed failed ({arm}/{seed}); see {stderr_path}") from error
    high_resolution_wall = time.perf_counter() - start_clock
    adjacency, dpi, sampling = validate_output(
        path=partial,
        stdout_path=stdout_path,
        source_commit=source_commit,
        driver=driver,
        panel_ids=panel_ids,
        expression_ids=expression_ids,
    )
    os.replace(partial, output)
    adjacency["full_sha256"] = sha256_file(output)

    anchor = None
    if count == driver.full_count and anchor_root is not None:
        expected_path = anchor_path(anchor_root, driver.key, seed)
        if not expected_path.is_file():
            raise RuntimeError(f"Missing full-size anchor: {expected_path}")
        expected = core.validate_adjacency(
            expected_path,
            stage=validation_stage(source_commit, driver.key),
            driver_ids=panel_ids,
            expression_ids=expression_ids,
        )
        anchor = {
            "path": str(expected_path),
            "sha256": expected["full_sha256"],
            "data_sha256": expected["data_sha256"],
            "data_match": expected["data_sha256"] == adjacency["data_sha256"],
            "edge_count_match": expected["edges"] == adjacency["edges"],
        }
        if not anchor["data_match"]:
            raise RuntimeError(f"Full-size anchor data mismatch: {arm}/{seed}")

    record = {
        "schema": "sjaracne-brca100-hub-size-dpi-seed-v1",
        "fingerprint": run_fingerprint,
        "source_commit": source_commit,
        "source_mode": build["source_mode"],
        "source_tree_fingerprint": build["source_provenance"]["source_tree_fingerprint"],
        "arm": arm,
        "driver": driver.key,
        "hub_count": count,
        "seed": seed,
        "command": command_for_job(
            build=build, driver=driver, panel=panel, seed=seed, output=output
        ),
        "binary_sha256": build["binary_sha256"],
        "config_sha256": build["config_sha256"],
        "null_model_sha256": build["null_model_sha256"],
        "started_at_utc": started,
        "finished_at_utc": utc_now(),
        "high_resolution_wall_s": high_resolution_wall,
        "gnu_time": core.parse_gnu_time(time_path),
        "stdout_sha256": sha256_file(stdout_path),
        "stderr_sha256": sha256_file(stderr_path),
        "stderr_bytes": stderr_path.stat().st_size,
        "adjacency": adjacency,
        "dpi": dpi,
        "sampling": sampling,
        "anchor": anchor,
    }
    atomic_json(marker, record)
    return f"{arm}/{seed:03d}", False, record


def aggregate_run_manifest(work_root: Path) -> None:
    rows: list[dict[str, object]] = []
    for marker in sorted((work_root / "results").glob("*/seed_metadata/TF_run_*.json")):
        record = load_json(marker)
        timing = record["gnu_time"]
        adjacency = record["adjacency"]
        dpi = record["dpi"]
        sampling = record["sampling"]
        rows.append(
            {
                "arm": record["arm"],
                "driver": record["driver"],
                "hub_count": record["hub_count"],
                "seed": record["seed"],
                "source_commit": record["source_commit"],
                "binary_sha256": record["binary_sha256"],
                "elapsed_s": timing["elapsed_s"],
                "user_s": timing["user_s"],
                "system_s": timing["system_s"],
                "max_rss_kib": timing["max_rss_kib"],
                "pre_edges": dpi["pre_edges"],
                "pruned_edges": dpi["pruned_edges"],
                "post_edges": dpi["post_edges"],
                "pruned_fraction": dpi["pruned_fraction"],
                "sampling_indices": ",".join(str(index) for index in sampling["indices"]),
                "sampling_sha256": sampling["sha256"],
                "source_rows": adjacency["source_rows"],
                "adjacency_bytes": adjacency["bytes"],
                "adjacency_sha256": adjacency["full_sha256"],
                "data_sha256": adjacency["data_sha256"],
                "anchor_data_match": (
                    "" if record.get("anchor") is None else record["anchor"]["data_match"]
                ),
            }
        )
    output = work_root / "results" / "run_manifest.tsv"
    temporary = output.with_name(output.name + ".partial")
    fieldnames = list(rows[0]) if rows else ["arm"]
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, output)


def validate_k6_output(path: Path, expected_rows: int) -> None:
    rows = 0
    seen: set[tuple[str, str]] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if header != ["source", "target", "support_count", "support_fraction"]:
            raise ValueError(f"Unexpected provisional K6 header in {path}")
        for line_number, line in enumerate(handle, 2):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 4:
                raise ValueError(f"Malformed K6 row {path}:{line_number}")
            edge = (fields[0], fields[1])
            support = int(fields[2])
            support_fraction = float(fields[3])
            if edge in seen or edge[0] == edge[1]:
                raise ValueError(f"Duplicate/self K6 edge at {path}:{line_number}")
            if not K_MINIMUM_RECURRENCE <= support <= len(SEEDS):
                raise ValueError(f"Invalid K6 support at {path}:{line_number}")
            if not math.isclose(support_fraction, support / len(SEEDS), abs_tol=1e-12):
                raise ValueError(f"Invalid K6 support fraction at {path}:{line_number}")
            seen.add(edge)
            rows += 1
    if rows != expected_rows:
        raise ValueError(f"K6 row mismatch in {path}: {rows} != {expected_rows}")


def aggregate_k6_arm(work_root: Path, driver, count: int) -> tuple[bool, dict]:
    """Disk-backed, direct K>=6 support aggregation for benchmark verification."""

    arm = arm_key(driver, count)
    arm_root = work_root / "results" / arm
    adjacencies = [
        arm_root / "adjacency" / f"TF_run_{seed:03d}.adj" for seed in SEEDS
    ]
    if any(not path.is_file() for path in adjacencies):
        missing = sum(not path.is_file() for path in adjacencies)
        raise RuntimeError(f"{arm} lacks {missing} of 100 adjacency inputs")
    adjacency_hashes = [sha256_file(path) for path in adjacencies]
    panel = work_root / "panels" / driver.key / f"h{count:05d}" / driver.filename
    panel_ids = read_nonempty_unique_ids(panel, expected_count=count)
    run_fingerprint = fingerprint(
        {
            "schema": "sjaracne-brca100-hub-size-dpi-provisional-k6-v1",
            "arm": arm,
            "k": K_MINIMUM_RECURRENCE,
            "seed_count": len(SEEDS),
            "panel_sha256": sha256_file(panel),
            "adjacency_sha256": adjacency_hashes,
        }
    )
    output_root = arm_root / "provisional_k6"
    output_root.mkdir(parents=True, exist_ok=True)
    output = output_root / "consensus_support_ge6.tsv"
    manifest_path = output_root / "manifest.json"
    if manifest_path.is_file() and output.is_file():
        existing = load_json(manifest_path)
        validate_k6_output(output, int(existing["k6_edges"]))
        if (
            existing.get("fingerprint") == run_fingerprint
            and existing.get("output_sha256") == sha256_file(output)
        ):
            return True, existing
        raise RuntimeError(f"Stale provisional K6 output: {output}")

    database = output_root / "support.sqlite.partial"
    if database.exists():
        database.unlink()
    connection = sqlite3.connect(database)
    try:
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute(
            "CREATE TABLE support (source TEXT NOT NULL, target TEXT NOT NULL, "
            "count INTEGER NOT NULL, PRIMARY KEY (source, target)) WITHOUT ROWID"
        )
        statement = (
            "INSERT INTO support(source,target,count) VALUES(?,?,1) "
            "ON CONFLICT(source,target) DO UPDATE SET count=count+1"
        )
        for path in adjacencies:
            connection.executemany(statement, iter_adjacency_edges(path))
            connection.commit()
        temporary = output.with_name(output.name + ".partial")
        target_counts = {accession: 0 for accession in panel_ids}
        k6_edges = 0
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("source\ttarget\tsupport_count\tsupport_fraction\n")
            cursor = connection.execute(
                "SELECT source,target,count FROM support WHERE count>=? "
                "ORDER BY source,target",
                (K_MINIMUM_RECURRENCE,),
            )
            for source, target, support in cursor:
                if source not in target_counts:
                    raise ValueError(f"K6 source {source} is outside panel {arm}")
                target_counts[source] += 1
                k6_edges += 1
                handle.write(
                    f"{source}\t{target}\t{support}\t{support / len(SEEDS):.12g}\n"
                )
        os.replace(temporary, output)
    finally:
        connection.close()
        if database.exists():
            database.unlink()

    active_hubs = sum(value > 0 for value in target_counts.values())
    manifest = {
        "schema": "sjaracne-brca100-hub-size-dpi-provisional-k6-v1",
        "warning": (
            "Benchmark-only direct recurrence aggregation; not the production "
            "SJARACNe minimum-recurrence implementation"
        ),
        "fingerprint": run_fingerprint,
        "arm": arm,
        "driver": driver.key,
        "hub_count": count,
        "k": K_MINIMUM_RECURRENCE,
        "seed_count": len(SEEDS),
        "k6_edges": k6_edges,
        "zero_filled_median_target_count": median(target_counts.values()),
        "active_hubs": active_hubs,
        "active_hub_fraction": active_hubs / count,
        "output_sha256": sha256_file(output),
        "output_bytes": output.stat().st_size,
        "adjacency_sha256": adjacency_hashes,
    }
    validate_k6_output(output, k6_edges)
    atomic_json(manifest_path, manifest)
    return False, manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("/mnt/d/GitHub/SJARACNe-hub-dpi"))
    parser.add_argument(
        "--work-root",
        type=Path,
        default=Path.home() / "sjaracne-benchmarks" / "brca100-hub-size-dpi-pilot",
    )
    parser.add_argument("--source-commit", default="HEAD")
    parser.add_argument(
        "--phase", choices=("prepare", "infer", "aggregate", "all"), default="all"
    )
    parser.add_argument("--drivers", default="all")
    parser.add_argument("--hub-counts", default="all")
    parser.add_argument("--seed-start", type=int, default=1)
    parser.add_argument("--seed-end", type=int, default=100)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument(
        "--anchor-root",
        type=Path,
        help="Existing threshold-sweep work root used to verify 100%% arms",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Prepare the frozen design and print planned inference commands only",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not 1 <= args.seed_start <= args.seed_end <= 100:
        raise ValueError("Seed range must be within 1..100")
    if args.workers < 1:
        raise ValueError("--workers must be positive")
    args.work_root.mkdir(parents=True, exist_ok=True)
    drivers = selected_drivers(args.drivers)
    source_commit = resolve_commit(args.repo, args.source_commit)
    build, panel_manifest, design = prepare(
        repo=args.repo,
        work_root=args.work_root,
        source_commit=source_commit,
        drivers=drivers,
        counts_spec=args.hub_counts,
    )
    core.console(f"[PREPARE] frozen source commit {source_commit}")

    seeds = list(range(args.seed_start, args.seed_end + 1))
    tasks = [
        (driver, count, seed)
        for seed in seeds
        for driver in drivers
        for count in selected_counts(driver, args.hub_counts)
    ]
    invocation_path = args.work_root / "invocations.json"
    invocations = (
        load_json(invocation_path)
        if invocation_path.is_file()
        else {"schema": "sjaracne-brca100-hub-size-dpi-invocations-v1", "invocations": []}
    )
    invocation = {
        "started_at_utc": utc_now(),
        "status": "running",
        "phase": args.phase,
        "source_commit": source_commit,
        "arms": sorted({arm_key(driver, count) for driver, count, _ in tasks}),
        "seed_start": args.seed_start,
        "seed_end": args.seed_end,
        "workers": args.workers,
        "anchor_root": None if args.anchor_root is None else str(args.anchor_root),
        "dry_run": args.dry_run,
    }
    invocations["invocations"].append(invocation)
    atomic_json(invocation_path, invocations)

    if args.dry_run:
        for driver, count, seed in tasks:
            panel = (
                args.work_root / "panels" / driver.key / f"h{count:05d}" / driver.filename
            )
            output = (
                args.work_root
                / "results"
                / arm_key(driver, count)
                / "adjacency"
                / f"TF_run_{seed:03d}.adj"
            )
            print(" ".join(command_for_job(build=build, driver=driver, panel=panel, seed=seed, output=output)))
        invocation.update(status="complete", finished_at_utc=utc_now(), planned_jobs=len(tasks))
        atomic_json(invocation_path, invocations)
        return 0

    if args.phase in ("infer", "all"):
        expression_ids = core.parse_expression_ids(
            args.work_root / "inputs" / EXPRESSION_FILENAME
        )
        panel_ids = {
            (driver.key, count): set(
                read_nonempty_unique_ids(
                    args.work_root
                    / "panels"
                    / driver.key
                    / f"h{count:05d}"
                    / driver.filename,
                    expected_count=count,
                )
            )
            for driver in drivers
            for count in selected_counts(driver, args.hub_counts)
        }
        core.console(f"[INFER] {len(tasks)} matched jobs with {args.workers} workers")
        completed = 0
        resumed = 0
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(
                    run_seed_job,
                    work_root=args.work_root,
                    build=build,
                    source_commit=source_commit,
                    driver=driver,
                    count=count,
                    seed=seed,
                    expression_ids=expression_ids,
                    panel_ids=panel_ids[(driver.key, count)],
                    anchor_root=args.anchor_root,
                ): (driver, count, seed)
                for driver, count, seed in tasks
            }
            for future in as_completed(futures):
                label, was_resumed, record = future.result()
                completed += 1
                resumed += int(was_resumed)
                core.console(
                    f"[INFER {completed}/{len(tasks)}] {label} "
                    f"{'resume' if was_resumed else 'done'}; "
                    f"pre={record['dpi']['pre_edges']} "
                    f"pruned={record['dpi']['pruned_edges']} "
                    f"post={record['dpi']['post_edges']}"
                )
        aggregate_run_manifest(args.work_root)
        invocation.update(
            inference_jobs=len(tasks),
            inference_resumed_jobs=resumed,
            inference_new_jobs=len(tasks) - resumed,
        )

    if args.phase in ("aggregate", "all"):
        if seeds != list(SEEDS):
            raise ValueError("Provisional K6 aggregation requires seeds 1..100")
        records = []
        for driver in drivers:
            for count in selected_counts(driver, args.hub_counts):
                resumed, record = aggregate_k6_arm(args.work_root, driver, count)
                records.append(record)
                core.console(
                    f"[K6] {record['arm']} {'resume' if resumed else 'done'}; "
                    f"edges={record['k6_edges']}"
                )
        atomic_json(
            args.work_root / "results" / "provisional_k6_manifest.json",
            {
                "schema": "sjaracne-brca100-hub-size-dpi-provisional-k6-aggregate-v1",
                "warning": "Benchmark-only; not the production consensus implementation",
                "records": records,
            },
        )

    invocation.update(status="complete", finished_at_utc=utc_now())
    atomic_json(invocation_path, invocations)
    core.console(f"[DONE] {args.phase} at {args.work_root}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
