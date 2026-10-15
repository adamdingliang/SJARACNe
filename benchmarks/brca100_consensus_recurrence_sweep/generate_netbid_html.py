#!/usr/bin/env python3
"""Generate and package four provenance-locked NetBID2 HTML reports.

This is deliberately an adjunct to ``run_recurrence_sweep.py``.  It reads the
frozen recurrence work root, writes HTML-mode NetBID2 products to a separate
overlay, and never modifies the frozen design, networks, summary-mode QC, or
compact result package.
"""

from __future__ import annotations

import argparse
import csv
from contextlib import contextmanager
import datetime as dt
from dataclasses import dataclass
import errno
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
from typing import Any, Iterable
import uuid


DESIGN_SCHEMA = "sjaracne-brca100-consensus-recurrence-sweep-v1"
NETWORK_SCHEMA = "sjaracne-brca100-consensus-recurrence-network-v1"
SUMMARY_SCHEMA = "sjaracne-brca100-consensus-recurrence-netbid2-v1"
HTML_RUN_SCHEMA = "sjaracne-brca100-consensus-recurrence-netbid2-html-v1"
HTML_AGGREGATE_SCHEMA = (
    "sjaracne-brca100-consensus-recurrence-netbid2-html-aggregate-v1"
)
HTML_PACKAGE_SCHEMA = (
    "sjaracne-brca100-consensus-recurrence-netbid2-html-package-v1"
)
OVERLAY_LOCK_SCHEMA = (
    "sjaracne-brca100-consensus-recurrence-netbid2-html-lock-v1"
)
PARTIAL_OWNER_SCHEMA = (
    "sjaracne-brca100-consensus-recurrence-netbid2-html-partial-owner-v1"
)

EXPECTED_DESIGN_SHA256 = (
    "1228de2bd5f2ae0e59bfaadefc43d8a9a080c68e6ae7313713f1bba762269216"
)
EXPECTED_DESIGN_FINGERPRINT = (
    "8948b3d2be1a9e1a69591131b0475b1c5265475b16665caca34be4f72d0a6e8e"
)
EXPECTED_FROZEN_PACKAGE_SHA256 = (
    "cec702cc601e115ab23d9d59e842d0cdb70b01f6e218987ae092941d15009abf"
)
EXPECTED_FROZEN_PACKAGE_FINGERPRINT = (
    "14020b969a9769ee2db99b2b665680f3b07b791992f56336acf9e5a6caeec889"
)
EXPECTED_WRAPPER_SHA256 = (
    "c706c7665fc2e589bf70030bb5dc0fc21935b93922b16fb256a5801752220095"
)
EXPECTED_FROZEN_SCRIPT_HASHES = {
    "aggregator_source_sha256": (
        "12fd3a0f629920365691aedd191f1ee06171b37a9bbf7663f88b784063bec3f0"
    ),
    "netbid_r_sha256": (
        "f9485b923e278ce7883d359181819cd728402dff11fcd4172b7fd6c0c51149a4"
    ),
    "runner_sha256": (
        "c276044854b590fcfb55780633d9cc1b749011905ede207a69cfe42b2f1f6d9a"
    ),
}
EXPECTED_ENVIRONMENT = {
    "R": "R version 4.4.3 (2025-02-28)",
    "NetBID2": "2.2.0",
    "NetBID2_remote_sha": "5defa454d600b94f5dd6d1f9f4428f99759a6821",
    "igraph": "2.3.3",
}


@dataclass(frozen=True)
class DriverSpec:
    filename: str
    prefix: str
    candidate_drivers: int
    sha256: str
    source_point: str
    per_subsample_p: float


@dataclass(frozen=True)
class RunSpec:
    edges: int
    network_bytes: int
    network_sha256: str
    network_manifest_sha256: str
    summary_manifest_sha256: str
    summary_outputs: dict[str, str]


DRIVER_SPECS = {
    "tf": DriverSpec(
        "BRCA100_TF.txt",
        "TF_",
        2608,
        "9b1219a489b99432175e4c4ad46add7b06f25aae388ee8dd3261fa91e4c43ffd",
        "p1e-03",
        1e-3,
    ),
    "sig": DriverSpec(
        "BRCA100_SIG.txt",
        "SIG_",
        10680,
        "16ca27df655f16684f880a4ad719c4e2ae3f8dc0d7e6b9eccdd24cd97c40797c",
        "p5e-04",
        5e-4,
    ),
}

RUN_ORDER = (("tf", 6), ("sig", 6), ("tf", 8), ("sig", 8))
RUN_SPECS = {
    ("tf", 6): RunSpec(
        416408,
        31191254,
        "9e55caaabce7c126ab165965452379be52e1ad8a481fbb955044a648879e1e55",
        "1acb007a5d1f30889a050451c8ab55efedce47f62ac71c69301e2c6f80d6f5b7",
        "ff86e127bbb203e75be57628608c0d605bf06bf1690402715aa93889cd7ef973",
        {
            "driver_target_sizes.tsv": (
                "a821872b52eb5eb136c233ad151d5b604667e29a89b349c5e30d83ee01cb251d"
            ),
            "netbid_environment.tsv": (
                "589fd14d9f877353eaf0b97e8dbc5a0b5e09bb4248ef6ee0114d709279facb9b"
            ),
            "network_summary.tsv": (
                "4a950f9d13df5cd74c501dc07059ed038b31f17b65871c67788234c61adcaf64"
            ),
        },
    ),
    ("sig", 6): RunSpec(
        739958,
        55365716,
        "47fef1d0d151e002c405fd68fcdeaca0eb838eaee9651c3eb0fc09cead88ce98",
        "04c1dcaaef2d3a8bc36f42379eac89a597d604d258243e61810f3b129b50e40c",
        "91ad1de904ad8ed7dd3ed16c0ce89cfc3d810c1d294082207408924f38e70474",
        {
            "driver_target_sizes.tsv": (
                "1e22aafb1c218e31f59cdd5c461625ac92180e8084481e2911a44432bb103af1"
            ),
            "netbid_environment.tsv": (
                "589fd14d9f877353eaf0b97e8dbc5a0b5e09bb4248ef6ee0114d709279facb9b"
            ),
            "network_summary.tsv": (
                "192274f640aeefa5020b99b95b12b5f6e4cbc3dde3bc35d753139ef0fb33c6d0"
            ),
        },
    ),
    ("tf", 8): RunSpec(
        269294,
        20148472,
        "da5d4b8d4e3777f32780198903ab1f8315c6034d84b6b93cd1b57544d8711272",
        "e7e490d8a6e5331c2c1343b59dcd7c3e17adc9e728bf9ad1f66c24752a361266",
        "831294e4824c03ca2a9afba70bb4539edc8396b3bb7ff94dc61684ef77e26baa",
        {
            "driver_target_sizes.tsv": (
                "1cbb25d56fe041ff0a5351010f2ad18317e4b1bfad27aaeb0326526d46deb1f6"
            ),
            "netbid_environment.tsv": (
                "589fd14d9f877353eaf0b97e8dbc5a0b5e09bb4248ef6ee0114d709279facb9b"
            ),
            "network_summary.tsv": (
                "f1122950f1aaffd2458ec0c418010fbbfa89a74a11064c6b235eedc87572bff3"
            ),
        },
    ),
    ("sig", 8): RunSpec(
        462099,
        34529971,
        "91136399f1608321dd1a9c9158c508805f27459c460ee18bd1aed376e9243a8a",
        "9f7ed9366eb11e8d516e308cf85c70186e2020ba8cd78b930b3b533abacb0b75",
        "eb40a0c24bd5ba944e5ac6cea815c3506f3530d97989af798844aa25ebac42ea",
        {
            "driver_target_sizes.tsv": (
                "df5ae55c2fc572c218e75eb916ab5e07710ef2b31c934d0dd4186648c78203a8"
            ),
            "netbid_environment.tsv": (
                "589fd14d9f877353eaf0b97e8dbc5a0b5e09bb4248ef6ee0114d709279facb9b"
            ),
            "network_summary.tsv": (
                "b117d7a47b55ef5c5063fba16ada17dfdadce50622844fe65ae79f3bef459c80"
            ),
        },
    ),
}

REQUIRED_METRICS = {
    "candidate_drivers",
    "active_drivers",
    "active_driver_fraction",
    "edges",
    "incident_nodes",
    "weak_components",
    "largest_weak_component",
    "largest_weak_component_fraction",
    "density",
    "target_size_zero_mean",
    "target_size_zero_median",
    "target_size_zero_q25",
    "target_size_zero_q75",
    "target_size_zero_max",
    "target_size_active_mean",
    "target_size_active_median",
    "target_size_active_q25",
    "target_size_active_q75",
    "target_size_active_max",
    "scale_free_adjusted_r2",
}
SUMMARY_FILENAMES = {
    "driver_target_sizes.tsv",
    "netbid_environment.tsv",
    "network_summary.tsv",
}

_ACTIVE_OVERLAY_LOCKS: dict[Path, dict[str, Any]] = {}


@dataclass(frozen=True)
class Context:
    repo_root: Path
    source_root: Path
    recurrence_root: Path
    html_root: Path
    design_path: Path
    design_sha256: str
    design_fingerprint: str
    frozen_package_manifest: Path
    wrapper: Path
    r_script: Path
    generator: Path
    frozen_runner: Path
    aggregator_source: Path
    environment: dict[str, str]
    generator_sha256: str
    driver_ids: dict[str, tuple[str, ...]]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_fingerprint(value: object) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def canonical_record(payload: dict[str, Any]) -> dict[str, Any]:
    result = dict(payload)
    result["fingerprint"] = canonical_fingerprint(result)
    return result


def validate_record_fingerprint(record: dict[str, Any], source: Path) -> None:
    payload = dict(record)
    observed = payload.pop("fingerprint", None)
    expected = canonical_fingerprint(payload)
    if observed != expected:
        raise ValueError(f"Canonical fingerprint mismatch: {source}")


def serialized_json(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def atomic_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.partial.{os.getpid()}")
    try:
        with temporary.open("xb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_json(path: Path, value: object) -> None:
    atomic_bytes(path, serialized_json(value))


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def require_file(path: Path, description: str) -> Path:
    if not path.is_file() or path.is_symlink():
        raise FileNotFoundError(f"Missing/non-regular {description}: {path}")
    return path


def require_hash(path: Path, expected: str, description: str) -> None:
    actual = sha256_file(require_file(path, description))
    if actual != expected:
        raise ValueError(
            f"{description} SHA-256 mismatch: {path}: {actual} != {expected}"
        )


def is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def secure_subdirectory(base: Path, *parts: str) -> Path:
    """Create a directory path without following nested symlinks."""
    base.mkdir(parents=True, exist_ok=True)
    if base.is_symlink() or not base.is_dir():
        raise RuntimeError(f"Unsafe managed directory root: {base}")
    resolved_base = base.resolve(strict=True)
    current = base
    for part in parts:
        if not part or part in {".", ".."} or "/" in part or "\\" in part:
            raise ValueError(f"Unsafe managed directory component: {part!r}")
        candidate = current / part
        if candidate.is_symlink():
            raise RuntimeError(f"Symlink in managed directory path: {candidate}")
        if candidate.exists():
            if not candidate.is_dir():
                raise RuntimeError(f"Non-directory in managed path: {candidate}")
        else:
            candidate.mkdir()
        resolved_candidate = candidate.resolve(strict=True)
        if not is_within(resolved_candidate, resolved_base):
            raise RuntimeError(f"Managed directory escaped its root: {candidate}")
        current = candidate
    return current


def assert_managed_path(base: Path, path: Path) -> None:
    """Fail closed if a managed path or any existing ancestor is a symlink."""
    resolved_base = base.resolve(strict=True)
    try:
        relative = path.relative_to(base)
    except ValueError as exc:
        raise RuntimeError(f"Managed path is outside its root: {path}") from exc
    current = base
    for part in relative.parts[:-1]:
        current = current / part
        if current.is_symlink() or not current.is_dir():
            raise RuntimeError(f"Unsafe managed path ancestor: {current}")
    parent = path.parent.resolve(strict=True)
    if not is_within(parent, resolved_base):
        raise RuntimeError(f"Managed path resolved outside its root: {path}")
    if path.is_symlink():
        raise RuntimeError(f"Managed path is a symlink: {path}")


def managed_atomic_json(context: Context, path: Path, value: object) -> None:
    assert_managed_path(context.html_root, path)
    atomic_json(path, value)


def managed_unlink(context: Context, path: Path) -> None:
    assert_managed_path(context.html_root, path)
    path.unlink()


def managed_rename(context: Context, source: Path, destination: Path) -> None:
    assert_managed_path(context.html_root, source)
    assert_managed_path(context.html_root, destination)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"Managed rename destination exists: {destination}")
    os.rename(source, destination)


def pid_is_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError as exc:
        if exc.errno in {errno.ESRCH, errno.EINVAL} or getattr(exc, "winerror", None) == 87:
            return False
        raise
    return True


def lock_path(context: Context) -> Path:
    return context.html_root / ".generate_netbid_html.lock"


def validate_lock_record(record: dict[str, Any], source: Path) -> None:
    expected_fields = {
        "schema",
        "pid",
        "token",
        "design_fingerprint",
        "generator_sha256",
        "started_at_utc",
        "fingerprint",
    }
    if set(record) != expected_fields or record.get("schema") != OVERLAY_LOCK_SCHEMA:
        raise RuntimeError(f"Malformed HTML overlay lock: {source}")
    validate_record_fingerprint(record, source)
    if not isinstance(record.get("pid"), int) or not isinstance(record.get("token"), str):
        raise RuntimeError(f"Malformed HTML overlay lock owner: {source}")


@contextmanager
def overlay_lock(context: Context):
    secure_subdirectory(context.html_root)
    path = lock_path(context)
    assert_managed_path(context.html_root, path)
    if context.html_root in _ACTIVE_OVERLAY_LOCKS:
        raise RuntimeError(f"HTML overlay lock is already active in this process: {path}")
    if path.exists() or path.is_symlink():
        require_file(path, "HTML overlay lock")
        existing = load_json(path)
        validate_lock_record(existing, path)
        if pid_is_alive(int(existing["pid"])):
            raise RuntimeError(
                f"HTML overlay is owned by live PID {existing['pid']}: {path}"
            )
        managed_unlink(context, path)

    record = canonical_record(
        {
            "schema": OVERLAY_LOCK_SCHEMA,
            "pid": os.getpid(),
            "token": uuid.uuid4().hex,
            "design_fingerprint": context.design_fingerprint,
            "generator_sha256": context.generator_sha256,
            "started_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        }
    )
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(serialized_json(record))
            handle.flush()
            os.fsync(handle.fileno())
    except Exception:
        if path.exists() and not path.is_symlink():
            path.unlink()
        raise
    _ACTIVE_OVERLAY_LOCKS[context.html_root] = record
    try:
        yield record
    finally:
        active = _ACTIVE_OVERLAY_LOCKS.pop(context.html_root, None)
        if active != record:
            raise RuntimeError(f"HTML overlay lock ownership changed: {path}")
        require_file(path, "HTML overlay lock")
        if load_json(path) != record:
            raise RuntimeError(f"HTML overlay lock bytes changed: {path}")
        managed_unlink(context, path)


def active_overlay_lock(context: Context) -> dict[str, Any]:
    record = _ACTIVE_OVERLAY_LOCKS.get(context.html_root)
    if record is None:
        raise RuntimeError("An active HTML overlay lock is required")
    path = lock_path(context)
    require_file(path, "HTML overlay lock")
    if load_json(path) != record:
        raise RuntimeError(f"HTML overlay lock no longer matches its owner: {path}")
    return record


def validate_path_separation(
    repo_root: Path, source_root: Path, recurrence_root: Path, html_root: Path
) -> None:
    roots = {
        "repository": repo_root,
        "source": source_root,
        "recurrence": recurrence_root,
    }
    for label, root in roots.items():
        if html_root == root or is_within(html_root, root) or is_within(root, html_root):
            raise ValueError(
                f"HTML overlay must be disjoint from the {label} root: {html_root}"
            )


def parse_environment_table(lines: Iterable[str], source: str) -> dict[str, str]:
    rows = list(csv.DictReader(lines, delimiter="\t"))
    if not rows or any(set(row) != {"component", "version"} for row in rows):
        raise ValueError(f"Malformed NetBID2 environment table: {source}")
    result: dict[str, str] = {}
    for row in rows:
        component = row["component"]
        if not component or component in result:
            raise ValueError(f"Duplicate/empty environment component: {source}")
        result[component] = row["version"]
    if result != EXPECTED_ENVIRONMENT:
        raise ValueError(f"Pinned NetBID2 environment mismatch: {result!r}")
    return result


def probe_environment(wrapper: Path, r_script: Path) -> dict[str, str]:
    completed = subprocess.run(
        [str(wrapper), "Rscript", str(r_script), "--probe"],
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "NetBID2 environment probe failed:\n" + (completed.stderr or "")
        )
    return parse_environment_table(completed.stdout.splitlines(), "environment probe")


def read_driver_ids(path: Path, expected_count: int) -> tuple[str, ...]:
    values = tuple(
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    )
    if len(values) != expected_count or len(set(values)) != len(values):
        raise ValueError(f"Unexpected driver list count/content: {path}")
    return values


def validate_design(
    repo_root: Path, source_root: Path, recurrence_root: Path
) -> tuple[Path, dict[str, Any]]:
    path = recurrence_root / "design.json"
    require_hash(path, EXPECTED_DESIGN_SHA256, "frozen recurrence design")
    design = load_json(path)
    if design.get("schema") != DESIGN_SCHEMA:
        raise ValueError(f"Unexpected recurrence design schema: {path}")
    if design.get("fingerprint") != EXPECTED_DESIGN_FINGERPRINT:
        raise ValueError(f"Unexpected recurrence design fingerprint: {path}")
    fingerprint_payload = dict(design)
    observed = fingerprint_payload.pop("fingerprint")
    if canonical_fingerprint(fingerprint_payload) != observed:
        raise ValueError(f"Internally inconsistent recurrence design: {path}")
    if Path(str(design.get("source_work_root"))).resolve() != source_root:
        raise ValueError("Recurrence design/source work-root mismatch")
    if design.get("minimum_supports") != list(range(6, 21)):
        raise ValueError("Frozen recurrence support grid changed")
    if design.get("benchmark_scripts") != EXPECTED_FROZEN_SCRIPT_HASHES:
        raise ValueError("Frozen recurrence script hashes changed")

    script_paths = {
        "aggregator_source_sha256": (
            repo_root
            / "benchmarks/brca100_consensus_recurrence_sweep/aggregate_recurrence.cpp"
        ),
        "runner_sha256": (
            repo_root
            / "benchmarks/brca100_consensus_recurrence_sweep/run_recurrence_sweep.py"
        ),
        "netbid_r_sha256": (
            repo_root
            / "benchmarks/brca100_pr67_threshold_sweep/run_netbid_qc.R"
        ),
    }
    for field, script_path in script_paths.items():
        require_hash(script_path, EXPECTED_FROZEN_SCRIPT_HASHES[field], field)
    return path, design


def validate_frozen_package(
    repo_root: Path, recurrence_root: Path
) -> Path:
    path = (
        repo_root
        / "benchmarks/brca100_consensus_recurrence_sweep"
        / "results_2026-08-20/package_manifest.json"
    )
    require_hash(path, EXPECTED_FROZEN_PACKAGE_SHA256, "frozen result package manifest")
    manifest = load_json(path)
    if (
        manifest.get("fingerprint") != EXPECTED_FROZEN_PACKAGE_FINGERPRINT
        or manifest.get("design_fingerprint") != EXPECTED_DESIGN_FINGERPRINT
        or Path(str(manifest.get("live_work_root"))).resolve() != recurrence_root
    ):
        raise ValueError(f"Frozen result package provenance mismatch: {path}")
    return path


def read_environment_file(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return parse_environment_table(handle, str(path))


def read_summary_metrics(path: Path) -> dict[str, float]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or any(set(row) != {"metric", "value"} for row in rows):
        raise ValueError(f"Malformed network summary: {path}")
    result: dict[str, float] = {}
    for row in rows:
        metric = row["metric"]
        if not metric or metric in result:
            raise ValueError(f"Duplicate/empty network metric: {path}")
        value = row["value"]
        result[metric] = math.nan if value == "NA" else float(value)
    if set(result) != REQUIRED_METRICS:
        raise ValueError(f"Unexpected network summary metrics: {path}")
    return result


def validate_summary_semantics(
    root: Path,
    driver_ids: tuple[str, ...],
    expected_edges: int,
    expected_environment: dict[str, str],
) -> dict[str, float]:
    metrics = read_summary_metrics(root / "network_summary.tsv")
    if int(metrics["candidate_drivers"]) != len(driver_ids):
        raise ValueError(f"Candidate-driver count mismatch: {root}")
    if int(metrics["edges"]) != expected_edges:
        raise ValueError(f"Network edge count mismatch: {root}")

    with (root / "driver_target_sizes.tsv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if any(set(row) != {"driver", "target_count"} for row in rows):
        raise ValueError(f"Malformed driver-target table: {root}")
    if tuple(row["driver"] for row in rows) != driver_ids:
        raise ValueError(f"Driver order/content mismatch: {root}")
    counts = [int(row["target_count"]) for row in rows]
    if any(value < 0 for value in counts) or sum(counts) != expected_edges:
        raise ValueError(f"Driver target counts disagree with edge count: {root}")
    if sum(value > 0 for value in counts) != int(metrics["active_drivers"]):
        raise ValueError(f"Active-driver count mismatch: {root}")
    if read_environment_file(root / "netbid_environment.tsv") != expected_environment:
        raise ValueError(f"NetBID2 environment artifact mismatch: {root}")
    return metrics


def validate_frozen_arm(
    context: Context, driver: str, threshold: int
) -> dict[str, Any]:
    driver_spec = DRIVER_SPECS[driver]
    run_spec = RUN_SPECS[(driver, threshold)]
    arm_root = context.recurrence_root / "results" / driver / f"k{threshold:03d}"
    network_path = arm_root / "consensus_network_ncol_.txt"
    network_manifest_path = arm_root / "network_manifest.json"
    summary_manifest_path = arm_root / "netbid2_manifest.json"
    summary_root = arm_root / "netbid2_qc"
    driver_path = context.source_root / "inputs" / driver_spec.filename

    require_hash(driver_path, driver_spec.sha256, f"{driver} driver list")
    if driver_path.stat().st_size == 0:
        raise ValueError(f"Empty driver list: {driver_path}")
    require_hash(network_path, run_spec.network_sha256, f"{driver} K={threshold} network")
    if network_path.stat().st_size != run_spec.network_bytes:
        raise ValueError(f"Network byte-count mismatch: {network_path}")
    require_hash(
        network_manifest_path,
        run_spec.network_manifest_sha256,
        f"{driver} K={threshold} network manifest",
    )
    network_manifest = load_json(network_manifest_path)
    expected_network_fields = {
        "schema": NETWORK_SCHEMA,
        "driver": driver,
        "minimum_support": threshold,
        "edges": run_spec.edges,
        "ncol_sha256": run_spec.network_sha256,
        "source_point": driver_spec.source_point,
        "per_subsample_p": driver_spec.per_subsample_p,
    }
    if any(network_manifest.get(key) != value for key, value in expected_network_fields.items()):
        raise ValueError(f"Frozen network manifest mismatch: {network_manifest_path}")

    require_hash(
        summary_manifest_path,
        run_spec.summary_manifest_sha256,
        f"{driver} K={threshold} summary manifest",
    )
    summary_manifest = load_json(summary_manifest_path)
    summary_input = summary_manifest.get("input")
    if not isinstance(summary_input, dict):
        raise ValueError(f"Malformed summary manifest input: {summary_manifest_path}")
    expected_summary_input = {
        "network_manifest_sha256": run_spec.network_manifest_sha256,
        "network_sha256": run_spec.network_sha256,
        "driver_sha256": driver_spec.sha256,
        "r_script_sha256": EXPECTED_FROZEN_SCRIPT_HASHES["netbid_r_sha256"],
        "wrapper_sha256": EXPECTED_WRAPPER_SHA256,
        "environment": context.environment,
        "generate_html": False,
    }
    if (
        summary_manifest.get("schema") != SUMMARY_SCHEMA
        or summary_manifest.get("driver") != driver
        or summary_manifest.get("minimum_support") != threshold
        or summary_manifest.get("outputs") != run_spec.summary_outputs
        or summary_input != expected_summary_input
    ):
        raise ValueError(f"Frozen summary manifest mismatch: {summary_manifest_path}")

    if not summary_root.is_dir() or summary_root.is_symlink():
        raise FileNotFoundError(f"Missing summary-mode NetBID2 output: {summary_root}")
    summary_children = list(summary_root.iterdir())
    if any(path.is_symlink() or not path.is_file() for path in summary_children):
        raise ValueError(f"Non-regular summary-mode output: {summary_root}")
    actual_summary_files = {path.name for path in summary_children}
    if actual_summary_files != SUMMARY_FILENAMES:
        raise ValueError(f"Unexpected summary-mode output inventory: {summary_root}")
    for filename, expected_hash in run_spec.summary_outputs.items():
        require_hash(summary_root / filename, expected_hash, f"summary output {filename}")
    metrics = validate_summary_semantics(
        summary_root,
        context.driver_ids[driver],
        run_spec.edges,
        context.environment,
    )
    return {
        "arm_root": arm_root,
        "network": network_path,
        "network_manifest": network_manifest_path,
        "summary_manifest": summary_manifest_path,
        "summary_root": summary_root,
        "driver_file": driver_path,
        "metrics": metrics,
    }


def build_context(
    repo_root: Path,
    source_root: Path,
    recurrence_root: Path,
    html_root: Path,
) -> Context:
    repo_root = repo_root.resolve(strict=True)
    source_root = source_root.resolve(strict=True)
    recurrence_root = recurrence_root.resolve(strict=True)
    html_root = html_root.resolve(strict=False)
    validate_path_separation(repo_root, source_root, recurrence_root, html_root)
    design_path, _ = validate_design(repo_root, source_root, recurrence_root)
    frozen_package = validate_frozen_package(repo_root, recurrence_root)
    wrapper = repo_root / "benchmarks/brca100_netbid_qc/netbid2-r"
    r_script = repo_root / "benchmarks/brca100_pr67_threshold_sweep/run_netbid_qc.R"
    generator = (
        repo_root
        / "benchmarks/brca100_consensus_recurrence_sweep/generate_netbid_html.py"
    )
    frozen_runner = (
        repo_root
        / "benchmarks/brca100_consensus_recurrence_sweep/run_recurrence_sweep.py"
    )
    aggregator = (
        repo_root
        / "benchmarks/brca100_consensus_recurrence_sweep/aggregate_recurrence.cpp"
    )
    require_hash(wrapper, EXPECTED_WRAPPER_SHA256, "NetBID2 environment wrapper")
    require_hash(
        r_script,
        EXPECTED_FROZEN_SCRIPT_HASHES["netbid_r_sha256"],
        "NetBID2 QC R script",
    )
    require_file(generator, "HTML adjunct generator")
    environment = probe_environment(wrapper, r_script)
    require_hash(wrapper, EXPECTED_WRAPPER_SHA256, "NetBID2 environment wrapper")
    require_hash(
        r_script,
        EXPECTED_FROZEN_SCRIPT_HASHES["netbid_r_sha256"],
        "NetBID2 QC R script after environment probe",
    )
    driver_ids = {
        driver: read_driver_ids(
            source_root / "inputs" / spec.filename, spec.candidate_drivers
        )
        for driver, spec in DRIVER_SPECS.items()
    }
    context = Context(
        repo_root=repo_root,
        source_root=source_root,
        recurrence_root=recurrence_root,
        html_root=html_root,
        design_path=design_path,
        design_sha256=EXPECTED_DESIGN_SHA256,
        design_fingerprint=EXPECTED_DESIGN_FINGERPRINT,
        frozen_package_manifest=frozen_package,
        wrapper=wrapper,
        r_script=r_script,
        generator=generator,
        frozen_runner=frozen_runner,
        aggregator_source=aggregator,
        environment=environment,
        generator_sha256=sha256_file(generator),
        driver_ids=driver_ids,
    )
    for driver, threshold in RUN_ORDER:
        validate_frozen_arm(context, driver, threshold)
    return context


def assert_runtime_sources(context: Context) -> None:
    expected = {
        context.wrapper: (EXPECTED_WRAPPER_SHA256, "NetBID2 environment wrapper"),
        context.r_script: (
            EXPECTED_FROZEN_SCRIPT_HASHES["netbid_r_sha256"],
            "NetBID2 QC R script",
        ),
        context.generator: (context.generator_sha256, "HTML adjunct generator"),
        context.frozen_runner: (
            EXPECTED_FROZEN_SCRIPT_HASHES["runner_sha256"],
            "frozen recurrence runner",
        ),
        context.aggregator_source: (
            EXPECTED_FROZEN_SCRIPT_HASHES["aggregator_source_sha256"],
            "frozen recurrence aggregator",
        ),
        context.design_path: (context.design_sha256, "frozen recurrence design"),
        context.frozen_package_manifest: (
            EXPECTED_FROZEN_PACKAGE_SHA256,
            "frozen result package manifest",
        ),
    }
    for path, (digest, description) in expected.items():
        require_hash(path, digest, description)


def runtime_r_script(context: Context) -> Path:
    return context.html_root / "runtime/run_netbid_qc.R"


def prepare_runtime(context: Context) -> Path:
    active_overlay_lock(context)
    assert_runtime_sources(context)
    runtime_root = secure_subdirectory(context.html_root, "runtime")
    snapshot = runtime_root / "run_netbid_qc.R"
    assert_managed_path(context.html_root, snapshot)
    expected_hash = EXPECTED_FROZEN_SCRIPT_HASHES["netbid_r_sha256"]
    if snapshot.exists() or snapshot.is_symlink():
        require_hash(snapshot, expected_hash, "immutable NetBID2 R-script snapshot")
    else:
        atomic_bytes(snapshot, context.r_script.read_bytes())
        require_hash(snapshot, expected_hash, "immutable NetBID2 R-script snapshot")
    environment = probe_environment(context.wrapper, snapshot)
    assert_runtime_sources(context)
    require_hash(snapshot, expected_hash, "NetBID2 R-script snapshot after probe")
    if environment != context.environment:
        raise ValueError("Runtime NetBID2 environment changed after preflight")
    return snapshot


def inventory(root: Path) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Symlink not allowed in artifact inventory: {path}")
        if path.is_dir():
            continue
        if not path.is_file():
            raise ValueError(f"Non-regular artifact: {path}")
        records.append(
            {
                "path": path.relative_to(root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    return records


def expected_html_names(prefix: str) -> set[str]:
    return SUMMARY_FILENAMES | {f"{prefix}netQC.Rmd", f"{prefix}netQC.html"}


def validate_html_output(
    root: Path,
    summary_root: Path,
    prefix: str,
    driver_ids: tuple[str, ...],
    expected_edges: int,
    expected_environment: dict[str, str],
) -> tuple[list[dict[str, object]], dict[str, float]]:
    if not root.is_dir() or root.is_symlink():
        raise FileNotFoundError(f"Missing HTML-mode output directory: {root}")
    actual_inventory = inventory(root)
    if {item["path"] for item in actual_inventory} != expected_html_names(prefix):
        raise ValueError(f"Unexpected HTML-mode output inventory: {root}")
    if any(int(item["bytes"]) <= 0 for item in actual_inventory):
        raise ValueError(f"Empty HTML-mode output artifact: {root}")
    html_path = root / f"{prefix}netQC.html"
    with html_path.open("rb") as handle:
        if b"<html" not in handle.read(4096).lower():
            raise ValueError(f"NetBID2 report is not recognizable HTML: {html_path}")
    for filename in SUMMARY_FILENAMES:
        generated = root / filename
        frozen = summary_root / filename
        hashes_differ = sha256_file(generated) != sha256_file(frozen)
        bytes_differ = generated.read_bytes() != frozen.read_bytes()
        if hashes_differ or bytes_differ:
            raise ValueError(
                f"HTML/summary NetBID2 TSV bytes differ: {generated} != {frozen}"
            )
    metrics = validate_summary_semantics(
        root, driver_ids, expected_edges, expected_environment
    )
    return actual_inventory, metrics


def log_record(path: Path) -> dict[str, object]:
    require_file(path, "NetBID2 log")
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def arm_paths(context: Context, driver: str, threshold: int) -> dict[str, Path]:
    root = context.html_root / "runs" / driver / f"k{threshold:03d}"
    return {
        "root": root,
        "final": root / "netbid2_qc_html",
        "partial": root / "netbid2_qc_html.partial",
        "manifest": root / "netbid2_qc_html_manifest.json",
        "pending": root / "netbid2_qc_html_manifest.pending.json",
        "owner": root / "netbid2_qc_html.partial.owner.json",
        "stdout": root / "logs/netbid2_qc_html.stdout.log",
        "stderr": root / "logs/netbid2_qc_html.stderr.log",
    }


def partial_owner_record(
    context: Context, driver: str, threshold: int, input_fingerprint: str
) -> dict[str, Any]:
    lock = active_overlay_lock(context)
    return canonical_record(
        {
            "schema": PARTIAL_OWNER_SCHEMA,
            "pid": os.getpid(),
            "lock_token": lock["token"],
            "driver": driver,
            "minimum_support": threshold,
            "input_fingerprint": input_fingerprint,
        }
    )


def validate_partial_owner(
    path: Path,
    *,
    driver: str,
    threshold: int,
    input_fingerprint: str,
) -> dict[str, Any]:
    require_file(path, "HTML partial owner record")
    record = load_json(path)
    expected_fields = {
        "schema",
        "pid",
        "lock_token",
        "driver",
        "minimum_support",
        "input_fingerprint",
        "fingerprint",
    }
    if set(record) != expected_fields:
        raise RuntimeError(f"Malformed HTML partial owner record: {path}")
    validate_record_fingerprint(record, path)
    if (
        record.get("schema") != PARTIAL_OWNER_SCHEMA
        or record.get("driver") != driver
        or record.get("minimum_support") != threshold
        or record.get("input_fingerprint") != input_fingerprint
        or not isinstance(record.get("pid"), int)
        or not isinstance(record.get("lock_token"), str)
    ):
        raise RuntimeError(f"HTML partial owner provenance mismatch: {path}")
    return record


def arm_input(
    context: Context,
    driver: str,
    threshold: int,
    frozen: dict[str, Any],
) -> dict[str, Any]:
    driver_spec = DRIVER_SPECS[driver]
    run_spec = RUN_SPECS[(driver, threshold)]
    return {
        "design_fingerprint": context.design_fingerprint,
        "design_sha256": context.design_sha256,
        "network_manifest_sha256": run_spec.network_manifest_sha256,
        "network_sha256": run_spec.network_sha256,
        "network_bytes": run_spec.network_bytes,
        "driver_sha256": driver_spec.sha256,
        "summary_manifest_sha256": run_spec.summary_manifest_sha256,
        "summary_outputs": run_spec.summary_outputs,
        "r_script_sha256": EXPECTED_FROZEN_SCRIPT_HASHES["netbid_r_sha256"],
        "wrapper_sha256": EXPECTED_WRAPPER_SHA256,
        "generator_sha256": context.generator_sha256,
        "environment": context.environment,
        "expected_edges": run_spec.edges,
        "expected_candidate_drivers": driver_spec.candidate_drivers,
        "prefix": driver_spec.prefix,
        "generate_html": True,
        "network": str(frozen["network"]),
        "driver_file": str(frozen["driver_file"]),
    }


def validate_html_run_record(
    context: Context,
    driver: str,
    threshold: int,
    record: dict[str, Any],
    source: Path,
    output_root: Path,
) -> dict[str, Any]:
    expected_fields = {
        "schema",
        "driver",
        "minimum_support",
        "input",
        "input_fingerprint",
        "command",
        "output",
        "output_inventory",
        "stdout",
        "stderr",
        "finished_at_utc",
        "fingerprint",
    }
    if set(record) != expected_fields:
        raise ValueError(f"Malformed HTML run manifest fields: {source}")
    validate_record_fingerprint(record, source)
    frozen = validate_frozen_arm(context, driver, threshold)
    expected_input = arm_input(context, driver, threshold, frozen)
    paths = arm_paths(context, driver, threshold)
    expected_command = [
        str(context.wrapper),
        "Rscript",
        str(runtime_r_script(context)),
        str(frozen["network"]),
        str(frozen["driver_file"]),
        str(paths["partial"]),
        DRIVER_SPECS[driver].prefix,
        "true",
    ]
    if (
        record.get("schema") != HTML_RUN_SCHEMA
        or record.get("driver") != driver
        or record.get("minimum_support") != threshold
        or record.get("input") != expected_input
        or record.get("input_fingerprint") != canonical_fingerprint(expected_input)
        or record.get("command") != expected_command
        or record.get("output") != str(paths["final"])
    ):
        raise ValueError(f"HTML run provenance mismatch: {source}")
    try:
        dt.datetime.fromisoformat(str(record["finished_at_utc"]))
    except ValueError as exc:
        raise ValueError(f"Invalid HTML run completion timestamp: {source}") from exc
    actual_inventory, metrics = validate_html_output(
        output_root,
        frozen["summary_root"],
        DRIVER_SPECS[driver].prefix,
        context.driver_ids[driver],
        RUN_SPECS[(driver, threshold)].edges,
        context.environment,
    )
    if record.get("output_inventory") != actual_inventory:
        raise ValueError(f"HTML output inventory changed: {source}")
    for key in ("stdout", "stderr"):
        expected_log = log_record(paths[key])
        if record.get(key) != expected_log:
            raise ValueError(f"HTML run log changed: {source}")
    return {"record": record, "metrics": metrics, "frozen": frozen}


def reclaim_partial_owner(
    context: Context,
    paths: dict[str, Path],
    driver: str,
    threshold: int,
    input_fingerprint: str,
) -> dict[str, Any]:
    owner = validate_partial_owner(
        paths["owner"],
        driver=driver,
        threshold=threshold,
        input_fingerprint=input_fingerprint,
    )
    active = active_overlay_lock(context)
    owned_by_this_run = (
        owner["pid"] == os.getpid() and owner["lock_token"] == active["token"]
    )
    if pid_is_alive(int(owner["pid"])) and not owned_by_this_run:
        raise RuntimeError(
            f"HTML partial is owned by live PID {owner['pid']}: {paths['partial']}"
        )
    return owner


def remove_managed_partial(
    context: Context, partial: Path, expected_parent: Path
) -> None:
    if partial.parent != expected_parent or partial.name != "netbid2_qc_html.partial":
        raise RuntimeError(f"Refusing to remove unexpected partial path: {partial}")
    assert_managed_path(context.html_root, partial)
    if partial.is_symlink() or not partial.is_dir():
        raise RuntimeError(f"Unverifiable HTML partial state: {partial}")
    resolved = partial.resolve(strict=True)
    if not is_within(resolved, context.html_root.resolve(strict=True)):
        raise RuntimeError(f"HTML partial resolved outside its overlay: {partial}")
    shutil.rmtree(partial)


def run_html_arm(context: Context, driver: str, threshold: int) -> dict[str, Any]:
    active_overlay_lock(context)
    paths = arm_paths(context, driver, threshold)
    frozen = validate_frozen_arm(context, driver, threshold)
    expected_input = arm_input(context, driver, threshold, frozen)
    input_fingerprint = canonical_fingerprint(expected_input)
    secure_subdirectory(
        context.html_root, "runs", driver, f"k{threshold:03d}"
    )
    secure_subdirectory(
        context.html_root, "runs", driver, f"k{threshold:03d}", "logs"
    )
    for path in paths.values():
        assert_managed_path(context.html_root, path)

    if paths["manifest"].is_file() and paths["final"].is_dir():
        if paths["partial"].exists() or paths["partial"].is_symlink():
            raise RuntimeError(f"Completed and partial HTML outputs coexist: {paths['root']}")
        record = load_json(paths["manifest"])
        if paths["pending"].exists() or paths["pending"].is_symlink():
            require_file(paths["pending"], "pending HTML run manifest")
            pending = load_json(paths["pending"])
            if pending != record:
                raise RuntimeError(f"Pending/completed HTML manifests differ: {paths['root']}")
        validate_html_run_record(
            context, driver, threshold, record, paths["manifest"], paths["final"]
        )
        if paths["pending"].exists():
            managed_unlink(context, paths["pending"])
        if paths["owner"].exists() or paths["owner"].is_symlink():
            reclaim_partial_owner(
                context, paths, driver, threshold, input_fingerprint
            )
            managed_unlink(context, paths["owner"])
        print(f"[NETBID2 HTML] {driver} K={threshold} resume", flush=True)
        return record

    if paths["pending"].is_file():
        reclaim_partial_owner(context, paths, driver, threshold, input_fingerprint)
        if paths["final"].is_dir() == paths["partial"].is_dir():
            raise RuntimeError(f"Ambiguous pending HTML state: {paths['root']}")
        recovery_root = paths["final"] if paths["final"].is_dir() else paths["partial"]
        record = load_json(paths["pending"])
        validate_html_run_record(
            context, driver, threshold, record, paths["pending"], recovery_root
        )
        if recovery_root == paths["partial"]:
            managed_rename(context, paths["partial"], paths["final"])
        managed_atomic_json(context, paths["manifest"], record)
        managed_unlink(context, paths["pending"])
        managed_unlink(context, paths["owner"])
        print(f"[NETBID2 HTML] {driver} K={threshold} recovered", flush=True)
        return record

    if any(
        path.exists() or path.is_symlink()
        for path in (paths["manifest"], paths["final"], paths["pending"])
    ):
        raise RuntimeError(f"Unverifiable completed HTML state: {paths['root']}")
    if paths["partial"].exists() or paths["partial"].is_symlink():
        reclaim_partial_owner(context, paths, driver, threshold, input_fingerprint)
        remove_managed_partial(context, paths["partial"], paths["root"])
        managed_unlink(context, paths["owner"])
    elif paths["owner"].exists() or paths["owner"].is_symlink():
        reclaim_partial_owner(context, paths, driver, threshold, input_fingerprint)
        managed_unlink(context, paths["owner"])

    owner = partial_owner_record(context, driver, threshold, input_fingerprint)
    managed_atomic_json(context, paths["owner"], owner)
    command = [
        str(context.wrapper),
        "Rscript",
        str(runtime_r_script(context)),
        str(frozen["network"]),
        str(frozen["driver_file"]),
        str(paths["partial"]),
        DRIVER_SPECS[driver].prefix,
        "true",
    ]
    print(f"[NETBID2 HTML] {driver} K={threshold}", flush=True)
    assert_runtime_sources(context)
    require_hash(
        runtime_r_script(context),
        EXPECTED_FROZEN_SCRIPT_HASHES["netbid_r_sha256"],
        "NetBID2 R-script runtime snapshot",
    )
    assert_managed_path(context.html_root, paths["stdout"])
    assert_managed_path(context.html_root, paths["stderr"])
    with paths["stdout"].open("w", encoding="utf-8", newline="\n") as stdout, paths[
        "stderr"
    ].open("w", encoding="utf-8", newline="\n") as stderr:
        completed = subprocess.run(command, stdout=stdout, stderr=stderr, text=True)
    assert_runtime_sources(context)
    require_hash(
        runtime_r_script(context),
        EXPECTED_FROZEN_SCRIPT_HASHES["netbid_r_sha256"],
        "NetBID2 R-script runtime snapshot after HTML run",
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"NetBID2 HTML failed for {driver} K={threshold}; see {paths['stderr']}"
        )

    output_inventory, _ = validate_html_output(
        paths["partial"],
        frozen["summary_root"],
        DRIVER_SPECS[driver].prefix,
        context.driver_ids[driver],
        RUN_SPECS[(driver, threshold)].edges,
        context.environment,
    )
    record = canonical_record(
        {
            "schema": HTML_RUN_SCHEMA,
            "driver": driver,
            "minimum_support": threshold,
            "input": expected_input,
            "input_fingerprint": input_fingerprint,
            "command": command,
            "output": str(paths["final"]),
            "output_inventory": output_inventory,
            "stdout": log_record(paths["stdout"]),
            "stderr": log_record(paths["stderr"]),
            "finished_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        }
    )
    managed_atomic_json(context, paths["pending"], record)
    managed_rename(context, paths["partial"], paths["final"])
    managed_atomic_json(context, paths["manifest"], record)
    managed_unlink(context, paths["pending"])
    managed_unlink(context, paths["owner"])
    validate_html_run_record(
        context, driver, threshold, record, paths["manifest"], paths["final"]
    )
    return record


def aggregate_payload(context: Context, records: list[dict[str, Any]]) -> dict[str, Any]:
    runs: list[dict[str, Any]] = []
    for (driver, threshold), record in zip(RUN_ORDER, records, strict=True):
        paths = arm_paths(context, driver, threshold)
        html_name = f"{DRIVER_SPECS[driver].prefix}netQC.html"
        runs.append(
            {
                "driver": driver,
                "minimum_support": threshold,
                "manifest": str(paths["manifest"]),
                "manifest_sha256": sha256_file(paths["manifest"]),
                "manifest_fingerprint": record["fingerprint"],
                "html_sha256": sha256_file(paths["final"] / html_name),
                "html_bytes": (paths["final"] / html_name).stat().st_size,
            }
        )
    return canonical_record(
        {
            "schema": HTML_AGGREGATE_SCHEMA,
            "design_sha256": context.design_sha256,
            "design_fingerprint": context.design_fingerprint,
            "generator_sha256": context.generator_sha256,
            "environment": context.environment,
            "selection": [
                {"driver": driver, "minimum_support": threshold}
                for driver, threshold in RUN_ORDER
            ],
            "runs": runs,
        }
    )


def _generate_html_locked(context: Context) -> dict[str, Any]:
    active_overlay_lock(context)
    records = [run_html_arm(context, driver, threshold) for driver, threshold in RUN_ORDER]
    aggregate = aggregate_payload(context, records)
    aggregate_path = context.html_root / "netbid2_html_manifest.json"
    assert_managed_path(context.html_root, aggregate_path)
    if aggregate_path.exists() or aggregate_path.is_symlink():
        require_file(aggregate_path, "HTML aggregate manifest")
        existing = load_json(aggregate_path)
        validate_record_fingerprint(existing, aggregate_path)
        if existing != aggregate:
            raise RuntimeError(f"Existing HTML aggregate is stale: {aggregate_path}")
    else:
        managed_atomic_json(context, aggregate_path, aggregate)
    return aggregate


def generate_html(context: Context) -> dict[str, Any]:
    with overlay_lock(context):
        prepare_runtime(context)
        return _generate_html_locked(context)


def _validate_html_overlay_locked(context: Context) -> dict[str, Any]:
    active_overlay_lock(context)
    records: list[dict[str, Any]] = []
    for driver, threshold in RUN_ORDER:
        paths = arm_paths(context, driver, threshold)
        for path in paths.values():
            assert_managed_path(context.html_root, path)
        if not paths["manifest"].is_file() or not paths["final"].is_dir():
            raise FileNotFoundError(f"Incomplete HTML arm: {paths['root']}")
        if any(
            path.exists() or path.is_symlink()
            for path in (paths["pending"], paths["partial"], paths["owner"])
        ):
            raise RuntimeError(f"Unresolved HTML state: {paths['root']}")
        record = load_json(paths["manifest"])
        validate_html_run_record(
            context, driver, threshold, record, paths["manifest"], paths["final"]
        )
        records.append(record)
    expected = aggregate_payload(context, records)
    aggregate_path = context.html_root / "netbid2_html_manifest.json"
    assert_managed_path(context.html_root, aggregate_path)
    require_file(aggregate_path, "HTML aggregate manifest")
    actual = load_json(aggregate_path)
    validate_record_fingerprint(actual, aggregate_path)
    if actual != expected:
        raise ValueError(f"HTML aggregate manifest mismatch: {aggregate_path}")
    return actual


def validate_html_overlay(context: Context) -> dict[str, Any]:
    with overlay_lock(context):
        prepare_runtime(context)
        return _validate_html_overlay_locked(context)


def copy_verified(
    source: Path, destination: Path, expected_sha256: str | None = None
) -> None:
    require_file(source, "package source")
    source_hash = sha256_file(source)
    if expected_sha256 is not None and source_hash != expected_sha256:
        raise ValueError(f"Package source changed before copy: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    if source_hash != sha256_file(destination):
        raise RuntimeError(f"Package copy hash mismatch: {source} -> {destination}")
    if sha256_file(source) != source_hash:
        raise RuntimeError(f"Package source changed during copy: {source}")


def package_report_name(driver: str, threshold: int) -> str:
    return f"k{threshold:03d}_{driver}_netbid2_qc.html"


def render_readme(reports: list[dict[str, Any]], context: Context) -> str:
    lines = [
        "# Representative NetBID2 recurrence-QC reports",
        "",
        "These four single-file reports compare the fixed BRCA100 TF and SIG",
        "networks at minimum recurrence counts `K=6` and `K=8`. Inputs are the",
        "post-DPI networks from the frozen recurrence sweep; HTML generation",
        "does not rerun SJARACNe inference or change the recurrence design.",
        "",
        "| K | Driver | Edges | Adjusted scale-free R2 | Report | Bytes | SHA-256 |",
        "|---:|---|---:|---:|---|---:|---|",
    ]
    for report in reports:
        lines.append(
            "| {minimum_support} | {driver} | {edges:,} | {scale_free_adjusted_r2:.6g} "
            "| [`{filename}`]({filename}) | {bytes:,} | `{sha256}` |".format(**report)
        )
    lines.extend(
        [
            "",
            "## Provenance",
            "",
            f"- Frozen design fingerprint: `{context.design_fingerprint}`",
            f"- Frozen design file SHA-256: `{context.design_sha256}`",
            "- R 4.4.3 (2025-02-28)",
            "- NetBID2 2.2.0, remote commit "
            "`5defa454d600b94f5dd6d1f9f4428f99759a6821`",
            "- igraph 2.3.3",
            "",
            "The `provenance/` directory contains the aggregate and per-arm run",
            "manifests, logs, R Markdown sources, exact TSV outputs, frozen input",
            "manifests, and scripts used for generation. HTML-mode TSV outputs",
            "were required to be byte-identical to the prior summary-mode outputs.",
            "",
            "## Viewing",
            "",
            "GitHub generally displays HTML as source. Download the raw file, verify",
            "the package with `sha256sum -c SHA256SUMS`, and open it locally. The",
            "reports are topology QC, not biological validation or edge-level FDR",
            "estimates.",
            "",
        ]
    )
    return "\n".join(lines)


def package_source_files(
    context: Context, partial: Path, reports: list[dict[str, Any]]
) -> None:
    active_overlay_lock(context)
    assert_runtime_sources(context)
    copy_verified(
        context.html_root / "netbid2_html_manifest.json",
        partial / "provenance/netbid2_html_manifest.json",
    )
    copy_verified(context.design_path, partial / "provenance/frozen/design.json")
    copy_verified(
        context.frozen_package_manifest,
        partial / "provenance/frozen/result_package_manifest.json",
    )
    script_sources = {
        "generate_netbid_html.py": (context.generator, context.generator_sha256),
        "run_netbid_qc.R": (
            runtime_r_script(context),
            EXPECTED_FROZEN_SCRIPT_HASHES["netbid_r_sha256"],
        ),
        "netbid2-r": (context.wrapper, EXPECTED_WRAPPER_SHA256),
        "run_recurrence_sweep.py": (
            context.frozen_runner,
            EXPECTED_FROZEN_SCRIPT_HASHES["runner_sha256"],
        ),
        "aggregate_recurrence.cpp": (
            context.aggregator_source,
            EXPECTED_FROZEN_SCRIPT_HASHES["aggregator_source_sha256"],
        ),
    }
    for name, (source, expected_hash) in script_sources.items():
        copy_verified(
            source, partial / "provenance/scripts" / name, expected_hash
        )

    for driver, threshold in RUN_ORDER:
        paths = arm_paths(context, driver, threshold)
        frozen = validate_frozen_arm(context, driver, threshold)
        arm_destination = partial / "provenance/arms" / driver / f"k{threshold:03d}"
        copy_verified(paths["manifest"], arm_destination / "netbid2_qc_html_manifest.json")
        copy_verified(paths["stdout"], arm_destination / "netbid2_qc_html.stdout.log")
        copy_verified(paths["stderr"], arm_destination / "netbid2_qc_html.stderr.log")
        copy_verified(frozen["network_manifest"], arm_destination / "network_manifest.json")
        copy_verified(frozen["summary_manifest"], arm_destination / "summary_manifest.json")
        prefix = DRIVER_SPECS[driver].prefix
        for filename in SUMMARY_FILENAMES | {f"{prefix}netQC.Rmd"}:
            copy_verified(paths["final"] / filename, arm_destination / filename)
        html_source = paths["final"] / f"{prefix}netQC.html"
        html_destination = partial / package_report_name(driver, threshold)
        copy_verified(html_source, html_destination)

    atomic_bytes(partial / ".gitattributes", b"* -text -whitespace\n")
    atomic_bytes(partial / "README.md", render_readme(reports, context).encode("utf-8"))
    assert_runtime_sources(context)


def build_report_records(context: Context) -> list[dict[str, Any]]:
    reports: list[dict[str, Any]] = []
    for driver, threshold in RUN_ORDER:
        paths = arm_paths(context, driver, threshold)
        prefix = DRIVER_SPECS[driver].prefix
        html_path = paths["final"] / f"{prefix}netQC.html"
        frozen = validate_frozen_arm(context, driver, threshold)
        reports.append(
            {
                "driver": driver,
                "minimum_support": threshold,
                "filename": package_report_name(driver, threshold),
                "bytes": html_path.stat().st_size,
                "sha256": sha256_file(html_path),
                "edges": RUN_SPECS[(driver, threshold)].edges,
                "scale_free_adjusted_r2": frozen["metrics"]["scale_free_adjusted_r2"],
                "network_sha256": RUN_SPECS[(driver, threshold)].network_sha256,
                "run_manifest_sha256": sha256_file(paths["manifest"]),
            }
        )
    return reports


def parse_sha256sums(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if len(line) < 67 or line[64:66] != "  ":
            raise ValueError(f"Malformed SHA256SUMS line: {line!r}")
        digest, relative = line[:64], line[66:]
        if relative in result or any(character not in "0123456789abcdef" for character in digest):
            raise ValueError(f"Malformed/duplicate SHA256SUMS entry: {line!r}")
        result[relative] = digest
    return result


def validate_package(
    package_root: Path, context: Context, aggregate: dict[str, Any]
) -> dict[str, Any]:
    active_overlay_lock(context)
    aggregate_path = context.html_root / "netbid2_html_manifest.json"
    assert_managed_path(context.html_root, aggregate_path)
    require_file(aggregate_path, "HTML aggregate manifest")
    canonical_aggregate = load_json(aggregate_path)
    validate_record_fingerprint(canonical_aggregate, aggregate_path)
    if aggregate != canonical_aggregate:
        raise ValueError("Caller-supplied aggregate is not the canonical overlay aggregate")
    if not package_root.is_dir() or package_root.is_symlink():
        raise FileNotFoundError(f"Missing HTML package: {package_root}")
    manifest_path = package_root / "package_manifest.json"
    checksum_path = package_root / "SHA256SUMS"
    manifest = load_json(manifest_path)
    expected_fields = {
        "schema",
        "source_aggregate_sha256",
        "source_aggregate_fingerprint",
        "design_sha256",
        "design_fingerprint",
        "generator_sha256",
        "reports",
        "files",
        "fingerprint",
    }
    if set(manifest) != expected_fields:
        raise ValueError(f"Malformed HTML package manifest: {manifest_path}")
    validate_record_fingerprint(manifest, manifest_path)
    if (
        manifest.get("schema") != HTML_PACKAGE_SCHEMA
        or manifest.get("source_aggregate_sha256") != sha256_file(aggregate_path)
        or manifest.get("source_aggregate_fingerprint") != aggregate["fingerprint"]
        or manifest.get("design_sha256") != context.design_sha256
        or manifest.get("design_fingerprint") != context.design_fingerprint
        or manifest.get("generator_sha256") != context.generator_sha256
        or manifest.get("reports") != build_report_records(context)
    ):
        raise ValueError(f"HTML package provenance mismatch: {manifest_path}")
    actual_without_control = [
        item
        for item in inventory(package_root)
        if item["path"] not in {"package_manifest.json", "SHA256SUMS"}
    ]
    if manifest.get("files") != actual_without_control:
        raise ValueError(f"HTML package content inventory mismatch: {package_root}")
    sums = parse_sha256sums(checksum_path)
    expected_paths = {
        item["path"]
        for item in inventory(package_root)
        if item["path"] != "SHA256SUMS"
    }
    if set(sums) != expected_paths:
        raise ValueError(f"SHA256SUMS inventory mismatch: {checksum_path}")
    for relative, expected_hash in sums.items():
        if sha256_file(package_root / relative) != expected_hash:
            raise ValueError(f"SHA256SUMS mismatch: {relative}")
    if (package_root / ".gitattributes").read_bytes() != b"* -text -whitespace\n":
        raise ValueError(f"Unexpected .gitattributes: {package_root}")
    return manifest


def _package_html_locked(
    context: Context, package_root: Path, aggregate: dict[str, Any]
) -> dict[str, Any]:
    active_overlay_lock(context)
    package_root = package_root.resolve(strict=False)
    benchmark_root = (
        context.repo_root / "benchmarks/brca100_consensus_recurrence_sweep"
    ).resolve(strict=True)
    if not is_within(package_root, benchmark_root) or package_root == benchmark_root:
        raise ValueError(f"Package root must be a child of {benchmark_root}")
    frozen_results = (benchmark_root / "results_2026-08-20").resolve(strict=True)
    if package_root == frozen_results or is_within(package_root, frozen_results):
        raise ValueError("Refusing to write inside the frozen recurrence result package")
    if package_root.exists():
        return validate_package(package_root, context, aggregate)

    partial = package_root.with_name(package_root.name + ".partial")
    package_root.parent.mkdir(parents=True, exist_ok=True)
    if partial.exists():
        try:
            validate_package(partial, context, aggregate)
        except Exception as exc:
            raise RuntimeError(
                f"Incomplete or unverifiable package partial: {partial}"
            ) from exc
        if package_root.exists():
            raise FileExistsError(f"Package appeared during recovery: {package_root}")
        os.rename(partial, package_root)
        return validate_package(package_root, context, aggregate)
    partial.mkdir(parents=False, exist_ok=False)
    reports = build_report_records(context)
    package_source_files(context, partial, reports)
    files = inventory(partial)
    package_manifest = canonical_record(
        {
            "schema": HTML_PACKAGE_SCHEMA,
            "source_aggregate_sha256": sha256_file(
                context.html_root / "netbid2_html_manifest.json"
            ),
            "source_aggregate_fingerprint": aggregate["fingerprint"],
            "design_sha256": context.design_sha256,
            "design_fingerprint": context.design_fingerprint,
            "generator_sha256": context.generator_sha256,
            "reports": reports,
            "files": files,
        }
    )
    atomic_json(partial / "package_manifest.json", package_manifest)
    checksum_records = [
        item for item in inventory(partial) if item["path"] != "SHA256SUMS"
    ]
    checksum_payload = "".join(
        f"{item['sha256']}  {item['path']}\n" for item in checksum_records
    ).encode("utf-8")
    atomic_bytes(partial / "SHA256SUMS", checksum_payload)
    validate_package(partial, context, aggregate)
    if package_root.exists():
        raise FileExistsError(f"Package appeared during build: {package_root}")
    os.rename(partial, package_root)
    return validate_package(package_root, context, aggregate)


def package_html(
    context: Context, package_root: Path, aggregate: dict[str, Any]
) -> dict[str, Any]:
    with overlay_lock(context):
        prepare_runtime(context)
        canonical_aggregate = _validate_html_overlay_locked(context)
        if aggregate != canonical_aggregate:
            raise ValueError(
                "Caller-supplied aggregate is not the canonical overlay aggregate"
            )
        return _package_html_locked(context, package_root, canonical_aggregate)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=("all", "generate", "package"), default="all")
    parser.add_argument(
        "--benchmark-repo",
        type=Path,
        default=Path("/mnt/d/GitHub/SJARACNe-brca100-netbid-qc"),
    )
    parser.add_argument(
        "--source-work-root",
        type=Path,
        default=Path(
            "/home/adam/sjaracne-benchmarks/brca100-pr67-threshold-sweep-20260818"
        ),
    )
    parser.add_argument(
        "--recurrence-work-root",
        type=Path,
        default=Path(
            "/home/adam/sjaracne-benchmarks/"
            "brca100-consensus-recurrence-sweep-20260820"
        ),
    )
    parser.add_argument(
        "--html-work-root",
        type=Path,
        default=Path(
            "/home/adam/sjaracne-benchmarks/"
            "brca100-consensus-recurrence-netbid2-html-20260820"
        ),
    )
    parser.add_argument("--package-root", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.phase in {"all", "package"} and args.package_root is None:
        raise ValueError("--package-root is required for the package/all phase")
    context = build_context(
        args.benchmark_repo,
        args.source_work_root,
        args.recurrence_work_root,
        args.html_work_root,
    )
    if args.phase in {"all", "generate"}:
        aggregate = generate_html(context)
    else:
        aggregate = validate_html_overlay(context)
    if args.phase in {"all", "package"}:
        package_html(context, args.package_root, aggregate)
        print(f"[NETBID2 HTML] package complete: {args.package_root}", flush=True)
    else:
        print(f"[NETBID2 HTML] overlay complete: {context.html_root}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
