#!/usr/bin/env python3
"""Dependency-free helpers for the BRCA100 SIG K_DPI witness screen."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import re
from typing import Iterable, Iterator


BASELINE_DEFAULT = Path(
    "/home/adam/sjaracne-benchmarks/brca100-hub-size-dpi-pilot-20260825"
)
WORK_ROOT_DEFAULT = Path(
    "/home/adam/sjaracne-benchmarks/brca100-kdpi-witness-screen-20260826"
)
REPO_DEFAULT = Path("/mnt/d/GitHub/SJARACNe-hub-dpi")
HUB_COUNTS = (1335, 5340, 10680)
SEEDS = tuple(range(1, 11))
K_VALUES = (1, 2, 3, 5, 10)
P_TOKEN = "5e-4"
P_VALUE = 5e-4
MI_CUTOFF = 0.1644671599536221
N_PAR = 40
EXPRESSION_FILENAME = "BRCA100.exp"
DRIVER_FILENAME = "BRCA100_SIG.txt"
SIDECAR_FIELDS = (
    "source_index",
    "pre_edges",
    "witnesses_ge_1",
    "witnesses_ge_2",
    "witnesses_ge_3",
    "witnesses_ge_5",
    "witnesses_ge_10",
)
SIDECAR_PROVENANCE_KEYS = (
    "schema",
    "graph_state",
    "count_unit",
    "count_semantics",
    "source_index_basis",
    "dpi_epsilon",
    "source_mode",
    "source_count",
    "annotated_gene_count",
    "input_file",
    "input_adjacency_file",
    "network_output_file",
    "subnetwork_file",
    "annotation_file",
    "k1_pruned_edges",
)

EXPECTED_BASELINE = {
    "schema": "sjaracne-brca100-hub-size-dpi-pilot-v1",
    "source_commit": "32fe12c168ef80291e487dbf4045f430b9c5d90a",
    "source_tree_fingerprint": (
        "3be6f0c75cf905f983d789410542adda388d07bc0231b4439371e0dd3f22169d"
    ),
    "binary_sha256": (
        "836c8edbfc75b517cefdde6426755f9e3b71021df09b06c6f3b9c262057cd47e"
    ),
    "config_sha256": (
        "dca25c4041800933e6d439fe8730d4e78039839f956932968401540bdf454dca"
    ),
    "null_model_sha256": (
        "e3a8522682a8ea239821aaa10b12db72d00e07bfdcad43599d8e76a06be80944"
    ),
    "expression_sha256": (
        "ad8a334f5f8cdf46a1000d3ee259b35258a18b3da2e314bb3a0cf7a421d98bc8"
    ),
    "panel_sha256": {
        1335: "1e082c351e5b57dcff968dbde4b52ebb7385b3cc1f6bb79cd3981ad6bba6b92c",
        5340: "530bba715d8a5bb3216d9b4cccccbf34a253bd4678b2fa88595899d1e34a6f93",
        10680: "16ca27df655f16684f880a4ad719c4e2ae3f8dc0d7e6b9eccdd24cd97c40797c",
    },
}

_DPI_STATS = re.compile(
    r"^\[DPI_STATS\] pre_edges=(\d+) pruned_edges=(\d+) "
    r"post_edges=(\d+) dpi_applied=([01])$"
)
_SAMPLING_PREFIX = "[SAMPLING] Selected original observation indices (0-based):"


def serialized_json(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode(
        "utf-8"
    )


def atomic_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("wb") as handle:
        handle.write(payload)
    os.replace(temporary, path)


def atomic_json(path: Path, value: object) -> None:
    atomic_bytes(path, serialized_json(value))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint(value: object) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return sha256_bytes(payload)


def utc_now() -> str:
    import datetime as dt

    return dt.datetime.now(dt.timezone.utc).isoformat()


def arm_key(count: int) -> str:
    if count not in HUB_COUNTS:
        raise ValueError(f"Unsupported SIG hub count: {count}")
    return f"sig_h{count:05d}"


def result_root(work_root: Path, count: int) -> Path:
    return work_root / "results" / arm_key(count)


def baseline_arm_root(baseline_root: Path, count: int) -> Path:
    return baseline_root / "results" / arm_key(count)


def panel_path(baseline_root: Path, count: int) -> Path:
    return (
        baseline_root
        / "panels"
        / "sig"
        / f"h{count:05d}"
        / DRIVER_FILENAME
    )


def read_unique_ids(path: Path, expected_count: int | None = None) -> list[str]:
    values = path.read_text(encoding="utf-8").splitlines()
    if any(not value for value in values) or len(values) != len(set(values)):
        raise ValueError(f"Blank or duplicate identifier in {path}")
    if expected_count is not None and len(values) != expected_count:
        raise ValueError(f"Expected {expected_count} identifiers in {path}, got {len(values)}")
    return values


def expression_index(path: Path) -> tuple[dict[str, int], set[str]]:
    """Map expression identifiers to zero-based data-row indices."""

    mapping: dict[str, int] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if header[:2] != ["isoformId", "geneSymbol"]:
            raise ValueError(f"Unexpected expression header in {path}")
        for index, line in enumerate(handle):
            identifier = line.split("\t", 1)[0]
            if not identifier or identifier in mapping:
                raise ValueError(f"Blank or duplicate expression identifier in {path}")
            mapping[identifier] = index
    return mapping, set(mapping)


def panel_source_indices(
    expression_mapping: dict[str, int], panel_ids: Iterable[str]
) -> set[int]:
    identifiers = set(panel_ids)
    missing = identifiers - set(expression_mapping)
    if missing:
        raise ValueError(f"Panel has {len(missing)} IDs absent from expression matrix")
    return {expression_mapping[identifier] for identifier in identifiers}


def parse_dpi_stats(path: Path) -> dict[str, int | bool | float]:
    matches = []
    with path.open("r", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if line.startswith("[DPI_STATS]"):
                match = _DPI_STATS.fullmatch(line.rstrip("\r\n"))
                if match is None:
                    raise ValueError(f"Malformed [DPI_STATS] record in {path}")
                matches.append(match)
    if len(matches) != 1:
        raise ValueError(f"Expected one [DPI_STATS] record in {path}, got {len(matches)}")
    pre, pruned, post, applied = (int(value) for value in matches[0].groups())
    if not applied or pre != pruned + post:
        raise ValueError(f"Invalid DPI accounting in {path}")
    return {
        "pre_edges": pre,
        "pruned_edges": pruned,
        "post_edges": post,
        "dpi_applied": True,
        "pruned_fraction": 0.0 if pre == 0 else pruned / pre,
    }


def parse_sampling(path: Path) -> dict[str, object]:
    payloads: list[str] = []
    with path.open("r", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if line.startswith(_SAMPLING_PREFIX):
                payloads.append(line[len(_SAMPLING_PREFIX) :].strip())
    if len(payloads) != 1:
        raise ValueError(f"Expected one sampling record in {path}, got {len(payloads)}")
    indices = [int(value) for value in payloads[0].split()]
    if (
        len(indices) != 80
        or len(set(indices)) != 80
        or indices != sorted(indices)
        or any(value < 0 or value >= 100 for value in indices)
    ):
        raise ValueError(f"Invalid 80-of-100 sample in {path}")
    canonical = (" ".join(str(value) for value in indices) + "\n").encode("ascii")
    return {"indices": indices, "sha256": sha256_bytes(canonical)}


def validate_adjacency(
    path: Path, *, panel_ids: set[str], all_expression_ids: set[str]
) -> dict[str, object]:
    """Validate the output and hash data rows independently of headers."""

    headers: list[str] = []
    sources: set[str] = set()
    edges = 0
    data_digest = hashlib.sha256()
    with path.open("rb") as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(b">"):
                headers.append(raw.decode("utf-8").rstrip("\r\n"))
                continue
            data_digest.update(raw)
            fields = raw.decode("utf-8").rstrip("\r\n").split("\t")
            if len(fields) < 1 or len(fields) % 2 == 0:
                raise ValueError(f"Malformed adjacency row {path}:{line_number}")
            source = fields[0]
            if source in sources or source not in panel_ids:
                raise ValueError(f"Invalid/duplicate source {source!r} in {path}")
            sources.add(source)
            seen: set[str] = set()
            for index in range(1, len(fields), 2):
                target = fields[index]
                if target in seen or target not in all_expression_ids or target == source:
                    raise ValueError(f"Invalid/duplicate target {source}->{target} in {path}")
                seen.add(target)
                mi = float(fields[index + 1])
                if not math.isfinite(mi) or mi <= 0.0:
                    raise ValueError(f"Invalid MI for {source}->{target} in {path}")
                edges += 1
    required = (
        ">  MI threshold method estimator-matched AP-MI permutation-null GPD tail",
        ">  AP-MI null model m 80",
        ">  AP-MI null model Npar 40",
        ">  AP-MI cutoff tail extrapolated no",
        ">  Sampling method fixed-size without replacement",
        ">  Sampling request 80%",
        ">  Eligible observations 100",
        ">  Sampled observations 80",
        ">  MI threshold    0.164467",
        ">  MI P-value      0.0005",
        ">  DPI tolerance   0",
    )
    header_text = "\n".join(headers)
    for item in required:
        if item not in header_text:
            raise ValueError(f"Missing header {item!r} in {path}")
    if sources != panel_ids:
        raise ValueError(
            f"Adjacency source rows differ from panel in {path}: "
            f"{len(sources)} != {len(panel_ids)}"
        )
    return {
        "full_sha256": sha256_file(path),
        "data_sha256": data_digest.hexdigest(),
        "bytes": path.stat().st_size,
        "header_lines": len(headers),
        "source_rows": len(sources),
        "edges": edges,
    }


def parse_witness_sidecar(
    path: Path,
    *,
    expected_source_indices: set[int],
    expected_provenance: dict[str, str] | None = None,
) -> dict[str, object]:
    """Parse and strictly validate one candidate witness-count sidecar."""

    rows: dict[int, dict[str, int]] = {}
    provenance: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        header_line = ""
        for line_number, line in enumerate(handle, 1):
            if not line.startswith("# "):
                header_line = line
                break
            fields = line[2:].rstrip("\r\n").split("\t")
            if len(fields) != 2 or not fields[0] or fields[0] in provenance:
                raise ValueError(f"Malformed/duplicate provenance in {path}:{line_number}")
            provenance[fields[0]] = fields[1]
        if not header_line:
            raise ValueError(f"Missing witness table header: {path}")
        if tuple(provenance) != SIDECAR_PROVENANCE_KEYS:
            raise ValueError(
                f"Unexpected witness provenance keys/order in {path}: {tuple(provenance)}"
            )
        expected_static = {
            "schema": "sjaracne.dpi_witness_threshold_counts.v1",
            "graph_state": "unchanged pre-DPI edges after native K=1 marking",
            "count_unit": "source-target edges",
            "count_semantics": (
                "edges having at least the indicated number of distinct eligible intermediates"
            ),
            "source_index_basis": "zero-based expression-row index",
            "dpi_epsilon": "0",
            "source_mode": "selected rows",
            "source_count": str(len(expected_source_indices)),
            "annotated_gene_count": str(len(expected_source_indices)),
        }
        for key, expected in expected_static.items():
            if provenance[key] != expected:
                raise ValueError(
                    f"Unexpected witness provenance {key}={provenance[key]!r} in {path}"
                )
        if expected_provenance is not None:
            unknown = set(expected_provenance) - set(provenance)
            if unknown:
                raise ValueError(
                    f"Unknown expected witness provenance keys for {path}: {sorted(unknown)}"
                )
            for key, expected in expected_provenance.items():
                if provenance[key] != expected:
                    raise ValueError(
                        f"Witness provenance mismatch for {key} in {path}: "
                        f"expected {expected!r}, got {provenance[key]!r}"
                    )
        header = tuple(header_line.rstrip("\r\n").split("\t"))
        if header != SIDECAR_FIELDS:
            raise ValueError(f"Unexpected witness sidecar header in {path}: {header}")
        for line_number, line in enumerate(handle, 2):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(SIDECAR_FIELDS):
                raise ValueError(f"Malformed witness row {path}:{line_number}")
            try:
                values = [int(value) for value in fields]
            except ValueError as error:
                raise ValueError(f"Non-integer witness row {path}:{line_number}") from error
            source_index = values[0]
            counts = values[1:]
            if source_index < 0 or source_index in rows or any(value < 0 for value in counts):
                raise ValueError(f"Invalid/duplicate witness row {path}:{line_number}")
            pre = counts[0]
            witness_counts = counts[1:]
            if any(value > pre for value in witness_counts) or any(
                right > left for left, right in zip(witness_counts, witness_counts[1:])
            ):
                raise ValueError(f"Nonmonotone witness counts in {path}:{line_number}")
            rows[source_index] = dict(zip(SIDECAR_FIELDS[1:], counts))
    if set(rows) != expected_source_indices:
        missing = len(expected_source_indices - set(rows))
        extra = len(set(rows) - expected_source_indices)
        raise ValueError(
            f"Witness source-index coverage mismatch in {path}: missing={missing}, extra={extra}"
        )
    totals = {
        field: sum(row[field] for row in rows.values()) for field in SIDECAR_FIELDS[1:]
    }
    return {
        "sha256": sha256_file(path),
        "bytes": path.stat().st_size,
        "source_rows": len(rows),
        "provenance": provenance,
        "rows": rows,
        "totals": totals,
    }


def witness_field(k: int) -> str:
    if k not in K_VALUES:
        raise ValueError(f"Unsupported K_DPI value: {k}")
    return f"witnesses_ge_{k}"


def aggregate_sidecar_group(
    parsed: dict[str, object], source_indices: set[int], k: int
) -> dict[str, int | float]:
    rows = parsed["rows"]
    if not isinstance(rows, dict) or not source_indices.issubset(rows):
        raise ValueError("Requested source group is absent from sidecar")
    pre = sum(int(rows[index]["pre_edges"]) for index in source_indices)
    pruned = sum(int(rows[index][witness_field(k)]) for index in source_indices)
    if pruned > pre:
        raise ValueError("Aggregated pruned count exceeds pre-DPI count")
    return {
        "pre_edges": pre,
        "pruned_edges": pruned,
        "post_edges": pre - pruned,
        "pruned_fraction": 0.0 if pre == 0 else pruned / pre,
    }


def median(values: Iterable[float | int]) -> float:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        raise ValueError("Median requires data")
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2.0


def quartiles(values: Iterable[float | int]) -> tuple[float, float, float]:
    ordered = sorted(float(value) for value in values)
    if len(ordered) < 2:
        raise ValueError("Quartiles require at least two values")
    middle = len(ordered) // 2
    lower = ordered[:middle]
    upper = ordered[middle:] if len(ordered) % 2 == 0 else ordered[middle + 1 :]
    return median(lower), median(ordered), median(upper)


def ensure_exact_json(path: Path, payload: dict) -> None:
    if path.is_file():
        if load_json(path) != payload:
            raise RuntimeError(f"Existing frozen design differs: {path}")
        return
    atomic_json(path, payload)
