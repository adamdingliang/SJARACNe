#!/usr/bin/env python3
"""Shared, dependency-free helpers for the BRCA100 hub-size DPI pilot."""

from __future__ import annotations

import dataclasses
import hashlib
import json
import math
import os
from pathlib import Path
import re
from typing import Iterable, Iterator


PANEL_SEED = "brca100-hub-size-dpi-pilot-panel-v1"
FRACTIONS = (0.125, 0.5, 1.0)
K_MINIMUM_RECURRENCE = 6
SEEDS = tuple(range(1, 101))
EXPRESSION_FILENAME = "BRCA100.exp"


@dataclasses.dataclass(frozen=True)
class DriverSpec:
    key: str
    filename: str
    full_count: int
    counts: tuple[int, int, int]
    p_token: str
    p_value: float
    mi_cutoff: float


DRIVERS = (
    DriverSpec(
        key="tf",
        filename="BRCA100_TF.txt",
        full_count=2608,
        counts=(326, 1304, 2608),
        p_token="1e-3",
        p_value=1e-3,
        mi_cutoff=0.14732247558240297,
    ),
    DriverSpec(
        key="sig",
        filename="BRCA100_SIG.txt",
        full_count=10680,
        counts=(1335, 5340, 10680),
        p_token="5e-4",
        p_value=5e-4,
        mi_cutoff=0.1644671599536221,
    ),
)
DRIVER_BY_KEY = {driver.key: driver for driver in DRIVERS}


def arm_key(driver: DriverSpec, count: int) -> str:
    return f"{driver.key}_h{count:05d}"


def utc_now() -> str:
    import datetime as dt

    return dt.datetime.now(dt.timezone.utc).isoformat()


def serialized_json(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


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
        raise ValueError(f"Expected JSON object in {path}")
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
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def read_nonempty_unique_ids(path: Path, *, expected_count: int | None = None) -> list[str]:
    ids = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    if any(not item for item in ids):
        raise ValueError(f"Blank identifier in {path}")
    if len(set(ids)) != len(ids):
        raise ValueError(f"Duplicate identifier in {path}")
    if expected_count is not None and len(ids) != expected_count:
        raise ValueError(f"Expected {expected_count} IDs in {path}, got {len(ids)}")
    return ids


def expression_variances(path: Path, selected_ids: set[str]) -> dict[str, float]:
    """Return population variance over all samples for requested expression IDs."""

    result: dict[str, float] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if header[:2] != ["isoformId", "geneSymbol"] or len(header) < 4:
            raise ValueError(f"Unexpected BRCA100 expression header in {path}")
        sample_count = len(header) - 2
        for line_number, line in enumerate(handle, 2):
            fields = line.rstrip("\r\n").split("\t")
            accession = fields[0] if fields else ""
            if accession not in selected_ids:
                continue
            if len(fields) != len(header):
                raise ValueError(f"Malformed expression row {path}:{line_number}")
            values = [float(value) for value in fields[2:]]
            if any(not math.isfinite(value) for value in values):
                raise ValueError(f"Non-finite expression value at {path}:{line_number}")
            mean = math.fsum(values) / sample_count
            variance = math.fsum((value - mean) ** 2 for value in values) / sample_count
            result[accession] = variance
    missing = selected_ids - set(result)
    if missing:
        example = sorted(missing)[0]
        raise ValueError(f"Expression matrix lacks {len(missing)} selected IDs; e.g. {example}")
    return result


def variance_balanced_order(
    ids: list[str], variances: dict[str, float], *, driver_key: str
) -> tuple[list[str], dict[str, int]]:
    """Create one deterministic order interleaving five variance quintiles.

    Membership is randomized reproducibly within each rank-based quintile using
    SHA-256.  Round-robin interleaving keeps every prefix balanced to within one
    item per non-exhausted quintile.  Serialization later restores input order.
    """

    if set(ids) != set(variances):
        raise ValueError("Variance map and driver identifiers do not agree")
    ranked = sorted(ids, key=lambda accession: (variances[accession], accession))
    buckets: list[list[str]] = [[] for _ in range(5)]
    membership: dict[str, int] = {}
    for rank, accession in enumerate(ranked):
        quintile = min(4, (rank * 5) // len(ranked))
        buckets[quintile].append(accession)
        membership[accession] = quintile
    for quintile, bucket in enumerate(buckets):
        bucket.sort(
            key=lambda accession: (
                hashlib.sha256(
                    f"{PANEL_SEED}\0{driver_key}\0q{quintile}\0{accession}".encode(
                        "utf-8"
                    )
                ).digest(),
                accession,
            )
        )
    order: list[str] = []
    next_index = [0] * 5
    while len(order) < len(ids):
        progressed = False
        for quintile, bucket in enumerate(buckets):
            index = next_index[quintile]
            if index < len(bucket):
                order.append(bucket[index])
                next_index[quintile] += 1
                progressed = True
        if not progressed:
            raise RuntimeError("Variance-balanced panel ordering stalled")
    if len(order) != len(ids) or set(order) != set(ids):
        raise RuntimeError("Variance-balanced ordering is not a permutation")
    return order, membership


def create_panel_files(input_root: Path, panel_root: Path) -> dict[str, object]:
    """Create/validate deterministic nested hub panels and their manifest."""

    all_ids: set[str] = set()
    driver_ids: dict[str, list[str]] = {}
    for driver in DRIVERS:
        ids = read_nonempty_unique_ids(
            input_root / driver.filename, expected_count=driver.full_count
        )
        driver_ids[driver.key] = ids
        all_ids.update(ids)
    variances = expression_variances(input_root / EXPRESSION_FILENAME, all_ids)

    driver_records: dict[str, object] = {}
    for driver in DRIVERS:
        ids = driver_ids[driver.key]
        order, membership = variance_balanced_order(
            ids, {item: variances[item] for item in ids}, driver_key=driver.key
        )
        original_position = {accession: index for index, accession in enumerate(ids)}
        panel_records: list[dict[str, object]] = []
        previous_set: set[str] = set()
        for fraction, count in zip(FRACTIONS, driver.counts):
            selected = set(order[:count])
            if not previous_set.issubset(selected):
                raise RuntimeError(f"Non-nested generated panel for {driver.key}/{count}")
            previous_set = selected
            serialized_ids = sorted(selected, key=original_position.__getitem__)
            payload = ("\n".join(serialized_ids) + "\n").encode("utf-8")
            if count == driver.full_count:
                # Preserve the staged input bytes exactly (including its line
                # ending convention) so the 100% arm is a true anchor.
                payload = (input_root / driver.filename).read_bytes()
            relative = Path(driver.key) / f"h{count:05d}" / driver.filename
            path = panel_root / relative
            if path.is_file() and path.read_bytes() != payload:
                raise RuntimeError(f"Existing panel differs from fixed design: {path}")
            if not path.is_file():
                atomic_bytes(path, payload)
            quintile_counts = [0] * 5
            for accession in selected:
                quintile_counts[membership[accession]] += 1
            panel_records.append(
                {
                    "arm": arm_key(driver, count),
                    "fraction": fraction,
                    "hub_count": count,
                    "path": relative.as_posix(),
                    "sha256": sha256_bytes(payload),
                    "bytes": len(payload),
                    "variance_quintile_counts": quintile_counts,
                    "nested_selection_prefix_sha256": sha256_bytes(
                        ("\n".join(order[:count]) + "\n").encode("utf-8")
                    ),
                }
            )
        full_input = input_root / driver.filename
        if (panel_root / Path(driver.key) / f"h{driver.full_count:05d}" / driver.filename).read_bytes() != full_input.read_bytes():
            raise RuntimeError(
                f"Full {driver.key} panel is not byte-identical to staged input"
            )
        driver_records[driver.key] = {
            "source_filename": driver.filename,
            "source_sha256": sha256_file(full_input),
            "full_count": driver.full_count,
            "p_value": driver.p_value,
            "mi_cutoff": driver.mi_cutoff,
            "selection_order_sha256": sha256_bytes(
                ("\n".join(order) + "\n").encode("utf-8")
            ),
            "panels": panel_records,
        }

    manifest = {
        "schema": "sjaracne-brca100-hub-size-dpi-panel-v1",
        "panel_seed": PANEL_SEED,
        "selection": (
            "rank population variance over 100 BRCA100 samples into five quintiles; "
            "SHA-256 order within quintile; round-robin quintile interleaving; "
            "serialize selected IDs in original input order"
        ),
        "expression_sha256": sha256_file(input_root / EXPRESSION_FILENAME),
        "drivers": driver_records,
    }
    manifest_path = panel_root / "panel_manifest.json"
    if manifest_path.is_file() and load_json(manifest_path) != manifest:
        raise RuntimeError(f"Existing panel manifest differs: {manifest_path}")
    if not manifest_path.is_file():
        atomic_json(manifest_path, manifest)
    return manifest


_DPI_STATS_PATTERN = re.compile(
    r"^\[DPI_STATS\] pre_edges=(\d+) pruned_edges=(\d+) "
    r"post_edges=(\d+) dpi_applied=([01])$"
)
_SAMPLING_PREFIX = "[SAMPLING] Selected original observation indices (0-based):"


def parse_dpi_stats(path: Path, *, require_applied: bool = True) -> dict[str, int | bool]:
    matches: list[re.Match[str]] = []
    with path.open("r", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if line.startswith("[DPI_STATS]"):
                match = _DPI_STATS_PATTERN.fullmatch(line.rstrip("\r\n"))
                if match is None:
                    raise ValueError(f"Malformed [DPI_STATS] record in {path}: {line.rstrip()}")
                matches.append(match)
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one [DPI_STATS] record in {path}, got {len(matches)}")
    pre_edges, pruned_edges, post_edges, applied = (
        int(value) for value in matches[0].groups()
    )
    if pre_edges != pruned_edges + post_edges:
        raise ValueError(
            f"DPI accounting mismatch in {path}: {pre_edges} != "
            f"{pruned_edges} + {post_edges}"
        )
    if require_applied and applied != 1:
        raise ValueError(f"DPI was not applied according to {path}")
    return {
        "pre_edges": pre_edges,
        "pruned_edges": pruned_edges,
        "post_edges": post_edges,
        "dpi_applied": bool(applied),
        "pruned_fraction": 0.0 if pre_edges == 0 else pruned_edges / pre_edges,
    }


def parse_sampling_indices(path: Path) -> dict[str, object]:
    payloads: list[str] = []
    with path.open("r", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if line.startswith(_SAMPLING_PREFIX):
                payloads.append(line.rstrip("\r\n")[len(_SAMPLING_PREFIX) :].strip())
    if len(payloads) != 1:
        raise ValueError(
            f"Expected exactly one verbose sampling-index record in {path}, "
            f"got {len(payloads)}"
        )
    try:
        indices = [int(value) for value in payloads[0].split()]
    except ValueError as error:
        raise ValueError(f"Non-integer sampling index in {path}") from error
    if (
        len(indices) != 80
        or len(set(indices)) != 80
        or indices != sorted(indices)
        or any(index < 0 or index >= 100 for index in indices)
    ):
        raise ValueError(
            f"Expected 80 unique sorted sampling indices in [0,99] in {path}"
        )
    canonical = (" ".join(str(index) for index in indices) + "\n").encode("ascii")
    return {"indices": indices, "sha256": sha256_bytes(canonical)}


def iter_adjacency_edges(path: Path) -> Iterator[tuple[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith(">"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3 or len(fields) % 2 == 0:
                raise ValueError(f"Malformed adjacency row {path}:{line_number}")
            source = fields[0]
            for index in range(1, len(fields), 2):
                yield source, fields[index]


def median(values: Iterable[float | int]) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("Median requires at least one value")
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return float(ordered[middle])
    return (float(ordered[middle - 1]) + float(ordered[middle])) / 2.0


def quartiles(values: Iterable[float | int]) -> tuple[float, float, float]:
    """Tukey hinges: median of each half, suitable for the fixed n=100 design."""

    ordered = sorted(float(value) for value in values)
    if len(ordered) < 2:
        raise ValueError("Quartiles require at least two values")
    middle = len(ordered) // 2
    lower = ordered[:middle]
    upper = ordered[middle:] if len(ordered) % 2 == 0 else ordered[middle + 1 :]
    return median(lower), median(ordered), median(upper)
