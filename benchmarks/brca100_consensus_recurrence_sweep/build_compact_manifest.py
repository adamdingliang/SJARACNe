#!/usr/bin/env python3
"""Build or verify the compact BRCA100 recurrence-evidence inventory."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
MANIFEST_PATH = ROOT / "PACKAGE_MANIFEST.json"
CHECKSUM_PATH = ROOT / "SHA256SUMS"
GENERATED_PATHS = {MANIFEST_PATH.name, CHECKSUM_PATH.name}
GENERATED_CACHE_DIRECTORIES = {"__pycache__", ".pytest_cache"}
GENERATED_CACHE_SUFFIXES = {".pyc", ".pyo"}
SCHEMA = "sjaracne-brca100-consensus-recurrence-compact-v1"
SOURCE_COMMIT = "7633ebb4a0d966dbda15a4e32d0efa492fb71aeb"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def payload_paths() -> list[Path]:
    paths: list[Path] = []
    for path in ROOT.rglob("*"):
        relative_path = path.relative_to(ROOT)
        relative = relative_path.as_posix()
        if path.is_symlink():
            raise ValueError(f"Symlinks are not allowed in the compact package: {relative}")
        if (
            GENERATED_CACHE_DIRECTORIES.intersection(relative_path.parts)
            or path.suffix in GENERATED_CACHE_SUFFIXES
        ):
            continue
        if path.is_file() and relative not in GENERATED_PATHS:
            paths.append(path)
    return sorted(paths, key=lambda path: path.relative_to(ROOT).as_posix())


def file_record(path: Path) -> dict[str, object]:
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def manifest_payload(records: list[dict[str, object]]) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "source_evidence_date": "2026-08-20",
        "source_sjaracne_commit": SOURCE_COMMIT,
        "scope": (
            "Compact review evidence for the BRCA100 minimum-recurrence "
            "consensus proposal; not a self-contained rerun bundle."
        ),
        "proposed_operating_points": {
            "tf": {
                "per_subsample_p_b": 0.001,
                "ap_mi_cutoff": 0.14732247558240297,
                "minimum_recurrence_k": 6,
                "subsamples_b": 100,
            },
            "sig": {
                "per_subsample_p_b": 0.0005,
                "ap_mi_cutoff": 0.1644671599536221,
                "minimum_recurrence_k": 6,
                "subsamples_b": 100,
            },
        },
        "selection": {
            "included": [
                "complete K=6..20 aggregate analysis and plot",
                "frozen matched-design and source-run provenance",
                (
                    "frozen aggregate exact-tail records and both sides of the "
                    "K=9 reproduction hash chain"
                ),
                (
                    "hash inventory for omitted materialized recurrence and "
                    "source-support artifacts"
                ),
                "representative K=6 TF and SIG NetBID2 HTML reports and provenance",
                "aggregation, packaging, report-generation scripts, and focused tests",
            ],
            "excluded": [
                "30 materialized consensus networks",
                "per-K NetBID2 tables and logs except the K=9 hash-chain records",
                "K=8 representative HTML reports",
                "source BRCA100 data and 200 post-DPI adjacency files",
                "generated caches and duplicate SVG output",
            ],
        },
        "file_count": len(records),
        "total_bytes": sum(int(record["bytes"]) for record in records),
        "files": records,
    }


def write_outputs() -> None:
    records = [file_record(path) for path in payload_paths()]
    manifest = manifest_payload(records)
    MANIFEST_PATH.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )

    checksummed_paths = payload_paths() + [MANIFEST_PATH]
    checksummed_paths.sort(key=lambda path: path.relative_to(ROOT).as_posix())
    lines = [
        f"{sha256_file(path)}  {path.relative_to(ROOT).as_posix()}\n"
        for path in checksummed_paths
    ]
    CHECKSUM_PATH.write_text("".join(lines), encoding="utf-8", newline="\n")


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def expected_manifest_records(manifest: dict[str, Any]) -> dict[str, dict[str, object]]:
    raw_records = manifest.get("files")
    if not isinstance(raw_records, list):
        raise ValueError("Manifest field 'files' must be a list")
    records: dict[str, dict[str, object]] = {}
    for record in raw_records:
        if not isinstance(record, dict):
            raise ValueError("Every manifest file record must be an object")
        relative = record.get("path")
        if not isinstance(relative, str) or not relative:
            raise ValueError("Every manifest file record needs a nonempty path")
        if relative in records:
            raise ValueError(f"Duplicate manifest path: {relative}")
        records[relative] = record
    return records


def verify_manifest() -> None:
    manifest = load_json(MANIFEST_PATH)
    if manifest.get("schema") != SCHEMA:
        raise ValueError(f"Unexpected compact-manifest schema: {manifest.get('schema')}")
    records = expected_manifest_records(manifest)
    actual_paths = {
        path.relative_to(ROOT).as_posix(): path for path in payload_paths()
    }
    if set(records) != set(actual_paths):
        missing = sorted(set(records) - set(actual_paths))
        unexpected = sorted(set(actual_paths) - set(records))
        raise ValueError(
            f"Manifest path mismatch; missing={missing}, unexpected={unexpected}"
        )
    for relative, path in actual_paths.items():
        record = records[relative]
        if record.get("bytes") != path.stat().st_size:
            raise ValueError(f"Manifest byte-size mismatch: {relative}")
        if record.get("sha256") != sha256_file(path):
            raise ValueError(f"Manifest SHA-256 mismatch: {relative}")
    if manifest.get("file_count") != len(records):
        raise ValueError("Manifest file_count does not match its inventory")
    total_bytes = sum(path.stat().st_size for path in actual_paths.values())
    if manifest.get("total_bytes") != total_bytes:
        raise ValueError("Manifest total_bytes does not match its inventory")


def parse_checksum_line(line: str, line_number: int) -> tuple[str, str]:
    parts = line.rstrip("\r\n").split("  ", 1)
    if len(parts) != 2 or len(parts[0]) != 64 or not parts[1]:
        raise ValueError(f"Malformed SHA256SUMS line {line_number}")
    try:
        int(parts[0], 16)
    except ValueError as error:
        raise ValueError(f"Non-hex digest on SHA256SUMS line {line_number}") from error
    return parts[0], parts[1]


def verify_checksums() -> None:
    expected_paths = payload_paths() + [MANIFEST_PATH]
    expected = {path.relative_to(ROOT).as_posix(): path for path in expected_paths}
    recorded: dict[str, str] = {}
    with CHECKSUM_PATH.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            digest, relative = parse_checksum_line(line, line_number)
            if relative in recorded:
                raise ValueError(f"Duplicate SHA256SUMS path: {relative}")
            recorded[relative] = digest
    if set(recorded) != set(expected):
        missing = sorted(set(expected) - set(recorded))
        unexpected = sorted(set(recorded) - set(expected))
        raise ValueError(
            f"SHA256SUMS path mismatch; missing={missing}, unexpected={unexpected}"
        )
    for relative, path in expected.items():
        if recorded[relative] != sha256_file(path):
            raise ValueError(f"SHA256SUMS mismatch: {relative}")


def verify() -> None:
    verify_manifest()
    verify_checksums()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--verify",
        action="store_true",
        help="verify the existing manifest and SHA256SUMS instead of rebuilding them",
    )
    args = parser.parse_args()
    if args.verify:
        verify()
        print("Compact manifest and SHA-256 checksums verified.")
    else:
        write_outputs()
        verify()
        print("Compact manifest and SHA-256 checksums written and verified.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
