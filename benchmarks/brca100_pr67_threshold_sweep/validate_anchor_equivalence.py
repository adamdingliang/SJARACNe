#!/usr/bin/env python3
"""Prove that the two sweep anchors reproduce their prior seed networks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path


PR67_COMMIT = "7633ebb4a0d966dbda15a4e32d0efa492fb71aeb"
PR66_COMMIT = "58091832848b2eaf2ae08f6f69482357b6b9b72c"
MODEL_SHA256 = "e3a8522682a8ea239821aaa10b12db72d00e07bfdcad43599d8e76a06be80944"
ANCHORS = (
    (
        "p1e-07",
        "pr67_7633ebb",
        PR67_COMMIT,
        True,
        1e-7,
        0.3224649956324025,
    ),
    (
        "p_pr66_cutoff_match",
        "pr66_5809183",
        PR66_COMMIT,
        False,
        0.000352804562601613,
        0.17280321515749669,
    ),
)
DRIVERS = ("tf", "sig")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def data_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for raw_line in handle:
            if not raw_line.startswith(b">"):
                digest.update(raw_line)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def option_value(command: object, option: str) -> str:
    if not isinstance(command, list) or command.count(option) != 1:
        raise ValueError(f"Expected exactly one {option} in command: {command!r}")
    index = command.index(option)
    if index + 1 >= len(command):
        raise ValueError(f"Missing value after {option}")
    return str(command[index + 1])


def atomic_json(path: Path, value: object) -> None:
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    os.replace(partial, path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sweep-work-root",
        type=Path,
        default=(
            Path.home()
            / "sjaracne-benchmarks"
            / "brca100-pr67-threshold-sweep-20260818"
        ),
    )
    parser.add_argument(
        "--prior-work-root",
        type=Path,
        default=(
            Path.home()
            / "sjaracne-benchmarks"
            / "brca100-netbid-qc-20260817-rerun"
        ),
    )
    parser.add_argument("--seed-start", type=int, default=1)
    parser.add_argument("--seed-end", type=int, default=100)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not 1 <= args.seed_start <= args.seed_end <= 100:
        raise ValueError("Seed range must be within 1..100")
    seeds = list(range(args.seed_start, args.seed_end + 1))
    for filename in ("BRCA100.exp", "BRCA100_TF.txt", "BRCA100_SIG.txt"):
        sweep_input = args.sweep_work_root / "inputs" / filename
        prior_input = args.prior_work_root / "inputs" / filename
        if sha256_file(sweep_input) != sha256_file(prior_input):
            raise RuntimeError(f"Sweep/prior input mismatch: {filename}")
    sweep_design_path = args.sweep_work_root / "sweep_design.json"
    sweep_design = load_json(sweep_design_path)
    if (
        sweep_design.get("schema") != "sjaracne-brca100-pr67-p-sweep-v2"
        or sweep_design.get("commit") != PR67_COMMIT
        or sweep_design.get("null_model_sha256") != MODEL_SHA256
        or sweep_design.get("fixed_parameters", {}).get("m") != 80
        or sweep_design.get("fixed_parameters", {}).get("npar") != 40
    ):
        raise RuntimeError(f"Unexpected sweep design provenance: {sweep_design_path}")
    design_points = sweep_design.get("all_points")
    if not isinstance(design_points, list):
        raise RuntimeError(f"Missing sweep-design point list: {sweep_design_path}")
    design_by_key = {
        str(record.get("key")): record
        for record in design_points
        if isinstance(record, dict)
    }
    if len(design_by_key) != len(design_points):
        raise RuntimeError(f"Invalid/duplicate sweep-design points: {sweep_design_path}")
    rows: list[dict[str, object]] = []
    point_manifest_hashes: dict[str, str] = {}
    for (
        sweep_point,
        prior_stage,
        prior_commit,
        prior_uses_model,
        expected_p,
        expected_cutoff,
    ) in ANCHORS:
        point_manifest_path = (
            args.sweep_work_root / "results" / sweep_point / "point_manifest.json"
        )
        point_manifest = load_json(point_manifest_path)
        design_point = design_by_key.get(sweep_point)
        if (
            point_manifest.get("key") != sweep_point
            or not isinstance(design_point, dict)
            or not math.isclose(
                float(point_manifest.get("p_value")), expected_p, rel_tol=0.0, abs_tol=0.0
            )
            or not math.isclose(
                float(point_manifest.get("mi_cutoff")),
                expected_cutoff,
                rel_tol=0.0,
                abs_tol=1e-15,
            )
            or point_manifest.get("p_value") != design_point.get("p_value")
            or point_manifest.get("mi_cutoff") != design_point.get("mi_cutoff")
        ):
            raise RuntimeError(f"Point manifest identity mismatch: {point_manifest_path}")
        point_manifest_hashes[sweep_point] = sha256_file(point_manifest_path)
        for driver in DRIVERS:
            for seed in seeds:
                stem = f"TF_run_{seed:03d}"
                sweep_root = args.sweep_work_root / "results" / sweep_point / driver
                prior_root = args.prior_work_root / "results" / prior_stage / driver
                sweep_metadata_path = sweep_root / "seed_metadata" / f"{stem}.json"
                prior_metadata_path = prior_root / "seed_metadata" / f"{stem}.json"
                sweep_adjacency = sweep_root / "adjacency" / f"{stem}.adj"
                prior_adjacency = prior_root / "adjacency" / f"{stem}.adj"
                for required in (
                    sweep_metadata_path,
                    prior_metadata_path,
                    sweep_adjacency,
                    prior_adjacency,
                ):
                    if not required.is_file():
                        raise FileNotFoundError(required)
                sweep_metadata = load_json(sweep_metadata_path)
                prior_metadata = load_json(prior_metadata_path)
                sweep_point_record = sweep_metadata.get("point")
                if (
                    not isinstance(sweep_point_record, dict)
                    or sweep_point_record.get("key") != sweep_point
                    or sweep_point_record.get("p_value") != point_manifest.get("p_value")
                    or sweep_point_record.get("mi_cutoff")
                    != point_manifest.get("mi_cutoff")
                    or sweep_metadata.get("driver") != driver
                    or sweep_metadata.get("seed") != seed
                    or prior_metadata.get("stage") != prior_stage
                    or prior_metadata.get("driver") != driver
                    or prior_metadata.get("seed") != seed
                ):
                    raise RuntimeError(
                        f"Seed metadata identity mismatch: "
                        f"{sweep_point}/{driver}/{seed:03d}"
                    )
                if (
                    sweep_metadata.get("commit") != PR67_COMMIT
                    or sweep_metadata.get("null_model_sha256") != MODEL_SHA256
                    or prior_metadata.get("commit") != prior_commit
                    or (prior_uses_model and prior_metadata.get("null_model_sha256") != MODEL_SHA256)
                    or (not prior_uses_model and prior_metadata.get("null_model_sha256") is not None)
                ):
                    raise RuntimeError(
                        f"Seed implementation provenance mismatch: "
                        f"{sweep_point}/{driver}/{seed:03d}"
                    )
                sweep_command = sweep_metadata.get("command")
                prior_command = prior_metadata.get("command")
                if (
                    float(option_value(sweep_command, "-p"))
                    != float(point_manifest["p_value"])
                    or option_value(prior_command, "-p") != "0.0000001"
                    or option_value(sweep_command, "-u") != "80%"
                    or option_value(prior_command, "-u") != "80%"
                    or option_value(sweep_command, "-N") != "40"
                    or option_value(prior_command, "-N") != "40"
                    or option_value(sweep_command, "-e") != "0"
                    or option_value(prior_command, "-e") != "0"
                    or option_value(sweep_command, "-a") != "adaptive_partitioning"
                    or option_value(prior_command, "-a") != "adaptive_partitioning"
                    or option_value(sweep_command, "-S") != str(seed)
                    or option_value(prior_command, "-S") != str(seed)
                    or not isinstance(sweep_command, list)
                    or sweep_command.count("-M") != 1
                    or not isinstance(prior_command, list)
                    or prior_command.count("-M") != int(prior_uses_model)
                ):
                    raise RuntimeError(
                        f"Seed command provenance mismatch: "
                        f"{sweep_point}/{driver}/{seed:03d}"
                    )
                sweep_data_hash = data_sha256(sweep_adjacency)
                prior_data_hash = data_sha256(prior_adjacency)
                for label, metadata, actual_hash in (
                    ("sweep", sweep_metadata, sweep_data_hash),
                    ("prior", prior_metadata, prior_data_hash),
                ):
                    recorded_hash = metadata.get("adjacency", {}).get("data_sha256")
                    if recorded_hash != actual_hash:
                        raise RuntimeError(
                            f"{label} metadata/data hash mismatch: "
                            f"{sweep_point}/{driver}/{seed:03d}"
                        )
                comparable_fields = ("edges", "source_rows", "mi_min", "mi_max")
                mismatches = [
                    field
                    for field in comparable_fields
                    if sweep_metadata["adjacency"].get(field)
                    != prior_metadata["adjacency"].get(field)
                ]
                if sweep_data_hash != prior_data_hash or mismatches:
                    raise RuntimeError(
                        f"Anchor mismatch {sweep_point}/{driver}/{seed:03d}: "
                        f"data_equal={sweep_data_hash == prior_data_hash}, "
                        f"fields={mismatches}"
                    )
                rows.append(
                    {
                        "sweep_point": sweep_point,
                        "prior_stage": prior_stage,
                        "driver": driver,
                        "seed": seed,
                        "data_sha256": sweep_data_hash,
                        "edges": sweep_metadata["adjacency"]["edges"],
                        "source_rows": sweep_metadata["adjacency"]["source_rows"],
                        "sweep_metadata_sha256": sha256_file(sweep_metadata_path),
                        "prior_metadata_sha256": sha256_file(prior_metadata_path),
                    }
                )

    output_root = args.sweep_work_root / "results" / "validation"
    output_root.mkdir(parents=True, exist_ok=True)
    table = output_root / "anchor_seed_equivalence.tsv"
    manifest_path = output_root / "anchor_seed_equivalence_manifest.json"
    # Never leave an older completion record beside a newly promoted table.
    if manifest_path.exists():
        manifest_path.unlink()
    temporary = table.with_name(table.name + ".partial")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, table)
    manifest = {
        "schema": "sjaracne-brca100-pr67-p-sweep-anchor-equivalence-v1",
        "anchors": [
            {
                "sweep_point": sweep,
                "prior_stage": prior,
                "expected_p": expected_p,
                "expected_mi_cutoff": expected_cutoff,
                "point_manifest_sha256": point_manifest_hashes[sweep],
            }
            for sweep, prior, _, _, expected_p, expected_cutoff in ANCHORS
        ],
        "drivers": list(DRIVERS),
        "seeds": seeds,
        "comparisons": len(rows),
        "all_data_sections_equal": True,
        "table": str(table),
        "table_sha256": sha256_file(table),
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "sweep_design_sha256": sha256_file(sweep_design_path),
        "sweep_work_root": str(args.sweep_work_root),
        "prior_work_root": str(args.prior_work_root),
    }
    atomic_json(manifest_path, manifest)
    print(f"[ANCHORS] {len(rows)} matched seed networks; all data sections equal")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
