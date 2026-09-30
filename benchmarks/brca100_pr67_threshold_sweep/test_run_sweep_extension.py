#!/usr/bin/env python3
"""Synthetic tests for the append-only sweep-design extension."""

from __future__ import annotations

import copy
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_sweep as sweep


class SweepDesignExtensionTest(unittest.TestCase):
    def synthetic_designs(
        self,
    ) -> tuple[dict[str, object], dict[str, object], dict[str, str], dict[str, object]]:
        build = {
            "binary_sha256": "b" * 64,
            "config_sha256": "c" * 64,
            "null_model_sha256": "d" * 64,
        }
        inputs: dict[str, object] = {
            "BRCA100.exp": {"sha256": "e" * 64, "bytes": 1},
            "BRCA100_TF.txt": {"sha256": "f" * 64, "bytes": 1},
            "BRCA100_SIG.txt": {"sha256": "a" * 64, "bytes": 1},
            "expression_id_count": {"count": 28278},
        }
        invariants: dict[str, object] = {
            "commit": sweep.PR67_COMMIT,
            **build,
            "fixed_parameters": {
                "sampling": "fixed 80% without replacement",
                "m": 80,
                "npar": 40,
                "dpi_epsilon": 0,
                "consensus_p": 1e-5,
                "seeds": list(range(1, 101)),
            },
            "inputs": inputs,
        }
        legacy_points = [
            {"key": point.key, "p_token": point.p_token, "p_value": point.p_value}
            for point in sweep.LEGACY_SWEEP_POINTS
        ]
        appended_points = [
            {"key": point.key, "p_token": point.p_token, "p_value": point.p_value}
            for point in sweep.EXTENDED_SWEEP_POINTS
        ]
        legacy = {
            "schema": sweep.LEGACY_SWEEP_DESIGN_SCHEMA,
            **invariants,
            "all_points": legacy_points,
        }
        extended = {
            "schema": sweep.SWEEP_DESIGN_SCHEMA,
            **invariants,
            "all_points": legacy_points + appended_points,
        }
        return legacy, extended, build, inputs

    def test_exact_v1_is_archived_and_migrated_once(self) -> None:
        legacy, extended, build, inputs = self.synthetic_designs()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design_path = root / "sweep_design.json"
            sweep.atomic_json(design_path, legacy)
            legacy_bytes = design_path.read_bytes()
            legacy_hash = sweep.sha256_bytes(legacy_bytes)
            extended_hash = sweep.sha256_bytes(sweep.serialized_json(extended))

            status = sweep.ensure_sweep_design(
                work_root=root,
                legacy_design=legacy,
                extended_design=extended,
                build=build,
                input_metadata=inputs,
            )

            self.assertEqual(status, "migrated-v1-to-v2")
            self.assertEqual(sweep.load_json(design_path), extended)
            history = root / sweep.SWEEP_DESIGN_HISTORY_DIRECTORY
            archive = history / f"{legacy_hash}.sweep_design.json"
            migration_path = history / f"{legacy_hash}_to_{extended_hash}.migration.json"
            self.assertEqual(archive.read_bytes(), legacy_bytes)
            migration = sweep.load_json(migration_path)
            fingerprint = migration.pop("fingerprint")
            self.assertEqual(fingerprint, sweep.json_fingerprint(migration))
            self.assertEqual(migration["from"]["sha256"], legacy_hash)
            self.assertEqual(migration["to"]["sha256"], extended_hash)
            self.assertEqual(
                [point["key"] for point in migration["appended_points"]],
                [point.key for point in sweep.EXTENDED_SWEEP_POINTS],
            )

            design_bytes = design_path.read_bytes()
            self.assertEqual(
                sweep.ensure_sweep_design(
                    work_root=root,
                    legacy_design=legacy,
                    extended_design=extended,
                    build=build,
                    input_metadata=inputs,
                ),
                "existing-v2",
            )
            self.assertEqual(design_path.read_bytes(), design_bytes)

    def test_mutated_v1_fails_without_history_or_design_change(self) -> None:
        legacy, extended, build, inputs = self.synthetic_designs()
        mutated = copy.deepcopy(legacy)
        mutated["fixed_parameters"]["consensus_p"] = 2e-5
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design_path = root / "sweep_design.json"
            sweep.atomic_json(design_path, mutated)
            original = design_path.read_bytes()

            with self.assertRaisesRegex(RuntimeError, "Incompatible existing"):
                sweep.ensure_sweep_design(
                    work_root=root,
                    legacy_design=legacy,
                    extended_design=extended,
                    build=build,
                    input_metadata=inputs,
                )

            self.assertEqual(design_path.read_bytes(), original)
            self.assertFalse((root / sweep.SWEEP_DESIGN_HISTORY_DIRECTORY).exists())
            self.assertFalse((root / "results").exists())

    def test_incompatible_legacy_point_manifest_blocks_migration(self) -> None:
        legacy, extended, build, inputs = self.synthetic_designs()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design_path = root / "sweep_design.json"
            sweep.atomic_json(design_path, legacy)
            original = design_path.read_bytes()
            point = legacy["all_points"][0]
            manifest = sweep.point_manifest_payload(
                design=point, build=build, input_metadata=inputs
            )
            manifest["consensus_p"] = 2e-5
            point_path = root / "results" / str(point["key"]) / "point_manifest.json"
            sweep.atomic_json(point_path, manifest)

            with self.assertRaisesRegex(RuntimeError, "legacy point manifest"):
                sweep.ensure_sweep_design(
                    work_root=root,
                    legacy_design=legacy,
                    extended_design=extended,
                    build=build,
                    input_metadata=inputs,
                )

            self.assertEqual(design_path.read_bytes(), original)
            self.assertFalse((root / sweep.SWEEP_DESIGN_HISTORY_DIRECTORY).exists())

    def test_fresh_v2_refuses_orphan_point_manifest(self) -> None:
        legacy, extended, build, inputs = self.synthetic_designs()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            orphan = root / "results" / "p4e-04" / "point_manifest.json"
            sweep.atomic_json(orphan, {"orphan": True})

            with self.assertRaisesRegex(RuntimeError, "without a sweep design"):
                sweep.ensure_sweep_design(
                    work_root=root,
                    legacy_design=legacy,
                    extended_design=extended,
                    build=build,
                    input_metadata=inputs,
                )

            self.assertFalse((root / "sweep_design.json").exists())
            self.assertFalse((root / sweep.SWEEP_DESIGN_HISTORY_DIRECTORY).exists())


if __name__ == "__main__":
    unittest.main()
