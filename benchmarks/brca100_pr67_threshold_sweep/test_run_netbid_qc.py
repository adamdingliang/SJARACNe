#!/usr/bin/env python3
"""Focused fail-closed tests for NetBID2 manifest provenance migration."""

from __future__ import annotations

import csv
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_netbid_qc as qc


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class NetbidManifestMigrationTest(unittest.TestCase):
    def make_fixture(self, root: Path) -> dict[str, object]:
        old_design = {
            "schema": qc.LEGACY_SWEEP_DESIGN_SCHEMA,
            "commit": "fixed-commit",
            "all_points": [{"key": "p_old", "p_value": 0.0001}],
            "fixed_parameters": {"m": 80},
        }
        new_design = {
            **old_design,
            "schema": qc.SWEEP_DESIGN_SCHEMA,
            "all_points": [
                *old_design["all_points"],
                {"key": "p_new", "p_value": 0.001},
            ],
        }
        old_bytes = qc.serialized_json(old_design)
        old_hash = qc.hashlib.sha256(old_bytes).hexdigest()
        active_path = root / "sweep_design.json"
        qc.atomic_json(active_path, new_design)
        new_hash = qc.sha256_file(active_path)
        history = root / qc.SWEEP_DESIGN_HISTORY_DIRECTORY
        archive_path = history / f"{old_hash}.sweep_design.json"
        archive_path.parent.mkdir(parents=True)
        archive_path.write_bytes(old_bytes)
        migration_path = history / f"{old_hash}_to_{new_hash}.migration.json"
        migration_payload = {
            "schema": qc.SWEEP_DESIGN_MIGRATION_SCHEMA,
            "operation": "append-only-point-extension",
            "from": {
                "schema": qc.LEGACY_SWEEP_DESIGN_SCHEMA,
                "sha256": old_hash,
                "archived_path": archive_path.relative_to(root).as_posix(),
                "point_keys": ["p_old"],
            },
            "to": {
                "schema": qc.SWEEP_DESIGN_SCHEMA,
                "sha256": new_hash,
                "active_path": "sweep_design.json",
                "point_keys": ["p_old", "p_new"],
            },
            "appended_points": [new_design["all_points"][1]],
            "manifest_path": migration_path.relative_to(root).as_posix(),
        }
        migration_payload["fingerprint"] = qc.fingerprint(migration_payload)
        qc.atomic_json(migration_path, migration_payload)
        migration = qc.load_sweep_design_migration(
            work_root=root,
            sweep_design=new_design,
            sweep_design_hash=new_hash,
        )
        self.assertIsNotNone(migration)

        arm = root / "results" / "p_old" / "tf"
        output = arm / "netbid2_qc"
        environment = {
            "R": "R synthetic",
            "NetBID2": "2.2.0",
            "NetBID2_remote_sha": "abc",
            "igraph": "2.3.3",
        }
        metrics: dict[str, object] = {metric: 0 for metric in qc.REQUIRED_METRICS}
        metrics.update(
            {
                "candidate_drivers": 2,
                "active_drivers": 1,
                "active_driver_fraction": 0.5,
                "edges": 1,
            }
        )
        write_tsv(
            output / "network_summary.tsv",
            ["metric", "value"],
            [{"metric": metric, "value": value} for metric, value in metrics.items()],
        )
        write_tsv(
            output / "driver_target_sizes.tsv",
            ["driver", "target_count"],
            [
                {"driver": "driver_a", "target_count": 1},
                {"driver": "driver_b", "target_count": 0},
            ],
        )
        write_tsv(
            output / "netbid_environment.tsv",
            ["component", "version"],
            [
                {"component": component, "version": version}
                for component, version in environment.items()
            ],
        )
        logs = arm / "logs"
        logs.mkdir(parents=True)
        stdout_path = logs / "netbid2_qc.stdout.log"
        stderr_path = logs / "netbid2_qc.stderr.log"
        stdout_path.write_text("synthetic output\n", encoding="utf-8")
        stderr_path.write_bytes(b"")
        manifest_path = arm / "netbid2_qc_manifest.json"
        command = [
            "wrapper",
            "Rscript",
            "script",
            "consensus",
            "drivers",
            "partial",
            "TF_",
            "false",
        ]
        current_payload = {
            "schema": qc.RUN_SCHEMA,
            "mode": "summary",
            "point": "p_old",
            "p_value": 0.0001,
            "mi_cutoff": 0.2,
            "point_manifest_sha256": "1" * 64,
            "sweep_design_sha256": new_hash,
            "driver": "tf",
            "driver_sha256": "2" * 64,
            "consensus_sha256": "3" * 64,
            "consensus_manifest_sha256": "4" * 64,
            "r_script_sha256": "5" * 64,
            "wrapper_sha256": "6" * 64,
            "environment": environment,
            "prefix": "TF_",
        }
        old_payload = dict(current_payload)
        old_payload["sweep_design_sha256"] = old_hash
        record = {
            **old_payload,
            "fingerprint": qc.fingerprint(old_payload),
            "command": command,
            "finished_at_utc": "2026-08-19T00:00:00+00:00",
            "output": str(output),
            "output_inventory": qc.inventory(output),
            "stdout_sha256": qc.sha256_file(stdout_path),
            "stderr_sha256": qc.sha256_file(stderr_path),
            "stderr_bytes": 0,
        }
        qc.atomic_json(manifest_path, record)
        return {
            "migration": migration,
            "manifest_path": manifest_path,
            "output": output,
            "record": record,
            "record_bytes": manifest_path.read_bytes(),
            "current_payload": current_payload,
            "command": command,
            "environment": environment,
            "stdout": stdout_path,
            "stderr": stderr_path,
            "old_hash": old_hash,
            "new_hash": new_hash,
        }

    def migrate(self, root: Path, fixture: dict[str, object]) -> dict[str, object]:
        return qc.migrate_completed_manifest_if_eligible(
            work_root=root,
            migration=fixture["migration"],
            point="p_old",
            driver="tf",
            mode="summary",
            manifest_path=fixture["manifest_path"],
            output_root=fixture["output"],
            record=qc.load_json(fixture["manifest_path"]),
            expected_payload=fixture["current_payload"],
            expected_command=fixture["command"],
            prefix="TF_",
            driver_ids=["driver_a", "driver_b"],
            expected_edges=1,
            stdout_path=fixture["stdout"],
            stderr_path=fixture["stderr"],
            expected_environment=fixture["environment"],
        )

    def test_exact_legacy_record_is_archived_and_refingerprinted(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = self.make_fixture(root)
            migrated = self.migrate(root, fixture)
            self.assertEqual(migrated["sweep_design_sha256"], fixture["new_hash"])
            self.assertEqual(
                migrated["fingerprint"], qc.fingerprint(fixture["current_payload"])
            )
            projected = dict(fixture["record"])
            projected["sweep_design_sha256"] = fixture["new_hash"]
            projected["fingerprint"] = qc.fingerprint(fixture["current_payload"])
            self.assertEqual(migrated, projected)

            pair = f"{fixture['old_hash']}_to_{fixture['new_hash']}"
            archive = (
                root
                / qc.NETBID_MANIFEST_HISTORY_DIRECTORY
                / pair
                / "arms/p_old/tf/netbid2_qc_manifest.json"
            )
            audit_path = archive.parents[3] / "migration.json"
            self.assertEqual(archive.read_bytes(), fixture["record_bytes"])
            audit = qc.load_json(audit_path)
            self.assertEqual(audit["fingerprint"], qc.fingerprint({
                key: value for key, value in audit.items() if key != "fingerprint"
            }))
            self.assertEqual(len(audit["migrated_runs"]), 1)

            # A crash after active rewrite but before audit promotion is
            # recoverable from the exact archive/current pair.
            audit_path.unlink()
            self.assertEqual(self.migrate(root, fixture), migrated)
            self.assertTrue(audit_path.is_file())

    def test_changed_log_rejects_before_any_migration_write(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = self.make_fixture(root)
            fixture["stdout"].write_text("changed\n", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "manifest validation failed"):
                self.migrate(root, fixture)
            self.assertEqual(fixture["manifest_path"].read_bytes(), fixture["record_bytes"])
            self.assertFalse((root / qc.NETBID_MANIFEST_HISTORY_DIRECTORY).exists())

    def test_changed_output_inventory_rejects_before_any_migration_write(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = self.make_fixture(root)
            summary = fixture["output"] / "network_summary.tsv"
            summary.write_bytes(summary.read_bytes() + b"\n")
            with self.assertRaisesRegex(RuntimeError, "manifest validation failed"):
                self.migrate(root, fixture)
            self.assertEqual(fixture["manifest_path"].read_bytes(), fixture["record_bytes"])
            self.assertFalse((root / qc.NETBID_MANIFEST_HISTORY_DIRECTORY).exists())

    def test_self_consistent_but_different_payload_is_not_adopted(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = self.make_fixture(root)
            record = qc.load_json(fixture["manifest_path"])
            record["driver_sha256"] = "9" * 64
            record["fingerprint"] = qc.fingerprint(qc.manifest_fingerprint_payload(record))
            qc.atomic_json(fixture["manifest_path"], record)
            with self.assertRaisesRegex(RuntimeError, "matches neither"):
                self.migrate(root, fixture)
            self.assertFalse((root / qc.NETBID_MANIFEST_HISTORY_DIRECTORY).exists())


if __name__ == "__main__":
    unittest.main()
