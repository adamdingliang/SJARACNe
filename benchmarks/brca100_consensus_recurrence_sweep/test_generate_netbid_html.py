#!/usr/bin/env python3
"""Focused tests for the provenance-locked NetBID2 HTML adjunct."""

from __future__ import annotations

import csv
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import generate_netbid_html as html


class HtmlFixture:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.repo = root / "repo"
        self.source = root / "source"
        self.recurrence = root / "recurrence"
        self.overlay = root / "html-overlay"
        self.package = (
            self.repo
            / "benchmarks/brca100_consensus_recurrence_sweep/representative_netbid2_qc"
        )
        self.script_root = (
            self.repo / "benchmarks/brca100_consensus_recurrence_sweep"
        )
        self.r_script = (
            self.repo / "benchmarks/brca100_pr67_threshold_sweep/run_netbid_qc.R"
        )
        self.wrapper = self.repo / "benchmarks/brca100_netbid_qc/netbid2-r"
        self.driver_ids = {
            "tf": ("TF1", "TF2"),
            "sig": ("SIG1", "SIG2", "SIG3"),
        }
        self.edge_counts = {
            ("tf", 6): 3,
            ("sig", 6): 4,
            ("tf", 8): 2,
            ("sig", 8): 3,
        }
        self.patchers: list[mock._patch] = []
        self._build()

    @staticmethod
    def write(path: Path, payload: bytes) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)

    @staticmethod
    def json(path: Path, value: object) -> None:
        HtmlFixture.write(path, html.serialized_json(value))

    def _summary_bytes(self, driver: str, threshold: int) -> dict[str, bytes]:
        ids = self.driver_ids[driver]
        edges = self.edge_counts[(driver, threshold)]
        counts = [1] * len(ids)
        counts[0] += edges - len(ids)
        metrics = {
            name: 0 for name in html.REQUIRED_METRICS
        }
        metrics.update(
            {
                "candidate_drivers": len(ids),
                "active_drivers": sum(value > 0 for value in counts),
                "active_driver_fraction": 1,
                "edges": edges,
                "incident_nodes": len(ids) + 1,
                "weak_components": 1,
                "largest_weak_component": len(ids) + 1,
                "largest_weak_component_fraction": 1,
                "density": 0.1,
                "target_size_zero_mean": edges / len(ids),
                "target_size_zero_median": 1,
                "target_size_active_mean": edges / len(ids),
                "target_size_active_median": 1,
                "scale_free_adjusted_r2": 0.8 + threshold / 100,
            }
        )
        summary = "metric\tvalue\n" + "".join(
            f"{name}\t{metrics[name]}\n" for name in sorted(metrics)
        )
        targets = "driver\ttarget_count\n" + "".join(
            f"{driver_id}\t{count}\n" for driver_id, count in zip(ids, counts, strict=True)
        )
        environment = "component\tversion\n" + "".join(
            f"{component}\t{version}\n"
            for component, version in html.EXPECTED_ENVIRONMENT.items()
        )
        return {
            "network_summary.tsv": summary.encode(),
            "driver_target_sizes.tsv": targets.encode(),
            "netbid_environment.tsv": environment.encode(),
        }

    def _build(self) -> None:
        generator = self.script_root / "generate_netbid_html.py"
        self.write(generator, Path(html.__file__).read_bytes())
        frozen_runner = self.script_root / "run_recurrence_sweep.py"
        aggregator = self.script_root / "aggregate_recurrence.cpp"
        self.write(frozen_runner, b"fixture recurrence runner\n")
        self.write(aggregator, b"fixture aggregator\n")
        self.write(self.r_script, b"fixture NetBID R script\n")
        self.write(self.wrapper, b"fixture NetBID wrapper\n")
        frozen_script_hashes = {
            "aggregator_source_sha256": html.sha256_file(aggregator),
            "netbid_r_sha256": html.sha256_file(self.r_script),
            "runner_sha256": html.sha256_file(frozen_runner),
        }
        wrapper_hash = html.sha256_file(self.wrapper)

        driver_specs: dict[str, html.DriverSpec] = {}
        for driver, ids in self.driver_ids.items():
            filename = f"fixture_{driver}.txt"
            driver_path = self.source / "inputs" / filename
            self.write(driver_path, ("\n".join(ids) + "\n").encode())
            driver_specs[driver] = html.DriverSpec(
                filename=filename,
                prefix="TF_" if driver == "tf" else "SIG_",
                candidate_drivers=len(ids),
                sha256=html.sha256_file(driver_path),
                source_point="tf-point" if driver == "tf" else "sig-point",
                per_subsample_p=0.01 if driver == "tf" else 0.02,
            )

        run_specs: dict[tuple[str, int], html.RunSpec] = {}
        for driver, threshold in html.RUN_ORDER:
            arm = self.recurrence / "results" / driver / f"k{threshold:03d}"
            edges = self.edge_counts[(driver, threshold)]
            network = arm / "consensus_network_ncol_.txt"
            network_rows = "".join(
                f"{driver.upper()}{index}\tTARGET{index}\t0.5\n"
                for index in range(edges)
            ).encode()
            self.write(network, network_rows)
            network_hash = html.sha256_file(network)
            network_manifest = {
                "schema": html.NETWORK_SCHEMA,
                "driver": driver,
                "minimum_support": threshold,
                "edges": edges,
                "ncol_sha256": network_hash,
                "source_point": driver_specs[driver].source_point,
                "per_subsample_p": driver_specs[driver].per_subsample_p,
            }
            network_manifest_path = arm / "network_manifest.json"
            self.json(network_manifest_path, network_manifest)
            network_manifest_hash = html.sha256_file(network_manifest_path)

            summary_root = arm / "netbid2_qc"
            summary_outputs: dict[str, str] = {}
            for filename, payload in self._summary_bytes(driver, threshold).items():
                self.write(summary_root / filename, payload)
                summary_outputs[filename] = html.sha256_file(summary_root / filename)
            summary_manifest = {
                "schema": html.SUMMARY_SCHEMA,
                "driver": driver,
                "minimum_support": threshold,
                "input": {
                    "network_manifest_sha256": network_manifest_hash,
                    "network_sha256": network_hash,
                    "driver_sha256": driver_specs[driver].sha256,
                    "r_script_sha256": frozen_script_hashes["netbid_r_sha256"],
                    "wrapper_sha256": wrapper_hash,
                    "environment": html.EXPECTED_ENVIRONMENT,
                    "generate_html": False,
                },
                "outputs": summary_outputs,
            }
            summary_manifest_path = arm / "netbid2_manifest.json"
            self.json(summary_manifest_path, summary_manifest)
            run_specs[(driver, threshold)] = html.RunSpec(
                edges=edges,
                network_bytes=len(network_rows),
                network_sha256=network_hash,
                network_manifest_sha256=network_manifest_hash,
                summary_manifest_sha256=html.sha256_file(summary_manifest_path),
                summary_outputs=summary_outputs,
            )

        design_payload = {
            "schema": html.DESIGN_SCHEMA,
            "source_work_root": str(self.source.resolve()),
            "minimum_supports": list(range(6, 21)),
            "benchmark_scripts": frozen_script_hashes,
        }
        design = html.canonical_record(design_payload)
        design_path = self.recurrence / "design.json"
        self.json(design_path, design)

        frozen_package_fingerprint = "fixture-frozen-package"
        frozen_package_path = (
            self.script_root / "results_2026-08-20/package_manifest.json"
        )
        self.json(
            frozen_package_path,
            {
                "fingerprint": frozen_package_fingerprint,
                "design_fingerprint": design["fingerprint"],
                "live_work_root": str(self.recurrence.resolve()),
            },
        )

        replacements = {
            "DRIVER_SPECS": driver_specs,
            "RUN_SPECS": run_specs,
            "EXPECTED_FROZEN_SCRIPT_HASHES": frozen_script_hashes,
            "EXPECTED_WRAPPER_SHA256": wrapper_hash,
            "EXPECTED_DESIGN_FINGERPRINT": design["fingerprint"],
            "EXPECTED_DESIGN_SHA256": html.sha256_file(design_path),
            "EXPECTED_FROZEN_PACKAGE_FINGERPRINT": frozen_package_fingerprint,
            "EXPECTED_FROZEN_PACKAGE_SHA256": html.sha256_file(frozen_package_path),
        }
        for name, value in replacements.items():
            patcher = mock.patch.object(html, name, value)
            patcher.start()
            self.patchers.append(patcher)

    def close(self) -> None:
        for patcher in reversed(self.patchers):
            patcher.stop()

    def fake_subprocess(self, *, corrupt_arm: tuple[str, int] | None = None):
        environment_stdout = "component\tversion\n" + "".join(
            f"{component}\t{version}\n"
            for component, version in html.EXPECTED_ENVIRONMENT.items()
        )

        def invoke(command, **kwargs):
            if command[-1] == "--probe":
                return subprocess.CompletedProcess(command, 0, environment_stdout, "")
            network = Path(command[3])
            driver = network.parents[1].name
            threshold = int(network.parent.name[1:])
            output = Path(command[5])
            prefix = command[6]
            self.assert_command(command, driver, threshold)
            output.mkdir(parents=True)
            summary = network.parent / "netbid2_qc"
            for filename in html.SUMMARY_FILENAMES:
                shutil.copyfile(summary / filename, output / filename)
            if corrupt_arm == (driver, threshold):
                with (output / "network_summary.tsv").open("ab") as handle:
                    handle.write(b"tampered\t1\n")
            self.write(output / f"{prefix}netQC.Rmd", b"fixture Rmd\n")
            self.write(
                output / f"{prefix}netQC.html",
                f"<!doctype html><html><body>{driver} K={threshold}</body></html>\n".encode(),
            )
            if kwargs.get("stdout") is not None:
                kwargs["stdout"].write("NetBID2 QC complete\n")
            return subprocess.CompletedProcess(command, 0)

        return invoke

    def assert_command(self, command: list[str], driver: str, threshold: int) -> None:
        if command[-1] != "true" or command[6] != html.DRIVER_SPECS[driver].prefix:
            raise AssertionError(f"Unexpected fake NetBID command: {command}")
        expected_network = (
            self.recurrence
            / "results"
            / driver
            / f"k{threshold:03d}"
            / "consensus_network_ncol_.txt"
        )
        if Path(command[3]) != expected_network:
            raise AssertionError(f"Unexpected fake network: {command[3]}")

    def context(self) -> html.Context:
        return html.build_context(
            self.repo, self.source, self.recurrence, self.overlay
        )


class GenerateNetbidHtmlTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.fixture = HtmlFixture(Path(self.temporary.name))

    def tearDown(self) -> None:
        self.fixture.close()
        self.temporary.cleanup()

    def test_generates_resumes_recovers_and_packages_exact_inventory(self) -> None:
        fake = self.fixture.fake_subprocess()
        with mock.patch.object(html.subprocess, "run", side_effect=fake) as runner:
            context = self.fixture.context()
            aggregate = html.generate_html(context)
            html_calls = [
                call
                for call in runner.call_args_list
                if call.args[0][-1] == "true"
            ]
            self.assertEqual(len(html_calls), 4)
            self.assertEqual(aggregate["schema"], html.HTML_AGGREGATE_SCHEMA)
            self.assertEqual(len(aggregate["runs"]), 4)

            resumed = html.generate_html(context)
            self.assertEqual(resumed, aggregate)
            self.assertEqual(
                sum(call.args[0][-1] == "true" for call in runner.call_args_list),
                4,
            )

            driver, threshold = html.RUN_ORDER[0]
            paths = html.arm_paths(context, driver, threshold)
            shutil.copyfile(paths["manifest"], paths["pending"])
            paths["manifest"].unlink()
            completed_record = html.load_json(paths["pending"])
            owner = html.canonical_record(
                {
                    "schema": html.PARTIAL_OWNER_SCHEMA,
                    "pid": 99_999_999,
                    "lock_token": "dead-fixture-owner",
                    "driver": driver,
                    "minimum_support": threshold,
                    "input_fingerprint": completed_record["input_fingerprint"],
                }
            )
            html.atomic_json(paths["owner"], owner)
            recovered = html.generate_html(context)
            self.assertEqual(recovered, aggregate)
            self.assertTrue(paths["manifest"].is_file())
            self.assertFalse(paths["pending"].exists())
            self.assertEqual(
                sum(call.args[0][-1] == "true" for call in runner.call_args_list),
                4,
            )

            manifest = html.package_html(context, self.fixture.package, aggregate)
            self.assertEqual(manifest["schema"], html.HTML_PACKAGE_SCHEMA)
            self.assertEqual(len(manifest["reports"]), 4)
            for run_driver, run_threshold in html.RUN_ORDER:
                self.assertTrue(
                    (
                        self.fixture.package
                        / html.package_report_name(run_driver, run_threshold)
                    ).is_file()
                )
            self.assertTrue((self.fixture.package / "SHA256SUMS").is_file())
            self.assertTrue(
                (
                    self.fixture.package
                    / "provenance/scripts/generate_netbid_html.py"
                ).is_file()
            )
            package_partial = self.fixture.package.with_name(
                self.fixture.package.name + ".partial"
            )
            self.fixture.package.rename(package_partial)
            self.assertEqual(
                html.package_html(context, self.fixture.package, aggregate), manifest
            )
            self.assertFalse(package_partial.exists())
            self.assertEqual(
                html.package_html(context, self.fixture.package, aggregate), manifest
            )

        for driver, threshold in html.RUN_ORDER:
            frozen_arm = (
                self.fixture.recurrence / "results" / driver / f"k{threshold:03d}"
            )
            self.assertFalse((frozen_arm / "netbid2_qc_html").exists())
            self.assertFalse((frozen_arm / "netbid2_qc_html_manifest.json").exists())

    def test_html_mode_tsv_difference_fails_before_manifest_or_finalization(self) -> None:
        fake = self.fixture.fake_subprocess(corrupt_arm=("tf", 6))
        with mock.patch.object(html.subprocess, "run", side_effect=fake):
            context = self.fixture.context()
            with self.assertRaisesRegex(ValueError, "TSV bytes differ"):
                html.generate_html(context)
        paths = html.arm_paths(context, "tf", 6)
        self.assertTrue(paths["partial"].is_dir())
        self.assertFalse(paths["final"].exists())
        self.assertFalse(paths["manifest"].exists())

    def test_frozen_network_mutation_is_rejected_before_any_html_run(self) -> None:
        network = (
            self.fixture.recurrence
            / "results/tf/k006/consensus_network_ncol_.txt"
        )
        network.write_bytes(network.read_bytes() + b"X\tY\t0.1\n")
        fake = self.fixture.fake_subprocess()
        with mock.patch.object(html.subprocess, "run", side_effect=fake) as runner:
            with self.assertRaisesRegex(ValueError, "network SHA-256 mismatch"):
                self.fixture.context()
        self.assertEqual(runner.call_count, 1)  # environment probe only
        self.assertEqual(runner.call_args.args[0][-1], "--probe")

    def test_completed_html_or_package_tampering_fails_closed(self) -> None:
        fake = self.fixture.fake_subprocess()
        with mock.patch.object(html.subprocess, "run", side_effect=fake):
            context = self.fixture.context()
            aggregate = html.generate_html(context)
            paths = html.arm_paths(context, "sig", 8)
            report = paths["final"] / "SIG_netQC.html"
            report.write_bytes(report.read_bytes() + b"tamper")
            with self.assertRaisesRegex(ValueError, "inventory changed"):
                html.generate_html(context)

            report.write_bytes(report.read_bytes()[:-6])
            aggregate = html.validate_html_overlay(context)
            html.package_html(context, self.fixture.package, aggregate)
            packaged = self.fixture.package / "k006_tf_netbid2_qc.html"
            packaged.write_bytes(packaged.read_bytes() + b"tamper")
            with self.assertRaisesRegex(ValueError, "content inventory mismatch"):
                html.package_html(context, self.fixture.package, aggregate)

    def test_nested_symlink_cannot_redirect_writes_into_frozen_arm(self) -> None:
        fake = self.fixture.fake_subprocess()
        with mock.patch.object(html.subprocess, "run", side_effect=fake) as runner:
            context = self.fixture.context()
            link_parent = self.fixture.overlay / "runs/tf"
            link_parent.mkdir(parents=True)
            frozen_arm = self.fixture.recurrence / "results/tf/k006"
            try:
                (link_parent / "k006").symlink_to(
                    frozen_arm, target_is_directory=True
                )
            except OSError as exc:
                self.skipTest(f"directory symlinks unavailable: {exc}")
            with self.assertRaisesRegex(RuntimeError, "Symlink in managed"):
                html.generate_html(context)
        self.assertFalse((frozen_arm / "logs").exists())
        self.assertFalse((frozen_arm / "netbid2_qc_html").exists())
        self.assertEqual(
            sum(call.args[0][-1] == "true" for call in runner.call_args_list),
            0,
        )

    def test_live_overlay_lock_and_unowned_partial_fail_without_deletion(self) -> None:
        fake = self.fixture.fake_subprocess()
        with mock.patch.object(html.subprocess, "run", side_effect=fake) as runner:
            context = self.fixture.context()
            self.fixture.overlay.mkdir()
            lock = html.canonical_record(
                {
                    "schema": html.OVERLAY_LOCK_SCHEMA,
                    "pid": os.getpid(),
                    "token": "other-live-runner",
                    "design_fingerprint": context.design_fingerprint,
                    "generator_sha256": context.generator_sha256,
                    "started_at_utc": "2026-08-20T00:00:00+00:00",
                }
            )
            html.atomic_json(html.lock_path(context), lock)
            with self.assertRaisesRegex(RuntimeError, "owned by live PID"):
                html.generate_html(context)
            html.lock_path(context).unlink()

            partial = html.arm_paths(context, "tf", 6)["partial"]
            partial.mkdir(parents=True)
            sentinel = partial / "do-not-delete.txt"
            sentinel.write_text("unowned\n", encoding="utf-8")
            with self.assertRaisesRegex(FileNotFoundError, "partial owner"):
                html.generate_html(context)
            self.assertTrue(sentinel.is_file())
        self.assertEqual(
            sum(call.args[0][-1] == "true" for call in runner.call_args_list),
            0,
        )

    def test_runtime_script_mutation_is_detected_after_first_html_call(self) -> None:
        base_fake = self.fixture.fake_subprocess()
        mutated = False

        def mutate_after_call(command, **kwargs):
            nonlocal mutated
            result = base_fake(command, **kwargs)
            if command[-1] == "true" and not mutated:
                self.fixture.r_script.write_text(
                    "concurrent mutation\n", encoding="utf-8"
                )
                mutated = True
            return result

        with mock.patch.object(
            html.subprocess, "run", side_effect=mutate_after_call
        ):
            context = self.fixture.context()
            with self.assertRaisesRegex(ValueError, "R script SHA-256 mismatch"):
                html.generate_html(context)
        paths = html.arm_paths(context, "tf", 6)
        self.assertTrue(paths["partial"].is_dir())
        self.assertTrue(paths["owner"].is_file())
        self.assertFalse(paths["final"].exists())
        self.assertFalse(paths["manifest"].exists())

    def test_package_rejects_a_canonically_fingerprinted_forged_aggregate(self) -> None:
        fake = self.fixture.fake_subprocess()
        with mock.patch.object(html.subprocess, "run", side_effect=fake):
            context = self.fixture.context()
            aggregate = html.generate_html(context)
            payload = dict(aggregate)
            payload.pop("fingerprint")
            payload["selection"] = list(reversed(payload["selection"]))
            forged = html.canonical_record(payload)
            with self.assertRaisesRegex(ValueError, "Caller-supplied aggregate"):
                html.package_html(context, self.fixture.package, forged)
        self.assertFalse(self.fixture.package.exists())


if __name__ == "__main__":
    unittest.main()
