"""Contracts for native consensus-before-DPI orchestration, without native mocks of science."""
import csv
import importlib
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest import mock

from test_wdl_workflow import work_directory
import test_wdl_workflow

WDL = importlib.import_module("SJARACNe.bin.wdl_workflow")
CLI = importlib.import_module("SJARACNe.sjaracne")
OPTIONS = importlib.import_module("SJARACNe.workflow_options")


def current(**extra):
    return OPTIONS.resolve_options({"network_type": "TF", **extra}, 100)


def native_output(final_edges=1):
    return (">  Consensus input runs B 2\n"
            ">  Consensus minimum recurrence K 1\n"
            ">  Consensus union edges 1\n"
            ">  Consensus candidate edges 1\n"
            ">  Consensus final edges {}\n"
            ">  Consensus source rows H 2\n"
            ">  DPI epsilon binary64 0.025000000000000001\n"
            "A\tB\t0.12345678901234567\nB\n").format(final_edges)


class TestConsensusFirstOptions(unittest.TestCase):
    def test_profile_defaults_and_previous_order_opt_out(self):
        first = current()
        self.assertEqual(first["consensus_order"], "consensus-first")
        self.assertEqual(first["dpi_offset"], 0.025)
        previous = current(consensus_order="dpi-first")
        self.assertEqual(previous["dpi_offset"], 0)
        legacy = OPTIONS.resolve_options({"profile": "legacy"}, 100)
        self.assertEqual(legacy["consensus_order"], "dpi-first")
        self.assertNotIn("dpi_offset", legacy)
        self.assertEqual(current(dpi_offset=0)["dpi_offset"], 0)

    def test_offset_and_order_validation(self):
        invalid = [dict(dpi_offset=v) for v in (True, "0.025", -0.01, 1, float("nan"), float("inf"))]
        invalid += [dict(consensus_order="unknown"), dict(p_value_consensus=0.01),
                    dict(dynamic_dpi=False, dpi_offset=0), dict(tolerance=0.1, dpi_offset=0.025)]
        for extra in invalid:
            with self.subTest(extra=extra), self.assertRaises(ValueError):
                current(**extra)
        with self.assertRaisesRegex(ValueError, "legacy profile"):
            OPTIONS.resolve_options({"profile": "legacy", "consensus_order": "consensus-first"}, 100)
        self.assertEqual(current(dpi_offset=0.1, dpi_max=0)["dpi_offset"], 0.1)
        self.assertNotIn("dpi_offset", current(tolerance=0.1))

    def test_standalone_cwl_scatter_count_mismatch_fails_before_staging(self):
        with work_directory() as root:
            expression, hubs = test_wdl_workflow.TestWdlTaskBoundaries.inputs(root)
            with mock.patch.object(WDL, "runtime_fingerprint") as fingerprint:
                with self.assertRaisesRegex(ValueError, "Scatter bootstrap_num differs"):
                    WDL.prepare({"exp_file": str(expression), "probe_file": str(hubs),
                        "bootstrap_num": 3,
                        "parameters": {"network_type": "TF", "bootstrap_num": 2, "min_recurrence": 1}})
                fingerprint.assert_not_called()
            self.assertFalse(Path("parameters.json").exists())
            self.assertFalse(Path("expression.exp").exists())

    def test_raw_scatter_defers_dpi_and_keeps_sampling(self):
        request = {"exp_file": "expression.exp", "probe_file": "hubs.txt", "seed": 1}
        for options in (current(), current(tolerance=0.1), current(dpi_offset=0)):
            with self.subTest(options=options):
                command = WDL.native_command(request, options)
                self.assertIn("--raw-mi", command)
                for flag in ("--dynamic-dpi", "--dpi-offset", "-e"):
                    self.assertNotIn(flag, command)
                self.assertEqual(command[command.index("-u") + 1], "80%")
        previous = WDL.native_command(request, current(consensus_order="dpi-first"))
        self.assertNotIn("--raw-mi", previous)
        self.assertEqual(previous[previous.index("--dpi-offset") + 1], "0.0")


class TestNativeConsensusAdapter(unittest.TestCase):
    def test_manifest_uses_existing_files_in_seed_order_without_copying(self):
        with work_directory() as root:
            files = [root / "input 'quoted'; 1.adj", root / "input space 2.adj"]
            for path in files:
                path.write_text("A\nB\n", encoding="utf-8")
            request = {"exp_file": "expression.exp", "probe_file": "hubs.txt", "adjacencies": list(map(str, files))}
            command = WDL.native_consensus_command(request, current(bootstrap_num=2, min_recurrence=1))
            self.assertEqual(Path("consensus_inputs.txt").read_text(encoding="utf-8").splitlines(), list(map(str, files)))
            self.assertFalse(Path("adjacencies").exists())
            self.assertEqual(command[command.index("--dpi-offset") + 1], "0.025")
            self.assertEqual(command[command.index("--min-recurrence") + 1], "1")
            self.assertEqual(command[command.index("-N") + 1], "40")
            for flag in ("-u", "-r", "-M", "-p", "--raw-mi"):
                self.assertNotIn(flag, command)

    def test_final_export_retains_precision_empty_hubs_and_metadata(self):
        with work_directory():
            Path("final.adj").write_text(native_output(), encoding="utf-8")
            network = WDL.export_native_consensus("final.adj", "results", current(bootstrap_num=2, min_recurrence=1))
            with network.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows, [{"source": "A", "target": "B", "MI": "0.12345678901234567"}])
            metadata = Path("results/parameter_info_.txt").read_text()
            self.assertIn("Consensus source rows H 2", metadata)
            self.assertIn("no MI rethresholding", metadata)
            self.assertIn('"p_value_bootstrap": 0.001', metadata)
            self.assertIn('"depth": 40', metadata)
            self.assertIn('"sampling_method": "without-replacement"', metadata)
            self.assertIn("MI aggregation: supporting-run arithmetic mean", Path("results/bootstrap_info_.txt").read_text())

    def test_export_rejects_inconsistent_native_completion_metadata(self):
        with work_directory():
            Path("final.adj").write_text(native_output(final_edges=2), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "metadata differs"):
                WDL.export_native_consensus("final.adj", "results", current(bootstrap_num=2, min_recurrence=1))

    def fixture(self, root):
        expression, hubs = test_wdl_workflow.TestWdlTaskBoundaries.inputs(root)
        fingerprint = {"native_sha256": "frozen", "python_runtime": {}}
        with mock.patch.object(WDL, "runtime_fingerprint", return_value=fingerprint):
            WDL.prepare({"exp_file": str(expression), "probe_file": str(hubs),
                         "parameters": {"network_type": "TF", "bootstrap_num": 2, "min_recurrence": 1}})
        settings = WDL.read_json("parameters.json")
        folder = root / "resamples"
        folder.mkdir()
        files, records = [], []
        for seed in (1, 2):
            path = folder / WDL.adjacency_name(seed)
            path.write_text("> raw\nA\tB\t0.12345678901234567\nB\n", encoding="utf-8")
            record = folder / ("seed%d.json" % seed)
            WDL.write_json(record, {"seed": seed, "returncode": 0, "runtime": fingerprint,
                "input_sha256": settings["input_sha256"], "prepared_settings_sha256": WDL.sha256("parameters.json"),
                "adjacency_sha256": WDL.sha256(path)})
            files.append(str(path)); records.append(str(record))
        return {"exp_file": str(root / "expression.exp"), "probe_file": str(root / "hubs.txt"),
                "settings": str(root / "parameters.json"), "adjacencies": files, "run_provenance": records}, fingerprint

    def test_consensus_native_completion_binds_final_output_hashes(self):
        with work_directory() as root:
            request, fingerprint = self.fixture(root)
            def native(command, stdout, stderr, check):
                self.assertFalse(check)
                Path(command[command.index("-o") + 1]).write_text(native_output(), encoding="utf-8")
                return SimpleNamespace(returncode=0)
            with mock.patch.object(WDL, "runtime_fingerprint", return_value=fingerprint), \
                    mock.patch.object(WDL.subprocess, "run", side_effect=native) as run:
                WDL.consensus(request)
            run.assert_called_once()
            self.assertFalse(Path("adjacencies").exists())
            value = WDL.read_json("consensus_provenance.json")
            self.assertEqual(value["returncode"], 0)
            self.assertEqual(value["parameters"]["dpi_offset"], 0.025)
            self.assertEqual(value["native_adjacency_sha256"], WDL.sha256("consensus.adj"))
            with Path("results/consensus_network_ncol_.txt").open(encoding="utf-8") as handle:
                enhanced = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(enhanced[0]["MI"], "0.12345678901234567")
            for filename, digest in value["output_sha256"].items():
                self.assertEqual(digest, WDL.sha256(Path("results") / filename))

    def test_failed_native_consensus_retains_evidence_and_does_not_enhance(self):
        with work_directory() as root:
            request, fingerprint = self.fixture(root)
            with mock.patch.object(WDL, "runtime_fingerprint", return_value=fingerprint), \
                    mock.patch.object(WDL.subprocess, "run", return_value=SimpleNamespace(returncode=9)):
                with self.assertRaisesRegex(RuntimeError, "Native consensus/DPI failed"):
                    WDL.consensus(request)
            self.assertEqual(WDL.read_json("consensus_provenance.json")["returncode"], 9)
            self.assertFalse(Path("results").exists())


class TestEnhancementMiPrecision(unittest.TestCase):
    def test_opt_in_preserves_token_and_default_keeps_legacy_rounding(self):
        with work_directory() as root:
            expression = root / "expression.exp"
            expression.write_text("isoformId\tgeneSymbol\ts1\ts2\ts3\nA\tAlpha\t1\t2\t4\nB\tBeta\t2\t3\t7\n", encoding="utf-8")
            network = root / "consensus_network_3col_.txt"
            token = "0.12345678901234567"
            network.write_text("source\ttarget\tMI\nA\tB\t" + token + "\n", encoding="utf-8")
            module = importlib.import_module("SJARACNe.bin.create_consensus_network")
            rows = []
            for name, options in (("old", {}), ("native", {"preserve_mi": True})):
                output = root / name
                module.create_enhanced_consensus_network(expression, network, output, **options)
                with (output / "consensus_network_ncol_.txt").open(encoding="utf-8") as handle:
                    rows.append(next(csv.DictReader(handle, delimiter="\t")))
            self.assertEqual(rows[0]["MI"], "0.1235")
            self.assertEqual(rows[1]["MI"], token)
            self.assertEqual({k: v for k, v in rows[0].items() if k != "MI"},
                             {k: v for k, v in rows[1].items() if k != "MI"})

    def test_preserved_mi_rejects_nonfinite_and_negative_tokens(self):
        for token in ("nan", "inf", "-0.01"):
            with self.subTest(token=token), work_directory() as root:
                expression = root / "expression.exp"
                expression.write_text("isoformId\tgeneSymbol\ts1\ts2\ts3\nA\tAlpha\t1\t2\t4\nB\tBeta\t2\t3\t7\n", encoding="utf-8")
                network = root / "consensus_network_3col_.txt"
                network.write_text("source\ttarget\tMI\nA\tB\t" + token + "\n", encoding="utf-8")
                module = importlib.import_module("SJARACNe.bin.create_consensus_network")
                with self.assertRaisesRegex(ValueError, "finite and nonnegative"):
                    module.create_enhanced_consensus_network(expression, network, root / "out", preserve_mi=True)


class TestNativeCwlPublication(unittest.TestCase):
    def fixture(self, root):
        engine, output = root / "engine", root / "published"
        engine.mkdir(); output.mkdir()
        parameters = current()
        for name in CLI.FINAL_OUTPUTS.values():
            (engine / name).write_text("output\n", encoding="utf-8")
        for name in ("consensus.native.stdout.txt", "consensus.native.stderr.txt"):
            (engine / name).write_text("", encoding="utf-8")
        hashes = {"exp_file": "expression-a", "probe_file": "hubs-a"}
        prepared = engine / "parameters.json"
        prepared.write_text(json.dumps({"schema": CLI.PARAMETER_SCHEMA,
            "runtime": {}, "parameters": parameters, "input_sha256": hashes,
            "original_input_sha256": hashes}), encoding="utf-8")
        receipt = {"schema": "sjaracne-consensus-v1", "returncode": 0,
                   "runtime": {}, "parameters": parameters,
                   "input_sha256": hashes, "prepared_settings_sha256": CLI.sha256(prepared),
                   "output_sha256": {name: CLI.sha256(engine / name) for name in CLI.FINAL_OUTPUTS.values()}}
        (engine / "consensus_provenance.json").write_text(json.dumps(receipt), encoding="utf-8")
        return engine, output, parameters, receipt, hashes

    def test_bad_output_hash_prevents_any_final_publication(self):
        with work_directory() as root:
            engine, output, parameters, receipt, hashes = self.fixture(root)
            receipt["output_sha256"] = {name: "wrong" for name in CLI.FINAL_OUTPUTS.values()}
            (engine / "consensus_provenance.json").write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "differs from consensus provenance"):
                CLI.publish_native_cwl_outputs(engine, output, {}, parameters, hashes)
            self.assertFalse(any(output.iterdir()))

    def test_foreign_prepared_receipt_or_original_inputs_prevent_publication(self):
        for case in ("prepared_sha", "receipt_inputs", "prepared_schema", "original_inputs"):
            with self.subTest(case=case), work_directory() as root:
                engine, output, parameters, receipt, hashes = self.fixture(root)
                if case == "prepared_sha":
                    receipt["prepared_settings_sha256"] = "foreign-dataset-prepared-settings"
                elif case == "receipt_inputs":
                    receipt["input_sha256"] = {"exp_file": "foreign-expression", "probe_file": "foreign-hubs"}
                elif case == "prepared_schema":
                    prepared = engine / "parameters.json"
                    value = json.loads(prepared.read_text())
                    value["schema"] = "foreign"
                    prepared.write_text(json.dumps(value), encoding="utf-8")
                    receipt["prepared_settings_sha256"] = CLI.sha256(prepared)
                else:
                    hashes = {"exp_file": "other-original-expression", "probe_file": "other-original-hubs"}
                (engine / "consensus_provenance.json").write_text(json.dumps(receipt), encoding="utf-8")
                with self.assertRaises(ValueError):
                    CLI.publish_native_cwl_outputs(engine, output, {}, parameters, hashes)
                self.assertFalse(any(output.iterdir()))


if __name__ == "__main__":
    unittest.main()
