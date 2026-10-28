"""Ordinary task-plumbing regressions; no WDL runner or installed binary needed."""
import contextlib
import importlib
import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock


WDL = importlib.import_module("SJARACNe.bin.wdl_workflow")


@contextlib.contextmanager
def work_directory():
    previous = Path.cwd()
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary).resolve()
        os.chdir(root)
        try:
            yield root
        finally:
            os.chdir(previous)


class TestWdlParameters(unittest.TestCase):
    def test_workflow_defaults_and_null_inputs(self):
        values = WDL.validate_parameters({}, 100)
        for name, expected in {
            "bootstrap_num": 100, "subsample_spec": "80%", "depth": 40,
            "p_value_bootstrap": 1e-7, "apmi_null_model": None,
            "dynamic_dpi": False, "tolerance": 0,
            "min_recurrence": 6, "p_value_consensus": None,
        }.items():
            self.assertEqual(values[name], expected, name)
        self.assertEqual(values, WDL.validate_parameters(
            {name: None for name in WDL.OPTION_NAMES}, 100))

    def test_pchip_requires_explicit_probability_without_tf_sig_default(self):
        for probability in ({}, {"p_value_bootstrap": None}):
            with self.subTest(probability=probability), self.assertRaisesRegex(
                    ValueError, "PCHIP requires explicit"):
                WDL.validate_parameters({"apmi_null_model": "pchip", **probability}, 100)
        for probability in (0.001, 0.0005):
            values = WDL.validate_parameters({"apmi_null_model": "pchip",
                                               "p_value_bootstrap": probability}, 100)
            self.assertEqual(values["p_value_bootstrap"], probability)
            self.assertEqual(values["subsample_spec"], "80%")

    def test_pchip_uses_resolved_m_and_native_ceil_operation(self):
        cases = [(100, "80%", 80), (517, "80%", 414),
                 (805, "80%", 644), (101, "73.5%", 75), (100, "80", 80)]
        for population, spec, expected in cases:
            with self.subTest(population=population, spec=spec):
                with mock.patch("SJARACNe.apmi_cutoff.get_mi_cutoff", return_value=0.1) as cutoff:
                    values = WDL.validate_parameters({"apmi_null_model": "pchip",
                        "p_value_bootstrap": 0.001, "subsample_spec": spec}, population)
                cutoff.assert_called_once_with(expected, 0.001)
                self.assertEqual(values["subsample_spec"], spec)

    def test_pchip_bounds_are_checked_by_shared_callable(self):
        for selected in (50, 80, 10000):
            with self.subTest(selected=selected):
                self.assertEqual(WDL.validate_parameters({"apmi_null_model": "pchip",
                    "p_value_bootstrap": 0.001, "subsample_spec": str(selected)}, selected)["depth"], 40)
        invalid = [({"subsample_spec": "49"}, 100),
                   ({"subsample_spec": "10001"}, 10001),
                   ({"depth": 20}, 100), ({"p_value_bootstrap": 1e-8}, 100),
                   ({"p_value_bootstrap": 0.1}, 100)]
        for extra, population in invalid:
            with self.subTest(extra=extra), self.assertRaises(ValueError):
                WDL.validate_parameters({"apmi_null_model": "pchip",
                    "p_value_bootstrap": 0.001, **extra}, population)

    def test_probability_one_bypasses_pchip_conversion(self):
        with mock.patch("SJARACNe.apmi_cutoff.get_mi_cutoff",
                        side_effect=AssertionError("p=1 must not request a cutoff")):
            values = WDL.validate_parameters({"apmi_null_model": "pchip",
                "p_value_bootstrap": 1, "depth": 1, "subsample_spec": "2"}, 2)
        self.assertEqual(values["p_value_bootstrap"], 1)
        self.assertEqual(values["depth"], 1)

    def test_consensus_default_and_explicit_legacy_selection(self):
        values = WDL.validate_parameters({"bootstrap_num": 2, "min_recurrence": 1}, 100)
        self.assertEqual(values["min_recurrence"], 1)
        values = WDL.validate_parameters({"bootstrap_num": 2, "p_value_consensus": 1e-5}, 100)
        self.assertIsNone(values["min_recurrence"])
        self.assertEqual(values["p_value_consensus"], 1e-5)
        for options in ({"bootstrap_num": 2},
                        {"bootstrap_num": 2, "min_recurrence": 3},
                        {"min_recurrence": 6, "p_value_consensus": 0.01}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                WDL.validate_parameters(options, 100)

    def test_dynamic_conflicts_and_tuning_requires_opt_in(self):
        for fixed in (0, 0.1, 1):
            with self.subTest(fixed=fixed), self.assertRaisesRegex(ValueError, "conflicts"):
                WDL.validate_parameters({"dynamic_dpi": True, "tolerance": fixed}, 100)
        for key, value in (("dpi_href", 479), ("dpi_slope", 0.05), ("dpi_max", 0.3)):
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "require dynamic_dpi"):
                WDL.validate_parameters({key: value}, 100)
        values = WDL.validate_parameters({"dynamic_dpi": True}, 100)
        self.assertNotIn("tolerance", values)
        for key in ("dpi_href", "dpi_slope", "dpi_max"):
            self.assertNotIn(key, values, "Native executable must own omitted tuning defaults")
        values = WDL.validate_parameters({"dynamic_dpi": True, "dpi_href": 2,
            "dpi_slope": 0.1, "dpi_max": 0}, 100)
        self.assertEqual(values["dpi_max"], 0)

    def test_invalid_parameter_types_ranges_and_nonfinite_values(self):
        cases = {
            "bootstrap_num": (True, 0, -1, 1.0, "100", 2147483648),
            "depth": (True, 0, -1, 40.0, "40"),
            "min_recurrence": (True, 0, -1, 1.0, 101),
            "p_value_bootstrap": (True, 0, -1, 1.1, "0.001", float("nan"), float("inf")),
            "p_value_consensus": (True, 0, -1, 1.1, float("nan"), float("inf")),
            "tolerance": (True, -0.1, 1.1, "0", float("nan"), float("inf")),
            "dynamic_dpi": (0, 1, "true"),
            "apmi_null_model": (True, "gpd", "legacy", "PCHIP"),
        }
        for name, invalid_values in cases.items():
            for value in invalid_values:
                with self.subTest(name=name, value=value), self.assertRaises(ValueError):
                    WDL.validate_parameters({name: value}, 100)
        for name, invalid_values in {
            "dpi_href": (True, 0, -1, 2.0, 2147483648),
            "dpi_slope": (True, 0, -1, float("nan"), float("inf")),
            "dpi_max": (True, -1, 1, float("nan"), float("inf")),
        }.items():
            for value in invalid_values:
                with self.subTest(name=name, value=value), self.assertRaises(ValueError):
                    WDL.validate_parameters({"dynamic_dpi": True, name: value}, 100)
        for parameters in ([], None, {"unknown": 1}):
            with self.subTest(parameters=parameters), self.assertRaises(ValueError):
                WDL.validate_parameters(parameters, 100)

    def test_subsampling_count_percentage_and_invalid_forms(self):
        for spec in ("2", "100", "100%", "80.0%", ".8%", "00100"):
            population = 1000 if spec == ".8%" else 100
            self.assertEqual(WDL.validate_parameters({"subsample_spec": spec}, population)["subsample_spec"], spec)
        for spec in (True, 80, "", "0", "1", "101", "0%", "100.1%", "nan%",
                     "inf%", "80 %", " 80%", "-80", "1e2", "1e2%", "80.0"):
            with self.subTest(spec=spec), self.assertRaises(ValueError):
                WDL.validate_parameters({"subsample_spec": spec}, 100)


class TestWdlNativeCommand(unittest.TestCase):
    @staticmethod
    def command(options=None, request=None):
        values = WDL.validate_parameters(options or {}, 100)
        request = {"exp_file": "/expression with space;$(literal).exp",
                   "probe_file": "/hubs 'quoted';literal.txt", "seed": 1, **(request or {})}
        return WDL.native_command(request, values)

    def test_sampling_and_scientific_options_map_to_native_flags(self):
        command = self.command()
        self.assertEqual(command[command.index("-s") + 1], command[command.index("-l") + 1])
        for flag, expected in (("-u", "80%"), ("-N", "40"), ("-p", "1e-07"),
                               ("-a", "adaptive_partitioning"), ("-S", "1"),
                               ("-o", "TF_run_001.adj")):
            self.assertEqual(command[command.index(flag) + 1], expected)
        self.assertNotIn("-r", command)
        self.assertNotIn("-M", command)
        self.assertEqual(Path(command[0]), WDL.PACKAGE / "bin/sjaracne.exe")

    def test_paths_with_shell_characters_are_literal_argv_elements(self):
        request = {"exp_file": "/input 'quoted';$(do_not_run).exp",
                   "probe_file": '/hubs `literal` "quoted".txt',
                   "aracne_config_dir": "/config ; $(literal)"}
        command = self.command(request=request)
        for flag, key in (("-i", "exp_file"), ("-s", "probe_file"),
                          ("-l", "probe_file"), ("-H", "aracne_config_dir")):
            self.assertEqual(command[command.index(flag) + 1], request[key])
        self.assertIsInstance(command, list)

    def test_fixed_zero_fixed_one_and_dynamic_flags(self):
        for tolerance in (0, 1):
            command = self.command({"tolerance": tolerance})
            self.assertEqual(float(command[command.index("-e") + 1]), tolerance)
            self.assertNotIn("--dynamic-dpi", command)
        command = self.command({"dynamic_dpi": True})
        self.assertIn("--dynamic-dpi", command)
        for flag in ("-e", "--dpi-href", "--dpi-slope", "--dpi-max"):
            self.assertNotIn(flag, command)
        command = self.command({"dynamic_dpi": True, "dpi_href": 2,
            "dpi_slope": 0.1, "dpi_max": 0, "apmi_null_model": "pchip",
            "p_value_bootstrap": 0.001, "subsample_spec": "80"})
        for flag, expected in (("--dpi-href", "2"), ("--dpi-slope", "0.1"),
                               ("--dpi-max", "0"), ("-M", "pchip"), ("-u", "80")):
            self.assertEqual(command[command.index(flag) + 1], expected)
        self.assertNotIn("-e", command)

    def test_seed_names_do_not_repeat_above_999(self):
        expected = {1: "TF_run_001.adj", 9: "TF_run_009.adj", 999: "TF_run_999.adj",
                    1000: "TF_run_1000.adj", 1001: "TF_run_1001.adj"}
        self.assertEqual({seed: WDL.adjacency_name(seed) for seed in expected}, expected)
        for seed in (True, 0, -1, 1.0, "1", 2147483648):
            with self.subTest(seed=seed), self.assertRaises(ValueError):
                WDL.adjacency_name(seed)
        for seed in (0, True, 101):
            with self.subTest(seed=seed), self.assertRaises(ValueError):
                self.command(request={"seed": seed})


class TestWdlFiles(unittest.TestCase):
    def test_lf_normalization_preserves_bytes_order_and_source(self):
        block = 1024 * 1024
        cases = [b"already\nLF\n", b"a\r\nb\rc\nlast\r", b"", b"\r\n\r\n",
                 b"x" * (block - 1) + b"\r\ny\rfinal",
                 b"x" * (block - 1) + b"\ry\r\nfinal\r",
                 b"\xef\xbb\xbfheader\r\n\xc3\xa9\t1\r"]
        with work_directory() as root:
            for number, content in enumerate(cases):
                with self.subTest(number=number):
                    source, destination = root / (str(number) + ".in"), root / (str(number) + ".out")
                    source.write_bytes(content)
                    WDL.normalize_lf(source, destination)
                    self.assertEqual(source.read_bytes(), content)
                    self.assertEqual(destination.read_bytes(), content.replace(b"\r\n", b"\n").replace(b"\r", b"\n"))

    def test_normalization_never_overwrites_destination_or_source(self):
        with work_directory() as root:
            source, destination = root / "input", root / "output"
            source.write_bytes(b"A\r\nB\r")
            destination.write_bytes(b"keep")
            for target in (source, destination):
                with self.subTest(target=target), self.assertRaises(FileExistsError):
                    WDL.normalize_lf(source, target)
            self.assertEqual(source.read_bytes(), b"A\r\nB\r")
            self.assertEqual(destination.read_bytes(), b"keep")

    def test_json_requires_object_and_rejects_nonfinite_constants(self):
        with work_directory() as root:
            path = root / "request.json"
            for text in ('[]', 'null', '1', '{"value": NaN}',
                         '{"value": Infinity}', '{"value": -Infinity}'):
                with self.subTest(text=text):
                    path.write_text(text, encoding="utf-8")
                    with self.assertRaises(ValueError):
                        WDL.read_json(path)
            path.write_text('{"path": "spaces ; $(literal)", "n": 1}', encoding="utf-8")
            self.assertEqual(WDL.read_json(path), {"path": "spaces ; $(literal)", "n": 1})


class TestWdlTaskBoundaries(unittest.TestCase):
    @staticmethod
    def inputs(root, observations=100):
        folder = root / "input files"
        folder.mkdir()
        expression, hubs = folder / "expression.exp", folder / "hubs.txt"
        expression.write_bytes(("isoformId\tgeneSymbol\t" + "\t".join("s%d" % i for i in range(observations))
            + "\r\nA\tSYMBOL_A\t" + "\t".join(str(i) for i in range(observations))
            + "\r\nB\tSYMBOL_B\t" + "\t".join(str(observations-i) for i in range(observations))
            + "\r\n").encode("ascii"))
        hubs.write_bytes(b"A\r\nA\r\nMISSING\r\n")
        return expression, hubs

    def prepare(self, root, options=None):
        expression, hubs = self.inputs(root)
        request = {"exp_file": str(expression), "probe_file": str(hubs), "parameters": options or {}}
        fingerprint = {"native_sha256": "frozen-binary", "python_runtime": {"module.py": "frozen"}}
        with mock.patch.object(WDL, "runtime_fingerprint", return_value=fingerprint):
            WDL.prepare(request)
        return request, fingerprint

    def test_prepare_validates_overlap_normalizes_and_preserves_originals(self):
        with work_directory() as root:
            request, fingerprint = self.prepare(root)
            settings = WDL.read_json(root / "parameters.json")
            self.assertEqual(settings["schema"], WDL.SCHEMA)
            self.assertEqual(settings["observation_count"], 100)
            self.assertEqual(settings["runtime"], fingerprint)
            for key, normalized in (("exp_file", "expression.exp"), ("probe_file", "hubs.txt")):
                original = Path(request[key]).read_bytes()
                self.assertIn(b"\r\n", original)
                self.assertEqual(Path(normalized).read_bytes(), original.replace(b"\r\n", b"\n"))
                self.assertEqual(settings["input_sha256"][key], WDL.sha256(normalized))
                self.assertEqual(settings["original_input_sha256"][key], WDL.sha256(request[key]))
            validation = Path("hub_overlap_validation.txt").read_text(encoding="utf-8")
            self.assertIn("status\tpassed\n", validation)
            self.assertIn("hub_genes_matched\t1\n", validation)
            self.assertIn("hub_genes_missing\t1\n", validation)
            self.assertIn("hub_duplicates_ignored\t1\n", validation)

    def test_invalid_options_fail_before_runtime_or_outputs(self):
        with work_directory() as root:
            expression, hubs = self.inputs(root)
            with mock.patch.object(WDL, "runtime_fingerprint") as runtime:
                with self.assertRaises(ValueError):
                    WDL.prepare({"exp_file": str(expression), "probe_file": str(hubs),
                                 "parameters": {"bootstrap_num": 2}})
                runtime.assert_not_called()
            for name in ("parameters.json", "expression.exp", "hubs.txt", "hub_overlap_validation.txt"):
                self.assertFalse(Path(name).exists(), name)

    def test_zero_overlap_fails_before_normalized_inputs_are_published(self):
        with work_directory() as root:
            expression, hubs = self.inputs(root)
            hubs.write_bytes(b"not_present\n")
            with mock.patch.object(WDL, "runtime_fingerprint", return_value={}):
                with self.assertRaises(SystemExit):
                    WDL.prepare({"exp_file": str(expression), "probe_file": str(hubs), "parameters": {}})
            for name in ("expression.exp", "hubs.txt", "parameters.json"):
                self.assertFalse(Path(name).exists(), name)

    def test_checked_settings_rejects_schema_runtime_and_input_drift(self):
        with work_directory() as root:
            _, fingerprint = self.prepare(root)
            request = {"settings": str(root / "parameters.json"),
                       "exp_file": str(root / "expression.exp"), "probe_file": str(root / "hubs.txt")}
            with mock.patch.object(WDL, "runtime_fingerprint", return_value=fingerprint):
                self.assertEqual(WDL.checked_settings(request)["schema"], WDL.SCHEMA)
                Path(request["probe_file"]).write_bytes(b"changed\n")
                with self.assertRaisesRegex(ValueError, "Prepared input changed"):
                    WDL.checked_settings(request)
            with mock.patch.object(WDL, "runtime_fingerprint", return_value={"changed": True}):
                with self.assertRaisesRegex(ValueError, "runtime/configuration differs"):
                    WDL.checked_settings(request)
            settings = WDL.read_json(request["settings"])
            settings["schema"] = "foreign"
            Path(request["settings"]).write_text(json.dumps(settings), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "schema"):
                WDL.checked_settings(request)

    def test_infer_records_exact_seed_argv_and_successful_output_hash(self):
        with work_directory() as root:
            _, fingerprint = self.prepare(root)
            request = {"settings": str(root / "parameters.json"), "seed": 1,
                       "exp_file": str(root / "expression.exp"), "probe_file": str(root / "hubs.txt")}

            def native(command, stdout, stderr, check):
                self.assertFalse(check)
                self.assertIsInstance(command, list)
                Path(command[command.index("-o") + 1]).write_bytes(b"> parameters\nA\tB\t0.5\n")
                stdout.write("native success")
                return SimpleNamespace(returncode=0)

            with mock.patch.object(WDL, "runtime_fingerprint", return_value=fingerprint):
                with mock.patch.object(WDL.subprocess, "run", side_effect=native) as launch:
                    WDL.infer(request)
                    launch.assert_called_once()
                provenance = WDL.read_json("run_provenance.json")
                self.assertEqual(provenance["seed"], 1)
                self.assertEqual(provenance["returncode"], 0)
                self.assertEqual(provenance["runtime"], fingerprint)
                self.assertEqual(provenance["adjacency_sha256"], WDL.sha256("TF_run_001.adj"))
                self.assertEqual(provenance["prepared_settings_sha256"], WDL.sha256(request["settings"]))
                self.assertEqual(provenance["argv"], launch.call_args.args[0])
                with mock.patch.object(WDL.subprocess, "run") as repeat:
                    with self.assertRaisesRegex(ValueError, "overwrite"):
                        WDL.infer(request)
                    repeat.assert_not_called()

    def test_infer_failure_records_provenance_and_raises(self):
        with work_directory() as root:
            _, fingerprint = self.prepare(root)
            request = {"settings": str(root / "parameters.json"), "seed": 1,
                       "exp_file": str(root / "expression.exp"), "probe_file": str(root / "hubs.txt")}
            with mock.patch.object(WDL, "runtime_fingerprint", return_value=fingerprint), \
                    mock.patch.object(WDL.subprocess, "run", return_value=SimpleNamespace(returncode=7)), \
                    mock.patch.object(WDL.sys, "stderr"):
                with self.assertRaisesRegex(RuntimeError, "Native inference failed"):
                    WDL.infer(request)
            provenance = WDL.read_json("run_provenance.json")
            self.assertEqual(provenance["returncode"], 7)
            self.assertNotIn("adjacency_sha256", provenance)
            self.assertFalse(Path("TF_run_001.adj").exists())

    def consensus_fixture(self, root, legacy=False):
        options = {"bootstrap_num": 2, **({"p_value_consensus": 0.01} if legacy else {"min_recurrence": 1})}
        _, fingerprint = self.prepare(root, options)
        settings = WDL.read_json("parameters.json")
        source = root / "source_networks"
        source.mkdir()
        (source / "native.log").write_text("must never enter consensus inputs", encoding="utf-8")
        adjacencies, records = [], []
        for seed in (1, 2):
            adjacency = source / WDL.adjacency_name(seed)
            adjacency.write_bytes(b"" if seed == 1 else b"> parameters\n")
            record = source / ("seed_%d.json" % seed)
            WDL.write_json(record, {"seed": seed, "returncode": 0, "runtime": fingerprint,
                "input_sha256": settings["input_sha256"], "adjacency_sha256": WDL.sha256(adjacency),
                "prepared_settings_sha256": WDL.sha256(root / "parameters.json")})
            adjacencies.append(str(adjacency))
            records.append(str(record))
        return {"settings": str(root / "parameters.json"), "exp_file": str(root / "expression.exp"),
                "probe_file": str(root / "hubs.txt"), "adjacencies": adjacencies,
                "run_provenance": records}, fingerprint

    def test_consensus_accepts_empty_networks_and_stages_only_exact_adjacencies(self):
        with work_directory() as root:
            request, fingerprint = self.consensus_fixture(root)
            with mock.patch.object(WDL, "runtime_fingerprint", return_value=fingerprint):
                WDL.consensus(request)
            self.assertEqual(sorted(path.name for path in Path("adjacencies").iterdir()),
                             ["TF_run_001.adj", "TF_run_002.adj"])
            self.assertEqual(Path("results/consensus_network_3col_.txt").read_text(), "source\ttarget\tMI\n")
            self.assertIn("Bootstrap networks: 2\n", Path("results/bootstrap_info_.txt").read_text())
            self.assertTrue(Path("results/consensus_network_ncol_.txt").is_file())

    def test_consensus_rejects_missing_extra_reordered_or_mismatched_evidence_before_staging(self):
        cases = ("missing_adjacency", "extra_adjacency", "missing_record", "reorder_adjacencies",
                 "reorder_records", "wrong_seed", "failed_seed", "runtime", "input_hash", "adjacency_hash", "settings_hash")
        for case in cases:
            with self.subTest(case=case), work_directory() as root:
                request, fingerprint = self.consensus_fixture(root)
                if case == "missing_adjacency":
                    request["adjacencies"].pop()
                elif case == "extra_adjacency":
                    request["adjacencies"].append(request["adjacencies"][0])
                elif case == "missing_record":
                    request["run_provenance"].pop()
                elif case.startswith("reorder_"):
                    key = "adjacencies" if case == "reorder_adjacencies" else "run_provenance"
                    request[key].reverse()
                elif case == "adjacency_hash":
                    Path(request["adjacencies"][0]).write_bytes(b"> altered\n")
                else:
                    path = Path(request["run_provenance"][0])
                    record = WDL.read_json(path)
                    key, value = {"wrong_seed": ("seed", 2), "failed_seed": ("returncode", 1),
                        "runtime": ("runtime", {}), "input_hash": ("input_sha256", {}),
                        "settings_hash": ("prepared_settings_sha256", "different-prepared-settings")}[case]
                    record[key] = value
                    path.write_text(json.dumps(record), encoding="utf-8")
                with mock.patch.object(WDL, "runtime_fingerprint", return_value=fingerprint):
                    with self.assertRaises(ValueError):
                        WDL.consensus(request)
                self.assertFalse(Path("adjacencies").exists())
                self.assertFalse(Path("results").exists())

    def test_consensus_passes_legacy_probability_without_implicit_k(self):
        with work_directory() as root:
            request, fingerprint = self.consensus_fixture(root, legacy=True)
            module = importlib.import_module("SJARACNe.bin.create_consensus_network")
            with mock.patch.object(WDL, "runtime_fingerprint", return_value=fingerprint), \
                    mock.patch.object(module, "create_consensus_network",
                               return_value="results/consensus_network_3col_.txt") as consensus, \
                    mock.patch.object(module, "create_enhanced_consensus_network") as enhance:
                WDL.consensus(request)
            consensus.assert_called_once_with("adjacencies", out_dir="results", p_value=0.01)
            enhance.assert_called_once_with(request["exp_file"], "results/consensus_network_3col_.txt", "results")


if __name__ == "__main__":
    unittest.main()
