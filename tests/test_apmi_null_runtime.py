#!/usr/bin/env python3
"""Native PCHIP integration and public workflow-profile regressions."""

import importlib
import importlib.util
import json
import math
import os
import re
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from native_test_support import find_native_executable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PRIMARY_SHA256 = "8b4104be3f3ffe0a476a4dd7c8f7f6d95f00f1a51fa467c7194edf3806d427f5"


SJARACNE_EXE = find_native_executable()


def python_cutoff(m, p):
    module = importlib.import_module("SJARACNe.apmi_cutoff")
    return module.get_mi_cutoff(m, p, npar=40)


@unittest.skipUnless(SJARACNE_EXE, "sjaracne.exe is not built for this OS")
class TestApmiNullRuntime(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.workdir = Path(self.folder.name)
        self.expression = self.write_expression(100)
        self.hubs = self.workdir / "hubs.txt"
        self.hubs.write_text("H\n", encoding="utf-8")

    def tearDown(self):
        self.folder.cleanup()

    def write_expression(self, m):
        path = self.workdir / f"input_{m}.exp"
        rows = (("H", list(range(m))), ("A", [(i * 7) % m for i in range(m)]),
                ("B", list(reversed(range(m)))))
        lines = ["isoformId\tgeneSymbol\t" + "\t".join(f"s{i}" for i in range(m))]
        lines.extend("\t".join([name, name] + [str(v) for v in values])
                     for name, values in rows)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path

    def run_sjaracne(self, *extra, output="out.adj", expression=None, npar=40):
        command = [SJARACNE_EXE, "-i", str(expression or self.expression),
                   "-s", str(self.hubs), "-e", "1", "-S", "9",
                   "-o", str(self.workdir / output)]
        if npar is not None:
            command += ["-N", str(npar)]
        return subprocess.run(command + list(extra), capture_output=True,
                              text=True, check=False)

    def assert_success(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def assert_failure(self, result, *markers):
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        for marker in markers:
            self.assertIn(marker.lower(), (result.stdout + result.stderr).lower())

    def header(self, output="out.adj"):
        return (self.workdir / output).read_text(encoding="utf-8")

    @staticmethod
    def header_value(header, key):
        match = re.search(r"^>\s+" + re.escape(key) + r"\s+([^\r\n]+)$", header, re.M | re.I)
        if match is None:
            raise AssertionError(f"Missing provenance field {key!r}:\n{header}")
        return match.group(1).strip()

    @staticmethod
    def data_rows(path):
        return [line for line in path.read_text(encoding="utf-8").splitlines()
                if line and not line.startswith(">")]

    def test_pchip_records_full_precision_cutoff_and_pinned_provenance(self):
        result = self.run_sjaracne("-M", "pchip", "-p", "0.001")
        self.assert_success(result)
        header = self.header()
        self.assertEqual(self.header_value(header, "MI threshold method"),
                         "experimental PCHIP empirical cutoff")
        self.assertEqual(self.header_value(header, "AP-MI PCHIP primary SHA256"), PRIMARY_SHA256)
        for key in ("AP-MI PCHIP data SHA256", "AP-MI PCHIP summary SHA256"):
            self.assertRegex(self.header_value(header, key), r"^[a-f0-9]{64}$")
        self.assertEqual(self.header_value(header, "AP-MI actual m"), "100")
        self.assertEqual(self.header_value(header, "AP-MI outside principal p range"), "no")
        self.assertEqual(self.header_value(header, "AP-MI extreme-tail p"), "no")
        self.assertEqual(self.header_value(header, "AP-MI off-grid m"), "no")
        actual = float(self.header_value(header, "MI threshold binary64"))
        self.assertEqual(actual.hex(), python_cutoff(100, 0.001).hex())
        self.assertEqual(float.fromhex(self.header_value(header, "MI threshold hex")), actual)
        self.assertNotIn("generalized-pareto", header.lower())

    def test_actual_subsample_size_not_full_input_size_selects_cutoff(self):
        for specification in ("80", "80%"):
            with self.subTest(specification=specification):
                result = self.run_sjaracne("-M", "pchip", "-p", "0.001", "-u", specification)
                self.assert_success(result)
                header = self.header()
                self.assertEqual(self.header_value(header, "AP-MI actual m"), "80")
                self.assertEqual(self.header_value(header, "AP-MI off-grid m"), "yes")
                value = float(self.header_value(header, "MI threshold binary64"))
                self.assertAlmostEqual(value, 0.14783093032516142, places=15)
                self.assertAlmostEqual(value, python_cutoff(80, 0.001), places=15)

    def test_native_python_parity_at_sensitive_knots_boundary_and_endpoints(self):
        cases = ((50, 0.05), (80, 0.0005), (350, 0.001),
                 (100, 0.0001), (10000, 0.05), (10000, 1e-7))
        for m, p in cases:
            with self.subTest(m=m, p=p):
                arguments = ["-M", "pchip", "-p", format(p, ".17g")]
                result = self.run_sjaracne(*arguments, expression=self.write_expression(m))
                self.assert_success(result)
                header = self.header()
                actual = float(self.header_value(header, "MI threshold binary64"))
                expected = python_cutoff(m, p)
                if m in (50, 100, 10000):
                    self.assertEqual(actual.hex(), expected.hex())
                else:
                    self.assertAlmostEqual(actual, expected, places=15)
                self.assertEqual(float.fromhex(self.header_value(header, "MI threshold hex")), actual)

    def test_rare_tail_is_computable_without_acknowledgement_flag(self):
        for p in ("0.0001", "0.000499999999", "1e-5", "1e-7"):
            with self.subTest(p=p):
                result = self.run_sjaracne("-M", "pchip", "-p", p)
                self.assert_success(result)
                header = self.header()
                self.assertEqual(self.header_value(header, "AP-MI outside principal p range"), "yes")
                self.assertEqual(self.header_value(header, "AP-MI extreme-tail p"),
                                 "yes" if float(p) <= 1e-5 else "no")
                self.assertNotIn("AP-MI rare-tail acknowledged", header)
                self.assertRegex(result.stdout + result.stderr, r"(?i)outside.*principal")

    def test_exact_moderate_boundary_does_not_require_acknowledgement(self):
        self.assert_success(self.run_sjaracne("-M", "pchip", "-p", "0.0005"))

    def test_pchip_fails_closed_outside_m_range(self):
        for m in (49, 10001):
            with self.subTest(m=m):
                result = self.run_sjaracne("-M", "pchip", "-p", "0.001",
                                           expression=self.write_expression(m))
                self.assert_failure(result, "m")
                self.assertFalse((self.workdir / "out.adj").exists())

    def test_pchip_fails_closed_for_nonmatching_npar(self):
        for npar in (20, 41):
            with self.subTest(npar=npar):
                self.assert_failure(self.run_sjaracne("-M", "pchip", "-p", "0.001",
                                                      npar=npar), "40")

    def test_bad_probabilities_are_rejected_without_legacy_fallback(self):
        for p in ("0", "-0.001", "0.05000001", "1e-8", "nan", "inf", "nonsense"):
            with self.subTest(p=p):
                result = self.run_sjaracne("-M", "pchip", "-p", p)
                self.assert_failure(result)
                self.assertNotIn("using the legacy affine", result.stdout.lower())

    def test_unknown_model_and_old_gpd_file_are_rejected_explicitly(self):
        self.assert_failure(self.run_sjaracne("-M", "unknown", "-p", "0.001"), "pchip")
        model = self.workdir / "historical-gpd.model"
        model.write_text("format=sjaracne-apmi-gpd-tail-v1\n", encoding="utf-8")
        result = self.run_sjaracne("-M", str(model), "-p", "0.001")
        self.assert_failure(result, "pchip")
        self.assertRegex(result.stdout + result.stderr, r"(?i)GPD|historical|file|retired")

    def test_active_pchip_rejects_replacement_sampling_and_adjacency_replay(self):
        self.assert_failure(self.run_sjaracne("-M", "pchip", "-p", "0.001", "-r", "1"),
                            "replacement")
        self.assert_failure(self.run_sjaracne("-M", "pchip", "-p", "0.001", "-j",
                                              str(self.workdir / "old.adj")), "adjacency")

    def test_removed_rare_acknowledgement_flag_is_an_unknown_option(self):
        for arguments in (("-R", "-p", "0.001"),
                          ("-M", "pchip", "-R", "-p", "1e-7"),
                          ("-M", "pchip", "-R", "-t", "0")):
            with self.subTest(arguments=arguments):
                self.assert_failure(self.run_sjaracne(*arguments), "unknown", "-R")

    def test_explicit_threshold_including_zero_overrides_model(self):
        for threshold in ("0", "0.14783093032516142"):
            with self.subTest(threshold=threshold):
                explicit = self.run_sjaracne("-t", threshold, "-p", "1e-7", output="explicit.adj")
                overridden = self.run_sjaracne("-M", str(self.workdir / "missing-gpd.model"),
                                               "-p", "1e-7", "-t", threshold,
                                               output="overridden.adj")
                self.assert_success(explicit)
                self.assert_success(overridden)
                self.assertIn("ignored", overridden.stdout.lower())
                self.assertEqual(self.data_rows(self.workdir / "explicit.adj"),
                                 self.data_rows(self.workdir / "overridden.adj"))
                self.assertNotIn("AP-MI PCHIP primary SHA256", self.header("overridden.adj"))

    def test_p_one_and_native_default_bypass_pchip_even_outside_model_domain(self):
        expression = self.write_expression(8)
        default = self.run_sjaracne(expression=expression, npar=None, output="default.adj")
        self.assert_success(default)
        for model in ("pchip", str(self.workdir / "missing-gpd.model")):
            with self.subTest(model=model):
                explicit = self.run_sjaracne("-M", model, "-p", "1",
                                             expression=expression, npar=20,
                                             output="preserve.adj")
                self.assert_success(explicit)
                self.assertEqual(self.data_rows(self.workdir / "default.adj"),
                                 self.data_rows(self.workdir / "preserve.adj"))
                self.assertNotIn("AP-MI PCHIP primary SHA256", self.header("preserve.adj"))

    def test_inactive_model_does_not_reject_legacy_sampling(self):
        for options in (("-p", "1"), ("-p", "0.001", "-t", "0")):
            with self.subTest(options=options):
                result = self.run_sjaracne("-M", "pchip", "-r", "1", *options)
                self.assert_success(result)
                self.assertNotIn("AP-MI PCHIP primary SHA256", self.header())

    def test_absent_model_retains_legacy_affine_formula(self):
        result = self.run_sjaracne("-p", "0.001", "-H", str(PROJECT_ROOT / "SJARACNe" / "config"))
        self.assert_success(result)
        self.assertIn("using the legacy affine threshold calibration", result.stdout)
        match = re.search(r"MI threshold determined for p=.*: ([^\s]+)", result.stdout)
        self.assertIsNotNone(match, result.stdout)
        expected = (1.062 - math.log(0.001)) / (48.7 + 0.634 * 100)
        self.assertAlmostEqual(float(match.group(1)), expected, places=6)
        self.assertNotIn("AP-MI PCHIP primary SHA256", self.header())

    def test_bare_model_option_and_malformed_explicit_threshold_fail_cleanly(self):
        self.assert_failure(self.run_sjaracne("-M"), "requires")
        for args in (("-t",), ("-t", "not-a-number"), ("-t", "nan"), ("-t", "inf")):
            with self.subTest(args=args):
                self.assert_failure(self.run_sjaracne(*args), "finite")


class TestApmiNullWorkflowWiring(unittest.TestCase):
    @staticmethod
    def load_wrapper():
        path = PROJECT_ROOT / "SJARACNe" / "sjaracne.py"
        spec = importlib.util.spec_from_file_location("sjaracne_wrapper_null", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def render_wrapper(self, *extra, hub_name="hubs.txt", network_type="TF"):
        module = self.load_wrapper()
        commands = []
        module.run_shell_command_call = commands.append
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            from test_sjaracne_acceptance import write_multihub_dpi_inputs
            expression, original_hubs = write_multihub_dpi_inputs(root)
            hubs = root / hub_name
            hubs.write_bytes(original_hubs.read_bytes())
            arguments = ["sjaracne", "local", "--engine", "cwl", "-e", str(expression),
                         "-g", str(hubs), "-o", str(root / "output"),
                         "-tmp", str(root / "tmp"), *extra]
            if network_type is not None:
                arguments += ['--network-type', network_type]
            with mock.patch.object(os.sys, "argv", arguments), \
                    mock.patch.object(module, 'runtime_fingerprint', return_value={}), \
                    mock.patch.object(module, 'validate_cwl_outputs', return_value={}), \
                    mock.patch.object(module, 'require_tool', side_effect=lambda name: name):
                module.main()
            rendered = json.loads((root / "output" / "sjaracne_workflow.yml").read_text())
        self.assertEqual(len(commands), 1)
        return rendered

    def test_cwl_model_is_optional_literal_string_not_an_external_file(self):
        workflow = (PROJECT_ROOT / "SJARACNe" / "cwl" / "sjaracne_workflow.cwl").read_text()
        command = (PROJECT_ROOT / "SJARACNe" / "cwl" / "sjaracne.cwl").read_text()
        self.assertIn("apmi_null_model: apmi_null_model", workflow)
        self.assertNotIn("allow_rare_tail", workflow)
        self.assertNotIn("allow_rare_tail", command)
        bootstrap_input = workflow.split("  p_value_bootstrap:", 1)[1].split("  apmi_null_model:", 1)[0]
        self.assertIn("type: double?", bootstrap_input)
        self.assertNotIn("default:", bootstrap_input)
        model_input = command.split("  apmi_null_model:", 1)[1].split("  p_value:", 1)[0]
        self.assertIn("type: string?", model_input)
        self.assertNotIn("type: File", model_input)
        self.assertIn("prefix: -M", model_input)
        self.assertNotIn("prefix: -R", command)

    def test_wrapper_forwards_pchip_without_file_staging(self):
        rendered = self.render_wrapper("-M", "pchip", "-pb", "0.001")
        self.assertEqual(rendered['apmi_null_model'], 'pchip')
        self.assertIsInstance(rendered['apmi_null_model'], str)
        self.assertEqual(rendered['p_value_bootstrap'], 0.001)

    def test_wrapper_accepts_explicit_rare_p_without_extra_flag(self):
        rendered = self.render_wrapper("-M", "pchip", "-pb", "1e-7")
        self.assertEqual(rendered['p_value_bootstrap'], 1e-7)
        self.assertNotIn("allow_rare_tail", rendered)

    def test_wrapper_requires_network_type_before_launch_even_with_explicit_probability(self):
        module = self.load_wrapper()
        commands = []
        module.run_shell_command_call = commands.append
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            from test_sjaracne_acceptance import write_multihub_dpi_inputs
            expression, hubs = write_multihub_dpi_inputs(root)
            arguments = ["sjaracne", "local", "--engine", "cwl", "-e", str(expression),
                         "-g", str(hubs), "-o", str(root / "output"),
                         "-tmp", str(root / "tmp"), "-M", "pchip", '-pb', '0.001']
            with mock.patch.object(os.sys, "argv", arguments), \
                    mock.patch.object(module, 'runtime_fingerprint', return_value={}), \
                    mock.patch.object(module, 'validate_cwl_outputs', return_value={}), \
                    mock.patch.object(module, 'require_tool', side_effect=lambda name: name), \
                    self.assertRaises(SystemExit):
                module.main()
            self.assertFalse((root / 'output').exists())
            self.assertFalse((root / 'tmp').exists())
        self.assertEqual(commands, [])

    def test_wrapper_accepts_explicit_preserve_all_probability(self):
        rendered = self.render_wrapper("-M", "pchip", "-pb", "1")
        self.assertEqual(rendered['p_value_bootstrap'], 1)

    def test_wrapper_does_not_select_probability_from_tf_or_sig_filename(self):
        for hub_name in ("TF_hubs.txt", "SIG_hubs.txt"):
            with self.subTest(hub_name=hub_name):
                current = self.render_wrapper(hub_name=hub_name)
                self.assertEqual(current['p_value_bootstrap'], 0.001)
                explicit = self.render_wrapper("-M", "pchip", "-pb", "0.007", hub_name=hub_name)
                self.assertEqual(explicit['p_value_bootstrap'], 0.007)
        sig = self.render_wrapper(network_type='SIG', hub_name='TF_named_but_explicit_SIG.txt')
        self.assertEqual(sig['p_value_bootstrap'], 0.0005)

    def test_wrapper_rejects_removed_acknowledgement_option(self):
        with self.assertRaises(SystemExit):
            self.render_wrapper("-M", "pchip", "-pb", "1e-7", "--allow-rare-tail")

    def test_wrapper_current_defaults_use_pchip_dynamic_and_unique_sampling(self):
        rendered = self.render_wrapper()
        self.assertEqual(rendered['p_value_bootstrap'], 0.001)
        self.assertEqual(rendered['subsample_spec'], '80%')
        self.assertEqual(rendered['apmi_null_model'], 'pchip')
        self.assertTrue(rendered['dynamic_dpi'])
        self.assertNotIn("allow_rare_tail", rendered)

    def test_wrapper_explicit_legacy_profile_retains_historical_scientific_defaults(self):
        rendered = self.render_wrapper('--profile', 'legacy', network_type=None)
        self.assertEqual(rendered['p_value_bootstrap'], 1e-7)
        self.assertIsNone(rendered.get('apmi_null_model'))
        self.assertFalse(rendered['dynamic_dpi'])
        self.assertEqual(rendered['tolerance'], 0)
        self.assertEqual(rendered['p_value_consensus'], 1e-5)
        self.assertIsNone(rendered.get('min_recurrence'))
        self.assertIsNone(rendered.get('subsample_spec'))
        self.assertEqual(rendered['sample_number'], 1)

    def test_wrapper_rejects_old_model_paths_and_unknown_selectors(self):
        for model in ("old-gpd.model", "unknown"):
            with self.subTest(model=model), self.assertRaises(SystemExit):
                self.render_wrapper("-M", model)


if __name__ == "__main__":
    unittest.main()
