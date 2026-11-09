#!/usr/bin/env python3
"""Native dynamic-DPI regression on tiny synthetic expression matrices only."""

import math
from pathlib import Path
import random
import re
import subprocess
import tempfile
import unittest

from native_test_support import find_native_executable


SJARACNE_EXE = find_native_executable()


@unittest.skipUnless(SJARACNE_EXE, "sjaracne.exe is not built")
class TestDynamicDpi(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.workdir = Path(self.folder.name)
        self.expression = self.workdir / "tiny.exp"
        generator = random.Random(291)
        genes = []
        for name, weight in (("H0", 100), ("H1", 80), ("H2", 60),
                             ("T0", 90), ("T1", 70), ("N", 0)):
            noise = list(range(100))
            generator.shuffle(noise)
            values = [weight * i + (100 - weight) * noise[i] + i / 1000.0
                      for i in range(100)]
            genes.append((name, values))
        genes.append(("ZERO", [0] * 100))
        lines = ["isoformId\tgeneSymbol\t" + "\t".join("s%d" % i for i in range(100))]
        lines.extend("\t".join([name, "symbol_" + name] + [str(v) for v in values])
                     for name, values in genes)
        self.expression.write_text("\n".join(lines) + "\n", encoding="ascii")
        self.hubs = self.workdir / "hubs.txt"
        self.hubs.write_text("H0\nH1\nH2\nT0\n", encoding="ascii")
        self.counter = 0

    def run_native(self, *extra, selected=True, hubs=None, annotation=None):
        self.counter += 1
        output = self.workdir / ("network_%d.adj" % self.counter)
        source = self.hubs if hubs is None else hubs
        command = [SJARACNE_EXE, "-i", str(self.expression), "-N", "40",
                   "-t", "0", "-S", "9", "-o", str(output)]
        if selected:
            command += ["-s", str(source)]
        command += ["-l", str(source if annotation is None else annotation)]
        result = subprocess.run(command + list(extra), capture_output=True,
                                text=True, check=False)
        return result, output

    def assert_success(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def assert_failure(self, result, output, *markers):
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertFalse(output.exists(), result.stdout + result.stderr)
        for marker in markers:
            self.assertIn(marker.lower(), (result.stdout + result.stderr).lower())

    @staticmethod
    def header_value(path, key):
        text = path.read_text(encoding="utf-8")
        match = re.search(r"^>\s+" + re.escape(key) + r"\s+([^\r\n]+)$", text, re.M)
        if match is None:
            raise AssertionError("Missing header %r:\n%s" % (key, text))
        return match.group(1).strip()

    @staticmethod
    def data_rows(path):
        return [line for line in path.read_text(encoding="utf-8").splitlines()
                if line and not line.startswith(">")]

    def epsilon(self, path):
        value = float(self.header_value(path, "DPI epsilon binary64"))
        self.assertEqual(float.fromhex(self.header_value(path, "DPI epsilon hex")), value)
        return value

    def test_default_opt_in_records_parameters_and_actual_source_rows(self):
        result, output = self.run_native("--dynamic-dpi")
        self.assert_success(result)
        self.assertEqual(self.header_value(output, "DPI mode"), "experimental-dynamic")
        self.assertEqual(self.header_value(output, "DPI source rows H"), "4")
        self.assertEqual(self.header_value(output, "DPI H_ref"), "479")
        self.assertEqual(float(self.header_value(output, "DPI slope")), 0.05)
        self.assertEqual(float(self.header_value(output, "DPI epsilon max")), 0.30)
        self.assertEqual(self.epsilon(output), 0.0)
        self.assertIn("[DPI] mode=experimental-dynamic H=4", result.stdout)
        self.assertIn("[NETWORK] Applying DPI", result.stdout)

    def test_formula_anchor_below_doubling_and_cap(self):
        cases = ((8, .05, .30, 0.0), (4, .05, .30, 0.0),
                 (2, .05, .30, .05), (1, .05, .30, .10),
                 (1, .25, .30, .30), (2, .05, 0.0, 0.0))
        for reference, slope, cap, expected in cases:
            with self.subTest(reference=reference, slope=slope, cap=cap):
                result, output = self.run_native(
                    "--dynamic-dpi", "--dpi-href", str(reference),
                    "--dpi-slope", str(slope), "--dpi-max", str(cap))
                self.assert_success(result)
                self.assertEqual(self.epsilon(output), expected)

    def test_non_power_of_two_formula_and_equals_syntax(self):
        result, output = self.run_native("--dynamic-dpi", "--dpi-href=3",
                                         "--dpi-slope=0.07", "--dpi-max=0.2")
        self.assert_success(result)
        expected = min(.2, .07 * max(0.0, math.log2(4 / 3)))
        self.assertAlmostEqual(self.epsilon(output), expected, places=15)

    def test_dynamic_graph_matches_exact_corresponding_fixed_epsilon(self):
        for arguments in (("--dpi-href", "2"),
                          ("--dpi-href", "3", "--dpi-slope", "0.07"),
                          ("--dpi-href", "1", "--dpi-slope", "0.25")):
            with self.subTest(arguments=arguments):
                dynamic, dynamic_path = self.run_native("--dynamic-dpi", *arguments)
                self.assert_success(dynamic)
                fixed, fixed_path = self.run_native("-e", format(self.epsilon(dynamic_path), ".17g"))
                self.assert_success(fixed)
                self.assertEqual(self.data_rows(dynamic_path), self.data_rows(fixed_path))

    def test_opt_out_retains_native_no_dpi_default(self):
        default, default_path = self.run_native()
        fixed, fixed_path = self.run_native("-e", "1")
        self.assert_success(default)
        self.assert_success(fixed)
        self.assertEqual(self.header_value(default_path, "DPI tolerance"), "1")
        self.assertEqual(self.data_rows(default_path), self.data_rows(fixed_path))
        self.assertNotIn("experimental-dynamic", default.stdout)
        self.assertNotIn("DPI source rows H", default_path.read_text(encoding="utf-8"))

    def test_explicit_fixed_epsilon_remains_supported_without_opt_in(self):
        result, output = self.run_native("-e", "0.1")
        self.assert_success(result)
        self.assertEqual(self.header_value(output, "DPI tolerance"), "0.1")
        self.assertNotIn("DPI H_ref", output.read_text(encoding="utf-8"))

    def test_fixed_dynamic_conflict_is_order_independent(self):
        for arguments in (("--dynamic-dpi", "-e", "0"),
                          ("-e", "0", "--dynamic-dpi"),
                          ("--dynamic-dpi", "-e", "1")):
            with self.subTest(arguments=arguments):
                result, output = self.run_native(*arguments)
                self.assert_failure(result, output, "-e", "dynamic")

    def test_tuning_without_opt_in_is_rejected(self):
        for flag, value in (("--dpi-href", "2"), ("--dpi-slope", "0.1"),
                            ("--dpi-max", "0.2")):
            with self.subTest(flag=flag):
                result, output = self.run_native(flag, value)
                self.assert_failure(result, output, "dynamic")

    def test_invalid_reference_counts_are_rejected(self):
        for value in ("0", "-1", "1.5", "nan", "inf", "abc", "2147483648"):
            with self.subTest(value=value):
                result, output = self.run_native("--dynamic-dpi", "--dpi-href", value)
                self.assert_failure(result, output, "dpi-href")

    def test_invalid_slopes_are_rejected(self):
        for value in ("0", "-0.01", "nan", "inf", "-inf", "1e309", "abc"):
            with self.subTest(value=value):
                result, output = self.run_native("--dynamic-dpi", "--dpi-slope", value)
                self.assert_failure(result, output, "dpi-slope")

    def test_invalid_caps_are_rejected(self):
        for value in ("-0.01", "1", "1.1", "nan", "inf", "-inf", "abc"):
            with self.subTest(value=value):
                result, output = self.run_native("--dynamic-dpi", "--dpi-max", value)
                self.assert_failure(result, output, "dpi-max")

    def test_missing_values_and_boolean_value_are_rejected(self):
        for flag in ("--dpi-href", "--dpi-slope", "--dpi-max"):
            with self.subTest(flag=flag):
                result, output = self.run_native("--dynamic-dpi", flag)
                self.assert_failure(result, output, flag)
        result, output = self.run_native("--dynamic-dpi=true")
        self.assert_failure(result, output, "dynamic-dpi")

    def test_positional_boolean_and_hidden_fixed_options_are_rejected(self):
        cases = [
            ("--dynamic-dpi", "false"),
            ("--dynamic-dpi", "--", "-e", "0"),
            ("stray", "--dynamic-dpi"),
            ("stray", "--dynamic-dpi=true"),
            ("--", "--dynamic-dpi"),
            ("--", "stray", "--dynamic-dpi"),
        ]
        for flag, value in (("--dpi-href", "2"), ("--dpi-slope", "0.1"),
                            ("--dpi-max", "0.2")):
            cases.extend((("stray", flag, value), ("stray", flag + "=" + value),
                          ("--", flag, value), ("--", "stray", flag + "=" + value)))
        for arguments in cases:
            with self.subTest(arguments=arguments):
                result, output = self.run_native(*arguments)
                self.assert_failure(result, output, "positional")

    def test_long_option_spelling_consumed_as_filename_is_not_a_mode_request(self):
        # Scan only arguments left after parsing, not option-like filenames
        # already consumed as values by the historical short-option parser.
        filename = self.workdir / "--dynamic-dpi"
        filename.write_bytes(self.expression.read_bytes())
        output = self.workdir / "literal_filename.adj"
        result = subprocess.run([
            SJARACNE_EXE, "-i", filename.name, "-s", str(self.hubs),
            "-e", "0", "-N", "40", "-t", "0", "-o", str(output),
        ], cwd=str(self.workdir), capture_output=True, text=True, check=False)
        self.assert_success(result)
        self.assertTrue(output.is_file())
        self.assertNotIn("experimental-dynamic", result.stdout)
        self.assertNotIn("DPI H_ref", output.read_text(encoding="utf-8"))

    def test_requires_selected_source_list_not_all_gene_or_single_hub_mode(self):
        result, output = self.run_native("--dynamic-dpi", selected=False)
        self.assert_failure(result, output, "-s")
        result, output = self.run_native("--dynamic-dpi", "-h", "H0", selected=False)
        self.assert_failure(result, output, "dynamic")

    def test_replay_is_rejected(self):
        baseline, adjacency = self.run_native("-e", "1")
        self.assert_success(baseline)
        result, output = self.run_native("--dynamic-dpi", "-j", str(adjacency))
        self.assert_failure(result, output, "-j", "dynamic")

    def test_duplicate_and_unresolved_ids_do_not_inflate_h(self):
        sources = self.workdir / "messy_hubs.txt"
        sources.write_bytes(b"\xef\xbb\xbf H0\r\nH1\n H0 \r\nMISSING\n\n")
        result, output = self.run_native("--dynamic-dpi", "--dpi-href", "1", hubs=sources)
        self.assert_success(result)
        self.assertEqual(self.header_value(output, "DPI source rows H"), "2")
        self.assertEqual(self.epsilon(output), .05)
        self.assertIn("matched: 2, missing: 1", result.stdout)

    def test_annotation_count_is_not_h(self):
        annotation = self.workdir / "all_annotation.txt"
        annotation.write_text("H0\nH1\nH2\nT0\nT1\nN\nZERO\n", encoding="ascii")
        result, output = self.run_native("--dynamic-dpi", "--dpi-href", "2", annotation=annotation)
        self.assert_success(result)
        self.assertEqual(self.header_value(output, "DPI source rows H"), "4")
        self.assertEqual(self.epsilon(output), .05)

    def test_empty_rows_after_thresholding_still_count_as_constructed_sources(self):
        result, output = self.run_native("--dynamic-dpi", "--dpi-href", "2", "-t", "1000000")
        self.assert_success(result)
        self.assertEqual(self.header_value(output, "DPI source rows H"), "4")
        self.assertEqual(self.epsilon(output), .05)
        self.assertIn("pre_edges=0", result.stdout)

    def test_inactive_explicit_hub_is_still_a_constructed_row(self):
        sources = self.workdir / "inactive_hub.txt"
        sources.write_text("H0\nH1\nZERO\n", encoding="ascii")
        result, output = self.run_native("--dynamic-dpi", hubs=sources)
        self.assert_success(result)
        self.assertEqual(self.header_value(output, "DPI source rows H"), "3")

    def test_conditional_mode_is_rejected_before_constructing_sources(self):
        # Conditional reconstruction is outside the first dynamic-DPI release.
        # Exercise both first and last requested rows to guard the legacy
        # unconstructed-control-row boundary without changing fixed mode.
        for control in ("+_H0", "+_T0"):
            with self.subTest(control=control):
                result, output = self.run_native("--dynamic-dpi", "-c", control, "0.8")
                self.assert_failure(result, output, "-c", "dynamic")

    def test_empty_and_all_unresolved_source_lists_fail_before_dpi(self):
        sources = self.workdir / "empty_or_missing_hubs.txt"
        for contents in ("", "MISSING\n"):
            with self.subTest(contents=contents):
                sources.write_text(contents, encoding="ascii")
                result, output = self.run_native("--dynamic-dpi", hubs=sources)
                self.assert_failure(result, output, "zero probe IDs matched")

    def test_fewer_than_two_constructed_rows_fail_clearly(self):
        sources = self.workdir / "one_hub.txt"
        sources.write_text("H0\nMISSING\n", encoding="ascii")
        result, output = self.run_native("--dynamic-dpi", hubs=sources)
        self.assert_failure(result, output, "two", "source")

    def test_witness_diagnostic_uses_resolved_dynamic_epsilon(self):
        witness = self.workdir / "witnesses.tsv"
        result, output = self.run_native("--dynamic-dpi", "--dpi-href", "2", "-W", str(witness))
        self.assert_success(result)
        self.assertTrue(witness.is_file())
        header = witness.read_text(encoding="utf-8")
        recorded = re.search(r"^# dpi_epsilon\t([^\r\n]+)$", header, re.M)
        self.assertIsNotNone(recorded)
        self.assertEqual(float(recorded.group(1)), self.epsilon(output))
        self.assertIn("# source_count\t4", header)


if __name__ == "__main__":
    unittest.main()
