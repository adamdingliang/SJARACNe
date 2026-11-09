#!/usr/bin/env python3

import subprocess
import tempfile
import unittest
from pathlib import Path

from native_test_support import find_native_executable


SJARACNE_EXE = find_native_executable()


@unittest.skipUnless(SJARACNE_EXE, "sjaracne.exe is not built")
class TestSubnetworkResolution(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.workdir = Path(self.folder.name)
        self.expression = self.workdir / "small.exp"
        self.expression.write_text(
            "isoformId\tgeneSymbol\ts1\ts2\ts3\ts4\ts5\ts6\n"
            "A\tA\t1\t2\t3\t4\t5\t6\n"
            "B\tB\t6\t5\t4\t3\t2\t1\n"
            "C\tC\t1\t3\t2\t6\t4\t5\n",
            encoding="utf-8",
        )

    def tearDown(self):
        self.folder.cleanup()

    def run_sjaracne(self, *extra_args):
        output = self.workdir / "network.adj"
        command = [
            SJARACNE_EXE,
            "-i",
            str(self.expression),
            "-p",
            "1",
            "-e",
            "1",
            "-o",
            str(output),
            *extra_args,
        ]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        return result, output

    def test_omitted_subnetwork_retains_all_gene_mode(self):
        result, output = self.run_sjaracne()

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Gene: 3", result.stdout)
        self.assertNotIn("[SUBNETWORK]", result.stdout)
        self.assertTrue(output.is_file())

    def test_partial_subnetwork_match_continues_with_summary(self):
        hubs = self.workdir / "partial_hubs.txt"
        hubs.write_text("A\nNOT_IN_EXPRESSION\n", encoding="utf-8")

        result, output = self.run_sjaracne("-s", str(hubs))

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(
            "[SUBNETWORK] Requested: 2, matched: 1, missing: 1", result.stdout
        )
        self.assertIn("Gene: 1", result.stdout)
        self.assertTrue(output.is_file())

    def test_single_computed_hub_skips_impossible_dpi(self):
        hubs = self.workdir / "single_hub.txt"
        hubs.write_text("A\n", encoding="utf-8")

        result, output = self.run_sjaracne("-s", str(hubs), "-e", "0")

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(
            "[NETWORK] Skipping DPI: fewer than two source rows are available",
            result.stdout,
        )
        self.assertNotIn("[NETWORK] Applying DPI", result.stdout)
        self.assertTrue(output.is_file())

    def test_two_computed_hubs_still_apply_dpi(self):
        hubs = self.workdir / "two_hubs.txt"
        hubs.write_text("A\nB\n", encoding="utf-8")

        result, output = self.run_sjaracne("-s", str(hubs), "-e", "0")

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("[NETWORK] Applying DPI", result.stdout)
        self.assertNotIn("[NETWORK] Skipping DPI", result.stdout)
        self.assertTrue(output.is_file())

    def test_subnetwork_list_normalizes_bom_whitespace_blank_lines_and_endings(self):
        hubs = self.workdir / "normalized_hubs.txt"
        hubs.write_bytes(b"\xef\xbb\xbf  A  \r\n\r\n\tB\t\r C \n")

        result, output = self.run_sjaracne("-s", str(hubs))

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(
            "[SUBNETWORK] Requested: 3, matched: 3, missing: 0", result.stdout
        )
        self.assertIn("Gene: 3", result.stdout)
        self.assertNotIn("Cannot find probe", result.stdout)
        self.assertTrue(output.is_file())

    def test_tf_annotation_list_uses_the_same_normalization(self):
        hubs = self.workdir / "normalized_tf_hubs.txt"
        hubs.write_bytes(b"\xef\xbb\xbf\tA \r\n\r\n B\t\r")

        result, output = self.run_sjaracne("-s", str(hubs), "-l", str(hubs))

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(f"[PARA] TF annotation list: {hubs} (2)", result.stdout)
        self.assertIn(
            "[SUBNETWORK] Requested: 2, matched: 2, missing: 0", result.stdout
        )
        self.assertNotIn("Cannot find probe", result.stdout)
        self.assertTrue(output.is_file())

    def test_duplicate_subnetwork_hub_matches_single_hub_output(self):
        hubs = self.workdir / "equivalent_hubs.txt"
        hubs.write_text("A\n", encoding="utf-8")

        single_result, output = self.run_sjaracne("-s", str(hubs))

        self.assertEqual(
            single_result.returncode,
            0,
            single_result.stdout + single_result.stderr,
        )
        self.assertTrue(output.is_file())
        single_hub_output = output.read_bytes()

        hubs.write_bytes(b"A\n A \r\n\tA\t\r")

        result, output = self.run_sjaracne("-s", str(hubs))

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(f"Subset of probes to reconstruct: {hubs} (1)", result.stdout)
        self.assertIn("Duplicate subnetwork probes ignored: 2", result.stdout)
        self.assertIn(
            "[SUBNETWORK] Requested: 1, matched: 1, missing: 0", result.stdout
        )
        self.assertIn("Gene: 1", result.stdout)
        self.assertEqual(output.read_bytes(), single_hub_output)

        source_rows = [
            line
            for line in output.read_text(encoding="utf-8").splitlines()
            if line.split("\t", 1)[0] == "A"
        ]
        self.assertEqual(len(source_rows), 1)

    def test_duplicate_tf_annotation_hubs_are_deduplicated(self):
        subnet = self.workdir / "single_hub.txt"
        subnet.write_text("A\n", encoding="utf-8")
        tf_hubs = self.workdir / "duplicate_tf_hubs.txt"
        tf_hubs.write_bytes(b"A\r\n A \nB\r\n\tB\t\r")

        result, output = self.run_sjaracne(
            "-s", str(subnet), "-l", str(tf_hubs)
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(f"[PARA] TF annotation list: {tf_hubs} (2)", result.stdout)
        self.assertIn("Duplicate TF annotation probes ignored: 2", result.stdout)
        self.assertNotIn("Cannot find probe", result.stdout)
        self.assertTrue(output.is_file())

    def test_unresolved_subnetwork_fails_instead_of_running_all_genes(self):
        hubs = self.workdir / "missing_hubs.txt"
        hubs.write_text("NOT_IN_EXPRESSION\n", encoding="utf-8")

        result, output = self.run_sjaracne("-s", str(hubs))

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("zero probe IDs matched", result.stderr)
        self.assertIn("(requested: 1)", result.stderr)
        self.assertIn("Refusing to construct an all-gene network", result.stderr)
        self.assertNotIn("Gene: 3", result.stdout)
        self.assertFalse(output.exists())

    def test_empty_subnetwork_file_fails_instead_of_running_all_genes(self):
        hubs = self.workdir / "empty_hubs.txt"
        hubs.write_text("", encoding="utf-8")

        result, output = self.run_sjaracne("-s", str(hubs))

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("zero probe IDs matched", result.stderr)
        self.assertIn("(requested: 0)", result.stderr)
        self.assertNotIn("Gene: 3", result.stdout)
        self.assertFalse(output.exists())

    def test_subnetwork_option_without_file_is_rejected(self):
        result, output = self.run_sjaracne("-s")

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Option '-s' requires a subnetwork file", result.stderr)
        self.assertNotIn("Gene: 3", result.stdout)
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
