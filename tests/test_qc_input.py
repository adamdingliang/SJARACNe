#!/usr/bin/env python3

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
QC_SCRIPT = PROJECT_ROOT / "SJARACNe" / "bin" / "QC_input.py"


class TestQCInput(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.workdir = Path(self.folder.name)
        self.expression = self.workdir / "small.exp"
        self.expression.write_text(
            "isoformId\tgeneSymbol\ts1\ts2\n"
            "A\tA\t1\t2\n"
            "B\tB\t2\t1\n"
            "C\tC\t3\t4\n",
            encoding="utf-8",
        )

    def tearDown(self):
        self.folder.cleanup()

    def run_qc(self, probe_file):
        report = self.workdir / "hub_overlap_validation.txt"
        result = subprocess.run(
            [
                sys.executable,
                str(QC_SCRIPT),
                "-e",
                str(self.expression),
                "-g",
                str(probe_file),
                "-o",
                str(report),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        return result, report

    @staticmethod
    def read_report(report):
        return dict(
            line.split("\t", 1)
            for line in report.read_text(encoding="utf-8").splitlines()
        )

    def test_partial_overlap_passes_and_reports_normalized_unique_hubs(self):
        probes = self.workdir / "partial_hubs.txt"
        probes.write_bytes(b"\xef\xbb\xbf A \r\n\r\n\tA\t\rMISSING\n")

        result, report = self.run_qc(probes)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(report.is_file())
        self.assertIn("Hub overlap: requested=2, matched=1, missing=1", result.stderr)
        self.assertEqual(
            self.read_report(report),
            {
                "status": "passed",
                "expression_genes": "3",
                "hub_genes_requested": "2",
                "hub_genes_matched": "1",
                "hub_genes_missing": "1",
                "hub_duplicates_ignored": "1",
                "matched_hubs": "A",
                "missing_hubs": "MISSING",
            },
        )

    def test_zero_overlap_fails_before_creating_validation_report(self):
        probes = self.workdir / "unresolved_hubs.txt"
        probes.write_text("NOT_PRESENT\nNOT_PRESENT\n", encoding="utf-8")

        result, report = self.run_qc(probes)

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("requested=1, matched=0, missing=1", result.stderr)
        self.assertIn("Refusing to start bootstrap jobs", result.stderr)
        self.assertFalse(report.exists())

    def test_blank_hub_file_fails_before_creating_validation_report(self):
        probes = self.workdir / "blank_hubs.txt"
        probes.write_text(" \n\t\r\n", encoding="utf-8")

        result, report = self.run_qc(probes)

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("requested=0, matched=0, missing=0", result.stderr)
        self.assertIn("Refusing to start bootstrap jobs", result.stderr)
        self.assertFalse(report.exists())

    def test_bootstrap_has_hard_dependency_on_validation_report(self):
        workflow = (
            PROJECT_ROOT / "SJARACNe" / "cwl" / "dpi_before_consensus.cwl"
        ).read_text(encoding="utf-8")
        qc_tool = (PROJECT_ROOT / "SJARACNe" / "cwl" / "QC_input.cwl").read_text(
            encoding="utf-8"
        )
        bootstrap_tool = (
            PROJECT_ROOT / "SJARACNe" / "cwl" / "sjaracne.cwl"
        ).read_text(encoding="utf-8")

        self.assertIn("out: [validation_report]", workflow)
        self.assertIn(
            "preflight_report: validate_inputs/validation_report", workflow
        )
        self.assertIn(
            "validation_report:\n    type: File\n    outputBinding:", qc_tool
        )
        self.assertIn(
            "preflight_report:\n    type: File?\n    doc:", bootstrap_tool
        )

    def test_active_python_cwl_tools_do_not_depend_on_script_shebangs(self):
        cwl_dir = PROJECT_ROOT / "SJARACNe" / "cwl"
        expected_commands = {
            "QC_input.cwl": "baseCommand: [python3, -m, SJARACNe.bin.QC_input]",
            "ch_line_ending.cwl": (
                "baseCommand: [python3, -m, SJARACNe.bin.ch_line_ending]"
            ),
            "create_consensus_network.cwl": (
                "baseCommand: [python3, -m, SJARACNe.bin.create_consensus_network]"
            ),
        }

        for filename, expected in expected_commands.items():
            with self.subTest(filename=filename):
                self.assertIn(
                    expected,
                    (cwl_dir / filename).read_text(encoding="utf-8"),
                )


if __name__ == "__main__":
    unittest.main()
