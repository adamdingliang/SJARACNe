#!/usr/bin/env python3

import subprocess
import tempfile
import unittest
from pathlib import Path

from native_test_support import find_native_executable


SJARACNE_EXE = find_native_executable()


@unittest.skipUnless(SJARACNE_EXE, "sjaracne.exe is not built")
class TestExpressionAndMIValidation(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.workdir = Path(self.folder.name)
        self.expression = self.workdir / "input.exp"
        self.adjacency = self.workdir / "input.adj"
        self.output = self.workdir / "network.adj"

    def tearDown(self):
        self.folder.cleanup()

    def run_sjaracne(self, expression_body, *extra_args):
        self.expression.write_text(expression_body, encoding="utf-8")
        if self.output.exists():
            self.output.unlink()

        result = subprocess.run(
            [
                SJARACNE_EXE,
                "-i",
                str(self.expression),
                "-p",
                "1",
                "-e",
                "1",
                "-o",
                str(self.output),
                *extra_args,
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        return result

    def test_first_expression_row_observation_count_is_validated(self):
        cases = {
            "expression_only": (
                "isoformId\tgeneSymbol\ts1\ts2\n"
                "A\tA\t1 2\t3\n"
                "B\tB\t4\t5\n"
            ),
            "expression_pvalue_pairs": (
                "isoformId\tgeneSymbol\ts1\ts2\n"
                "A\tA\t1 0.1 2\t0.2\t3\t0.3\n"
                "B\tB\t4\t0.1\t5\t0.2\n"
            ),
        }

        for name, expression_body in cases.items():
            with self.subTest(format=name):
                result = self.run_sjaracne(expression_body)

                self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                self.assertIn(
                    "Incorrect expression dimensions at line 2: "
                    "expected 2 observations, found 3",
                    result.stderr,
                )
                self.assertFalse(self.output.exists())

    def test_at_least_two_observations_are_required(self):
        result = self.run_sjaracne(
            "isoformId\tgeneSymbol\ts1\n"
            "A\tA\t1\n"
            "B\tB\t2\n"
        )

        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn(
            "must contain at least 2 observation columns; found 1",
            result.stderr,
        )
        self.assertFalse(self.output.exists())

    def test_conditional_selection_must_retain_two_observations(self):
        result = self.run_sjaracne(
            "isoformId\tgeneSymbol\ts1\ts2\ts3\ts4\ts5\ts6\n"
            "A\tA\t1\t2\t3\t4\t5\t6\n"
            "B\tB\t6\t5\t4\t3\t2\t1\n",
            "-c",
            "+_A",
            "0.2",
        )

        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn(
            "At least 2 observations are required for MI calculation; "
            "found 1 after conditional selection",
            result.stderr,
        )
        self.assertFalse(self.output.exists())

    def test_adjacency_input_rejects_nonfinite_mi(self):
        expression_body = (
            "isoformId\tgeneSymbol\ts1\ts2\ts3\n"
            "A\tA\t1\t2\t3\n"
            "B\tB\t3\t2\t1\n"
        )

        for value in ("nan", "inf", "-inf"):
            with self.subTest(value=value):
                self.adjacency.write_text(
                    "A\tB\t{}\n".format(value), encoding="utf-8"
                )
                result = self.run_sjaracne(
                    expression_body, "-j", str(self.adjacency)
                )

                self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                self.assertIn("Non-finite MI value", result.stderr)
                self.assertIn("A -> B", result.stderr)
                self.assertFalse(self.output.exists())

    def test_adjacency_input_rejects_malformed_mi(self):
        expression_body = (
            "isoformId\tgeneSymbol\ts1\ts2\ts3\n"
            "A\tA\t1\t2\t3\n"
            "B\tB\t3\t2\t1\n"
        )

        for value in ("not-a-number", "   "):
            with self.subTest(value=repr(value)):
                self.adjacency.write_text(
                    "A\tB\t{}\n".format(value), encoding="utf-8"
                )
                result = self.run_sjaracne(
                    expression_body, "-j", str(self.adjacency)
                )

                self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                self.assertIn("Invalid MI value", result.stderr)
                self.assertIn("A -> B", result.stderr)
                self.assertFalse(self.output.exists())

    def test_nonfinite_noise_corrected_mi_is_rejected(self):
        result = self.run_sjaracne(
            "isoformId\tgeneSymbol\ts1\ts2\n"
            "A\tA\t0\t2\n"
            "B\tB\t0\t2\n",
            "-n",
            "10",
        )

        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("Non-finite noise-corrected MI calculated", result.stderr)
        self.assertIn("A -> B", result.stderr)
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()
