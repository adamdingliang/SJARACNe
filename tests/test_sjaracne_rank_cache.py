#!/usr/bin/env python3

import csv
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = PROJECT_ROOT / "tests" / "fixtures" / "rank_cache"


def find_sjaracne_executable():
    configured = os.environ.get("SJARACNE_TEST_EXE")
    candidates = [
        configured,
        PROJECT_ROOT / "SJARACNe" / "bin" / "sjaracne.exe",
        shutil.which("sjaracne.exe"),
    ]
    for candidate in candidates:
        if not candidate or not Path(candidate).is_file():
            continue
        if os.name == "nt":
            with Path(candidate).open("rb") as handle:
                if handle.read(2) != b"MZ":
                    continue
        return str(candidate)
    return None


SJARACNE_EXE = find_sjaracne_executable()
REFERENCE_RUNTIME_SUPPORTED = (
    sys.platform.startswith("linux") and platform.libc_ver()[0] == "glibc"
)


@unittest.skipUnless(
    SJARACNE_EXE and REFERENCE_RUNTIME_SUPPORTED,
    "fixed bootstrap references require a Linux glibc build",
)
class TestRankBenchmarkCorrectnessFixture(unittest.TestCase):
    def assert_fixture_matches_reference(self, reference_name, *extra_args):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "network.adj"
            result = subprocess.run(
                [
                    SJARACNE_EXE,
                    "-i",
                    str(FIXTURES / "tied_counts.exp"),
                    "-s",
                    str(FIXTURES / "tied_hubs.txt"),
                    "-S",
                    "17",
                    "-t",
                    "0",
                    "-e",
                    "1",
                    *extra_args,
                    "-o",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            actual = {}
            for line in output.read_text(encoding="utf-8").splitlines():
                if not line or line.startswith(">"):
                    continue
                fields = line.split("\t")
                source = fields[0]
                self.assertEqual(len(fields[1:]) % 2, 0)
                for index in range(1, len(fields), 2):
                    actual[(source, fields[index])] = float(fields[index + 1])

            expected = {}
            with (FIXTURES / reference_name).open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    expected[(row["source"], row["target"])] = float(row["mi"])

            self.assertEqual(set(actual), set(expected))
            for edge, expected_mi in expected.items():
                self.assertAlmostEqual(actual[edge], expected_mi, places=6, msg=edge)

    def test_tied_bootstrap_seed_17_matches_reference_mi(self):
        self.assert_fixture_matches_reference(
            "tied_seed17_reference.tsv", "-r", "1"
        )

    def test_conditional_bootstrap_ranks_match_reference_mi(self):
        self.assert_fixture_matches_reference(
            "tied_conditional_bootstrap_seed17_reference.tsv",
            "-c",
            "-_FOXP1",
            "0.75",
            "-r",
            "1",
        )


if __name__ == "__main__":
    unittest.main()
