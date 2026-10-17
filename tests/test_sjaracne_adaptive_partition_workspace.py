#!/usr/bin/env python3

import os
import shutil
import subprocess
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


@unittest.skipUnless(SJARACNE_EXE, "SJARACNe executable is not available")
class TestAdaptivePartitionWorkspace(unittest.TestCase):
    @staticmethod
    def read_network(path):
        network = {}
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            if not line or line.startswith(">"):
                continue
            fields = line.split("\t")
            source = fields[0]
            if len(fields[1:]) % 2:
                raise AssertionError("adjacency row has an incomplete edge")
            for index in range(1, len(fields), 2):
                network[(source, fields[index])] = fields[index + 1]
        return network

    def test_deep_shallow_deep_pairs_do_not_leak_workspace_state(self):
        with tempfile.TemporaryDirectory() as temporary:
            workdir = Path(temporary)
            expression = workdir / "deep_shallow_deep.exp"
            hubs = workdir / "hubs.txt"
            output = workdir / "network.adj"

            expression.write_text(
                "isoformId\tgeneSymbol\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\n"
                "H\tH\t1\t2\t3\t4\t5\t6\t7\t8\n"
                "DEEP_A1\tDEEP_A1\t1\t2\t3\t4\t5\t6\t7\t8\n"
                "SHALLOW_B\tSHALLOW_B\t1\t5\t3\t7\t2\t6\t4\t8\n"
                "DEEP_A2\tDEEP_A2\t1\t2\t3\t4\t5\t6\t7\t8\n",
                encoding="utf-8",
            )
            hubs.write_text("H\n", encoding="utf-8")

            result = subprocess.run(
                [
                    SJARACNE_EXE,
                    "-i",
                    str(expression),
                    "-s",
                    str(hubs),
                    "-S",
                    "1",
                    "-t",
                    "0",
                    "-e",
                    "1",
                    "-N",
                    "20",
                    "-o",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            body = [
                line
                for line in output.read_text(encoding="utf-8").splitlines()
                if line and not line.startswith(">")
            ]
            self.assertEqual(
                body,
                ["H\tDEEP_A1\t0.693147\tDEEP_A2\t0.693147"],
            )
            self.assertNotIn("SHALLOW_B", body[0])

    def test_reused_workspace_is_independent_of_hub_order(self):
        with tempfile.TemporaryDirectory() as temporary:
            workdir = Path(temporary)
            forward_hubs = FIXTURES / "tied_hubs.txt"
            reverse_hubs = workdir / "reversed_hubs.txt"
            reverse_hubs.write_text(
                "\n".join(
                    reversed(forward_hubs.read_text(encoding="utf-8").splitlines())
                )
                + "\n",
                encoding="utf-8",
            )

            networks = []
            for name, hubs in (("forward", forward_hubs), ("reverse", reverse_hubs)):
                output = workdir / "{}.adj".format(name)
                result = subprocess.run(
                    [
                        SJARACNE_EXE,
                        "-i",
                        str(FIXTURES / "tied_counts.exp"),
                        "-s",
                        str(hubs),
                        "-S",
                        "17",
                        "-r",
                        "1",
                        "-t",
                        "0",
                        "-e",
                        "1",
                        "-o",
                        str(output),
                    ],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                networks.append(self.read_network(output))

            self.assertEqual(networks[0], networks[1])

    def test_mixed_leaf_and_recursive_quadrants_match_legacy_mi(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "mixed_quadrants.adj"
            result = subprocess.run(
                [
                    SJARACNE_EXE,
                    "-i",
                    str(FIXTURES / "mixed_quadrant_counts.exp"),
                    "-s",
                    str(FIXTURES / "mixed_quadrant_hub.txt"),
                    "-S",
                    "1",
                    "-t",
                    "0",
                    "-e",
                    "1",
                    "-N",
                    "20",
                    "-o",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            body = [
                line
                for line in output.read_text(encoding="utf-8").splitlines()
                if line and not line.startswith(">")
            ]
            self.assertEqual(body, ["H\tMIXED_Q\t0.056633"])


if __name__ == "__main__":
    unittest.main()
