#!/usr/bin/env python3

import math
import os
import shutil
import struct
import subprocess
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def find_generator():
    configured = os.environ.get("SJARACNE_APMI_NULL_EXE")
    candidates = [
        configured,
        PROJECT_ROOT / "SJARACNe" / "bin" / "apmi_null_generator.exe",
        shutil.which("apmi_null_generator.exe"),
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


GENERATOR = find_generator()


@unittest.skipUnless(GENERATOR, "apmi_null_generator.exe is not built")
class TestApmiNullGenerator(unittest.TestCase):
    COMMON = ("--m", "80", "--draws", "8", "--npar", "40", "--seed", "17")

    @staticmethod
    def run_generator(*arguments):
        return subprocess.run(
            [GENERATOR, *arguments],
            capture_output=True,
            check=False,
        )

    @staticmethod
    def parse_tsv(payload):
        lines = payload.decode("ascii").splitlines()
        rows = [line for line in lines if line and not line.startswith("#")]
        if not rows or rows[0] != "draw\tmi":
            raise AssertionError("missing generator TSV header")
        parsed = []
        for expected_draw, row in enumerate(rows[1:]):
            draw, mi = row.split("\t")
            if int(draw) != expected_draw:
                raise AssertionError("draw indices are not consecutive")
            parsed.append(float(mi))
        return parsed

    def test_seeded_tsv_stream_is_reproducible_and_matches_golden_values(self):
        first = self.run_generator(*self.COMMON)
        second = self.run_generator(*self.COMMON)
        self.assertEqual(first.returncode, 0, first.stderr.decode(errors="replace"))
        self.assertEqual(second.returncode, 0, second.stderr.decode(errors="replace"))
        self.assertEqual(first.stdout, second.stdout)
        self.assertIn(b"# kernel_schema=sjaracne-apmi-v1\n", first.stdout)
        self.assertIn(b"# rng=mt19937-rejection-fisher-yates-v1\n", first.stdout)
        expected = [
            0.0,
            0.0,
            0.0050083668463560826,
            0.0050083668463560826,
            0.020135513550687989,
            0.0050083668463560826,
            0.0012505213548648086,
            0.0,
        ]
        for observed, reference in zip(self.parse_tsv(first.stdout), expected):
            self.assertAlmostEqual(observed, reference, delta=1e-14)

    def test_binary_stream_bitwise_matches_tsv_and_has_metadata(self):
        tsv = self.run_generator(*self.COMMON)
        self.assertEqual(tsv.returncode, 0, tsv.stderr.decode(errors="replace"))
        expected = self.parse_tsv(tsv.stdout)

        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "null-mi.bin"
            binary = self.run_generator(
                *self.COMMON, "--format", "binary", "--output", str(output)
            )
            self.assertEqual(
                binary.returncode, 0, binary.stderr.decode(errors="replace")
            )
            self.assertEqual(binary.stdout, b"")
            payload = output.read_bytes()
            self.assertEqual(len(payload), 8 * len(expected))
            self.assertEqual(
                list(struct.unpack("<{}d".format(len(expected)), payload)), expected
            )

            metadata = Path(str(output) + ".meta").read_text(encoding="ascii")
            self.assertIn("format=sjaracne-apmi-null-binary-v1\n", metadata)
            self.assertIn("m=80\n", metadata)
            self.assertIn("draws=8\n", metadata)
            self.assertIn("npar_limit=40\n", metadata)
            self.assertIn("seed=17\n", metadata)
            self.assertIn("dtype=float64\n", metadata)
            self.assertIn("byte_order=little\n", metadata)

    def test_invalid_or_incomplete_requests_fail_clearly(self):
        cases = (
            (("--m", "1", "--draws", "2"), b"--m must be at least 2"),
            (("--m", "8", "--draws", "0"), b"--draws must be positive"),
            (
                ("--m", "8", "--draws", "2", "--format", "binary"),
                b"--format binary requires --output FILE",
            ),
        )
        for arguments, message in cases:
            with self.subTest(arguments=arguments):
                result = self.run_generator(*arguments)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stderr)

    def test_complete_m4_enumeration_matches_exact_distribution(self):
        result = self.run_generator(
            "--m", "4", "--enumerate", "--npar", "40"
        )
        self.assertEqual(result.returncode, 0, result.stderr.decode(errors="replace"))
        observed = self.parse_tsv(result.stdout)
        self.assertEqual(len(observed), math.factorial(4))
        self.assertEqual(sum(abs(value) < 1e-14 for value in observed), 16)
        self.assertEqual(sum(abs(value - math.log(2.0)) < 1e-14 for value in observed), 8)
        self.assertIn(
            b"# rng=complete-lexicographic-permutation-enumeration-v1\n",
            result.stdout,
        )


if __name__ == "__main__":
    unittest.main()
