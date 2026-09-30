#!/usr/bin/env python3

import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "SJARACNe" / "calibrate_apmi_null.py"


def load_calibration_module():
    specification = importlib.util.spec_from_file_location("apmi_null_calibration", SCRIPT)
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


CALIBRATION = load_calibration_module()


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


class TestApmiNullTailMath(unittest.TestCase):
    def test_exponential_and_bounded_gpd_cutoffs(self):
        exponential = CALIBRATION.TailFit(
            0.95, 0.1, 5000, 0.05, 0.0, 0.02, None, 0.0
        )
        expected = 0.1 + 0.02 * math.log(0.05 / 1e-4)
        self.assertAlmostEqual(CALIBRATION.gpd_cutoff(exponential, 1e-4), expected)

        bounded = CALIBRATION.TailFit(
            0.95, 0.1, 5000, 0.05, -0.25, 0.02, 0.18, 0.0
        )
        self.assertLess(CALIBRATION.gpd_cutoff(bounded, 1e-7), 0.18)
        self.assertGreater(CALIBRATION.gpd_cutoff(bounded, 1e-7), 0.1)

    def test_validated_range_stops_at_first_failed_level(self):
        rows = [
            {"nominal_p": 1e-2, "nominal_in_interval": True},
            {"nominal_p": 1e-3, "nominal_in_interval": True},
            {"nominal_p": 1e-4, "nominal_in_interval": False},
            {"nominal_p": 1e-5, "nominal_in_interval": True},
        ]
        self.assertEqual(CALIBRATION.contiguous_validated_minimum(rows), 1e-3)
        rows[0]["nominal_in_interval"] = False
        self.assertIsNone(CALIBRATION.contiguous_validated_minimum(rows))

    @staticmethod
    def valid_arguments(**overrides):
        values = {
            "generator": SCRIPT,
            "m": [80],
            "npar": 40,
            "fit_draws": 1000,
            "validation_draws": 1000,
            "tail_quantile": 0.8,
            "min_tail_count": 100,
            "p_min": 1e-7,
            "p_max": 0.01,
            "default_p": 1e-7,
            "validation_p": [0.01, 0.001],
            "fit_seed": 17,
            "validation_seed": 23,
            "stability_relative_tolerance": 0.1,
            "min_stability_fits": 2,
            "validation_family_confidence": 0.99,
            "validation_family_comparisons": 2,
        }
        values.update(overrides)
        return SimpleNamespace(**values)

    def test_v1_requires_default_p_to_equal_supported_minimum(self):
        arguments = self.valid_arguments(p_min=1e-8)
        with self.assertRaisesRegex(ValueError, "p-min to equal --default-p"):
            CALIBRATION.validate_args(arguments)

    def test_cross_m_fit_and_validation_seed_collision_is_rejected(self):
        arguments = self.valid_arguments(fit_seed=1_000_000_000, validation_seed=0)
        with self.assertRaisesRegex(ValueError, "seed streams overlap"):
            CALIBRATION.validate_args(arguments)


@unittest.skipUnless(GENERATOR, "apmi_null_generator.exe is not built")
class TestApmiNullCalibrationSmoke(unittest.TestCase):
    def test_small_end_to_end_run_writes_finite_strict_artifacts(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "calibration"
            output.mkdir()
            stale_model = output / "apmi_null_m00032_npar040.model"
            stale_model.write_text("stale model\n", encoding="utf-8")
            command = [
                sys.executable,
                str(SCRIPT),
                "--generator",
                GENERATOR,
                "--output-dir",
                str(output),
                "--m",
                "32",
                "--npar",
                "40",
                "--fit-draws",
                "5000",
                "--validation-draws",
                "5000",
                "--tail-quantile",
                "0.8",
                "--min-tail-count",
                "100",
                "--min-stability-fits",
                "2",
                "--stability-relative-tolerance",
                "0.99",
                "--p-max",
                "0.02",
                "--validation-p",
                "0.02",
                "--allow-rejected-calibration",
            ]
            completed = subprocess.run(command, text=True, capture_output=True, check=False)
            self.assertEqual(completed.returncode, 0, completed.stderr)

            report = json.loads(
                (output / "calibration_report.json").read_text(encoding="utf-8")
            )
            self.assertEqual(report["format"], "sjaracne-apmi-null-calibration-report-v1")
            result = report["results"][0]
            self.assertEqual(result["m"], 32)
            self.assertEqual(result["model"]["npar_limit"], 40)
            self.assertEqual(result["model"]["default_p"], 1e-7)
            self.assertTrue(math.isfinite(result["model"]["default_p_cutoff"]))
            self.assertEqual(
                result["model"]["calibration_status"],
                "accepted" if result["accepted"] else "rejected",
            )
            self.assertIn("validated_p_max", result["model"])

            model = output / "apmi_null_m00032_npar040.model"
            self.assertEqual(model.is_file(), bool(result["accepted"]))
            if model.is_file():
                model_text = model.read_text(encoding="utf-8")
                self.assertIn("format=sjaracne-apmi-gpd-tail-v1\n", model_text)
                self.assertIn("kernel_schema=sjaracne-apmi-v1\n", model_text)
                self.assertIn("m=32\n", model_text)
                self.assertIn("npar_limit=40\n", model_text)
                self.assertNotIn("nan", model_text.lower())
                self.assertNotIn("inf", model_text.lower())
                self.assertNotEqual(model_text, "stale model\n")

    @unittest.skipIf(os.name == "nt", "POSIX executable fixture")
    def test_allow_rejected_does_not_hide_generator_failure(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            failing = root / "failing-generator"
            failing.write_text("#!/bin/sh\nexit 9\n", encoding="utf-8")
            failing.chmod(0o755)
            command = [
                sys.executable,
                str(SCRIPT),
                "--generator", str(failing),
                "--output-dir", str(root / "output"),
                "--m", "32",
                "--fit-draws", "1000",
                "--validation-draws", "1000",
                "--min-tail-count", "100",
                "--allow-rejected-calibration",
            ]
            completed = subprocess.run(command, text=True, capture_output=True)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("null generator failed (9)", completed.stderr)


if __name__ == "__main__":
    unittest.main()
