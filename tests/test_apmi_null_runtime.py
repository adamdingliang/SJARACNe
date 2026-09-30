#!/usr/bin/env python3

import importlib.util
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def find_sjaracne_executable():
    configured = os.environ.get("SJARACNE_TEST_EXE")
    for candidate in (
        configured,
        PROJECT_ROOT / "SJARACNe" / "bin" / "sjaracne.exe",
        shutil.which("sjaracne.exe"),
    ):
        if not candidate or not Path(candidate).is_file():
            continue
        if os.name == "nt":
            with Path(candidate).open("rb") as handle:
                if handle.read(2) != b"MZ":
                    continue
        return str(candidate)
    return None


SJARACNE_EXE = find_sjaracne_executable()


def cutoff(threshold, probability, shape, scale, pvalue):
    log_ratio = math.log(pvalue / probability)
    if abs(shape) < 1e-10:
        return threshold - scale * log_ratio
    return threshold + scale / shape * math.expm1(-shape * log_ratio)


def model_fields(*, m=8, npar=20, shape=0.0):
    threshold = 0.1
    probability = 0.05
    scale = 0.02
    default_p = 1e-7
    endpoint = "none" if shape >= 0.0 else f"{threshold - scale / shape:.17g}"
    return {
        "format": "sjaracne-apmi-gpd-tail-v1",
        "kernel_schema": "sjaracne-apmi-v1",
        "estimator": "sjaracne-adaptive-partitioning",
        "sampling_null": "independent-uniform-rank-permutation",
        "rank_policy": "unique-ordinal-ranks",
        "m": str(m),
        "npar_limit": str(npar),
        "tail_model": "generalized-pareto-mle-floc0",
        "tail_threshold_quantile": "0.95",
        "tail_threshold": str(threshold),
        "tail_probability": str(probability),
        "tail_shape": str(shape),
        "tail_scale": str(scale),
        "tail_endpoint": endpoint,
        "calibration_status": "accepted",
        "calibrator_schema": "sjaracne-apmi-gpd-calibrator-v1",
        "calibrator_sha256": "d" * 64,
        "validation_method": "independent-rank-permutation-stream",
        "stability_probability": "1e-7",
        "stability_relative_range": "0.02",
        "stability_relative_tolerance": "0.1",
        "validation_family_confidence": "0.99",
        "validation_point_confidence": "0.999",
        "supported_p_min": "1e-7",
        "supported_p_max": "0.01",
        "validated_p_min": "1e-5",
        "validated_p_max": "0.01",
        "default_p": str(default_p),
        "default_p_cutoff": format(
            cutoff(threshold, probability, shape, scale, default_p), ".17g"
        ),
        "fit_draws": "1000000",
        "validation_draws": "1000000",
        "fit_seed": "17",
        "validation_seed": "23",
        "rng": "mt19937-rejection-fisher-yates-v1",
        "generator_sha256": "a" * 64,
        "fit_values_sha256": "b" * 64,
        "validation_values_sha256": "c" * 64,
        "scipy_version": "1.17.1",
    }


@unittest.skipUnless(SJARACNE_EXE, "sjaracne.exe is not built")
class TestApmiNullRuntime(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.workdir = Path(self.folder.name)
        self.expression = self.workdir / "input.exp"
        rows = [
            ("H", "H", [0, 1, 2, 3, 4, 5, 6, 7]),
            ("A", "A", [0, 2, 1, 4, 3, 6, 5, 7]),
            ("B", "B", [7, 5, 6, 3, 4, 1, 2, 0]),
        ]
        lines = ["isoformId\tgeneSymbol\t" + "\t".join(f"s{i}" for i in range(8))]
        lines.extend(
            "\t".join([accession, symbol] + [str(value) for value in values])
            for accession, symbol, values in rows
        )
        self.expression.write_text("\n".join(lines) + "\n", encoding="utf-8")
        self.hubs = self.workdir / "hubs.txt"
        self.hubs.write_text("H\n", encoding="utf-8")

    def tearDown(self):
        self.folder.cleanup()

    def write_model(self, name="null.model", fields=None, extra_lines=()):
        path = self.workdir / name
        values = model_fields() if fields is None else fields
        path.write_text(
            "\n".join(f"{key}={value}" for key, value in values.items())
            + "\n"
            + "\n".join(extra_lines)
            + ("\n" if extra_lines else ""),
            encoding="utf-8",
        )
        return path

    def run_sjaracne(self, *extra, output="out.adj"):
        command = [
            SJARACNE_EXE,
            "-i", str(self.expression),
            "-s", str(self.hubs),
            "-e", "1",
            "-S", "9",
            "-o", str(self.workdir / output),
            *extra,
        ]
        return subprocess.run(command, capture_output=True, text=True, check=False)

    @staticmethod
    def threshold_from_stdout(stdout):
        match = re.search(r"MI threshold determined for p=.*: ([^\s]+)", stdout)
        if match is None:
            raise AssertionError(stdout)
        return float(match.group(1))

    @staticmethod
    def data_rows(path):
        return [
            line
            for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith(">")
        ]

    def test_exponential_gpd_cutoff_and_header_provenance(self):
        model = self.write_model()
        result = self.run_sjaracne("-M", str(model), "-p", "0.001")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        expected = cutoff(0.1, 0.05, 0.0, 0.02, 0.001)
        self.assertAlmostEqual(self.threshold_from_stdout(result.stdout), expected, places=6)

        header = (self.workdir / "out.adj").read_text(encoding="utf-8")
        self.assertIn(
            ">  MI threshold method estimator-matched AP-MI permutation-null GPD tail\n",
            header,
        )
        self.assertIn(">  AP-MI null model format sjaracne-apmi-gpd-tail-v1\n", header)
        self.assertIn(">  AP-MI kernel schema sjaracne-apmi-v1\n", header)
        self.assertIn(">  AP-MI null model m 8\n", header)
        self.assertIn(">  AP-MI null model Npar 20\n", header)
        self.assertIn(">  AP-MI validated p min 1e-05\n", header)
        self.assertIn(">  AP-MI validated p max 0.01\n", header)
        self.assertIn(">  AP-MI cutoff tail extrapolated no\n", header)
        self.assertIn(">  AP-MI calibrator schema sjaracne-apmi-gpd-calibrator-v1\n", header)
        self.assertIn(">  AP-MI calibrator SHA256 " + "d" * 64 + "\n", header)
        self.assertIn(">  AP-MI generator SHA256 " + "a" * 64 + "\n", header)
        self.assertIn(">  AP-MI fit values SHA256 " + "b" * 64 + "\n", header)
        self.assertIn(">  AP-MI validation values SHA256 " + "c" * 64 + "\n", header)

    def test_nonzero_shape_uses_general_gpd_formula(self):
        for shape in (0.1, -0.1):
            with self.subTest(shape=shape):
                fields = model_fields(shape=shape)
                model = self.write_model(f"shape_{shape}.model", fields=fields)
                result = self.run_sjaracne("-M", str(model), "-p", "0.001")
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                expected = cutoff(0.1, 0.05, shape, 0.02, 0.001)
                self.assertAlmostEqual(
                    self.threshold_from_stdout(result.stdout), expected, places=6
                )

    def test_explicit_threshold_bypasses_model_and_matches_model_cutoff(self):
        model = self.write_model()
        modeled = self.run_sjaracne("-M", str(model), "-p", "0.001", output="modeled.adj")
        self.assertEqual(modeled.returncode, 0, modeled.stdout + modeled.stderr)
        threshold = cutoff(0.1, 0.05, 0.0, 0.02, 0.001)

        explicit = self.run_sjaracne(
            "-M", str(self.workdir / "does-not-exist.model"),
            "-p", "0.001",
            "-t", format(threshold, ".17g"),
            output="explicit.adj",
        )
        self.assertEqual(explicit.returncode, 0, explicit.stdout + explicit.stderr)
        self.assertIn("'-M' is ignored", explicit.stdout)
        self.assertEqual(
            self.data_rows(self.workdir / "modeled.adj"),
            self.data_rows(self.workdir / "explicit.adj"),
        )

    def test_parser_and_estimator_mismatches_fail_closed(self):
        cases = []
        wrong_kernel = model_fields()
        wrong_kernel["kernel_schema"] = "wrong-kernel"
        cases.append((self.write_model("kernel.model", wrong_kernel), (), "kernel_schema"))

        wrong_m = model_fields(m=7)
        cases.append((self.write_model("m.model", wrong_m), (), "calibrated for m=7"))

        wrong_npar = model_fields(npar=40)
        cases.append((self.write_model("npar.model", wrong_npar), (), "npar_limit=40"))

        rejected = model_fields()
        rejected["calibration_status"] = "rejected"
        cases.append(
            (self.write_model("rejected.model", rejected), (), "calibration_status")
        )

        invalid_calibrator_sha = model_fields()
        invalid_calibrator_sha["calibrator_sha256"] = "not-a-sha"
        cases.append(
            (
                self.write_model("calibrator_sha.model", invalid_calibrator_sha),
                (),
                "SHA-256 provenance fields",
            )
        )

        duplicate = self.write_model(
            "duplicate.model", model_fields(), extra_lines=("m=8",)
        )
        cases.append((duplicate, (), "Duplicate AP-MI null-model field 'm'"))

        invalid_default = model_fields()
        invalid_default["default_p_cutoff"] = "0.2"
        cases.append((self.write_model("cutoff.model", invalid_default), (), "inconsistent"))

        unsupported_below_default = model_fields()
        unsupported_below_default["supported_p_min"] = "1e-8"
        cases.append(
            (
                self.write_model("supported_min.model", unsupported_below_default),
                (),
                "supported_p_min must equal default_p",
            )
        )

        endpoint_beyond_support = model_fields(shape=-0.005)
        cases.append(
            (
                self.write_model("endpoint.model", endpoint_beyond_support),
                (),
                "tail_endpoint exceeds the theoretical AP-MI maximum log(m)",
            )
        )

        cutoff_beyond_support = model_fields()
        cutoff_beyond_support["tail_scale"] = "1"
        cutoff_beyond_support["default_p_cutoff"] = format(
            cutoff(0.1, 0.05, 0.0, 1.0, 1e-7), ".17g"
        )
        cases.append(
            (
                self.write_model("support.model", cutoff_beyond_support),
                (),
                "cutoff beyond the theoretical AP-MI maximum log(m)",
            )
        )

        model = self.write_model("range.model")
        cases.append((model, ("-p", "0.02"), "outside the AP-MI null model's supported range"))

        for path, extra, message in cases:
            with self.subTest(path=path.name):
                result = self.run_sjaracne("-M", str(path), *(extra or ("-p", "0.001")))
                self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                self.assertIn(message, result.stderr)

    def test_model_rejects_legacy_replacement_and_adjacency_replay(self):
        model = self.write_model()
        legacy = self.run_sjaracne("-M", str(model), "-p", "0.001", "-r", "1")
        self.assertEqual(legacy.returncode, 1, legacy.stdout + legacy.stderr)
        self.assertIn("cannot be used with legacy replacement sampling", legacy.stderr)

        replay = self.run_sjaracne(
            "-M", str(model), "-p", "0.001", "-j", str(self.workdir / "old.adj")
        )
        self.assertEqual(replay.returncode, 1, replay.stdout + replay.stderr)
        self.assertIn("replaying an existing adjacency matrix", replay.stderr)

    def test_missing_model_uses_legacy_calibration_with_explicit_warning(self):
        result = self.run_sjaracne(
            "-p", "0.001", "-H", str(PROJECT_ROOT / "SJARACNe" / "config")
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("using the legacy affine threshold calibration", result.stdout)
        header = (self.workdir / "out.adj").read_text(encoding="utf-8")
        self.assertNotIn("AP-MI null model", header)

    def test_below_validated_range_warns_and_records_extrapolation(self):
        model = self.write_model()
        result = self.run_sjaracne("-M", str(model), "-p", "1e-7")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("below the model's held-out validated_p_min=1e-05", result.stdout)
        header = (self.workdir / "out.adj").read_text(encoding="utf-8")
        self.assertIn(">  AP-MI validated p min 1e-05\n", header)
        self.assertIn(">  AP-MI cutoff tail extrapolated yes\n", header)

    def test_above_validated_range_warns_and_records_extrapolation(self):
        fields = model_fields()
        fields["validated_p_max"] = "0.002"
        model = self.write_model(fields=fields)
        result = self.run_sjaracne("-M", str(model), "-p", "0.003")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("above the model's held-out validated_p_max=0.002", result.stdout)
        header = (self.workdir / "out.adj").read_text(encoding="utf-8")
        self.assertIn(">  AP-MI validated p max 0.002\n", header)
        self.assertIn(">  AP-MI cutoff tail extrapolated yes\n", header)

    def test_absent_validated_range_is_rejected(self):
        fields = model_fields()
        fields["validated_p_min"] = "none"
        model = self.write_model(fields=fields)
        result = self.run_sjaracne("-M", str(model), "-p", "0.001")
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("must record a held-out validated_p_min", result.stderr)

    def test_bare_model_option_fails_cleanly(self):
        result = self.run_sjaracne("-M")
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("Option '-M' requires a model file", result.stderr)

    def test_explicit_threshold_requires_a_complete_finite_number(self):
        for arguments in (("-t",), ("-t", "not-a-number"), ("-t", "nan"), ("-t", "inf")):
            with self.subTest(arguments=arguments):
                result = self.run_sjaracne(*arguments)
                self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                self.assertIn("requires a finite numeric value", result.stderr)


class TestApmiNullWorkflowWiring(unittest.TestCase):
    @staticmethod
    def load_wrapper():
        path = PROJECT_ROOT / "SJARACNe" / "sjaracne.py"
        spec = importlib.util.spec_from_file_location("sjaracne_wrapper_null", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def test_cwl_threads_optional_model_without_changing_default_p(self):
        workflow = (PROJECT_ROOT / "SJARACNe" / "cwl" / "sjaracne_workflow.cwl").read_text()
        command = (PROJECT_ROOT / "SJARACNe" / "cwl" / "sjaracne.cwl").read_text()
        self.assertIn("apmi_null_model: apmi_null_model", workflow)
        self.assertIn("default: 1e-7", workflow)
        self.assertIn("prefix: -M", command)
        self.assertIn("type: File?", command.split("  apmi_null_model:", 1)[1].split("  p_value:", 1)[0])

    def test_wrapper_serializes_optional_model_file(self):
        module = self.load_wrapper()
        commands = []
        module.run_shell_command_call = commands.append
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            arguments = [
                "sjaracne", "local",
                "-e", str(root / "input.exp"),
                "-g", str(root / "hubs.txt"),
                "-M", str(root / "null.model"),
                "-o", str(root / "output"),
                "-tmp", str(root / "tmp"),
            ]
            with mock.patch.object(os.sys, "argv", arguments):
                module.main()
            rendered = (root / "output" / "sjaracne_workflow.yml").read_text()
        self.assertIn("apmi_null_model:\n  class: File\n", rendered)
        self.assertIn(
            '  path: {}'.format(json.dumps(str((root / "null.model").resolve()))),
            rendered,
        )
        self.assertEqual(len(commands), 1)


if __name__ == "__main__":
    unittest.main()
