"""Tiny real regressions for the new default engine and historical profile.

Unlike optional developer tests, these fail for a missing Sprocket/runtime.
The complete small runner and CI explicitly provision and pin that dependency.
"""
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

from native_test_support import find_native_executable
from test_sjaracne_acceptance import write_multihub_dpi_inputs


ROOT = Path(__file__).resolve().parents[1]
FINAL_FILES = (
    "consensus_network_3col_.txt", "consensus_network_ncol_.txt",
    "bootstrap_info_.txt", "parameter_info_.txt",
)
PATH_HEADER = re.compile(r"^(>\s+(?:Input file\s+|Subnetwork file\s+|Output network:\s+)).*$")


def header(text, name):
    match = re.search(r"^>\s+" + re.escape(name) + r"\s+([^\r\n]+)$", text, re.M)
    if match is None:
        raise AssertionError("Missing native metadata: " + name)
    return match.group(1).strip()


def normalized_parameters(path):
    return "\n".join(PATH_HEADER.sub(r"\1<localized-path>", line)
                     for line in path.read_text(encoding="utf-8").splitlines())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class TestDefaultWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not sys.platform.startswith("linux"):
            raise RuntimeError("Real default-engine regression requires Linux/WSL, not a skipped test")
        cls.native = find_native_executable()
        if cls.native is None:
            raise RuntimeError("Build the native executable before the default-engine regression")
        configured = os.environ.get("SJARACNE_TEST_SPROCKET")
        found = configured if configured is not None else shutil.which("sprocket")
        if not found or not Path(found).is_file() or not os.access(found, os.X_OK):
            raise RuntimeError("Sprocket 0.32.0 is required; set SJARACNE_TEST_SPROCKET or PATH")
        cls.sprocket = str(Path(found).resolve())
        version = subprocess.check_output([cls.sprocket, "--version"], text=True)
        if not re.search(r"^sprocket v?0\.32\.0(?:\s|$)", version):
            raise RuntimeError("Default-engine tests pin Sprocket 0.32.0: " + version.strip())
        cls.env = os.environ.copy()
        cls.env.update(PYTHONPATH=str(ROOT), PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1",
                       OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        cls.env["PATH"] = os.pathsep.join((str(Path(sys.executable).parent),
                                          str(ROOT / "SJARACNe/bin"), cls.env.get("PATH", "")))

    def run_cli(self, folder, expression, hubs, *options, engine=None):
        output = folder / "output"
        command = [sys.executable, "-m", "SJARACNe.sjaracne", "local",
                   "-e", str(expression), "-g", str(hubs), "-n", "6",
                   "-o", str(output), "-tmp", str(folder / "scratch")]
        if engine is not None:
            command += ["--engine", engine]
        if engine in (None, "sprocket"):
            command += ["--sprocket", self.sprocket, "--jobs", "2"]
        command += list(options)
        process = subprocess.run(command, cwd=ROOT, env=self.env, capture_output=True, text=True)
        return process, output

    def successful_run(self, folder, expression, hubs, *options, engine=None):
        process, output = self.run_cli(folder, expression, hubs, *options, engine=engine)
        self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
        manifest_path = output / "sjaracne_run.json"
        self.assertTrue(manifest_path.is_file(), process.stdout + process.stderr)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(manifest["state"], "complete", manifest)
        self.assertEqual(manifest["engine"], "sprocket" if engine is None else engine)
        for name in FINAL_FILES:
            self.assertTrue((output / name).is_file(), name)
        if manifest["engine"] == "sprocket":
            engine_root = Path(manifest["engine_output_dir"]).resolve()
            outputs_path = Path(manifest["engine_outputs_json"]).resolve()
            self.assertTrue(outputs_path.is_file())
            self.assertIn(engine_root, outputs_path.parents)
            self.assertEqual(outputs_path.parent.parent.name, "SJARACNe")
            self.assertNotEqual(outputs_path.parent.name, "_latest")
            outputs = json.loads(outputs_path.read_text(encoding="utf-8"))
            for name, key in (("consensus_network_3col_.txt", "consensus_network"),
                              ("consensus_network_ncol_.txt", "out_dir"),
                              ("bootstrap_info_.txt", "bootstrap_info"),
                              ("parameter_info_.txt", "parameter_info")):
                self.assertEqual((output / name).read_bytes(),
                                 Path(outputs["SJARACNe." + key]).read_bytes())
        else:
            outputs = None
        return output, manifest, outputs

    def test_default_tf_and_sig_use_sprocket_pchip_dynamic_and_k6(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = root / "input files"
            fixture.mkdir()
            expression, hubs = write_multihub_dpi_inputs(fixture)
            original_hashes = (sha256(expression), sha256(hubs))
            for network, probability in (("TF", 0.001), ("SIG", 0.0005)):
                with self.subTest(network=network):
                    output, manifest, outputs = self.successful_run(
                        root / (network + " with spaces"), expression, hubs, "--network-type", network)
                    options = manifest["scientific_options"]
                    self.assertEqual(options["profile"], "current")
                    self.assertEqual(options["network_type"], network)
                    self.assertEqual(options["apmi_null_model"], "pchip")
                    self.assertEqual(options["p_value_bootstrap"], probability)
                    self.assertIs(options["dynamic_dpi"], True)
                    self.assertEqual(options["sampling_method"], "without-replacement")
                    self.assertEqual(options["subsample_spec"], "80%")
                    self.assertEqual(options["depth"], 40)
                    self.assertEqual(options["min_recurrence"], 6)
                    self.assertIsNone(options["p_value_consensus"])
                    adjacencies = [Path(value) for value in outputs["SJARACNe.adjacencies"]]
                    records = [Path(value) for value in outputs["SJARACNe.run_provenance"]]
                    self.assertEqual([path.name for path in adjacencies],
                                     ["TF_run_%03d.adj" % seed for seed in range(1, 7)])
                    self.assertEqual(len(records), 6)
                    for seed, (adjacency, record) in enumerate(zip(adjacencies, records), 1):
                        text = adjacency.read_text(encoding="utf-8")
                        provenance = json.loads(record.read_text(encoding="utf-8"))
                        argv = provenance["argv"]
                        self.assertEqual(provenance["seed"], seed)
                        self.assertEqual(provenance["adjacency_sha256"], sha256(adjacency))
                        self.assertEqual(argv[argv.index("-S") + 1], str(seed))
                        self.assertEqual(argv[argv.index("-M") + 1], "pchip")
                        self.assertEqual(argv[argv.index("-u") + 1], "80%")
                        self.assertIn("--dynamic-dpi", argv)
                        self.assertNotIn("-e", argv)
                        self.assertNotIn("-r", argv)
                        self.assertEqual(header(text, "MI threshold method"), "experimental PCHIP empirical cutoff")
                        self.assertEqual(float(header(text, "MI P-value")), probability)
                        self.assertEqual(header(text, "AP-MI actual m"), "80")
                        self.assertEqual(header(text, "DPI mode"), "experimental-dynamic")
                        self.assertEqual(header(text, "DPI H_ref"), "479")
                        self.assertEqual(header(text, "DPI source rows H"), "3")
                        self.assertEqual(float(header(text, "DPI epsilon binary64")), 0)
                    metadata = (output / "bootstrap_info_.txt").read_text(encoding="utf-8")
                    self.assertIn("Consensus selection: minimum recurrence\n", metadata)
                    self.assertIn("Minimum recurrence: 6\n", metadata)
                    self.assertGreater(len((output / "consensus_network_ncol_.txt").read_bytes().splitlines()), 1)
            self.assertEqual((sha256(expression), sha256(hubs)), original_hashes)

    def test_legacy_profile_restores_replacement_affine_strict_and_probability_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            expression, hubs = write_multihub_dpi_inputs(root)
            wdl_output, manifest, outputs = self.successful_run(
                root / "wdl", expression, hubs, "--profile", "legacy")
            options = manifest["scientific_options"]
            self.assertEqual(options["profile"], "legacy")
            self.assertEqual(options["sampling_method"], "with-replacement")
            self.assertIsNone(options["subsample_spec"])
            self.assertEqual(options["apmi_null_model"], "affine")
            self.assertEqual(options["p_value_bootstrap"], 1e-7)
            self.assertIs(options["dynamic_dpi"], False)
            self.assertEqual(options["tolerance"], 0)
            self.assertEqual(options["p_value_consensus"], 1e-5)
            self.assertIsNone(options["min_recurrence"])
            expected_cutoff = (1.062 - math.log(1e-7)) / (48.7 + 0.634 * 100)
            for adjacency, record in zip(outputs["SJARACNe.adjacencies"], outputs["SJARACNe.run_provenance"]):
                text = Path(adjacency).read_text(encoding="utf-8")
                argv = json.loads(Path(record).read_text(encoding="utf-8"))["argv"]
                self.assertEqual(argv[argv.index("-r") + 1], "1")
                self.assertEqual(float(argv[argv.index("-e") + 1]), 0)
                self.assertEqual(argv[argv.index("-N") + 1], "40")
                for flag in ("-u", "-M", "--dynamic-dpi"):
                    self.assertNotIn(flag, argv)
                self.assertEqual(header(text, "Sampling method"), "legacy bootstrap with replacement")
                self.assertEqual(header(text, "Resample draws"), "100")
                self.assertTrue(math.isclose(float(header(text, "MI threshold")), expected_cutoff,
                                             rel_tol=0, abs_tol=5e-7))
            cwl_output, cwl_manifest, _ = self.successful_run(
                root / "cwl", expression, hubs, "--profile", "legacy", engine="cwl")
            self.assertEqual(cwl_manifest["scientific_options"], options)
            for name in FINAL_FILES[:3]:
                self.assertEqual((wdl_output / name).read_bytes(), (cwl_output / name).read_bytes(), name)
            self.assertEqual(normalized_parameters(wdl_output / "parameter_info_.txt"),
                             normalized_parameters(cwl_output / "parameter_info_.txt"))
            metadata = (wdl_output / "bootstrap_info_.txt").read_text(encoding="utf-8")
            self.assertNotIn("Minimum recurrence:", metadata)
            self.assertIn("Bonferroni corrected (0.05) alpha:", metadata)
            self.assertIn("mu:", metadata)
            self.assertIn("sigma:", metadata)

    def test_missing_network_type_and_legacy_conflicts_fail_before_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            expression, hubs = write_multihub_dpi_inputs(root)
            for name, options, diagnostic in (
                    ("missing", (), "network_type"),
                    ("explicit_pb_without_type", ("-pb", "0.001"), "network_type"),
                    ("legacy_dynamic", ("--profile", "legacy", "--dynamic-dpi"), "legacy")):
                with self.subTest(name=name):
                    process, output = self.run_cli(root / name, expression, hubs, *options)
                    self.assertNotEqual(process.returncode, 0)
                    self.assertIn(diagnostic, process.stdout + process.stderr)
                    self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
