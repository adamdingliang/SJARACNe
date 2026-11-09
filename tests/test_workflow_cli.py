"""Default CLI adapter contracts without requiring either external engine."""
from contextlib import ExitStack
import importlib
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock


CLI = importlib.import_module("SJARACNe.sjaracne")


class TestWorkflowCli(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.expression = self.root / "expression with space;literal.exp"
        self.hubs = self.root / "SIG_named_hubs.txt"
        header = "isoformId\tgeneSymbol\t" + "\t".join("s%d" % n for n in range(100))
        self.expression.write_text(header + "\nA\tAlpha\t" + "\t".join(str(n) for n in range(100))
                                   + "\nB\tBeta\t" + "\t".join(str(99-n) for n in range(100)) + "\n", encoding="utf-8")
        self.hubs.write_text("A\n", encoding="utf-8")
        self.output = self.root / "output with space"
        self.base = ["local", "--network-type", "TF", "-e", str(self.expression), "-g", str(self.hubs), "-o", str(self.output)]

    def invoke(self, arguments=None, callback=None, tool=None):
        with ExitStack() as stack:
            stack.enter_context(mock.patch.object(CLI, "runtime_fingerprint", return_value={"frozen": "runtime"}))
            stack.enter_context(mock.patch.object(CLI, "require_tool", side_effect=tool or (lambda value: value)))
            stack.enter_context(mock.patch.object(CLI, "check_sprocket_version", return_value="sprocket v0.32.0"))
            launch = stack.enter_context(mock.patch.object(CLI, "run_shell_command_call", side_effect=callback or self.fake_engine))
            result = CLI.main(self.base if arguments is None else arguments)
        return result, launch

    def fake_engine(self, command):
        if "--output-dir" not in command:
            output = Path(command[command.index("--outdir") + 1])
            output.mkdir(parents=True, exist_ok=True)
            for key, filename in CLI.FINAL_OUTPUTS.items():
                (output / filename).write_text("scientific output: " + key + "\n", encoding="utf-8")
            manifest = json.loads((self.output / "sjaracne_run.json").read_text())
            if manifest['scientific_options']['consensus_order'] == 'consensus-first':
                for name in ('consensus.native.stdout.txt', 'consensus.native.stderr.txt'):
                    (output / name).write_text('', encoding='utf-8')
                prepared = output / 'parameters.json'
                prepared.write_text(json.dumps({'schema': CLI.PARAMETER_SCHEMA,
                    'runtime': manifest['runtime'], 'parameters': manifest['scientific_options'],
                    'input_sha256': manifest['original_input_sha256'],
                    'original_input_sha256': manifest['original_input_sha256']}), encoding='utf-8')
                (output / 'consensus_provenance.json').write_text(json.dumps({
                    'schema': 'sjaracne-consensus-v1', 'returncode': 0,
                    'runtime': manifest['runtime'], 'parameters': manifest['scientific_options'],
                    'prepared_settings_sha256': CLI.sha256(prepared),
                    'input_sha256': manifest['original_input_sha256'],
                    'raw_hash_verification': {'requested': manifest['verify_raw_hashes'],
                        'performed': manifest['verify_raw_hashes'],
                        'method': 'sha256' if manifest['verify_raw_hashes'] else 'not-reverified'},
                    'output_sha256': {name: CLI.sha256(output / name) for name in CLI.FINAL_OUTPUTS.values()}
                }), encoding='utf-8')
            return 0
        engine_root = Path(command[command.index("--output-dir") + 1])
        run = engine_root / "runs/SJARACNe/2026-test-run"
        run.mkdir(parents=True)
        manifest = json.loads((self.output / "sjaracne_run.json").read_text())
        values = {}
        for key, filename in CLI.FINAL_OUTPUTS.items():
            path = run / filename
            path.write_text("scientific output: " + key + "\n", encoding="utf-8")
            values["SJARACNe." + key] = str(path)
        prepared = run / "parameters.json"
        prepared.write_text(json.dumps({"schema": CLI.PARAMETER_SCHEMA,
            "runtime": manifest["runtime"], "parameters": manifest["scientific_options"],
            "input_sha256": manifest["original_input_sha256"],
            "original_input_sha256": manifest["original_input_sha256"]}), encoding="utf-8")
        values["SJARACNe.parameters"] = str(prepared)
        evidence = run / 'consensus_provenance.json'
        evidence.write_text(json.dumps({'schema': 'sjaracne-consensus-v1', 'returncode': 0,
            'runtime': manifest['runtime'], 'parameters': manifest['scientific_options'],
            'prepared_settings_sha256': CLI.sha256(prepared),
            'input_sha256': manifest['original_input_sha256'],
            'raw_hash_verification': {'requested': manifest['verify_raw_hashes'],
                'performed': manifest['verify_raw_hashes'],
                'method': 'sha256' if manifest['verify_raw_hashes'] else 'not-reverified'},
            'output_sha256': {name: CLI.sha256(run / name) for name in CLI.FINAL_OUTPUTS.values()}}), encoding='utf-8')
        values['SJARACNe.consensus_provenance'] = str(evidence)
        (run / "outputs.json").write_text(json.dumps(values), encoding="utf-8")
        return 0

    def read_job(self):
        return json.loads((self.output / "sjaracne_workflow.yml").read_text(encoding="utf-8"))

    def test_default_routes_sprocket_with_pchip_dynamic_and_tf_probability(self):
        result, launch = self.invoke()
        self.assertEqual(result, 0)
        command = launch.call_args.args[0]
        self.assertIsInstance(command, list)
        self.assertEqual(command[0], "sprocket")
        self.assertNotIn("cwltool", command)
        job = self.read_job()
        for key, value in {"profile": "current", "network_type": "TF", "apmi_null_model": "pchip",
                           "p_value_bootstrap": 0.001, "dynamic_dpi": True,
                           "subsample_spec": "80%", "min_recurrence": 20}.items():
            self.assertEqual(job["SJARACNe." + key], value)
        self.assertEqual(job["SJARACNe.exp_file"], str(self.expression))
        self.assertIn("cpu = 4", (self.output / "sprocket.toml").read_text())
        record = json.loads((self.output / "sjaracne_run.json").read_text())
        self.assertEqual(record["state"], "complete")
        self.assertEqual(record["engine"], "sprocket")
        for filename in CLI.FINAL_OUTPUTS.values():
            self.assertTrue((self.output / filename).is_file())
        self.assertEqual(len(record["published_outputs"]), 4)

    def test_raw_hash_verification_defaults_off_outside_scientific_options(self):
        self.invoke()
        job = self.read_job()
        record = json.loads((self.output / "sjaracne_run.json").read_text())
        self.assertIs(job["SJARACNe.verify_raw_hashes"], False)
        self.assertIs(record["verify_raw_hashes"], False)
        self.assertNotIn("verify_raw_hashes", record["scientific_options"])
        self.assertNotIn("verify_raw_hashes", CLI.SCIENTIFIC_FIELDS)

    def test_raw_hash_verification_opt_in_reaches_both_native_workflow_engines(self):
        for engine in ("sprocket", "cwl"):
            with self.subTest(engine=engine):
                self.output = self.root / ("verified_" + engine)
                arguments = self.base.copy()
                arguments[arguments.index("-o") + 1] = str(self.output)
                self.invoke(arguments + ["--engine", engine, "--verify-raw-hashes"])
                job = self.read_job()
                key = "SJARACNe.verify_raw_hashes" if engine == "sprocket" else "verify_raw_hashes"
                self.assertIs(job[key], True)
                record = json.loads((self.output / "sjaracne_run.json").read_text())
                self.assertIs(record["verify_raw_hashes"], True)
                self.assertNotIn("verify_raw_hashes", record["scientific_options"])

    def test_wdl_dpi_first_accepts_hash_verification_but_legacy_cwl_rejects_it(self):
        self.invoke(self.base + ["--consensus-order", "dpi-first", "--verify-raw-hashes"])
        self.assertIs(self.read_job()["SJARACNe.verify_raw_hashes"], True)
        self.output = self.root / "unsupported_cwl"
        arguments = self.base.copy()
        arguments[arguments.index("-o") + 1] = str(self.output)
        with self.assertRaises(SystemExit):
            self.invoke(arguments + ["--engine", "cwl", "--consensus-order", "dpi-first",
                                     "--verify-raw-hashes"])
        self.assertFalse(self.output.exists())

    def test_verification_flag_does_not_accept_a_string_boolean(self):
        with self.assertRaises(SystemExit):
            self.invoke(self.base + ["--verify-raw-hashes", "false"])
        self.assertFalse(self.output.exists())

    def test_raw_hash_policy_must_match_before_any_final_publication(self):
        for engine in ("sprocket", "cwl"):
            for enabled in (False, True):
                for mutation in ("missing", "requested", "performed", "method", "integer_boolean"):
                    with self.subTest(engine=engine, enabled=enabled, mutation=mutation):
                        self.output = self.root / ("policy_%s_%s_%s" % (engine, enabled, mutation))
                        arguments = self.base.copy()
                        arguments[arguments.index("-o") + 1] = str(self.output)
                        arguments += ["--engine", engine]
                        if enabled:
                            arguments += ["--verify-raw-hashes"]

                        def mutate(command):
                            self.fake_engine(command)
                            if "--output-dir" in command:
                                folder = Path(command[command.index("--output-dir") + 1]) / "runs/SJARACNe/2026-test-run"
                            else:
                                folder = Path(command[command.index("--outdir") + 1])
                            path = folder / "consensus_provenance.json"
                            receipt = json.loads(path.read_text())
                            if mutation == "missing":
                                del receipt["raw_hash_verification"]
                            elif mutation in ("requested", "performed"):
                                receipt["raw_hash_verification"][mutation] = not enabled
                            elif mutation == "method":
                                receipt["raw_hash_verification"]["method"] = "different"
                            else:
                                receipt["raw_hash_verification"]["requested"] = int(enabled)
                            path.write_text(json.dumps(receipt))

                        with self.assertRaises(SystemExit):
                            self.invoke(arguments, callback=mutate)
                        self.assertFalse((self.output / "consensus_network_ncol_.txt").exists())
                        self.assertEqual(json.loads((self.output / "sjaracne_run.json").read_text())["state"], "failed")

    def test_non_boolean_expected_verification_policy_is_rejected(self):
        for invalid in (0, 1, "false", {}, []):
            with self.subTest(invalid=invalid), self.assertRaisesRegex(ValueError, "Boolean"):
                CLI.checked_consensus_provenance("absent.json", {}, {}, "absent-settings.json",
                                                 expected_verify_raw_hashes=invalid)

    def test_explicit_smaller_k_remains_supported(self):
        self.invoke(self.base + ["-n", "6", "-k", "6"])
        self.assertEqual(self.read_job()["SJARACNe.min_recurrence"], 6)
        self.assertEqual(self.read_job()["SJARACNe.bootstrap_num"], 6)

    def test_default_k20_rejects_fewer_resamples_before_launch(self):
        with mock.patch.object(CLI, "run_shell_command_call") as launch:
            with self.assertRaises(SystemExit):
                CLI.main(self.base + ["-n", "19"])
        launch.assert_not_called()
        self.assertFalse(self.output.exists())

    def test_sig_probability_and_explicit_override(self):
        arguments = self.base.copy()
        arguments[arguments.index("TF")] = "SIG"
        self.invoke(arguments + ["-pb", "0.002"])
        self.assertEqual(self.read_job()["SJARACNe.p_value_bootstrap"], 0.002)
        self.assertEqual(self.read_job()["SJARACNe.network_type"], "SIG")

    def test_network_type_remains_required_with_explicit_probability(self):
        arguments = [item for n, item in enumerate(self.base) if n not in (1, 2)] + ["-pb", "0.001"]
        with self.assertRaises(SystemExit):
            self.invoke(arguments)
        self.assertFalse(self.output.exists())

    def test_omitted_subcommand_is_local_and_serial_means_one_task(self):
        self.invoke(self.base[1:] + ["--serial"])
        text = (self.output / "sprocket.toml").read_text()
        self.assertIn("cpu = 1", text)
        self.assertIn("concurrency = 1", text)

    def test_explicit_cwl_is_execution_choice_not_legacy_science(self):
        _, launch = self.invoke(self.base + ["--engine", "cwl"])
        self.assertEqual(launch.call_args.args[0][0], "cwltool")
        job = self.read_job()
        self.assertEqual(job["scientific_options"]["apmi_null_model"], "pchip")
        self.assertTrue(job["scientific_options"]["dynamic_dpi"])
        self.assertEqual(job["scientific_options"]["p_value_bootstrap"], 0.001)
        self.assertEqual(job["scientific_options"]["consensus_order"], "consensus-first")
        self.assertIn("consensus_before_dpi.cwl", launch.call_args.args[0][-2])
        self.assertIs(job["verify_raw_hashes"], False)

    def test_legacy_profile_defaults_restore_replacement_affine_strict_and_pc(self):
        arguments = self.base.copy()
        del arguments[1:3]
        self.invoke(arguments + ["--profile", "legacy", "--engine", "cwl"])
        job = self.read_job()
        self.assertEqual(job["sample_number"], 1)
        self.assertIsNone(job["subsample_spec"])
        self.assertIsNone(job["apmi_null_model"])
        self.assertFalse(job["dynamic_dpi"])
        self.assertEqual(job["tolerance"], 0)
        self.assertEqual(job["p_value_bootstrap"], 1e-7)
        self.assertEqual(job["p_value_consensus"], 1e-5)
        self.assertNotIn("min_recurrence", job)
        self.assertNotIn("verify_raw_hashes", job)

    def test_native_cwl_reuses_scratch_parent_with_fresh_owned_children(self):
        scratch = self.root / "shared scratch"
        scratch.mkdir()
        sibling = scratch / "keep.txt"
        sibling.write_text("user data", encoding="utf-8")
        self.invoke(self.base + ["--engine", "cwl", "-tmp", str(scratch)])
        first = json.loads((self.output / "sjaracne_run.json").read_text())
        self.output = self.root / "second output"
        self.base[-1] = str(self.output)
        self.invoke(self.base + ["--engine", "cwl", "-tmp", str(scratch)])
        second = json.loads((self.output / "sjaracne_run.json").read_text())
        self.assertNotEqual(first["engine_work_dir"], second["engine_work_dir"])
        for record in (first, second):
            self.assertEqual(Path(record["engine_work_dir"]).parent, scratch)
            self.assertEqual(Path(record["engine_output_dir"]).parent, Path(record["engine_work_dir"]))
        self.assertEqual(sibling.read_text(), "user data")

    def test_fixed_tolerance_is_explicit_default_override(self):
        self.invoke(self.base + ["--engine", "cwl", "--tolerance", "0.2"])
        self.assertFalse(self.read_job()["scientific_options"]["dynamic_dpi"])
        self.assertEqual(self.read_job()["scientific_options"]["tolerance"], 0.2)

    def test_missing_engine_fails_before_output_without_fallback(self):
        def missing(value):
            raise ValueError("Required executable not found: " + value)
        with self.assertRaises(SystemExit):
            self.invoke(tool=missing)
        self.assertFalse(self.output.exists())

    def test_help_advertises_only_local_workflow_routes(self):
        help_text = CLI.build_parser().format_help()
        self.assertIn('local', help_text)
        self.assertIn('run', help_text)
        self.assertNotRegex(help_text, r'\blsf\b')

    def test_retired_lsf_route_fails_before_dependencies_or_outputs(self):
        diagnostic = io.StringIO()
        with mock.patch.object(CLI, 'require_tool') as require_tool, \
                mock.patch.object(CLI, 'runtime_fingerprint') as runtime, \
                mock.patch.object(CLI, 'run_shell_command_call') as launch, \
                mock.patch.object(CLI.sys, 'stderr', diagnostic):
            with self.assertRaises(SystemExit) as result:
                CLI.main(['lsf'] + self.base[1:] + ['--profile', 'legacy'])
        self.assertEqual(result.exception.code, 2)
        self.assertIn('removed', diagnostic.getvalue().lower())
        require_tool.assert_not_called()
        runtime.assert_not_called()
        launch.assert_not_called()
        self.assertFalse(self.output.exists())

    def test_retired_scheduler_flags_are_not_local_options(self):
        for flag in ('-j', '--config-json'):
            with self.subTest(flag=flag), mock.patch.object(CLI, 'require_tool') as require_tool, \
                    mock.patch.object(CLI, 'run_shell_command_call') as launch, \
                    mock.patch.object(CLI.sys, 'stderr', io.StringIO()):
                with self.assertRaises(SystemExit):
                    CLI.main(self.base + [flag, 'scheduler.json'])
                require_tool.assert_not_called()
                launch.assert_not_called()
                self.assertFalse(self.output.exists())

    def test_cwl_rejects_explicit_job_limit_instead_of_ignoring_it(self):
        with self.assertRaises(SystemExit):
            self.invoke(self.base + ["--engine", "cwl", "--jobs", "2"])
        self.assertFalse(self.output.exists())

    def test_nonempty_output_is_never_clobbered(self):
        self.output.mkdir()
        keep = self.output / "keep.txt"
        keep.write_text("user data")
        with self.assertRaises(SystemExit):
            self.invoke()
        self.assertEqual(keep.read_text(), "user data")
        self.assertEqual(list(self.output.iterdir()), [keep])

    def test_conflicting_dynamic_fixed_or_legacy_modern_inputs_fail_before_output(self):
        cases = (["--dynamic-dpi", "--tolerance", "0"], ["--profile", "legacy", "-M", "pchip"],
                 ["--fixed-dpi", "--dpi-slope", "0.05"], ["--memory", "1 GiB"],
                 ["--sprocket-config", "anything.toml", "--jobs", "2"])
        for extra in cases:
            with self.subTest(extra=extra), self.assertRaises(SystemExit):
                self.invoke(self.base + extra)
            self.assertFalse(self.output.exists())

    def test_failed_engine_records_failure_and_retains_evidence(self):
        def failure(command):
            raise subprocess.CalledProcessError(7, command)
        with self.assertRaises(SystemExit):
            self.invoke(callback=failure)
        record = json.loads((self.output / "sjaracne_run.json").read_text())
        self.assertEqual(record["state"], "failed")
        self.assertIn("returned non-zero exit status 7", record["error"])
        self.assertTrue((self.output / "sjaracne_workflow.yml").is_file())

    def test_success_exit_without_all_cwl_outputs_is_a_failed_workflow(self):
        with self.assertRaises(SystemExit):
            self.invoke(self.base + ["--engine", "cwl"], callback=lambda command: 0)
        record = json.loads((self.output / "sjaracne_run.json").read_text())
        self.assertEqual(record["state"], "failed")
        self.assertIn("required final output", record["error"])

    def test_execution_path_prefers_fingerprinted_package_native_binary(self):
        environment = CLI.execution_env()
        entries = environment["PATH"].split(CLI.os.pathsep)
        self.assertEqual(entries[0], str(CLI.PACKAGE / "bin"))
        self.assertEqual(entries[1], str(Path(CLI.sys.executable).parent))

    def test_output_evidence_cannot_point_outside_fresh_engine_tree(self):
        def escape(command):
            self.fake_engine(command)
            engine_root = Path(command[command.index("--output-dir") + 1])
            outputs = engine_root / "runs/SJARACNe/2026-test-run/outputs.json"
            values = json.loads(outputs.read_text())
            values["SJARACNe.out_dir"] = str(self.expression)
            outputs.write_text(json.dumps(values))
        with self.assertRaises(SystemExit):
            self.invoke(callback=escape)
        self.assertFalse((self.output / "consensus_network_ncol_.txt").exists())
        self.assertEqual(json.loads((self.output / "sjaracne_run.json").read_text())["state"], "failed")

    def test_foreign_runtime_or_different_options_are_rejected_before_publication(self):
        for field in ("runtime", "parameters"):
            with self.subTest(field=field):
                output = self.root / ("out_" + field)
                self.output = output
                arguments = self.base.copy()
                arguments[arguments.index("-o") + 1] = str(output)
                def mutate(command):
                    self.fake_engine(command)
                    engine_root = Path(command[command.index("--output-dir") + 1])
                    prepared = engine_root / "runs/SJARACNe/2026-test-run/parameters.json"
                    values = json.loads(prepared.read_text())
                    values[field] = {"foreign": True}
                    prepared.write_text(json.dumps(values))
                with self.assertRaises(SystemExit):
                    self.invoke(arguments, callback=mutate)
                self.assertFalse((output / "consensus_network_ncol_.txt").exists())

    def test_sprocket_rejects_crossdataset_prepared_or_receipt_before_publication(self):
        for field in ("prepared_settings_sha256", "input_sha256", "original_input_sha256"):
            with self.subTest(field=field):
                self.output = self.root / ("foreign_" + field)
                arguments = self.base.copy()
                arguments[arguments.index("-o") + 1] = str(self.output)
                def mutate(command):
                    self.fake_engine(command)
                    folder = Path(command[command.index("--output-dir") + 1]) / "runs/SJARACNe/2026-test-run"
                    path = folder / "consensus_provenance.json"
                    receipt = json.loads(path.read_text())
                    if field == "original_input_sha256":
                        prepared = folder / "parameters.json"
                        value = json.loads(prepared.read_text())
                        value[field] = {"foreign": "dataset"}
                        prepared.write_text(json.dumps(value))
                        receipt["prepared_settings_sha256"] = CLI.sha256(prepared)
                    else:
                        receipt[field] = {"foreign": "dataset"} if field == "input_sha256" else "foreign-settings"
                    path.write_text(json.dumps(receipt))
                with self.assertRaises(SystemExit):
                    self.invoke(arguments, callback=mutate)
                self.assertFalse((self.output / "consensus_network_ncol_.txt").exists())

    def test_scratch_parent_contents_survive_and_execution_uses_new_child(self):
        scratch = self.root / "scratch"
        scratch.mkdir()
        keep = scratch / "keep"
        keep.write_text("unchanged")
        self.invoke(self.base + ["-tmp", str(scratch)])
        record = json.loads((self.output / "sjaracne_run.json").read_text())
        self.assertEqual(Path(record["engine_output_dir"]).parent, scratch)
        self.assertTrue(Path(record["engine_output_dir"]).name.startswith("sjaracne-sprocket-"))
        self.assertEqual(keep.read_text(), "unchanged")

    def test_version_mismatch_and_missing_tool_are_clear_errors(self):
        with mock.patch.object(CLI.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "sprocket v0.31.0", "")):
            with self.assertRaisesRegex(ValueError, "tested with Sprocket 0.32.0"):
                CLI.check_sprocket_version("sprocket")
        with mock.patch.object(CLI.shutil, "which", return_value=None):
            with self.assertRaisesRegex(ValueError, "no engine fallback"):
                CLI.require_tool("missing-sprocket")


if __name__ == "__main__":
    unittest.main()
