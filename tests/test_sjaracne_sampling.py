#!/usr/bin/env python3

import math
import importlib.util
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


@unittest.skipUnless(SJARACNE_EXE, "sjaracne.exe is not built")
class TestWithoutReplacementSampling(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.workdir = Path(self.folder.name)
        self.hubs = self.workdir / "hubs.txt"
        self.hubs.write_text("H\n", encoding="utf-8")

    def tearDown(self):
        self.folder.cleanup()

    def write_expression(self, observation_count=10):
        expression = self.workdir / "input.exp"
        names = [f"s{i}" for i in range(observation_count)]
        control = list(range(observation_count))
        if observation_count == 10:
            control = [9, 0, 8, 1, 7, 2, 6, 3, 5, 4]
        rows = [
            ("CONTROL", "CONTROL", control),
            ("H", "H", list(range(observation_count))),
            ("A", "A", [(i * 3) % observation_count for i in range(observation_count)]),
            ("B", "B", [(i * 7 + 1) % observation_count for i in range(observation_count)]),
            ("C", "C", [(i * 9 + 2) % observation_count for i in range(observation_count)]),
        ]
        lines = ["isoformId\tgeneSymbol\t" + "\t".join(names)]
        for accession, symbol, values in rows:
            lines.append(
                "\t".join([accession, symbol] + [str(value) for value in values])
            )
        expression.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return expression

    def run_sjaracne(self, *extra_args, observation_count=10, output_name="out.adj"):
        expression = self.write_expression(observation_count)
        output = self.workdir / output_name
        command = [
            SJARACNE_EXE,
            "-i",
            str(expression),
            "-s",
            str(self.hubs),
            "-S",
            "17",
            "-t",
            "0",
            "-e",
            "1",
            "-v",
            "on",
            "-o",
            str(output),
            *extra_args,
        ]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        return result, output

    @staticmethod
    def selected_indices(stdout):
        match = re.search(
            r"^\[SAMPLING\] Selected original observation indices \(0-based\):(.*)$",
            stdout,
            flags=re.MULTILINE,
        )
        if match is None:
            raise AssertionError("sampling index trace was not reported")
        return [int(value) for value in match.group(1).split()]

    @staticmethod
    def data_rows(path):
        return [
            line
            for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith(">")
        ]

    def test_percentage_uses_ceiling_and_selects_unique_observations(self):
        cases = (
            (10, [0, 2, 4, 5, 6, 7, 8, 9]),
            (11, [0, 1, 3, 4, 5, 6, 7, 8, 9]),
        )
        for observation_count, expected in cases:
            with self.subTest(observations=observation_count):
                result, _ = self.run_sjaracne(
                    "-u", "80%", observation_count=observation_count
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                selected = self.selected_indices(result.stdout)
                self.assertEqual(selected, expected)
                self.assertEqual(len(set(selected)), len(expected))
                self.assertEqual(selected, sorted(selected))
                self.assertTrue(all(0 <= value < observation_count for value in selected))
                self.assertIn(
                    f"selected {len(expected)} of {observation_count} eligible observations",
                    result.stdout,
                )

    def test_exact_size_is_reproducible_and_seed_dependent(self):
        first, first_output = self.run_sjaracne("-u", "6", output_name="first.adj")
        second, second_output = self.run_sjaracne("-u", "6", output_name="second.adj")
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
        expected = [0, 2, 5, 6, 7, 9]
        self.assertEqual(self.selected_indices(first.stdout), expected)
        self.assertEqual(self.selected_indices(second.stdout), expected)
        self.assertEqual(self.data_rows(first_output), self.data_rows(second_output))

        selections = {tuple(self.selected_indices(first.stdout))}
        included = set(expected)
        for seed in range(1, 9):
            result, _ = self.run_sjaracne(
                "-u", "6", "-S", str(seed), output_name=f"seed_{seed}.adj"
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            selected = self.selected_indices(result.stdout)
            self.assertEqual(len(selected), len(set(selected)))
            self.assertTrue(all(0 <= value < 10 for value in selected))
            selections.add(tuple(selected))
            included.update(selected)
        self.assertGreater(len(selections), 1)
        self.assertEqual(included, set(range(10)))

    def test_sampling_entire_population_matches_unsampled_network(self):
        unsampled, unsampled_output = self.run_sjaracne(output_name="full.adj")
        sampled, sampled_output = self.run_sjaracne(
            "-u", "100%", output_name="sampled_full.adj"
        )
        self.assertEqual(unsampled.returncode, 0, unsampled.stdout + unsampled.stderr)
        self.assertEqual(sampled.returncode, 0, sampled.stdout + sampled.stderr)
        self.assertEqual(self.data_rows(unsampled_output), self.data_rows(sampled_output))

        exact, exact_output = self.run_sjaracne(
            "-u", "10", output_name="sampled_exact_full.adj"
        )
        self.assertEqual(exact.returncode, 0, exact.stdout + exact.stderr)
        self.assertEqual(self.data_rows(unsampled_output), self.data_rows(exact_output))

    def test_minimum_exact_sample_size_runs_end_to_end(self):
        result, output = self.run_sjaracne("-u", "2", output_name="minimum.adj")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        selected = self.selected_indices(result.stdout)
        self.assertEqual(len(selected), 2)
        self.assertEqual(len(set(selected)), 2)
        self.assertIn("x 2 observations", result.stdout)
        self.assertTrue(output.is_file())

    def test_automatic_filename_distinguishes_sampling_seeds(self):
        expression = self.write_expression(10)
        for seed in (17, 18):
            result = subprocess.run(
                [
                    SJARACNE_EXE,
                    "-i",
                    str(expression),
                    "-s",
                    str(self.hubs),
                    "-u",
                    "80%",
                    "-S",
                    str(seed),
                    "-t",
                    "100",
                    "-e",
                    "1",
                ],
                cwd=str(self.workdir),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

        self.assertEqual(
            {path.name for path in self.workdir.glob("input*_u80pct_S*.adj")},
            {
                "input_k99_t1e+02_u80pct_S17.adj",
                "input_k99_t1e+02_u80pct_S18.adj",
            },
        )

    def test_sampled_network_matches_materialized_selected_columns(self):
        sampled, sampled_output = self.run_sjaracne(
            "-u", "80%", "-n", "0.01", output_name="sampled.adj"
        )
        self.assertEqual(sampled.returncode, 0, sampled.stdout + sampled.stderr)
        selected = self.selected_indices(sampled.stdout)

        full_fields = [
            line.split("\t")
            for line in (self.workdir / "input.exp")
            .read_text(encoding="utf-8")
            .splitlines()
        ]
        materialized = self.workdir / "materialized.exp"
        materialized.write_text(
            "\n".join(
                "\t".join(row[:2] + [row[index + 2] for index in selected])
                for row in full_fields
            )
            + "\n",
            encoding="utf-8",
        )
        materialized_output = self.workdir / "materialized.adj"
        result = subprocess.run(
            [
                SJARACNE_EXE,
                "-i",
                str(materialized),
                "-s",
                str(self.hubs),
                "-S",
                "17",
                "-t",
                "0",
                "-e",
                "1",
                "-n",
                "0.01",
                "-o",
                str(materialized_output),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(
            self.data_rows(sampled_output), self.data_rows(materialized_output)
        )

    def test_conditional_population_is_sampled_before_mi(self):
        result, _ = self.run_sjaracne("-c", "-_CONTROL", "0.6", "-u", "80%")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        selected = self.selected_indices(result.stdout)
        self.assertEqual(selected, [1, 3, 7, 9, 8])
        self.assertEqual(len(set(selected)), 5)
        self.assertIn("selected 5 of 6 eligible observations", result.stdout)

        too_large, output = self.run_sjaracne(
            "-c", "-_CONTROL", "0.6", "-u", "7", output_name="too_large.adj"
        )
        self.assertEqual(too_large.returncode, 1, too_large.stdout + too_large.stderr)
        self.assertIn("only 6 are eligible", too_large.stderr)
        self.assertFalse(output.exists())

    def test_invalid_sampling_requests_fail_closed(self):
        cases = {
            "not-a-size": "invalid observation count",
            "0%": "must be within (0%,100%]",
            "101%": "must be within (0%,100%]",
            "1": "at least 2 observations",
            "11": "only 10 are eligible",
            "nan%": "invalid percentage",
        }
        for request, message in cases.items():
            with self.subTest(request=request):
                result, output = self.run_sjaracne(
                    "-u", request, output_name=f"invalid_{request.replace('%', 'pct')}.adj"
                )
                self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                self.assertIn(message, result.stderr)
                self.assertFalse(output.exists())

        missing, missing_output = self.run_sjaracne(
            "-u", output_name="missing_request.adj"
        )
        self.assertEqual(missing.returncode, 1, missing.stdout + missing.stderr)
        self.assertIn("Option '-u' requires", missing.stderr)
        self.assertFalse(missing_output.exists())

    def test_legacy_bootstrap_and_unique_sampling_are_mutually_exclusive(self):
        result, output = self.run_sjaracne("-r", "1", "-u", "80%")
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("cannot be used together", result.stderr)
        self.assertFalse(output.exists())

        adjacency = self.workdir / "input.adj"
        adjacency.write_text("H\tA\t0.5\n", encoding="utf-8")
        replay, replay_output = self.run_sjaracne(
            "-j", str(adjacency), "-u", "80%", output_name="replay.adj"
        )
        self.assertEqual(replay.returncode, 1, replay.stdout + replay.stderr)
        self.assertIn("cannot be used with an existing adjacency matrix", replay.stderr)
        self.assertFalse(replay_output.exists())

    def test_p_value_threshold_is_calibrated_for_sampled_m(self):
        expression = self.write_expression(10)
        output = self.workdir / "threshold.adj"
        p_value = 1e-7
        result = subprocess.run(
            [
                SJARACNE_EXE,
                "-i",
                str(expression),
                "-s",
                str(self.hubs),
                "-u",
                "80%",
                "-p",
                str(p_value),
                "-e",
                "1",
                "-H",
                str(PROJECT_ROOT / "SJARACNe" / "config"),
                "-o",
                str(output),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        alpha, beta, gamma = 1.062, -48.7, -0.634
        expected = (alpha - math.log(p_value)) / (-beta - gamma * 8)
        match = re.search(r"MI threshold determined for p=.*: ([^\s]+)", result.stdout)
        self.assertIsNotNone(match, result.stdout)
        self.assertAlmostEqual(float(match.group(1)), expected, places=5)
        self.assertNotAlmostEqual(
            float(match.group(1)),
            (alpha - math.log(p_value)) / (-beta - gamma * 10),
            places=5,
        )

        header = output.read_text(encoding="utf-8")
        self.assertIn(">  Sampling method fixed-size without replacement\n", header)
        self.assertIn(">  Sampling request 80%\n", header)
        self.assertIn(">  Eligible observations 10\n", header)
        self.assertIn(">  Sampled observations 8\n", header)


class TestWorkflowSamplingConfiguration(unittest.TestCase):
    @staticmethod
    def load_wrapper():
        path = PROJECT_ROOT / "SJARACNe" / "sjaracne.py"
        spec = importlib.util.spec_from_file_location("sjaracne_wrapper", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def run_wrapper(self, *sampling_args):
        module = self.load_wrapper()
        commands = []
        module.run_shell_command_call = commands.append

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "output"
            temp_dir = root / "tmp"
            arguments = [
                "sjaracne",
                "local",
                "-e",
                str(root / "input.exp"),
                "-g",
                str(root / "hubs.txt"),
                "-o",
                str(output),
                "-tmp",
                str(temp_dir),
                *sampling_args,
            ]
            with mock.patch.object(os.sys, "argv", arguments):
                module.main()
            workflow = (output / "sjaracne_workflow.yml").read_text(
                encoding="utf-8"
            )
        self.assertEqual(len(commands), 1)
        return workflow

    def test_wrapper_defaults_to_eighty_percent(self):
        workflow = self.run_wrapper()
        self.assertIn('subsample_spec: "80%"\n', workflow)

    def test_wrapper_defaults_to_minimum_recurrence_six(self):
        workflow = self.run_wrapper()
        self.assertIn('min_recurrence: 6\n', workflow)
        self.assertNotIn('p_value_consensus:', workflow)

    def test_wrapper_serializes_recurrence_or_explicit_legacy_probability(self):
        recurrence = self.run_wrapper('--min-recurrence', '8')
        legacy = self.run_wrapper('--p-value-consensus', '0.01')
        self.assertIn('min_recurrence: 8\n', recurrence)
        self.assertNotIn('p_value_consensus:', recurrence)
        self.assertIn('p_value_consensus: 0.01\n', legacy)
        self.assertNotIn('min_recurrence:', legacy)

        with self.assertRaises(SystemExit):
            self.run_wrapper(
                '--min-recurrence',
                '6',
                '--p-value-consensus',
                '0.01',
            )

    def test_wrapper_preserves_explicit_brca100_bootstrap_probabilities(self):
        tf = self.run_wrapper('--p-value-bootstrap', '1e-3')
        sig = self.run_wrapper('--p-value-bootstrap', '5e-4')
        self.assertIn('p_value_bootstrap: 1e-3\n', tf)
        self.assertIn('p_value_bootstrap: 5e-4\n', sig)

    def test_wrapper_rejects_recurrence_above_bootstrap_count_early(self):
        with self.assertRaises(SystemExit):
            self.run_wrapper('--bootstrap-num', '5')
        with self.assertRaises(SystemExit):
            self.run_wrapper(
                '--bootstrap-num', '6', '--min-recurrence', '7'
            )

        legacy = self.run_wrapper(
            '--bootstrap-num', '5', '--p-value-consensus', '0.01'
        )
        self.assertIn('bootstrap_num: 5\n', legacy)
        self.assertIn('p_value_consensus: 0.01\n', legacy)

    def test_wrapper_serializes_fraction_or_exact_size(self):
        fraction = self.run_wrapper("--subsample-fraction", "0.64")
        exact = self.run_wrapper("--subsample-size", "7")
        self.assertIn('subsample_spec: "64%"\n', fraction)
        self.assertIn('subsample_spec: "7"\n', exact)

    def test_cwl_defaults_and_wiring_use_unique_sampling(self):
        workflow = (PROJECT_ROOT / "SJARACNe" / "cwl" / "sjaracne_workflow.cwl").read_text(
            encoding="utf-8"
        )
        command = (PROJECT_ROOT / "SJARACNe" / "cwl" / "sjaracne.cwl").read_text(
            encoding="utf-8"
        )
        self.assertIn('default: "80%"', workflow)
        self.assertIn("subsample_spec: subsample_spec", workflow)
        self.assertIn("prefix: -u", command)
        sample_section = command.split("  sample_number:", 1)[1].split(
            "  subsample_spec:", 1
        )[0]
        self.assertNotIn("default:", sample_section)

    def test_cwl_wires_both_consensus_selection_modes(self):
        workflow = (
            PROJECT_ROOT / "SJARACNe" / "cwl" / "sjaracne_workflow.cwl"
        ).read_text(encoding="utf-8")
        command = (
            PROJECT_ROOT / "SJARACNe" / "cwl" / "create_consensus_network.cwl"
        ).read_text(encoding="utf-8")

        self.assertIn("min_recurrence:\n    type: int?", workflow)
        self.assertIn("p_value_consensus:\n    type: float?", workflow)
        self.assertNotIn("p_value_consensus:\n    type: float\n    default:", workflow)
        self.assertIn("min_recurrence: min_recurrence", workflow)
        self.assertIn("p_thresh_arg: p_value_consensus", workflow)
        self.assertIn("min_recurrence:\n    type: int?", command)
        self.assertIn("prefix: -k", command)
        self.assertIn("p_thresh_arg:\n    type: float?", command)
        self.assertIn(
            "glob: $(inputs.output_dir)/bootstrap_info_.txt",
            command,
        )
        self.assertIn(
            "glob: $(inputs.output_dir)/parameter_info_.txt",
            command,
        )
        self.assertIn("outputSource: consensus/bootstrap_info", workflow)
        self.assertIn("outputSource: consensus/parameter_info", workflow)
        self.assertIn(
            "out: [out_dir, bootstrap_info, parameter_info]",
            workflow,
        )


if __name__ == "__main__":
    unittest.main()
