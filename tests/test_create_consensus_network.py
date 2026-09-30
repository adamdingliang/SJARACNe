#!/usr/bin/env python3

import argparse
import importlib
import unittest
import tempfile
from pathlib import Path
from unittest import mock
from SJARACNe.bin.create_consensus_network import create_consensus_network as cn
from SJARACNe.bin.create_consensus_network import create_enhanced_consensus_network as ecn
from SJARACNe.bin.create_consensus_network import consensus_probability
from SJARACNe.bin.create_consensus_network import minimum_recurrence
from SJARACNe.bin.create_consensus_network import uprob


CONSENSUS_MODULE = importlib.import_module(
    "SJARACNe.bin.create_consensus_network"
)


class TestConsensusNetwork(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.folder = tempfile.TemporaryDirectory()
        cn('./tests/inputs/adjmat_dir', 0.05, self.folder.name)

    def test_consensus_network_3col(self):
        expected = Path('./tests/answerkey/consensus_network_3col_.txt')
        actual = Path(self.folder.name) / 'consensus_network_3col_.txt'
        self.assertEqual(expected.read_text(encoding='utf-8'), actual.read_text(encoding='utf-8'))
        
    def test_bootstrap_info(self):
        expected = Path('./tests/answerkey/bootstrap_info_.txt')
        actual = Path(self.folder.name) / 'bootstrap_info_.txt'
        self.assertEqual(expected.read_text(encoding='utf-8'), actual.read_text(encoding='utf-8'))

    def test_parameter_info(self):
        expected = Path('./tests/answerkey/parameter_info_.txt').read_text(
            encoding='utf-8'
        ).splitlines()
        output_network = Path(self.folder.name) / 'consensus_network_3col_.txt'
        actual = (Path(self.folder.name) / 'parameter_info_.txt').read_text(
            encoding='utf-8'
        ).splitlines()

        self.assertEqual(
            actual[-1],
            f'>  Output network: {output_network}',
        )
        self.assertEqual(expected, actual[:-1])

    '''    
    def test_enhanced_consensus_network(self):
        ecn('./inputs/BRCA100.exp', './answerkey/consensus_network_3col_.txt', self.folder.name + '/consensus_network_ncol_.txt')
        self.assertTrue(filecmp.cmp('./answerkey/consensus_network_ncol_.txt', self.folder.name + '/consensus_network_ncol_.txt'))
    '''

    #Testing different values on the function uprob(z-score)
    def test_uprob_100(self):
        self.assertAlmostEqual(0, uprob(100))

    def test_uprob_2(self):
        self.assertAlmostEqual(0.022750150253, uprob(2))

    def test_uprob_half(self):
        self.assertAlmostEqual(0.308518356555, uprob(0.5))

    def test_uprob_neg1(self):
        self.assertAlmostEqual(0.841344680778, uprob(-1))


class TestMinimumRecurrenceConsensus(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.workdir = Path(self.folder.name)
        self.adjacency_dir = self.workdir / "bootstraps"
        self.adjacency_dir.mkdir()

    def tearDown(self):
        self.folder.cleanup()

    def write_bootstrap(self, run_number, edges):
        body = "A"
        for target, mi in edges:
            body += "\t{}\t{}".format(target, mi)
        (self.adjacency_dir / "run_{:03d}.adj".format(run_number)).write_text(
            ">  Input file test.exp\n" + body + "\n",
            encoding="utf-8",
        )

    def write_boundary_fixture(self):
        for run_number in range(1, 8):
            edges = [("B", 0.1 + run_number * 0.01)]
            if run_number <= 6:
                edges.append(("C", 0.2))
            if run_number <= 5:
                edges.append(("D", 0.3))
            self.write_bootstrap(run_number, edges)

    def test_default_k6_is_inclusive_and_skips_probability_calculation(self):
        self.write_boundary_fixture()
        output = self.workdir / "output"

        with mock.patch.object(
            CONSENSUS_MODULE,
            "uprob",
            side_effect=AssertionError("K mode must not calculate a probability"),
        ):
            cn(self.adjacency_dir, out_dir=output)

        network = (output / "consensus_network_3col_.txt").read_text(
            encoding="utf-8"
        )
        self.assertEqual(
            network,
            "source\ttarget\tMI\nA\tB\t0.1400\nA\tC\t0.2000\n",
        )

        bootstrap_info = (output / "bootstrap_info_.txt").read_text(
            encoding="utf-8"
        )
        self.assertEqual(
            bootstrap_info,
            "Total edge tested: 3\n"
            "Consensus selection: minimum recurrence\n"
            "Bootstrap networks: 7\n"
            "Minimum recurrence: 6\n"
            "Minimum recurrence fraction: 0.857142857143\n",
        )
        self.assertNotIn("mu:", bootstrap_info)
        self.assertNotIn("sigma:", bootstrap_info)
        self.assertNotIn("Bonferroni", bootstrap_info)

        parameters = (output / "parameter_info_.txt").read_text(
            encoding="utf-8"
        )
        self.assertIn(">  Bootstrap No: 7\n", parameters)
        self.assertIn(">  Consensus selection minimum recurrence\n", parameters)
        self.assertIn(">  Minimum recurrence 6 of 7\n", parameters)
        self.assertIn(">  Minimum recurrence fraction 0.857142857143\n", parameters)

    def test_explicit_k1_retains_every_observed_ordered_edge(self):
        self.write_bootstrap(1, [("B", 0.5), ("C", 0.2)])
        output = self.workdir / "output"
        cn(self.adjacency_dir, out_dir=output, min_recurrence=1)
        self.assertEqual(
            (output / "consensus_network_3col_.txt").read_text(encoding="utf-8"),
            "source\ttarget\tMI\nA\tB\t0.5000\nA\tC\t0.2000\n",
        )

    def test_duplicate_occurrences_do_not_satisfy_recurrence(self):
        (self.adjacency_dir / "run_001.adj").write_text(
            ">  Input file test.exp\nA\tB\t0.5\tB\t0.5\nA\tB\t0.5\n",
            encoding="utf-8",
        )
        self.write_bootstrap(2, [("C", 0.2)])
        output = self.workdir / "output"
        with self.assertLogs(level="WARNING"):
            cn(self.adjacency_dir, out_dir=output, min_recurrence=2)
        self.assertEqual(
            (output / "consensus_network_3col_.txt").read_text(encoding="utf-8"),
            "source\ttarget\tMI\n",
        )

    def test_invalid_or_impossible_recurrence_fails_closed(self):
        self.write_bootstrap(1, [("B", 0.5)])
        invalid_values = (0, -1, 1.5, True, "1")
        for index, value in enumerate(invalid_values):
            output = self.workdir / "invalid_{}".format(index)
            with self.subTest(value=value):
                with self.assertRaisesRegex(
                    ValueError,
                    "minimum recurrence must",
                ):
                    cn(
                        self.adjacency_dir,
                        out_dir=output,
                        min_recurrence=value,
                    )
                self.assertFalse(
                    (output / "consensus_network_3col_.txt").exists()
                )

        output = self.workdir / "too_large"
        with self.assertRaisesRegex(
            ValueError,
            "minimum recurrence 2 exceeds the number of bootstrap networks 1",
        ):
            cn(self.adjacency_dir, out_dir=output, min_recurrence=2)
        self.assertFalse(output.exists())

        default_output = self.workdir / "default_too_large"
        with self.assertRaisesRegex(
            ValueError,
            "minimum recurrence 6 exceeds the number of bootstrap networks 1",
        ):
            cn(self.adjacency_dir, out_dir=default_output)
        self.assertFalse(default_output.exists())

    def test_explicit_none_preserves_legacy_bonferroni_mode(self):
        self.write_bootstrap(1, [("B", 0.5)])
        bonferroni_output = self.workdir / "bonferroni"
        explicit_output = self.workdir / "explicit"

        cn(self.adjacency_dir, None, bonferroni_output)
        cn(self.adjacency_dir, 0.05, explicit_output)

        for filename in (
            "consensus_network_3col_.txt",
            "bootstrap_info_.txt",
        ):
            self.assertEqual(
                (bonferroni_output / filename).read_bytes(),
                (explicit_output / filename).read_bytes(),
            )
        bonferroni_parameters = (
            bonferroni_output / "parameter_info_.txt"
        ).read_text(encoding="utf-8").splitlines()
        explicit_parameters = (
            explicit_output / "parameter_info_.txt"
        ).read_text(encoding="utf-8").splitlines()
        self.assertEqual(
            [
                line
                for line in bonferroni_parameters
                if not line.startswith(">  Output network:")
            ],
            [
                line
                for line in explicit_parameters
                if not line.startswith(">  Output network:")
            ],
        )
        bootstrap_info = (
            bonferroni_output / "bootstrap_info_.txt"
        ).read_text(encoding="utf-8")
        self.assertIn("Bonferroni corrected", bootstrap_info)
        self.assertNotIn("Minimum recurrence", bootstrap_info)

    def test_probability_and_recurrence_modes_are_mutually_exclusive(self):
        self.write_bootstrap(1, [("B", 0.5)])
        output = self.workdir / "output"
        with self.assertRaisesRegex(ValueError, "mutually exclusive"):
            cn(
                self.adjacency_dir,
                0.05,
                output,
                min_recurrence=1,
            )
        self.assertFalse(output.exists())

    def test_cli_value_parsers_fail_closed(self):
        self.assertEqual(minimum_recurrence("6"), 6)
        self.assertEqual(consensus_probability("0.01"), 0.01)
        for value in ("0", "-1", "6.0", "nan"):
            with self.subTest(minimum_recurrence=value):
                with self.assertRaises(argparse.ArgumentTypeError):
                    minimum_recurrence(value)
        for value in ("0", "-1", "1.1", "nan", "inf"):
            with self.subTest(consensus_probability=value):
                with self.assertRaises(argparse.ArgumentTypeError):
                    consensus_probability(value)


class TestConsensusDuplicateEdges(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.workdir = Path(self.folder.name)

    def tearDown(self):
        self.folder.cleanup()

    @staticmethod
    def write_bootstrap(directory, filename, body):
        directory.mkdir(parents=True, exist_ok=True)
        (directory / filename).write_text(
            ">  Input file test.exp\n" + body,
            encoding="utf-8",
        )

    @staticmethod
    def read_network(output_dir):
        return (output_dir / "consensus_network_3col_.txt").read_text(
            encoding="utf-8"
        )

    def test_duplicate_edge_occurrences_count_once_per_bootstrap(self):
        clean_dir = self.workdir / "clean"
        duplicate_dir = self.workdir / "duplicate"
        clean_output = self.workdir / "clean_output"
        duplicate_output = self.workdir / "duplicate_output"

        self.write_bootstrap(clean_dir, "run_1.adj", "A\tB\t0.5\tC\t0.2\n")
        self.write_bootstrap(clean_dir, "run_2.adj", "A\tB\t0.7\n")
        self.write_bootstrap(
            duplicate_dir,
            "run_1.adj",
            "A\tB\t0.5\tB\t0.5\tC\t0.2\nA\tB\t0.5\n",
        )
        self.write_bootstrap(duplicate_dir, "run_2.adj", "A\tB\t0.7\n")

        cn(clean_dir, 1.0, clean_output)
        with self.assertLogs(level="WARNING") as logs:
            cn(duplicate_dir, 1.0, duplicate_output)

        self.assertIn("Ignored 2 duplicate edge occurrence(s)", "\n".join(logs.output))
        self.assertEqual(
            self.read_network(duplicate_output), self.read_network(clean_output)
        )
        self.assertIn("A\tB\t0.6000\n", self.read_network(duplicate_output))
        self.assertEqual(
            (duplicate_output / "bootstrap_info_.txt").read_text(encoding="utf-8"),
            (clean_output / "bootstrap_info_.txt").read_text(encoding="utf-8"),
        )

    def test_conflicting_duplicate_mi_values_fail_closed(self):
        adjacency_dir = self.workdir / "conflicting"
        output_dir = self.workdir / "output"
        self.write_bootstrap(
            adjacency_dir,
            "run_1.adj",
            "A\tB\t0.5\nA\tB\t0.6\n",
        )

        with self.assertRaisesRegex(
            ValueError,
            "Conflicting MI values for duplicate edge A----B.*0.5 versus 0.6",
        ):
            cn(adjacency_dir, 1.0, output_dir)

        self.assertFalse((output_dir / "consensus_network_3col_.txt").exists())

    def test_equivalent_mi_text_is_an_exact_numeric_duplicate(self):
        adjacency_dir = self.workdir / "equivalent"
        output_dir = self.workdir / "output"
        self.write_bootstrap(
            adjacency_dir,
            "run_1.adj",
            "A\tB\t0.5\nA\tB\t0.5000\n",
        )

        with self.assertLogs(level="WARNING"):
            cn(adjacency_dir, 1.0, output_dir)

        self.assertIn("A\tB\t0.5000\n", self.read_network(output_dir))

    def test_nonfinite_mi_fails_closed(self):
        adjacency_dir = self.workdir / "nonfinite"
        output_dir = self.workdir / "output"
        self.write_bootstrap(adjacency_dir, "run_1.adj", "A\tB\tnan\n")

        with self.assertRaisesRegex(
            ValueError,
            "Non-finite MI value for edge A----B.*nan",
        ):
            cn(adjacency_dir, 1.0, output_dir)

        self.assertFalse((output_dir / "consensus_network_3col_.txt").exists())

    def test_reverse_ordered_edge_is_not_treated_as_a_duplicate(self):
        adjacency_dir = self.workdir / "directed"
        output_dir = self.workdir / "output"
        self.write_bootstrap(
            adjacency_dir,
            "run_1.adj",
            "A\tB\t0.5\nB\tA\t0.7\n",
        )

        cn(adjacency_dir, 1.0, output_dir)

        network = self.read_network(output_dir)
        self.assertIn("A\tB\t0.5000\n", network)
        self.assertIn("B\tA\t0.7000\n", network)


class TestEmptyConsensusNetwork(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.workdir = Path(self.folder.name)
        self.adjacency_dir = self.workdir / "bootstraps"
        self.output_dir = self.workdir / "output"
        self.expression = self.workdir / "small.exp"
        self.expression.write_text(
            "isoformId\tgeneSymbol\ts1\ts2\ts3\n"
            "A\tA\t1\t2\t3\n",
            encoding="utf-8",
        )

    def tearDown(self):
        self.folder.cleanup()

    def write_empty_bootstrap(self, filename):
        self.adjacency_dir.mkdir(parents=True, exist_ok=True)
        (self.adjacency_dir / filename).write_text(
            ">  Input file test.exp\n>  Subnetwork file hubs.txt\n",
            encoding="utf-8",
        )

    def test_empty_bootstrap_networks_produce_headers_and_metadata(self):
        self.write_empty_bootstrap("run_1.adj")
        self.write_empty_bootstrap("run_2.adj")

        network = cn(self.adjacency_dir, 0.05, self.output_dir)
        ecn(self.expression, network, self.output_dir)

        self.assertEqual(
            (self.output_dir / "consensus_network_3col_.txt").read_text(
                encoding="utf-8"
            ),
            "source\ttarget\tMI\n",
        )
        self.assertEqual(
            (self.output_dir / "consensus_network_ncol_.txt").read_text(
                encoding="utf-8"
            ),
            "source\ttarget\tsource.symbol\ttarget.symbol\tMI\tpearson\t"
            "spearman\tslope\tp-value\n",
        )
        self.assertEqual(
            (self.output_dir / "bootstrap_info_.txt").read_text(encoding="utf-8"),
            "Total edge tested: 0\n"
            "Bonferroni corrected (0.05) alpha: N/A (no edges tested)\n"
            "mu: 0.0\n"
            "sigma: 0.0\n",
        )

        parameter_info = (self.output_dir / "parameter_info_.txt").read_text(
            encoding="utf-8"
        )
        self.assertIn(">  Input file test.exp\n", parameter_info)
        self.assertIn(">  Subnetwork file hubs.txt\n", parameter_info)
        self.assertIn(">  Bootstrap No: 2\n", parameter_info)
        self.assertIn(">  Source: sjaracne2\n", parameter_info)

    def test_missing_bootstrap_files_fail_clearly(self):
        self.adjacency_dir.mkdir()

        with self.assertRaisesRegex(ValueError, "No bootstrap adjacency files found"):
            cn(self.adjacency_dir, 0.05, self.output_dir)

        self.assertFalse(
            (self.output_dir / "consensus_network_3col_.txt").exists()
        )

if __name__ == '__main__':
    unittest.main()
