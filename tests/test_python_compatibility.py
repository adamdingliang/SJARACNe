"""Contracts for the standalone cross-Python output comparator."""
import math
from pathlib import Path
import tempfile
import unittest

from python_compatibility import (ATOL, ENHANCED_COLUMNS, close_number,
                                  compare_enhanced, compare_metadata, independent_selection)


class TestPythonCompatibilityComparator(unittest.TestCase):
    def test_float_tolerance_does_not_accept_nonfinite_mismatch(self):
        self.assertTrue(close_number(0.5, 0.5 + 1e-14))
        self.assertFalse(close_number(0.5, 0.5001))
        self.assertTrue(close_number(math.nan, math.nan))
        self.assertFalse(close_number(math.nan, 0))
        self.assertTrue(close_number(math.inf, math.inf))
        self.assertFalse(close_number(math.inf, -math.inf))

    @staticmethod
    def enhanced(path, source="A", mi="0.1234", pearson="0.5"):
        path.write_text("\t".join(ENHANCED_COLUMNS) + "\n" +
                        "\t".join((source, "B", "symbol_A", "symbol_B", mi,
                                   pearson, "0.6", "1.2", "0.01")) + "\n", encoding="utf-8")

    def test_enhancement_roundoff_is_reported_but_identity_and_mi_stay_exact(self):
        with tempfile.TemporaryDirectory() as temporary:
            a, b = Path(temporary) / "a.tsv", Path(temporary) / "b.tsv"
            self.enhanced(a)
            self.enhanced(b, pearson=str(0.5 + ATOL))
            result = compare_enhanced(a, b)
            self.assertEqual(result["edges"], 1)
            self.assertGreater(result["statistics"]["pearson"]["max_abs_difference"], 0)
            for kwargs in ({"source": "C"}, {"mi": "0.1235"}, {"pearson": "0.51"}):
                self.enhanced(b, **kwargs)
                with self.subTest(kwargs=kwargs), self.assertRaises(AssertionError):
                    compare_enhanced(a, b)

    def test_empty_outputs_do_not_count_as_scientific_agreement(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "empty.tsv"
            path.write_text("\t".join(ENHANCED_COLUMNS) + "\n", encoding="utf-8")
            with self.assertRaises(AssertionError):
                compare_enhanced(path, path)

    def test_consensus_metadata_only_tolerates_numeric_roundoff(self):
        with tempfile.TemporaryDirectory() as temporary:
            a, b = Path(temporary) / "a.txt", Path(temporary) / "b.txt"
            a.write_text("Total edge tested: 2\nmu: 1.5\nsigma: 0.5\n", encoding="utf-8")
            b.write_text("Total edge tested: 2\nmu: 1.50000000000001\nsigma: 0.5\n", encoding="utf-8")
            compare_metadata(a, b)
            b.write_text("Total edge tested: 3\nmu: 1.5\nsigma: 0.5\n", encoding="utf-8")
            with self.assertRaises(AssertionError):
                compare_metadata(a, b)

    def test_independent_recurrence_checks_identity_support_and_mean_mi(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            adjacencies = [root / "seed1.adj", root / "seed2.adj"]
            adjacencies[0].write_text("> header\nA\tB\t0.2\tC\t0.1\n", encoding="utf-8")
            adjacencies[1].write_text("> header\nA\tB\t0.4\n", encoding="utf-8")
            consensus = root / "consensus.tsv"
            consensus.write_text("source\ttarget\tMI\nA\tB\t0.3000\n", encoding="utf-8")
            options = {"p_value_consensus": None, "min_recurrence": 2}
            result = independent_selection(adjacencies, options, consensus)
            self.assertEqual(result["union_edges"], 2)
            self.assertEqual(result["retained_edges"], 1)
            consensus.write_text("source\ttarget\tMI\nA\tC\t0.1000\n", encoding="utf-8")
            with self.assertRaises(AssertionError):
                independent_selection(adjacencies, options, consensus)

    def test_independent_legacy_probability_checks_selection_away_from_boundary(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            adjacencies = [root / "seed1.adj", root / "seed2.adj"]
            adjacencies[0].write_text("> header\nA\tB\t0.2\tC\t0.1\n", encoding="utf-8")
            adjacencies[1].write_text("> header\nA\tB\t0.4\n", encoding="utf-8")
            consensus = root / "consensus.tsv"
            consensus.write_text("source\ttarget\tMI\nA\tB\t0.3000\n", encoding="utf-8")
            result = independent_selection(adjacencies, {"p_value_consensus": 0.2}, consensus)
            self.assertEqual(result["mode"], "legacy normal probability")
            self.assertEqual(result["retained_edges"], 1)
            boundary = 0.5 * math.erfc(1 / math.sqrt(2))
            with self.assertRaisesRegex(AssertionError, "decision boundary"):
                independent_selection(adjacencies, {"p_value_consensus": boundary}, consensus)
