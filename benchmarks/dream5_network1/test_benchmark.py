"""Scorer tests and bundled-input integrity check; no native binary required."""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from benchmark import adjacency, digest, gold_labels, load_inputs, ranking_metrics, score


class BenchmarkTests(unittest.TestCase):
    def test_bundled_network1_inputs(self):
        root = Path(__file__).resolve().parent / "data" / "Network1"
        inputs = load_inputs(data_root=root)
        provenance = json.loads((root / "provenance.json").read_text(encoding="utf-8"))
        self.assertEqual(
            digest((root / "dream5_net1.exp").read_bytes()),
            provenance["prepared_expression_sha256"],
        )
        self.assertEqual(
            (root / "all195_tfs.txt").read_text(encoding="utf-8").splitlines(),
            inputs["tfs"].splitlines(),
        )
        self.assertEqual(len(gold_labels(inputs["gold"])), 278392)

    def test_ties_and_average_precision(self):
        auc, ap = ranking_metrics([0.9, 0.8, 0.8, 0], [1, 0, 1, 0])
        self.assertAlmostEqual(auc, 0.875)
        self.assertAlmostEqual(ap, 5 / 6)

    def test_explicit_gold_universe_only(self):
        labels = {("T", "A"): 1, ("T", "B"): 0,
                  ("T", "C"): 1, ("T", "D"): 0}
        raw = {("T", "A"): 0.9, ("T", "B"): 0.8,
               ("T", "C"): 0.8, ("T", "X"): 0.7}
        final = {("T", "A"): 0.9, ("T", "X"): 0.7}
        got = score(raw, final, labels, ["T"])
        self.assertEqual((got["tp"], got["fp"], got["fn"], got["tn"]), (1, 0, 1, 2))
        self.assertEqual(got["final_edges_unlabeled"], 1)
        self.assertEqual(got["precision"], 1)
        self.assertEqual(got["sensitivity"], 0.5)
        self.assertAlmostEqual(got["auroc_all_unthresholded_mi"], 0.875)

    def test_duplicate_gold_pair_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            gold_labels("T\tA\t1\nT\tA\t0\n")

    def test_final_mi_must_match_raw(self):
        with self.assertRaisesRegex(ValueError, "disagrees"):
            score({("T", "A"): 0.5}, {("T", "A"): 0.4},
                  {("T", "A"): 1, ("T", "B"): 0}, ["T"])

    def test_native_adjacency_parser(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "tiny.adj"
            path.write_text(">  MI threshold 0.1\nT\tA\t0.5\tB\t0.2\n", encoding="utf-8")
            self.assertEqual(adjacency(path), {("T", "A"): 0.5, ("T", "B"): 0.2})
            path.write_text("T\tA\t0.5\nT\tA\t0.6\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                adjacency(path)


if __name__ == "__main__":
    unittest.main()
