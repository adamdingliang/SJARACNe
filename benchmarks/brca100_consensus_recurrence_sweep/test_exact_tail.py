#!/usr/bin/env python3
"""Independent tests for exact Poisson-binomial consensus tails."""

from __future__ import annotations

import csv
import itertools
import math
from pathlib import Path
import sys
import unittest


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_recurrence_sweep as sweep


def brute_force_tails(probabilities: list[float]) -> list[float]:
    """Enumerate every Bernoulli outcome, independently of the implementation."""

    mass_by_count = [0.0] * (len(probabilities) + 1)
    for outcome in itertools.product((0, 1), repeat=len(probabilities)):
        mass = math.prod(
            probability if observed else 1.0 - probability
            for probability, observed in zip(probabilities, outcome, strict=True)
        )
        mass_by_count[sum(outcome)] += mass
    return [sum(mass_by_count[k:]) for k in range(len(mass_by_count))]


class PoissonBinomialTailTest(unittest.TestCase):
    def assert_tails_close(
        self, actual: list[float], expected: list[float], *, rel_tol: float = 1e-13
    ) -> None:
        self.assertEqual(len(actual), len(expected))
        for support, (observed, reference) in enumerate(zip(actual, expected, strict=True)):
            self.assertTrue(
                math.isclose(observed, reference, rel_tol=rel_tol, abs_tol=1e-15),
                msg=(
                    f"tail mismatch at support {support}: "
                    f"observed={observed!r}, expected={reference!r}"
                ),
            )

    def test_empty_distribution_has_certain_zero_tail(self) -> None:
        self.assert_tails_close(sweep.poisson_binomial_tails([]), [1.0])

    def test_matches_exhaustive_enumeration_for_small_cases(self) -> None:
        cases = [
            [0.0],
            [1.0],
            [0.2, 0.7],
            [0.1, 0.35, 0.8],
            [0.0, 0.15, 0.5, 0.9, 1.0],
        ]
        for probabilities in cases:
            with self.subTest(probabilities=probabilities):
                self.assert_tails_close(
                    sweep.poisson_binomial_tails(probabilities),
                    brute_force_tails(probabilities),
                )

    def test_tail_sequence_is_monotone_and_has_expected_endpoints(self) -> None:
        probabilities = [0.03, 0.12, 0.41, 0.77, 0.98]
        tails = sweep.poisson_binomial_tails(probabilities)
        self.assertEqual(len(tails), len(probabilities) + 1)
        self.assertTrue(math.isclose(tails[0], 1.0, abs_tol=1e-15))
        self.assertTrue(
            math.isclose(
                tails[-1], math.prod(probabilities), rel_tol=1e-13, abs_tol=1e-15
            )
        )
        self.assertTrue(all(left >= right for left, right in zip(tails, tails[1:])))
        self.assertTrue(all(0.0 <= tail <= 1.0 for tail in tails))

    def test_reproduces_published_brca100_anchors(self) -> None:
        manifest_path = (
            HERE
            / "results_2026-08-20"
            / "provenance"
            / "source_sweep"
            / "run_manifest.tsv"
        )
        with manifest_path.open(newline="", encoding="utf-8") as handle:
            records = list(csv.DictReader(handle, delimiter="\t"))

        anchors = {
            ("p1e-03", "tf"): {
                "union_edges": 5_666_377,
                "tails": {
                    6: 0.03270220879190407,
                    9: 0.0006664162246044295,
                    12: 4.998826469828169e-06,
                    15: 1.680784531297306e-08,
                    16: 2.156159767372582e-09,
                },
            },
            ("p5e-04", "sig"): {
                "union_edges": 13_243_781,
                "tails": {
                    6: 0.017702567659849355,
                    9: 0.00023668891255159368,
                    11: 7.499485238057152e-06,
                    14: 2.0978977420788344e-08,
                    15: 2.5074146505000815e-09,
                },
            },
        }

        for (point, driver), anchor in anchors.items():
            with self.subTest(point=point, driver=driver):
                edge_counts = [
                    int(record["edges"])
                    for record in records
                    if record["point"] == point and record["driver"] == driver
                ]
                self.assertEqual(len(edge_counts), 100)
                union_edges = int(anchor["union_edges"])
                probabilities = [count / union_edges for count in edge_counts]
                tails = sweep.poisson_binomial_tails(probabilities)
                for support, expected in anchor["tails"].items():
                    self.assertTrue(
                        math.isclose(
                            tails[support], expected, rel_tol=1e-12, abs_tol=1e-15
                        ),
                        msg=(
                            f"{point}/{driver} tail mismatch at K={support}: "
                            f"observed={tails[support]!r}, expected={expected!r}"
                        ),
                    )

    def test_legacy_uprob_matches_sjaracne_polynomial_branches(self) -> None:
        expected = {
            0.0: 0.5,
            1.0: 0.1586553192214073,
            1.899999: 0.02907575719600064,
            1.9: 0.028716605046543553,
            5.0: 2.8665157187920085e-07,
            -1.0: 0.8413446807785927,
        }
        for z_score, reference in expected.items():
            with self.subTest(z_score=z_score):
                self.assertTrue(
                    math.isclose(
                        sweep.legacy_uprob(z_score),
                        reference,
                        rel_tol=1e-13,
                        abs_tol=1e-15,
                    )
                )


if __name__ == "__main__":
    unittest.main()
