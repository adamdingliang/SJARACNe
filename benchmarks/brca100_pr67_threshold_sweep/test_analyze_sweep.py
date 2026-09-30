#!/usr/bin/env python3
"""Focused unit tests for threshold-sweep target-size reporting."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_sweep as analysis


class TargetSizeReportingTest(unittest.TestCase):
    def make_network(self) -> tuple[pd.DataFrame, list[str], set[str]]:
        candidates = ["driver10", "driver30", "driver100", "driver0"]
        records: list[dict[str, object]] = []
        expression_ids = set(candidates)
        for driver, target_count, support_count in (
            ("driver10", 10, 20),
            ("driver30", 30, 50),
            ("driver100", 100, 80),
        ):
            for index in range(target_count):
                target = f"{driver}_target_{index:03d}"
                expression_ids.add(target)
                records.append(
                    {
                        "source": driver,
                        "target": target,
                        "MI": 0.25 + index / 1000.0,
                        "support_count": support_count,
                        "support_fraction": support_count / 100.0,
                    }
                )
        frame = (
            pd.DataFrame.from_records(records)
            .set_index(["source", "target"], drop=False)
            .sort_index()
        )
        return frame, candidates, expression_ids

    def test_total_and_recurrence_qualified_driver_fractions(self) -> None:
        frame, candidates, expression_ids = self.make_network()
        row, target_sizes = analysis.summarize_network(
            frame,
            candidates,
            expression_ids,
            prefix={"p_key": "synthetic", "p_value": 1e-3},
            seed_runs=100,
            candidate_pair_tests=1000,
        )

        self.assertEqual(
            target_sizes.to_dict(),
            {"driver10": 10, "driver30": 30, "driver100": 100, "driver0": 0},
        )
        expected_total_coverage = (
            (10, 3),
            (20, 2),
            (30, 2),
            (50, 1),
            (100, 1),
        )
        for minimum, expected_count in expected_total_coverage:
            with self.subTest(minimum=minimum):
                self.assertEqual(
                    row[f"target_size_ge_{minimum}_driver_count"], expected_count
                )
                self.assertEqual(
                    row[f"target_size_ge_{minimum}_driver_fraction"],
                    expected_count / len(candidates),
                )

        # Threshold comparisons are inclusive: support 0.20, 0.50, and 0.80
        # all qualify at their corresponding boundary.
        self.assertEqual(row["support_ge_20pct_edges"], 140)
        self.assertEqual(row["support_ge_50pct_edges"], 130)
        self.assertEqual(row["support_ge_80pct_edges"], 100)
        self.assertEqual(row["support_ge_20pct_target_size_median_zero_filled"], 20)
        self.assertEqual(row["support_ge_50pct_target_size_median_zero_filled"], 15)
        self.assertEqual(row["support_ge_80pct_target_size_median_zero_filled"], 0)
        self.assertEqual(
            row["support_ge_50pct_target_size_ge_30_driver_count"], 2
        )
        self.assertEqual(
            row["support_ge_50pct_target_size_ge_30_driver_fraction"], 0.5
        )
        self.assertEqual(
            row["support_ge_80pct_target_size_ge_100_driver_count"], 1
        )
        self.assertEqual(
            row["support_ge_80pct_target_size_ge_100_driver_fraction"], 0.25
        )
        interpretation = row["recurrence_qualified_target_interpretation"]
        self.assertIn("stability", interpretation)
        self.assertIn("not biological truth", interpretation)
        self.assertIn("FDR", interpretation)

    def test_empty_consensus_zero_fills_all_recurrence_metrics(self) -> None:
        columns = ["source", "target", "MI", "support_count", "support_fraction"]
        empty = pd.DataFrame(columns=columns).set_index(
            ["source", "target"], drop=False
        )
        candidates = ["driver1", "driver2"]
        row, target_sizes = analysis.summarize_network(
            empty,
            candidates,
            set(candidates),
            prefix={"p_key": "empty", "p_value": 1e-7},
            seed_runs=100,
            candidate_pair_tests=10,
        )

        self.assertEqual(target_sizes.tolist(), [0, 0])
        for support in (20, 50, 80):
            self.assertEqual(row[f"support_ge_{support:02d}pct_edges"], 0)
            self.assertEqual(row[f"support_ge_{support:02d}pct_edge_fraction"], 0.0)
            self.assertEqual(
                row[f"support_ge_{support:02d}pct_target_size_ge_10_driver_count"],
                0,
            )

    def test_topology_screen_does_not_skip_stricter_interpolated_point(self) -> None:
        rows: list[dict[str, object]] = []
        for p_key, p_value, point_class in (
            ("p3e-04", 3e-4, "interpolation-within-accepted-range"),
            ("p5e-04", 5e-4, "exact-gate2-grid"),
        ):
            for driver in analysis.DRIVERS:
                rows.append(
                    {
                        "p_key": p_key,
                        "p_value": p_value,
                        "p_label": str(p_value),
                        "driver_class": driver,
                        "calibration_point_class": point_class,
                        "active_driver_fraction": 0.95,
                        "largest_weak_component_fraction_incident": 0.98,
                        "incident_node_fraction_expression": 0.75,
                        "candidate_pair_tests": 1000,
                        "candidate_pair_tests_interpretation": "synthetic",
                        "nominal_null_exceedances_before_DPI": p_value * 1000,
                        "nominal_null_exceedances_interpretation": "synthetic",
                    }
                )
        screen, selection = analysis.operating_point_screen(pd.DataFrame(rows))

        self.assertTrue(screen["joint_tf_sig_held_out_pass"].all())
        self.assertEqual(selection["selected_p_key"], "p3e-04")
        self.assertEqual(
            selection["selection_status"],
            "provisional-interpolated-operating-point",
        )
        self.assertIn("not a regulon-density selection rule", selection["scope"])


if __name__ == "__main__":
    unittest.main()
