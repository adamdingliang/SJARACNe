#!/usr/bin/env python3
"""Focused tests for the one-pass consensus-recurrence aggregator."""

from __future__ import annotations

import csv
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


COMPILER = os.environ.get("CXX", "g++")


@unittest.skipUnless(shutil.which(COMPILER), f"{COMPILER} is required")
class AggregateRecurrenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.build_directory = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.build_directory.name) / "aggregate_recurrence"
        source = Path(__file__).with_name("aggregate_recurrence.cpp")
        subprocess.run(
            [
                COMPILER,
                "-O2",
                "-std=c++11",
                "-Wall",
                "-Wextra",
                "-Wpedantic",
                "-Wconversion",
                "-Wshadow",
                "-Werror",
                "-o",
                str(cls.binary),
                str(source),
            ],
            check=True,
            capture_output=True,
            text=True,
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls.build_directory.cleanup()

    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    @staticmethod
    def write_adjacencies(
        directory: Path,
        count: int = 100,
        body_for_run=None,
    ) -> None:
        directory.mkdir(parents=True)
        for run in range(1, count + 1):
            body = "X\tY\t0.5\n" if body_for_run is None else body_for_run(run)
            (directory / f"TF_run_{run}.adj").write_text(
                ">  Input file BRCA100.exp\n" + body,
                encoding="utf-8",
                newline="\n",
            )

    def invoke(self, adjacency: Path, label: str = "output") -> tuple:
        output = self.root / label
        output.mkdir()
        edges = output / "eligible_edges.tsv"
        runs = output / "run_edge_counts.tsv"
        summary = output / "aggregate_summary.tsv"
        completed = subprocess.run(
            [
                str(self.binary),
                str(adjacency),
                str(edges),
                str(runs),
                str(summary),
            ],
            capture_output=True,
            text=True,
        )
        return completed, edges, runs, summary

    @staticmethod
    def read_rows(path: Path) -> list[dict[str, str]]:
        with path.open(encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_aggregates_ordered_edges_once_and_keeps_support_six_or_more(self) -> None:
        adjacency = self.root / "adjacency"

        def body(run: int) -> str:
            rows = ["X\tY\t0.5"]
            if run <= 6:
                rows.append(f"A\tB\t{run / 10.0}")
            if run <= 5:
                rows.append("A\tC\t0.2")
            if run <= 7:
                rows.append("B\tA\t0.123456")
            return "\n".join(rows) + "\n"

        self.write_adjacencies(adjacency, body_for_run=body)
        completed, edges, runs, summary = self.invoke(adjacency)
        self.assertEqual(completed.returncode, 0, completed.stderr)

        edge_rows = self.read_rows(edges)
        self.assertEqual(
            list(edge_rows[0]),
            [
                "source",
                "target",
                "mean_observed_MI",
                "consensus_MI",
                "support_count",
                "support_fraction",
            ],
        )
        self.assertEqual(
            [(row["source"], row["target"]) for row in edge_rows],
            [("A", "B"), ("B", "A"), ("X", "Y")],
        )
        by_edge = {(row["source"], row["target"]): row for row in edge_rows}
        self.assertNotIn(("A", "C"), by_edge)
        self.assertEqual(by_edge[("A", "B")]["support_count"], "6")
        self.assertAlmostEqual(
            float(by_edge[("A", "B")]["support_fraction"]), 0.06
        )
        self.assertAlmostEqual(
            float(by_edge[("A", "B")]["mean_observed_MI"]), 0.35
        )
        self.assertEqual(by_edge[("A", "B")]["consensus_MI"], "0.3500")
        self.assertEqual(by_edge[("B", "A")]["support_count"], "7")
        self.assertEqual(by_edge[("B", "A")]["consensus_MI"], "0.1235")
        self.assertEqual(by_edge[("X", "Y")]["support_count"], "100")
        self.assertEqual(by_edge[("X", "Y")]["consensus_MI"], "0.5000")

        run_rows = self.read_rows(runs)
        self.assertEqual(len(run_rows), 100)
        self.assertEqual(
            [row["adjacency_file"] for row in run_rows[:4]],
            ["TF_run_1.adj", "TF_run_10.adj", "TF_run_100.adj", "TF_run_11.adj"],
        )
        run_counts = {
            row["adjacency_file"]: int(row["edge_count"]) for row in run_rows
        }
        self.assertEqual(run_counts["TF_run_1.adj"], 4)
        self.assertEqual(run_counts["TF_run_6.adj"], 3)
        self.assertEqual(run_counts["TF_run_7.adj"], 2)
        self.assertEqual(run_counts["TF_run_8.adj"], 1)

        summary_values = {
            row["metric"]: row["value"] for row in self.read_rows(summary)
        }
        self.assertEqual(
            summary_values,
            {
                "bootstrap_runs": "100",
                "minimum_support": "6",
                "union_edges": "4",
                "retained_edges": "3",
            },
        )

    def test_outputs_are_byte_deterministic(self) -> None:
        adjacency = self.root / "adjacency"
        self.write_adjacencies(adjacency)
        first = self.invoke(adjacency, "first")
        second = self.invoke(adjacency, "second")
        self.assertEqual(first[0].returncode, 0, first[0].stderr)
        self.assertEqual(second[0].returncode, 0, second[0].stderr)
        for first_path, second_path in zip(first[1:], second[1:]):
            self.assertEqual(first_path.read_bytes(), second_path.read_bytes())

    def test_duplicate_ordered_edge_within_one_run_fails_closed(self) -> None:
        adjacency = self.root / "adjacency"

        def body(run: int) -> str:
            if run == 1:
                return "A\tB\t0.5\tB\t0.5\n"
            return "A\tB\t0.5\n"

        self.write_adjacencies(adjacency, body_for_run=body)
        completed, edges, runs, summary = self.invoke(adjacency)
        self.assertEqual(completed.returncode, 1)
        self.assertIn(
            "Duplicate ordered edge A----B in adjacency file TF_run_1.adj",
            completed.stderr,
        )
        self.assertFalse(edges.exists())
        self.assertFalse(runs.exists())
        self.assertFalse(summary.exists())

    def test_nonfinite_mi_fails_closed(self) -> None:
        adjacency = self.root / "adjacency"

        def body(run: int) -> str:
            return "A\tB\tnan\n" if run == 1 else "A\tB\t0.5\n"

        self.write_adjacencies(adjacency, body_for_run=body)
        completed, edges, runs, summary = self.invoke(adjacency)
        self.assertEqual(completed.returncode, 1)
        self.assertIn("Invalid numeric value", completed.stderr)
        self.assertFalse(edges.exists())
        self.assertFalse(runs.exists())
        self.assertFalse(summary.exists())

    def test_requires_exactly_one_hundred_adjacency_files(self) -> None:
        for count in (99, 101):
            with self.subTest(count=count):
                adjacency = self.root / f"adjacency_{count}"
                self.write_adjacencies(adjacency, count=count)
                completed, edges, runs, summary = self.invoke(
                    adjacency, f"output_{count}"
                )
                self.assertEqual(completed.returncode, 1)
                self.assertIn(
                    f"Expected exactly 100 .adj files, got {count}",
                    completed.stderr,
                )
                self.assertFalse(edges.exists())
                self.assertFalse(runs.exists())
                self.assertFalse(summary.exists())


if __name__ == "__main__":
    unittest.main()
