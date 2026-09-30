#!/usr/bin/env python3
"""Synthetic unit tests for the K_DPI witness screen harness."""

from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

import analyze_screen
import run_screen
from screen_common import (
    K_VALUES,
    SIDECAR_FIELDS,
    aggregate_sidecar_group,
    parse_witness_sidecar,
)


def sidecar_text(rows: list[tuple[int, ...]], *, source_count: int | None = None) -> str:
    count = len(rows) if source_count is None else source_count
    provenance = (
        "# schema\tsjaracne.dpi_witness_threshold_counts.v1\n"
        "# graph_state\tunchanged pre-DPI edges after native K=1 marking\n"
        "# count_unit\tsource-target edges\n"
        "# count_semantics\tedges having at least the indicated number of distinct eligible intermediates\n"
        "# source_index_basis\tzero-based expression-row index\n"
        "# dpi_epsilon\t0\n"
        "# source_mode\tselected rows\n"
        f"# source_count\t{count}\n"
        f"# annotated_gene_count\t{count}\n"
        "# input_file\tinput.exp\n"
        "# input_adjacency_file\t\n"
        "# network_output_file\tnetwork.adj\n"
        "# subnetwork_file\tpanel.txt\n"
        "# annotation_file\tpanel.txt\n"
        f"# k1_pruned_edges\t{sum(row[2] for row in rows)}\n"
    )
    header = "\t".join(SIDECAR_FIELDS) + "\n"
    body = "".join("\t".join(str(value) for value in row) + "\n" for row in rows)
    return provenance + header + body


class SidecarTests(unittest.TestCase):
    def write(self, root: Path, payload: str) -> Path:
        path = root / "witness.tsv"
        path.write_text(payload, encoding="utf-8", newline="\n")
        return path

    def test_valid_sidecar_and_group_aggregation(self) -> None:
        rows = [(2, 10, 6, 4, 3, 1, 0), (5, 20, 8, 6, 4, 2, 1)]
        with tempfile.TemporaryDirectory() as directory:
            parsed = parse_witness_sidecar(
                self.write(Path(directory), sidecar_text(rows)),
                expected_source_indices={2, 5},
            )
        self.assertEqual(parsed["totals"]["pre_edges"], 30)
        self.assertEqual(parsed["totals"]["witnesses_ge_1"], 14)
        common = aggregate_sidecar_group(parsed, {2}, 3)
        self.assertEqual(common["pre_edges"], 10)
        self.assertEqual(common["pruned_edges"], 3)
        self.assertAlmostEqual(common["pruned_fraction"], 0.3)

    def test_rejects_nonmonotone_counts(self) -> None:
        rows = [(2, 10, 6, 7, 3, 1, 0)]
        with tempfile.TemporaryDirectory() as directory:
            path = self.write(Path(directory), sidecar_text(rows))
            with self.assertRaisesRegex(ValueError, "Nonmonotone"):
                parse_witness_sidecar(path, expected_source_indices={2})

    def test_rejects_incomplete_source_coverage(self) -> None:
        rows = [(2, 10, 6, 4, 3, 1, 0)]
        with tempfile.TemporaryDirectory() as directory:
            path = self.write(Path(directory), sidecar_text(rows, source_count=2))
            with self.assertRaisesRegex(ValueError, "coverage mismatch"):
                parse_witness_sidecar(path, expected_source_indices={2, 5})

    def test_rejects_provenance_drift(self) -> None:
        rows = [(2, 10, 6, 4, 3, 1, 0)]
        payload = sidecar_text(rows).replace("# dpi_epsilon\t0", "# dpi_epsilon\t0.1")
        with tempfile.TemporaryDirectory() as directory:
            path = self.write(Path(directory), payload)
            with self.assertRaisesRegex(ValueError, "dpi_epsilon"):
                parse_witness_sidecar(path, expected_source_indices={2})

    def test_dynamic_path_provenance_must_match(self) -> None:
        rows = [(2, 10, 6, 4, 3, 1, 0)]
        expected = {
            "input_file": "input.exp",
            "input_adjacency_file": "",
            "network_output_file": "network.adj",
            "subnetwork_file": "panel.txt",
            "annotation_file": "panel.txt",
        }
        with tempfile.TemporaryDirectory() as directory:
            path = self.write(Path(directory), sidecar_text(rows))
            parse_witness_sidecar(
                path,
                expected_source_indices={2},
                expected_provenance=expected,
            )
            mismatched = dict(expected, network_output_file="published.adj")
            with self.assertRaisesRegex(ValueError, "network_output_file"):
                parse_witness_sidecar(
                    path,
                    expected_source_indices={2},
                    expected_provenance=mismatched,
                )


class HarnessProvenanceTests(unittest.TestCase):
    def test_frozen_harness_hashes_detect_file_drift(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "harness.py"
            path.write_text("version = 1\n", encoding="utf-8")
            paths = {"synthetic/harness.py": path}
            frozen = run_screen.harness_hashes(paths)
            self.assertEqual(
                run_screen.verify_frozen_harness_hashes(frozen, paths), frozen
            )
            path.write_text("version = 2\n", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "changed=.*synthetic/harness.py"):
                run_screen.verify_frozen_harness_hashes(frozen, paths)


class AnalysisTests(unittest.TestCase):
    def test_common_source_pre_edges_are_compared_per_row(self) -> None:
        common = {2, 5}
        records = {
            (10, 1): {2: 10, 5: 20},
            (20, 1): {2: 10, 5: 20},
            (30, 1): {2: 10, 5: 21},
            (10, 2): {2: 11, 5: 19},
            (20, 2): {2: 11, 5: 19},
            (30, 2): {2: 11, 5: 19},
        }
        result = analyze_screen.compare_common_source_pre_edges(
            records,
            common,
            hub_counts=(10, 20, 30),
            seeds=(1, 2),
        )
        self.assertEqual(result["common_pre_edge_expected_comparisons"], 8)
        self.assertEqual(result["common_pre_edge_comparisons"], 8)
        self.assertEqual(result["common_pre_edge_mismatches"], 1)
        self.assertEqual(result["common_pre_edge_missing_panel_pairs"], 0)

    def test_common_source_pre_edge_comparison_detects_missing_panel(self) -> None:
        result = analyze_screen.compare_common_source_pre_edges(
            {
                (10, 1): {2: 10},
                (20, 1): {2: 10},
            },
            {2},
            hub_counts=(10, 20, 30),
            seeds=(1,),
        )
        self.assertEqual(result["common_pre_edge_expected_comparisons"], 2)
        self.assertEqual(result["common_pre_edge_comparisons"], 1)
        self.assertEqual(result["common_pre_edge_mismatches"], 0)
        self.assertEqual(result["common_pre_edge_missing_panel_pairs"], 1)

    def test_paired_effect_metrics_are_seed_matched(self) -> None:
        metrics: list[dict] = []
        fractions = {
            1: (0.20, 0.40),
            2: (0.15, 0.25),
            3: (0.10, 0.18),
            5: (0.05, 0.09),
            10: (0.01, 0.02),
        }
        for group in analyze_screen.SOURCE_GROUPS:
            for seed in (1, 2):
                for k in K_VALUES:
                    small, full = fractions[k]
                    for hub_count, fraction in ((1335, small), (10680, full)):
                        metrics.append(
                            {
                                "source_group": group,
                                "hub_count": hub_count,
                                "seed": seed,
                                "k_dpi": k,
                                "pruned_fraction": fraction,
                                "pruning_retained_vs_k1": fraction / fractions[1][
                                    0 if hub_count == 1335 else 1
                                ],
                            }
                        )
        paired = analyze_screen.paired_effects(metrics)
        row = next(
            item
            for item in paired
            if item["source_group"] == "common_1335_sources" and item["k_dpi"] == 2
        )
        self.assertEqual(row["paired_seeds"], 2)
        self.assertAlmostEqual(row["delta_pruned_fraction_median"], 0.10)
        self.assertAlmostEqual(row["normalized_gap_vs_k1_median"], 0.50)
        self.assertEqual(row["abs_gap_reduced_vs_k1_seeds"], 2)
        self.assertEqual(row["full_gt_small_seeds"], 2)
        self.assertAlmostEqual(row["small_pruning_retained_vs_k1_median"], 0.75)
        self.assertAlmostEqual(row["full_pruning_retained_vs_k1_median"], 0.625)

        report = analyze_screen.render_report(
            [],
            [row],
            {"gates": []},
        )
        self.assertIn("small-panel pruning retained", report)
        self.assertIn("full-panel pruning retained", report)


if __name__ == "__main__":
    unittest.main()
