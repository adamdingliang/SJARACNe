#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

import analyze_pilot
from pilot_common import (
    DRIVERS,
    arm_key,
    create_panel_files,
    parse_dpi_stats,
    parse_sampling_indices,
    sha256_bytes,
    sha256_file,
    variance_balanced_order,
)
import run_pilot


THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parents[1]


class PanelTests(unittest.TestCase):
    def test_variance_balanced_order_is_deterministic_and_prefix_balanced(self) -> None:
        ids = [f"ID{i:02d}" for i in range(23)]
        variances = {accession: float(index) for index, accession in enumerate(ids)}
        first, membership = variance_balanced_order(ids, variances, driver_key="toy")
        second, _ = variance_balanced_order(ids, variances, driver_key="toy")
        self.assertEqual(first, second)
        self.assertEqual(set(first), set(ids))
        for size in (5, 10, 15, 20):
            counts = [0] * 5
            for accession in first[:size]:
                counts[membership[accession]] += 1
            self.assertLessEqual(max(counts) - min(counts), 1)

    def test_real_brca100_panels_are_nested_balanced_and_full_exact(self) -> None:
        input_root = REPO_ROOT / "tests" / "inputs"
        with tempfile.TemporaryDirectory() as temporary:
            panel_root = Path(temporary) / "panels"
            manifest = create_panel_files(input_root, panel_root)
            self.assertEqual(set(manifest["drivers"]), {"tf", "sig"})
            for driver in DRIVERS:
                records = manifest["drivers"][driver.key]["panels"]
                previous: set[str] = set()
                for expected_count, record in zip(driver.counts, records):
                    path = panel_root / record["path"]
                    selected = set(path.read_text(encoding="utf-8").splitlines())
                    self.assertEqual(len(selected), expected_count)
                    self.assertTrue(previous.issubset(selected))
                    self.assertLessEqual(
                        max(record["variance_quintile_counts"])
                        - min(record["variance_quintile_counts"]),
                        1,
                    )
                    previous = selected
                full = panel_root / records[-1]["path"]
                self.assertEqual(full.read_bytes(), (input_root / driver.filename).read_bytes())


class DpiStatsTests(unittest.TestCase):
    def write(self, root: Path, text: str) -> Path:
        path = root / "stdout.log"
        path.write_text(text, encoding="utf-8", newline="\n")
        return path

    def test_exact_record(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            indices = " ".join(str(index) for index in range(80))
            stats = parse_dpi_stats(
                self.write(
                    Path(temporary),
                    "[SAMPLING] Selected original observation indices (0-based): "
                    + indices
                    + "\nnoise\n[DPI_STATS] pre_edges=100 pruned_edges=40 "
                    "post_edges=60 dpi_applied=1\n",
                )
            )
            sampling = parse_sampling_indices(Path(temporary) / "stdout.log")
        self.assertEqual(stats["pruned_edges"], 40)
        self.assertAlmostEqual(stats["pruned_fraction"], 0.4)
        self.assertEqual(sampling["indices"], list(range(80)))

    def test_rejects_malformed_duplicate_and_unbalanced_records(self) -> None:
        cases = (
            "[DPI_STATS] pre_edges=100 pruned_edges=40 post_edges=60\n",
            "[DPI_STATS] pre_edges=100 pruned_edges=40 post_edges=60 dpi_applied=1\n"
            "[DPI_STATS] pre_edges=100 pruned_edges=40 post_edges=60 dpi_applied=1\n",
            "[DPI_STATS] pre_edges=100 pruned_edges=41 post_edges=60 dpi_applied=1\n",
            "[DPI_STATS] pre_edges=100 pruned_edges=40 post_edges=60 dpi_applied=0\n",
        )
        for text in cases:
            with self.subTest(text=text), tempfile.TemporaryDirectory() as temporary:
                with self.assertRaises(ValueError):
                    parse_dpi_stats(self.write(Path(temporary), text))


class AggregateTests(unittest.TestCase):
    def test_disk_backed_k6_counts_direct_recurrence(self) -> None:
        driver = DRIVERS[0]
        count = driver.counts[0]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            panel = root / "panels" / driver.key / f"h{count:05d}" / driver.filename
            panel.parent.mkdir(parents=True)
            ids = [f"H{i:03d}" for i in range(count)]
            panel.write_text("\n".join(ids) + "\n", encoding="utf-8", newline="\n")
            adjacency_root = root / "results" / arm_key(driver, count) / "adjacency"
            adjacency_root.mkdir(parents=True)
            for seed in range(1, 101):
                rows = [f"{ids[0]}\tTARGET_A\t0.5"]
                if seed <= 6:
                    rows.append(f"{ids[1]}\tTARGET_B\t0.4")
                if seed <= 5:
                    rows.append(f"{ids[2]}\tTARGET_C\t0.3")
                (adjacency_root / f"TF_run_{seed:03d}.adj").write_text(
                    "> test\n" + "\n".join(rows) + "\n",
                    encoding="utf-8",
                    newline="\n",
                )
            resumed, record = run_pilot.aggregate_k6_arm(root, driver, count)
            self.assertFalse(resumed)
            self.assertEqual(record["k6_edges"], 2)
            self.assertEqual(record["active_hubs"], 2)
            self.assertEqual(record["zero_filled_median_target_count"], 0.0)
            resumed, repeated = run_pilot.aggregate_k6_arm(root, driver, count)
            self.assertTrue(resumed)
            self.assertEqual(repeated, record)

    def test_analysis_gates_and_summary_on_complete_synthetic_run(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            results = root / "results"
            results.mkdir()
            panels = root / "panels"
            panels.mkdir()
            panel_manifest = panels / "panel_manifest.json"
            panel_manifest.write_text("{}\n", encoding="utf-8", newline="\n")
            design = {
                "source_commit": "a" * 40,
                "build": {"binary_sha256": "b" * 64},
                "panel_manifest_sha256": sha256_file(panel_manifest),
            }
            (root / "pilot_design.json").write_text(
                json.dumps(design) + "\n", encoding="utf-8", newline="\n"
            )
            fields = [
                "arm",
                "driver",
                "hub_count",
                "seed",
                "source_commit",
                "binary_sha256",
                "pre_edges",
                "pruned_edges",
                "post_edges",
                "pruned_fraction",
                "sampling_indices",
                "sampling_sha256",
                "adjacency_sha256",
                "data_sha256",
                "anchor_data_match",
            ]
            with (results / "run_manifest.tsv").open(
                "w", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
                writer.writeheader()
                for driver in DRIVERS:
                    for count in driver.counts:
                        for seed in range(1, 101):
                            pre = count * 10 + seed
                            pruned = count + seed
                            post = pre - pruned
                            writer.writerow(
                                {
                                    "arm": arm_key(driver, count),
                                    "driver": driver.key,
                                    "hub_count": count,
                                    "seed": seed,
                                    "source_commit": "a" * 40,
                                    "binary_sha256": "b" * 64,
                                    "pre_edges": pre,
                                    "pruned_edges": pruned,
                                    "post_edges": post,
                                    "pruned_fraction": pruned / pre,
                                    "sampling_indices": ",".join(
                                        str(index) for index in range(80)
                                    ),
                                    "sampling_sha256": sha256_bytes(
                                        (" ".join(str(index) for index in range(80)) + "\n").encode(
                                            "ascii"
                                        )
                                    ),
                                    "adjacency_sha256": "c" * 64,
                                    "data_sha256": "d" * 64,
                                    "anchor_data_match": (
                                        "True" if count == driver.full_count else ""
                                    ),
                                }
                            )
                        output_root = (
                            results / arm_key(driver, count) / "provisional_k6"
                        )
                        output_root.mkdir(parents=True)
                        output = output_root / "consensus_support_ge6.tsv"
                        output.write_text(
                            "source\ttarget\tsupport_count\tsupport_fraction\n",
                            encoding="utf-8",
                            newline="\n",
                        )
                        expected_full_edges = {"tf": 416408, "sig": 739958}
                        manifest = {
                            "k": 6,
                            "seed_count": 100,
                            "hub_count": count,
                            "k6_edges": (
                                expected_full_edges[driver.key]
                                if count == driver.full_count
                                else 0
                            ),
                            "zero_filled_median_target_count": 0.0,
                            "active_hubs": 0,
                            "active_hub_fraction": 0.0,
                            "output_sha256": sha256_file(output),
                        }
                        (output_root / "manifest.json").write_text(
                            json.dumps(manifest) + "\n", encoding="utf-8", newline="\n"
                        )
            summaries, validation = analyze_pilot.analyze(root, require_anchor=True)
            self.assertEqual(len(summaries), 6)
            self.assertTrue(validation["all_required_gates_pass"])
            self.assertTrue(all(item["status"] == "pass" for item in validation["gates"]))


if __name__ == "__main__":
    unittest.main()
