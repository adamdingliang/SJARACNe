#!/usr/bin/env python3
"""Create matched machine-readable comparisons and summary plots."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "sjaracne-brca100-netbid-qc-v1"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


STAGES = ("baseline_12113fb", "pr66_5809183", "pr67_7633ebb")
STAGE_LABELS = {
    "baseline_12113fb": "Baseline",
    "pr66_5809183": "PR66: 80% WOR",
    "pr67_7633ebb": "PR67: matched null",
}
DRIVERS = {"tf": "BRCA100_TF.txt", "sig": "BRCA100_SIG.txt"}
PAIRS = (
    ("baseline_12113fb", "pr66_5809183"),
    ("pr66_5809183", "pr67_7633ebb"),
    ("baseline_12113fb", "pr67_7633ebb"),
)
COLORS = {
    "baseline_12113fb": "#D55E00",
    "pr66_5809183": "#009E73",
    "pr67_7633ebb": "#0072B2",
}
MAX_MI_PLOT_POINTS = 100_000
MI_PLOT_SEED = 20260817
SVG_METADATA = {
    "Creator": "SJARACNe BRCA100 NetBID2 QC benchmark",
    "Date": None,
}
GZIP_COMPRESSION = {"method": "gzip", "compresslevel": 9, "mtime": 0}
EDGE_MEMBERSHIP_PATTERNS = (
    ("baseline_only", "Baseline\nonly", (1, 0, 0)),
    ("pr66_only", "PR66\nonly", (0, 1, 0)),
    ("pr67_only", "PR67\nonly", (0, 0, 1)),
    ("baseline_pr66_only", "Baseline + PR66\nonly", (1, 1, 0)),
    ("baseline_pr67_only", "Baseline + PR67\nonly", (1, 0, 1)),
    ("pr66_pr67_only", "PR66 + PR67\nonly", (0, 1, 1)),
    ("all_three", "All three", (1, 1, 1)),
)


def safe_correlation(x: pd.Series, y: pd.Series, method: str) -> tuple[float, str]:
    paired = pd.concat([x, y], axis=1).dropna()
    if len(paired) < 3:
        return math.nan, "fewer than 3 paired observations"
    if paired.iloc[:, 0].nunique() < 2 or paired.iloc[:, 1].nunique() < 2:
        return math.nan, "zero variance"
    if method == "pearson":
        return float(stats.pearsonr(paired.iloc[:, 0], paired.iloc[:, 1]).statistic), ""
    return float(stats.spearmanr(paired.iloc[:, 0], paired.iloc[:, 1]).statistic), ""


def save_plot(figure: plt.Figure, plots_root: Path, stem: str) -> None:
    figure.savefig(plots_root / f"{stem}.png", dpi=180)
    figure.savefig(plots_root / f"{stem}.svg", metadata=SVG_METADATA)


def deterministic_pair_sample(
    x: pd.Series,
    y: pd.Series,
    *,
    maximum: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return matching, reproducibly sampled positions from two numeric series."""
    x_values = x.to_numpy(dtype=float, copy=False)
    y_values = y.to_numpy(dtype=float, copy=False)
    if x_values.shape != y_values.shape:
        raise ValueError("Paired plot vectors have different shapes")
    if len(x_values) <= maximum:
        return x_values, y_values
    # Systematic positions cover the sorted common-edge list without the O(n)
    # permutation allocation that random sampling without replacement may use.
    positions = (
        np.arange(maximum, dtype=np.int64) * len(x_values) // maximum
        + seed % len(x_values)
    ) % len(x_values)
    positions.sort()
    return x_values[positions], y_values[positions]


def membership_pattern_counts(
    baseline: set[tuple[str, str]],
    pr66: set[tuple[str, str]],
    pr67: set[tuple[str, str]],
) -> dict[str, int]:
    """Count all seven nonempty membership regions without materializing a union."""
    all_three = len(baseline & pr66 & pr67)
    baseline_pr66_total = len(baseline & pr66)
    baseline_pr67_total = len(baseline & pr67)
    pr66_pr67_total = len(pr66 & pr67)
    baseline_pr66_only = baseline_pr66_total - all_three
    baseline_pr67_only = baseline_pr67_total - all_three
    pr66_pr67_only = pr66_pr67_total - all_three
    counts = {
        "baseline_only": (
            len(baseline) - baseline_pr66_only - baseline_pr67_only - all_three
        ),
        "pr66_only": len(pr66) - baseline_pr66_only - pr66_pr67_only - all_three,
        "pr67_only": len(pr67) - baseline_pr67_only - pr66_pr67_only - all_three,
        "baseline_pr66_only": baseline_pr66_only,
        "baseline_pr67_only": baseline_pr67_only,
        "pr66_pr67_only": pr66_pr67_only,
        "all_three": all_three,
    }
    if any(value < 0 for value in counts.values()):
        raise ValueError("Negative directed-edge membership count")
    expected_union = (
        len(baseline)
        + len(pr66)
        + len(pr67)
        - baseline_pr66_total
        - baseline_pr67_total
        - pr66_pr67_total
        + all_three
    )
    if sum(counts.values()) != expected_union:
        raise ValueError("Directed-edge membership counts do not sum to the union")
    return counts


def load_network(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype={"source": str, "target": str})
    required = {
        "source", "target", "source.symbol", "target.symbol", "MI",
        "pearson", "spearman", "slope", "p-value",
    }
    if set(frame.columns) != required:
        raise ValueError(f"Unexpected columns in {path}: {list(frame.columns)}")
    if frame.empty:
        raise ValueError(f"Consensus network contains no edges: {path}")
    if frame[["source", "target"]].isna().any().any():
        raise ValueError(f"Missing endpoints in {path}")
    if frame.duplicated(["source", "target"]).any():
        raise ValueError(f"Duplicate directed edges in {path}")
    if not np.isfinite(frame["MI"]).all() or (frame["MI"] <= 0).any():
        raise ValueError(f"Invalid MI values in {path}")
    return frame


def load_support(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype={"source": str, "target": str})
    required = {
        "source", "target", "consensus_MI", "support_count",
        "support_fraction", "mean_observed_MI", "consensus_MI_roundtrip_match",
    }
    if set(frame.columns) != required or frame.empty:
        raise ValueError(f"Unexpected or empty support table {path}")
    if frame.duplicated(["source", "target"]).any():
        raise ValueError(f"Duplicate directed edges in {path}")
    if (
        not frame["support_count"].between(1, 100).all()
        or not np.allclose(
            frame["support_fraction"], frame["support_count"] / 100.0,
            rtol=0.0, atol=1e-12,
        )
        or not (frame["consensus_MI_roundtrip_match"] == 1).all()
    ):
        raise ValueError(f"Invalid support values in {path}")
    return frame


def load_and_validate_netbid_metrics(
    qc_root: Path,
    candidates: list[str],
    degree: pd.Series,
    edge_count: int,
    incident_nodes: int,
) -> tuple[dict[str, float], pd.DataFrame]:
    summary = pd.read_csv(qc_root / "network_summary.tsv", sep="\t")
    required_metrics = {
        "candidate_drivers", "active_drivers", "edges", "incident_nodes",
        "weak_components", "largest_weak_component", "density",
        "target_size_zero_mean", "target_size_zero_median",
        "target_size_zero_q25", "target_size_zero_q75", "target_size_zero_max",
        "target_size_active_mean", "target_size_active_median",
        "target_size_active_q25", "target_size_active_q75",
        "target_size_active_max", "scale_free_adjusted_r2",
    }
    if (
        set(summary.columns) != {"metric", "value"}
        or summary["metric"].duplicated().any()
        or set(summary["metric"]) != required_metrics
        or not np.isfinite(summary["value"]).all()
    ):
        raise ValueError(f"Invalid NetBID2 summary in {qc_root}")
    metrics = dict(zip(summary["metric"], summary["value"]))

    target_sizes = pd.read_csv(
        qc_root / "driver_target_sizes.tsv",
        sep="\t",
        dtype={"driver": str},
    )
    if (
        set(target_sizes.columns) != {"driver", "target_count"}
        or target_sizes["driver"].duplicated().any()
        or set(target_sizes["driver"]) != set(candidates)
    ):
        raise ValueError(f"Invalid NetBID2 target-size table in {qc_root}")
    aligned_target_sizes = target_sizes.set_index("driver").loc[candidates, "target_count"]
    if not np.array_equal(aligned_target_sizes.to_numpy(), degree.to_numpy()):
        raise ValueError(f"NetBID2 target sizes disagree with directed edges in {qc_root}")

    exact_expectations = {
        "candidate_drivers": len(candidates),
        "active_drivers": int((degree > 0).sum()),
        "edges": edge_count,
        "incident_nodes": incident_nodes,
        "target_size_zero_mean": float(degree.mean()),
        "target_size_zero_median": float(degree.median()),
        "target_size_zero_q25": float(degree.quantile(0.25)),
        "target_size_zero_q75": float(degree.quantile(0.75)),
        "target_size_zero_max": int(degree.max()),
    }
    for name, expected in exact_expectations.items():
        if not math.isclose(metrics[name], expected, rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError(
                f"NetBID2 metric {name}={metrics[name]} disagrees with "
                f"edge-derived value {expected} in {qc_root}"
            )
    return metrics, summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--work-root",
        type=Path,
        default=(
            Path.home()
            / "sjaracne-benchmarks"
            / "brca100-netbid-qc-20260817-rerun"
        ),
    )
    args = parser.parse_args()
    results_root = args.work_root / "results"
    comparison_root = results_root / "comparison"
    plots_root = comparison_root / "plots"
    comparison_root.mkdir(exist_ok=True)
    plots_root.mkdir(exist_ok=True)

    inference = pd.read_csv(results_root / "run_manifest.tsv", sep="\t")
    if (
        len(inference) != 600
        or set(zip(inference["stage"], inference["driver"]))
        != {(stage, driver) for stage in STAGES for driver in DRIVERS}
        or not (inference.groupby(["stage", "driver"]).size() == 100).all()
    ):
        raise ValueError("Inference manifest does not contain 100 runs for all six arms")
    inference_summary = (
        inference.groupby(["stage", "driver"], sort=False)
        .agg(
            runs=("seed", "size"),
            edges_median=("edges", "median"),
            edges_min=("edges", "min"),
            edges_max=("edges", "max"),
            user_s_median=("user_s", "median"),
            system_s_median=("system_s", "median"),
            elapsed_s_median_mixed_load=("elapsed_s", "median"),
            max_rss_mib_median=("max_rss_kib", lambda values: values.median() / 1024.0),
            max_rss_mib_max=("max_rss_kib", lambda values: values.max() / 1024.0),
            adjacency_mib_median=("adjacency_bytes", lambda values: values.median() / 2**20),
        )
        .reset_index()
    )
    inference_summary.to_csv(
        comparison_root / "inference_summary.tsv", sep="\t", index=False
    )

    networks: dict[tuple[str, str], pd.DataFrame] = {}
    sizes: dict[tuple[str, str], pd.Series] = {}
    summary_rows: list[dict[str, object]] = []
    size_rows: list[pd.DataFrame] = []
    netbid_summary_rows: list[pd.DataFrame] = []

    for driver, driver_filename in DRIVERS.items():
        candidates = [
            line.strip()
            for line in (args.work_root / "inputs" / driver_filename)
            .read_text(encoding="utf-8")
            .splitlines()
            if line.strip()
        ]
        if len(candidates) != len(set(candidates)):
            raise ValueError(f"Duplicate candidate IDs in {driver_filename}")
        for stage in STAGES:
            arm_root = results_root / stage / driver
            path = arm_root / "consensus" / "consensus_network_ncol_.txt"
            frame = load_network(path)
            support_path = arm_root / "consensus" / "consensus_support.tsv"
            support = load_support(support_path)
            frame = frame.merge(
                support,
                on=["source", "target"],
                how="left",
                validate="one_to_one",
                sort=False,
            )
            if frame["support_count"].isna().any() or len(frame) != len(support):
                raise ValueError(f"Support edges do not match consensus edges: {path}")
            if not np.allclose(
                frame["MI"], frame["consensus_MI"], rtol=0.0, atol=5.00001e-5
            ):
                raise ValueError(f"Support MI does not match consensus MI: {path}")
            if not set(frame["source"]).issubset(candidates):
                raise ValueError(f"Non-candidate source in {path}")
            frame = frame.set_index(["source", "target"], drop=False)
            networks[(stage, driver)] = frame
            degree = (
                frame.groupby(level="source", sort=False)
                .size()
                .reindex(candidates, fill_value=0)
            )
            degree.name = "target_count"
            sizes[(stage, driver)] = degree
            incident_nodes = len(set(frame["source"]) | set(frame["target"]))
            netbid_metrics, netbid_summary = load_and_validate_netbid_metrics(
                arm_root / "netbid2_qc",
                candidates,
                degree,
                len(frame),
                incident_nodes,
            )
            netbid_summary_rows.append(
                netbid_summary.assign(stage=stage, driver_class=driver)
            )
            size_rows.append(
                pd.DataFrame(
                    {
                        "stage": stage,
                        "driver_class": driver,
                        "driver": degree.index,
                        "target_count": degree.to_numpy(),
                    }
                )
            )
            summary_rows.append(
                {
                    "stage": stage,
                    "driver_class": driver,
                    "candidate_drivers": len(candidates),
                    "active_drivers": int((degree > 0).sum()),
                    "edges": len(frame),
                    "incident_nodes": incident_nodes,
                    "target_size_mean_zero_filled": float(degree.mean()),
                    "target_size_median_zero_filled": float(degree.median()),
                    "target_size_q25_zero_filled": float(degree.quantile(0.25)),
                    "target_size_q75_zero_filled": float(degree.quantile(0.75)),
                    "target_size_max": int(degree.max()),
                    "MI_mean": float(frame["MI"].mean()),
                    "MI_median": float(frame["MI"].median()),
                    "MI_q25": float(frame["MI"].quantile(0.25)),
                    "MI_q75": float(frame["MI"].quantile(0.75)),
                    "support_mean": float(frame["support_fraction"].mean()),
                    "support_median": float(frame["support_fraction"].median()),
                    "support_q25": float(frame["support_fraction"].quantile(0.25)),
                    "support_q75": float(frame["support_fraction"].quantile(0.75)),
                    "support_min": float(frame["support_fraction"].min()),
                    "support_max": float(frame["support_fraction"].max()),
                    "netbid2_weak_components": int(
                        netbid_metrics["weak_components"]
                    ),
                    "netbid2_largest_weak_component": int(
                        netbid_metrics["largest_weak_component"]
                    ),
                    "netbid2_density": netbid_metrics["density"],
                    "netbid2_target_size_active_mean": netbid_metrics[
                        "target_size_active_mean"
                    ],
                    "netbid2_target_size_active_median": netbid_metrics[
                        "target_size_active_median"
                    ],
                    "netbid2_target_size_active_q25": netbid_metrics[
                        "target_size_active_q25"
                    ],
                    "netbid2_target_size_active_q75": netbid_metrics[
                        "target_size_active_q75"
                    ],
                    "netbid2_target_size_active_max": int(
                        netbid_metrics["target_size_active_max"]
                    ),
                    "netbid2_scale_free_total_degree_adjusted_r2": netbid_metrics[
                        "scale_free_adjusted_r2"
                    ],
                }
            )
            edge_output = arm_root / "edge_metrics.tsv.gz"
            frame.reset_index(drop=True).assign(
                stage=stage, driver_class=driver
            ).to_csv(
                edge_output,
                sep="\t",
                index=False,
                compression=GZIP_COMPRESSION,
            )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(comparison_root / "network_summary.tsv", sep="\t", index=False)
    pd.concat(netbid_summary_rows, ignore_index=True).to_csv(
        comparison_root / "netbid2_network_summary_long.tsv", sep="\t", index=False
    )
    all_sizes = pd.concat(size_rows, ignore_index=True)
    all_sizes.to_csv(
        comparison_root / "driver_target_sizes.tsv.gz",
        sep="\t",
        index=False,
        compression=GZIP_COMPRESSION,
    )

    target_figure, target_axes = plt.subplots(
        len(DRIVERS),
        len(PAIRS),
        figsize=(14.5, 8.2),
        constrained_layout=True,
        squeeze=False,
    )
    mi_figure, mi_axes = plt.subplots(
        len(DRIVERS),
        len(PAIRS),
        figsize=(14.5, 8.2),
        constrained_layout=True,
        squeeze=False,
    )

    pairwise_rows: list[dict[str, object]] = []
    for driver_index, driver in enumerate(DRIVERS):
        for pair_index, (first, second) in enumerate(PAIRS):
            left = networks[(first, driver)]
            right = networks[(second, driver)]
            left_edges = set(left.index)
            right_edges = set(right.index)
            intersection = left_edges & right_edges
            union_size = len(left_edges) + len(right_edges) - len(intersection)
            left_sizes = sizes[(first, driver)]
            right_sizes = sizes[(second, driver)]
            active_union = (left_sizes > 0) | (right_sizes > 0)

            target_pearson, target_pearson_reason = safe_correlation(
                left_sizes, right_sizes, "pearson"
            )
            target_spearman, target_spearman_reason = safe_correlation(
                left_sizes, right_sizes, "spearman"
            )
            active_pearson, active_pearson_reason = safe_correlation(
                left_sizes[active_union], right_sizes[active_union], "pearson"
            )
            active_spearman, active_spearman_reason = safe_correlation(
                left_sizes[active_union], right_sizes[active_union], "spearman"
            )
            common_index = pd.MultiIndex.from_tuples(
                sorted(intersection), names=["source", "target"]
            )
            left_mi = left.loc[common_index, "MI"] if intersection else pd.Series(dtype=float)
            right_mi = right.loc[common_index, "MI"] if intersection else pd.Series(dtype=float)
            mi_pearson, mi_pearson_reason = safe_correlation(
                left_mi.reset_index(drop=True), right_mi.reset_index(drop=True), "pearson"
            )
            mi_spearman, mi_spearman_reason = safe_correlation(
                left_mi.reset_index(drop=True), right_mi.reset_index(drop=True), "spearman"
            )

            target_axis = target_axes[driver_index, pair_index]
            target_x = np.log1p(left_sizes.to_numpy(dtype=float, copy=False))
            target_y = np.log1p(right_sizes.to_numpy(dtype=float, copy=False))
            target_limit = max(float(target_x.max()), float(target_y.max()), 1.0)
            target_axis.scatter(
                target_x,
                target_y,
                s=7,
                alpha=0.28,
                color="#3C5488",
                edgecolors="none",
                rasterized=True,
            )
            target_axis.plot(
                [0.0, target_limit],
                [0.0, target_limit],
                linestyle="--",
                linewidth=1.0,
                color="#666666",
            )
            target_axis.set_xlim(0.0, target_limit * 1.02)
            target_axis.set_ylim(0.0, target_limit * 1.02)
            target_axis.set_aspect("equal", adjustable="box")
            target_axis.set_xlabel(f"log1p targets: {STAGE_LABELS[first]}")
            target_axis.set_ylabel(f"log1p targets: {STAGE_LABELS[second]}")
            target_axis.set_title(
                f"{driver.upper()} | all n={len(left_sizes):,}, "
                f"Spearman={target_spearman:.3f}\n"
                f"union-active n={int(active_union.sum()):,}, "
                f"Spearman={active_spearman:.3f}"
            )
            target_axis.grid(alpha=0.2)

            mi_axis = mi_axes[driver_index, pair_index]
            sampled_left_mi, sampled_right_mi = deterministic_pair_sample(
                left_mi,
                right_mi,
                maximum=MAX_MI_PLOT_POINTS,
                seed=MI_PLOT_SEED + 10 * driver_index + pair_index,
            )
            if len(sampled_left_mi):
                mi_limit = max(
                    float(left_mi.max()),
                    float(right_mi.max()),
                    1e-6,
                )
                collection = mi_axis.hexbin(
                    sampled_left_mi,
                    sampled_right_mi,
                    gridsize=55,
                    mincnt=1,
                    bins="log",
                    cmap="viridis",
                    extent=(0.0, mi_limit, 0.0, mi_limit),
                )
                collection.set_rasterized(True)
                mi_axis.plot(
                    [0.0, mi_limit],
                    [0.0, mi_limit],
                    linestyle="--",
                    linewidth=1.0,
                    color="#777777",
                )
                mi_axis.set_xlim(0.0, mi_limit)
                mi_axis.set_ylim(0.0, mi_limit)
                mi_axis.set_aspect("equal", adjustable="box")
            else:
                mi_axis.text(
                    0.5,
                    0.5,
                    "No common directed edges",
                    ha="center",
                    va="center",
                    transform=mi_axis.transAxes,
                )
            mi_axis.set_xlabel(f"Consensus MI: {STAGE_LABELS[first]}")
            mi_axis.set_ylabel(f"Consensus MI: {STAGE_LABELS[second]}")
            mi_axis.set_title(
                f"{driver.upper()} | common={len(intersection):,} | "
                f"shown={len(sampled_left_mi):,}\nSpearman={mi_spearman:.3f}"
            )
            mi_axis.grid(alpha=0.15)

            left_support = (
                left.loc[common_index, "support_fraction"]
                if intersection else pd.Series(dtype=float)
            )
            right_support = (
                right.loc[common_index, "support_fraction"]
                if intersection else pd.Series(dtype=float)
            )
            support_pearson, support_pearson_reason = safe_correlation(
                left_support.reset_index(drop=True),
                right_support.reset_index(drop=True),
                "pearson",
            )
            support_spearman, support_spearman_reason = safe_correlation(
                left_support.reset_index(drop=True),
                right_support.reset_index(drop=True),
                "spearman",
            )
            left_active = set(left["source"])
            right_active = set(right["source"])
            pairwise_rows.append(
                {
                    "driver_class": driver,
                    "first_stage": first,
                    "second_stage": second,
                    "first_edges": len(left_edges),
                    "second_edges": len(right_edges),
                    "intersection_edges": len(intersection),
                    "union_edges": union_size,
                    "edge_jaccard": len(intersection) / union_size,
                    "lost_from_first": len(left_edges) - len(intersection),
                    "gained_in_second": len(right_edges) - len(intersection),
                    "active_driver_jaccard": (
                        len(left_active & right_active) / len(left_active | right_active)
                    ),
                    "target_size_n_all": len(left_sizes),
                    "target_size_pearson_all": target_pearson,
                    "target_size_pearson_all_reason": target_pearson_reason,
                    "target_size_spearman_all": target_spearman,
                    "target_size_spearman_all_reason": target_spearman_reason,
                    "target_size_n_union_active": int(active_union.sum()),
                    "target_size_pearson_union_active": active_pearson,
                    "target_size_pearson_union_active_reason": active_pearson_reason,
                    "target_size_spearman_union_active": active_spearman,
                    "target_size_spearman_union_active_reason": active_spearman_reason,
                    "common_edge_mi_n": len(intersection),
                    "common_edge_mi_pearson": mi_pearson,
                    "common_edge_mi_pearson_reason": mi_pearson_reason,
                    "common_edge_mi_spearman": mi_spearman,
                    "common_edge_mi_spearman_reason": mi_spearman_reason,
                    "common_edge_support_n": len(intersection),
                    "common_edge_support_pearson": support_pearson,
                    "common_edge_support_pearson_reason": support_pearson_reason,
                    "common_edge_support_spearman": support_spearman,
                    "common_edge_support_spearman_reason": support_spearman_reason,
                }
            )

    del (
        left_edges,
        right_edges,
        intersection,
        common_index,
        left_mi,
        right_mi,
        left_support,
        right_support,
        sampled_left_mi,
        sampled_right_mi,
    )
    pairwise = pd.DataFrame(pairwise_rows)
    pairwise.to_csv(comparison_root / "pairwise_comparison.tsv", sep="\t", index=False)

    target_figure.suptitle(
        "Paired candidate-driver target sizes (all candidates, zeros included)"
    )
    save_plot(target_figure, plots_root, "target_size_pairwise_log1p")
    plt.close(target_figure)

    mi_figure.suptitle(
        "Common directed-edge consensus MI; hexbin color is log bin count; "
        f"at most {MAX_MI_PLOT_POINTS:,} deterministic samples per panel"
    )
    save_plot(mi_figure, plots_root, "common_edge_mi_pairwise_hexbin")
    plt.close(mi_figure)

    membership_rows: list[dict[str, object]] = []
    membership_figure, membership_axes = plt.subplots(
        1, len(DRIVERS), figsize=(13.5, 4.8), constrained_layout=True, squeeze=False
    )
    membership_colors = (
        COLORS["baseline_12113fb"],
        COLORS["pr66_5809183"],
        COLORS["pr67_7633ebb"],
        "#B07A00",
        "#8B5FBF",
        "#008C95",
        "#555555",
    )
    for driver_index, driver in enumerate(DRIVERS):
        baseline_edges = set(networks[("baseline_12113fb", driver)].index)
        pr66_edges = set(networks[("pr66_5809183", driver)].index)
        pr67_edges = set(networks[("pr67_7633ebb", driver)].index)
        counts = membership_pattern_counts(baseline_edges, pr66_edges, pr67_edges)
        values = [counts[key] for key, _, _ in EDGE_MEMBERSHIP_PATTERNS]
        labels = [label for _, label, _ in EDGE_MEMBERSHIP_PATTERNS]
        for key, _, flags in EDGE_MEMBERSHIP_PATTERNS:
            membership_rows.append(
                {
                    "driver_class": driver,
                    "membership_pattern": key,
                    "in_baseline": flags[0],
                    "in_pr66": flags[1],
                    "in_pr67": flags[2],
                    "directed_edges": counts[key],
                }
            )

        axis = membership_axes[0, driver_index]
        positions = np.arange(len(EDGE_MEMBERSHIP_PATTERNS))
        axis.bar(positions, values, color=membership_colors)
        axis.set_xticks(positions, labels, rotation=25, ha="right")
        axis.set_yscale("symlog", linthresh=1)
        axis.set_ylabel("Directed consensus edges")
        axis.set_title(f"{driver.upper()} | union={sum(values):,}")
        axis.grid(axis="y", alpha=0.25)

    pd.DataFrame(membership_rows).to_csv(
        comparison_root / "directed_edge_membership_patterns.tsv",
        sep="\t",
        index=False,
    )
    membership_figure.suptitle("Three-stage directed-edge membership patterns")
    save_plot(
        membership_figure,
        plots_root,
        "directed_edge_membership_patterns",
    )
    plt.close(membership_figure)
    del baseline_edges, pr66_edges, pr67_edges

    figure, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    for axis, driver in zip(axes, ("tf", "sig")):
        for stage in STAGES:
            values = np.sort(sizes[(stage, driver)].to_numpy())
            probabilities = np.arange(1, len(values) + 1) / len(values)
            axis.step(
                values,
                probabilities,
                where="post",
                label=STAGE_LABELS[stage],
                color=COLORS[stage],
                linewidth=1.8,
            )
        axis.set_xscale("symlog", linthresh=1)
        axis.set_xlabel("Targets per candidate driver (zero included)")
        axis.set_ylabel("Cumulative fraction")
        axis.set_title(driver.upper())
        axis.grid(alpha=0.25)
    axes[0].legend(frameon=False)
    save_plot(figure, plots_root, "target_size_ecdf")
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    for axis, driver in zip(axes, ("tf", "sig")):
        for stage in STAGES:
            values = np.sort(networks[(stage, driver)]["support_fraction"].to_numpy())
            probabilities = np.arange(1, len(values) + 1) / len(values)
            axis.step(
                values,
                probabilities,
                where="post",
                label=STAGE_LABELS[stage],
                color=COLORS[stage],
                linewidth=1.8,
            )
        axis.set_xlabel("Consensus support fraction")
        axis.set_ylabel("Cumulative fraction of retained edges")
        axis.set_title(driver.upper())
        axis.grid(alpha=0.25)
    axes[0].legend(frameon=False)
    save_plot(figure, plots_root, "consensus_support_ecdf")
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
    for axis, metric, label in (
        (axes[0], "edges", "Consensus edges"),
        (axes[1], "active_drivers", "Active drivers"),
    ):
        x = np.arange(len(STAGES))
        width = 0.36
        for offset, driver in ((-width / 2, "tf"), (width / 2, "sig")):
            values = [
                summary.loc[
                    (summary["stage"] == stage)
                    & (summary["driver_class"] == driver),
                    metric,
                ].iloc[0]
                for stage in STAGES
            ]
            axis.bar(x + offset, values, width, label=driver.upper())
        axis.set_xticks(x, [STAGE_LABELS[stage] for stage in STAGES], rotation=15)
        axis.set_ylabel(label)
        axis.grid(axis="y", alpha=0.25)
    axes[0].legend(frameon=False)
    save_plot(figure, plots_root, "network_size")
    plt.close(figure)

    print(f"Wrote comparison outputs to {comparison_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
