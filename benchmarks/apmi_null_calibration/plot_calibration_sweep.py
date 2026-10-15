#!/usr/bin/env python3
"""Plot the exact-m AP-MI null-tail sweep from its tracked CSV summary."""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    return parser.parse_args()


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("calibration summary is empty")
    for row in rows:
        row["m"] = int(row["m"])
        row["accepted"] = row["accepted"] == "True"
        row["cutoff"] = float(row["cutoff_p1e-7"])
        row["stability"] = float(row["relative_stability_range"])
    return rows


def main():
    args = parse_args()
    rows = read_rows(args.summary)
    accepted = [row for row in rows if row["accepted"]]
    rejected = [row for row in rows if not row["accepted"]]

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False})
    figure, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), constrained_layout=True)

    for subset, color, marker, label in (
        (rejected, "#D55E00", "x", "Rejected diagnostic fit"),
        (accepted, "#0072B2", "o", "Accepted model"),
    ):
        axes[0].scatter([row["m"] for row in subset],
                        [row["cutoff"] for row in subset],
                        color=color, marker=marker, s=65, linewidth=2, label=label)
        axes[1].scatter([row["m"] for row in subset],
                        [100.0 * row["stability"] for row in subset],
                        color=color, marker=marker, s=65, linewidth=2, label=label)

    axes[0].set_xlabel("Exact subsample size m")
    axes[0].set_ylabel("Candidate AP-MI cutoff at p = 1e-7")
    axes[0].set_title("Extreme-tail cutoff")
    axes[1].axhline(10.0, color="#444444", linestyle="--", linewidth=1.2,
                    label="10% stability limit")
    axes[1].set_xlabel("Exact subsample size m")
    axes[1].set_ylabel("Cutoff range across tail thresholds (%)")
    axes[1].set_title("p = 1e-7 threshold stability")
    axes[1].legend(frameon=False, fontsize=9)
    for axis in axes:
        axis.grid(axis="y", color="#dddddd", linewidth=0.7)

    figure.suptitle("SJARACNe estimator-matched AP-MI null sweep (Npar = 40)")
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(str(args.output_prefix) + ".png", dpi=180)
    figure.savefig(str(args.output_prefix) + ".svg")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
