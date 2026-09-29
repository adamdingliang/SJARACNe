#!/usr/bin/env python3
"""Opt-in DREAM5 Network 1 benchmark for PR #70 sampling (stdlib only)."""

import argparse
import csv
import hashlib
import io
import json
import math
import subprocess
import sys
import zipfile
from pathlib import Path


ROOT = "DREAM5_network_inference_challenge/Network1/"
FILES = {
    "expression": "input data/net1_expression_data.tsv",
    "tfs": "input data/net1_transcription_factors.tsv",
    "gold": "gold standard/DREAM5_NetworkInference_GoldStandard_Network1.tsv",
}
SHA256 = {
    "expression": "fb39ea0bfbede740e28aef21df768dd14c1cb5063d6ed8a384aec1e2fc417663",
    "tfs": "ce9ab43c396fbdb0a02091127b895b71689295bff8bc167f25a1b55d92b77478",
    "gold": "2b5b808846e553e229723cb50f0be324f0250e722965665978aeda6c7fc12703",
}
ARCHIVE_MD5 = "fabfb9cbba5761d29d681e752a45df75"


def digest(data, algorithm="sha256"):
    return hashlib.new(algorithm, data).hexdigest()


def load_inputs(archive=None, data_root=None):
    """Read only three exact members; never extract a ZIP onto the filesystem."""
    if (archive is None) == (data_root is None):
        raise ValueError("Supply exactly one of --archive or --data-root")
    if archive is not None:
        archive = Path(archive)
        if digest(archive.read_bytes(), "md5") != ARCHIVE_MD5:
            raise ValueError("DREAM5 archive MD5 mismatch")
        with zipfile.ZipFile(archive) as zf:
            raw = {key: zf.read(ROOT + name) for key, name in FILES.items()}
    else:
        root = Path(data_root)
        raw = {key: (root / name).read_bytes() for key, name in FILES.items()}
    for key, data in raw.items():
        if digest(data) != SHA256[key]:
            raise ValueError("DREAM5 {} SHA256 mismatch".format(key))
    return {key: data.decode("utf-8") for key, data in raw.items()}


def gold_labels(text):
    labels = {}
    for row in csv.reader(io.StringIO(text), delimiter="\t"):
        if len(row) != 3 or row[2] not in ("0", "1"):
            raise ValueError("Malformed gold row: {}".format(row))
        pair = (row[0], row[1])
        if pair in labels:
            raise ValueError("Duplicate gold pair: {}".format(pair))
        labels[pair] = int(row[2])
    if not labels:
        raise ValueError("Empty gold standard")
    return labels


def prepare(inputs, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    exp = out / "dream5_net1.exp"
    hubs = out / "all195_tfs.txt"
    if exp.exists() or hubs.exists():
        raise FileExistsError("Prepared input already exists; choose a fresh --out")
    reader = csv.reader(io.StringIO(inputs["expression"]), delimiter="\t")
    genes = next(reader)
    rows = list(reader)
    tfs = inputs["tfs"].splitlines()
    labels = gold_labels(inputs["gold"])
    if (len(genes), len(rows), len(tfs), len(labels), sum(labels.values())) != (1643, 805, 195, 278392, 4012):
        raise ValueError("Unexpected DREAM5 Network 1 dimensions")
    if len(set(genes)) != len(genes) or len(set(tfs)) != len(tfs) or not set(tfs) <= set(genes):
        raise ValueError("Duplicate or unknown gene/TF")
    if any(len(row) != len(genes) for row in rows):
        raise ValueError("Ragged expression matrix")
    with exp.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["isoformId", "geneSymbol"] + ["S{}".format(i + 1) for i in range(len(rows))])
        for gene, values in zip(genes, zip(*rows)):
            writer.writerow([gene, gene] + list(values))
    hubs.write_text("\n".join(tfs) + "\n", encoding="utf-8")
    meta = {
        "source_sha256": SHA256,
        "prepared_expression_sha256": digest(exp.read_bytes()),
        "profiles": len(rows), "genes": len(genes), "tfs": len(tfs),
        "labeled_pairs": len(labels), "gold_positives": sum(labels.values()),
        "design": "all supplied TFs, all genes and profiles; gold used only for scoring",
    }
    (out / "provenance.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    return meta


def adjacency(path):
    edges = {}
    with Path(path).open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith(">") or not line.strip():
                continue
            cells = line.rstrip("\r\n").split("\t")
            if not cells[0] or len(cells) % 2 != 1:
                raise ValueError("Malformed adjacency row in {}".format(path))
            for i in range(1, len(cells), 2):
                pair = (cells[0], cells[i])
                if pair in edges:
                    raise ValueError("Duplicate adjacency pair: {}".format(pair))
                score = float(cells[i + 1])
                if not math.isfinite(score) or score < 0:
                    raise ValueError("Invalid MI for {}".format(pair))
                edges[pair] = score
    return edges


def ranking_metrics(scores, labels):
    """Tie-aware Mann-Whitney AUROC and grouped-threshold average precision."""
    if len(scores) != len(labels) or not scores:
        raise ValueError("Mismatched or empty ranking")
    positives = sum(labels)
    negatives = len(labels) - positives
    if not positives or not negatives:
        raise ValueError("Both classes are required")
    groups = {}
    for score, label in zip(scores, labels):
        if not math.isfinite(score):
            raise ValueError("Nonfinite score")
        counts = groups.setdefault(score, [0, 0])
        counts[label] += 1
    negative_below = 0
    favorable = 0.0
    for score in sorted(groups):
        neg, pos = groups[score]
        favorable += pos * (negative_below + neg / 2)
        negative_below += neg
    tp = fp = 0
    ap = 0.0
    for score in sorted(groups, reverse=True):
        neg, pos = groups[score]
        tp += pos
        fp += neg
        ap += (pos / positives) * (tp / (tp + fp))
    return favorable / (positives * negatives), ap


def score(edges_raw, edges_final, labels, hubs):
    hubs = set(hubs)
    for kind, edges in (("raw", edges_raw), ("final", edges_final)):
        if any(source not in hubs for source, _ in edges):
            raise ValueError("{} adjacency has an unselected source".format(kind))
    if any(edges_raw.get(pair) != mi for pair, mi in edges_final.items()):
        raise ValueError("Final MI disagrees with unthresholded MI")
    # The released gold file is the evaluation universe; omitted output pairs
    # outside it are NOT negatives. A missing labeled pair receives score zero.
    pairs = list(labels)
    gold = list(labels.values())
    p = sum(gold)
    n = len(gold) - p
    called = [pair for pair in edges_final if pair in labels]
    tp = sum(labels[pair] for pair in called)
    fp = len(called) - tp
    raw_auc, raw_ap = ranking_metrics([edges_raw.get(pair, 0) for pair in pairs], gold)
    final_auc, final_ap = ranking_metrics([edges_final.get(pair, 0) for pair in pairs], gold)
    return {
        "pre_dpi_edges_all": len(edges_raw), "final_edges_all": len(edges_final),
        "final_edges_labeled": len(called), "final_edges_unlabeled": len(edges_final) - len(called),
        "tp": tp, "fp": fp, "fn": p - tp, "tn": n - fp,
        "precision": tp / len(called) if called else None,
        "sensitivity": tp / p, "specificity": (n - fp) / n,
        "auroc_all_unthresholded_mi": raw_auc,
        "average_precision_all_unthresholded_mi": raw_ap,
        "auroc_final_mi": final_auc,
        "average_precision_final_mi": final_ap,
        "auroc_final_binary": (tp / p + (n - fp) / n) / 2,
    }


def run(args):
    out = args.out
    exp, hubs = out / "dream5_net1.exp", out / "all195_tfs.txt"
    if not exp.is_file() or not hubs.is_file():
        raise FileNotFoundError("Run prepare first")
    if not args.binary.is_file() or not args.config.is_dir():
        raise FileNotFoundError("Supply existing --binary and --config")
    runs = out / "runs"
    runs.mkdir(parents=True, exist_ok=True)
    arms = ["legacy"] + ["u{}".format(f) for f in args.fractions]
    for seed in args.seeds:
        for arm in arms:
            sampler = ["-r", "1"] if arm == "legacy" else ["-u", arm[1:] + "%"]
            # p=1/t=0 gives the entire pre-DPI MI ranking. Nominal p gives
            # the final post-DPI operating point; never confuse the two.
            stages = ["all_mi", "nominal_raw", "final"]
            if args.fixed_mi is not None:
                stages += ["fixed_raw", "fixed_final"]
            for stage in stages:
                target = runs / "{}_seed{}_{}.adj".format(arm, seed, stage)
                log = target.with_suffix(".log")
                if target.exists() or log.exists():
                    raise FileExistsError("Refusing to overwrite {}".format(target))
                if stage == "all_mi":
                    cutoff = ["-p", "1", "-t", "0", "-e", "1"]
                elif stage.startswith("fixed_"):
                    cutoff = ["-p", "1", "-t", str(args.fixed_mi),
                              "-e", "1" if stage == "fixed_raw" else "0"]
                else:
                    cutoff = ["-p", str(args.pvalue),
                              "-e", "1" if stage == "nominal_raw" else "0"]
                cmd = [str(args.binary.resolve()), "-i", str(exp.resolve()), "-l", str(hubs.resolve()), "-s", str(hubs.resolve()),
                       *cutoff, "-N", str(args.npar), "-H", str(args.config.resolve()), "-S", str(seed), *sampler,
                       "-o", str(target.resolve())]
                with log.open("w", encoding="utf-8") as stream:
                    result = subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT, check=False)
                if result.returncode:
                    raise RuntimeError("Command failed ({}); see {}".format(result.returncode, log))
                print(target)
    manifest = {"binary_sha256": digest(args.binary.read_bytes()), "pvalue": args.pvalue,
                "fixed_mi": args.fixed_mi,
                "npar": args.npar, "seeds": args.seeds, "arms": arms,
                "note": "legacy uses PR70 compatibility path -r 1; u arms use -u N%"}
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def score_runs(args, inputs):
    out = args.out
    labels = gold_labels(inputs["gold"])
    hubs = inputs["tfs"].splitlines()
    manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    results = []
    for seed in manifest["seeds"]:
        for arm in manifest["arms"]:
            stem = "{}_seed{}_".format(arm, seed)
            raw = adjacency(out / "runs" / (stem + "all_mi.adj"))
            nominal_raw = adjacency(out / "runs" / (stem + "nominal_raw.adj"))
            final = adjacency(out / "runs" / (stem + "final.adj"))
            if any(raw.get(pair) != mi for pair, mi in nominal_raw.items()):
                raise ValueError("Nominal-p raw MI disagrees with unthresholded MI")
            if any(nominal_raw.get(pair) != mi for pair, mi in final.items()):
                raise ValueError("Final edge missing from nominal-p pre-DPI network")
            row = {"seed": seed, "arm": arm, **score(raw, final, labels, hubs)}
            row["nominal_pre_dpi_edges_all"] = len(nominal_raw)
            row["nominal_dpi_pruning_fraction"] = (
                1 - len(final) / len(nominal_raw) if nominal_raw else None)
            if manifest["fixed_mi"] is not None:
                fixed_raw = adjacency(out / "runs" / (stem + "fixed_raw.adj"))
                fixed_final = adjacency(out / "runs" / (stem + "fixed_final.adj"))
                if any(raw.get(pair) != mi for pair, mi in fixed_raw.items()):
                    raise ValueError("Fixed-cutoff MI disagrees with unthresholded MI")
                if any(fixed_raw.get(pair) != mi for pair, mi in fixed_final.items()):
                    raise ValueError("Fixed final edge missing from fixed pre-DPI network")
                fixed = score(raw, fixed_final, labels, hubs)
                row["fixed_mi"] = {
                    "threshold": manifest["fixed_mi"],
                    "pre_dpi_edges_all": len(fixed_raw),
                    "final_edges_all": len(fixed_final),
                    "precision": fixed["precision"],
                    "sensitivity": fixed["sensitivity"],
                    "specificity": fixed["specificity"],
                }
            results.append(row)
    report = {"universe": {"pairs": len(labels), "positives": sum(labels.values()),
                            "negatives": len(labels) - sum(labels.values())},
              "definitions": {"unthresholded_mi": "pre-DPI -p 1 -t 0 -e 1; absent labeled pairs score zero",
                              "nominal_pre_dpi": "-p nominal -e 1; DPI pruning fraction uses all output edges",
                              "final": "post-DPI nominal-p edges; absent labeled pairs score zero",
                              "average_precision": "stepwise PR area, grouping tied scores",
                              "auroc": "Mann-Whitney, half credit for positive-negative ties"},
              "per_seed": results}
    (out / "metrics.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for row in results:
        print("{arm} seed {seed}: P={precision:.4f} recall={sensitivity:.4f} specificity={specificity:.4f} "
              "raw AUROC={auroc_all_unthresholded_mi:.4f} raw AP={average_precision_all_unthresholded_mi:.4f}".format(**row))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "score"))
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--archive", type=Path, help="Original Zenodo challenge ZIP")
    source.add_argument("--data-root", type=Path, help="Extracted Network1 directory")
    parser.add_argument("--out", type=Path, required=True, help="Ignored working/output directory")
    parser.add_argument("--binary", type=Path, help="PR70 sjaracne executable")
    parser.add_argument("--config", type=Path, help="Matching SJARACNe/config directory")
    parser.add_argument("--fractions", type=int, nargs="+", default=[64, 80, 90])
    parser.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    parser.add_argument("--pvalue", type=float, default=1e-7)
    parser.add_argument("--fixed-mi", type=float, help="Optional common explicit MI threshold control")
    parser.add_argument("--npar", type=int, default=40)
    args = parser.parse_args()
    if args.action in ("prepare", "score") and not (args.archive or args.data_root):
        parser.error("prepare/score requires --archive or --data-root")
    if args.action == "run" and (not args.binary or not args.config):
        parser.error("run requires --binary and --config")
    if any(not 1 <= f <= 100 for f in args.fractions) or any(seed < 0 for seed in args.seeds):
        parser.error("invalid fractions or seeds")
    if not 0 < args.pvalue <= 1 or args.npar < 1:
        parser.error("invalid p-value or npar")
    if args.fixed_mi is not None and (not math.isfinite(args.fixed_mi) or args.fixed_mi < 0):
        parser.error("invalid fixed MI threshold")
    if args.action == "prepare":
        print(json.dumps(prepare(load_inputs(args.archive, args.data_root), args.out), indent=2))
    elif args.action == "run":
        run(args)
    else:
        score_runs(args, load_inputs(args.archive, args.data_root))


if __name__ == "__main__":
    main()
