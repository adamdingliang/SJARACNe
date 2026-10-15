# BRCA100 hub-size DPI pilot: compact evidence package

This directory is the compact, reviewable evidence package for the completed
BRCA100 hub-size DPI pilot. The complete run remains at:

`/home/adam/sjaracne-benchmarks/brca100-hub-size-dpi-pilot-20260825`

The run used source commit `32fe12c168ef80291e487dbf4045f430b9c5d90a`,
frozen source-tree fingerprint
`3be6f0c75cf905f983d789410542adda388d07bc0231b4439371e0dd3f22169d`,
and binary SHA-256
`836c8edbfc75b517cefdde6426755f9e3b71021df09b06c6f3b9c262057cd47e`.
See `pilot_design.json` and `invocations.json` for the complete parameters and
provenance.

## Result

The descriptive trend is clear in this one matched panel: the median fraction
of pre-DPI edges pruned rose with the hub-list size in both networks.

| Network | 12.5% hubs | 50% hubs | 100% hubs | Full minus 12.5% |
|---|---:|---:|---:|---:|
| TF | 0.1079 | 0.1917 | 0.2420 | +0.1341 |
| SIG | 0.2249 | 0.3432 | 0.4166 | +0.1917 |

The paired result is stronger than the arm medians alone: every adjacent-size
pruning-fraction difference was positive in all 100 matched seeds. The median
paired increases were TF 326->1,304: +8.243 percentage points; TF
1,304->2,608: +5.115 points; SIG 1,335->5,340: +11.749 points; and SIG
5,340->10,680: +7.555 points.

This supports the narrow conclusion that hub-list size materially affects DPI
pruning under the current workflow. It passes the pilot gate for a fixed-`-l`
control and replicated panels if a causal DPI-mechanism claim is needed. It
does **not** identify a biologically optimal hub count.

All seven recorded validation gates passed. The full-size K>=6 edge counts
exactly reproduced the independent recurrence anchors: TF 416,408 and SIG
739,958. In an additional post-hoc check, the complete per-edge recurrence
support tables for the full-size TF and SIG arms were compared with the prior
independent artifacts; both comparisons were exact. That extra comparison was
performed after `validation_gates.json` was generated, so it is documented
here rather than represented as an eighth JSON gate.

## Limitations

- One deterministic panel confounds hub count with panel membership.
- The workflow-faithful arms change `-s` and DPI annotation `-l` together, so
  this pilot does not isolate the DPI mechanism.
- Edge topology and pruning fractions do not establish biological accuracy.
- The full BRCA100 target universe was held fixed; these results do not
  validate Xenium 5k.
- K>=6 is a benchmark-only direct support-count calculation, not the production
  minimum-recurrence implementation.

## Runtime and storage

- End-to-end wall time, from `prepare` start through `aggregate` finish:
  1:34:47.9 (5,687.9 seconds).
- Main 600-job inference invocation: 1:27:06.0 wall time with 12 workers; 598
  jobs were new and the two full-size seed-1 smoke jobs were resumed.
- Sum of the 600 per-job elapsed times: 61,286.99 seconds (17:01:26.99).
- Aggregation wall time: 2:23.98.
- Maximum recorded per-job RSS: 78,664 KiB.
- Complete WSL run tree: 1,680,003,537 bytes by `du -sb` (`du -sh`: 1.6G).
  Its `results/` subtree accounts for 1,616,875,798 bytes.
- Compact package: 475105 bytes across 15 files, including this README and its
  checksum manifest.

Timestamps in `invocations.json` are UTC; the local run date was 2026-08-25.

## Package contents

Included:

- `REPORT.md`: human-readable result table and validation summary.
- `dpi_summary.tsv`: six-arm DPI and provisional K>=6 summary.
- `validation_gates.json`: machine-readable gate results and limitations.
- `pilot_design.json`, `panel_manifest.json`, and `invocations.json`: frozen
  design, panel construction, and execution provenance.
- `run_manifest.tsv`: all 600 per-run DPI counts, hashes, resource metrics, and
  full-size anchor checks.
- `provisional_k6/*_manifest.json`: six per-arm K>=6 manifests.
- `MANIFEST.sha256`: SHA-256 for every package file except itself.

Deliberately excluded to keep this package compact:

- the 600 adjacency files and their logs;
- the six full provisional K>=6 consensus/support edge tables;
- raw BRCA100 input data, panel-list payloads, null model, binary, and build
  directory.

The WSL source tree remains the authoritative location for excluded artifacts.
The copied evidence files were checked byte-for-byte by comparing their SHA-256
hashes with the WSL originals. Run `sha256sum -c MANIFEST.sha256` from this
directory to verify the package.
