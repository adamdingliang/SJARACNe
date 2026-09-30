#!/usr/bin/env bash
set -euo pipefail

benchmark_repo=/mnt/d/GitHub/SJARACNe-brca100-netbid-qc
work_root=/home/adam/sjaracne-benchmarks/brca100-pr67-threshold-sweep-20260818
prior_work_root=/home/adam/sjaracne-benchmarks/brca100-netbid-qc-20260817-rerun
environment_runner="$benchmark_repo/benchmarks/brca100_netbid_qc/netbid2-r"
sweep_root="$benchmark_repo/benchmarks/brca100_pr67_threshold_sweep"

"$environment_runner" python "$sweep_root/validate_anchor_equivalence.py" \
  --sweep-work-root "$work_root" \
  --prior-work-root "$prior_work_root"

"$environment_runner" python "$sweep_root/run_sweep.py" \
  --phase consensus \
  --points all \
  --drivers all \
  --seed-start 1 \
  --seed-end 100 \
  --workers 1 \
  --work-root "$work_root"

"$environment_runner" python "$sweep_root/run_support_summaries.py" \
  --benchmark-repo "$benchmark_repo" \
  --work-root "$work_root"

"$environment_runner" python "$sweep_root/run_netbid_qc.py" \
  --points all \
  --drivers all \
  --html-points none \
  --work-root "$work_root"

"$environment_runner" python "$sweep_root/analyze_sweep.py" \
  --work-root "$work_root" \
  --pr66-work-root "$prior_work_root"
