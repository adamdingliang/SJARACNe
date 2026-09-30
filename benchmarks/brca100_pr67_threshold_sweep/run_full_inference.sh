#!/usr/bin/env bash
set -euo pipefail

benchmark_repo=/mnt/d/GitHub/SJARACNe-brca100-netbid-qc
work_root=/home/adam/sjaracne-benchmarks/brca100-pr67-threshold-sweep-20260818
environment_runner="$benchmark_repo/benchmarks/brca100_netbid_qc/netbid2-r"
sweep_runner="$benchmark_repo/benchmarks/brca100_pr67_threshold_sweep/run_sweep.py"

exec "$environment_runner" python "$sweep_runner" \
  --phase infer \
  --points all \
  --drivers all \
  --seed-start 1 \
  --seed-end 100 \
  --workers 12 \
  --work-root "$work_root"
