#!/usr/bin/env bash
set -euo pipefail

# Run from the repository root after configuring the LSF queue and resources.
# This five-network smoke test uses K=1, not scientific-analysis settings.
# Keep legacy affine MI cutoffs and the default 80% subsampling.
sjaracne lsf -j ./SJARACNe/config/config_cwlexec.json \
  -e ./tests/inputs/BRCA100.exp -g ./tests/inputs/BRCA100_TF.txt \
  -n 5 -k 1 -pb 1e-7 --subsample-fraction 0.8 \
  -o ./results/tf_lsf_smoke -tmp ./results/tmp_tf_lsf_smoke/
