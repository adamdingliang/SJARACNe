#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."

# Run from the repository root after configuring the LSF queue and resources.
# Set both input variables before running; output paths may be overridden.
: "${SJARACNE_EXPRESSION_FILE:?Set SJARACNE_EXPRESSION_FILE to the expression matrix}"
: "${SJARACNE_HUB_FILE:?Set SJARACNE_HUB_FILE to the hub ID list}"
SJARACNE_OUTPUT_DIR="${SJARACNE_OUTPUT_DIR:-./results/lsf_real_data}"
SJARACNE_TMP_DIR="${SJARACNE_TMP_DIR:-./results/tmp_lsf_real_data/}"

# Legacy affine MI cutoffs are retained; no dataset-specific PCHIP choice is made.
sjaracne lsf -j ./SJARACNe/config/config_cwlexec.json \
  -e "$SJARACNE_EXPRESSION_FILE" -g "$SJARACNE_HUB_FILE" \
  -n 100 -k 6 -pb 1e-7 --subsample-fraction 0.8 \
  -o "$SJARACNE_OUTPUT_DIR" -tmp "$SJARACNE_TMP_DIR"
