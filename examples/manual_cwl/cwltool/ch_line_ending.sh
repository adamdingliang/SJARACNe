#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
export PATH="$(pwd)/SJARACNe/bin:$PATH"
mkdir -p ./results/cwl/ch_line_ending
cwltool --outdir ./results/cwl/ch_line_ending ./SJARACNe/cwl/ch_line_ending.cwl ./examples/manual_cwl/jobs/brca_ch_line_ending.yml
