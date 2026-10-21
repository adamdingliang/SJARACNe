#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
export PATH="$(pwd)/SJARACNe/bin:$PATH"
mkdir -p ./results/cwl/cwlexec/ch_line_ending
cwlexec -pe PATH -c ./examples/manual_cwl/cwlexec/config.json --outdir ./results/cwl/cwlexec/ch_line_ending ./SJARACNe/cwl/ch_line_ending.cwl ./examples/manual_cwl/jobs/brca_ch_line_ending.yml
