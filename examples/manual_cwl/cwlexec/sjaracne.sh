#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
export PATH="$(pwd)/SJARACNe/bin:$PATH"
mkdir -p ./results/cwl/cwlexec/sjaracne
cwlexec -pe PATH -c ./examples/manual_cwl/cwlexec/config.json --outdir ./results/cwl/cwlexec/sjaracne ./SJARACNe/cwl/sjaracne.cwl ./examples/manual_cwl/jobs/brca_sjaracne.yml
