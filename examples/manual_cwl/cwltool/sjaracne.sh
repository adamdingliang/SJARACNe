#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
export PATH="$(pwd)/SJARACNe/bin:$PATH"
mkdir -p ./results/cwl/sjaracne
cwltool --outdir ./results/cwl/sjaracne ./SJARACNe/cwl/sjaracne.cwl ./examples/manual_cwl/jobs/brca_sjaracne.yml
