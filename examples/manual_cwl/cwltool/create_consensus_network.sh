#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
export PATH="$(pwd)/SJARACNe/bin:$PATH"
mkdir -p ./results/cwl/create_consensus_network
cwltool --outdir ./results/cwl/create_consensus_network ./SJARACNe/cwl/create_consensus_network.cwl \
./examples/manual_cwl/jobs/brca_create_consensus_network.yml
