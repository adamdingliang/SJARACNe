#!/usr/bin/env cwl-runner
# Compatibility alias for the consensus and DPI tool.
# Prefer consensus_and_dpi.cwl; both paths load the same main process.
cwlVersion: v1.0
$graph:
  - $import: consensus_and_dpi.cwl
