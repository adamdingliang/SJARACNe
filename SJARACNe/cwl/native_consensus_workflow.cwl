#!/usr/bin/env cwl-runner
# Compatibility alias for the consensus-first workflow.
# Prefer consensus_before_dpi.cwl; both paths load the same main process.
cwlVersion: v1.0
$graph:
  - $import: consensus_before_dpi.cwl
