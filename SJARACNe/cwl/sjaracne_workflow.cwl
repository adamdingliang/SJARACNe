#!/usr/bin/env cwl-runner
# Compatibility alias for the historical DPI-first workflow.
# Prefer dpi_before_consensus.cwl; both paths load the same main process.
cwlVersion: v1.0
$graph:
  - $import: dpi_before_consensus.cwl
