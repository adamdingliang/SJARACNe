#!/usr/bin/env cwl-runner
# Compatibility alias for the input preparation tool.
# Prefer prepare_inputs.cwl; both paths load the same main process.
cwlVersion: v1.0
$graph:
  - $import: prepare_inputs.cwl
