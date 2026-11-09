#!/usr/bin/env cwl-runner
# Compatibility alias for the resampled MI estimation tool.
# Prefer resample_and_estimate_mi.cwl; both paths load the same main process.
cwlVersion: v1.0
$graph:
  - $import: resample_and_estimate_mi.cwl
