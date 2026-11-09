#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: ExpressionTool
label: Collect resample networks into a directory
doc: |
  Describe a CWL Directory containing every supplied resample graph. In the
  DPI-first workflow these graphs have already undergone per-resample DPI;
  the downstream Python consensus consumes the collected directory.

requirements:
 - class: InlineJavascriptRequirement

inputs:
  input_files:
    type: File[]
  dirname:
    type: string

outputs:
  out_dir:
    type: Directory

expression: |
  ${
    return {"out_dir": { "class": "Directory", "basename": inputs.dirname, "listing": inputs.input_files} };
  }
