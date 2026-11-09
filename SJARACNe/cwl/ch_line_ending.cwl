#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: Run the historical line-ending check or convert to an explicit output
doc: |
  Inspect the localized expression or hub-list file. An explicit output_file
  selects the converted LF file. For compatibility, an omitted output_file
  returns the staged input filename; CRLF or CR conversion may leave a .tmp
  sidecar without replacing that returned file. The DPI-first workflow keeps
  this historical omitted-output behavior. Original user files are unchanged.

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input_file)

baseCommand: [python3, -m, SJARACNe.bin.ch_line_ending]

inputs:
  input_file:
    type: File
    inputBinding:
      position: 1
      prefix: -i
      valueFrom: $(self.basename)
  output_file:
    type: string?
    inputBinding:
      position: 2
      prefix: -o

outputs:
  out_file:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.output_file == null) {
            return inputs.input_file.basename;
          }
          else {
            return inputs.output_file;
          }
        }
