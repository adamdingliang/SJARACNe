#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: ExpressionTool
label: Create resample adjacency filenames
doc: |
  Return the historical TF_run_*.adj output filenames in resample order.
  The DPI-first workflow pairs this array with the seed array by dotproduct;
  the filename prefix is retained for both TF and SIG network roles.

requirements:
 - class: InlineJavascriptRequirement

inputs:
  number:
    type: int
    label: a positive integer

outputs:
  str_array:
    type: string[]

expression: |
  ${ var str_arr = [], str_i = '';
     for (var i = 1; i < inputs.number+1; i++) {
       if (i <= 999) {
         str_i = ("00" + i).slice(-3);
       }
       str_arr.push('TF_run_' + str_i + '.adj');
     }
     return { "str_array": str_arr };
   }
