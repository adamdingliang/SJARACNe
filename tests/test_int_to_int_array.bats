#!/usr/bin/env bats

load "bats_cwl_helper"

setup() {
    setup_cwl_test
}

@test "int=0 to int[]" {
    printf '{"number":0}\n' > "$CWL_JOB"
    run cwltool_to_json "$PROJECT_ROOT/SJARACNe/cwl/int_to_int_array.cwl" "$CWL_JOB"
    [ "$status" -eq 0 ]

    run assert_array_json "int_array" '[]'
    [ "$status" -eq 0 ]
}

@test "int=3 to int[1, 2, 3]" {
    printf '{"number":3}\n' > "$CWL_JOB"
    run cwltool_to_json "$PROJECT_ROOT/SJARACNe/cwl/int_to_int_array.cwl" "$CWL_JOB"
    [ "$status" -eq 0 ]

    run assert_array_json "int_array" '[1,2,3]'
    [ "$status" -eq 0 ]
}
