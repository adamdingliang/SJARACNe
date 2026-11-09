#!/usr/bin/env bats

load "bats_cwl_helper"

setup() {
    setup_cwl_test
}

@test "int to string array" {
    printf '{"number":3}\n' > "$CWL_JOB"
    run cwltool_to_json "$PROJECT_ROOT/SJARACNe/cwl/int_to_str_array.cwl" "$CWL_JOB"
    [ "$status" -eq 0 ]

    run assert_array_json "str_array" '["TF_run_001.adj","TF_run_002.adj","TF_run_003.adj"]'
    [ "$status" -eq 0 ]
}
