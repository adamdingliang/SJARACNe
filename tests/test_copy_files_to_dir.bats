#!/usr/bin/env bats

load "bats_cwl_helper"

setup() {
    setup_cwl_test
}

@test "copy files from A to B" {
    input_dir="$BATS_TEST_TMPDIR/source files"
    mkdir -p "$input_dir"
    printf '%s\n' 'small_word' > "$input_dir/small"
    printf '%s\n' 'medium sentence, in one line' > "$input_dir/medium"
    printf '%s\n' 'large file has multiple' 'lines, and can still be copied' \
        'properly' > "$input_dir/large"
    printf '%s\n' \
        '{"dirname":"copydir","input_files":[{"class":"File","path":"source files/small"},{"class":"File","path":"source files/medium"},{"class":"File","path":"source files/large"}]}' \
        > "$CWL_JOB"

    run cwltool_to_json "$PROJECT_ROOT/SJARACNe/cwl/copy_files_to_dir.cwl" "$CWL_JOB"
    [ "$status" -eq 0 ]

    run python3 - "$CWL_RESULT" <<'PY'
import json
import sys
from pathlib import Path

actual = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
assert type(actual) is dict and set(actual) == {"out_dir"}, actual
directory = actual["out_dir"]
assert type(directory) is dict, directory
assert directory["class"] == "Directory", directory
assert directory["basename"] == "copydir", directory
listing = directory["listing"]
assert type(listing) is list, listing
assert all(type(item) is dict for item in listing), listing
assert [(item["class"], item["basename"]) for item in listing] == [
    ("File", "small"), ("File", "medium"), ("File", "large")
], listing
PY
    [ "$status" -eq 0 ]

    run cmp -- "$input_dir/small" "$BATS_TEST_TMPDIR/output/copydir/small"
    [ "$status" -eq 0 ]
    run cmp -- "$input_dir/medium" "$BATS_TEST_TMPDIR/output/copydir/medium"
    [ "$status" -eq 0 ]
    run cmp -- "$input_dir/large" "$BATS_TEST_TMPDIR/output/copydir/large"
    [ "$status" -eq 0 ]
}
