setup_cwl_test() {
    PROJECT_ROOT="$(cd "$BATS_TEST_DIRNAME/.." && pwd)"
    CWL_JOB="$BATS_TEST_TMPDIR/job.json"
    CWL_RESULT="$BATS_TEST_TMPDIR/result.json"
    cd "$BATS_TEST_TMPDIR"
}

cwltool_to_json() {
    cwltool --outdir "$BATS_TEST_TMPDIR/output" \
        --tmpdir-prefix "$BATS_TEST_TMPDIR/cwl-tmp-" \
        --tmp-outdir-prefix "$BATS_TEST_TMPDIR/cwl-out-" \
        "$@" > "$CWL_RESULT"
}

assert_array_json() {
    python3 - "$CWL_RESULT" "$1" "$2" <<'PY'
import json
import sys
from pathlib import Path

actual = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
key = sys.argv[2]
expected = json.loads(sys.argv[3])
assert type(actual) is dict, actual
assert set(actual) == {key}, actual
assert type(actual[key]) is list, actual
assert actual[key] == expected, actual
assert all(type(item) is type(reference)
           for item, reference in zip(actual[key], expected)), actual
PY
}
