# Acceptance fixtures

`cnn_5.txt` is the active end-to-end acceptance result for five fixed-size
80% subsamples without replacement. Its test explicitly selects DPI-first,
affine MI and probability consensus (`-pc 0.01`); it is a compatibility
fixture, not the current K20 consensus-first default.

[Acceptance tests](../../test_sjaracne_acceptance.py) cover these compatibility
cases. [Default-workflow tests](../../test_default_workflow.py) cover current
consensus-first behavior with generated inputs. Historical unused snapshots
remain in Git history; `cnn_5.txt` was not regenerated.
