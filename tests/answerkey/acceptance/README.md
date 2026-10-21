# Acceptance fixtures

`cnn_5.txt` is the active end-to-end acceptance result. It is generated from
five fixed-size 80% subsamples without replacement.

The test explicitly uses legacy consensus `-pc 0.01` and affine MI cutoffs;
it is a supported-compatibility check, not the current K=6/PCHIP workflow.
Current-mode coverage uses tiny generated inputs in `test_sjaracne_acceptance.py`.

Unused 100-run historical snapshots and absolute-path cluster job files were
removed during test maintenance. They remain recoverable from Git history;
`cnn_5.txt` was not regenerated or changed.
