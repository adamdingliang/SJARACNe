# Native rank-cache correctness fixtures

These six small, synthetic files protect native MI-ranking and adaptive-partition
behavior. They are test inputs and expected outputs, not biological benchmark
datasets. The engine tests do not need `sjaracne_benchmark`.

The files were moved without byte changes from
`benchmarks/rank_cache/fixtures/` at SJARACNe commit
`edf7a2f` when the historical performance benchmarks were separated from the
engine repository. Historical benchmark tools, timing results and their
provenance are maintained in
[sjaracne_benchmark](https://github.com/adamdingliang/sjaracne_benchmark)
(currently private).

- `tied_counts.exp` and `tied_hubs.txt` exercise tied expression values and hub
  ordering.
- `tied_seed17_reference.tsv` and
  `tied_conditional_bootstrap_seed17_reference.tsv` contain expected MI values
  for seed 17, with and without conditional bootstrap selection. These fixed
  bootstrap references require a Linux glibc build.
- `mixed_quadrant_counts.exp` and `mixed_quadrant_hub.txt` exercise mixed leaf
  and recursive quadrants in adaptive partitioning.

The regression tests are in `tests/test_sjaracne_rank_cache.py` and
`tests/test_sjaracne_adaptive_partition_workspace.py`. Run them from the repository
root with:

```sh
python -m unittest discover -s tests -p 'test_sjaracne_rank_cache.py' -v
python -m unittest discover -s tests -p 'test_sjaracne_adaptive_partition_workspace.py' -v
```

Set `SJARACNE_TEST_EXE` to a built native executable if it is not available at
`SJARACNe/bin/sjaracne.exe` or on `PATH`. Tests skip if the required runtime is
unavailable.
