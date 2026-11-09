# SJARACNe2 CWL workflows and steps

SJARACNe2 keeps the `sjaracne` command and the `SJARACNe/` source directory
for compatibility.

Use workflow names that state the processing order. The public CLI chooses the
entry point after resolving the scientific profile and `--consensus-order`.

| Entry point | Processing order | Intended use |
| --- | --- | --- |
| [`consensus_before_dpi.cwl`](consensus_before_dpi.cwl) | Resample and estimate significant MI edges -> consensus -> DPI | Current CLI default |
| [`dpi_before_consensus.cwl`](dpi_before_consensus.cwl) | Resample and estimate significant MI edges -> DPI in each resample -> consensus | Explicit `dpi-first` or the historical `legacy` profile |

`dpi-first` is a processing order, not a synonym for the `legacy` profile.
The current profile can explicitly select `dpi-first` while retaining current
sampling and MI settings. The `legacy` profile also restores historical sampling,
MI thresholds and probability-based consensus.

## Consensus before DPI

```text
prepare_inputs ---------+
                        +-> resample_and_estimate_mi (one job per seed)
create_resample_seeds --+      -> consensus_and_dpi
```

`prepare_inputs` validates expression and hub IDs, normalizes line endings in
staged copies and records resolved settings and input/runtime fingerprints.
The seed array is a scheduling input; a resample job also depends on preparation.

`resample_and_estimate_mi` selects observations using its seed, estimates MI and
retains edges above the significance cutoff. For consensus-first settings it
passes `--raw-mi`, so each output graph has not undergone DPI.

`consensus_and_dpi` waits for all successful resamples and verifies their
provenance. One C++ invocation retains edges supported by at least K runs,
averages each retained edge's MI over supporting runs only and applies DPI to
the consensus graph in memory. It then annotates only the final network using
the full expression matrix. Here, C++ describes the implementation; the step
name describes its scientific operations.

The CLI supplies `scientific_options` with `consensus_order: consensus-first`.
Direct callers must supply settings for that same order. For the current
profile, defaults are 100 seeded 80% resamples without replacement, K20 and
dynamic DPI with an additive offset of 0.025. The hub reference count remains
479, the slope 0.05 and the cap 0.30. Use the public CLI to resolve omitted
settings consistently with WDL; the shared adapter validates them again.

## Optional raw graph hash verification

For `consensus_before_dpi.cwl`, the operational input `verify_raw_hashes` is a
Boolean defaulting to `false`, separate from `scientific_options`. Enable it
with the public CLI's `--verify-raw-hashes`, or set `verify_raw_hashes: true` in a
direct CWL job, to reread each resample graph and compare its SHA-256 with the
producer receipt before consensus.

The producer hashes remain computed and recorded. With the default, the adapter
checks successful ordered seeds and matching runtime, input and settings
fingerprints, and records the producer graph hashes without freshly verifying
the graph bytes. Input/runtime/configuration, settings and final-output hash
checks remain in place, as do native row, ID, duplicate-edge and MI checks.
For unchanged inputs this option changes no scientific calculation or output.

The retained `dpi_before_consensus.cwl` compatibility workflow does not implement
this verification option. The public CLI rejects `--verify-raw-hashes` with CWL
DPI-first or CWL legacy runs rather than ignoring it. Use WDL/Sprocket when
fresh graph-hash verification is needed with DPI-first; WDL supports both orders.

## DPI before consensus

The compatibility workflow makes every operation visible in its step names:
`validate_inputs`, `create_resample_seeds`, `create_resample_graph_names`,
`check_expression_line_endings`, `check_hub_line_endings`,
`resample_estimate_mi_and_apply_dpi`, `collect_resample_networks` and
`build_consensus_network`.

The line-ending steps retain the historical utility behavior: with no explicit
output filename they expose the staged input and do not guarantee that the
returned file was converted to LF. Input preparation in the consensus-first
workflow uses its validated adapter instead.

The per-seed step applies DPI before writing each graph. The final Python step
merges those pruned graphs and annotates the result; it does not run DPI again.
Explicit current-profile `dpi-first` uses a zero additive dynamic-DPI offset by
default. The legacy profile instead uses strict fixed DPI and its historical
probability consensus.

## Updating older CWL paths

Use the descriptive filenames directly. The following import-only aliases have
been removed; update saved commands and custom CWL `run` references accordingly.

| Removed path | Use directly |
| --- | --- |
| `native_consensus_workflow.cwl` | `consensus_before_dpi.cwl` |
| `sjaracne_workflow.cwl` | `dpi_before_consensus.cwl` |
| `native_prepare.cwl` | `prepare_inputs.cwl` |
| `native_infer.cwl` | `resample_and_estimate_mi.cwl` |
| `native_consensus.cwl` | `consensus_and_dpi.cwl` |

The unused `get_consensus_network.cwl` and `get_enhanced_consensus_network.cwl`
shims were also removed: they called scripts no longer supplied by the package.
`create_consensus_network.cwl` performs consensus merging and annotation for
the DPI-first workflow through the supported Python module.

Public input/output keys and native filenames remain unchanged. In particular,
`out_dir` is the historical output key for an enhanced network **file**. Internal
Python adapter subcommands remain `prepare`, `infer` and `consensus`; changing
workflow labels does not change their interface or the native algorithms.
