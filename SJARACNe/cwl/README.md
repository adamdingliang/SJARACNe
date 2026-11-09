# SJARACNe2 CWL workflows and steps

Select CWL with `--engine cwl`; the default CLI engine is WDL/Sprocket. The CLI
resolves the scientific profile and selects a document by `--consensus-order`.

| Entry point | Processing order | Intended use |
| --- | --- | --- |
| [`consensus_before_dpi.cwl`](consensus_before_dpi.cwl) | Resample and estimate MI -> consensus -> DPI | Default for `--engine cwl` |
| [`dpi_before_consensus.cwl`](dpi_before_consensus.cwl) | Resample and estimate MI -> per-resample DPI -> consensus | Explicit `dpi-first` or `legacy` |

`dpi-first` changes processing order only. The `legacy` profile also restores
historical sampling, MI thresholds and probability consensus; see the
[compatibility guide](../../docs/BACKWARD_COMPATIBILITY.md).

## Consensus before DPI

```text
prepare_inputs ---------+
                        +-> resample_and_estimate_mi (one job per seed)
create_resample_seeds --+      -> consensus_and_dpi
```

`prepare_inputs` validates IDs/options, normalizes staged line endings and records
input/runtime fingerprints. Each resample depends on both preparation and its
seed. `resample_and_estimate_mi` selects observations, estimates MI and applies
the cutoff; `--raw-mi` disables per-resample DPI.

`consensus_and_dpi` checks successful resample receipts, then uses one C++ call
to select edges with support >= K, average MI over supporting runs and prune the
in-memory graph. Python annotates the final network using the full matrix.

Direct callers supply resolved `scientific_options` with
`consensus_order: consensus-first` and a matching `bootstrap_num`. Current profile
defaults are 100 seeded 80% resamples without replacement, K20, and dynamic DPI
with Href=479, slope=0.05, cap=0.30 and additive offset=0.025. The public CLI
resolves omissions consistently with WDL; the shared adapter validates inputs.

## Optional raw graph hash verification

`consensus_before_dpi.cwl` accepts operational Boolean `verify_raw_hashes`,
default `false`, separately from scientific settings. Use CLI
`--verify-raw-hashes` or direct input `verify_raw_hashes: true` to compare current
graph bytes with each producer's SHA-256 before consensus.

Producers always compute hashes. With re-verification off, consensus checks
ordered successful seeds and matching runtime/input/settings fingerprints, and
records producer hashes without freshly checking graph bytes. Native format
checks and input/runtime/settings/final-output hash checks still apply. This
option changes no scientific calculation for unchanged inputs; see
[verification details](../../docs/native_consensus_before_dpi.md#optional-raw-graph-hash-verification).

The retained DPI-first CWL workflow does not support this option. The CLI rejects
it for CWL DPI-first and legacy runs. Use WDL/Sprocket for opt-in re-verification
with DPI-first.

## DPI before consensus

The retained workflow stages are `validate_inputs`, `create_resample_seeds`,
`create_resample_graph_names`, `check_expression_line_endings`,
`check_hub_line_endings`, `resample_estimate_mi_and_apply_dpi`,
`collect_resample_networks` and `build_consensus_network`.

The historical line-ending utilities expose their staged input when no output
filename is given; they do not guarantee LF conversion. Consensus-first uses
the validated preparation adapter instead.

Each resample applies DPI before writing its graph. Python merges and annotates
those pruned graphs without a second DPI pass. Current-profile DPI-first uses
zero dynamic offset by default; legacy uses strict fixed DPI and historical
probability consensus.

## Updating older CWL paths

Replace removed import-only aliases in saved commands and CWL `run` references:

| Removed path | Use directly |
| --- | --- |
| `native_consensus_workflow.cwl` | `consensus_before_dpi.cwl` |
| `sjaracne_workflow.cwl` | `dpi_before_consensus.cwl` |
| `native_prepare.cwl` | `prepare_inputs.cwl` |
| `native_infer.cwl` | `resample_and_estimate_mi.cwl` |
| `native_consensus.cwl` | `consensus_and_dpi.cwl` |

The removed `get_consensus_network.cwl` and `get_enhanced_consensus_network.cwl`
shims called scripts no longer supplied. DPI-first uses
`create_consensus_network.cwl` for supported merging and annotation.

Public input/output keys and native filenames are unchanged. `out_dir` is the
historical key for an enhanced network **file**. Internal Python adapter
subcommands remain `prepare`, `infer` and `consensus`.
