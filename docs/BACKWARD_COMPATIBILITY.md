# Backward-compatible workflows

Use `--profile legacy` for the scientific settings used before sampling without
replacement. It selects retained behavior in the current executable, not a
historical checkout.

## Current and historical profiles

| Setting | Current workflow profile | `--profile legacy` |
| --- | --- | --- |
| Local workflow/engine | WDL 1.3 / Sprocket; CWL selectable | Same |
| Network type | Required `--network-type TF` or `SIG`, even with explicit `pb` | Optional; does not change legacy probability |
| Per-resample sampling | `ceil(0.8*N)` observations without replacement | N draws with replacement, native `-r 1` |
| Native seeds | `-S 1` through `-S 100` by default | Same; `-r 1` is not a seed |
| MI conversion | Experimental PCHIP | Affine; omit native `-M` |
| Per-resample `pb` | TF `0.001`; SIG `0.0005` | `1e-7` |
| AP-MI Npar limit | 40 | 40 |
| Order | Raw MI -> K20/supporting-run mean MI -> DPI | Per-run DPI -> probability consensus |
| DPI | Dynamic; Href=479, slope=0.05, offset=0.025, cap=0.30 | Strict fixed tolerance, native `-e 0` |
| Consensus | Distinct-network recurrence K=20 | Normal-approximation `pc=1e-5` |
| Reconstruction/annotation | Same full hub list for native `-s` and `-l` | Same |

Legacy affine conversion uses the full observation count N; PCHIP calibration
is not used. Compatibility was checked against pre-sampling revision
`e59973a1891672a14f06d7847545c53368f3066c`; that revision is not downloaded or run
automatically.

### Run the historical profile

```bash
sjaracne local --profile legacy \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_TF.txt \
  -o legacy_brca100_tf
```

To use local CWL instead of Sprocket:

```bash
uv sync --locked --extra cwl
source .venv/bin/activate
sjaracne local --engine cwl --profile legacy \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_TF.txt \
  -o legacy_cwl_brca100_tf -tmp /absolute/path/cwl-scratch
```

CWL also requires Node.js. The retired CWLEXEC runner, `sjaracne lsf`, `-j` and
`--config-json` are unsupported. Old scheduler JSON is not converted to Sprocket
configuration. Configure other Sprocket backends explicitly with
`--sprocket-config`; cluster execution remains unvalidated here.

## Explicit overrides and saved inputs

Legacy permits historical `pb`, `pc`, depth, network-count and fixed-tolerance
overrides. It rejects PCHIP, dynamic DPI, unique subsampling and K consensus.
Use the current profile for mixed configurations.

Current-profile `-M affine`, `--fixed-dpi`, `--tolerance` or `-pb` changes only
that setting. `-pc` requires `--consensus-order dpi-first`. Likewise, `-k 6`
restores only the previous recurrence cutoff, not the legacy profile. Without
probability consensus, fewer than 20 networks require a smaller explicit K.

Direct saved `cwl/dpi_before_consensus.cwl` jobs retain affine `pb=1e-7`, strict
DPI and 80% sampling without replacement, but an omitted consensus cutoff selects
K20. Set `min_recurrence: 6` to preserve the prior cutoff. The public CLI instead
passes a complete resolved profile.

For a direct legacy CWL job, set these scientific inputs explicitly; replacement
and unique-sampling settings cannot be combined:

```yaml
sample_number: 1
subsample_spec: null
apmi_null_model: null
p_value_bootstrap: 1e-7
depth: 40
dynamic_dpi: false
tolerance: 0
min_recurrence: null
p_value_consensus: 1e-5
```

Here null `apmi_null_model` means affine for the retained CWL document. The CLI
accepts `-M affine` and translates it to that representation.

Final network/metadata files retain their names:

- `consensus_network_3col_.txt`
- `consensus_network_ncol_.txt`
- `bootstrap_info_.txt`
- `parameter_info_.txt`

The CLI also writes `sjaracne_run.json` with resolved settings and engine/output
paths. Sprocket retains per-seed graphs, provenance and logs. Choose an empty
output directory; existing results are not overwritten.

## Direct native executable

Bare native inference defaults differ from workflow defaults:

- `-N` defaults to 20; workflows pass 40.
- Omitted `-M` uses affine conversion when a probability cutoff is requested;
  workflow PCHIP is explicit `-M pchip`.
- `-r 1` requests full-size replacement sampling; unique subsampling uses `-u`.
  They cannot be combined.
- Fixed `-e` remains available; `-e 1` disables DPI.
- `-p 1` preserves all MI; explicit `-t`, including zero, overrides conversion.
- Single-hub, all-gene, conditional and adjacency-replay interfaces retain their
  constraints. Unsupported PCHIP/replacement or dynamic/replay combinations fail.

## What compatibility does not promise

`get_mi_cutoff(m, p)` still requires explicit numeric inputs; workflow defaults
do not change that Python API.

The legacy profile provides configuration/workflow compatibility, not universal
byte-identical historical reproduction. Input validation, duplicate handling,
provenance, compilers and numerical libraries can differ. Exact reproduction
requires the pinned historical runtime, inputs, configuration and environment.

A local six-gene/100-observation check compared the pre-sampling executable with
current legacy settings. Seeds 1 and 17 produced byte-identical adjacency
bodies, not necessarily metadata. This is narrow evidence, recorded in the
[historical validation report](default_profiles_validation.md).

Current calibration and dynamic-DPI choices remain experimental. See
[interpretation and calibration limits](../SJARACNe/wdl/README.md#interpretation-and-calibration-limits)
for their scope.

## Retaining the previous current workflow

Keep the preceding order and zero-offset schedule with:

```bash
sjaracne local --network-type TF --consensus-order dpi-first \
  -e input.exp -g hubs.txt -o previous_current
```

`--dpi-offset 0` instead keeps consensus first with the unshifted dynamic
schedule. Neither request selects legacy. Native `--raw-mi` and
`--consensus-inputs` are explicit modes; native offset defaults to zero and
requires `--dynamic-dpi`. The workflow supplies its resolved offset.

The CLI selects descriptive [CWL entry points](../SJARACNe/cwl/README.md#updating-older-cwl-paths)
by order; update saved import-only aliases accordingly. WDL retains
`wdl/sjaracne.wdl`, workflow ID `SJARACNe` and its public keys. Its task labels are
`prepare_inputs`, `resample_and_estimate_mi` and `build_consensus_network`.
Renamed tasks can invalidate old call-cache identities: resume already-started
work with its frozen workflow document. Adapter subcommands remain `prepare`,
`infer` and `consensus`.
