# Default WDL 1.3 / Sprocket workflow

The public local workflow now uses **WDL 1.3 and Sprocket** by default. It
validates and stages inputs, scatters seeded raw-MI reconstructions, then runs
one native C++ consensus/DPI invocation and enhances the final graph. WDL
schedules these stages; it does not reimplement their scientific calculations.

This branch also makes experimental **PCHIP MI cutoffs and dynamic DPI** the
current-profile defaults. The retained CWL path and an explicit pre-sampling
legacy profile are documented in [BACKWARD_COMPATIBILITY.md](../../docs/BACKWARD_COMPATIBILITY.md).

## Local use

Use Linux or WSL. Build and install this source checkout in the Python environment
that will execute its tasks, and install **Sprocket 0.32.0** separately using the
[official installation guide](https://sprocket.bio/getting-started/installation).
For the Linux x86_64 archive from the
[0.32.0 release](https://github.com/stjude-rust-labs/sprocket/releases/tag/v0.32.0),
the SHA256 is
`27dbf45e3ceaad9e13d416fdf44e8741af1a98c6a516723b27057e695cab0756`.
The CLI checks the engine; it does not silently install Sprocket or fall back to
CWL when it is missing. Sprocket is not a Python package dependency.

Choose the network type explicitly:

```bash
sjaracne local --network-type TF \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_TF.txt \
  -o brca100_tf

sjaracne local --network-type SIG \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_SIG.txt \
  -o brca100_sig
```

TF selects `pb=0.001`; SIG selects `pb=0.0005`. These probabilities can be
overridden with `-pb`. The network type is required for the current profile even
with an explicit probability; no filename or hub-list biology is inferred.
Inputs retain the [main README format and ID-matching rules](../../README.md#inputs).
The same full hub list supplies native `-s` and `-l`.

| Setting | Current workflow default |
| --- | --- |
| Resampled networks | 100, seeds 1 through 100 |
| Sampling | 80% without replacement; native ceiling rule for the selected count |
| AP-MI Npar limit | 40 |
| MI conversion | Experimental PCHIP |
| Per-resample probability | TF `0.001`; SIG `0.0005` |
| Order | Raw MI resamples, K-support consensus, supporting-run mean MI, then DPI |
| DPI | Dynamic epsilon, Href=479, slope=0.05, cap=0.30, additive offset=0.025 |
| Consensus | Distinct-network recurrence K=20, inclusive; 20 of 100 by default |

Dynamic tolerance is one scalar for the final consensus graph:

```text
epsilon(H) = min(dpi_max, dpi_offset + dpi_slope * max(0, log2(H / dpi_href)))
```

H counts unique reconstructed source rows, including empty rows. It is not
annotation-list length or emitted-edge count. Larger epsilon relaxes pruning;
H at or below Href resolves to the offset, subject to the cap. The native
executable calculates and records the actual tolerance. Missing edges do not
enter the supporting-run MI mean denominator. Raw MI files retain 17 significant
digits, including every constructed source row even when it has no edges.
Native-order final three-column and enhanced outputs preserve the computed MI
tokens; the prior order keeps its historical four-decimal MI formatting.

Default status does not remove scientific limitations. The TF/SIG probabilities
are provisional operating points, not universal optima. PCHIP calibration does
not establish biological FDR, and dynamic tolerance does not guarantee more
accurate direct edges. Native principal-range and extreme-tail warnings remain
visible in the retained logs. Unsupported model inputs fail rather than silently
extrapolating or selecting another model.

## Explicit alternatives

`--dpi-offset 0` removes the added tolerance. `--consensus-order dpi-first
--dpi-offset 0` retains the previous current-profile order and tolerance schedule.
Consensus-first requires minimum recurrence; probability consensus requires
`--consensus-order dpi-first` or the legacy profile. Explicit fixed tolerance
bypasses the dynamic schedule and cannot be combined with an explicit offset.


```bash
# Current sampling and K=20, but strict fixed DPI:
sjaracne local --network-type TF --fixed-dpi \
  -e /absolute/path/input.exp -g /absolute/path/hubs.txt -o strict_tf

# Historical pre-sampling scientific profile, still executed with Sprocket:
sjaracne local --profile legacy \
  -e /absolute/path/input.exp -g /absolute/path/hubs.txt -o legacy_network

# Same historical profile through the retained CWL runner:
sjaracne local --engine cwl --profile legacy \
  -e /absolute/path/input.exp -g /absolute/path/hubs.txt -o legacy_cwl
```

`--tolerance 0` also selects strict fixed DPI unless dynamic DPI is explicitly
requested, in which case the conflict is rejected. `-M affine` chooses affine
conversion without implicitly changing other current-profile settings. Legacy
defaults, supported overrides and low-level native behavior are detailed in the
[compatibility guide](../../docs/BACKWARD_COMPATIBILITY.md).

## Direct WDL inputs

The standalone document is [`sjaracne.wdl`](sjaracne.wdl). Its display name is
"Reconstruct a regulatory network". The public workflow identifier remains
`SJARACNe`, preserving existing `SJARACNe.*` input and output keys.

| Step | Operation |
| --- | --- |
| `prepare_inputs` | Validate expression/hub IDs and options, normalize line endings and record input/runtime hashes. |
| `resample_and_estimate_mi` | For each seed, resample observations, estimate MI and apply its significance cutoff. The default consensus-first order emits raw MI without DPI; explicit DPI-first applies DPI here. |
| `build_consensus_network` | Verify all resample receipts, combine their edges and enhance the final graph. Consensus-first performs K-support selection, supporting-run mean MI and DPI in one C++ invocation; DPI-first combines already-pruned graphs. |

The WDL document supports both processing orders, so the final step has an
order-neutral name. The selected `consensus_order` determines where DPI runs.
Task names and descriptions change how the engine labels its logs; scientific
options, native commands and final output schemas are unchanged.

Example JSON for a six-resample functional check:

```json
{
  "SJARACNe.exp_file": "/absolute/path/input.exp",
  "SJARACNe.probe_file": "/absolute/path/hubs.txt",
  "SJARACNe.network_type": "TF",
  "SJARACNe.bootstrap_num": 6,
  "SJARACNe.min_recurrence": 6
}
```

Save it as `wdl-inputs.json`, then from the repository root:

```bash
sprocket --skip-config-search --config SJARACNe/wdl/sprocket.local.toml \
  check SJARACNe/wdl/sjaracne.wdl
sprocket --skip-config-search --config SJARACNe/wdl/sprocket.local.toml \
  validate SJARACNe/wdl/sjaracne.wdl @wdl-inputs.json
sprocket --skip-config-search --config SJARACNe/wdl/sprocket.local.toml \
  run SJARACNe/wdl/sjaracne.wdl @wdl-inputs.json --output-dir /absolute/path/wdl-output
```

Omitted scientific options resolve through the same shared profile code used by
the CLI. The reduced example explicitly uses K6; otherwise the current default
is K20 and fewer than 20 resamples require a smaller explicit K. This is not an
automatic support-fraction rule. Set `profile` to `legacy` to use that profile; its network type is
optional. Optional inputs include `aracne_config_dir`, `subsample_spec`, `depth`,
`p_value_bootstrap`, `apmi_null_model`, `min_recurrence`, `p_value_consensus`,
`tolerance`, `dynamic_dpi`, `dpi_href`, `dpi_slope`, `dpi_max`, `dpi_offset` and
`consensus_order`. Explicit sampling
mode uses `sampling_method`; conflicting replacement/unique inputs are rejected.

Do not copy the old standalone document's defaults from historical reports.
Omissions now select the current profile, not affine/strict settings. Raw CWL
keeps its older MI/sampling/DPI defaults, but an omitted consensus cutoff now
uses K20; specify `min_recurrence: 6` to retain its prior recurrence cutoff.
Engine selection through the public CLI always forwards a fully resolved
scientific profile. For consensus-first, the CLI uses
`cwl/consensus_before_dpi.cwl`, whose tasks invoke the same validated JSON
adapter as WDL. Its stages are `prepare_inputs`, `create_resample_seeds`,
`resample_and_estimate_mi` and `consensus_and_dpi`. Its scatter adjacencies, provenance and logs remain in the
recorded owned engine output directory. The retained low-level CWL workflow
is named `cwl/dpi_before_consensus.cwl`; `cwl/sjaracne_workflow.cwl` is a
forwarding compatibility alias. `cwl/native_consensus_workflow.cwl` forwards to
`cwl/consensus_before_dpi.cwl`. See the [CWL index](../cwl/README.md) for the
complete entry-point and step mapping. A smaller K may be appropriate for scRNA-seq, but neither
a scRNA-seq default nor dynamic K has been validated; inspect target coverage.

## Resources, outputs and provenance

The local host backend requires no Docker but provides **no container isolation**.
Run trusted workflows only. CLI local scheduling defaults to four jobs and a
16-GiB memory budget, configurable with `--jobs` and `--memory`; `--serial` uses
one job. Task memory defaults are 2 GiB for preparation, 4 GiB for each inference
and 8 GiB for consensus. Override them through `--prepare-memory`,
`--inference-memory` and `--consensus-memory`, or the corresponding WDL inputs.
These are scheduling requests, not a claim that every dataset fits those sizes.
The bundled direct-WDL configuration uses two CPUs/two scatter slots and a
16-GiB budget. An explicit `--sprocket-config` selects custom engine settings.

The CLI publishes the familiar four final network/metadata files in its output
directory and saves `sjaracne_run.json` with resolved options and engine output
locations. It refuses a nonempty output directory rather than overwriting a
previous run. Engine work, logs and task attempts are retained. `-tmp` remains
accepted for choosing a scratch parent; existing files there are not deleted.

Direct WDL outputs include `out_dir` (the enhanced network **file**, a historical
name), `consensus_network` (three columns), `bootstrap_info`, `parameter_info`,
`adjacencies`, `run_provenance`, `native_stdout`, `native_stderr`,
`validation_report`, `parameters`, `consensus_provenance`,
`consensus_native_stdout` and `consensus_native_stderr`. Sprocket's workflow
`outputs.json` gives the
actual file paths; do not assume they are at the engine output root.

Every seed has a separate adjacency and provenance record. The preparation step
freezes runtime, native binary, configuration and normalized input hashes;
reconstructions also record the prepared-settings hash. Consensus checks the
expected seed set and hashes. Consensus-first passes a manifest of the existing
localized adjacency paths to native C++ without copying raw graphs or writing
an intermediate mean-MI candidate graph. DPI-first retains Python consensus.
Empty networks still count as runs. Retained MI is averaged across runs
containing an edge; enhancement uses the **full original expression matrix**,
not a sampled subset. Its regression p-value is not an MI or consensus-edge
probability. Final publication checks successful consensus provenance and all
four final hashes before copying native-order outputs.

Raw scatter files still need to be written and transferred/localized across
nodes. Native integration removes Python aggregation caches, adjacency staging
copies and the candidate-graph write/read; it does not eliminate distributed
scatter I/O or provenance hashing. Raw files remain available according to the
engine's retention policy; no automatic deletion is introduced.

No expression transformation, gene filtering or hub reordering is introduced.
Line endings are normalized on staged copies without changing original files.
The old CWL line-ending helper's staging defect is not copied into WDL.

## Regression and validation limits

The complete small regression now requires pinned Sprocket and executes real
default CLI and CWL/WDL comparison tests. CI verifies the release archive SHA256
before installing the engine. A focused paired-engine run remains available:

```bash
python tests/wdl/run_parity.py --sprocket /absolute/path/sprocket \
  --output-dir /absolute/path/new-parity-evidence
```

The runner installs nothing and preserves commands, versions, hashes and logs.
It compares paired per-seed edge/MI bytes, scientific metadata, consensus and
enhancement, with independent recurrence/mean-MI checks. Current TF/SIG defaults,
explicit fixed/dynamic settings and the legacy profile are covered. CRLF and
invalid-input checks are separate from biological accuracy.

See the [historical default-profile validation report](../../docs/default_profiles_validation.md).
The [initial standalone parity report](../../docs/wdl13_parity.md) records the
earlier opt-in defaults and is historical evidence, not a description of the new
default profile. Cluster/cloud, container, production-scale memory, cache and
failure-recovery validation remain separate work; the local configuration keeps
call caching and retries disabled.
