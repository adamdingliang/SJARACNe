# Standalone WDL workflow (experimental)

This is Step 1 of the workflow migration: a **WDL 1.3** implementation of the
existing SJARACNe scientific workflow, intended for local execution with
[Sprocket](https://sprocket.bio/). The current Python CLI and CWL workflow remain
unchanged. Run this WDL directly; `sjaracne local` still invokes CWL.

The workflow is [`sjaracne.wdl`](sjaracne.wdl), named `SJARACNe`. It validates and
stages inputs, scatters independent seeded reconstructions, then runs the
existing consensus and network-enhancement code. It uses the installed
SJARACNe package and its native executable, not a second MI or DPI implementation.

## Inputs and scientific defaults

Expression and hub files have the same format and ID matching rules as the
[CWL inputs](../../README.md#inputs). The same hub list is passed as both native
`-s` (reconstruction sources) and `-l` (DPI annotations).

| Input | Default or requirement | Meaning |
| --- | --- | --- |
| `exp_file` | Required `File` | Tab-separated matrix beginning with `isoformId`, `geneSymbol`, then observations. |
| `probe_file` | Required `File` | Hub IDs matching the matrix's first column. |
| `aracne_config_dir` | Optional `Directory` | Omission uses the installed package's bundled configuration. |
| `bootstrap_num` | `100` | Number of resampled networks; legacy input name, not sampling with replacement. |
| `subsample_spec` | `"80%"` | Native `-u`: fixed-size sampling without replacement; a count such as `"80"` is also accepted. Percentages resolve using the native ceiling rule. |
| `depth` | `40` | Native AP-MI Npar limit, `-N`. |
| `apmi_null_model` | Omitted | Legacy affine conversion; set `"pchip"` to opt in to experimental PCHIP. |
| `p_value_bootstrap` | `1e-7` only in legacy mode | Must be supplied explicitly with `"pchip"`; no automatic TF/SIG probability choice. |
| `min_recurrence` | Effective default `6` | Keep an ordered edge appearing in at least K distinct input networks. K cannot exceed `bootstrap_num`. |
| `p_value_consensus` | Omitted | Explicit deprecated normal-approximation consensus mode, mutually exclusive with `min_recurrence`. |
| `tolerance` | Effective default `0` | Fixed DPI epsilon, native `-e`; `1` disables DPI. Conflicts with `dynamic_dpi`, including explicit zero. |
| `dynamic_dpi` | `false` | Opt in to experimental hub-count-dependent epsilon instead of fixed epsilon. |
| `dpi_href`, `dpi_slope`, `dpi_max` | Native defaults `479`, `0.05`, `0.30` when enabled | Optional dynamic tuning; requires `dynamic_dpi=true`. |

PCHIP and dynamic DPI retain the existing native validation, numerical behavior,
warnings and limitations. In particular, PCHIP is not a biological FDR guarantee.
The provisional TF/SIG probabilities in the main README are examples, not WDL
defaults.

Dynamic epsilon is one scalar per reconstruction:

```text
epsilon(H) = min(dpi_max, dpi_slope * max(0, log2(H / dpi_href)))
```

`H` is the number of unique reconstructed source rows, including empty rows;
it is not the annotation-list length or the number of emitted edges. Larger
epsilon relaxes pruning. The native executable computes and records this value.

## Local example and runner verification

Use Linux or WSL, with this checkout built and installed into the Python
environment used by the tasks. The native executable is resolved from the
installed package's `bin` directory. No separate CWL invocation is needed for
the WDL run; the retained CWL workflow is needed only for parity comparisons.

The local validation target is **Sprocket 0.32.0** with the bundled
[`sprocket.local.toml`](sprocket.local.toml). Its host backend requires no Docker,
but provides **no container isolation**: run only trusted workflows. The example
limits scheduling to two CPUs/two concurrent resamples and 2 GiB total task
memory. Its 1-GiB-per-task requirements are for small parity tests, not a
production bulk-network resource recommendation.

Install Sprocket separately from its official
[0.32.0 release](https://github.com/stjude-rust-labs/sprocket/releases/tag/v0.32.0).
The Linux x86_64 release archive SHA256 is
`27dbf45e3ceaad9e13d416fdf44e8741af1a98c6a516723b27057e695cab0756`.
This is optional tooling, not a new SJARACNe Python dependency. See the
[Sprocket configuration reference](https://sprocket.bio/reference/configuration)
and [WDL 1.3 specification](https://github.com/openwdl/wdl/blob/wdl-1.3/SPEC.md).

Example input JSON for a two-resample functional smoke test:

```json
{
  "SJARACNe.exp_file": "/absolute/path/BRCA100.exp",
  "SJARACNe.probe_file": "/absolute/path/BRCA100_TF.txt",
  "SJARACNe.bootstrap_num": 2,
  "SJARACNe.min_recurrence": 1,
  "SJARACNe.apmi_null_model": "pchip",
  "SJARACNe.p_value_bootstrap": 0.001
}
```

Save the edited JSON as `wdl-inputs.json`. From the repository root, with the
installed `python3` and `sprocket` on PATH:

```bash
sprocket --skip-config-search --config SJARACNe/wdl/sprocket.local.toml \
  check SJARACNe/wdl/sjaracne.wdl
sprocket --skip-config-search --config SJARACNe/wdl/sprocket.local.toml \
  validate SJARACNe/wdl/sjaracne.wdl @wdl-inputs.json
sprocket --skip-config-search --config SJARACNe/wdl/sprocket.local.toml \
  run SJARACNe/wdl/sjaracne.wdl @wdl-inputs.json --output-dir /absolute/path/wdl-output
```

The JSON keys above already include the workflow name. Sprocket retains outputs,
logs and task attempts under the selected output directory; the workflow's
`outputs.json` gives the actual file paths. Call caching and task retries are
disabled in this initial parity configuration.

This preserves strict DPI. To request dynamic DPI, add
`"SJARACNe.dynamic_dpi": true`; optional overrides use `SJARACNe.dpi_href`,
`SJARACNe.dpi_slope` and `SJARACNe.dpi_max`. Do not also supply `tolerance`.
For a full run, omit the smoke-test `bootstrap_num` and `min_recurrence`
overrides to restore 100 networks and K=6. To use the deprecated consensus
probability mode, supply `p_value_consensus` and omit `min_recurrence`.

## Outputs and provenance

| Workflow output | Content |
| --- | --- |
| `out_dir` | Enhanced network file, `consensus_network_ncol_.txt`; this legacy output name denotes a file. |
| `consensus_network` | Three-column `source`, `target`, `MI` consensus file. |
| `bootstrap_info` | Consensus selection and resample metadata. |
| `parameter_info` | Native header metadata from the first adjacency plus consensus settings. |
| `adjacencies` | Per-seed native networks in seed order. |
| `run_provenance` | Per-run command and execution provenance, including seeds. |
| `native_stdout`, `native_stderr` | Per-seed native logs, preserving cutoff warnings and DPI statistics. |
| `validation_report` | Hub-overlap validation result. |
| `parameters` | Resolved workflow parameters. |

The internal consensus output directory is fixed as `results`; engine output
placement is separate. The CWL input `final_out_dir_name` is not exposed here.
Native seeds are `1..bootstrap_num`, with names `TF_run_001.adj`, etc. The WDL
uses minimum-width three-digit formatting without repeating names above 999,
correcting the old CWL naming limitation. The `TF_run_` prefix is historical
and is also used for SIG networks.

Only adjacency files are staged into the consensus input directory. Empty
networks still count as resamples. Support counts distinct ordered pairs in
each run; retained MI is the mean across runs containing the edge, rounded by
the existing consensus code. Correlations, slopes and regression p-values in
the enhanced output use **all original observations**, not a sampled subset.
The regression p-value is not an MI or consensus-edge probability.

No expression transformation, gene filtering or hub reordering is introduced.
Input files are validated and copied with safe LF line-ending normalization;
the originals are not modified. The WDL does not copy the existing CWL helper's
failure to replace CRLF input when its optional output name is omitted. LF
inputs define the primary direct CWL/WDL parity comparison.

Scientific parity means the same per-seed edge identities and MI values,
consensus network, and relevant scientific metadata when using the **same
native executable, inputs and numerical environment**. Staged input paths and
output paths can differ in metadata. It is not independent validation of the
algorithm or a claim that experimental PCHIP/dynamic settings improve every
dataset. Native headers omit the random seed, so per-run provenance must be
retained rather than relying only on `parameter_info`.
The preparation step freezes input, native binary, configuration and Python
source hashes; each reconstruction records the prepared-settings hash too.
Consensus rejects mixed or altered evidence rather than accepting a directory
of unrelated networks.

## Small parity regression

Ordinary regression discovery includes the task-boundary unit tests without
requiring Sprocket. To additionally execute the actual WDL and CWL engines:

```bash
python tests/wdl/run_parity.py --sprocket /absolute/path/sprocket \
  --output-dir /absolute/path/new-wdl-parity-evidence
```

The runner installs nothing and refuses an existing evidence directory. It
checks PCHIP/strict DPI/default K=6, positive dynamic epsilon, omitted-probability
legacy affine mode and explicit legacy consensus filtering. It compares every
seed's edge/MI bytes and scientific headers, final consensus/enhancement and
metadata; only localized input/output path values may differ. Independent
recurrence/mean-MI checks prevent two matching engines from being the only
oracle. CRLF normalization and invalid-input failures are checked separately.
Commands, versions, hashes, logs and `summary.json` are retained.

See the [local validation report](../../docs/wdl13_parity.md) for measured results
and the conditional conclusion. Cluster execution and large-data workflow
validation are not established by these tiny local checks.

## Scope boundary

Step 1 does not replace the Python CLI, change CWL, add a scheduler backend, or
claim LSF/Slurm/cloud portability has been tested. The later CLI adapter, cache
and recovery behavior, resource policies and cluster tests remain separate
work. Local scientific parity must be established first; large benchmark
datasets are not required for that regression.
