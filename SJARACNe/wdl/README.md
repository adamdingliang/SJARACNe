# WDL 1.3 / Sprocket workflow

This guide documents the retained SJARACNe2 `master` snapshot, where DPI
runs before recurrence consensus. For the current native consensus-before-DPI
workflow, see the [public `develop` guide](https://github.com/adamdingliang/SJARACNe2/blob/develop/SJARACNe/wdl/README.md).

The workflow prepares inputs, scatters native MI/DPI reconstructions, then
combines and enhances their networks in Python. The current profile uses
experimental PCHIP cutoffs and dynamic DPI. See the
[compatibility guide](../../docs/BACKWARD_COMPATIBILITY.md) for legacy settings and CWL.

## Local use

Use Linux or WSL. Build and install this source checkout in the Python environment
that will execute its tasks, and install **Sprocket 0.32.0** separately using the
[official installation guide](https://sprocket.bio/getting-started/installation).
For the Linux x86_64 archive from the
[0.32.0 release](https://github.com/stjude-rust-labs/sprocket/releases/tag/v0.32.0),
the SHA256 is
`27dbf45e3ceaad9e13d416fdf44e8741af1a98c6a516723b27057e695cab0756`.
Sprocket is installed separately. The CLI rejects a missing or incompatible engine.

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

## Default parameter values

| Setting | Current workflow default |
| --- | --- |
| Resampled networks | 100, seeds 1 through 100 |
| Sampling | 80% without replacement; native ceiling rule for the selected count |
| AP-MI Npar limit | 40 |
| MI conversion | Experimental PCHIP |
| Per-resample probability | TF `0.001`; SIG `0.0005` |
| DPI | Experimental dynamic epsilon, Href=479, slope=0.05, cap=0.30 |
| Consensus | Distinct-network recurrence K=20, inclusive; 20 of 100 by default |

### Recurrence cutoff

K20 is a practical operating point, not biological FDR control or a universal
optimum. Inspect target coverage: the bundled BRCA100 example becomes sparse.
Use `-k COUNT` to adjust it; `-k 6` restores only the prior cutoff. A smaller K
may suit scRNA-seq, but no single-cell default or dynamic K rule is validated. See the
[benchmark summary](https://github.com/adamdingliang/sjaracne_benchmark/pull/4).

### Dynamic DPI tolerance

Dynamic tolerance is one scalar for each reconstruction:

```text
epsilon(H) = min(dpi_max, dpi_slope * max(0, log2(H / dpi_href)))
```

`H` counts unique reconstructed source rows, including empty rows. Larger
epsilon relaxes pruning; `H <= Href` gives zero. The native executable records
the actual tolerance.

### Interpretation and calibration limits

Inferred edges are candidate associations, not verified direct or causal
regulation. The provisional TF/SIG probabilities are operating choices, not
universal optima or biological false-discovery-rate guarantees.

PCHIP lookup uses a bundled empirical null surface with Npar=40,
`m=50..10000` and probabilities `1e-7..0.05`. The special `pb=1` value bypasses
lookup to preserve all MI candidates. Unsupported lookup inputs fail rather
than falling back or extrapolating. Probabilities below `0.0005` warn about
limited assessment; extreme-tail Monte Carlo precision is limited. Warnings
remain in the retained logs.

Dynamic DPI relaxes pruning as the reconstructed hub count increases above its
reference. It addresses hub-context sensitivity as an engineering choice;
DREAM5 tests did not establish a dynamic-specific accuracy advantage. See the
[DPI benchmark summaries](https://github.com/adamdingliang/sjaracne_benchmark/blob/main/docs/benchmark_summaries/README.md)
for the evidence and its limitations.

## Explicit alternatives

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

The standalone document is [`sjaracne.wdl`](sjaracne.wdl), workflow `SJARACNe`.
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

Omitted options use the same profile resolver as the CLI. The example sets K6
explicitly; fewer than 20 resamples require a smaller K. For `profile: legacy`,
network type is optional. Other inputs include `aracne_config_dir`, `subsample_spec`, `depth`,
`p_value_bootstrap`, `apmi_null_model`, `min_recurrence`, `p_value_consensus`,
`tolerance`, `dynamic_dpi`, `dpi_href`, `dpi_slope` and `dpi_max`. Explicit sampling
mode uses `sampling_method`; conflicting replacement/unique inputs are rejected.

Direct CWL inputs retain affine/strict-DPI defaults but use K20 when selection
is omitted; set `min_recurrence: 6` for the prior cutoff. The public CLI always
forwards its fully resolved scientific profile to either engine.

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

The CLI publishes four final network/metadata files and `sjaracne_run.json`.
Use a new or empty output directory; engine work, logs and attempts are retained.
`-tmp` selects a scratch parent whose existing contents are preserved.

Direct WDL outputs include `out_dir` (the enhanced network **file**, a historical
name), `consensus_network` (three columns), `bootstrap_info`, `parameter_info`,
`adjacencies`, `run_provenance`, `native_stdout`, `native_stderr`,
`validation_report` and `parameters`. Sprocket's workflow `outputs.json` gives the
actual file paths; do not assume they are at the engine output root.

Every seed has a separate adjacency and provenance record. The preparation step
freezes runtime, native binary, configuration and normalized input hashes;
reconstructions also record the prepared-settings hash. Consensus checks the
expected seed set and hashes and stages **only** the adjacency files. Empty
networks still count as runs. Retained MI is averaged across runs containing an
edge; enhancement uses the **full original expression matrix**, not a sampled
subset. Its regression p-value is not an MI or consensus-edge probability.

No expression transformation, gene filtering or hub reordering is introduced.
Line endings are normalized on staged copies without changing original files.
WDL normalizes its own staged copies independently of the retained CWL helper.

## Regression and validation limits

The complete regression requires pinned Sprocket and runs real CLI and
CWL/WDL comparisons. CI verifies the engine archive checksum. For a focused run:

```bash
python tests/wdl/run_parity.py --sprocket /absolute/path/sprocket \
  --output-dir /absolute/path/new-parity-evidence
```

The runner installs nothing and preserves commands, versions, hashes and logs.
It compares paired per-seed edge/MI bytes, scientific metadata, consensus and
enhancement, with independent recurrence/mean-MI checks. Current TF/SIG defaults,
explicit fixed/dynamic settings and the legacy profile are covered. CRLF and
invalid-input checks are separate from biological accuracy.

The [profile](../../docs/default_profiles_validation.md) and
[initial parity](../../docs/wdl13_parity.md) reports preserve dated evidence.
They do not validate production-scale resources or remote/container backends.
The local configuration disables call caching and retries.
