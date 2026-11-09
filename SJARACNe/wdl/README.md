# SJARACNe2 WDL 1.3 / Sprocket workflow

The default local engine is **Sprocket 0.32.0** with WDL 1.3. It prepares inputs,
scatters seeded MI reconstructions, then runs native consensus/DPI and annotates
the final graph. The `sjaracne` command and `SJARACNe/` source paths are retained.

## Local use

Use Linux or WSL. Install this checkout in the Python environment used by its
tasks, and install Sprocket separately using the
[installation guide](https://sprocket.bio/getting-started/installation).
The CLI requires version 0.32.0. CI checks the Linux x86_64
[release archive](https://github.com/stjude-rust-labs/sprocket/releases/tag/v0.32.0)
against SHA256 `27dbf45e3ceaad9e13d416fdf44e8741af1a98c6a516723b27057e695cab0756`.

Choose the network type explicitly:

```bash
sjaracne local --network-type TF \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_TF.txt \
  -o brca100_tf

sjaracne local --network-type SIG \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_SIG.txt \
  -o brca100_sig
```

TF selects `pb=0.001`; SIG selects `pb=0.0005`. Override with `-pb`. The current
profile requires a network type even with an explicit probability. Supply inputs
using the [matrix and ID rules](../../README.md#inputs); the same full hub list
supplies native `-s` and `-l`.

## Default parameter values

The current profile uses 100 seeded 80% resamples, K20 candidate selection,
supporting-run mean MI and one final DPI pass. These are workflow defaults;
[bare native defaults differ](../../docs/BACKWARD_COMPATIBILITY.md#direct-native-executable).

| Setting | Current workflow default |
| --- | --- |
| Resampled networks | 100, seeds 1 through 100 |
| Sampling | 80% without replacement; `ceil(0.8 * N)` observations |
| AP-MI Npar limit | 40 |
| MI conversion | Experimental PCHIP |
| Per-resample probability | TF `0.001`; SIG `0.0005` |
| Order | Raw MI resamples, K-support consensus, supporting-run mean MI, then DPI |
| DPI | Dynamic epsilon, Href=479, slope=0.05, cap=0.30, additive offset=0.025 |
| Consensus | Edge present in at least 20 distinct resampled networks |
| Raw-graph hash re-verification | Off; operational option `--verify-raw-hashes` |

### Recurrence cutoff

K20 is an operating choice, not a biological FDR guarantee or a universal optimum.
Inspect target coverage: the bundled BRCA100 example becomes sparse at K20.
A smaller K may suit scRNA-seq, but no single-cell default or dynamic-K rule has
been validated. Override with `-k COUNT`; `-k 6` restores only the previous
recurrence cutoff. See the [benchmark summary](https://github.com/adamdingliang/sjaracne_benchmark/pull/4).

### Dynamic DPI tolerance

One scalar tolerance applies to the final consensus graph:

```text
epsilon(H) = min(dpi_max, dpi_offset + dpi_slope * max(0, log2(H / dpi_href)))
```

With the workflow defaults:

```text
epsilon(H) = min(0.30, 0.025 + 0.05 * max(0, log2(H / 479)))
```

`H` counts unique reconstructed source rows, including zero-target rows. Larger
epsilon relaxes pruning. At or below Href, epsilon equals the offset, subject to
the cap. The executable records the resolved value. Set `--dpi-offset 0` to test
the unshifted schedule; fixed DPI does not use the offset.

Retained MI is averaged over supporting runs only. Raw MI uses 17 significant
digits, and consensus-first final outputs preserve those tokens. DPI-first keeps
its historical four-decimal output format. See the
[native interface](../../docs/native_consensus_before_dpi.md).

### Interpretation and calibration limits

Edges are candidate associations, not verified direct or causal regulation.
TF/SIG probabilities are provisional operating choices, not biological FDR
controls or universal optima.

For a PCHIP cutoff (`pb != 1`), the empirical null surface requires Npar=40,
selected observation count `m=50..10000`, probabilities `1e-7..0.05` and sampling
without replacement. Unsupported inputs fail without fallback or extrapolation.
`pb=1` bypasses cutoff lookup and retains all MI candidates. Probabilities below
`0.0005` warn that they are outside the main assessment range; extreme-tail
assessment has limited Monte Carlo precision. Warnings remain in the logs.

Dynamic DPI adjusts pruning for reconstructed hub count; DREAM5 tests did not
establish a dynamic-specific accuracy advantage. See the
[DPI benchmark summaries](https://github.com/adamdingliang/sjaracne_benchmark/blob/main/docs/benchmark_summaries/README.md).

## Explicit alternatives

Use `--dpi-offset 0` to keep consensus first with the unshifted dynamic schedule,
or `--consensus-order dpi-first` to prune each resample first. DPI-first defaults
to zero offset. These choices do not select the full legacy profile:

```bash
sjaracne local --network-type TF --dpi-offset 0 \
  -e /absolute/path/input.exp -g /absolute/path/hubs.txt -o unshifted_tf

sjaracne local --network-type TF --consensus-order dpi-first \
  -e /absolute/path/input.exp -g /absolute/path/hubs.txt -o dpi_first_tf
```

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

`--tolerance 0` also selects strict fixed DPI. Fixed tolerance conflicts with
explicit dynamic-DPI tuning. `-M affine` changes MI conversion only. Probability
consensus (`-pc`) requires DPI-first; consensus-first requires recurrence K.
See the [compatibility guide](../../docs/BACKWARD_COMPATIBILITY.md) for legacy
settings and supported overrides.

## Direct WDL inputs

[`sjaracne.wdl`](sjaracne.wdl) supports both processing orders. Its workflow ID
remains `SJARACNe`, preserving `SJARACNe.*` input/output keys; its display name is
"Reconstruct a regulatory network".

| Step | Operation |
| --- | --- |
| `prepare_inputs` | Validate matrix/hub IDs and options, normalize line endings and record input/runtime hashes. |
| `resample_and_estimate_mi` | Resample by seed, estimate MI and apply the cutoff. Consensus-first emits raw MI; DPI-first prunes each graph here. |
| `build_consensus_network` | Check resample receipts, combine edges and annotate the final graph. Consensus-first performs K-support selection, supporting-run mean MI and DPI in C++; DPI-first combines pruned graphs in Python. |

Example inputs for a six-resample functional check:

```json
{
  "SJARACNe.exp_file": "/absolute/path/input.exp",
  "SJARACNe.probe_file": "/absolute/path/hubs.txt",
  "SJARACNe.network_type": "TF",
  "SJARACNe.bootstrap_num": 6,
  "SJARACNe.min_recurrence": 6
}
```

Save as `wdl-inputs.json`, then run from the repository root:

```bash
sprocket --skip-config-search --config SJARACNe/wdl/sprocket.local.toml \
  check SJARACNe/wdl/sjaracne.wdl
sprocket --skip-config-search --config SJARACNe/wdl/sprocket.local.toml \
  validate SJARACNe/wdl/sjaracne.wdl @wdl-inputs.json
sprocket --skip-config-search --config SJARACNe/wdl/sprocket.local.toml \
  run SJARACNe/wdl/sjaracne.wdl @wdl-inputs.json --output-dir /absolute/path/wdl-output
```

Omitted scientific inputs use the CLI's shared profile resolver. This example
sets K6 explicitly: fewer than 20 resamples require a smaller K, with no automatic
support-fraction rescaling. `profile: legacy` selects historical settings and
makes `network_type` optional.

Optional scientific inputs include `sampling_method`, `subsample_spec`, `depth`,
`p_value_bootstrap`, `apmi_null_model`, `min_recurrence`, `p_value_consensus`,
`tolerance`, `dynamic_dpi`, `dpi_href`, `dpi_slope`, `dpi_max`, `dpi_offset` and
`consensus_order`. Replacement and unique-sampling inputs cannot be combined.
`aracne_config_dir` selects a configuration directory.

The operational input `SJARACNe.verify_raw_hashes` defaults to `false`. To add
fresh graph-byte verification to either order, including legacy, add this entry
to the input JSON:

```json
{
  "SJARACNe.verify_raw_hashes": true
}
```

The CLI equivalent is `--verify-raw-hashes`. Direct CWL documents have different
input schemas and retained defaults; see the [CWL guide](../cwl/README.md).

## Resources, outputs and provenance

The local backend runs trusted tasks directly on the host, without container
isolation. CLI defaults are four task slots (`--jobs`) and a 16-GiB scheduling
budget (`--memory`); `--serial` uses one slot. Task requests are 2 GiB for
preparation, 4 GiB per reconstruction and 8 GiB for consensus, adjustable with
`--prepare-memory`, `--inference-memory` and `--consensus-memory`. Requests are
scheduling reservations, not measured peaks or hard limits.

The direct-WDL configuration uses two CPUs/two scatter slots and 16 GiB.
`--sprocket-config` selects a trusted custom backend configuration, which owns
concurrency and the overall budget; do not combine it with `--jobs` or `--serial`.
Per-task memory requests still enter WDL.

The CLI publishes four network/metadata files and `sjaracne_run.json` with resolved
options and engine paths. It refuses nonempty output directories. Logs, engine
work and task attempts are retained. `-tmp` chooses a scratch parent without
deleting existing contents.

Direct outputs include `out_dir` (the enhanced network **file**),
`consensus_network` (three columns), `bootstrap_info`, `parameter_info`,
`adjacencies`, `run_provenance`, `native_stdout`, `native_stderr`,
`validation_report`, `parameters`, `consensus_provenance`,
`consensus_native_stdout` and `consensus_native_stderr`. Use Sprocket's
`outputs.json` for their paths. See the [output formats](../../README.md#outputs).

Preparation records runtime, native, configuration and normalized input hashes.
Each seed records its settings and graph SHA-256. Consensus always checks ordered
successful seeds and matching fingerprints; by default it records producer
hashes without rereading graphs to recompute them. `--verify-raw-hashes` enables
that extra comparison. Input/settings/runtime and final-output checks remain in
place, as do native format checks. The receipt records the verification policy.
See [hash verification](../../docs/native_consensus_before_dpi.md#optional-raw-graph-hash-verification).

Consensus-first passes existing graph paths directly to C++, without aggregation
caches, adjacency staging copies or an intermediate candidate graph. Scatter
still writes/localizes raw graphs. Empty graphs count as resamples. Enhancement
uses the full original expression matrix; its regression p-value is not an MI or
consensus-edge probability. Publication checks final provenance and all four
output hashes. Input transformation, filtering and hub reordering are not added;
line endings are normalized on staged copies.

## Regression and validation limits

The small regression executes real default-CLI and paired CWL/WDL runs with
pinned Sprocket. A focused paired-engine run is also available:

```bash
python tests/wdl/run_parity.py --sprocket /absolute/path/sprocket \
  --output-dir /absolute/path/new-parity-evidence
```

The runner installs nothing and records versions, commands, hashes and logs.
Checks cover per-seed MI/edge bytes, consensus, enhancement, recurrence/mean-MI,
current TF/SIG settings, explicit DPI choices, legacy, CRLF and invalid inputs.
These are implementation checks, not biological-accuracy validation.

The [default-profile report](../../docs/default_profiles_validation.md) and
[initial parity report](../../docs/wdl13_parity.md) are historical snapshots.
The local configuration disables call caching and retries. Production scheduler,
cloud, container and recovery behavior require separate validation.
