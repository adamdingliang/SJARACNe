# Backward-compatible workflows

This branch changes the public workflow defaults. An unchanged command is not a
promise of an unchanged inferred network. Use the explicit `legacy` profile when
you want the scientific settings used before the sampling-without-replacement
change.

## Current and historical profiles

| Setting | Current profile, now the workflow default | `--profile legacy` |
| --- | --- | --- |
| Local workflow/engine | WDL 1.3 / Sprocket | WDL 1.3 / Sprocket; CWL remains selectable |
| Network type | Explicit `--network-type TF` or `SIG`, even when supplying `pb` | Optional; does not select a different legacy probability |
| Per-resample sampling | 80% without replacement, using `ceil(0.8*N)` observations | N draws with replacement, native `-r 1` |
| Native random seeds | `-S 1` through `-S 100` by default | Same; `-r 1` is not the random seed |
| MI conversion | Experimental PCHIP | Historical affine conversion, no native `-M` |
| Per-resample `pb` | TF `0.001`; SIG `0.0005` | `1e-7` |
| AP-MI Npar limit | 40 | 40 |
| DPI | Experimental dynamic epsilon; Href=479, slope=0.05, cap=0.30 | Fixed strict tolerance, native `-e 0` |
| Consensus | Distinct-network recurrence K=6 | Explicit legacy normal-approximation `pc=1e-5` |
| Reconstruction/annotation | Same full hub list supplied to native `-s` and `-l` | Same |

The legacy settings were checked against pre-sampling revision
`e59973a1891672a14f06d7847545c53368f3066c`. The profile uses the retained
replacement-sampling implementation in the **current** native executable. It
does not download, check out or run that historical revision automatically.
For replacement sampling, the affine conversion uses the full observation count
N; no PCHIP calibration is used.

### Run the historical profile

```bash
sjaracne local --profile legacy \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_TF.txt \
  -o legacy_brca100_tf
```

This restores the complete legacy scientific profile while using Sprocket to
orchestrate it. To use the retained CWL path instead:

```bash
python -m pip install -e '.[cwl]'
sjaracne local --engine cwl --profile legacy \
  -e tests/inputs/BRCA100.exp -g tests/inputs/BRCA100_TF.txt \
  -o legacy_cwl_brca100_tf -tmp /absolute/path/cwl-scratch
```

CWL additionally needs Node.js. CWLEXEC and the old `sjaracne lsf` command have
been removed, including their configuration templates and launchers. They are
not installed, invoked or retained as an optional execution dependency.
Old `lsf`, `-j` and `--config-json` invocations now fail instead of submitting
jobs or silently switching engines.

Use WDL/Sprocket with an explicit `--sprocket-config` for separately configured
scheduler backends. Cluster execution remains unvalidated here; old CWLEXEC
JSON files are not Sprocket configurations and are not converted automatically.
The `legacy` scientific profile remains available with either Sprocket or local
`cwltool`; it does not depend on the retired scheduler runner.

## Explicit overrides and saved inputs

Historical `pb`, `pc`, depth, resample-count and fixed-tolerance tuning remain
possible in the legacy profile. PCHIP, dynamic DPI, unique subsampling and K-based
consensus conflict with the definition of that profile and are rejected rather
than silently ignored. Use the current profile for mixed modern configurations.

The current profile also supports explicit `-M affine`, `--fixed-dpi`,
`--tolerance`, `-pb` and `-pc` settings. Selecting **one** old option does not
restore every historical setting. For example, `-M affine` still leaves the
current TF/SIG probability, without-replacement sampling, dynamic DPI and K=6 in
place unless you explicitly change them.

Raw saved CWL jobs retain the previous `develop` defaults: affine `pb=1e-7`,
strict DPI, 80% **without** replacement and K=6. Those defaults are not the
pre-sampling legacy profile. The public CLI injects its complete resolved profile
into CWL. A raw legacy CWL job must explicitly set `sample_number: 1`,
`subsample_spec: null`, affine mode, fixed tolerance and the legacy consensus
probability. Do not set native replacement and unique-sampling inputs together.

The scientific portion of such a raw CWL input is:

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

Here null `apmi_null_model` means affine for the retained CWL document; the
public CLI uses the explicit spelling `-M affine` and performs that translation.

The familiar final files remain:

- `consensus_network_3col_.txt`
- `consensus_network_ncol_.txt`
- `bootstrap_info_.txt`
- `parameter_info_.txt`

The new CLI additionally saves its resolved scientific options and engine/output
locations in `sjaracne_run.json`. Sprocket retains per-seed adjacencies, command
provenance and native logs. A nonempty output directory is not silently reused or
overwritten.

## Direct native executable

Low-level `sjaracne.exe` defaults remain unchanged. This is deliberate: blindly
applying workflow defaults would break legitimate all-gene, single-hub,
conditional and adjacency-replay uses. In particular:

- Omitting native `-M` still selects affine conversion when a probability cutoff
  is requested; workflow PCHIP is passed explicitly as `-M pchip`.
- `-r 1` retains full-size replacement sampling; `-u` is the separate unique
  subsampling option. They cannot be combined.
- Fixed `-e` remains available; `-e 1` disables DPI.
- Native `-p 1` preserves all MI. An explicit `-t`, including zero, overrides
  probability conversion.
- Existing single-hub, all-gene, conditional and replay interfaces retain their
  constraints. Unsupported PCHIP/replacement or dynamic/replay combinations fail
  explicitly; there is no silent model fallback.

## What compatibility does not promise

The existing Python `get_mi_cutoff(m, p)` callable is also retained. It requires
explicit numeric inputs; workflow role defaults do not change that API.

This is **configuration and scientific-workflow compatibility**, not universal
byte-for-byte historical reproduction. The current executable includes newer
input validation, safety checks, duplicate handling and provenance. Numerical
libraries, compiler versions and platforms may also differ. Exact historical
reproduction requires the pinned historical source/executable, inputs,
configuration and numerical environment.

A narrow local compatibility check compared the verified pre-sampling native
executable at `e59973a` with the current executable using the legacy settings.
On the six-gene, 100-observation synthetic fixture, paired seeds **1 and 17**
produced byte-identical adjacency **bodies**. This is useful evidence that the
retained replacement/affine/strict-DPI path reproduces that small example; it
does not establish identical metadata or universal fidelity across datasets,
compilers and platforms. Artifact provenance is recorded in the
[current validation report](default_profiles_validation.md).

The new PCHIP and dynamic-DPI defaults remain experimental scientific operating
choices. Default status does not establish biological FDR control, direct-edge
accuracy, universal calibration or an optimal TF/SIG probability for every
dataset. See the [workflow guide](../SJARACNe/wdl/README.md) and
[current validation report](default_profiles_validation.md).
