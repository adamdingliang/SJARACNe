# Current defaults and pre-sampling compatibility: validation

Validation date: 2026-10-28 (America/New_York). Work is on the existing
`feat/wdl13-scientific-parity` branch, following Step 1 commit `d824d7a`.

The complete small regression passed: **284 Python tests**, **ten CWL
validations**, **four Bats tests**, and **422 assertions** in the expanded
engine matrix, with no failed checks or skips.

## Implemented defaults

| Setting | Current profile | Explicit legacy profile |
| --- | --- | --- |
| Default workflow/engine | WDL 1.3 / Sprocket 0.32.0 | Same; CWL also remains selectable |
| Role | Explicit TF or SIG required, including when pb is supplied | Optional; does not change legacy pb |
| MI conversion and pb | PCHIP; TF 0.001, SIG 0.0005, explicit pb overrides | Affine; pb=1e-7 unless overridden |
| Sampling | 80% without replacement | N draws with replacement, native `-r 1` |
| DPI | Dynamic Href=479, slope=0.05, cap=0.30 | Fixed epsilon=0 unless explicitly tuned |
| Consensus | Inclusive K=6 | Normal-approximation pc=1e-5 unless overridden |
| Networks/Npar | 100 / 40 | 100 / 40 |

The defaults apply to public workflows/CLI, not every bare native invocation.
The C++ estimator, cutoff implementation and DPI algorithm were not modified.
PCHIP and dynamic DPI remain experimental scientific operating choices.

## Checks

| Validation | Result |
| --- | --- |
| Full ordinary Python/native/CLI suite | 284 tests passed, zero skips/errors/failures |
| Production CWL documents | Ten validated |
| CWL utility tests | Four passed |
| Expanded real CWL/WDL matrix | Seven paired configurations plus CRLF/fail-fast checks; 422 assertions passed |
| Historical/current native legacy check | Seeds 1 and 17 had identical adjacency bodies on the tiny fixture |
| Source package | Built successfully; required WDL/profile/configuration/documentation assets and `version.py` verified; CWL dependency remains optional |
| Hosted GitHub Actions | Workflow updated; no hosted result is claimed by these local checks |

The ordinary regression includes real default CLI runs for TF and SIG, legacy
WDL/CWL engine equivalence, missing-role rejection and incompatible-option
failures. It also retains the existing native correctness, input, probability,
consensus and explicit older-CWL acceptance tests rather than regenerating
their answerkeys to fit the new defaults.

The expanded direct-engine matrix uses the existing six-gene, 100-observation,
three-source fixture and paired seeds:

| Configuration | Resamples per engine | Selection |
| --- | ---: | --- |
| PCHIP TF pb=0.001, explicitly fixed strict DPI | 6 | K=6 |
| Same, positive dynamic epsilon with Href=2, slope=0.1, cap=0.3 | 6 | K=6 |
| Explicit older-develop affine pb=1e-7, strict DPI, 80% unique sampling | 6 | K=6 |
| PCHIP with explicit legacy consensus probability | 3 | pc=0.01 |
| Current TF defaults, omitted MI/DPI options | 6 | K=6 |
| Current SIG defaults, omitted MI/DPI options | 6 | K=6 |
| Pre-sampling legacy profile, full-size replacement sampling | 6 | pc=1e-5 |

Each paired run compares adjacency bodies and scientific metadata exactly.
Only localized input/output path values are normalized. Consensus/enhanced
networks and bootstrap metadata must agree; independent recurrence/mean-MI
and legacy-tail calculations check selection rather than treating engine
agreement alone as an oracle. Current defaults are inspected in native headers,
commands and per-seed provenance. CRLF/input preservation and invalid-input
fail-fast checks are included separately.

The fixture's H=3 makes the default dynamic epsilon zero. The separate Href=2
configuration exercises a positive epsilon, approximately 0.05849625. Neither
is a biological benchmark or a test of large-panel memory requirements.

## Direct historical-binary check

For seeds 1 and 17 on the same tiny fixture, the verified pre-sampling
`e59973a1891672a14f06d7847545c53368f3066c` binary and the current binary produced
**identical adjacency bodies**, including printed MI values, with `-r 1`,
affine pb=1e-7, Npar=40 and epsilon=0. This is a narrow valid-input check, not
universal historical bitwise reproducibility or header identity.

The historical binary SHA256 was
`eda4f07465160f9d2d536715023d439b36628b1e094dab8b7afed511efa24e11`.
The current binary SHA256 was
`dbdbe074f1e6263cddf296dd0f9ebbacc437edbb47e0904ffc79f34a119a4fde`.
No obsolete executable is bundled or required by ordinary regression/CI.

## Evidence and reproduction

Local artifacts, not committed to the repository:

- Full regression: `D:\GitHub\local_runs\wdl13_defaults_regression_20261028_r3`.
  Its `manifest.json`, `python_tests.json`, logs and nested `wdl_parity/summary.json`
  record commands, versions, results and tested source/binary hashes.
- Historical binary comparison:
  `D:\GitHub\local_runs\wdl13_defaults_legacy_native_20261028\summary.json`.
- Source-package check:
  `D:\GitHub\local_runs\wdl13_defaults_package_20261028_r2\sjaracne-0.2.1.tar.gz`.

```bash
python tests/run_regression.py --jobs 2 --output-dir /absolute/path/new-evidence
```

CI installs checksum-pinned Sprocket 0.32.0, optional CWL dependencies and runs
the full small suite. It does not launch BRCA100/LUAD/DREAM5 benchmark batches.
Local success does not by itself claim a completed hosted GitHub Actions run.

## Conditional conclusion

These checks support changing the public workflow defaults while retaining
explicit configuration compatibility with the pre-sampling workflow on the
tested local cases. They do not show improved biological edge accuracy, an optimal probability
for every dataset, or universal PCHIP calibration. Direct-native interfaces
remain compatible; current safety checks and metadata are not rolled back.
Production-scale resources, remote LSF/Slurm/cloud execution, containers,
call caching and failure recovery remain outside this local validation.

See the [compatibility guide](../BACKWARD_COMPATIBILITY.md) and
[WDL guide](../SJARACNe/wdl/README.md). The earlier
[Step 1 report](wdl13_parity.md) is preserved as historical evidence.
