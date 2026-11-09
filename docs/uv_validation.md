# Historical uv migration validation

This report records **October 29, 2026** validation on
`feat/wdl13-scientific-parity`, following `8a4f9f1`. Its six-resample/K6 results
are historical. Use the [environment guide](uv.md) and
[current workflow defaults](../SJARACNe/wdl/README.md#default-parameter-values)
for the maintained checkout.

The migration changed environment management and build freshness, preserving
sampling, MI estimation, PCHIP coefficients, DPI and consensus at that revision.

## Implementation

- Moved package metadata, dependencies, optional CWL support and the entry point
  to `pyproject.toml`, retaining setuptools and the Make build hook.
- Pinned uv 0.11.16 and transitive dependencies/artifact hashes. Numerical/CWL
  versions and build-tool pins matched the preceding Python 3.11/3.14 checks.
- Selected Python 3.14 by default and restricted uv resolution to 3.11-3.14;
  only 3.11 and 3.14 were tested. Historical package metadata remained `>=3.7.6`.
- Added native source/header/Makefile and runtime-asset cache keys, with Make
  dependencies to prevent stale objects.
- Prepared each CI lane once, then ran independent pip source-install,
  offline/cache and strict regression checks with separately pinned Sprocket.

Retained `requirements-*.txt` files are historical top-level snapshots, not the
active CI dependency source. Ordinary pip installation does not require uv.

## Validation evidence

The source-archive and deployment checks passed on CPython 3.11.15 and
3.14.5. Each clean source archive was compiled into a platform-specific wheel
and installed in a fresh environment through ordinary pip, using the full
locked export as constraints. Both fresh environments passed `pip check`,
installed-asset checks, native/CLI help and public PCHIP lookup checks.

The separate uv deployment probe passed all six checks in each Python lane:

1. A fresh non-editable core installation with no precompiled executables.
2. Unchanged source reused the installed/cached build without recompiling.
3. A second, fresh environment installed offline from the warmed cache.
4. A header-only edit to the disposable source invalidated the package cache
   and rebuilt all seven native objects.
5. A regulator-list edit refreshed the installed asset without unnecessary
   native compilation.
6. Imports, CLI/native help and PCHIP lookup worked with uv absent from `PATH`.

The restricted-PATH probe is not itself a complete network run: Sprocket is
deliberately absent from that probe. Full workflows are checked separately by
the regression and engine-parity suites.

An initial 3.14 probe changed the hash-protected PCHIP JSON and failed its
integrity check (`uv-install314`). The revised disposable regulator-list probe
passed in `uv-install314-resume1`; scientific assets remained unchanged.

The complete strict suite passed independently on each interpreter:

| Check | Python 3.11.15 | Python 3.14.5 |
| --- | --- | --- |
| Python tests | 314 passed; no failures, errors or skips | 314 passed; no failures, errors or skips |
| CWL document validation | 10 passed | 10 passed |
| Bats utility tests | 4 passed | 4 passed |
| Real WDL/CWL parity | 422 checks passed; zero failures | 422 checks passed; zero failures |
| Fresh source archive / native wheel / pip install | Passed | Passed |
| Locked install / offline reuse / native and asset cache probes | Passed | Passed |

The paired six-gene/100-observation scientific-output comparison also passed
all three configurations, using the same native executable, six seeds and
Sprocket version in both environments:

| Configuration | Consensus edges per environment | Native edges/MI and selected consensus | Maximum enhancement-statistic difference |
| --- | ---: | --- | ---: |
| Current TF defaults | 6 | Exact match | 0 |
| Current SIG defaults | 6 | Exact match | 0 |
| Controlled legacy probability case | 5 | Exact match | 0 |

The legacy fixture uses explicit `pb=1` and `pc=0.1` to produce informative
nonempty outputs, not the legacy default operating probabilities. Enhancement
statistics cover Pearson/Spearman correlations, slope and regression p-value.
The comparator permits `rtol=1e-12`, `atol=1e-14`, but observed no differences
here. Default dynamic epsilon is zero for this three-hub fixture; the full
regression separately covers positive dynamic tolerance.

Both completed strict runs preserved all 107 frozen source hashes. The
controller completed all seven stages successfully; separate uv-probe manifests
record the successful offline/cache checks. Evidence is retained under
`D:\GitHub\local_runs\uv_migration_20261029\`; environments/caches are not
committed or uploaded as source files.

## Conditional conclusion

At this snapshot, clean installations, warmed-cache offline reuse, native cache
refresh and complete small workflows passed on the recorded Python 3.11/3.14
stacks. Paired scientific outputs agreed with unchanged inference settings.

These measurements do not establish cold offline installation, universal
binary/numerical portability, biological accuracy or remote scheduler/cloud
support. They report local validation, not later hosted CI or full biological
benchmarks. See the [current regression guide](../tests/README.md) for maintained
checks and commands.
