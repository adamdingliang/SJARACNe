# uv migration validation

Scope: `feat/wdl13-scientific-parity`, following `8a4f9f1`, on October 29,
2026 (America/New_York). This changes package/environment management and native
build freshness, not sampling, MI estimation, PCHIP coefficients, DPI or
consensus settings. No full biological benchmark batches are rerun.

## Implementation

- Standard metadata, dependencies, optional CWL support and the entry point
  move to `pyproject.toml`. Setuptools and its existing Make build hook remain.
- uv 0.11.16 is pinned. The committed lock includes transitive dependencies
  and artifact hashes. Numerical/CWL versions match the previous Python
  3.11/3.14 validation; build tools retain their prior pins.
- `.python-version` selects 3.14. uv resolution is limited to 3.11 through
  3.14; only 3.11 and 3.14 are tested here. The historical published Python
  requirement remains `>=3.7.6`, not a claim of new support for old versions.
- Native sources, headers, Makefile and installed configuration/workflow
  assets participate in package cache invalidation. Make additionally tracks
  native headers and recipe changes so a rebuild cannot reuse stale objects.
- CI prepares a dedicated locked environment once per Python lane. Independent
  pip source-archive installation, offline/cache probes and the strict small
  regression remain required. Sprocket is installed separately with its
  existing version/checksum pin.

The retained `tests/requirements-*.txt` files describe the preceding top-level
validation/pip alternative; they are not the active CI dependency source.
Ordinary pip installation does not require uv. See the
[environment and deployment guide](uv.md) for commands and HPC/cloud limits.

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

An initial 3.14 probe edited whitespace in the hash-protected frozen PCHIP
JSON and correctly failed its integrity check. That failed fixture is retained
in `uv-install314`; the corrected probe changes only the disposable regulator
list, and passed in `uv-install314-resume1`. No scientific assets in the
checkout were changed to accommodate the probe.

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

The uv migration is ready for review and routine local/CI use with these tested
Python 3.11/3.14 environments: clean installations, offline reuse of a prepared
cache, native cache invalidation and complete small workflows passed. The
tested paired scientific outputs agreed without changing inference settings.

This does not validate every interpreter/platform in package metadata, a cold
offline installation without prepared artifacts, universal binary/numerical
portability, biological accuracy, or LSF/Slurm/AWS/Azure execution. Ordinary pip
installation remains available, and uv is not required solely to execute an
already-installed runtime. No full biological benchmark rerun or hosted GitHub
Actions result is claimed; the updated CI workflow has been checked locally.
