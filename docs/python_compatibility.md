# Python 3.11 and 3.14 compatibility

This report records the **pre-uv** validation at `8a4f9f1`; its measurements and
requirements-file provenance are retained unchanged. The subsequent
[uv migration](uv.md) centralizes active dependencies in `pyproject.toml` and
locks transitive versions. See the [migration validation](uv_validation.md)
for checks of that later environment-management change.

Validation scope: Linux/WSL, on `feat/wdl13-scientific-parity`, following
`5188c52`. This is installation and small-workflow validation, not a new
biological benchmark. Full commands, source hashes and dependency inventories
are retained with the local evidence.

Checks completed October 28, 2026 (America/New_York); report finalized October
29. The README now recommends Python 3.14 with its tested constraints while
retaining Python 3.11 as a compatibility option.

## Tested environments

| Component | Compatibility lane | Newer-Python lane |
| --- | --- | --- |
| CPython | 3.11.15 | 3.14.5 |
| NumPy | 1.26.4 | 2.4.6 |
| pandas | 2.2.3 | 2.3.3 |
| SciPy | 1.13.1 | 1.17.1 |
| cwltool | 3.1.20240508115724 | 3.3.20260925135507 |
| Sprocket | 0.32.0 | 0.32.0 |

The 3.11 numerical/CWL pins are unchanged. The 3.14 lane deliberately uses
versions with compatible binary distributions, rather than trying to build the
old numerical stack on an unsupported interpreter. Both use pip 26.2.2,
setuptools 84.2.1, wheel 0.48.0 and build 1.6.1.

Exact top-level pins are in
[`requirements-ci.txt`](../tests/requirements-ci.txt),
[`requirements-ci-py314.txt`](../tests/requirements-ci-py314.txt) and
[`requirements-build.txt`](../tests/requirements-build.txt). These do not lock
every transitive dependency; the resolved inventories are retained in the
evidence. Python patch versions above describe the actual local checks, not a
requirement to use those patches indefinitely.

## Packaging and CI changes

- Replace legacy `distutils` imports with the public setuptools build command
  and declare the PEP 517 backend in `pyproject.toml`.
- Compile the standalone native programs before collecting package data, even
  when building from a clean source archive without cached executables.
- Include native programs, configuration, the frozen PCHIP model and WDL/CWL
  assets explicitly. Wheels have conservative interpreter/platform tags rather
  than incorrectly claiming platform independence.
- Run both Python versions in independent CI matrix lanes, with no allowed
  failures and version-specific evidence artifacts.
- Include packaging controls in the strict regression's source hashes.

The native estimator, cutoff coefficients, sampling, DPI and consensus rules
were not changed. The optional CWL dependency remains optional; the default
workflow still uses the separately installed Sprocket executable. Nothing is
published to PyPI by these checks.

The existing package metadata's historical `>=3.7.6` declaration is unchanged.
This work tests **3.11 and 3.14**, not every version admitted by that declaration;
it does not establish support for end-of-life interpreters. A minimum-version
policy change is separate from recommending a tested newer environment. The
`legacy` scientific profile does not require an obsolete Python interpreter.

## Validation

| Check | Python 3.11 | Python 3.14 |
| --- | --- | --- |
| Pinned environment and `pip check` | Passed | Passed |
| Clean source archive, compiled native wheel, fresh isolated install | Passed | Passed |
| Complete Python suite | 305 passed; no failures, errors or skips | 305 passed; no failures, errors or skips |
| Production CWL validation and Bats utilities | Ten documents and four utility tests passed | Ten documents and four utility tests passed |
| Real WDL/CWL parity matrix | 422 checks passed; zero failures | 422 checks passed; zero failures |

The source-install check builds in a temporary copy with prebuilt executables,
objects and build caches excluded. It inspects archive contents and wheel
metadata, installs the resulting wheel into a fresh environment outside the
checkout, and checks imports, executable permissions, native/CLI help, PCHIP
lookup and installed asset hashes. Its temporary environments are not included
in uploaded artifacts. This packaging check is not itself a full workflow run.

The strict regression separately exercises real current TF/SIG workflows,
positive dynamic tolerance, fixed tolerance, sampling, both consensus modes,
native correctness and paired WDL/CWL execution. It rejects failures, skipped
tests, empty discovery and source/fixture changes. No BRCA100, LUAD or DREAM5
resampling benchmark batches are launched.

### Cross-version scientific output comparison

The additional comparator holds the six-gene/100-observation fixture, three
hub sources, native executable, six seeds and Sprocket fixed. It runs three
paired configurations:

1. Current TF defaults, with six resamples and K=6.
2. Current SIG defaults, with six resamples and K=6.
3. The legacy scientific profile explicitly tuned to `pb=1`, `pc=0.1` so the
   small fixture has informative nonempty probability-selected outputs. These
   are not the legacy default probabilities.

Each run's consensus is checked independently from per-resample edges and
support counts. Empty output agreement is not sufficient. Native edge bodies
and scientific headers, selected edge identities and MI must match exactly;
only correlation/regression statistics and floating-point metadata permit
`rtol=1e-12`, `atol=1e-14`. Actual differences are recorded. With three hubs the
default dynamic tolerance resolves to zero; positive-tolerance validation comes
from the full regression, not these three paired cases.

All three cross-version comparisons passed:

| Configuration | Consensus edges in each environment | Native edges/MI and selected consensus identities/MI | Maximum absolute difference in emitted enhancement statistics |
| --- | ---: | --- | ---: |
| Current TF | 6 | Exact match | 0 |
| Current SIG | 6 | Exact match | 0 |
| Controlled legacy probability case | 5 | Exact match | 0 |

The enhancement comparison covers Pearson and Spearman correlations, slope and
regression p-value. The comparator allowed the declared roundoff tolerance,
but observed no numeric differences in these emitted statistics on this fixture.
This does not establish universal bitwise identity for arbitrary datasets.

## Reproduce

Use the matching manifest in an isolated environment, then install this checkout:

```bash
python -m pip install -r tests/requirements-build.txt
# Python 3.14; replace with tests/requirements-ci.txt for Python 3.11.
python -m pip install -r tests/requirements-ci-py314.txt
python -m pip install --no-deps --no-build-isolation -e .
python tests/check_source_install.py --output-dir /absolute/new/source-install \
  --requirements tests/requirements-ci-py314.txt
python tests/run_regression.py --jobs 4 --output-dir /absolute/new/regression
```

Install Sprocket 0.32.0, the compiler/Make, Node.js and Bats separately. See the
[regression guide](../tests/README.md) for the cross-version comparison command.
Do not run two full regressions concurrently against the same checkout: each
force-rebuilds its native executables.

Local evidence root: `D:\GitHub\local_runs\python_compatibility_20261028\`.
It contains per-version source-install manifests, archives, installed-asset
fingerprints, dependency inventories and regression results. Large execution
artifacts and temporary environments are not committed to the repository.

The original Python 3.14 controller stopped during its final parity stage,
after all 305 Python tests, ten CWL validations and four Bats tests had passed;
no failing test was recorded. Its incomplete manifest is preserved, not relabeled
as a completed controller run. The detached continuation verified unchanged
source/native hashes, reran the whole parity matrix in
`regression314/wdl_parity-resume1`, then completed the full Python 3.11 regression
and cross-version comparison. `controller_state.json` records the successful
aggregate completion; `cross-version/summary.json` records the output comparison.

## Conditional conclusion

Python 3.14 is a defensible recommended environment **with the tested dependency
constraints**: fresh installation and all small workflow/regression checks
passed, and the paired scientific outputs agreed. Python 3.11 remains a tested
compatibility option with its previous numerical/CWL pins.

These results do not guarantee arbitrary dependency upgrades, universal bitwise
reproducibility, biological accuracy, free-threaded Python, native Windows
execution or remote scheduler/cloud behavior. The native scientific algorithms
and default parameters are unchanged. No full dataset benchmark rerun or hosted
GitHub Actions result is claimed. A uv project/lockfile migration is not part of
this change.
