# Automatic native builds during installation

Scope: `feat/wdl13-scientific-parity`, following `6920ffe`, October 29, 2026
(America/New_York). This is an installation/build improvement, not a change to
network inference, dependency versions or scientific defaults.

## User-facing change

On Ubuntu, including Ubuntu through WSL on Windows, install system prerequisites
once:

```bash
sudo apt update
sudo apt install --yes build-essential git curl
```

After installing the pinned uv and separately supplied Sprocket versions, the
normal source installation is:

```bash
uv sync --locked
source .venv/bin/activate
```

No separate `make -C SJARACNe` command is needed. The installer builds both
native executables through Make. Other Linux distributions need their
equivalent compiler/Make packages. Native Windows/MSYS2 source installation is
not supported; use an Ubuntu/WSL terminal rather than Git Bash/UCRT64.

## Build behavior and safety

The setuptools build command now starts with a dedicated public `build_native`
subcommand. Both regular wheels and PEP 660 editable installs execute it before
Python package data is collected. Metadata/source discovery does not compile.
The subcommand exposes source/output mappings for editable installs; it does
not rely on an editable-specific private setuptools API.

Missing Make or the compiler produces an actionable message identifying the
missing tools and Ubuntu/WSL package commands. A failed native build fails the
installation instead of leaving a successful-looking Python-only package. The
installer never runs sudo or installs operating-system packages automatically.

Unchanged uv synchronization can reuse its existing cached installation; it
does not compile on every invocation. Native source/header cache keys and
Make dependencies remain in place. Advanced users may still call Make
directly, and compiler/toolchain changes need an explicit rebuild as described
in the [environment guide](uv.md).

CI no longer manually builds the checkout before `uv sync`, so it exercises
the same automatic build path as a new user. Independent clean-source probes
cover default editable and non-editable uv installs, plus pip editable and
source-archive/wheel installs. No prebuilt native programs are supplied to
those probes.

## Validation

Checks use CPython 3.11.15 and 3.14.5 with the unchanged locked dependency
versions. The clean-install tests deliberately exclude prebuilt executables
and object files; they never manually invoke Make before installation.

| Check | Python 3.11 | Python 3.14 |
| --- | --- | --- |
| Fresh default editable `uv sync --locked` builds both native programs | Passed | Passed |
| Fresh non-editable uv install and warmed-cache offline install | Passed | Passed |
| Missing Make/compiler errors identify the tools and setup commands | Passed: two expected installation failures | Passed: two expected installation failures |
| Editable header-only refresh recompiles all seven objects | Passed | Passed |
| Editable C++ source refresh recompiles only the changed object | Passed | Passed |
| Unchanged sync reuses its installation without recompilation | Passed | Passed |
| Configured assets refresh and runtime import/help/PCHIP probes | Passed | Passed |
| Source archive -> native wheel -> fresh ordinary pip install | Passed | Passed |
| Fresh ordinary pip editable install without a manual Make step | Passed | Passed |
| Complete Python regression | 322 passed; zero failures, errors or skips | 322 passed; zero failures, errors or skips |
| CWL document validation / Bats utilities | 10 / 4 passed | 10 / 4 passed |
| Real WDL/CWL parity | 422 checks passed; zero failures | 422 checks passed; zero failures |

Each uv install probe records 19 stages, including two intentionally failing
missing-tool installs; all stages satisfied their expected outcome. The
ordinary pip fallback was also independently checked on 3.11 using the retained
top-level requirement snapshots rather than the full uv export, so a missing
runtime dependency cannot be masked by preinstalled regression dependencies.

An initial full-suite attempt was interrupted after the pip test helper's
default dependency handling was corrected. Its incomplete regression manifest
is retained and is not counted as a completed run. The final passing full
suite evidence is in `validation-resume1`, against a fresh source snapshot.

Evidence is retained under
`D:\GitHub\local_runs\native_autobuild_20261029\`: `uv-install311` and
`uv-install314` contain the deployment probes; `pip-fallback311` records the
separate pip fallback; `validation-resume1` contains fresh pip installs and
both complete regressions. Source copies, environments and caches are not
committed. The preceding [uv migration report](uv_validation.md) remains
historical evidence for its tested revision.

The final controller completed all nine stages successfully. Both strict
runs retained identical frozen source hashes throughout. The additional paired
six-gene/100-observation comparison passed for current TF, current SIG and an
explicitly tuned legacy probability case: 6/6/5 consensus edges per
environment, exact native edges/MI and selected consensus identities, and zero
observed differences in emitted correlation/regression statistics. This is
small-fixture compatibility evidence, not universal numerical identity.

## Conditional conclusion

These checks support simplifying normal Linux/WSL source installation to
`uv sync --locked`, after installing system prerequisites once. A manual Make
command is no longer part of the normal install path. Source builds still
require Make and a C++ compiler; no compiler-free binary release is published
by this change. Sprocket remains separately installed.

No inference algorithms, model coefficients, dependency versions or scientific
defaults changed. The tests do not establish native Windows/MSYS2, production
HPC/cloud, arbitrary toolchain compatibility or biological accuracy. No full
dataset benchmark rerun or hosted GitHub Actions result is claimed.
