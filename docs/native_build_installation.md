# Automatic native builds during installation

Source installation builds native programs automatically. For current setup,
see the [environment guide](uv.md).

The validation below is a **historical snapshot from October 29, 2026**,
following `6920ffe` on `feat/wdl13-scientific-parity`. Its counts and measurements
apply to that revision; [current CI](../.github/workflows/regression.yml) tests
the maintained checkout on Python 3.11 and 3.14.

## Source installation

On Ubuntu, including Ubuntu through WSL on Windows, install system prerequisites
once:

```bash
sudo apt update
sudo apt install --yes build-essential git curl
```

After installing uv, run the source installation below. Install Sprocket
separately before launching a WDL workflow:

```bash
uv sync --locked
source .venv/bin/activate
```

No separate `make -C SJARACNe` command is needed. The installer builds both
native executables through Make. Other Linux distributions need their
equivalent compiler/Make packages. Native Windows/MSYS2 source installation is
not supported; use an Ubuntu/WSL terminal rather than Git Bash/UCRT64.

## Build behavior and safety

Setuptools runs `build_native` for regular wheels and PEP 660 editable installs
before collecting package data. Metadata/source discovery does not compile.
Missing tools produce setup instructions; compilation failures fail installation.
The installer does not install system packages or run `sudo`.

Unchanged uv synchronization can reuse a cached installation. Native source,
header and Makefile changes participate in cache/build freshness. Compiler or
toolchain changes need an explicit rebuild; see the [environment guide](uv.md).

CI exercises the automatic build path through disposable uv editable/non-editable
and pip source/wheel/editable installations without prebuilt native programs.

## Validation

The October 29 checks used CPython 3.11.15 and 3.14.5 with the locked dependency
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

Each probe recorded 19 stages, including two expected missing-tool failures.
The Python 3.11 pip fallback also passed with the retained top-level requirement
snapshots, independently of a preinstalled regression environment.

An interrupted initial attempt remains marked incomplete. Final full-suite
results in `validation-resume1` use a fresh snapshot after correcting the pip
test helper's dependency handling.

Evidence is retained under
`D:\GitHub\local_runs\native_autobuild_20261029\`: `uv-install311` and
`uv-install314` contain the deployment probes; `pip-fallback311` records the
separate pip fallback; `validation-resume1` contains fresh pip installs and
both complete regressions. Source copies, environments and caches are not
committed. The preceding [uv migration report](uv_validation.md) remains
historical evidence for its tested revision.

The controller completed nine stages with unchanged frozen source hashes.
The paired six-gene/100-observation TF, SIG and tuned legacy cases produced
6/6/5 edges per environment, exact selected identities/MI, and zero observed
correlation/regression differences. These are small-fixture measurements.

## Conditional conclusion

The recorded Linux/WSL checks passed automatic source builds with Make and a
C++ compiler. They did not publish a compiler-free binary release or change
inference algorithms, coefficients, dependencies or scientific defaults.
Sprocket remains a separate installation.

This snapshot does not validate native Windows/MSYS2, arbitrary toolchains,
production HPC/cloud or biological accuracy. Hosted CI results and later
workflow changes are outside these dated measurements.
