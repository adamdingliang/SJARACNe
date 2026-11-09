# Environments and deployment with uv

Use uv **0.11.16** and the committed lockfile to install SJARACNe2. Source
installation builds the C++ programs automatically. Install system build tools
and **Sprocket 0.32.0** separately; CLI help does not check engine readiness.

For inference settings, see the [workflow guide](../SJARACNe/wdl/README.md).
The [installation validation](uv_validation.md) records the original uv migration.

## Contents

- [One-time setup on Ubuntu or Ubuntu in WSL](#one-time-setup-on-ubuntu-or-ubuntu-in-wsl)
- [Install SJARACNe2 with uv](#install-sjaracne2-with-uv)
- [Prepare Sprocket for workflows](#prepare-sprocket-for-workflows)
- [Choose a Python version](#choose-a-python-version)
- [Dependency sources and optional tools](#dependency-sources-and-optional-tools)
- [Native builds and cache invalidation](#native-builds-and-cache-invalidation)
- [HPC and cloud precautions](#hpc-and-cloud-precautions)
- [Ordinary pip alternative](#ordinary-pip-alternative)

## One-time setup on Ubuntu or Ubuntu in WSL

Use Linux or **Ubuntu in WSL** on Windows. Source installation requires Make
and a C++ compiler. Git Bash and MSYS2/UCRT64 are not supported substitutes.
On Ubuntu:

```bash
sudo apt update
sudo apt install --yes build-essential git curl
```

Use equivalent packages on other Linux distributions. The installer does not
install system packages or run `sudo`.

Install the pinned uv version using its
[official instructions](https://docs.astral.sh/uv/getting-started/installation/#installing-a-specific-version):

```bash
curl -LsSf https://astral.sh/uv/0.11.16/install.sh | sh
```

Open a new shell if needed, then check `uv --version`.

## Install SJARACNe2 with uv

From the repository root:

```bash
uv sync --locked
uv run --locked --no-sync sjaracne local --help
```

This creates `.venv`, installs Python dependencies and builds native programs.
The distribution name remains `SJARACNe`; no manual Make step is needed.
Missing build tools stop installation with setup instructions; see the
[native installation guide](native_build_installation.md).

Activate `.venv` with `source .venv/bin/activate` to use `sjaracne` directly.
Alternatively, use `uv run --locked --no-sync sjaracne ...`. `--no-sync` requires
an already-prepared environment; use `uv sync` for installation or updates.

Install Sprocket before launching the default WDL workflow.

## Prepare Sprocket for workflows

Install **Sprocket 0.32.0** using its
[official guide](https://sprocket.bio/getting-started/installation).
Choose the matching platform binary and verify its checksum; the
[WDL guide](../SJARACNe/wdl/README.md#local-use) records the Linux x86_64 checksum.

Put the executable on `PATH` and check:

```bash
sprocket --version
```

It must report `0.32.0`. For a default Cargo installation, add
`export PATH="$HOME/.cargo/bin:$PATH"` to `~/.bashrc` if needed. Alternatively,
pass `--sprocket /absolute/path/to/sprocket` to `sjaracne local`.
The CLI does not install the engine or silently fall back to CWL.

## Choose a Python version

CI tests **Python 3.11 and 3.14**. `.python-version` selects 3.14; uv can
download it when needed. To use an existing interpreter:

```bash
uv sync --locked --python /path/to/python3.14 --no-python-downloads
```

The package requires Python **3.11 or newer**. The lock resolves 3.11 through
3.14; intermediate versions are resolvable but not CI-tested. The locked
environment excludes Python versions outside that range.

For a separate Python 3.11 environment:

```bash
export UV_PROJECT_ENVIRONMENT=/path/to/sjaracne-py311
export UV_PYTHON=3.11
uv sync --locked --python 3.11
uv run --locked --no-sync sjaracne local --help
```

Keep both variables set for later `uv run` commands, or invoke the environment's
`bin/sjaracne` directly. `UV_PYTHON` avoids the checkout's 3.14 preference when
using `--no-sync`. Do not replace an environment used by running jobs.

## Dependency sources and optional tools

| File or setting | Purpose |
| --- | --- |
| [`pyproject.toml`](../pyproject.toml), `[project]` | Package metadata, runtime dependencies and optional `cwl` extra; also used by pip. |
| `[dependency-groups]` | Local build, regression and coverage tools. |
| `[tool.uv]` | Python-specific constraints, build constraints and native-source cache keys. |
| [`uv.lock`](../uv.lock) | Exact direct/transitive resolutions and artifact hashes. |
| [`.python-version`](../.python-version) | Preferred interpreter family. |

No dependency groups or CWL tools are installed by default. Add CWL support:

```bash
uv sync --locked --extra cwl
```

For the complete small regression suite:

```bash
uv sync --locked --group regression --extra cwl
uv run --locked --no-sync python tests/run_regression.py --jobs 4
```

The `regression` group includes pinned build tools. Install Make, a C++ compiler,
Node.js, Bats and Sprocket separately; see the [regression guide](../tests/README.md).
Keep the same extras/groups when resynchronizing: exact `uv sync` can remove
packages no longer selected.

`--locked` rejects a stale lockfile instead of updating it. Dependency upgrades
require an explicit lock update and validation in both CI lanes. A lockfile does
not guarantee identical numerical output across platforms or toolchains.

## Native builds and cache invalidation

Setuptools invokes Make during source installation. For development,
`make -C SJARACNe` remains available. Native sources, headers, the Makefile,
build metadata and runtime assets are cache keys; unchanged synchronization
can reuse a cached build.

After a compiler/toolchain change, rebuild in a dedicated environment with no
running jobs:

```bash
make -B -C SJARACNe
uv cache clean sjaracne
uv sync --locked --reinstall-package sjaracne
```

Keep your usual extras/groups in the last command. Clearing the package cache
prevents reuse of an older native wheel.

## HPC and cloud precautions

Prepare an environment or container **before submitting jobs** and leave it
unchanged during execution. Invoke `/path/to/prepared-env/bin/sjaracne` directly;
uv is not required to execute an installed package.

- **Isolation:** use a dedicated environment. Exact synchronization can remove
  unselected packages from a shared environment.
- **Storage:** follow site policy. `UV_PROJECT_ENVIRONMENT` and `UV_CACHE_DIR`
  relocate the environment and cache. Keep them on one filesystem for linking,
  or use `UV_LINK_MODE=copy`. Node-local storage must be prepared per node.
- **Offline nodes:** prepare the interpreter, dependencies and native programs.
  `uv sync --locked --offline` needs cached or local artifacts; a lockfile alone
  is not an offline bundle. A prepared environment avoids synchronization.
- **Portability:** build for the target OS, architecture and system libraries.
  Do not copy workstation virtual environments to unrelated cluster platforms.
- **Containers:** build native programs during image construction and use a
  non-editable installation. Pin and test the image alongside the lockfile.

uv does not schedule jobs or validate LSF, Slurm or cloud backends. The documented
workflow uses host execution without container isolation; container deployment
is guidance, not a validated backend.

Official references: [cache behavior](https://docs.astral.sh/uv/concepts/cache/),
[Python management](https://docs.astral.sh/uv/concepts/python-versions/) and
[containers](https://docs.astral.sh/uv/guides/integration/docker/).

## Ordinary pip alternative

uv is not a runtime dependency. For a basic pip source installation:

```bash
python3.14 -m venv .venv-pip
source .venv-pip/bin/activate
python -m pip install -r tests/requirements-build.txt
python -m pip install --no-build-isolation -c tests/requirements-ci-py314.txt -e .
sjaracne local --help
```

Pip also builds native programs; the same compiler/Make prerequisites apply.
Install Sprocket separately for WDL execution.

For Python 3.11, use that interpreter and `tests/requirements-ci.txt`. These
requirements files are historical top-level snapshots, not full transitive
locks. Pip does not read `uv.lock` and may resolve different transitive versions.

To supply the lock's complete runtime/build constraints, export them with uv:

```bash
uv export --locked --no-emit-project --no-hashes --group regression --extra cwl \
  --format requirements-txt --output-file locked-requirements.txt
```

Transfer the export with the package if needed, and use it as pip constraints
for build tools and installation. It includes conditional markers and excludes
the checkout. `--no-hashes` permits constraint-file use; it does not download
artifacts. Retain `uv.lock` for artifact hashes. The
[fresh-source check](../tests/README.md#ordinary-pip-and-fresh-source-installation)
verifies pip installation independently of uv.
