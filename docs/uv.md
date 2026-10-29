# Environments and deployment with uv

Use uv **0.11.16** to prepare reusable SJARACNe environments from the committed
lockfile. uv installs the Python package and triggers its native C++ build;
it **does not install Sprocket** or system build tools. Package installation
and CLI help do not require Sprocket, but the default WDL workflow does.

This guide covers environment setup and deployment, not inference settings.
See the [validation report](uv_validation.md) for the tested installation,
cache and workflow checks.

## Contents

- [One-time setup on Ubuntu or Ubuntu in WSL](#one-time-setup-on-ubuntu-or-ubuntu-in-wsl)
- [Install SJARACNe with uv](#install-sjaracne-with-uv)
- [Prepare Sprocket for workflows](#prepare-sprocket-for-workflows)
- [Choose a Python version](#choose-a-python-version)
- [Dependency sources and optional tools](#dependency-sources-and-optional-tools)
- [Native builds and cache invalidation](#native-builds-and-cache-invalidation)
- [HPC and cloud precautions](#hpc-and-cloud-precautions)
- [Ordinary pip alternative](#ordinary-pip-alternative)

## One-time setup on Ubuntu or Ubuntu in WSL

Use a Linux terminal, or **Ubuntu in WSL** on Windows. Git Bash and
MSYS2/UCRT64 are not substitutes for WSL. A source installation needs Make and
a C++ compiler, even though the Python installer invokes them automatically.
For Ubuntu or Ubuntu in WSL, install those prerequisites once:

```bash
sudo apt update
sudo apt install --yes build-essential git curl
```

Other Linux distributions should use their package manager for the equivalent
compiler and Make packages. The installer does not run `sudo` or install
system packages automatically.

Install the pinned uv version using the
[official installation instructions](https://docs.astral.sh/uv/getting-started/installation/#installing-a-specific-version).
For the standalone installer on Linux/WSL:

```bash
curl -LsSf https://astral.sh/uv/0.11.16/install.sh | sh
```

Open a new Ubuntu shell if needed, then check `uv --version`.

## Install SJARACNe with uv

Run these commands from the repository root (the directory containing
`pyproject.toml` and `uv.lock`):

```bash
uv sync --locked
uv run --locked --no-sync sjaracne local --help
```

`uv sync --locked` installs SJARACNe and its Python dependencies into `.venv`,
and automatically builds the native programs. No separate Make step is needed.
Missing build tools stop installation with an actionable error; see the
[native installation validation](native_build_installation.md).

**This does not install Sprocket.** The `--help` command checks the CLI, not
whether the default workflow is ready to run. Complete the next section before
launching a WDL workflow.

Activate `.venv` with `source .venv/bin/activate` to use `sjaracne` directly.
Alternatively, `uv run --locked --no-sync sjaracne ...` uses the prepared
environment without activation or resynchronization. `--no-sync` is for an
already-prepared environment, not initial installation or dependency updates.

## Prepare Sprocket for workflows

Install **Sprocket 0.32.0 separately** using the
[official release](https://github.com/stjude-rust-labs/sprocket/releases/tag/v0.32.0).
A prebuilt binary avoids compiling its Rust dependencies. Choose the correct
platform and verify the checksum; the
[WDL guide](../SJARACNe/wdl/README.md#local-use) records the Linux x86_64 checksum.

Ensure its installation directory is on `PATH`, then check:

```bash
sprocket --version
```

The output must report version `0.32.0`. If a default Cargo installation is
not found, add `export PATH="$HOME/.cargo/bin:$PATH"` to `~/.bashrc` and open a
new shell. For a downloaded binary, use its actual installation directory.
Alternatively, pass `--sprocket /absolute/path/to/sprocket` to `sjaracne local`.

SJARACNe neither installs the engine automatically nor silently falls back to
CWL. Once both the package and engine are ready, follow the
[workflow examples](../SJARACNe/wdl/README.md#local-use).

## Choose a Python version

`.python-version` selects Python 3.14 by default. uv can download a compatible
interpreter if one is missing. To use an existing interpreter explicitly:

```bash
uv sync --locked --python /path/to/python3.14 --no-python-downloads
```

The lock resolves Python **3.11 through 3.14**. Only 3.11 and 3.14 have the
documented regression coverage; intermediate versions are not newly validated
by being resolvable. The package's historical `requires-python >=3.7.6`
metadata is retained for ordinary installers, not as a current support claim.
`uv sync` rejects interpreters outside its narrower locked range.

To maintain a separate Python 3.11 environment, choose a dedicated absolute
path rather than replacing an environment used by running jobs:

```bash
export UV_PROJECT_ENVIRONMENT=/path/to/sjaracne-py311
export UV_PYTHON=3.11
uv sync --locked --python 3.11
uv run --locked --no-sync sjaracne local --help
```

Keep both variables set for subsequent `uv run` commands, or call the prepared
environment's `bin/sjaracne` directly. Environment directories are local
artifacts; do not commit them. `UV_PYTHON` also avoids a version-mismatch
warning from the checkout's default 3.14 preference when using `--no-sync`.

## Dependency sources and optional tools

| File or setting | Purpose |
| --- | --- |
| [`pyproject.toml`](../pyproject.toml), `[project]` | Standard package metadata, runtime dependencies and optional `cwl` extra; usable by pip without uv. |
| `pyproject.toml`, `[dependency-groups]` | Local build and regression tools, not published runtime requirements. |
| `pyproject.toml`, `[tool.uv]` | Python-specific validated dependency constraints, build constraints and native-source cache keys. |
| [`uv.lock`](../uv.lock) | Resolved package versions, sources and artifact hashes, including transitive dependencies. |
| [`.python-version`](../.python-version) | Preferred interpreter family; not an installed Python binary. |

The lock fixes direct and transitive dependency versions. The
[Python compatibility report](python_compatibility.md) lists the validated
numerical stacks and optional cwltool versions; cwltool is not installed by
the default synchronization.

No dependency groups are selected by default. Add optional CWL support with:

```bash
uv sync --locked --extra cwl
```

For the complete small regression suite:

```bash
uv sync --locked --group regression --extra cwl
uv run --locked --no-sync python tests/run_regression.py --jobs 4
```

The `regression` group includes the pinned `build` group. Neither installs
Make, a C++ compiler, Node.js, Bats or Sprocket. Prepare these external tools
before running the complete suite; see the [regression guide](../tests/README.md).
Use the same extras/groups when resynchronizing: an ordinary exact `uv sync`
can remove packages that are no longer selected.

`--locked` rejects a stale lockfile rather than updating it. Do not regenerate
the lock as part of ordinary execution. Dependency upgrades should be explicit
reviewed changes, followed by both regression lanes. A lockfile does not prove
identical numerical output across operating systems, CPUs or compiler builds.

## Native builds and cache invalidation

Normal source installation builds C++ through setuptools and Make. For C++
development, `make -C SJARACNe` remains available as an explicit rebuild.
An unchanged cached `uv sync` does not necessarily recompile native programs.

Cache keys cover native sources, headers, the Makefile, build metadata and
runtime assets. Edits invalidate the cached package build; a changed compiler
or toolchain requires an explicit rebuild and reinstallation. In a dedicated
environment with no running jobs:

```bash
make -B -C SJARACNe
uv cache clean sjaracne
uv sync --locked --reinstall-package sjaracne
```

Include the same extras/groups you normally use in that last command. Cleaning
the package's build cache avoids reinstalling an older cached native wheel.

## HPC and cloud precautions

Prepare an environment or container **once before submitting jobs**, then
execute it without modification. Do not install or synchronize dependencies in
every resample task. For example, a job can invoke
`/path/to/prepared-env/bin/sjaracne` directly; uv need not be installed in the
runtime image solely to execute an already-installed package.

- **Isolation:** use a dedicated environment, not global Python, an unrelated
  Conda environment or an environment being modified by other analyses.
  Exact synchronization can remove unselected packages.
- **Storage:** on shared filesystems, follow the site's policy.
  `UV_PROJECT_ENVIRONMENT` and `UV_CACHE_DIR` relocate the environment and
  cache; placing them on the same filesystem permits efficient linking.
  `UV_LINK_MODE=copy` avoids linking across separate filesystems at the cost of
  extra copying. Node-local scratch must be prepared per node or supplied by a
  container; it is not automatically visible to other compute nodes.
- **Offline nodes:** prepare the interpreter, all runtime/build
  dependencies and native executables in advance. `uv sync --locked --offline`
  works only when all required artifacts are already cached or available
  locally. A lockfile alone is not an offline package bundle. Using the
  prepared environment directly avoids synchronization entirely.
- **Portability:** build for the target operating system, CPU architecture and
  compatible system libraries. Do not copy a workstation virtual environment to an
  unrelated cluster or an x86_64 environment to ARM cloud instances.
- **Containers:** install dependencies and compile the native
  programs during image construction; use a non-editable installation so the
  runtime does not depend on a separate source checkout. Pin and test the
  image as well as the lockfile.

uv is not a scheduler or workflow engine. It does not add or validate LSF,
Slurm, AWS or Azure execution; that remains Sprocket/backend/container work.
The current documented workflow uses host execution without container
isolation. Container deployment is a recommendation, not a newly implemented
or validated backend.

Official references: [cache behavior](https://docs.astral.sh/uv/concepts/cache/),
[Python management](https://docs.astral.sh/uv/concepts/python-versions/) and
[container integration](https://docs.astral.sh/uv/guides/integration/docker/).

## Ordinary pip alternative

uv is not a runtime dependency of SJARACNe. An ordinary source installation
still uses the standard setuptools backend:

```bash
python3.14 -m venv .venv-pip
source .venv-pip/bin/activate
python -m pip install -r tests/requirements-build.txt
python -m pip install --no-build-isolation -c tests/requirements-ci-py314.txt -e .
sjaracne local --help
```

Pip also builds the native executables automatically. The same compiler/Make
prerequisites apply; Sprocket still requires separate installation for WDL
execution, and `--help` alone does not verify it.

For Python 3.11, select that interpreter and use `tests/requirements-ci.txt`
instead. These retained requirements files are **historical top-level
validation snapshots**, not the authoritative dependency definitions or full
transitive locks. Pip does not read `uv.lock`; the simple fallback above may
resolve different transitive versions.

To supply pip with the lock's complete runtime/build constraints, export them
on a machine with uv, then transfer the export along with the package if needed:

```bash
uv export --locked --no-emit-project --no-hashes --group regression --extra cwl \
  --format requirements-txt --output-file locked-requirements.txt
```

Use that export as pip constraints for installing the build tools and the
package. It contains conditional Python/platform markers and excludes the
checkout itself. `--no-hashes` makes it usable as a pip constraint file; it does
not download package artifacts or create an offline bundle. Retain the
original lockfile for its recorded artifact hashes. The fresh source-install
test checks the pip path independently of an installed uv executable.
