# Environments and deployment with uv

Use uv **0.11.16** and the committed lockfile to install SJARACNe2 and build
its native programs. Install Make, a compiler and **Sprocket 0.32.0** separately.
CLI help does not require Sprocket; WDL execution does.

See the [workflow guide](../SJARACNe/wdl/README.md) for inference settings and
the [historical uv report](uv_validation.md) for migration evidence.

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

Use Linux or **Ubuntu in WSL** on Windows. Install the Ubuntu build tools once:

```bash
sudo apt update
sudo apt install --yes build-essential git curl
```

Other Linux distributions need equivalent compiler/Make packages. The installer
does not install system packages; native Windows/MSYS2 builds are unsupported.

Install the pinned uv version using the
[official installation instructions](https://docs.astral.sh/uv/getting-started/installation/#installing-a-specific-version).
For the standalone installer on Linux/WSL:

```bash
curl -LsSf https://astral.sh/uv/0.11.16/install.sh | sh
```

Open a new Ubuntu shell if needed, then check `uv --version`.

## Install SJARACNe2 with uv

Run these commands from the repository root (the directory containing
`pyproject.toml` and `uv.lock`):

```bash
uv sync --locked
uv run --locked --no-sync sjaracne local --help
```

`uv sync --locked` installs into `.venv` and builds the native programs. No
separate Make step is needed; see [native installation](native_build_installation.md).

The help command checks CLI loading. Complete the Sprocket setup below before
running a WDL workflow.

Activate `.venv` to use `sjaracne` directly, or use `uv run --locked --no-sync`
with the prepared environment. Synchronize explicitly after dependency changes.

## Prepare Sprocket for workflows

Install **Sprocket 0.32.0 separately** following the
[official installation guide](https://sprocket.bio/getting-started/installation).
Select version `0.32.0`, rather than an unpinned latest version.
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

A missing or incompatible engine fails explicitly. Then follow the
[workflow examples](../SJARACNe/wdl/README.md#local-use).

## Choose a Python version

`.python-version` selects Python 3.14 by default. uv can download a compatible
interpreter if one is missing. To use an existing interpreter explicitly:

```bash
uv sync --locked --python /path/to/python3.14 --no-python-downloads
```

The lock covers Python **3.11–3.14**; CI tests 3.11 and 3.14. Intermediate
versions are resolvable but unvalidated. Historical package metadata allows
`>=3.7.6`; that is not the maintained validation range.

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

The lock fixes direct and transitive versions. The
[Python compatibility report](python_compatibility.md) records the original
validated stacks. Optional cwltool is excluded from a default sync.

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

`--locked` rejects stale metadata instead of updating the lock. Review dependency
upgrades separately and run both regression lanes; a lock does not guarantee
identical numerics across platforms or compilers.

## Native builds and cache invalidation

Source installation builds C++ through setuptools/Make. An unchanged cached
sync can reuse its native build. For source development, `make -C SJARACNe`
remains available.

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

Prepare the environment or container **before submitting jobs** and keep it
unchanged. Tasks can call `/path/to/prepared-env/bin/sjaracne` directly without
uv or per-resample dependency synchronization.

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

The documented backend runs on the local host without container isolation.
These deployment recommendations do not validate a container, scheduler or
cloud backend.

Official references: [cache behavior](https://docs.astral.sh/uv/concepts/cache/),
[Python management](https://docs.astral.sh/uv/concepts/python-versions/) and
[container integration](https://docs.astral.sh/uv/guides/integration/docker/).

## Ordinary pip alternative

uv is not a runtime dependency of SJARACNe2. An ordinary source installation
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

For Python 3.11, use that interpreter and `tests/requirements-ci.txt`. These
requirements files are historical top-level pins, not full transitive locks.
Pip does not read `uv.lock`; export its constraints when exact pins are needed.

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
