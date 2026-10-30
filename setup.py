from pathlib import Path
import platform
import shutil
import subprocess

from setuptools import Command, Distribution, find_packages, setup
from setuptools.command.build import build
from setuptools.errors import SetupError


ROOT = Path(__file__).resolve().parent
NATIVE_PROGRAMS = (
    "SJARACNe/bin/sjaracne.exe", "SJARACNe/bin/apmi_null_generator.exe",
)


class Build(build):
    # Editable wheels execute build subcommands directly rather than Build.run.
    # A dedicated command is intentional: setuptools may silence failures from
    # customized build_py commands during editable installation.
    sub_commands = [("build_native", None)] + build.sub_commands


class BuildNative(Command):
    """Compile standalone native programs before package data is collected.

    Both regular and PEP 660 editable builds execute this public build
    subcommand. Native files are always generated in place; build_py copies
    them into regular wheels, and editable output mappings point to the source.
    Metadata and source-file discovery do not execute this command.
    """

    description = "build SJARACNe native executables using Make"
    user_options = []
    editable_mode = False

    def initialize_options(self):
        self.build_lib = None

    def finalize_options(self):
        self.set_undefined_options("build", ("build_lib", "build_lib"))

    def run(self):
        host = platform.system()
        if host not in ("Linux", "Darwin"):
            raise SetupError(
                "SJARACNe native source builds require Linux (including Ubuntu "
                "through WSL on Windows). Native Windows/MSYS2 source builds "
                "are not supported. Open an Ubuntu/WSL terminal and retry "
                "the installation there."
            )
        compiler = "clang++" if host == "Darwin" else "g++"
        missing = [name for name in ("make", compiler) if not shutil.which(name)]
        if missing:
            advice = (
                "On Ubuntu/Debian, including Ubuntu through WSL, run:\n"
                "  sudo apt update\n"
                "  sudo apt install build-essential\n"
                "Then retry uv sync --locked (or your pip install command)."
                if host == "Linux" else
                "Install Make and the Xcode command-line tools, then retry "
                "the installation."
            )
            raise SetupError(
                "SJARACNe cannot compile its native executables: missing "
                + ", ".join(missing) + " on PATH.\n" + advice
                + "\nThe installer does not install system packages or run sudo."
            )
        self.announce("Building SJARACNe native executables automatically", level=2)
        try:
            subprocess.check_call(["make", "-C", str(ROOT / "SJARACNe")])
        except (OSError, subprocess.CalledProcessError) as error:
            raise SetupError(
                "SJARACNe native compilation failed; the source installation "
                "cannot complete without its executables. See the compiler "
                "output above, fix the build error, and retry the installation."
            ) from error

    def get_outputs(self):
        self.ensure_finalized()
        return [str(Path(self.build_lib) / name) for name in NATIVE_PROGRAMS]

    def get_output_mapping(self):
        # Supported by setuptools' editable build protocol, including strict
        # editable mode. These methods work before run() creates the binaries.
        return dict(zip(self.get_outputs(), NATIVE_PROGRAMS))

    def get_source_files(self):
        sources = list((ROOT / "SJARACNe/src").glob("*.cpp"))
        sources.extend((ROOT / "SJARACNe/src").glob("*.h"))
        sources.extend((ROOT / name for name in (
            "SJARACNe/Makefile", "tools/apmi_null_generator.cpp",
        )))
        return [str(path.relative_to(ROOT)) for path in sorted(sources)]


class NativeDistribution(Distribution):
    def has_ext_modules(self):
        # Conservatively tag wheels for their build interpreter/platform rather
        # than advertising a platform-independent py3-none-any archive.
        return True


setup(
    # Metadata and dependencies live in pyproject.toml; keep only the custom
    # native build integration and package discovery/data configuration here.
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    include_package_data=True,
    package_data={
        "SJARACNe": [
            "bin/*", "config/*.txt", "config/*.json", "cwl/*.cwl",
            "wdl/*.wdl", "wdl/*.md", "wdl/*.toml",
        ],
    },
    distclass=NativeDistribution,
    cmdclass={
        'build': Build,
        'build_native': BuildNative,
    }
)
