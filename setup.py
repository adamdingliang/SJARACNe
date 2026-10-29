from pathlib import Path
import subprocess

from setuptools import Distribution, find_packages, setup
from setuptools.command.build import build


ROOT = Path(__file__).resolve().parent


class Build(build):
    def run(self):
        # The native programs are standalone executables, not Python extensions.
        # Build them before package data is collected, including from a clean
        # source archive which has no precompiled binaries.
        subprocess.check_call(["make", "-C", str(ROOT / "SJARACNe")])
        super().run()


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
    }
)
