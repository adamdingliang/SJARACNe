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


version = {}
with (ROOT / "version.py").open(encoding="utf-8") as fp:
    exec(fp.read(), version)

setup(
    name="SJARACNe",
    version=version["__version__"],
    description="Gene network reverse engineering from big data",
    url="https://github.com/jyyulab/SJARACNe",
    author="Liang Ding, Jiyang Yu",
    author_email="liang.ding@stjude.org, jiyang.yu@stjude.org",
    license="See LICENSE.md",
    install_requires=[
        "pandas >= 1.2.3",
        "numpy >= 1.20.1",
        "scipy >= 1.7.3",
    ],
    extras_require={"cwl": ["cwltool >= 3.0.20201117141248"]},
    python_requires=">=3.7.6",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    include_package_data=True,
    package_data={
        "SJARACNe": [
            "bin/*", "config/*.txt", "config/*.json", "cwl/*.cwl",
            "wdl/*.wdl", "wdl/*.md", "wdl/*.toml",
        ],
    },
    distclass=NativeDistribution,
    # test_suite="tests",
    entry_points={
        "console_scripts": [
            "sjaracne=SJARACNe.sjaracne:main",
        ]
    },
    cmdclass={
        'build': Build,
    }
)
