import sys
from setuptools import setup, find_packages
from distutils.command.install import install
from distutils.command.build import build
from subprocess import call

class Build(build):
    def run(self):
        cmd = ["make", "-C", "./SJARACNe"]
        if call(cmd) != 0:
            sys.exit(-1)
        build.run(self)

version = {}
with open("./version.py") as fp:
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
        "cwltool >= 3.0.20201117141248",
    ],
    python_requires=">=3.7.6",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    include_package_data=True,
    # test_suite="tests",
    entry_points={
        "console_scripts": [
            "sjaracne=SJARACNe.sjaracne:main",
            "sjaracne-calibrate-apmi-null=SJARACNe.calibrate_apmi_null:main",
        ]
    },
    cmdclass={
        'build': Build,
    }
)
