#!/usr/bin/env python3
"""Strict, small Linux/WSL regression; never run dataset resampling benchmarks.

Dependencies must already be installed. Build both native programs, run every
Python test, validate all CWL documents and execute the Bats utility tests.
Missing dependencies, empty discovery and skipped tests are failures.
"""

import argparse
import datetime
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
import unittest


ROOT = Path(__file__).resolve().parents[1]
PACKAGES = ("numpy", "pandas", "scipy", "cwltool")


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def collect_test_ids(suite):
    for item in suite:
        if isinstance(item, unittest.TestSuite):
            yield from collect_test_ids(item)
        else:
            yield item.id()


def strict_success(result, discovered):
    return bool(discovered and result.testsRun == discovered and result.wasSuccessful()
                and not result.skipped and not result.expectedFailures)


def build_command(jobs):
    # Force a clean compilation decision independently of Make's dependency
    # tracking and uv's cached package rebuild decisions.
    return ["make", "-B", "-C", ROOT / "SJARACNe", "-j", jobs, "all"]


def source_hashes():
    paths = [ROOT / name for name in ("setup.py", "pyproject.toml", "uv.lock",
                                     ".python-version", "version.py", "MANIFEST.in",
                                     "requirements.txt")]
    paths.append(ROOT / "SJARACNe/Makefile")
    for folder in ("SJARACNe/src", "SJARACNe/config", "SJARACNe/cwl", "SJARACNe/wdl", "tests", "tools", ".github/workflows"):
        paths.extend(path for path in (ROOT / folder).rglob("*")
                     if path.is_file() and "__pycache__" not in path.parts
                     and path.suffix not in (".pyc", ".pyo"))
    paths.extend((ROOT / "SJARACNe").glob("*.py"))
    paths.extend((ROOT / "SJARACNe/bin").glob("*.py"))
    return {str(path.relative_to(ROOT)): sha256(path) for path in sorted(set(paths))}


def python_suite(output):
    os.chdir(ROOT)
    sys.path.insert(0, str(ROOT))
    loader = unittest.TestLoader()
    suite = loader.discover(str(ROOT / "tests"), pattern="test_*.py")
    identities = list(collect_test_ids(suite))
    if loader.errors or not identities or len(identities) != len(set(identities)):
        raise RuntimeError("Empty, duplicate or failed test discovery: " + repr(loader.errors))
    result = unittest.TextTestRunner(stream=sys.stdout, verbosity=2).run(suite)
    record = {
        "discovered": len(identities), "tests_run": result.testsRun, "test_ids": identities,
        "failures": [{"test": test.id(), "detail": detail} for test, detail in result.failures],
        "errors": [{"test": test.id(), "detail": detail} for test, detail in result.errors],
        "skipped": [{"test": test.id(), "reason": reason} for test, reason in result.skipped],
        "expected_failures": [test.id() for test, _ in result.expectedFailures],
        "passed": strict_success(result, len(identities)),
    }
    (output / "python_tests.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return 0 if record["passed"] else 1


def verify_installed_runtime(installed, package_root=None):
    # CWL subprocesses normally scrub PYTHONPATH. Verify the package they will
    # import outside this checkout, not only the parent process's source import.
    package_root = ROOT / "SJARACNe" if package_root is None else Path(package_root)
    paths = list(package_root.glob("*.py")) + list((package_root / "bin").glob("*.py"))
    paths += list((package_root / "cwl").glob("*.cwl"))
    paths += [path for path in (package_root / "config").glob("*") if path.is_file()]
    hashes = {}
    for path in sorted(paths):
        relative = path.relative_to(package_root)
        imported = Path(installed) / relative
        if not imported.is_file() or sha256(imported) != sha256(path):
            raise RuntimeError("Installed SJARACNe differs at {}; run python -m pip install -e . "
                               "in this environment first".format(relative))
        hashes[str(relative)] = sha256(path)
    if not hashes:
        raise RuntimeError("No Python/CWL/configuration runtime files found")
    return {"path": str(installed), "sha256": hashes}


def preflight(env):
    if platform.system() != "Linux":
        raise RuntimeError("Use Linux or WSL for the complete native/CWL/Bats regression.")
    tools = {}
    for command in ("make", "g++", "cwltool", "bats", "sprocket"):
        path = shutil.which(command, path=env["PATH"])
        if not path:
            raise RuntimeError("Required command not found: " + command)
        tools[command] = path
    if not any(shutil.which(command, path=env["PATH"]) for command in ("node", "nodejs")):
        raise RuntimeError("Node.js is required by the CWL JavaScript expressions.")
    sprocket_version = subprocess.check_output([tools["sprocket"], "--version"], env=env, text=True).strip()
    if not re.search(r"^sprocket v?0\.32\.0(?:\s|$)", sprocket_version):
        raise RuntimeError("The WDL regression requires Sprocket 0.32.0, got " + sprocket_version)
    imported = subprocess.check_output([
        sys.executable, "-I", "-c",
        "import pathlib,SJARACNe; print(pathlib.Path(SJARACNe.__file__).resolve().parent)"
    ], cwd=tempfile.gettempdir(), env=env, text=True).strip()
    return {"tools": tools, "sprocket_version": sprocket_version,
            "packages": {name: importlib.metadata.version(name) for name in PACKAGES},
            "installed_python_runtime": verify_installed_runtime(Path(imported))}


def run(output, jobs):
    output.mkdir(parents=True, exist_ok=False)
    env = os.environ.copy()
    env.update({"PYTHONDONTWRITEBYTECODE": "1", "PYTHONNOUSERSITE": "1",
                "SJARACNE_REQUIRE_NATIVE_TESTS": "1", "OPENBLAS_NUM_THREADS": "1",
                "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
    native = {"SJARACNE_TEST_EXE": ROOT / "SJARACNe/bin/sjaracne.exe",
              "SJARACNE_APMI_NULL_EXE": ROOT / "SJARACNe/bin/apmi_null_generator.exe"}
    for variable, expected in native.items():
        if variable in env and Path(env[variable]).resolve() != expected:
            raise RuntimeError("Unset foreign native override before running this source tree: " + variable)
        env[variable] = str(expected)
    env["PATH"] = os.pathsep.join((str(ROOT / "SJARACNe/bin"),
                                  str(Path(sys.executable).absolute().parent), env.get("PATH", "")))
    env["PYTHONPATH"] = str(ROOT)
    manifest = {"schema": "sjaracne-small-regression-v1", "complete": False,
                "source_root": str(ROOT), "started_at_utc": now(), "python": sys.version,
                "controller_sha256": sha256(__file__), "steps": [],
                "scope": "small functional/native/CWL regression; no dataset benchmark resampling"}
    if os.environ.get("GITHUB_ACTIONS") == "true":
        manifest["github"] = {key: os.environ.get(key) for key in (
            "GITHUB_REPOSITORY", "GITHUB_SHA", "GITHUB_REF", "GITHUB_EVENT_NAME",
            "GITHUB_RUN_ID", "GITHUB_RUN_ATTEMPT", "GITHUB_WORKFLOW", "RUNNER_OS")}

    def persist():
        (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    def command(name, arguments):
        entry = {"name": name, "argv": [str(value) for value in arguments], "started_at_utc": now()}
        manifest["steps"].append(entry)
        persist()
        print("START " + name, flush=True)
        started = time.monotonic()
        log = output / (name + ".log")
        with log.open("w", encoding="utf-8") as handle:
            process = subprocess.run(entry["argv"], cwd=ROOT, env=env, stdout=handle, stderr=subprocess.STDOUT)
        entry.update(returncode=process.returncode, elapsed_seconds=time.monotonic() - started,
                     log=log.name, log_sha256=sha256(log))
        persist()
        if process.returncode:
            raise RuntimeError("{} failed; inspect {}".format(name, log))
        print("PASS " + name, flush=True)

    try:
        manifest["environment"] = preflight(env)
        env["SJARACNE_TEST_SPROCKET"] = manifest["environment"]["tools"]["sprocket"]
        manifest["sprocket_sha256"] = sha256(env["SJARACNE_TEST_SPROCKET"])
        manifest["source_sha256"] = source_hashes()
        command("build", build_command(jobs))
        manifest["native_sha256"] = {key: sha256(path) for key, path in native.items()}
        command("python_tests", [sys.executable, __file__, "--suite-only", "--output-dir", output])
        documents = sorted((ROOT / "SJARACNe/cwl").glob("*.cwl"))
        if not documents:
            raise RuntimeError("No CWL documents found")
        for document in documents:
            command("validate_" + document.stem, ["cwltool", "--validate", document])
        command("bats_tests", ["bats", "--formatter", "tap", ROOT / "tests"])
        tap = (output / "bats_tests.log").read_text(encoding="utf-8")
        plans = re.findall(r"^1\.\.(\d+)\s*$", tap, re.M)
        passed = len(re.findall(r"^ok \d+(?: |$)", tap, re.M))
        if len(plans) != 1 or not passed or passed != int(plans[0]) or re.search(r"#\s*skip\b", tap, re.I):
            raise RuntimeError("Bats must run a nonempty plan with every test passing and no skips")
        manifest["bats_passed"] = passed
        command("wdl_parity", [sys.executable, ROOT / "tests/wdl/run_parity.py",
                               "--sprocket", env["SJARACNE_TEST_SPROCKET"],
                               "--output-dir", output / "wdl_parity"])
        if source_hashes() != manifest["source_sha256"]:
            raise RuntimeError("Source or fixtures changed during regression")
        manifest["complete"] = True
    except Exception as error:
        manifest["error"] = str(error)
        print("FAIL: " + str(error), file=sys.stderr)
    finally:
        manifest["finished_at_utc"] = now()
        persist()
        print("Evidence: " + str(output), flush=True)
    return 0 if manifest["complete"] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=min(4, os.cpu_count() or 1), help="Native build workers")
    parser.add_argument("--output-dir", type=Path, help="New evidence directory; default is under the system temp directory")
    parser.add_argument("--suite-only", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.suite_only:
        if args.output_dir is None:
            parser.error("Internal suite mode requires --output-dir")
        return python_suite(args.output_dir)
    output = args.output_dir
    if output is None:
        output = Path(tempfile.mkdtemp(prefix="sjaracne-regression-")) / "evidence"
    try:
        return run(output.resolve(), args.jobs)
    except (OSError, RuntimeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    sys.exit(main())
