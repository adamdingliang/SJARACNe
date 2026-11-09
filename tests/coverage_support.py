"""Coverage instrumentation for the small regression, never package runtime.

Only an explicitly selected, disposable virtual environment is modified. A
temporary startup hook captures engine-launched Python even when CWL scrubs
environment variables. Ordinary regression and installed SJARACNe are unchanged.
"""

import importlib.metadata
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import sysconfig


NATIVE_CFLAGS = "-c -pipe -g -Wall -W -fPIC -O0 --coverage -fprofile-abs-path -fprofile-update=atomic"
NATIVE_OVERRIDES = ("CC=g++", "CFLAGS=" + NATIVE_CFLAGS, "LIB=-L /usr/local/lib --coverage")
HOOK_MODULE = "_sjaracne_regression_coverage"
HOOK_NAME = HOOK_MODULE + ".pth"


def remove_stale_counters(build):
    """Remove only generated counter files within the exact native build tree."""
    build = Path(build).absolute()
    if build.is_symlink():
        raise RuntimeError("Coverage build directory must not be a symlink")
    resolved = build.resolve()
    counters = list(build.rglob("*.gcda"))
    for path in counters:
        if path.is_symlink() or not path.resolve().is_relative_to(resolved):
            raise RuntimeError("Refusing a coverage counter outside the native build tree")
    for path in counters:
        path.unlink()
    return len(counters)


def write_hook(directory, config, audit):
    """Create two exclusively owned startup files; refuse preexisting hooks.

    Returns:
        dict: File paths and their original content, for ownership-checked cleanup.
    """
    directory, config, audit = map(Path, (directory, config, audit))
    module = directory / (HOOK_MODULE + ".py")
    hook = directory / HOOK_NAME
    if module.exists() or hook.exists() or module.is_symlink() or hook.is_symlink():
        raise RuntimeError("Coverage startup files already exist; use a fresh environment")
    code = '''"""Generated temporary regression hook, not part of installed SJARACNe."""
import json
import os
from pathlib import Path
import sys
import coverage
os.environ.pop("COVERAGE_PROCESS_CONFIG", None)
os.environ.pop("COVERAGE_FILE", None)
os.environ["COVERAGE_PROCESS_START"] = __SJARACNE_CONFIG_PATH__
coverage.process_startup()
(Path(__SJARACNE_AUDIT_PATH__) / (str(os.getpid()) + ".json")).write_text(
    json.dumps({"pid": os.getpid(), "argv": sys.orig_argv}) + "\\n", encoding="utf-8")
'''.replace("__SJARACNE_CONFIG_PATH__", repr(str(config))).replace("__SJARACNE_AUDIT_PATH__", repr(str(audit)))
    owned = {}
    try:
        for path, content in ((module, code), (hook, "import " + HOOK_MODULE + "\n")):
            with path.open("x", encoding="utf-8", newline="\n") as handle:
                handle.write(content)
            owned[path] = content
    except Exception:
        remove_hook(owned)
        raise
    return owned


def remove_hook(owned):
    """Remove only unchanged files created by this instrumentation session."""
    for path, content in owned.items():
        if path.is_symlink() or not path.is_file() or path.read_text(encoding="utf-8") != content:
            raise RuntimeError("Coverage startup file changed; refusing to remove " + str(path))
    for path in reversed(list(owned)):
        path.unlink()
    owned.clear()


def persisted_task_coverage(data_dir, audit_dir, package):
    """Verify saved data from actual CWL and WDL Python task processes.

    Startup markers alone do not prove the process saved execution data. Inspect
    per-process coverage databases before combining tools can consume them.
    """
    from coverage import CoverageData

    expected = {
        "wdl:prepare": "bin/wdl_workflow.py", "wdl:infer": "bin/wdl_workflow.py",
        "wdl:consensus": "bin/wdl_workflow.py",
        "SJARACNe.bin.QC_input": "bin/QC_input.py",
        "SJARACNe.bin.create_consensus_network": "bin/create_consensus_network.py",
    }
    saved = {key: 0 for key in expected}
    datasets = list(Path(data_dir).glob(".coverage.*"))
    for marker in Path(audit_dir).glob("*.json"):
        entry = json.loads(marker.read_text(encoding="utf-8"))
        argv = entry["argv"]
        labels = [key for key in expected if key in argv]
        if "SJARACNe.bin.wdl_workflow" in argv:
            index = argv.index("SJARACNe.bin.wdl_workflow")
            if index + 1 < len(argv):
                labels.append("wdl:" + argv[index + 1])
        for label in labels:
            if label not in expected:
                continue
            pattern = r"\.(?:pid)?" + str(entry["pid"]) + r"\."
            for path in datasets:
                if not re.search(pattern, path.name):
                    continue
                data = CoverageData(basename=str(path))
                data.read()
                if data.lines(str(Path(package) / expected[label])):
                    saved[label] += 1
                    break
    missing = [key for key, count in saved.items() if not count]
    if missing:
        raise RuntimeError("Missing persisted task coverage: " + ", ".join(missing))
    return saved


def summarize(python_report, cpp_report, root):
    """Return line-only and branch counts, rejecting empty or out-of-scope data."""
    root = Path(root).resolve()
    package = root / "SJARACNe"
    expected = set(package.glob("*.py")) | set((package / "bin").glob("*.py"))
    actual = set()
    for name, record in python_report["files"].items():
        path = (root / name).resolve()
        if path not in expected:
            raise RuntimeError("Unexpected Python coverage source: " + name)
        actual.add(path)
        if record["summary"]["num_statements"] and not record["summary"]["covered_lines"]:
            raise RuntimeError("Runtime Python module was never executed: " + name)
    if actual != expected:
        raise RuntimeError("Python coverage omitted runtime modules")
    totals = python_report["totals"]
    if not totals["num_statements"] or not totals["covered_lines"]:
        raise RuntimeError("Empty Python coverage denominator")
    if not cpp_report["line_total"] or not cpp_report["line_covered"] or not cpp_report["files"]:
        raise RuntimeError("Empty C++ coverage denominator")
    for item in cpp_report["files"]:
        path = (root / item["filename"]).resolve()
        if not path.is_relative_to(package / "src") or path.name == "apmi_pchip_data.h":
            raise RuntimeError("Unexpected C++ coverage source: " + item["filename"])
    return {
        "python": {"covered_lines": totals["covered_lines"], "lines": totals["num_statements"],
                   "line_percent": 100 * totals["covered_lines"] / totals["num_statements"],
                   "covered_branches": totals["covered_branches"], "branches": totals["num_branches"]},
        "cpp": {"covered_lines": cpp_report["line_covered"], "lines": cpp_report["line_total"],
                "line_percent": 100 * cpp_report["line_covered"] / cpp_report["line_total"],
                "covered_branches": cpp_report["branch_covered"], "branches": cpp_report["branch_total"]},
    }


class CoverageSession:
    """Own instrumentation files and report generation for one strict test run."""

    def __init__(self, root, output):
        self.root = Path(root).resolve()
        self.output = Path(output).resolve() / "coverage"
        self.owned = {}
        self.details = {}

    def prepare(self):
        """Require a disposable external venv, then install a temporary hook."""
        prefix = Path(sys.prefix).resolve()
        site_packages = Path(sysconfig.get_path("purelib")).resolve()
        if sys.prefix == sys.base_prefix or prefix.is_relative_to(self.root):
            raise RuntimeError("Coverage requires a disposable virtual environment outside the checkout")
        if not site_packages.is_relative_to(prefix):
            raise RuntimeError("Coverage startup location is outside the virtual environment")
        versions = {name: importlib.metadata.version(name) for name in ("coverage", "gcovr")}
        compiler = subprocess.check_output(["g++", "-dumpfullversion"], text=True).strip()
        gcov_banner = subprocess.check_output(["gcov", "--version"], text=True).splitlines()[0]
        version_match = re.search(r"\b(\d+\.\d+(?:\.\d+)?)\s*$", gcov_banner)
        gcov_version = version_match.group(1) if version_match else None
        if compiler != gcov_version:
            raise RuntimeError("GCC and gcov versions must match for coverage")
        self.output.mkdir(exist_ok=False)
        self.data = self.output / "data"
        self.audit = self.output / "process_startup"
        self.data.mkdir()
        self.audit.mkdir()
        self.config = self.output / "coverage.rc"
        self.config.write_text(
            "[run]\nbranch = True\nparallel = True\nsource = " + str(self.root / "SJARACNe")
            + "\ndata_file = " + str(self.data / ".coverage")
            + "\ndisable_warnings = no-data-collected\n\n[report]\nprecision = 2\n",
            encoding="utf-8",
        )
        build = self.root / "SJARACNe/build"
        if not build.resolve().is_relative_to(self.root):
            raise RuntimeError("Coverage build directory is outside the checkout")
        removed = remove_stale_counters(build)
        self.owned = write_hook(site_packages, self.config, self.audit)
        self.details = {"versions": versions, "gcc_gcov_version": compiler,
                        "native_build_overrides": list(NATIVE_OVERRIDES), "stale_counters_removed": removed,
                        "scope": "runtime Python and C++ core; shared AP-MI includes null-generator tests"}

    def close(self):
        """Remove the temporary startup hook, including after test failure."""
        remove_hook(self.owned)

    def report(self):
        """Create separate nonempty reports only after all strict tests finish."""
        if self.owned:
            raise RuntimeError("Remove coverage startup hooks before reporting")
        previous = Path.cwd()
        try:
            # coverage.py anchors relative filenames during initialization.
            os.chdir(self.root)
            return self._report_in_root()
        finally:
            os.chdir(previous)

    def _report_in_root(self):
        """Generate reports with canonical repository-relative filenames."""
        from coverage import Coverage

        self.details["persisted_task_coverage"] = persisted_task_coverage(
            self.data, self.audit, self.root / "SJARACNe")
        coverage = Coverage(config_file=str(self.config), data_file=str(self.data / ".coverage"))
        coverage.combine(data_paths=[str(self.data)], strict=True, keep=True)
        coverage.save()
        coverage.xml_report(outfile=str(self.output / "python.xml"))
        coverage.json_report(outfile=str(self.output / "python.json"), pretty_print=True)
        args = [sys.executable, "-m", "gcovr", "--root", str(self.root),
                "--gcov-executable", "gcov", "--filter", r"^SJARACNe/src/.*\.(cpp|h)$",
                "--exclude", r"^SJARACNe/src/apmi_pchip_data\.h$",
                "--cobertura", str(self.output / "cpp.xml"),
                "--json-summary", str(self.output / "cpp.json"), "--json-summary-pretty",
                "--txt", str(self.output / "cpp.txt"), str(self.root / "SJARACNe/build")]
        with (self.output / "gcovr.log").open("w", encoding="utf-8") as handle:
            subprocess.run(args, cwd=self.root, stdout=handle, stderr=subprocess.STDOUT, check=True)
        self.details.update(summarize(
            json.loads((self.output / "python.json").read_text(encoding="utf-8")),
            json.loads((self.output / "cpp.json").read_text(encoding="utf-8")), self.root))
        self.details["complete"] = True
        (self.output / "summary.json").write_text(json.dumps(self.details, indent=2) + "\n", encoding="utf-8")
        return self.details
