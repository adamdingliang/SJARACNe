"""Stdlib-only contracts for opt-in regression coverage instrumentation."""

from contextlib import ExitStack
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest import mock

import coverage_support as support


class TestCounterCleanup(unittest.TestCase):
    def test_removes_only_counters_in_the_selected_build_tree(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            build = root / "SJARACNe" / "build"
            nested = build / "tools"
            nested.mkdir(parents=True)
            counters = (build / "apmi.gcda", nested / "generator.gcda")
            preserved = (build / "apmi.gcno", build / "apmi.o", root / "outside.gcda")
            for path in counters + preserved:
                path.write_text("evidence", encoding="utf-8")
            self.assertEqual(support.remove_stale_counters(build), 2)
            self.assertTrue(all(not path.exists() for path in counters))
            self.assertTrue(all(path.read_text(encoding="utf-8") == "evidence" for path in preserved))
            self.assertEqual(support.remove_stale_counters(build), 0)
            self.assertEqual(support.remove_stale_counters(root / "absent"), 0)

    def test_refuses_a_symlinked_build_directory(self):
        with mock.patch.object(Path, "is_symlink", return_value=True):
            with self.assertRaisesRegex(RuntimeError, "build directory must not be a symlink"):
                support.remove_stale_counters(Path("build"))

    def test_symlinked_counter_prevents_all_deletion(self):
        with tempfile.TemporaryDirectory() as temporary:
            build = Path(temporary)
            ordinary, linked = build / "ordinary.gcda", build / "linked.gcda"
            for path in (ordinary, linked):
                path.write_text("evidence", encoding="utf-8")
            # Mock the link predicate so this contract also runs without Windows
            # symlink privileges, while exercising the real traversal and refusal.
            with mock.patch.object(Path, "is_symlink", lambda path: path == linked):
                with self.assertRaisesRegex(RuntimeError, "outside the native build tree"):
                    support.remove_stale_counters(build)
            self.assertTrue(ordinary.is_file())
            self.assertTrue(linked.is_file())

    def test_resolved_escape_prevents_all_deletion(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            build = root / "build"
            build.mkdir()
            counter = build / "escaped.gcda"
            outside = root / "outside.gcda"
            for path in (counter, outside):
                path.write_text("evidence", encoding="utf-8")
            resolve = Path.resolve
            with mock.patch.object(Path, "resolve", lambda path: outside if path == counter else resolve(path)):
                with self.assertRaisesRegex(RuntimeError, "outside the native build tree"):
                    support.remove_stale_counters(build)
            self.assertTrue(counter.is_file())
            self.assertTrue(outside.is_file())


class TestStartupHook(unittest.TestCase):
    def test_refuses_either_existing_startup_file_without_overwriting(self):
        for filename in (support.HOOK_MODULE + ".py", support.HOOK_NAME):
            with self.subTest(filename=filename), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                existing = root / filename
                existing.write_text("preexisting", encoding="utf-8")
                with self.assertRaisesRegex(RuntimeError, "already exist"):
                    support.write_hook(root, root / "coverage.rc", root / "audit")
                self.assertEqual(existing.read_text(encoding="utf-8"), "preexisting")
                self.assertEqual(list(root.iterdir()), [existing])

    def test_refuses_a_dangling_startup_symlink(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            linked = root / support.HOOK_NAME
            with mock.patch.object(Path, "is_symlink", lambda path: path == linked):
                with self.assertRaisesRegex(RuntimeError, "already exist"):
                    support.write_hook(root, root / "coverage.rc", root / "audit")
            self.assertEqual(list(root.iterdir()), [])

    def test_child_uses_the_same_absolute_config_with_missing_or_foreign_environment(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            startup = root / "site-packages"
            startup.mkdir()
            audit = root / "audit"
            audit.mkdir()
            config = root / "coverage.rc"
            config.write_text("[run]\nparallel = True\n", encoding="utf-8")
            owned = support.write_hook(startup, config, audit)
            self.addCleanup(support.remove_hook, owned)
            # -I -S avoids any real site startup, installed coverage package, or
            # outer instrumentation. Only addsitedir executes the generated .pth.
            child = """
import json
import os
from pathlib import Path
import site
import sys
import types
coverage = types.ModuleType("coverage")
def process_startup():
    Path(sys.argv[2]).write_text(json.dumps({"config": os.environ["COVERAGE_PROCESS_START"]}), encoding="utf-8")
coverage.process_startup = process_startup
sys.modules["coverage"] = coverage
site.addsitedir(sys.argv[1])
"""
            try:
                for index, inherited in enumerate((None, str(root / "foreign.rc"))):
                    with self.subTest(inherited=inherited):
                        proof = root / ("proof-%d.json" % index)
                        environment = os.environ.copy()
                        environment.pop("COVERAGE_PROCESS_START", None)
                        environment.pop("PYTHONPATH", None)
                        if inherited is not None:
                            environment["COVERAGE_PROCESS_START"] = inherited
                        result = subprocess.run(
                            [sys.executable, "-I", "-S", "-c", child, str(startup), str(proof)],
                            cwd=root, env=environment, capture_output=True, text=True, check=False,
                        )
                        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                        self.assertEqual(json.loads(proof.read_text(encoding="utf-8")), {"config": str(config)})
                markers = [json.loads(path.read_text(encoding="utf-8")) for path in audit.glob("*.json")]
                self.assertEqual(len(markers), 2)
                self.assertTrue(all(isinstance(marker["pid"], int) and marker["argv"] for marker in markers))
            finally:
                support.remove_hook(owned)
            self.assertFalse(owned)

    def test_cleanup_refuses_changed_files_before_removing_any_owned_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            owned = support.write_hook(root, root / "coverage.rc", root / "audit")
            hook = root / support.HOOK_NAME
            original = owned[hook]
            hook.write_text("changed by someone else\n", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "changed; refusing to remove"):
                support.remove_hook(owned)
            self.assertEqual(len(owned), 2)
            self.assertTrue(all(path.is_file() for path in owned))
            self.assertEqual(hook.read_text(encoding="utf-8"), "changed by someone else\n")
            hook.write_text(original, encoding="utf-8")
            paths = list(owned)
            support.remove_hook(owned)
            self.assertFalse(owned)
            self.assertTrue(all(not path.exists() for path in paths))
            support.remove_hook(owned)

    def test_cleanup_refuses_a_replaced_symlink(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            owned = support.write_hook(root, root / "coverage.rc", root / "audit")
            hook = root / support.HOOK_NAME
            with mock.patch.object(Path, "is_symlink", lambda path: path == hook):
                with self.assertRaisesRegex(RuntimeError, "changed; refusing to remove"):
                    support.remove_hook(owned)
            self.assertTrue(all(path.is_file() for path in owned))
            support.remove_hook(owned)


class TestCoverageSummary(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        package = self.root / "SJARACNe"
        (package / "bin").mkdir(parents=True)
        (package / "src").mkdir()
        self.python = {"files": {}, "totals": {
            "num_statements": 8, "covered_lines": 3,
            "num_branches": 10, "covered_branches": 7, "percent_covered": 99.99,
        }}
        for name, statements, covered in (
                ("__init__.py", 0, 0), ("workflow.py", 4, 2),
                ("bin/__init__.py", 0, 0), ("bin/worker.py", 4, 1)):
            path = package / name
            path.write_text("# fixture\n", encoding="utf-8")
            self.python["files"][path.relative_to(self.root).as_posix()] = {
                "summary": {"num_statements": statements, "covered_lines": covered},
            }
        self.cpp = {"line_total": 10, "line_covered": 4, "branch_total": 8, "branch_covered": 3,
                    "files": [{"filename": "SJARACNe/src/main.cpp"}]}

    def test_percentages_use_lines_not_combined_line_and_branch_percent(self):
        result = support.summarize(self.python, self.cpp, self.root)
        self.assertEqual(result["python"], {
            "covered_lines": 3, "lines": 8, "line_percent": 37.5,
            "covered_branches": 7, "branches": 10,
        })
        self.assertEqual(result["cpp"], {
            "covered_lines": 4, "lines": 10, "line_percent": 40,
            "covered_branches": 3, "branches": 8,
        })

    def test_source_names_are_compared_by_canonical_path(self):
        record = self.python["files"].pop("SJARACNe/workflow.py")
        self.python["files"]["SJARACNe/bin/../workflow.py"] = record
        self.cpp["files"][0]["filename"] = "SJARACNe/src/../src/main.cpp"
        self.assertEqual(support.summarize(self.python, self.cpp, self.root)["python"]["line_percent"], 37.5)

    def test_rejects_omitted_unexecuted_or_foreign_python_modules(self):
        for failure in ("omitted", "unexecuted", "foreign"):
            with self.subTest(failure=failure):
                report = copy.deepcopy(self.python)
                if failure == "omitted":
                    report["files"].pop("SJARACNe/workflow.py")
                elif failure == "unexecuted":
                    report["files"]["SJARACNe/workflow.py"]["summary"]["covered_lines"] = 0
                else:
                    report["files"]["tests/not_runtime.py"] = {"summary": {"num_statements": 0, "covered_lines": 0}}
                with self.assertRaises(RuntimeError):
                    support.summarize(report, self.cpp, self.root)

    def test_rejects_empty_denominators_or_no_executed_lines(self):
        for language, field, value in (
                ("python", "num_statements", 0), ("python", "covered_lines", 0),
                ("cpp", "line_total", 0), ("cpp", "line_covered", 0), ("cpp", "files", [])):
            with self.subTest(language=language, field=field):
                python, cpp = copy.deepcopy(self.python), copy.deepcopy(self.cpp)
                if language == "python":
                    python["totals"][field] = value
                else:
                    cpp[field] = value
                with self.assertRaisesRegex(RuntimeError, "Empty .* coverage denominator"):
                    support.summarize(python, cpp, self.root)

    def test_rejects_generated_or_out_of_scope_cpp_sources(self):
        for filename in ("SJARACNe/src/apmi_pchip_data.h", "tools/apmi_null_generator.cpp",
                         "SJARACNe/src/../../outside.cpp"):
            with self.subTest(filename=filename):
                self.cpp["files"] = [{"filename": filename}]
                with self.assertRaisesRegex(RuntimeError, "Unexpected C[+][+] coverage source"):
                    support.summarize(self.python, self.cpp, self.root)


class TestPersistedTaskCoverage(unittest.TestCase):
    def test_requires_saved_lines_for_each_cwl_and_wdl_task_process(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            data, audit, package = root / "data", root / "audit", root / "SJARACNe"
            data.mkdir()
            audit.mkdir()
            records = {}
            expected = {}
            loaded = []
            tasks = [("wdl:" + phase, "SJARACNe.bin.wdl_workflow", phase, "bin/wdl_workflow.py")
                     for phase in ("prepare", "infer", "consensus")]
            tasks += [(name, name, None, filename) for name, filename in (
                ("SJARACNe.bin.QC_input", "bin/QC_input.py"),
                ("SJARACNe.bin.create_consensus_network", "bin/create_consensus_network.py"))]
            for pid, (label, module, phase, filename) in enumerate(tasks, 100):
                argv = ["python", "-m", module] + ([phase] if phase is not None else [])
                (audit / (str(pid) + ".json")).write_text(json.dumps({"pid": pid, "argv": argv}), encoding="utf-8")
                dataset = data / (".coverage.host.%s%d.token" % ("pid" if pid % 2 else "", pid))
                dataset.write_text("mock database", encoding="utf-8")
                records[str(dataset)] = str(package / filename)
                expected[label] = 1

            def make_data(basename):
                result = mock.Mock()
                result.lines.side_effect = lambda path: [1] if records[basename] == path else []
                loaded.append(result)
                return result

            fake = types.ModuleType("coverage")
            fake.CoverageData = mock.Mock(side_effect=make_data)
            with mock.patch.dict(sys.modules, {"coverage": fake}):
                self.assertEqual(support.persisted_task_coverage(data, audit, package), expected)
                self.assertEqual(fake.CoverageData.call_count, len(tasks))
                for result in loaded:
                    result.read.assert_called_once_with()
                    self.assertEqual(result.lines.call_count, 1)
                records.clear()
                fake.CoverageData.side_effect = lambda basename: types.SimpleNamespace(read=lambda: None, lines=lambda path: [])
                with self.assertRaisesRegex(RuntimeError, "Missing persisted task coverage:.*wdl:prepare"):
                    support.persisted_task_coverage(data, audit, package)

    def test_startup_markers_without_matching_data_do_not_count(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "audit").mkdir()
            (root / "data").mkdir()
            (root / "audit" / "123.json").write_text(
                json.dumps({"pid": 123, "argv": ["python", "-m", "SJARACNe.bin.wdl_workflow", "prepare"]}),
                encoding="utf-8",
            )
            (root / "data" / ".coverage.host.999.token").write_text("foreign process", encoding="utf-8")
            fake = types.ModuleType("coverage")
            fake.CoverageData = mock.Mock()
            with mock.patch.dict(sys.modules, {"coverage": fake}):
                with self.assertRaisesRegex(RuntimeError, "Missing persisted task coverage"):
                    support.persisted_task_coverage(root / "data", root / "audit", root / "SJARACNe")
            fake.CoverageData.assert_not_called()


class TestCoverageSession(unittest.TestCase):
    def test_prepare_rejects_non_disposable_or_out_of_venv_startup_locations(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            checkout, output, external = root / "source", root / "evidence", root / "environment"
            checkout.mkdir()
            output.mkdir()
            for prefix, base, purelib, diagnostic in (
                    (external, external, external / "site-packages", "disposable virtual environment"),
                    (checkout / ".venv", root / "base", checkout / ".venv/site-packages", "outside the checkout"),
                    (external, root / "base", root / "foreign-site-packages", "outside the virtual environment")):
                with self.subTest(prefix=prefix, purelib=purelib), ExitStack() as stack:
                    stack.enter_context(mock.patch.object(support.sys, "prefix", str(prefix)))
                    stack.enter_context(mock.patch.object(support.sys, "base_prefix", str(base)))
                    stack.enter_context(mock.patch.object(support.sysconfig, "get_path", return_value=str(purelib)))
                    versions = stack.enter_context(mock.patch.object(support.importlib.metadata, "version"))
                    with self.assertRaisesRegex(RuntimeError, diagnostic):
                        support.CoverageSession(checkout, output).prepare()
                    versions.assert_not_called()
                    self.assertFalse((output / "coverage").exists())

    def test_prepare_requires_matching_gcc_and_gcov_before_creating_artifacts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            checkout, output, external = root / "source", root / "evidence", root / "environment"
            checkout.mkdir()
            output.mkdir()
            with ExitStack() as stack:
                stack.enter_context(mock.patch.object(support.sys, "prefix", str(external)))
                stack.enter_context(mock.patch.object(support.sys, "base_prefix", str(root / "base")))
                stack.enter_context(mock.patch.object(support.sysconfig, "get_path", return_value=str(external / "site-packages")))
                stack.enter_context(mock.patch.object(support.importlib.metadata, "version", return_value="mock version"))
                compiler = stack.enter_context(mock.patch.object(support.subprocess, "check_output", side_effect=[
                    "12.3.0\n", "gcov (Ubuntu 11.4.0-1ubuntu1) 11.4.0\nCopyright notice\n"]))
                hook = stack.enter_context(mock.patch.object(support, "write_hook"))
                with self.assertRaisesRegex(RuntimeError, "GCC and gcov versions must match"):
                    support.CoverageSession(checkout, output).prepare()
                self.assertEqual(compiler.call_args_list, [
                    mock.call(["g++", "-dumpfullversion"], text=True),
                    mock.call(["gcov", "--version"], text=True),
                ])
                hook.assert_not_called()
                self.assertFalse((output / "coverage").exists())

    def test_prepare_uses_an_absolute_config_and_records_its_instrumentation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            checkout, output, external = root / "source", root / "evidence", root / "environment"
            checkout.mkdir()
            output.mkdir()
            with ExitStack() as stack:
                stack.enter_context(mock.patch.object(support.sys, "prefix", str(external)))
                stack.enter_context(mock.patch.object(support.sys, "base_prefix", str(root / "base")))
                purelib = external / "site-packages"
                stack.enter_context(mock.patch.object(support.sysconfig, "get_path", return_value=str(purelib)))
                stack.enter_context(mock.patch.object(support.importlib.metadata, "version", side_effect=lambda name: "mock " + name))
                stack.enter_context(mock.patch.object(support.subprocess, "check_output", side_effect=[
                    "12.3.0\n", "gcov (Ubuntu 12.3.0-1ubuntu1) 12.3.0\nCopyright notice\n"]))
                cleanup = stack.enter_context(mock.patch.object(support, "remove_stale_counters", return_value=3))
                hook = stack.enter_context(mock.patch.object(support, "write_hook", return_value={}))
                session = support.CoverageSession(checkout, output)
                session.prepare()
                self.assertTrue(session.config.is_absolute())
                self.assertIn("data_file = " + str(session.data / ".coverage"), session.config.read_text(encoding="utf-8"))
                self.assertIn("source = " + str(checkout / "SJARACNe"), session.config.read_text(encoding="utf-8"))
                cleanup.assert_called_once_with(checkout / "SJARACNe/build")
                hook.assert_called_once_with(purelib, session.config, session.audit)
                self.assertEqual(session.details["stale_counters_removed"], 3)
                self.assertEqual(session.details["native_build_overrides"], list(support.NATIVE_OVERRIDES))
                self.assertEqual(session.details["versions"], {"coverage": "mock coverage", "gcovr": "mock gcovr"})
                session.close()

    def test_report_requires_the_owned_hook_to_be_closed(self):
        fake = types.ModuleType("coverage")
        fake.Coverage = mock.Mock()
        with tempfile.TemporaryDirectory() as temporary, mock.patch.dict(sys.modules, {"coverage": fake}):
            session = support.CoverageSession(temporary, temporary)
            session.owned = {Path(temporary) / "hook.pth": "owned"}
            with self.assertRaisesRegex(RuntimeError, "Remove coverage startup hooks before reporting"):
                session.report()
        fake.Coverage.assert_not_called()


if __name__ == "__main__":
    unittest.main()
