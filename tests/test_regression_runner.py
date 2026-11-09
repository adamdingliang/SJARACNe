import unittest
from pathlib import Path
import tempfile

from run_regression import (strict_success, collect_test_ids, build_command,
                            verify_installed_runtime, source_hashes, ROOT)
from coverage_support import NATIVE_OVERRIDES


class TestStrictRunner(unittest.TestCase):
    def test_packaging_controls_are_part_of_the_frozen_sources(self):
        hashes = source_hashes()
        for name in ('setup.py', 'pyproject.toml', 'uv.lock', '.python-version',
                     'version.py', 'MANIFEST.in'):
            self.assertIn(name, hashes)

    def test_native_build_is_forced_including_both_targets(self):
        command = build_command(2)
        self.assertEqual(command, ['make', '-B', '-C', ROOT / 'SJARACNe', '-j', 2, 'all'])
        self.assertEqual(build_command(2, coverage=False), command)

    def test_coverage_build_preserves_forced_targets_and_compile_and_link_flags(self):
        ordinary = build_command(2)
        instrumented = build_command(2, coverage=True)
        self.assertEqual(instrumented, ordinary[:-1] + list(NATIVE_OVERRIDES) + ['all'])
        self.assertIn('CC=g++', instrumented)
        compiler_flags = next(value.split('=', 1)[1] for value in instrumented
                              if isinstance(value, str) and value.startswith('CFLAGS='))
        linker_flags = next(value.split('=', 1)[1] for value in instrumented
                            if isinstance(value, str) and value.startswith('LIB='))
        self.assertIn('-O0', compiler_flags.split())
        self.assertNotIn('-O3', compiler_flags.split())
        self.assertIn('--coverage', compiler_flags.split())
        self.assertIn('-fprofile-abs-path', compiler_flags.split())
        self.assertIn('-fprofile-update=atomic', compiler_flags.split())
        self.assertIn('--coverage', linker_flags.split())
        self.assertEqual(build_command(2), ordinary)

    def test_foreign_or_missing_installed_runtime_cannot_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, installed = root / 'source', root / 'installed'
            source.mkdir()
            installed.mkdir()
            (source / '__init__.py').write_text('current', encoding='utf-8')
            with self.assertRaisesRegex(RuntimeError, 'Installed SJARACNe differs'):
                verify_installed_runtime(installed, source)
            (installed / '__init__.py').write_text('older', encoding='utf-8')
            with self.assertRaisesRegex(RuntimeError, 'Installed SJARACNe differs'):
                verify_installed_runtime(installed, source)
            (installed / '__init__.py').write_text('current', encoding='utf-8')
            self.assertEqual(set(verify_installed_runtime(installed, source)['sha256']), {'__init__.py'})

    def test_rejects_empty_and_incomplete_runs(self):
        result = unittest.TestResult()
        self.assertFalse(strict_success(result, 0))
        self.assertFalse(strict_success(result, 1))

    def test_rejects_skips_and_expected_failures(self):
        result = unittest.TestResult()
        result.testsRun = 1
        self.assertTrue(strict_success(result, 1))
        result.skipped = [("example", "native missing")]
        self.assertFalse(strict_success(result, 1))
        result.skipped = []
        result.expectedFailures = [("example", "not a full pass")]
        self.assertFalse(strict_success(result, 1))

    def test_nested_discovery_keeps_every_test_id(self):
        case = unittest.FunctionTestCase(lambda: None, description="example")
        suite = unittest.TestSuite([unittest.TestSuite([case])])
        self.assertEqual(list(collect_test_ids(suite)), [case.id()])
