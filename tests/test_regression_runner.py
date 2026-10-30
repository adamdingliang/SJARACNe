import unittest
from pathlib import Path
import tempfile

from run_regression import (strict_success, collect_test_ids, build_command,
                            verify_installed_runtime, source_hashes)


class TestStrictRunner(unittest.TestCase):
    def test_packaging_controls_are_part_of_the_frozen_sources(self):
        hashes = source_hashes()
        for name in ('setup.py', 'pyproject.toml', 'uv.lock', '.python-version',
                     'version.py', 'MANIFEST.in'):
            self.assertIn(name, hashes)

    def test_native_build_is_forced_including_both_targets(self):
        command = build_command(2)
        self.assertEqual(command[0], 'make')
        self.assertIn('-B', command)
        self.assertEqual(command[-1], 'all')
        self.assertEqual(command[command.index('-j') + 1], 2)

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
