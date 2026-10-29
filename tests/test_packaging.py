"""Cheap packaging contracts; actual clean archive installs run separately."""

import ast
import importlib.util
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock
import zipfile

from setuptools import Distribution

from check_source_install import PACKAGE_ASSETS, PROBE, check_members, inspect_wheel


ROOT = Path(__file__).resolve().parents[1]


def load_setup():
    specification = importlib.util.spec_from_file_location(
        'sjaracne_packaging_contract', ROOT / 'setup.py')
    module = importlib.util.module_from_spec(specification)
    with mock.patch('setuptools.setup') as configure:
        specification.loader.exec_module(module)
    return module, configure.call_args.kwargs


class TestPackaging(unittest.TestCase):
    def test_public_setuptools_build_without_distutils_imports(self):
        tree = ast.parse((ROOT / 'setup.py').read_text(encoding='utf-8'))
        imports = [node.module for node in ast.walk(tree)
                   if isinstance(node, ast.ImportFrom)]
        self.assertIn('setuptools.command.build', imports)
        self.assertFalse(any(name and name.startswith('distutils') for name in imports))

    def test_native_build_precedes_python_package_collection(self):
        module, settings = load_setup()
        command = settings['cmdclass']['build'](Distribution())
        sequence = mock.Mock()
        with mock.patch.object(module.subprocess, 'check_call') as native, \
                mock.patch('setuptools.command.build.build.run', autospec=True) as package:
            sequence.attach_mock(native, 'native')
            sequence.attach_mock(package, 'package')
            command.run()
        self.assertEqual(sequence.mock_calls, [
            mock.call.native(['make', '-C', str(ROOT / 'SJARACNe')]),
            mock.call.package(command),
        ])

    def test_failed_native_build_cannot_publish_python_only_package(self):
        module, settings = load_setup()
        command = settings['cmdclass']['build'](Distribution())
        with mock.patch.object(module.subprocess, 'check_call',
                               side_effect=subprocess.CalledProcessError(2, ['make'])), \
                mock.patch('setuptools.command.build.build.run') as package:
            with self.assertRaises(subprocess.CalledProcessError):
                command.run()
            package.assert_not_called()

    def test_wheel_is_not_advertised_as_platform_independent(self):
        _, settings = load_setup()
        self.assertTrue(settings['distclass']().has_ext_modules())

    def test_native_and_workflow_assets_are_explicit_package_data(self):
        _, settings = load_setup()
        self.assertTrue(settings['include_package_data'])
        patterns = settings['package_data']['SJARACNe']
        for pattern in ('bin/*', 'config/*.json', 'config/*.txt', 'cwl/*.cwl',
                        'wdl/*.wdl', 'wdl/*.toml'):
            self.assertIn(pattern, patterns)
        self.assertEqual(settings['python_requires'], '>=3.7.6')

    def test_pep517_backend_and_source_build_inputs_are_declared(self):
        project = (ROOT / 'pyproject.toml').read_text(encoding='utf-8')
        self.assertIn('build-backend = "setuptools.build_meta"', project)
        self.assertIn('setuptools>=67.8', project)
        self.assertIn('wheel>=0.42', project)
        manifest = (ROOT / 'MANIFEST.in').read_text(encoding='utf-8')
        for source in ('include version.py', 'include SJARACNe/src/*',
                       'include SJARACNe/Makefile',
                       'include tools/apmi_null_generator.cpp'):
            self.assertIn(source, manifest)

    def test_archive_inspection_rejects_missing_or_retired_assets(self):
        check_members(PACKAGE_ASSETS, PACKAGE_ASSETS)
        with self.assertRaisesRegex(RuntimeError, 'missing required assets'):
            check_members(PACKAGE_ASSETS[1:], PACKAGE_ASSETS)
        with self.assertRaisesRegex(RuntimeError, 'Retired runner assets'):
            check_members(PACKAGE_ASSETS + ('old/config_cwlexec.json',), PACKAGE_ASSETS)

    def test_wheel_inspection_requires_native_platform_metadata(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'package.whl'
            for pure, tag, valid in ((False, 'cp314-cp314-linux_x86_64', True),
                                     (True, 'py3-none-any', False),
                                     (False, 'py3-none-any', False)):
                with zipfile.ZipFile(path, 'w') as wheel:
                    for name in PACKAGE_ASSETS:
                        wheel.writestr(name, '')
                    wheel.writestr('sjaracne-0.2.1.dist-info/WHEEL',
                                   'Root-Is-Purelib: {}\nTag: {}\n'.format(
                                       str(pure).lower(), tag))
                if valid:
                    self.assertEqual(inspect_wheel(path), [tag])
                else:
                    with self.assertRaises(RuntimeError):
                        inspect_wheel(path)

    def test_installed_runtime_probe_is_valid_python_and_writes_json_newline(self):
        compile(PROBE, '<installed-runtime-probe>', 'exec')
        tree = ast.parse(PROBE)
        writes = [node for node in ast.walk(tree) if isinstance(node, ast.Call)
                  and isinstance(node.func, ast.Attribute) and node.func.attr == 'write_text']
        self.assertEqual(len(writes), 1)
        self.assertEqual(writes[0].args[0].right.value, '\n')


if __name__ == '__main__':
    unittest.main()
