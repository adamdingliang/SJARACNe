"""Cheap packaging contracts; actual clean archive installs run separately."""

import ast
import importlib.util
from pathlib import Path
import subprocess
import tempfile
import tomllib
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
        distribution = Distribution({'packages': ['SJARACNe'],
                                     'cmdclass': settings['cmdclass']})
        command = settings['cmdclass']['build'](distribution)
        self.assertEqual(command.get_sub_commands()[0], 'build_native')
        self.assertLess(command.get_sub_commands().index('build_native'),
                        command.get_sub_commands().index('build_py'))
        sequence = []
        with mock.patch.object(command, 'run_command', side_effect=sequence.append):
            command.run()
        self.assertEqual(sequence[0], 'build_native')

    def test_regular_and_editable_native_commands_invoke_make_once(self):
        for editable in (False, True):
            module, settings = load_setup()
            command = settings['cmdclass']['build_native'](Distribution())
            command.editable_mode = editable
            with mock.patch.object(module.platform, 'system', return_value='Linux'), \
                    mock.patch.object(module.shutil, 'which', side_effect=lambda name: '/usr/bin/' + name), \
                    mock.patch.object(module.subprocess, 'check_call') as native:
                command.run()
            native.assert_called_once_with(['make', '-C', str(ROOT / 'SJARACNe')])

    def test_failed_native_build_cannot_publish_python_only_package(self):
        module, settings = load_setup()
        for editable in (False, True):
            command = settings['cmdclass']['build_native'](Distribution())
            command.editable_mode = editable
            with mock.patch.object(module.platform, 'system', return_value='Linux'), \
                    mock.patch.object(module.shutil, 'which', return_value='/available'), \
                    mock.patch.object(module.subprocess, 'check_call',
                                      side_effect=subprocess.CalledProcessError(2, ['make'])):
                with self.assertRaisesRegex(module.SetupError, 'native compilation failed'):
                    command.run()

    def test_missing_make_or_compiler_has_actionable_ubuntu_wsl_guidance(self):
        for absent in ('make', 'g++'):
            module, settings = load_setup()
            command = settings['cmdclass']['build_native'](Distribution())
            with mock.patch.object(module.platform, 'system', return_value='Linux'), \
                    mock.patch.object(module.shutil, 'which',
                                      side_effect=lambda name: None if name == absent else '/available'), \
                    mock.patch.object(module.subprocess, 'check_call') as native:
                with self.assertRaises(module.SetupError) as caught:
                    command.run()
            message = str(caught.exception)
            for expected in (absent, 'build-essential', 'WSL', 'sudo apt update',
                             'uv sync --locked', 'does not install system packages'):
                self.assertIn(expected, message)
            native.assert_not_called()

    def test_windows_and_msys2_fail_with_wsl_guidance_before_tool_lookup(self):
        for platform in ('Windows', 'MSYS_NT-10.0-26100'):
            module, settings = load_setup()
            command = settings['cmdclass']['build_native'](Distribution())
            with mock.patch.object(module.platform, 'system', return_value=platform), \
                    mock.patch.object(module.shutil, 'which') as which, \
                    mock.patch.object(module.subprocess, 'check_call') as native:
                with self.assertRaisesRegex(module.SetupError, 'Ubuntu/WSL'):
                    command.run()
            which.assert_not_called()
            native.assert_not_called()

    def test_metadata_loading_and_source_discovery_do_not_compile(self):
        with mock.patch('subprocess.check_call') as native, \
                mock.patch('shutil.which') as which:
            module, settings = load_setup()
            distribution = Distribution({'cmdclass': settings['cmdclass']})
            command = settings['cmdclass']['build_native'](distribution)
            outputs = command.get_outputs()
            mapping = command.get_output_mapping()
            sources = command.get_source_files()
        native.assert_not_called()
        which.assert_not_called()
        self.assertEqual(len(outputs), 2)
        self.assertEqual(set(mapping), set(outputs))
        self.assertEqual(set(mapping.values()), set(module.NATIVE_PROGRAMS))
        for source in ('SJARACNe/Makefile', 'tools/apmi_null_generator.cpp',
                       'SJARACNe/src/main.cpp', 'SJARACNe/src/param.h'):
            self.assertIn(str(Path(source)), sources)

    def test_build_launch_failure_is_not_silently_ignored(self):
        module, settings = load_setup()
        command = settings['cmdclass']['build_native'](Distribution())
        with mock.patch.object(module.platform, 'system', return_value='Linux'), \
                mock.patch.object(module.shutil, 'which', return_value='/available'), \
                mock.patch.object(module.subprocess, 'check_call',
                                  side_effect=FileNotFoundError('make disappeared')):
            with self.assertRaisesRegex(module.SetupError, 'cannot complete'):
                command.run()

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
        project = tomllib.loads((ROOT / 'pyproject.toml').read_text(encoding='utf-8'))
        self.assertEqual(project['project']['requires-python'], '>=3.7.6')

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

    def test_fresh_pip_editable_install_constrains_and_installs_its_runtime(self):
        tree = ast.parse((ROOT / 'tests/check_source_install.py').read_text(encoding='utf-8'))
        commands = [node for node in ast.walk(tree) if isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Name) and node.func.id == 'run'
                    and node.args and isinstance(node.args[0], ast.Constant)
                    and node.args[0].value == 'fresh-editable-install']
        self.assertEqual(len(commands), 1)
        arguments = commands[0].args[1].elts
        literals = [item.value for item in arguments if isinstance(item, ast.Constant)]
        self.assertNotIn('--no-deps', literals)
        constraint = next(index for index, item in enumerate(arguments)
                          if isinstance(item, ast.Constant) and item.value == '--constraint')
        self.assertIsInstance(arguments[constraint + 1], ast.Name)
        self.assertEqual(arguments[constraint + 1].id, 'requirements')
        self.assertTrue(any(isinstance(item, ast.BinOp) and isinstance(item.op, ast.Add)
                            and isinstance(item.right, ast.Constant)
                            and item.right.value == '[cwl]' for item in arguments))


if __name__ == '__main__':
    unittest.main()
