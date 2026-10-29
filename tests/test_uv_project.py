"""Offline contracts for the locked project and native-build integration.

These tests inspect committed configuration. They do not install packages,
contact an index, or require uv in the runtime used by SJARACNe.
"""

import ast
from pathlib import Path
import tomllib
import unittest

from packaging.markers import default_environment
from packaging.requirements import Requirement

from check_uv_install import RUNTIME_PROBE, package_was_built


ROOT = Path(__file__).resolve().parents[1]


def project_configuration():
    return tomllib.loads((ROOT / 'pyproject.toml').read_text(encoding='utf-8'))


class TestUvProject(unittest.TestCase):
    def test_standard_metadata_is_centralized_without_changing_pip_compatibility(self):
        configuration = project_configuration()
        project = configuration['project']
        self.assertEqual(project['name'].lower(), 'sjaracne')
        self.assertEqual(project['requires-python'], '>=3.7.6')
        self.assertEqual(project['scripts']['sjaracne'], 'SJARACNe.sjaracne:main')
        self.assertEqual(project['dynamic'], ['version'])
        self.assertEqual(configuration['tool']['setuptools']['dynamic']['version'],
                         {'attr': 'version.__version__'})
        dependencies = {Requirement(item).name.lower() for item in project['dependencies']}
        self.assertEqual(dependencies, {'numpy', 'pandas', 'scipy'})
        self.assertEqual({Requirement(item).name.lower()
                          for item in project['optional-dependencies']['cwl']}, {'cwltool'})
        setup = ast.parse((ROOT / 'setup.py').read_text(encoding='utf-8'))
        calls = [node for node in ast.walk(setup) if isinstance(node, ast.Call)
                 and isinstance(node.func, ast.Name) and node.func.id == 'setup']
        self.assertEqual(len(calls), 1)
        metadata = {'name', 'version', 'install_requires', 'extras_require',
                    'python_requires', 'entry_points'}
        self.assertFalse(metadata.intersection(keyword.arg for keyword in calls[0].keywords))

    def test_preferred_python_and_uv_are_pinned_without_becoming_runtime_dependencies(self):
        self.assertEqual((ROOT / '.python-version').read_text(encoding='utf-8').strip(), '3.14')
        configuration = project_configuration()
        self.assertEqual(configuration['tool']['uv']['required-version'], '==0.11.16')
        for dependency in configuration['project']['dependencies']:
            self.assertNotEqual(Requirement(dependency).name.lower(), 'uv')
        ignored = (ROOT / '.gitignore').read_text(encoding='utf-8').splitlines()
        self.assertNotIn('.python-version', ignored)
        self.assertNotIn('uv.lock', ignored)

    def test_locked_constraints_retain_both_validated_numerical_stacks(self):
        constraints = [Requirement(item) for item in
                       project_configuration()['tool']['uv']['constraint-dependencies']]
        expected = {
            '3.11': {'numpy': '1.26.4', 'pandas': '2.2.3', 'scipy': '1.13.1',
                     'cwltool': '3.1.20240508115724'},
            '3.14': {'numpy': '2.4.6', 'pandas': '2.3.3', 'scipy': '1.17.1',
                     'cwltool': '3.3.20260925135507'},
        }
        for python, versions in expected.items():
            environment = default_environment()
            environment.update(python_version=python, python_full_version=python + '.0')
            for package, version in versions.items():
                matches = [item for item in constraints if item.name.lower() == package
                           and (item.marker is None or item.marker.evaluate(environment))]
                self.assertEqual(len(matches), 1, (python, package, matches))
                self.assertEqual(str(matches[0].specifier), '==' + version)

    def test_lock_contains_the_validated_versions_and_no_editable_external_source(self):
        lock = tomllib.loads((ROOT / 'uv.lock').read_text(encoding='utf-8'))
        packages = lock['package']
        versions = {}
        for package in packages:
            # uv need not resolve the dynamic version of the editable root.
            if 'version' in package:
                versions.setdefault(package['name'], set()).add(package['version'])
        for name, expected in {
            'numpy': {'1.26.4', '2.4.6'}, 'pandas': {'2.2.3', '2.3.3'},
            'scipy': {'1.13.1', '1.17.1'},
            'cwltool': {'3.1.20240508115724', '3.3.20260925135507'},
            'pip': {'26.2.2'}, 'setuptools': {'84.2.1'},
            'wheel': {'0.48.0'}, 'build': {'1.6.1'},
        }.items():
            self.assertEqual(versions[name], expected)
        local = [item for item in packages if 'editable' in item.get('source', {})]
        self.assertEqual(len(local), 1)
        self.assertEqual(local[0]['name'], 'sjaracne')
        self.assertEqual(local[0]['source'], {'editable': '.'})

    def test_native_build_inputs_are_covered_by_file_content_cache_keys(self):
        keys = project_configuration()['tool']['uv']['cache-keys']
        covered = set()
        for key in keys:
            if 'file' in key:
                covered.update(path.resolve() for path in ROOT.glob(key['file']) if path.is_file())
        expected = {ROOT / name for name in (
            'pyproject.toml', 'setup.py', 'version.py', 'MANIFEST.in',
            'SJARACNe/Makefile', 'tools/apmi_null_generator.cpp',
        )}
        expected.update(path for path in (ROOT / 'SJARACNe/src').rglob('*')
                        if path.is_file() and path.suffix in ('.cpp', '.h'))
        expected.update(path for path in (ROOT / 'SJARACNe/config').glob('*') if path.is_file())
        missing = sorted(str(path.relative_to(ROOT)) for path in expected
                         if path.resolve() not in covered)
        self.assertEqual(missing, [], 'Native edits must invalidate the cached package build')

    def test_make_recompiles_all_native_objects_after_header_or_recipe_changes(self):
        makefile = (ROOT / 'SJARACNe/Makefile').read_text(encoding='utf-8')
        self.assertRegex(makefile, r'(?m)^HEADERS\s*:=.*find.*\$\(SRCDIR\)')
        self.assertRegex(makefile,
                         r'(?m)^\$\(OBJECTS\)\s+\$\(NULL_OBJECT\):\s+\$\(HEADERS\)\s+Makefile\s*$')

    def test_uv_runtime_probe_is_valid_and_build_detection_is_package_specific(self):
        compile(RUNTIME_PROBE, '<uv-runtime-probe>', 'exec')
        self.assertTrue(package_was_built('   Building sjaracne @ file:///tmp/source\n'))
        self.assertTrue(package_was_built('DEBUG Built `SJARACNe @ file:///tmp/source`\n'))
        self.assertFalse(package_was_built('Building setuptools\nInstalled 1 package\n'))


if __name__ == '__main__':
    unittest.main()
