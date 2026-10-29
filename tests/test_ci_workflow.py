"""Cheap contract checks for the ordinary regression workflow configuration."""

from pathlib import Path
import unittest

from ruamel.yaml import YAML


ROOT = Path(__file__).resolve().parents[1]


class TestCIWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workflow = YAML(typ='safe').load(
            (ROOT / '.github/workflows/regression.yml').read_text(encoding='utf-8'))

    def test_checks_pr_targets_and_pushes_on_both_branches(self):
        events = self.workflow['on']
        self.assertEqual(set(events), {'pull_request', 'push', 'workflow_dispatch'})
        for event in ('pull_request', 'push'):
            self.assertEqual(events[event], {'branches': ['develop', 'master']})

    def test_uses_read_only_standard_runner_and_pinned_actions(self):
        self.assertEqual(self.workflow['permissions'], {'contents': 'read'})
        job = self.workflow['jobs']['regression']
        self.assertEqual(job['runs-on'], 'ubuntu-24.04')
        self.assertLessEqual(job['timeout-minutes'], 20)
        actions = [step for step in job['steps'] if 'uses' in step]
        self.assertEqual(len(actions), 3)
        for step in actions:
            self.assertRegex(step['uses'], r'^actions/[a-z-]+@[0-9a-f]{40}$')
        checkout = next(step for step in actions if step['uses'].startswith('actions/checkout@'))
        self.assertIs(checkout['with']['persist-credentials'], False)

    def test_checks_python_311_compatibility_and_314_independently(self):
        job = self.workflow['jobs']['regression']
        strategy = job['strategy']
        self.assertIs(strategy['fail-fast'], False)
        self.assertNotIn('continue-on-error', job)
        self.assertEqual(strategy['matrix']['include'], [
            {'python-version': '3.11', 'requirements-file': 'tests/requirements-ci.txt'},
            {'python-version': '3.14', 'requirements-file': 'tests/requirements-ci-py314.txt'},
        ])
        setup = next(step for step in job['steps']
                     if step.get('uses', '').startswith('actions/setup-python@'))
        self.assertEqual(setup['with']['python-version'], '${{ matrix.python-version }}')
        for requirements in ('tests/requirements-ci.txt', 'tests/requirements-ci-py314.txt'):
            self.assertTrue((ROOT / requirements).is_file())
            self.assertIn(requirements, setup['with']['cache-dependency-path'])
        self.assertIn('tests/requirements-build.txt', setup['with']['cache-dependency-path'])
        self.assertIn('pyproject.toml', setup['with']['cache-dependency-path'])
        install = next(step for step in job['steps']
                       if step.get('name') == 'Install pinned CI baseline and this checkout')
        self.assertIn("-r '${{ matrix.requirements-file }}'", install['run'])
        self.assertIn('python -m pip install --upgrade -r tests/requirements-build.txt', install['run'])
        baseline = (ROOT / 'tests/requirements-ci.txt').read_text(encoding='utf-8')
        for pin in ('numpy==1.26.4', 'pandas==2.2.3', 'scipy==1.13.1',
                    'cwltool==3.1.20240508115724'):
            self.assertIn(pin, baseline)

    def test_checks_fresh_source_install_with_pinned_build_tools(self):
        steps = self.workflow['jobs']['regression']['steps']
        build = next(step for step in steps
                     if step.get('name') == 'Verify fresh installation from the source archive')
        self.assertEqual(build['run'],
                         'python tests/check_source_install.py '
                         '--output-dir "$BOOTSTRAP_OUTPUT/source-install" '
                         "--requirements '${{ matrix.requirements-file }}'")
        self.assertNotIn('continue-on-error', build)
        regression = next(step for step in steps
                          if step.get('run', '').startswith('python tests/run_regression.py'))
        self.assertLess(steps.index(build), steps.index(regression))
        requirements = (ROOT / 'tests/requirements-build.txt').read_text(encoding='utf-8')
        for package in ('pip', 'setuptools', 'wheel', 'build'):
            self.assertRegex(requirements, r'(?m)^' + package + r'==[0-9][^\s]*$')

    def test_requires_real_runner_and_preserves_short_lived_evidence(self):
        job = self.workflow['jobs']['regression']
        # runner context is unavailable in job-level env. Paths are initialized
        # by a shell step, where RUNNER_TEMP exists, before subsequent steps.
        self.assertTrue(all('${{ runner.' not in str(value) for value in job['env'].values()))
        steps = job['steps']
        self.assertIn('$GITHUB_ENV', steps[0]['run'])
        self.assertEqual(steps[0]['run'].count('${{ matrix.python-version }}'), 2)
        runs = [step for step in steps if step.get('run', '').startswith('python tests/run_regression.py')]
        self.assertEqual(len(runs), 1)
        self.assertIn('--output-dir "$REGRESSION_OUTPUT"', runs[0]['run'])
        self.assertLessEqual(runs[0]['timeout-minutes'], 15)
        self.assertNotIn('continue-on-error', runs[0])
        artifact = next(step for step in steps if step.get('uses', '').startswith('actions/upload-artifact@'))
        self.assertEqual(artifact['if'], 'always()')
        self.assertEqual(artifact['with']['retention-days'], 7)
        self.assertIn('${{ env.REGRESSION_OUTPUT }}', artifact['with']['path'])
        self.assertIn('${{ matrix.python-version }}', artifact['with']['name'])

    def test_installs_checksum_pinned_sprocket_before_regression(self):
        steps = self.workflow['jobs']['regression']['steps']
        candidates = [step for step in steps if 'sprocket/releases/download/' in step.get('run', '')]
        self.assertEqual(len(candidates), 1)
        install = candidates[0]
        command = install['run']
        self.assertIn('/v0.32.0/sprocket-v0.32.0-x86_64-unknown-linux-gnu.tar.gz', command)
        self.assertIn('27dbf45e3ceaad9e13d416fdf44e8741af1a98c6a516723b27057e695cab0756', command)
        self.assertIn('sha256sum --check --strict', command)
        self.assertIn('$GITHUB_PATH', command)
        self.assertIn('sprocket-version.txt', command)
        self.assertNotIn('continue-on-error', install)
        regression = next(step for step in steps
                          if step.get('run', '').startswith('python tests/run_regression.py'))
        self.assertLess(steps.index(install), steps.index(regression))
