"""Cheap contract checks for the ordinary regression workflow configuration."""

from pathlib import Path
import tomllib
import unittest

from ruamel.yaml import YAML


ROOT = Path(__file__).resolve().parents[1]


class TestCIWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workflow = YAML(typ='safe').load(
            (ROOT / '.github/workflows/regression.yml').read_text(encoding='utf-8'))
        cls.project = tomllib.loads((ROOT / 'pyproject.toml').read_text(encoding='utf-8'))

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
        self.assertEqual(len(actions), 4)
        for step in actions:
            self.assertRegex(step['uses'], r'^(?:actions|astral-sh)/[a-z-]+@[0-9a-f]{40}$')
        checkout = next(step for step in actions if step['uses'].startswith('actions/checkout@'))
        self.assertIs(checkout['with']['persist-credentials'], False)

    def test_checks_python_311_compatibility_and_314_independently(self):
        job = self.workflow['jobs']['regression']
        strategy = job['strategy']
        self.assertIs(strategy['fail-fast'], False)
        self.assertNotIn('continue-on-error', job)
        self.assertEqual(strategy['matrix']['include'], [
            {'python-version': '3.11'},
            {'python-version': '3.14'},
        ])
        setup = next(step for step in job['steps']
                     if step.get('uses', '').startswith('actions/setup-python@'))
        self.assertEqual(setup['with']['python-version'], '${{ matrix.python-version }}')
        self.assertNotIn('cache', setup['with'])
        install = next(step for step in job['steps']
                       if step.get('name') == 'Install locked regression environment and this checkout')
        self.assertIn("uv sync --locked --python '${{ matrix.python-version }}' --group regression --extra cwl", install['run'])
        self.assertNotIn('make -C SJARACNe', install['run'])
        system = next(step for step in job['steps']
                      if step.get('name') == 'Install system prerequisites')
        self.assertIn('build-essential', system['run'])
        baseline = '\n'.join(self.project['tool']['uv']['constraint-dependencies'])
        for pin in ('numpy==1.26.4', 'pandas==2.2.3', 'scipy==1.13.1',
                    'cwltool==3.1.20240508115724', 'numpy==2.4.6',
                    'pandas==2.3.3', 'scipy==1.17.1', 'cwltool==3.3.20260925135507'):
            self.assertIn(pin, baseline)

    def test_pins_uv_and_prepares_environment_once_outside_checkout(self):
        job = self.workflow['jobs']['regression']
        steps = job['steps']
        setup = next(step for step in steps
                     if step.get('uses', '').startswith('astral-sh/setup-uv@'))
        self.assertEqual(setup['with']['version'], '0.11.16')
        self.assertIs(setup['with']['enable-cache'], False)
        self.assertEqual(job['env']['UV_PYTHON_DOWNLOADS'], '0')
        self.assertEqual(job['env']['UV_PYTHON'], '${{ matrix.python-version }}')
        self.assertIn('UV_PROJECT_ENVIRONMENT=%s/sjaracne-env-py%s', steps[0]['run'])
        self.assertIn('UV_CACHE_DIR=%s/sjaracne-uv-cache', steps[0]['run'])
        self.assertIn('$RUNNER_TEMP', steps[0]['run'])
        commands = '\n'.join(step.get('run', '') for step in steps)
        self.assertEqual(commands.count('uv sync '), 1)
        self.assertNotIn('uv lock --upgrade', commands)
        self.assertIn('uv-version.txt', commands)
        self.assertIn('sha256sum uv.lock pyproject.toml', commands)

    def test_checks_fresh_source_install_with_pinned_build_tools(self):
        steps = self.workflow['jobs']['regression']['steps']
        build = next(step for step in steps
                     if step.get('name') == 'Verify fresh installation from the source archive')
        self.assertEqual(build['run'],
                         'uv run --locked --no-sync python tests/check_source_install.py '
                         '--output-dir "$BOOTSTRAP_OUTPUT/source-install" '
                         '--requirements "$BOOTSTRAP_OUTPUT/locked-requirements.txt" '
                         '--build-requirements "$BOOTSTRAP_OUTPUT/locked-requirements.txt"')
        self.assertNotIn('continue-on-error', build)
        regression = next(step for step in steps
                          if 'python tests/run_regression.py' in step.get('run', ''))
        self.assertLess(steps.index(build), steps.index(regression))
        install = next(step for step in steps
                       if step.get('name') == 'Install locked regression environment and this checkout')
        self.assertIn('uv export --locked --group regression --extra cwl --no-emit-project --no-hashes',
                      install['run'])
        requirements = '\n'.join(self.project['dependency-groups']['build'])
        for package in ('pip', 'setuptools', 'wheel', 'build'):
            self.assertRegex(requirements, r'(?m)^' + package + r'==[0-9][^\s]*$')

    def test_requires_real_runner_and_preserves_short_lived_evidence(self):
        job = self.workflow['jobs']['regression']
        # runner context is unavailable in job-level env. Paths are initialized
        # by a shell step, where RUNNER_TEMP exists, before subsequent steps.
        self.assertTrue(all('${{ runner.' not in str(value) for value in job['env'].values()))
        steps = job['steps']
        self.assertIn('$GITHUB_ENV', steps[0]['run'])
        self.assertEqual(steps[0]['run'].count('${{ matrix.python-version }}'), 3)
        runs = [step for step in steps if 'python tests/run_regression.py' in step.get('run', '')]
        self.assertEqual(len(runs), 1)
        self.assertTrue(runs[0]['run'].startswith('uv run --locked --no-sync '))
        self.assertIn('--output-dir "$REGRESSION_OUTPUT"', runs[0]['run'])
        self.assertLessEqual(runs[0]['timeout-minutes'], 15)
        self.assertNotIn('continue-on-error', runs[0])
        artifact = next(step for step in steps if step.get('uses', '').startswith('actions/upload-artifact@'))
        self.assertEqual(artifact['if'], 'always()')
        self.assertEqual(artifact['with']['retention-days'], 7)
        self.assertIn('${{ env.REGRESSION_OUTPUT }}', artifact['with']['path'])
        self.assertIn('${{ matrix.python-version }}', artifact['with']['name'])

    def test_checks_uv_offline_install_and_cache_invalidation_before_regression(self):
        steps = self.workflow['jobs']['regression']['steps']
        probe = next(step for step in steps if step.get('name') ==
                     'Verify uv offline installation and native cache invalidation')
        self.assertEqual(probe['run'],
                         'uv run --locked --no-sync python tests/check_uv_install.py '
                         '--uv "$(command -v uv)" '
                         '--python "$UV_PROJECT_ENVIRONMENT/bin/python" '
                         '--cache-dir "$RUNNER_TEMP/uv-install-cache-py${{ matrix.python-version }}" '
                         '--output-dir "$BOOTSTRAP_OUTPUT/uv-install"')
        # The probe must warm its own new cache before testing offline installs,
        # rather than silently inheriting CI's pre-existing cache contents.
        self.assertNotIn('--cache-dir "$UV_CACHE_DIR"', probe['run'])
        self.assertNotIn('continue-on-error', probe)
        source = next(step for step in steps if step.get('name') ==
                      'Verify fresh installation from the source archive')
        sprocket = next(step for step in steps if step.get('name') ==
                        'Install pinned Sprocket for WDL 1.3 regression')
        self.assertLess(steps.index(source), steps.index(probe))
        self.assertLess(steps.index(probe), steps.index(sprocket))

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
                          if 'python tests/run_regression.py' in step.get('run', ''))
        self.assertLess(steps.index(install), steps.index(regression))
