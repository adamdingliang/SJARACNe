"""Dynamic-DPI wrapper, CWL wiring and small real-workflow regressions."""
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

from native_test_support import find_native_executable
from test_sjaracne_acceptance import write_multihub_dpi_inputs


ROOT = Path(__file__).resolve().parents[1]
NATIVE = find_native_executable()


def load_wrapper():
    spec = importlib.util.spec_from_file_location('dynamic_dpi_wrapper', ROOT / 'SJARACNe/sjaracne.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestDynamicDpiWrapper(unittest.TestCase):
    def run_wrapper(self, *extra, platform='local'):
        module = load_wrapper()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output, scratch = root / 'output', root / 'scratch'
            expression, hubs = write_multihub_dpi_inputs(root)
            args = ['sjaracne', platform, '--engine', 'cwl', '--network-type', 'TF',
                    '-e', str(expression), '-g', str(hubs), '-o', str(output), '-tmp', str(scratch)]
            if platform == 'lsf':
                (root / 'lsf.json').write_text('{}', encoding='ascii')
                args += ['-j', str(root / 'lsf.json')]
            args += list(extra)
            with mock.patch.object(module, 'run_shell_command_call') as launch, \
                    mock.patch.object(module, 'runtime_fingerprint', return_value={}), \
                    mock.patch.object(module, 'validate_cwl_outputs', return_value={}), \
                    mock.patch.object(module, 'require_tool', side_effect=lambda name: name):
                with mock.patch.object(sys, 'argv', args):
                    module.main()
                self.assertEqual(launch.call_count, 1)
                command = launch.call_args.args[0]
            return json.loads((output / 'sjaracne_workflow.yml').read_text(encoding='utf-8')), command

    def test_explicit_fixed_dpi_disables_new_dynamic_default(self):
        job, _ = self.run_wrapper('--fixed-dpi')
        self.assertFalse(job['dynamic_dpi'])
        self.assertEqual(job['tolerance'], 0)

    def test_current_default_and_explicit_dynamic_use_same_native_tuning(self):
        default, _ = self.run_wrapper()
        explicit, _ = self.run_wrapper('--dynamic-dpi')
        self.assertTrue(default['dynamic_dpi'])
        self.assertEqual(default['dynamic_dpi'], explicit['dynamic_dpi'])
        self.assertIsNone(default.get('tolerance'))
        for key, value in (('dpi_href', 479), ('dpi_slope', 0.05), ('dpi_max', 0.3)):
            self.assertEqual(default.get(key) if default.get(key) is not None else value, value)

    def test_all_parameters_pass_through_for_local_and_lsf(self):
        for platform in ('local', 'lsf'):
            with self.subTest(platform=platform):
                job, command = self.run_wrapper('--dynamic-dpi', '--dpi-href', '2',
                                                  '--dpi-slope', '0.1', '--dpi-max', '0.25',
                                                  platform=platform)
                for key, value in (('dynamic_dpi', True), ('dpi_href', 2), ('dpi_slope', 0.1), ('dpi_max', 0.25)):
                    self.assertEqual(job[key], value)
                self.assertIn('cwlexec' if platform == 'lsf' else 'cwltool', command)
                self.assertIsNone(job.get('tolerance'))

    def test_zero_cap_is_valid(self):
        job, _ = self.run_wrapper('--dynamic-dpi', '--dpi-max', '0')
        self.assertEqual(job['dpi_max'], 0)

    def test_invalid_parameters_fail_before_launch_or_output(self):
        cases = [(flag, value) for flag, values in (
            ('--dpi-href', ('0', '-1', '1.5', 'nan', '2147483648')),
            ('--dpi-slope', ('0', '-1', 'nan', 'inf', 'invalid')),
            ('--dpi-max', ('-0.1', '1', 'nan', 'inf', 'invalid')),
        ) for value in values]
        for flag, value in cases:
            with self.subTest(flag=flag, value=value):
                self.assert_rejected('--dynamic-dpi', flag, value)

    def test_tuning_is_rejected_when_dynamic_is_explicitly_disabled(self):
        for flag, value in (('--dpi-href', '479'), ('--dpi-slope', '0.05'), ('--dpi-max', '0.3')):
            with self.subTest(flag=flag):
                self.assert_rejected('--fixed-dpi', flag, value)

    def assert_rejected(self, *extra):
        module = load_wrapper()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output, scratch = root / 'output', root / 'scratch'
            expression, hubs = write_multihub_dpi_inputs(root)
            args = ['sjaracne', 'local', '--engine', 'cwl', '--network-type', 'TF',
                    '-e', str(expression), '-g', str(hubs), '-o', str(output), '-tmp', str(scratch), *extra]
            with mock.patch.object(module, 'run_shell_command_call') as launch, \
                    mock.patch.object(module, 'runtime_fingerprint', return_value={}), \
                    mock.patch.object(module, 'validate_cwl_outputs', return_value={}), \
                    mock.patch.object(module, 'require_tool', side_effect=lambda name: name):
                with mock.patch.object(sys, 'argv', args), self.assertRaises(SystemExit):
                    module.main()
                launch.assert_not_called()
            self.assertFalse(output.exists())
            self.assertFalse(scratch.exists())

    def test_cwl_mode_and_parameters_are_forwarded_to_every_resample(self):
        workflow = (ROOT / 'SJARACNe/cwl/sjaracne_workflow.cwl').read_text(encoding='utf-8')
        command = (ROOT / 'SJARACNe/cwl/sjaracne.cwl').read_text(encoding='utf-8')
        for name in ('dynamic_dpi', 'dpi_href', 'dpi_slope', 'dpi_max', 'tolerance'):
            self.assertIn(name + ': ' + name, workflow)
        for flag in ('--dynamic-dpi', '--dpi-href', '--dpi-slope', '--dpi-max'):
            self.assertIn('prefix: ' + flag, command)
        self.assertIn('tolerance:\n    type: double?', command)
        self.assertIn("hasValue(fixed) ? fixed : 0", command)
        self.assertIn('dynamic_dpi conflicts with an explicit fixed tolerance', command)


@unittest.skipUnless(NATIVE, 'sjaracne.exe is not built')
class TestDynamicDpiRealCwl(unittest.TestCase):
    @staticmethod
    def header(text, key):
        match = re.search(r'^>\s+' + re.escape(key) + r'\s+([^\r\n]+)$', text, re.M)
        if match is None:
            raise AssertionError('Missing header ' + key + '\n' + text)
        return match.group(1).strip()

    @staticmethod
    def edges(path):
        return '\n'.join(line for line in path.read_text(encoding='utf-8').splitlines()
                         if line and not line.startswith('>'))

    def run_command_cwl(self, workdir, expression, hubs, **options):
        job = {'exp_file': {'class': 'File', 'path': str(expression)},
               'probe_file_tf': {'class': 'File', 'path': str(hubs)},
               'probe_file_subnetwork': {'class': 'File', 'path': str(hubs)},
               'aracne_config_dir': {'class': 'Directory', 'path': str(ROOT / 'SJARACNe/config')},
               'p_value': 1, 'npar_limit': 40, 'seed': 19,
               'subsample_spec': '80%', 'output_file_name': 'network.adj', **options}
        workdir.mkdir()
        inputs = workdir / 'job.json'
        inputs.write_text(json.dumps(job), encoding='ascii')
        result = subprocess.run(['cwltool', '--outdir', str(workdir / 'out'),
                                 '--tmpdir-prefix', str(workdir / 'scratch'),
                                 str(ROOT / 'SJARACNe/cwl/sjaracne.cwl'), str(inputs)],
                                cwd=str(ROOT), capture_output=True, text=True)
        return result, workdir / 'out/network.adj'

    def test_real_command_cwl_matches_exact_effective_fixed_epsilon(self):
        with tempfile.TemporaryDirectory() as temporary:
            workdir = Path(temporary)
            expression, hubs = write_multihub_dpi_inputs(workdir)
            dynamic, path = self.run_command_cwl(workdir / 'dynamic', expression, hubs,
                                                dynamic_dpi=True, dpi_href=2,
                                                dpi_slope=0.1, dpi_max=0.3)
            self.assertEqual(dynamic.returncode, 0, dynamic.stdout + dynamic.stderr)
            text = path.read_text(encoding='utf-8')
            self.assertEqual(self.header(text, 'DPI mode'), 'experimental-dynamic')
            self.assertEqual(self.header(text, 'DPI source rows H'), '3')
            epsilon = float.fromhex(self.header(text, 'DPI epsilon hex'))
            self.assertGreater(epsilon, 0)
            fixed, reference = self.run_command_cwl(workdir / 'fixed', expression, hubs,
                                                    tolerance=epsilon)
            self.assertEqual(fixed.returncode, 0, fixed.stdout + fixed.stderr)
            self.assertEqual(self.edges(path), self.edges(reference))

    def test_direct_cwl_rejects_fixed_dynamic_conflict_or_ignored_tuning(self):
        with tempfile.TemporaryDirectory() as temporary:
            workdir = Path(temporary)
            expression, hubs = write_multihub_dpi_inputs(workdir)
            for name, options, message in (
                ('conflict', {'dynamic_dpi': True, 'tolerance': 0}, 'conflicts with an explicit fixed tolerance'),
                ('no_opt_in', {'dpi_slope': 0.1}, 'require dynamic_dpi'),
            ):
                with self.subTest(name=name):
                    result, path = self.run_command_cwl(workdir / name, expression, hubs, **options)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn(message, result.stdout + result.stderr)
                    self.assertFalse(path.exists())

    def test_wrapper_executes_two_real_dynamic_resamples_and_consensus(self):
        with tempfile.TemporaryDirectory() as temporary:
            workdir = Path(temporary)
            expression, hubs = write_multihub_dpi_inputs(workdir)
            output = workdir / 'output'
            result = subprocess.run([
                sys.executable, '-m', 'SJARACNe.sjaracne', 'local',
                '--engine', 'cwl', '--network-type', 'TF',
                '-e', str(expression), '-g', str(hubs), '-n', '2', '--min-recurrence', '1',
                '-pb', '1', '--dynamic-dpi', '--dpi-href', '2',
                '--dpi-slope', '0.1', '--dpi-max', '0.3',
                '-o', str(output), '-tmp', str(workdir / 'scratch'),
            ], cwd=str(ROOT), capture_output=True, text=True)
            log = result.stdout + result.stderr
            self.assertEqual(result.returncode, 0, log)
            self.assertEqual(log.count('[DPI] mode=experimental-dynamic'), 2, log)
            self.assertNotIn('Skipping DPI', log)
            self.assertTrue((output / 'consensus_network_ncol_.txt').is_file())
            parameters = (output / 'parameter_info_.txt').read_text(encoding='utf-8')
            self.assertEqual(self.header(parameters, 'DPI mode'), 'experimental-dynamic')
            self.assertEqual(self.header(parameters, 'DPI H_ref'), '2')
            self.assertEqual(self.header(parameters, 'DPI source rows H'), '3')


if __name__ == '__main__':
    unittest.main()
