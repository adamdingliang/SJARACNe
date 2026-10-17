#!/usr/bin/env python3

import unittest
import tempfile
import subprocess
import sys
import json
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def write_tiny_inputs(workdir):
    expression = workdir / 'input.exp'
    lines = ['isoformId\tgeneSymbol\t' + '\t'.join('s%d' % i for i in range(100))]
    for name, values in [('H', range(100)), ('A', range(99, -1, -1)),
                         ('B', [(i * 31) % 100 for i in range(100)])]:
        lines.append('\t'.join([name, name] + [str(value) for value in values]))
    expression.write_text('\n'.join(lines) + '\n', encoding='ascii')
    hubs = workdir / 'hubs.txt'
    hubs.write_text('H\n', encoding='ascii')
    return expression, hubs


class TestSJARACNe(unittest.TestCase):
    def test_pchip_workflow_on_tiny_synthetic_input(self):
        for pvalue in ('0.001', '0.0004', '1e-7'):
            with self.subTest(p=pvalue), tempfile.TemporaryDirectory() as tmpdir:
                workdir = Path(tmpdir)
                expression, hubs = write_tiny_inputs(workdir)
                output = workdir / 'output'
                result = subprocess.run([
                    sys.executable, '-m', 'SJARACNe.sjaracne', 'local',
                    '-e', str(expression), '-g', str(hubs), '-n', '2', '-pc', '1',
                    '-pb', pvalue, '-M', 'pchip', '--subsample-size', '80',
                    '-o', str(output), '-tmp', str(workdir / 'tmp'),
                ], cwd=str(PROJECT_ROOT), capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                log = result.stdout + result.stderr
                self.assertEqual(log.count('[THRESHOLD] PCHIP cutoff binary64='), 2)
                self.assertTrue((output / 'consensus_network_ncol_.txt').is_file())
                self.assertNotIn('acknowledged', log)
                if float(pvalue) < 0.0005:
                    self.assertIn('not a demonstrated failure threshold', log)
                if float(pvalue) <= 1e-5:
                    self.assertIn('extreme-tail Monte Carlo precision', log)
                else:
                    self.assertNotIn('extreme-tail Monte Carlo precision', log)

    def test_direct_cwl_probability_omission_policy(self):
        for model in (None, 'pchip'):
            with self.subTest(model=model), tempfile.TemporaryDirectory() as tmpdir:
                workdir = Path(tmpdir)
                expression, hubs = write_tiny_inputs(workdir)
                job = {
                    'exp_file': {'class': 'File', 'path': str(expression)},
                    'probe_file': {'class': 'File', 'path': str(hubs)},
                    'aracne_config_dir': {'class': 'Directory', 'path': str(PROJECT_ROOT / 'SJARACNe/config')},
                    'p_value_consensus': 1, 'bootstrap_num': 1,
                    'subsample_spec': '80', 'final_out_dir_name': 'results',
                }
                if model:
                    job['apmi_null_model'] = model
                job_file = workdir / 'job.json'
                job_file.write_text(json.dumps(job), encoding='ascii')
                output = workdir / 'output'
                result = subprocess.run([
                    'cwltool', '--outdir', str(output), '--tmpdir-prefix', str(workdir / 'tmp'),
                    str(PROJECT_ROOT / 'SJARACNe/cwl/sjaracne_workflow.cwl'), str(job_file),
                ], cwd=str(PROJECT_ROOT), capture_output=True, text=True)
                log = result.stdout + result.stderr
                if model:
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn('PCHIP requires explicit p_value_bootstrap', log)
                    self.assertNotIn('PCHIP cutoff binary64=', log)
                else:
                    self.assertEqual(result.returncode, 0, log)
                    self.assertIn('legacy affine threshold', log)
                    self.assertTrue((output / 'consensus_network_ncol_.txt').is_file())

    def test_acceptance(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            workdir = Path(tmpdir)
            output_dir = workdir / 'output'
            temp_dir = workdir / 'tmp'
            subprocess.check_call([
                sys.executable,
                '-m',
                'SJARACNe.sjaracne',
                'local',
                '-e',
                str(PROJECT_ROOT / 'tests' / 'inputs' / 'Tcell1170.exp'),
                '-g',
                str(PROJECT_ROOT / 'tests' / 'inputs' / 'TcellTF.txt'),
                '-n',
                '5',
                '-pc',
                '0.01',
                '-o',
                str(output_dir),
                '-tmp',
                str(temp_dir),
            ], cwd=str(PROJECT_ROOT))

            expected = PROJECT_ROOT / 'tests' / 'answerkey' / 'acceptance' / 'cnn_5.txt'
            actual = output_dir / 'consensus_network_ncol_.txt'
            self.assertEqual(
                expected.read_text(encoding='utf-8'),
                actual.read_text(encoding='utf-8'),
            )


if __name__ == '__main__':
    unittest.main()

