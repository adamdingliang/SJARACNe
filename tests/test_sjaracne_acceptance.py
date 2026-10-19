#!/usr/bin/env python3

import unittest
import tempfile
import subprocess
import sys
import json
import csv
import math
import random
import re
import warnings
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


def write_recurrence_inputs(workdir):
    """Small deterministic fixture with stable and sample-sensitive edges."""
    expression, hubs = write_tiny_inputs(workdir)
    noise = list(range(100))
    random.Random(73).shuffle(noise)
    with expression.open('a', encoding='ascii') as handle:
        for weight in (42, 49):
            values = [weight * i + (100 - weight) * noise[i] for i in range(100)]
            name = 'W%d' % weight
            handle.write('\t'.join([name, name] + [str(value) for value in values]) + '\n')
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
            self.assertTrue((output_dir / 'bootstrap_info_.txt').is_file())
            self.assertTrue((output_dir / 'parameter_info_.txt').is_file())

    @staticmethod
    def header_value(text, key):
        match = re.search(r'^>\s+' + re.escape(key) + r'\s+([^\r\n]+)$', text, re.M)
        if match is None:
            raise AssertionError('Missing native header %r:\n%s' % (key, text))
        return match.group(1).strip()

    def read_retained_bootstraps(self, workdir):
        # The unmodified workflow does not publish its bootstrap files. Keep its
        # real temporary outputs, including copies made by copy_to_dir, and
        # count each named resample exactly once rather than counting copies.
        contents = {}
        for path in sorted(workdir.rglob('TF_run_*.adj')):
            text = path.read_text(encoding='utf-8')
            if path.name in contents:
                self.assertEqual(contents[path.name], text, str(path))
            contents[path.name] = text
        self.assertEqual(set(contents), {'TF_run_%03d.adj' % i for i in range(1, 7)})
        networks = []
        for name in sorted(contents):
            edges = {}
            for line in contents[name].splitlines():
                if not line or line.startswith('>'):
                    continue
                fields = line.split('\t')
                self.assertEqual(len(fields[1:]) % 2, 0, line)
                for index in range(1, len(fields), 2):
                    edge = (fields[0], fields[index])
                    self.assertNotIn(edge, edges, name)
                    edges[edge] = float(fields[index + 1])
                    self.assertTrue(math.isfinite(edges[edge]))
            networks.append(edges)
        return contents, networks

    def assert_consensus_rows(self, path, expected):
        with path.open(encoding='utf-8', newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        keys = [(row['source'], row['target']) for row in rows]
        self.assertEqual(keys, sorted(expected), str(path))
        self.assertEqual({key: row['MI'] for key, row in zip(keys, rows)}, expected)

    def test_consensus_selection_and_null_model_real_cwl_matrix(self):
        from SJARACNe import get_mi_cutoff

        native_networks = {}
        retained_by_mode = {}
        probability = 0.001
        config = PROJECT_ROOT / 'SJARACNe' / 'config'
        coefficient_rows = [line for line in (config / 'config_threshold.txt').read_text().splitlines()
                            if line and not line.startswith('>')]
        alpha, beta, gamma = map(float, coefficient_rows[0].split())
        affine_cutoff = (alpha - math.log(probability)) / (-beta - gamma * 80)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pchip_cutoff = get_mi_cutoff(80, probability)

        for model in (None, 'pchip'):
            for selection in ('default-k6', 'explicit-pc'):
                with self.subTest(model=model or 'affine', selection=selection), \
                        tempfile.TemporaryDirectory() as tmpdir:
                    workdir = Path(tmpdir)
                    expression, hubs = write_recurrence_inputs(workdir)
                    output = workdir / 'output'
                    job = {
                        'exp_file': {'class': 'File', 'path': str(expression)},
                        'probe_file': {'class': 'File', 'path': str(hubs)},
                        'aracne_config_dir': {'class': 'Directory', 'path': str(config)},
                        'bootstrap_num': 6, 'p_value_bootstrap': probability,
                        'subsample_spec': '80', 'final_out_dir_name': 'results',
                    }
                    if model:
                        job['apmi_null_model'] = model
                    if selection == 'explicit-pc':
                        job['p_value_consensus'] = 1
                    job_file = workdir / 'job.json'
                    job_file.write_text(json.dumps(job), encoding='ascii')
                    result = subprocess.run([
                        'cwltool', '--leave-tmpdir', '--outdir', str(output),
                        '--tmpdir-prefix', str(workdir / 'tmp'),
                        str(PROJECT_ROOT / 'SJARACNe/cwl/sjaracne_workflow.cwl'), str(job_file),
                    ], cwd=str(PROJECT_ROOT), capture_output=True, text=True)
                    log = result.stdout + result.stderr
                    self.assertEqual(result.returncode, 0, log)
                    headers, networks = self.read_retained_bootstraps(workdir)
                    support = {}
                    total_mi = {}
                    for network in networks:
                        for edge, mi in network.items():
                            support[edge] = support.get(edge, 0) + 1
                            total_mi[edge] = total_mi.get(edge, 0.0) + mi
                    self.assertEqual(support.get(('H', 'A')), 6)
                    self.assertTrue(any(0 < count < 6 for count in support.values()), support)
                    if model in native_networks:
                        self.assertEqual(native_networks[model], networks)
                    native_networks[model] = networks

                    # Derive the legacy moments from actual emitted edges. In
                    # this small fixture every finite normal tail is below pc=1,
                    # so explicit legacy mode retains the union, independently
                    # of the production consensus function and its uprob helper.
                    mu = sum(float(len(network)) / len(support) for network in networks)
                    variance = sum((float(len(network)) / len(support)) *
                                   (1.0 - float(len(network)) / len(support))
                                   for network in networks)
                    if selection == 'default-k6':
                        selected = {edge for edge, count in support.items() if count >= 6}
                    else:
                        if variance:
                            self.assertTrue(all(
                                0.5 * math.erfc((count - mu) / math.sqrt(2 * variance)) < 1
                                for count in support.values()
                            ))
                        selected = set(support)
                    expected = {edge: '{:.4f}'.format(total_mi[edge] / support[edge])
                                for edge in selected}
                    self.assertTrue(expected)
                    retained_by_mode[(model, selection)] = selected
                    self.assert_consensus_rows(output / 'consensus_network_ncol_.txt', expected)
                    three_column = list(workdir.rglob('consensus_network_3col_.txt'))
                    self.assertEqual(len(three_column), 1)
                    self.assert_consensus_rows(three_column[0], expected)

                    bootstrap_info = (output / 'bootstrap_info_.txt').read_text(encoding='utf-8')
                    metadata = dict(line.split(': ', 1) for line in bootstrap_info.splitlines())
                    self.assertEqual(metadata['Total edge tested'], str(len(support)))
                    parameter_info = (output / 'parameter_info_.txt').read_text(encoding='utf-8')
                    self.assertIn('>  Bootstrap No: 6\n', parameter_info)
                    if selection == 'default-k6':
                        self.assertEqual(metadata, {
                            'Total edge tested': str(len(support)),
                            'Consensus selection': 'minimum recurrence',
                            'Bootstrap networks': '6', 'Minimum recurrence': '6',
                            'Minimum recurrence fraction': '1',
                        })
                        self.assertIn('>  Consensus selection minimum recurrence\n', parameter_info)
                        self.assertIn('>  Minimum recurrence 6 of 6\n', parameter_info)
                        self.assertIn('>  Minimum recurrence fraction 1\n', parameter_info)
                        self.assertNotIn('normal approximation', log)
                    else:
                        self.assertAlmostEqual(float(metadata['mu']), mu, places=14)
                        self.assertAlmostEqual(float(metadata['sigma']), math.sqrt(variance), places=14)
                        self.assertNotIn('Minimum recurrence', bootstrap_info + parameter_info)
                        self.assertIn('Legacy consensus p-value filtering uses a normal approximation', log)

                    for header in headers.values():
                        self.assertEqual(self.header_value(header, 'Eligible observations'), '100')
                        self.assertEqual(self.header_value(header, 'Sampled observations'), '80')
                        self.assertEqual(self.header_value(header, 'Sampling request'), '80')
                        self.assertEqual(self.header_value(header, 'Sampling method'),
                                         'fixed-size without replacement')
                        self.assertEqual(float(self.header_value(header, 'MI P-value')), probability)
                        if model:
                            self.assertEqual(self.header_value(header, 'MI threshold method'),
                                             'experimental PCHIP empirical cutoff')
                            self.assertEqual(self.header_value(header, 'AP-MI actual m'), '80')
                            self.assertEqual(self.header_value(header, 'AP-MI Npar'), '40')
                            self.assertEqual(float(self.header_value(header, 'MI threshold binary64')).hex(),
                                             pchip_cutoff.hex())
                            self.assertEqual(float.fromhex(self.header_value(header, 'MI threshold hex')),
                                             pchip_cutoff)
                        else:
                            self.assertNotIn('AP-MI PCHIP', header)
                            self.assertNotIn('MI threshold binary64', header)
                            self.assertAlmostEqual(float(self.header_value(header, 'MI threshold')),
                                                   affine_cutoff, delta=1e-6)
                    if model:
                        self.assertEqual(log.count('[THRESHOLD] PCHIP cutoff binary64='), 6)
                        self.assertNotIn('legacy affine threshold', log)
                    else:
                        self.assertEqual(log.count('using the legacy affine threshold'), 6)
                        self.assertNotIn('PCHIP cutoff binary64=', log)
                    if selection == 'explicit-pc':
                        self.assertNotEqual(retained_by_mode[(model, 'default-k6')], selected)
                    print('\nRECURRENCE_MATRIX_PASS=' + json.dumps({
                        'model': model or 'affine',
                        'consensus': 'default_k6' if selection == 'default-k6' else 'legacy_pc',
                        'B': len(networks), 'K': 6 if selection == 'default-k6' else None,
                        'K_over_B': 1 if selection == 'default-k6' else None,
                        'p_value_consensus': 1 if selection == 'explicit-pc' else None,
                        'actual_m': 80, 'cutoff': pchip_cutoff if model else affine_cutoff,
                        'observed_edges': len(support), 'selected_edges': len(selected),
                        'support_counts': {'%s->%s' % edge: count for edge, count in support.items()},
                    }, sort_keys=True), flush=True)


if __name__ == '__main__':
    unittest.main()

