#!/usr/bin/env python3

import unittest
import tempfile
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class TestSJARACNe(unittest.TestCase):
    def test_pchip_workflow_on_tiny_synthetic_input(self):
        for pvalue, extra in [('0.001', []), ('1e-7', ['--allow-rare-tail'])]:
            with self.subTest(p=pvalue), tempfile.TemporaryDirectory() as tmpdir:
                workdir = Path(tmpdir)
                expression = workdir / 'input.exp'
                lines = ['isoformId\tgeneSymbol\t' + '\t'.join('s%d' % i for i in range(100))]
                for name, values in [('H', range(100)), ('A', range(99, -1, -1)),
                                     ('B', [(i * 31) % 100 for i in range(100)])]:
                    lines.append('\t'.join([name, name] + [str(value) for value in values]))
                expression.write_text('\n'.join(lines) + '\n', encoding='ascii')
                hubs = workdir / 'hubs.txt'
                hubs.write_text('H\n', encoding='ascii')
                output = workdir / 'output'
                result = subprocess.run([
                    sys.executable, '-m', 'SJARACNe.sjaracne', 'local',
                    '-e', str(expression), '-g', str(hubs), '-n', '2', '-pc', '1',
                    '-pb', pvalue, '-M', 'pchip', '--subsample-size', '80',
                    '-o', str(output), '-tmp', str(workdir / 'tmp'), *extra,
                ], cwd=str(PROJECT_ROOT), capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                log = result.stdout + result.stderr
                self.assertEqual(log.count('[THRESHOLD] PCHIP cutoff binary64='), 2)
                self.assertTrue((output / 'consensus_network_ncol_.txt').is_file())
                if extra:
                    self.assertIn('rare-tail use explicitly acknowledged', log)

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

