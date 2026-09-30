#!/usr/bin/env python3

import unittest
import tempfile
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class TestSJARACNe(unittest.TestCase):
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

    def test_default_minimum_recurrence_acceptance(self):
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
                '6',
                '-o',
                str(output_dir),
                '-tmp',
                str(temp_dir),
            ], cwd=str(PROJECT_ROOT))

            self.assertTrue(
                (output_dir / 'consensus_network_ncol_.txt').is_file()
            )
            bootstrap_info = (output_dir / 'bootstrap_info_.txt').read_text(
                encoding='utf-8'
            )
            self.assertIn('Consensus selection: minimum recurrence\n', bootstrap_info)
            self.assertIn('Bootstrap networks: 6\n', bootstrap_info)
            self.assertIn('Minimum recurrence: 6\n', bootstrap_info)
            self.assertIn('Minimum recurrence fraction: 1\n', bootstrap_info)

            parameter_info = (output_dir / 'parameter_info_.txt').read_text(
                encoding='utf-8'
            )
            self.assertIn('>  Bootstrap No: 6\n', parameter_info)
            self.assertIn('>  Consensus selection minimum recurrence\n', parameter_info)
            self.assertIn('>  Minimum recurrence 6 of 6\n', parameter_info)


if __name__ == '__main__':
    unittest.main()

