"""Behavioral checks for native pre-DPI recurrence consensus and raw MI output."""
import math
import os
from pathlib import Path
import random
import signal
import subprocess
import tempfile
import unittest

from native_test_support import find_native_executable


NATIVE = find_native_executable()


@unittest.skipUnless(NATIVE, 'sjaracne.exe is not built')
class TestNativeConsensus(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.work = Path(self.temp.name)
        self.exp = self.work / 'input.exp'
        self.hubs = self.work / 'hubs.txt'
        self.annotation = self.work / 'annotation.txt'
        self.manifest = self.work / 'runs.txt'
        self.output = self.work / 'result.adj'
        shuffled = list(range(1, 101))
        random.Random(4).shuffle(shuffled)
        values = [list(range(1, 101)), list(range(100, 0, -1)), shuffled,
                  [i * i % 97 for i in range(100)]]
        lines = ['isoformId\tgeneSymbol\t' + '\t'.join('s%d' % i for i in range(100))]
        lines += [gene + '\t' + gene + '\t' + '\t'.join(map(str, data))
                  for gene, data in zip('ABCD', values)]
        self.exp.write_text('\n'.join(lines) + '\n')
        self.hubs.write_text('A\nB\nC\nD\n')
        self.annotation.write_text('')
        self.runs = []

    def tearDown(self):
        self.temp.cleanup()

    def write_runs(self, bodies, header='>  DPI tolerance   1\n'):
        self.runs = []
        for i, body in enumerate(bodies):
            path = self.work / ('raw run %03d.adj' % i)
            path.write_text(header + body)
            self.runs.append(path)
        self.manifest.write_text(''.join(path.name + '\n' for path in self.runs))

    def command(self, k=1, epsilon='1', extra=None, output=None):
        command = [NATIVE, '-i', str(self.exp), '-s', str(self.hubs), '-l', str(self.annotation),
                   '--consensus-inputs', str(self.manifest), '--min-recurrence', str(k),
                   '-p', '1', '-t', '0', '-o', str(output or self.output)]
        if epsilon is not None:
            command += ['-e', str(epsilon)]
        return command + (extra or [])

    def invoke(self, command, **kwargs):
        return subprocess.run(command, cwd=self.work, capture_output=True, text=True, **kwargs)

    @staticmethod
    def rows(path):
        return [line for line in path.read_text().splitlines() if line and not line.startswith('>')]

    @classmethod
    def edges(cls, path):
        result = {}
        for line in cls.rows(path):
            fields = line.split('\t')
            for i in range(1, len(fields), 2):
                result[(fields[0], fields[i])] = float(fields[i + 1])
        return result

    @staticmethod
    def header(path, name):
        prefix = '>  ' + name + ' '
        found = [line[len(prefix):] for line in path.read_text().splitlines() if line.startswith(prefix)]
        if len(found) != 1:
            raise AssertionError('Expected one %r header, got %r' % (name, found))
        return found[0].strip()

    def assert_success(self, command):
        result = self.invoke(command)
        self.assertEqual(result.returncode, 0, result.stderr + '\n' + result.stdout)
        self.assertFalse(list(self.work.glob('result.adj.tmp.*')))
        return result

    def test_inclusive_k_and_supporting_run_mean_preserve_direction_and_zero(self):
        self.write_runs([
            'A\tB\t0.12345678901234566\tC\t0.9\tD\t0\nB\tA\t0.7\nC\nD\n',
            'A\tB\t0.32345678901234567\tD\t0\nB\nC\nD\n',
            'A\nB\tA\t0.6\nC\nD\n',
        ])
        self.assert_success(self.command(k=2))
        edges = self.edges(self.output)
        self.assertEqual(set(edges), {('A', 'B'), ('A', 'D'), ('B', 'A')})
        self.assertEqual(edges[('A', 'B')], (0.12345678901234566 + 0.32345678901234567) / 2)
        self.assertEqual(edges[('A', 'D')], 0.0)
        self.assertEqual(edges[('B', 'A')], (0.7 + 0.6) / 2)
        self.assertEqual(self.header(self.output, 'Consensus union edges'), '4')
        self.assertEqual(self.header(self.output, 'Consensus candidate edges'), '3')
        self.assertEqual(self.header(self.output, 'Consensus final edges'), '3')
        self.assertEqual(self.header(self.output, 'Consensus input runs B'), '3')
        self.assertEqual(self.header(self.output, 'Consensus minimum recurrence K'), '2')
        self.assertEqual([row.split('\t')[0] for row in self.rows(self.output)], list('ABCD'))

    def test_subset_sources_keep_requested_order_and_independent_support(self):
        # Source ordinals differ from their expression-gene indices. Shared
        # targets must contribute separately to each source, including zeros.
        self.hubs.write_text('C\nA\n')
        self.write_runs([
            'C\tB\t0.25\tD\t0\nA\tB\t0.5\tD\t0\n',
            'C\tB\t0.75\tD\t0\nA\tB\t1.5\tD\t0\n',
            'C\tB\t1.25\nA\tD\t0\n',
        ])
        self.assert_success(self.command(k=2))
        self.assertEqual(self.edges(self.output), {
            ('C', 'B'): 0.75, ('C', 'D'): 0.0,
            ('A', 'B'): 1.0, ('A', 'D'): 0.0,
        })
        self.assertEqual([row.split('\t')[0] for row in self.rows(self.output)], ['C', 'A'])
        self.assertEqual(self.header(self.output, 'Consensus source rows H'), '2')
        self.assertEqual(self.header(self.output, 'Consensus union edges'), '4')
        self.assertEqual(self.header(self.output, 'Consensus candidate edges'), '4')

    def test_finite_mi_sum_overflow_preserves_existing_output(self):
        self.write_runs(['A\tB\t1e308\nB\nC\nD\n'] * 2)
        self.output.write_text('previous good network\n')
        result = self.invoke(self.command(k=2))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('MI sum overflow', result.stderr)
        self.assertEqual(self.output.read_text(), 'previous good network\n')
        self.assertFalse(list(self.work.glob('result.adj.tmp.*')))

    def test_dpi_matches_existing_replay_on_the_same_consensus_triangle(self):
        body = 'A\tB\t0.1\tC\t0.9\nB\tC\t0.8\nC\nD\n'
        self.write_runs([body, body])
        self.assert_success(self.command(k=2, epsilon='0'))
        native_edges = self.edges(self.output)
        self.assertEqual(native_edges, {('A', 'C'): 0.9, ('B', 'C'): 0.8})
        replay = self.work / 'replay.adj'
        result = self.invoke([NATIVE, '-i', str(self.exp), '-s', str(self.hubs), '-l', str(self.annotation),
                              '-j', str(self.runs[0]), '-p', '1', '-t', '0', '-e', '0', '-o', str(replay)])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.edges(replay), native_edges)

    def test_consensus_is_deterministic_in_explicit_manifest_order(self):
        values = [1.0, 2.0 ** -53, 2.0 ** -53]
        self.write_runs(['A\tB\t%.17g\nB\nC\nD\n' % value for value in values])
        self.assert_success(self.command())
        first = self.output.read_bytes()
        total = 0.0
        for value in values:
            total += value
        self.assertEqual(self.edges(self.output)[('A', 'B')], total / 3)
        self.assert_success(self.command(extra=['-S', '99']))
        self.assertEqual(self.output.read_bytes(), first)
        self.manifest.write_text(''.join(path.name + '\n' for path in reversed(self.runs)))
        self.assert_success(self.command())
        total = 0.0
        for value in reversed(values):
            total += value
        self.assertEqual(self.edges(self.output)[('A', 'B')], total / 3)

    def test_support_counter_handles_more_than_255_distinct_runs(self):
        self.write_runs(['A\tB\t0\nB\nC\nD\n'] * 300)
        self.assert_success(self.command(k=300))
        self.assertEqual(self.edges(self.output), {('A', 'B'): 0.0})
        self.assertEqual(self.header(self.output, 'Consensus input runs B'), '300')

    @unittest.skipUnless(os.name == 'posix', 'uses POSIX open-file limits')
    def test_many_resamples_do_not_require_many_open_files(self):
        import resource
        self.write_runs(['A\tB\t0.25\nB\nC\nD\n'] * 80)

        def limit_open_files():
            soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
            limited = 32 if soft == resource.RLIM_INFINITY else min(32, soft)
            resource.setrlimit(resource.RLIMIT_NOFILE, (limited, hard))

        result = self.invoke(self.command(k=80), preexec_fn=limit_open_files)
        self.assertEqual(result.returncode, 0, result.stderr + '\n' + result.stdout)
        self.assertEqual(self.edges(self.output), {('A', 'B'): 0.25})
        self.assertEqual(self.header(self.output, 'Consensus input runs B'), '80')

    @unittest.skipUnless(os.name == 'posix', 'uses POSIX address-space limits')
    def test_accumulator_allocation_failure_preserves_existing_output(self):
        import resource
        if not hasattr(resource, 'RLIMIT_AS'):
            self.skipTest('address-space limits are unavailable')
        address_limit = 128 * 1024 * 1024
        try:
            _, hard = resource.getrlimit(resource.RLIMIT_AS)
        except (OSError, ValueError):
            self.skipTest('address-space limits are unsupported')
        if hard != resource.RLIM_INFINITY and hard < address_limit:
            self.skipTest('existing hard limit is too low for the loading fixture')

        # Loading 5000 genes and 100 observations takes much less than 128 MiB,
        # while the full 12*H*G consensus accumulator requires 300,000,000 bytes.
        # Limit only the native child, leaving the test runner unaffected.
        genes = ['GENE_%04d' % i for i in range(5000)]
        observations = '\t'.join(str(i + 1) for i in range(100))
        header = 'isoformId\tgeneSymbol\t' + '\t'.join('s%d' % i for i in range(100))
        self.exp.write_text(header + '\n' + ''.join(
            gene + '\t' + gene + '\t' + observations + '\n' for gene in genes))
        self.hubs.write_text('\n'.join(genes) + '\n')
        self.write_runs(['\n'.join(genes) + '\n'])
        self.output.write_text('previous good network\n')

        def limit_address_space():
            resource.setrlimit(resource.RLIMIT_AS, (address_limit, hard))

        result = self.invoke(self.command(), preexec_fn=limit_address_space)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Cannot allocate sequential consensus accumulators', result.stderr)
        self.assertIn('H=5000', result.stderr)
        self.assertIn('G=5000', result.stderr)
        self.assertIn('300000000 bytes', result.stderr)
        self.assertEqual(self.output.read_text(), 'previous good network\n')
        self.assertFalse(list(self.work.glob('result.adj.tmp.*')))

    def test_each_raw_file_is_fully_validated_before_reading_the_next(self):
        # The first graph fails at its last source; the second fails at its
        # first. File-major reading must report the first graph's error.
        self.write_runs([
            'A\nB\nC\nD\tA\tnan\n',
            'A\tB\tinf\nB\nC\nD\n',
        ])
        result = self.invoke(self.command())
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(self.runs[0].name, result.stderr)
        self.assertNotIn(self.runs[1].name, result.stderr)
        self.assertIn('nan', result.stderr)
        self.assertFalse(self.output.exists())

    def test_all_empty_consensus_preserves_sources_and_dynamic_h(self):
        self.write_runs(['A\nB\nC\nD\n'] * 2)
        self.assert_success(self.command(k=2, epsilon=None, extra=[
            '--dynamic-dpi', '--dpi-href', '2', '--dpi-offset', '0.025']))
        self.assertEqual(self.rows(self.output), list('ABCD'))
        self.assertEqual(self.header(self.output, 'Consensus source rows H'), '4')
        self.assertEqual(self.header(self.output, 'DPI source rows H'), '4')
        self.assertEqual(float(self.header(self.output, 'DPI epsilon binary64')), 0.025 + 0.05)

    def test_dynamic_offset_zero_baseline_cap_and_parameter_bounds(self):
        self.write_runs(['A\tB\t0.2\nB\nC\nD\n'])
        cases = [('8', '0', '0.3', 0.0), ('8', '0.025', '0.3', 0.025),
                 ('2', '0.025', '0.03', 0.03)]
        for href, offset, cap, expected in cases:
            with self.subTest(href=href, offset=offset, cap=cap):
                self.assert_success(self.command(epsilon=None, extra=[
                    '--dynamic-dpi', '--dpi-href', href, '--dpi-offset', offset, '--dpi-max', cap]))
                self.assertEqual(float(self.header(self.output, 'DPI epsilon binary64')), expected)
        for bad in ('nan', 'inf', '-0.01', '1', '0.1junk'):
            with self.subTest(bad_offset=bad):
                self.assertNotEqual(self.invoke(self.command(epsilon=None, extra=[
                    '--dynamic-dpi', '--dpi-offset', bad])).returncode, 0)
        self.assertNotEqual(self.invoke(self.command(extra=['--dpi-offset', '0'])).returncode, 0)

    def test_default_k20_and_invalid_k_fail_before_publication(self):
        self.write_runs(['A\nB\nC\nD\n'])
        base = self.command()
        start = base.index('--min-recurrence')
        del base[start:start + 2]
        self.assertNotEqual(self.invoke(base).returncode, 0)
        self.assertFalse(self.output.exists())
        for k in ('0', '-1', '1.5', '2', '2147483648', 'nan'):
            with self.subTest(k=k):
                self.assertNotEqual(self.invoke(self.command(k=k)).returncode, 0)
                self.assertFalse(self.output.exists())

    def test_missing_reordered_duplicate_and_extra_source_rows_fail(self):
        cases = ['A\nB\nC\n', 'B\nA\nC\nD\n', 'A\nA\nC\nD\n',
                 'A\nB\nC\nD\nA\n', 'A\nB\n\nC\nD\n', 'A\t\nB\nC\nD\n']
        for body in cases:
            with self.subTest(body=body):
                self.write_runs([body])
                result = self.invoke(self.command())
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(self.output.exists())

    def test_unknown_self_duplicate_malformed_and_invalid_mi_fail(self):
        for target, token in [('UNKNOWN', '0.1'), ('A', '0.1'), ('B', 'nan'),
                              ('B', 'inf'), ('B', '-0.1'), ('B', '0.2junk'), ('B', '')]:
            with self.subTest(target=target, token=token):
                self.write_runs(['A\t%s\t%s\nB\nC\nD\n' % (target, token)])
                self.assertNotEqual(self.invoke(self.command()).returncode, 0)
                self.assertFalse(self.output.exists())
        self.write_runs(['A\tB\t0.1\tB\t0.1\nB\nC\nD\n'])
        self.assertIn('duplicate target', self.invoke(self.command()).stderr)
        self.assertFalse(self.output.exists())

    def test_unpruned_dpi_header_is_mandatory_unique_and_equal_one(self):
        headers = ['', '>  DPI tolerance   0\n', '>  DPI tolerance   1\n> DPI tolerance 1\n',
                   '> DPI tolerance nan\n', '> DPI tolerance\n']
        for header in headers:
            with self.subTest(header=header):
                self.write_runs(['A\nB\nC\nD\n'], header=header)
                self.assertNotEqual(self.invoke(self.command()).returncode, 0)
                self.assertFalse(self.output.exists())

    def test_manifest_duplicates_and_hardlinks_are_rejected_but_equal_content_allowed(self):
        self.write_runs(['A\nB\nC\nD\n'] * 2)
        self.assert_success(self.command(k=2))
        self.output.unlink()
        self.manifest.write_text(self.runs[0].name + '\n./' + self.runs[0].name + '\n')
        self.assertIn('Duplicate consensus input', self.invoke(self.command()).stderr)
        alias = self.work / 'hardlink.adj'
        os.link(self.runs[0], alias)
        self.manifest.write_text(self.runs[0].name + '\n' + alias.name + '\n')
        self.assertIn('Duplicate consensus input', self.invoke(self.command()).stderr)
        self.assertFalse(self.output.exists())

    def test_consensus_cannot_rethreshold_resample_or_estimate_mi(self):
        self.write_runs(['A\nB\nC\nD\n'])
        for flags in (['-p', '0.001'], ['-t', '0.1'], ['-u', '80%'], ['-r', '1'],
                      ['-M', 'pchip'], ['-f', '1', '0'], ['-n', '0.1'], ['--raw-mi']):
            with self.subTest(flags=flags):
                self.assertNotEqual(self.invoke(self.command(extra=flags)).returncode, 0)
                self.assertFalse(self.output.exists())

    def test_output_aliases_inputs_and_hardlinks_are_rejected_unchanged(self):
        self.write_runs(['A\nB\nC\nD\n'])
        for target in (self.exp, self.hubs, self.annotation, self.manifest, self.runs[0]):
            with self.subTest(target=target):
                before = target.read_bytes()
                result = self.invoke(self.command(output=target))
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(target.read_bytes(), before)
        alias = self.work / 'alias.adj'
        os.link(self.runs[0], alias)
        before = self.runs[0].read_bytes()
        self.assertNotEqual(self.invoke(self.command(output=alias)).returncode, 0)
        self.assertEqual(self.runs[0].read_bytes(), before)

    def test_failed_input_does_not_replace_an_existing_valid_output(self):
        self.write_runs(['A\nB\nC\n'])
        self.output.write_text('previous good network\n')
        self.assertNotEqual(self.invoke(self.command()).returncode, 0)
        self.assertEqual(self.output.read_text(), 'previous good network\n')
        self.assertFalse(list(self.work.glob('result.adj.tmp.*')))

    @unittest.skipUnless(os.name == 'posix', 'uses POSIX file-size limits')
    def test_failed_output_write_does_not_publish_or_report_success(self):
        import resource
        self.write_runs(['A\nB\nC\nD\n'])
        self.output.write_text('previous good network\n')
        def limit():
            signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
            resource.setrlimit(resource.RLIMIT_FSIZE, (64, 64))
        # Instrumented native programs flush gcov counters on exit. Keep this
        # size-limited process from truncating the shared build's .gcda files.
        environment = dict(os.environ, GCOV_PREFIX=str(self.work / 'limited-gcov'),
                           GCOV_PREFIX_STRIP='0')
        result = self.invoke(self.command(), preexec_fn=limit, env=environment)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Unable to finish writing matrix', result.stderr)
        self.assertEqual(self.output.read_text(), 'previous good network\n')
        self.assertFalse(list(self.work.glob('result.adj.tmp.*')))

    def test_raw_mode_preserves_full_precision_and_hub_order_and_empty_rows(self):
        self.hubs.write_text('B\nA\nB\nMISSING\nC\nD\n')
        base = [NATIVE, '-i', str(self.exp), '-s', str(self.hubs), '-e', '1', '-t', '0', '-N', '40', '-S', '3']
        ordinary = self.work / 'ordinary.adj'
        result = self.invoke(base + ['-o', str(ordinary)])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assert_success(base + ['--raw-mi', '-o', str(self.output)])
        self.assertEqual([row.split('\t')[0] for row in self.rows(self.output)], list('BACD'))
        self.assertEqual(self.header(self.output, 'DPI mode'), 'disabled-raw-mi')
        self.assertEqual(self.header(self.output, 'MI output precision'), '17')
        raw_edges, ordinary_edges = self.edges(self.output), self.edges(ordinary)
        self.assertEqual(set(raw_edges), set(ordinary_edges))
        self.assertTrue(any(raw_edges[key] != ordinary_edges[key] for key in raw_edges))
        for key, value in raw_edges.items():
            self.assertEqual(float(format(value, '.6g')), ordinary_edges[key])
        self.assert_success(base + ['--raw-mi', '-t', '1000', '-o', str(self.output)])
        self.assertEqual(self.rows(self.output), list('BACD'))

    def test_raw_mode_rejects_enabled_dpi_adjacency_and_consensus_settings(self):
        self.write_runs(['A\nB\nC\nD\n'])
        base = [NATIVE, '-i', str(self.exp), '-s', str(self.hubs), '--raw-mi', '-o', str(self.output)]
        for flags in (['-e', '0'], ['--dynamic-dpi'], ['-j', str(self.runs[0])],
                      ['--min-recurrence', '1'], ['--raw-mi=true']):
            with self.subTest(flags=flags):
                self.assertNotEqual(self.invoke(base + flags).returncode, 0)
                self.assertFalse(self.output.exists())


if __name__ == '__main__':
    unittest.main()
