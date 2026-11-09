#!/usr/bin/env python3
import unittest
import tempfile
from pathlib import Path
from SJARACNe.bin.ch_line_ending import ch_line_ending as ch


class TestLineEnding(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.input = Path(self.folder.name) / 'input.txt'
        self.output = Path(self.folder.name) / 'output.txt'

    def test_infile_same_as_outfile(self):
        content = b'Hello world!\r\nThis is an example file\r\n'
        self.input.write_bytes(content)
        with self.assertRaises(SystemExit) as err:
            ch(str(self.input), str(self.input))
        self.assertEqual(err.exception.code, 'Error - you must omit output file argument if it is identical to input file')
        self.assertEqual(self.input.read_bytes(), content)

    def test_first_line_is_invalid(self):
        # Write and close the actual invalid content before the function reads it.
        # The former buffered fixture accidentally tested an empty input instead.
        content = b'Hello world!This is an example file'
        self.input.write_bytes(content)
        with self.assertRaises(SystemExit) as err:
            ch(str(self.input), str(self.output))
        self.assertEqual(err.exception.code, 'Error - invalid line ending in the first line: {}'.format(content))
        self.assertFalse(self.output.exists())


    def assert_conversion(self, ending):
        self.input.write_bytes(b'Hello world!' + ending + b'This is an example file')
        output_file = ch(str(self.input), str(self.output))
        self.assertEqual(Path(output_file).read_bytes(), b'Hello world!\nThis is an example file')

    def test_unix_ending_file(self):
        self.assert_conversion(b'\n')

    def test_windows_ending_file(self):
        self.assert_conversion(b'\r\n')

    def test_mac_ending(self):
        self.assert_conversion(b'\r')

if __name__ == "__main__":
    unittest.main()
