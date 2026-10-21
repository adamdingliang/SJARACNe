"""Synthetic discovery checks: these tests never launch native inference."""

import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import native_test_support as support


class TestNativeTestSupport(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.root = Path(self.folder.name)
        self.environment = mock.patch.dict(os.environ, {}, clear=True)
        self.environment.start()
        self.addCleanup(self.environment.stop)
        self.which = mock.patch.object(support.shutil, "which", return_value=None)
        self.mock_which = self.which.start()
        self.addCleanup(self.which.stop)

    def write_executable(self, path, payload=b"\x7fELFsynthetic"):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
        path.chmod(0o755)
        return path

    def discover(self, **kwargs):
        return support.find_native_executable(project_root=self.root, **kwargs)

    def test_optional_missing_binary_returns_none(self):
        self.assertIsNone(self.discover())

    def test_required_missing_binary_fails_instead_of_skipping(self):
        os.environ[support.REQUIRE_NATIVE_ENV] = "1"
        with self.assertRaisesRegex(RuntimeError, "requires sjaracne.exe"):
            self.discover()

    def test_explicit_missing_override_never_falls_back(self):
        os.environ["SJARACNE_TEST_EXE"] = str(self.root / "missing.exe")
        with self.assertRaisesRegex(RuntimeError, "SJARACNE_TEST_EXE=.*not a file"):
            self.discover()
        self.mock_which.assert_not_called()

    def test_empty_and_directory_overrides_are_errors(self):
        for configured in ("", str(self.root)):
            with self.subTest(configured=configured):
                os.environ["SJARACNE_TEST_EXE"] = configured
                with self.assertRaises(RuntimeError):
                    self.discover()

    def test_explicit_valid_override_is_absolute_and_authoritative(self):
        path = self.write_executable(self.root / "chosen.exe")
        os.environ["SJARACNE_TEST_EXE"] = str(path)
        with mock.patch.object(support, "_runtime_problem", return_value=None):
            self.assertEqual(self.discover(), str(path.resolve()))
        self.mock_which.assert_not_called()

    def test_repository_binary_precedes_path_binary(self):
        local = self.root / "SJARACNe" / "bin" / "sjaracne.exe"
        self.mock_which.return_value = str(self.root / "path.exe")
        with mock.patch.object(support, "_runtime_problem", return_value=None):
            self.assertEqual(self.discover(), str(local.resolve()))

    def test_unusable_local_binary_can_fall_back_to_path(self):
        found = self.root / "path.exe"
        self.mock_which.return_value = str(found)
        with mock.patch.object(support, "_runtime_problem",
                               side_effect=["wrong runtime", None]):
            self.assertEqual(self.discover(), str(found.resolve()))

    def test_generator_uses_its_own_override_and_filename(self):
        os.environ["SJARACNE_APMI_NULL_EXE"] = str(self.root / "generator.exe")
        with mock.patch.object(support, "_runtime_problem", return_value=None):
            self.assertEqual(self.discover(
                name="apmi_null_generator.exe", override_env="SJARACNE_APMI_NULL_EXE"),
                str((self.root / "generator.exe").resolve()))

    def test_runtime_formats_are_checked_consistently(self):
        path = self.root / "binary.exe"
        cases = (
            ("win32", b"MZsynthetic", None),
            ("win32", b"\x7fELFsynthetic", "Windows executable"),
            ("linux", b"\x7fELFsynthetic", None),
            ("linux", b"MZsynthetic", "Linux executable"),
            ("darwin", bytes.fromhex("cffaedfe") + b"synthetic", None),
            ("darwin", b"\x7fELFsynthetic", "macOS executable"),
            ("linux", b"#!/bin/sh\n", None),
        )
        for platform_name, payload, expected in cases:
            with self.subTest(platform=platform_name, payload=payload):
                self.write_executable(path, payload)
                with mock.patch.object(support.sys, "platform", platform_name), \
                        mock.patch.object(support.os, "access", return_value=True):
                    problem = support._runtime_problem(path)
                if expected is None:
                    self.assertIsNone(problem)
                else:
                    self.assertIn(expected, problem)

    def test_non_executable_posix_file_is_unusable(self):
        path = self.write_executable(self.root / "binary.exe")
        with mock.patch.object(support.sys, "platform", "linux"), \
                mock.patch.object(support.os, "access", return_value=False):
            self.assertIn("not executable", support._runtime_problem(path))

    def test_wrong_runtime_explicit_override_fails_even_when_optional(self):
        os.environ["SJARACNE_TEST_EXE"] = str(self.root / "wrong.exe")
        with mock.patch.object(support, "_runtime_problem", return_value="wrong runtime"):
            with self.assertRaisesRegex(RuntimeError, "wrong runtime"):
                self.discover()


if __name__ == "__main__":
    unittest.main()
