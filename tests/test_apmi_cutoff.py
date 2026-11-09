"""Standalone PCHIP runtime tests; no benchmark checkout or raw data required."""
from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import warnings

import SJARACNe.apmi_cutoff as cutoff
from SJARACNe import get_mi_cutoff


ROOT = Path(__file__).resolve().parents[1]
# Literal outputs of the already frozen benchmark callable, not fitted values.
# Source primary SHA256: 8b4104be3f3ffe0a476a4dd7c8f7f6d95f00f1a51fa467c7194edf3806d427f5
# NumPy/SciPy log/exp/PPoly and stdlib/native arithmetic can differ a few ULPs
# away from literal empirical knots; the snapshots measure that explicitly.
GOLDEN = (
    (50, 0.05, "0x1.0f793acaf4b60p-4"),
    (100, 0.01, "0x1.2336cdb370380p-4"),
    (350, 0.05, "0x1.6d0f073279deap-7"),
    (350, 0.0034641016151378, "0x1.dec9329a3efafp-6"),
    (450, 0.001, "0x1.ef97cc8c9b600p-6"),
    (644, 0.001, "0x1.60f6c092d59f0p-6"),
    (725, 0.0005, "0x1.64abb29dffb29p-6"),
    (1000, 0.0005, "0x1.0509973b2fa00p-6"),
    (5000, 0.0034641016151378, "0x1.0aaba3c3fff2cp-9"),
    (10000, 0.05, "0x1.a6ed352c28000p-12"),
    (80, 0.00075, "0x1.3d67265e603cep-3"),
    (90, 0.002, "0x1.e2b2e26fcad11p-4"),
    (516, 0.0002, "0x1.1e6cd4a25d336p-5"),
    (644, 0.0001, "0x1.045182e05b36ap-5"),
    (725, 0.00005, "0x1.f82c25bbaf175p-6"),
    (1000, 0.000005, "0x1.dbeb4c304a673p-6"),
    (10000, 1e-7, "0x1.231d452bbe800p-8"),
    (50, 1e-7, "0x1.19e1d4d279964p-1"),
    (350, 1e-7, "0x1.a94d4abf500bap-4"),
    (9850, 0.0005, "0x1.abc66ae609e74p-10"),
)


def bits(value):
    return struct.unpack(">Q", struct.pack(">d", value))[0]


class TestStandaloneApmiCutoff(unittest.TestCase):
    def test_public_package_export(self):
        self.assertIs(get_mi_cutoff, cutoff.get_mi_cutoff)

    def test_legacy_star_exports_remain_available(self):
        # Test the historical exports in a clean interpreter. Importing a new
        # submodule during test discovery adds incidental attributes to bin;
        # those are not new public star-import compatibility requirements.
        script = '''
scope = {}
exec('from SJARACNe import *', scope)
import SJARACNe.bin as legacy
from SJARACNe import get_mi_cutoff
for name in vars(legacy):
    if not name.startswith('_'):
        assert scope[name] is getattr(legacy, name), name
assert scope['get_mi_cutoff'] is get_mi_cutoff
'''
        result = subprocess.run([sys.executable, '-c', script], cwd=str(ROOT),
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_bundled_provenance_and_shared_native_data_hash(self):
        data = cutoff._load_data()
        self.assertEqual(data["NM"], 180)
        self.assertEqual(data["NP"], 25)
        self.assertEqual(data["npar"], 40)
        self.assertEqual(data["status"], "experimental-not-uniformly-validated")
        self.assertEqual(data["PRIMARY_SHA256"], "8b4104be3f3ffe0a476a4dd7c8f7f6d95f00f1a51fa467c7194edf3806d427f5")
        self.assertTrue(data["no_new_draws_fitting_selection_or_recalibration"])
        header = (ROOT / "SJARACNe" / "src" / "apmi_pchip_data.h").read_text(encoding="ascii")
        self.assertIn('DATA_SHA256[] = "' + cutoff._DATA_SHA256 + '"', header)
        self.assertIn('PRIMARY_SHA256[] = "' + data["PRIMARY_SHA256"] + '"', header)
        self.assertEqual(hashlib.sha256(cutoff._DATA_PATH.read_bytes()).hexdigest(), cutoff._DATA_SHA256)
        self.assertEqual(hashlib.sha256((ROOT / "tools" / "export_apmi_pchip.py").read_bytes()).hexdigest(), data["exporter_sha256"])

    def test_every_literal_empirical_knot_has_exact_binary64_value(self):
        data = cutoff._load_data()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for i, m in enumerate(data["M"]):
                for j, p in enumerate(data["P"]):
                    value = get_mi_cutoff(m, p)
                    self.assertEqual(value.hex(), data["OBSERVED"][j][i].hex(), (m, p))

    def test_frozen_benchmark_golden_snapshots(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for m, p, hexadecimal in GOLDEN:
                actual = get_mi_cutoff(m, p)
                expected = float.fromhex(hexadecimal)
                self.assertLessEqual(abs(bits(actual) - bits(expected)), 8, (m, p, actual.hex(), hexadecimal))

    def test_all200_core_predictions_finite_positive_and_in_support(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for m in range(50, 10001, 50):
                values = [get_mi_cutoff(m, p) for p in (0.05, 0.01, 0.001, 0.0005)]
                self.assertTrue(all(0 < value <= math.log(m) for value in values))
                self.assertEqual(values, sorted(values))

    def test_strict_m_domain_and_integer_semantics(self):
        for m in (True, False, 49, 10001, 80.0, "80", float("nan"), float("inf")):
            with self.assertRaises(ValueError, msg=str(m)):
                get_mi_cutoff(m, 0.001)

    def test_strict_probability_and_npar(self):
        for p in (True, 0, -0.001, 0.05001, 0.1, 1e-8, float("nan"), float("inf"), "0.001"):
            with self.assertRaises(ValueError, msg=str(p)):
                get_mi_cutoff(100, p)
        for npar in (20, 39, 41, True, 40.0, "40"):
            with self.assertRaises(ValueError, msg=str(npar)):
                get_mi_cutoff(100, 0.001, npar=npar)

    def test_full_probability_range_requires_no_acknowledgment_flag(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for p in (0.0004, 0.0002, 0.0001, 1e-5, 1e-6, 1e-7):
                self.assertGreater(get_mi_cutoff(100, p), 0)

    def test_probability_required_and_removed_flag_rejected(self):
        with self.assertRaises(TypeError):
            get_mi_cutoff(100)
        with self.assertRaises(TypeError):
            get_mi_cutoff(100, 0.001, allow_rare_tail=True)

    def test_experimental_off_grid_and_rare_warnings(self):
        with warnings.catch_warnings(record=True) as recorded:
            warnings.simplefilter("always")
            get_mi_cutoff(644, 0.0001)
        messages = [str(item.message) for item in recorded]
        self.assertEqual(len(messages), 3)
        self.assertTrue(any("not uniformly validated" in message for message in messages))
        self.assertTrue(any("off-grid m" in message for message in messages))
        self.assertTrue(any("outside the principal assessment range" in message for message in messages))
        self.assertFalse(any("Extreme-tail Monte Carlo precision" in message for message in messages))

    def test_probability_0004_has_range_warning_not_extreme_tail_warning(self):
        with warnings.catch_warnings(record=True) as recorded:
            warnings.simplefilter("always")
            value = get_mi_cutoff(100, 0.0004)
        self.assertGreater(value, 0)
        messages = [str(item.message) for item in recorded]
        self.assertEqual(len(messages), 2)
        self.assertTrue(any("outside the principal assessment range" in message for message in messages))
        self.assertFalse(any("Extreme-tail Monte Carlo precision" in message for message in messages))

    def test_extreme_tail_warning_at_and_below_1e5(self):
        for p in (1e-5, 1e-6, 1e-7):
            with warnings.catch_warnings(record=True) as recorded:
                warnings.simplefilter("always")
                value = get_mi_cutoff(100, p)
            self.assertGreater(value, 0)
            messages = [str(item.message) for item in recorded]
            self.assertEqual(len(messages), 3)
            self.assertTrue(any("outside the principal assessment range" in message for message in messages))
            self.assertTrue(any("Extreme-tail Monte Carlo precision" in message and "100, 10 and 1" in message for message in messages))

    def test_pchip_linear_plateau_and_no_extrapolation(self):
        self.assertEqual(cutoff._pchip_value([0, 1], [0, 2], 0.25), 0.5)
        self.assertEqual(cutoff._pchip_value([0, 1, 2], [0, 1, 2], 0.5), 0.5)
        self.assertEqual(cutoff._pchip_value([0, 1, 2, 3], [0, 1, 1, 2], 1.5), 1)
        for query in (-0.1, 3.1):
            with self.assertRaises(ValueError):
                cutoff._pchip_value([0, 1, 2, 3], [0, 1, 1, 2], query)

    def test_body_crossings_rejected_without_sorting(self):
        data = copy.deepcopy(cutoff._load_data())
        data["OBSERVED"][1][0] = data["OBSERVED"][0][0] / 2
        with self.assertRaisesRegex(ValueError, "cross"):
            cutoff._predict(data, 50, 0.001)

    def test_rare_only_failure_does_not_break_moderate_scope(self):
        data = copy.deepcopy(cutoff._load_data())
        data["OBSERVED"][-1][0] = -1.0
        self.assertGreater(cutoff._predict(data, 50, 0.001), 0)
        with self.assertRaises(ValueError):
            cutoff._predict(data, 50, 1e-7)

    def test_shared_boundary_failure_invalidates_both_scopes(self):
        data = copy.deepcopy(cutoff._load_data())
        data["OBSERVED"][21][0] = 0.0
        for p in (0.001, 1e-7):
            with self.assertRaises(ValueError):
                cutoff._predict(data, 50, p)

    def test_corrupt_data_hash_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tampered.json"
            path.write_bytes(b"{}")
            cutoff._load_data.cache_clear()
            try:
                with patch.object(cutoff, "_DATA_PATH", path):
                    with self.assertRaisesRegex(ValueError, "hash mismatch"):
                        cutoff._load_data()
            finally:
                cutoff._load_data.cache_clear()

    def test_standalone_import_and_call_without_numpy_scipy_or_benchmark(self):
        script = """
import builtins, sys, warnings
original = builtins.__import__
blocked = {'numpy', 'scipy', 'pandas', 'sjaracne_benchmark'}
def checked(name, *args, **kwargs):
    if name.split('.')[0] in blocked:
        raise RuntimeError('forbidden runtime dependency: ' + name)
    return original(name, *args, **kwargs)
builtins.__import__ = checked
from SJARACNe import get_mi_cutoff
warnings.simplefilter('ignore')
assert get_mi_cutoff(80, 0.001) > 0
assert not any(name.split('.')[0] in blocked for name in sys.modules)
"""
        result = subprocess.run([sys.executable, "-c", script], cwd=str(ROOT), capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
