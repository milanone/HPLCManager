"""Test of the DAD spectra reader (.uv version 31) and of the extraction by wavelength.

Two levels: a synthetic .uv (no real data or Tk needed) and, if present, the real folder
'example data/009-0201.D' (gitignored): the chromatogram extracted from the .uv at 280/295/320/350 nm
(4 nm band, 360/100 nm reference) must reproduce the dad1A..D.ch signals computed by ChemStation.

    py -m unittest discover -s tests -v
"""
import importlib.util
import os
import struct
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "hplc", os.path.join(HERE, "..", "HPLCManager.pyw"))
hplc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hplc)

SAMPLE = os.path.join(HERE, "..", "example data", "009-0201.D")


def write_uv(path, matrix, t0=240, dt=400, w0=4000, w1=12000, step=40):
    """Write a minimal v31 .uv: zero-filled header + one little-endian record per spectrum."""
    b = bytearray(0x200)
    b[0:3] = b"\x0231"
    for row, ms in zip(matrix, (t0 + dt * np.arange(len(matrix)))):
        entries, prev = [], 0
        for v in row:
            d = int(v) - prev
            if -0x7FFF <= d <= 0x7FFF:
                entries.append(struct.pack("<h", d))
            else:
                entries.append(b"\x00\x80" + struct.pack("<i", int(v)))
            prev = int(v)
        data = b"".join(entries)
        length = 22 + len(data)
        b += struct.pack("<HHI", 0x43, length, int(ms)) + struct.pack("<HHH", w0, w1, step) + bytes(8) + data
    b += bytes(100)    # tail after the last record
    with open(path, "wb") as f:
        f.write(b)


class TestSyntheticReadUv(unittest.TestCase):
    def setUp(self):
        self.app = hplc.HPLCManager.__new__(hplc.HPLCManager)

    def test_deltas_and_absolutes(self):
        rng = np.random.default_rng(1)
        m = rng.integers(-400, 400, size=(5, 201))
        m[2, 50:] += 600000          # jump that needs the absolute value (0x8000 + int32)
        m[3] = np.arange(201) * 3000  # ramp: all deltas fit in int16
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.uv")
            write_uv(p, m)
            t, wl, S, info = self.app.read_uv(p)
        np.testing.assert_allclose(S, m / 2000.0)
        np.testing.assert_allclose(wl, 200 + 2 * np.arange(201))
        np.testing.assert_allclose(t, (240 + 400 * np.arange(5)) / 60000.0)
        self.assertEqual(info["Spectra"], "5 spectra, 200-600 nm, step 2.0 nm")

    def test_unsupported_version(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.uv")
            with open(p, "wb") as f:
                f.write(b"\x02131" + bytes(0x600))
            with self.assertRaises(ValueError):
                self.app.read_uv(p)

    def test_extract_wavelength(self):
        wl = 200 + 2 * np.arange(201, dtype=float)
        S = np.tile(wl, (3, 1)) * 0.0
        S[:, :] = 1.0                   # flat spectrum = 1 mAU
        S[1, :] += 10.0                 # the second spectrum is 11 mAU
        app = self.app
        y = app.extract_wavelength(wl, S, 280, 4)
        np.testing.assert_allclose(y, [1.0, 11.0, 1.0])
        y = app.extract_wavelength(wl, S, 280, 4, 360, 100)    # band minus reference: all flat -> 0
        np.testing.assert_allclose(y, 0.0, atol=1e-12)
        with self.assertRaises(ValueError):
            app.extract_wavelength(wl, S, 900, 4)               # outside the acquired spectrum


@unittest.skipUnless(os.path.isdir(SAMPLE), "sample 009-0201.D not present")
class TestUvWithRealCh(unittest.TestCase):
    CHANNELS = {"dad1A.ch": 280, "dad1B.ch": 295, "dad1C.ch": 320, "dad1D.ch": 350}

    @classmethod
    def setUpClass(cls):
        cls.app = hplc.HPLCManager.__new__(hplc.HPLCManager)
        cls.t, cls.wl, cls.S, cls.info = cls.app.read_uv(os.path.join(SAMPLE, "dad1.uv"))

    def test_structure(self):
        self.assertEqual(self.S.shape, (13493, 201))
        self.assertEqual(self.wl[0], 200.0)
        self.assertEqual(self.wl[-1], 600.0)
        self.assertAlmostEqual(self.t[0], 240 / 60000.0)
        self.assertAlmostEqual(self.t[1] - self.t[0], 400 / 60000.0)

    def test_extraction_like_chemstation(self):
        # the .ch starts at -2560 ms, the .uv at 240 ms: same time scale, 7 samples of offset
        for ch, lam in self.CHANNELS.items():
            with self.subTest(channel=ch):
                tc, yc, _ = self.app.read_ch(os.path.join(SAMPLE, ch))
                y = self.app.extract_wavelength(self.wl, self.S, lam, 4, 360, 100)
                idx = np.round((self.t - tc[0]) * 60000.0 / 400.0).astype(int)
                np.testing.assert_allclose(tc[idx], self.t, atol=1e-6)
                deviation = y - yc[idx]
                # 4 channels, main peak up to 2100 mAU: RMS deviation below 0.3 mAU
                self.assertLess(deviation.std(), 0.3)
                self.assertLess(np.abs(deviation).max(), 8.0)


if __name__ == "__main__":
    unittest.main()
