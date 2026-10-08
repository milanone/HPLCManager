"""Tests for the ChemStation reader (.ch version 30 and REPORTnn.CSV) of HPLCManager.

Two levels: a synthetic .ch file built with the same delta/absolute scheme (no real
data or Tk needed) and, if present, the real folder 'example data/009-0201.D'
(gitignored) compared with the peak tables written by ChemStation.

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


def write_ch(path, values, t0=-2560, dt=400, sample=b"Test sample"):
    """Write a minimal v30 .ch: zero-filled header + data in segments of 25 entries."""
    b = bytearray(0x400)
    b[0:3] = b"\x0230"
    b[0x18:0x19 + len(sample)] = bytes([len(sample)]) + sample
    b[0x11A:0x122] = struct.pack(">ii", t0, t0 + dt * (len(values) - 1))
    entries, prev = [], 0
    for v in values:
        d = v - prev
        if -0x7FFF <= d <= 0x7FFF:
            entries.append(struct.pack(">h", d))
        else:
            entries.append(b"\x80\x00" + struct.pack(">i", v))
        prev = v
    for i in range(0, len(entries), 25):
        block = entries[i:i + 25]
        b += bytes([0x10, len(block)]) + b"".join(block)
    with open(path, "wb") as f:
        f.write(b)


class TestSyntheticReadCh(unittest.TestCase):
    def setUp(self):
        self.app = hplc.HPLCManager.__new__(hplc.HPLCManager)

    def test_deltas_and_absolutes(self):
        # the jump from 10 to 200000 does not fit in an int16: it must use the 0x8000 escape
        values = [0, 5, 10, 200000, 200010, 199990, 3, 0] * 7
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.ch")
            write_ch(p, values)
            t, y, info = self.app.read_ch(p)
        np.testing.assert_allclose(y, np.array(values) / 2000.0)
        self.assertEqual(len(t), len(values))
        self.assertAlmostEqual(t[0], -2560 / 60000.0)
        self.assertAlmostEqual(t[1] - t[0], 400 / 60000.0)
        self.assertEqual(info["Sample"], "Test sample")

    def test_unsupported_version(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.ch")
            with open(p, "wb") as f:
                f.write(b"\x02130" + bytes(0x500))
            with self.assertRaises(ValueError):
                self.app.read_ch(p)


@unittest.skipUnless(os.path.isdir(SAMPLE), "sample 009-0201.D not present")
class TestWithChemStationReport(unittest.TestCase):
    """Every peak of REPORTnn.CSV must fall where the .ch signal has its apex."""

    CHANNELS = {"dad1A.ch": "REPORT01.CSV", "dad1B.ch": "REPORT02.CSV", "dad1C.ch": "REPORT03.CSV",
              "dad1D.ch": "REPORT04.CSV", "dad1E.ch": "REPORT05.CSV"}

    def setUp(self):
        self.app = hplc.HPLCManager.__new__(hplc.HPLCManager)

    def test_header(self):
        t, y, info = self.app.read_ch(os.path.join(SAMPLE, "dad1A.ch"))
        self.assertEqual(len(y), 13501)
        self.assertAlmostEqual(t[1] - t[0], 400 / 60000.0)
        self.assertAlmostEqual(t[-1], 5397440 / 60000.0)   # about 89.96 min
        self.assertEqual(info["Sample"], "Estr acq gambe")
        self.assertEqual(info["Module"], "G1315B")
        self.assertEqual(info["Method"], "POLIFENB.M")

    def test_tall_report_peaks(self):
        # Only the tall peaks and before 30 min: there the ChemStation baseline lies a few mAU
        # below the signal. Small peaks (integration baseline below the local minimum) and
        # those beyond 80 min (negative drift at the tail of the gradient) cannot be compared
        # without replicating the integrator.
        checked = 0
        for ch, rep in self.CHANNELS.items():
            t, y, _ = self.app.read_ch(os.path.join(SAMPLE, ch))
            for p in self.app.read_report_csv(os.path.join(SAMPLE, rep)):
                if p["height"] < 50 or p["rt"] > 30:
                    continue
                with self.subTest(channel=ch, peak=p["n"], rt=p["rt"]):
                    checked += 1
                    m = np.abs(t - p["rt"]) <= 0.05
                    i = np.flatnonzero(m)[np.argmax(y[m])]
                    self.assertAlmostEqual(t[i], p["rt"], delta=0.01)   # apex at the report RT
                    ratio = y[i] / p["height"]
                    self.assertTrue(1.0 <= ratio <= 1.15, ratio)  # signal = height + baseline
        self.assertGreater(checked, 20)

    def test_main_peak(self):
        t, y, _ = self.app.read_ch(os.path.join(SAMPLE, "dad1A.ch"))
        p = max(self.app.read_report_csv(os.path.join(SAMPLE, "REPORT01.CSV")),
                key=lambda q: q["height"])
        self.assertAlmostEqual(p["rt"], 4.354, places=3)
        self.assertAlmostEqual(y.max(), 2100.5, delta=1.0)
        self.assertAlmostEqual(t[np.argmax(y)], p["rt"], delta=0.01)


class TestFolderMetadata(unittest.TestCase):
    def setUp(self):
        self.app = hplc.HPLCManager.__new__(hplc.HPLCManager)

    def test_synthetic(self):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "Report00.CSV"), "w", encoding="latin-1", newline="") as f:
                f.write('"Sample Name","Sample X",""\r\n"Injection Date","30-Sep-26, 11:08:12",""\r\n'
                        '"Acq. Operator","30-Sep-26, 11:08:12",""\r\n"Location",9,""\r\n'
                        '"Inj Volume",20,"\xb5l"\r\n"Start Flow",1,"ml/min"\r\n'
                        '"Solvent 1","PMP1, Solvent A","H2O"\r\n"Column 1","Peak Number",""\r\n'
                        '"Signal 1","DAD1 A, Sig=280,4 Ref=360,100",""\r\n')
            with open(os.path.join(d, "RUN.LOG"), "w", encoding="latin-1") as f:
                f.write("1100 THM   1 Column temperature = 29.9 \xb0C   11:08:16 09/30/26\r\n"
                        "1100 THM   1 Column temperature = 30.1 \xb0C   12:38:16 09/30/26\r\n")
            m = self.app.read_folder_metadata(d)
        self.assertEqual(m["Sample Name"], "Sample X")
        self.assertEqual(m["Location"], "9")
        self.assertEqual(m["Inj Volume"], "20 \xb5l")
        self.assertEqual(m["Start Flow"], "1 ml/min")
        self.assertEqual(m["Solvent 1"], "PMP1, Solvent A, H2O")
        self.assertEqual(m["Signal 1"], "DAD1 A, Sig=280,4 Ref=360,100")
        self.assertNotIn("Column 1", m)
        self.assertNotIn("Acq. Operator", m)       # it holds the injection date, not the operator
        self.assertEqual(m["Column temperature"], "29.9-30.1 \xb0C")

    def test_empty_folder(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(self.app.read_folder_metadata(d), {})

    @unittest.skipUnless(os.path.isdir(SAMPLE), "sample 009-0201.D not present")
    def test_real_sample(self):
        m = self.app.read_folder_metadata(SAMPLE)
        self.assertEqual(m["Sample Name"], "Estr acq gambe")
        self.assertEqual(m["Location"], "9")
        self.assertEqual(m["Inj Volume"], "20 \xb5l")
        self.assertEqual(m["Start Pressure"], "96.449997 bar")
        self.assertEqual(m["Solvent 2"], "PMP1, Solvent B, ACN")
        self.assertEqual(m["Signal 5"], "DAD1 E, Sig=520,4 Ref=off")
        self.assertTrue(m["Data File"].endswith("009-0201.D"))
        self.assertEqual(m["Column temperature"], "29.9-30.0 \xb0C")


class TestReadReportCsv(unittest.TestCase):
    def test_row(self):
        app = hplc.HPLCManager.__new__(hplc.HPLCManager)
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "REPORT01.CSV")
            with open(p, "w", encoding="latin-1", newline="") as f:
                f.write('1,2.78,"VV  ",0.085,100.7,16.9,0.18\r\n')
            r = app.read_report_csv(p)
        self.assertEqual(r, [{"n": 1, "rt": 2.78, "kind": "VV", "width": 0.085,
                              "area": 100.7, "height": 16.9, "area_pct": 0.18}])


if __name__ == "__main__":
    unittest.main()
