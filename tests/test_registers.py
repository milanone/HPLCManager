"""Tests of the ChemStation register readers: LCDIAG.REG (instrument curves) and ACQRES.REG
(modules and column), with the curves window.

The curves are tried on a synthetic file (same structure) and on the real sample
'example data/009-0201.D' (gitignored), comparing with the values written by ChemStation in
Report00.CSV / RUN.LOG; without the sample the related tests are skipped.

    py -m unittest discover -s tests -v
"""
import importlib.util
import os
import struct
import tempfile
import tkinter as tk
import unittest
from unittest import mock

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "hplc", os.path.join(HERE, "..", "HPLCManager.pyw"))
hplc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hplc)

SAMPLE = os.path.join(HERE, "..", "example data", "009-0201.D")


def write_diag(path, signals):
    """Write a synthetic LCDIAG.REG: for each signal (title, unit, values, dt_min, factor)
    prefix with N and dt, 'min\\0unit\\0', uint32 array and header with the title."""
    b = bytearray(b"\x02" + b"32\x00" + b"\x0dREGISTER FILE" + bytes(40))
    for title, unit, values, dt, factor in signals:
        pre = bytearray(400)                         # 400-byte prefix: N and dt at fixed offsets
        base = len(b) + 400                          # position of 'min' in the file
        b += pre
        b[base - 168:base - 164] = struct.pack("<I", len(values))
        b[base - 116:base - 108] = struct.pack("<d", dt)
        b += b"min\x00" + unit.encode("latin-1") + b"\x00"
        b += np.round(np.asarray(values) * factor).astype("<u4").tobytes()
        b += b"ObjClass" + bytes(30) + b"Title" + bytes(20) + title.encode() + b"\x00arial\x00" + bytes(60)
    with open(path, "wb") as f:
        f.write(b)


class TestSyntheticDiagnostics(unittest.TestCase):
    def setUp(self):
        self.app = hplc.HPLCManager.__new__(hplc.HPLCManager)
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_scales_and_times(self):
        p = os.path.join(self.tmp.name, "LCDIAG.REG")
        write_diag(p, [("PMP1, Pressure", "bar", [94.42, 95.0, 100.55], 0.005, 100),
                        ("PMP1, Flow", "ml/min", [1.0, 1.0, 0.95], 0.005, 1000),
                        ("PMP1, Solvent B", "%", [0.0, 12.5, 99.9], 0.005, 10),
                        ("THM1, Temperature (Left)", "\xb0C", [29.91, 30.0], 1 / 60, 100)])
        seg = self.app.read_diagnostics(p)
        self.assertEqual(list(seg), ["PMP1, Pressure", "PMP1, Flow", "PMP1, Solvent B", "THM1, Temperature (Left)"])
        np.testing.assert_allclose(seg["PMP1, Pressure"]["y"], [94.42, 95.0, 100.55])
        np.testing.assert_allclose(seg["PMP1, Flow"]["y"], [1.0, 1.0, 0.95])
        np.testing.assert_allclose(seg["PMP1, Solvent B"]["y"], [0.0, 12.5, 99.9])
        np.testing.assert_allclose(seg["THM1, Temperature (Left)"]["y"], [29.91, 30.0])
        np.testing.assert_allclose(seg["PMP1, Pressure"]["t"], [0, 0.005, 0.01])
        self.assertAlmostEqual(seg["THM1, Temperature (Left)"]["t"][1], 1 / 60)
        self.assertEqual(seg["PMP1, Pressure"]["unit"], "bar")

    def test_not_a_register_file(self):
        p = os.path.join(self.tmp.name, "x.REG")
        with open(p, "wb") as f:
            f.write(b"niente" * 100)
        with self.assertRaises(ValueError):
            self.app.read_diagnostics(p)
        with self.assertRaises(ValueError):
            self.app.read_acqres(p)

    def test_instrument_table(self):
        t_fast = np.arange(0, 3) * 0.005
        self.app.instrument = {'D': {'signals': {
            'PMP1, Pressure': {'t': t_fast, 'y': np.array([1.0, 2.0, 3.0]), 'unit': 'bar'},
            'THM1, Temperature (Left)': {'t': np.array([0.0, 0.01]), 'y': np.array([30.0, 31.0]), 'unit': '\xb0C'}}}}
        df = self.app.instrument_table('D')
        self.assertEqual(df.index.name, 'Time (min)')
        np.testing.assert_allclose(df['PMP1, Pressure [bar]'], [1, 2, 3])
        # the temperature (slower) is interpolated onto the time scale of the fast signal
        np.testing.assert_allclose(df['THM1, Temperature (Left) [\xb0C]'], [30.0, 30.5, 31.0])


@unittest.skipUnless(os.path.isdir(SAMPLE), "sample 009-0201.D not present")
class TestRealSample(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = hplc.HPLCManager.__new__(hplc.HPLCManager)
        cls.seg = cls.app.read_diagnostics(os.path.join(SAMPLE, "LCDIAG.REG"))
        cls.acq = cls.app.read_acqres(os.path.join(SAMPLE, "ACQRES.REG"))

    def test_signals_present(self):
        self.assertEqual(list(self.seg), ["PMP1, Pressure", "PMP1, Flow", "PMP1, Solvent A", "PMP1, Solvent B",
                                          "PMP1, Solvent C", "PMP1, Solvent D", "THM1, Temperature (Left)"])
        for t in ("PMP1, Pressure", "PMP1, Flow", "PMP1, Solvent A"):
            self.assertEqual(len(self.seg[t]['y']), 18000)
            self.assertAlmostEqual(self.seg[t]['t'][1], 0.005)
        self.assertEqual(len(self.seg["THM1, Temperature (Left)"]['y']), 5400)
        self.assertAlmostEqual(self.seg["PMP1, Pressure"]['t'][-1], 89.995)    # 90 min run

    def test_consistent_with_report00_and_run_log(self):
        m = self.app.read_folder_metadata(SAMPLE)
        p = self.seg["PMP1, Pressure"]
        # Report00.CSV: Start Pressure 96.449997 bar (a few seconds after injection), Stop 104.720001 bar
        self.assertAlmostEqual(np.interp(4 / 60, p['t'], p['y']), 96.45, delta=1.0)
        self.assertAlmostEqual(p['y'][-1], 104.72, delta=0.5)
        self.assertTrue(m['Start Pressure'].startswith('96.449997'))
        # flow 1 ml/min from start to end
        np.testing.assert_allclose(self.seg["PMP1, Flow"]['y'], 1.0)
        # RUN.LOG: column temperature 29.9-30.0 C; the curve is in the same range
        temp = self.seg["THM1, Temperature (Left)"]['y']
        self.assertGreaterEqual(temp.min(), 29.8)
        self.assertLessEqual(temp.max(), 30.2)
        self.assertEqual(self.seg["THM1, Temperature (Left)"]['unit'], '\xb0C')

    def test_gradient_sums_to_100(self):
        tot = sum(self.seg["PMP1, Solvent %s" % c]['y'] for c in "ABCD")
        self.assertGreater(tot.min(), 99.0)
        self.assertLess(tot.max(), 101.0)
        self.assertEqual(self.seg["PMP1, Solvent A"]['y'][0], 100.0)       # starts at 100% A
        self.assertEqual(self.seg["PMP1, Solvent B"]['y'][0], 0.0)
        self.assertGreater(self.seg["PMP1, Solvent B"]['y'].max(), 60.0)   # gradient towards B

    def test_modules_and_column(self):
        names = [m['name'] for m in self.acq['modules']]
        self.assertEqual(names, ['1100 Quaternary Pump', '1100 Autosampler', '1100 Diode Array Detector',
                                '1100 Column Thermostat'])
        det = self.acq['modules'][2]
        self.assertEqual((det['part'], det['firmware'], det['build']), ('G1315B', 'A.05.09', '007'))
        self.assertRegex(det['serial'], r'^DE\d{8}$')          # serial number: not reported in the repo
        self.assertEqual(self.acq['modules'][0]['part'], 'G1311A')
        c = self.acq['column']
        self.assertEqual(c['description'], 'ODS Hypersil')
        self.assertEqual((c['length_mm'], c['diameter_mm'], c['particle_um']), (100.0, 2.1, 5.0))

    def test_metadata(self):
        m = self.app.read_folder_metadata(SAMPLE)
        self.assertEqual(m['Column'], 'ODS Hypersil')
        self.assertEqual(m['Column size'], '100 x 2.1 mm, 5 um')
        self.assertRegex(m['1100 Diode Array Detector'], r'^G1315B, S/N DE\d{8}, FW A\.05\.09$')


@unittest.skipUnless(os.path.isdir(SAMPLE), "sample 009-0201.D not present")
class TestApplication(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.root = tk.Tk()
        except tk.TclError as e:
            raise unittest.SkipTest("Tk non disponibile: %s" % e)
        cls.root.withdraw()

    @classmethod
    def tearDownClass(cls):
        cls.root.destroy()

    def setUp(self):
        self.app = hplc.HPLCManager(self.root)
        self.app.process_file(SAMPLE)
        self.app.refresh_view()

    def tearDown(self):
        for w in self.root.winfo_children():
            w.destroy()

    def test_loading(self):
        self.assertIn('009-0201', self.app.instrument)
        d = self.app.instrument['009-0201']
        self.assertEqual(len(d['signals']), 7)
        self.assertEqual(len(d['modules']), 4)
        # the metadata of the extracted trace show the column
        n = list(self.app.chromatograms)[0]
        self.assertEqual(self.app.chromatograms[n]['info']['Column'], 'ODS Hypersil')

    def test_curves_window(self):
        w = self.app.open_instrument_curves()
        self.root.update()
        labels = [ax.get_ylabel() for ax in w.fig.axes]
        self.assertEqual(labels, ['Pressure (bar)', 'Solvent composition (%)', 'Flow (ml/min)',
                                     'Temperature (\xb0C)'])
        # solvents C and D stay at zero: the composition panel has only A and B
        self.assertEqual([l.get_label() for l in w.fig.axes[1].lines], ['Solvent A', 'Solvent B'])
        self.assertEqual(len(w.fig.axes[0].lines[0].get_xdata()), 18000)

    def test_export_curves(self):
        import pandas as pd
        p = os.path.join(tempfile.mkdtemp(), 'curve.csv')
        with mock.patch.object(hplc.filedialog, 'asksaveasfilename', return_value=p), \
                mock.patch.object(hplc.messagebox, 'showinfo'), \
                mock.patch.object(hplc.messagebox, 'showerror') as err:
            self.app.export_instrument_curves()
        self.assertFalse(err.called)
        df = pd.read_csv(p, sep=';', encoding='latin-1', index_col=0)
        self.assertEqual(df.shape, (18000, 7))
        self.assertAlmostEqual(df['PMP1, Pressure [bar]'].iloc[0], 94.42)
        self.assertAlmostEqual(df['PMP1, Flow [ml/min]'].iloc[100], 1.0)

    def test_session_carries_the_curves(self):
        import pickle
        data = pickle.loads(pickle.dumps(self.app._session_data()))   # come passando da file
        self.assertIn('009-0201', data['instrument'])
        self.app.clear_all()
        self.assertEqual(self.app.instrument, {})
        self.app._apply_session(data)
        self.assertIn('009-0201', self.app.instrument)

    def test_without_curves(self):
        self.app.instrument.clear()
        with mock.patch.object(hplc.messagebox, 'showinfo') as info:
            self.assertIsNone(self.app.open_instrument_curves())
        self.assertTrue(info.called)


if __name__ == "__main__":
    unittest.main()
