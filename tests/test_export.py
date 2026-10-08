"""Export and session tests of HPLCManager.

The tables (peaks, spectra) are tried without a window. The real saving (CSV, session, figure)
runs on the application with a hidden Tk window, replacing only the dialogs
(file and message boxes); so they need a screen and the sample 'example data/009-0201.D'
(gitignored): without it, those tests are skipped.

    py -m unittest discover -s tests -v
"""
import gzip
import importlib.util
import os
import pickle
import tempfile
import tkinter as tk
import unittest
from unittest import mock

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "hplc", os.path.join(HERE, "..", "HPLCManager.pyw"))
hplc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hplc)

SAMPLE = os.path.join(HERE, "..", "example data", "009-0201.D")


def trace(peaks=(), cs_peaks=()):
    df = pd.DataFrame({'mAU': np.arange(5.0)}, index=pd.Index(np.arange(5.0), name='Time (min)'))
    return {'df': df, 'info': {}, 'peaks': list(peaks), 'cs_peaks': list(cs_peaks), 'hidden': False}


class TestTables(unittest.TestCase):
    def setUp(self):
        self.app = hplc.HPLCManager.__new__(hplc.HPLCManager)

    def test_peaks_table(self):
        mine = {'n': 1, 'rt': 4.35, 'height': 100.0, 'area': 2000.0, 'pct': 80.0, 'start': 4.2, 'end': 4.5}
        cs = {'n': 9, 'rt': 4.354, 'kind': 'VV', 'width': 0.18, 'area': 23183.0, 'height': 2003.0,
              'area_pct': 41.9}
        self.app.chromatograms = {'A': trace([mine], [cs]), 'B': trace()}
        df = self.app.peaks_table()
        self.assertEqual(list(df['Source']), ['HPLCManager', 'ChemStation'])
        self.assertEqual(list(df['Trace']), ['A', 'A'])
        self.assertEqual(df.loc[0, 'Area (mAU*s)'], 2000.0)
        self.assertEqual(df.loc[1, 'Type'], 'VV')
        self.assertEqual(df.loc[1, 'Width (min)'], 0.18)
        self.assertTrue(np.isnan(df.loc[0, 'Width (min)']))

    def test_peaks_table_empty(self):
        self.app.chromatograms = {'A': trace()}
        self.assertTrue(self.app.peaks_table().empty)

    def test_spectra_table(self):
        wl = np.array([200.0, 202.0, 204.0])
        self.app._chosen_spectra = [
            {'dataset': 'D1', 't': 4.357, 'wl': wl, 'y': np.array([1.0, 2.0, 3.0])},
            {'dataset': 'D1', 't': 8.8, 'wl': wl, 'y': np.array([4.0, 5.0, 6.0])}]
        df = self.app.spectra_table()
        self.assertEqual(df.index.name, 'Wavelength (nm)')
        self.assertEqual(list(df.columns), ['D1 t=4.357 min', 'D1 t=8.800 min'])
        np.testing.assert_allclose(df['D1 t=8.800 min'], [4.0, 5.0, 6.0])

    def test_spectra_table_empty(self):
        self.app._chosen_spectra = []
        self.assertTrue(self.app.spectra_table().empty)


@unittest.skipUnless(os.path.isdir(SAMPLE), "sample 009-0201.D not present")
class TestWithApplication(unittest.TestCase):
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
        self.app.var_include_ch.set(True)                 # also the channels with the ChemStation peaks
        self.app.process_file(SAMPLE)
        self.app.refresh_view()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def tearDown(self):
        for w in self.root.winfo_children():
            w.destroy()

    def path(self, name):
        return os.path.join(self.tmp.name, name)

    def save_with(self, func, name):
        """Run `func` with the 'save as' dialog and the message boxes simulated."""
        p = self.path(name)
        with mock.patch.object(hplc.filedialog, 'asksaveasfilename', return_value=p), \
                mock.patch.object(hplc.messagebox, 'showinfo') as info, \
                mock.patch.object(hplc.messagebox, 'showerror') as err, \
                mock.patch.object(hplc.messagebox, 'showwarning') as warn:
            func()
        self.assertFalse(err.called, err.call_args)
        self.assertFalse(warn.called, warn.call_args)
        return p, info

    def test_export_peaks(self):
        n = [k for k in self.app.chromatograms if k.endswith('dad1A')][0]
        self.app.listbox.selection_clear(0, tk.END)       # a click replaces the selection
        self.app.listbox.selection_set(list(self.app.chromatograms).index(n))
        self.app.open_peaks()
        self.app.var_prom.set('5')
        self.app._pk_apply()
        p, _ = self.save_with(self.app.export_peaks, 'picchi.csv')
        df = pd.read_csv(p, sep=';', encoding='latin-1')
        self.assertEqual(set(df['Source']), {'HPLCManager', 'ChemStation'})
        cs = df[(df['Source'] == 'ChemStation') & (df['Trace'] == n)]
        self.assertEqual(len(cs), 60)                     # peaks of REPORT01.CSV
        self.assertAlmostEqual(cs['Area (mAU*s)'].max(), 23183.0, delta=5)

    def test_export_peaks_without_peaks(self):
        self.app.chromatograms = {'X': trace()}
        with mock.patch.object(hplc.messagebox, 'showinfo') as info:
            self.app.export_peaks()
        self.assertTrue(info.called)

    def test_export_spectra(self):
        self.app._sp_dataset = lambda: self.app.spectra[next(iter(self.app.spectra))]
        self.app.var_sp_set = tk.StringVar(value=next(iter(self.app.spectra)))
        self.app._sp_spectrum_at(4.354)
        self.app._sp_spectrum_at(8.8)
        p, _ = self.save_with(self.app.export_spectra, 'spettri.csv')
        df = pd.read_csv(p, sep=';', encoding='latin-1', index_col=0)
        self.assertEqual(df.shape, (201, 2))
        self.assertEqual(df.index[0], 200.0)
        # the spectrum at 4.35 min has its maximum at ~272 nm
        self.assertAlmostEqual(df.iloc[:, 0].idxmax(), 272, delta=6)

    def test_export_full_dad(self):
        p, _ = self.save_with(self.app.export_full_dad, 'dad.csv')
        df = pd.read_csv(p, sep=';', encoding='latin-1', index_col=0, nrows=3)
        self.assertEqual(df.shape, (3, 201))
        self.assertEqual(list(df.columns[:2]), ['200', '202'])

    def test_session_round_trip(self):
        self.app.var_sp_l.set('254')
        self.app._sp_extract(True)
        before = list(self.app.chromatograms)
        peak = self.app.chromatograms[before[0]]['cs_peaks'][:1]
        p, _ = self.save_with(self.app.save_session, 's.hplcsession')
        with gzip.open(p, 'rb') as f:
            self.assertEqual(pickle.load(f)['version'], 1)
        self.app.clear_all()
        self.assertEqual(self.app.chromatograms, {})
        with mock.patch.object(hplc.filedialog, 'askopenfilename', return_value=p):
            self.app.open_session()
        self.assertEqual(list(self.app.chromatograms), before)
        self.assertEqual(self.app.chromatograms[before[0]]['cs_peaks'][:1], peak)
        self.assertIn(next(iter(self.app.spectra)), self.app.spectra)
        self.assertEqual(self.app.spectra[next(iter(self.app.spectra))]['S'].shape, (13493, 201))
        self.assertTrue(self.app.var_include_ch.get())
        self.assertFalse(self.app._dirty)

    def test_open_invalid_session(self):
        p = self.path('x.hplcsession')
        with gzip.open(p, 'wb') as f:
            pickle.dump({'other': 1}, f)
        with mock.patch.object(hplc.filedialog, 'askopenfilename', return_value=p), \
                mock.patch.object(hplc.messagebox, 'showerror') as err:
            self.app.open_session()
        self.assertTrue(err.called)
        self.assertGreater(len(self.app.chromatograms), 0)   # the current state is not touched

    def test_figure(self):
        self.app._sp_dataset = lambda: self.app.spectra[next(iter(self.app.spectra))]
        self.app.var_sp_set = tk.StringVar(value=next(iter(self.app.spectra)))
        self.app._sp_spectrum_at(4.354)                      # adds a marker line
        self.assertTrue(any(a.get_gid() == 'marker' for a in self.app.ax.lines))
        p, _ = self.save_with(self.app.save_figure_image, 'f.png')
        self.assertGreater(os.path.getsize(p), 5000)
        p, _ = self.save_with(self.app.save_figure_pickle, 'f.fig.pickle')
        with open(p, 'rb') as f:
            fig = pickle.load(f)
        self.assertEqual(len(fig.axes), 1)
        # the spectra markers and the cursor do not end up in the saved figure...
        self.assertFalse(any(a.get_gid() == 'marker' for a in fig.axes[0].lines))
        # ...and instead stay on the live view
        self.assertTrue(any(a.get_gid() == 'marker' for a in self.app.ax.lines))

    def test_figure_editor(self):
        if not os.path.isfile(os.path.join(HERE, "..", "..", "PlotStyleKit", "plot_editor.pyw")):
            self.skipTest("PlotStyleKit not present")
        with mock.patch.object(hplc.messagebox, 'showerror') as err:
            self.app.open_figure_editor()
        self.assertFalse(err.called, err.call_args)
        self.assertTrue(any(isinstance(w, tk.Toplevel) for w in self.root.winfo_children()))


if __name__ == "__main__":
    unittest.main()
