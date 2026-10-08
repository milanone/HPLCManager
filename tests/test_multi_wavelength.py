"""Test of the extraction of several wavelengths in one go (the 'Several wavelengths' dialog).

`extract_series` is tried without a window on synthetic spectra; the real dialog (number of rows,
wavelengths, reference per row, errors) runs with a hidden Tk window on the sample
'example data/009-0201.D' (gitignored), if present.

    py -m unittest discover -s tests -v
"""
import importlib.util
import os
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


def synthetic_app():
    app = hplc.HPLCManager.__new__(hplc.HPLCManager)
    wl = 200 + 2 * np.arange(201, dtype=float)
    t = np.arange(10) / 60.0
    S = np.outer(np.arange(10) + 1.0, np.ones(201))      # flat spectrum, scales with time
    S[:, wl > 400] = 0.0                                    # above 400 nm it is 0
    app.chromatograms = {}
    app.spectra = {'D': {'t': t, 'wl': wl, 'S': S, 'info': {'Sample': 'X', 'Spectra': 's', 'Version': '31'}}}
    return app


class TestExtractSeries(unittest.TestCase):
    def test_several_traces(self):
        app = synthetic_app()
        names = app.extract_series('D', [(250, 4, None, None), (300, 8, 500, 20), (350, 4, 500, 20)])
        self.assertEqual(names, ['D 250,4 no ref', 'D 300,8 ref 500,20', 'D 350,4 ref 500,20'])
        self.assertEqual(list(app.chromatograms), names)
        # flat spectrum = (i+1) up to 400 nm, 0 beyond: with a reference at 500 nm the signal does not change
        for n in names:
            np.testing.assert_allclose(app.chromatograms[n]['df']['mAU'], np.arange(10) + 1.0)
        self.assertEqual(app.chromatograms[names[0]]['info']['Sample'], 'X')
        self.assertNotIn('Spectra', app.chromatograms[names[0]]['info'])

    def test_all_or_nothing(self):
        app = synthetic_app()
        with self.assertRaisesRegex(ValueError, "Row 2"):
            app.extract_series('D', [(250, 4, None, None), (900, 4, None, None)])   # outside the spectrum
        self.assertEqual(app.chromatograms, {})
        with self.assertRaisesRegex(ValueError, "Row 1"):
            app.extract_series('D', [(250, 0, None, None)])                         # zero bandwidth

    def test_same_wavelength_twice(self):
        app = synthetic_app()
        names = app.extract_series('D', [(250, 4, None, None), (250, 4, None, None)])
        self.assertEqual(len(set(names)), 2)       # the names stay distinct


@unittest.skipUnless(os.path.isdir(SAMPLE), "sample 009-0201.D not present")
class TestDialog(unittest.TestCase):
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

    def fill_in(self, top, values):
        top._var_n.set(len(values))
        top._set_n()
        self.assertEqual(len(top._rows), len(values))
        for r, (l, b, ref) in zip(top._rows, values):
            r['l'].set(l)
            r['b'].set(b)
            r['use'].set(ref is not None)
            if ref is not None:
                r['rl'].set(ref[0])
                r['rb'].set(ref[1])

    def test_number_and_values_chosen_by_user(self):
        before = len(self.app.chromatograms)
        top = self.app.open_multi_extraction()
        self.fill_in(top, [('254', '4', None), ('280', '4', ('360', '100')), ('310,5', '8', None), ('350', '4', ('400', '20'))])
        with mock.patch.object(hplc.messagebox, 'showerror') as err:
            top._extract()
        self.assertFalse(err.called, err.call_args)
        names = list(self.app.chromatograms)
        self.assertEqual(len(names) - before, 4)
        self.assertEqual(names[before:], ['009-0201 254,4 no ref', '009-0201 280,4 ref 360,100',
                                        '009-0201 310.5,8 no ref', '009-0201 350,4 ref 400,20'])
        # the 280/4 trace with reference 360/100 coincides with the recorded channel dad1A.ch
        t, y, _ = self.app.read_ch(os.path.join(SAMPLE, 'dad1A.ch'))
        mine = self.app.chromatograms['009-0201 280,4 ref 360,100']['df']['mAU'].to_numpy()
        self.assertLess(np.std(mine - y[7:7 + len(mine)]), 0.3)
        self.assertTrue(self.app._dirty)

    def test_error_adds_nothing(self):
        before = list(self.app.chromatograms)
        top = self.app.open_multi_extraction()
        self.fill_in(top, [('254', '4', None), ('', '4', None)])      # second row empty
        with mock.patch.object(hplc.messagebox, 'showerror') as err:
            top._extract()
        self.assertTrue(err.called)
        self.assertIn("Row 2", err.call_args[0][1])
        self.assertEqual(list(self.app.chromatograms), before)
        self.assertTrue(top.winfo_exists())                         # the dialog stays open so it can be corrected

    def test_reopens_with_last_values(self):
        top = self.app.open_multi_extraction()
        self.fill_in(top, [('254', '4', None), ('280', '6', ('360', '100'))])
        with mock.patch.object(hplc.messagebox, 'showerror'):
            top._extract()
        top2 = self.app.open_multi_extraction()
        self.assertEqual(len(top2._rows), 2)
        self.assertEqual([r['l'].get() for r in top2._rows], ['254', '280'])
        self.assertTrue(top2._rows[1]['use'].get())
        self.assertEqual(top2._rows[1]['b'].get(), '6')

    def test_without_spectra(self):
        self.app.spectra.clear()
        with mock.patch.object(hplc.messagebox, 'showinfo') as info:
            self.assertIsNone(self.app.open_multi_extraction())
        self.assertTrue(info.called)


if __name__ == "__main__":
    unittest.main()
