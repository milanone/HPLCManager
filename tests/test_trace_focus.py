"""Test: the extracted trace selected in the list is the one that 'Update trace' modifies.

Reproduces the usage sequence: extract a trace, add another with 'Add as new', change
wavelength or bandwidth and press Update: the active trace must change, not always the first.
The selection in the list is the real one (<<ListboxSelect>> event).

    py -m unittest discover -s tests -v
"""
import importlib.util
import os
import pickle
import tkinter as tk
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "hplc", os.path.join(HERE, "..", "HPLCManager.pyw"))
hplc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hplc)

SAMPLE = os.path.join(HERE, "..", "example data", "009-0201.D")


@unittest.skipUnless(os.path.isdir(SAMPLE), "sample 009-0201.D not present")
class TestActiveTrace(unittest.TestCase):
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
        self.root.update()

    def tearDown(self):
        for w in self.root.winfo_children():
            w.destroy()

    def names(self):
        return list(self.app.chromatograms)

    def select_row(self, i):
        """Real selection in the list: like a click on row i."""
        self.app.listbox.selection_clear(0, tk.END)
        self.app.listbox.selection_set(i)
        self.app.listbox.event_generate('<<ListboxSelect>>')
        self.root.update()

    def set_params(self, l=None, b=None, ref=None):
        if l is not None:
            self.app.var_sp_l.set(str(l))
        if b is not None:
            self.app.var_sp_b.set(str(b))
        if ref is not None:
            self.app.var_sp_use_ref.set(ref is not False)
            if ref is not False:
                self.app.var_sp_rl.set(str(ref[0]))
                self.app.var_sp_rb.set(str(ref[1]))
            self.app._sp_ref_state()

    def test_user_sequence(self):
        # first trace: 280/4, no reference (the one created on opening)
        self.assertEqual(self.names(), ['009-0201 280,4 no ref'])
        # second: 254/4 with Add as new
        self.set_params(l=254)
        self.app._sp_extract(True)
        self.assertEqual(self.names(), ['009-0201 280,4 no ref', '009-0201 254,4 no ref'])
        # the new one is selected: Update changes its wavelength, the first is not touched
        self.assertEqual(self.app._selected(), ['009-0201 254,4 no ref'])
        before = self.app.chromatograms['009-0201 280,4 no ref']['df']['mAU'].to_numpy().copy()
        self.set_params(l=300)
        self.app._sp_extract(False)
        self.assertEqual(self.names(), ['009-0201 280,4 no ref', '009-0201 300,4 no ref'])
        np.testing.assert_array_equal(self.app.chromatograms['009-0201 280,4 no ref']['df']['mAU'], before)
        # the modified one stays selected: a second Update (bandwidth) still acts on it
        self.assertEqual(self.app._selected(), ['009-0201 300,4 no ref'])
        self.set_params(b=8)
        self.app._sp_extract(False)
        self.assertEqual(self.names(), ['009-0201 280,4 no ref', '009-0201 300,8 no ref'])

    def test_selecting_a_trace_reloads_the_parameters(self):
        self.set_params(l=254, b=8, ref=(400, 20))
        self.app._sp_extract(True)                         # 254/8 with reference 400/20
        self.set_params(l=320, b=4, ref=False)
        self.app._sp_extract(True)                         # 320/4 without reference
        self.assertEqual(len(self.names()), 3)
        self.select_row(1)                                 # click on the second trace
        self.assertEqual((self.app.var_sp_l.get(), self.app.var_sp_b.get()), ('254', '8'))
        self.assertTrue(self.app.var_sp_use_ref.get())
        self.assertEqual((self.app.var_sp_rl.get(), self.app.var_sp_rb.get()), ('400', '20'))
        self.assertEqual(str(self.app.e_rl.cget('state')), 'normal')
        self.select_row(2)                                 # third: without reference
        self.assertEqual(self.app.var_sp_l.get(), '320')
        self.assertFalse(self.app.var_sp_use_ref.get())
        self.assertEqual(str(self.app.e_rl.cget('state')), 'disabled')
        self.select_row(0)                                 # first
        self.assertEqual((self.app.var_sp_l.get(), self.app.var_sp_b.get()), ('280', '4'))

    def test_modifies_the_clicked_trace_not_the_first(self):
        self.set_params(l=254)
        self.app._sp_extract(True)
        self.set_params(l=320)
        self.app._sp_extract(True)                         # three traces: 280, 254, 320
        self.select_row(1)                                 # back to the second (254)
        self.set_params(l=260)
        self.app._sp_extract(False)
        self.assertEqual(self.names(), ['009-0201 280,4 no ref', '009-0201 260,4 no ref',
                                       '009-0201 320,4 no ref'])
        self.assertEqual(self.app._selected(), ['009-0201 260,4 no ref'])

    def test_wavelength_series(self):
        names = self.app.extract_series('009-0201', [(254, 4, None, None), (280, 4, 360, 100)])
        self.assertEqual(self.app.chromatograms[names[1]]['estr'],
                         {'ds': '009-0201', 'l': 280, 'b': 4, 'ref': 360, 'rb': 100})
        self.app.refresh_view()
        self.select_row(self.names().index(names[1]))
        self.assertTrue(self.app.var_sp_use_ref.get())
        self.assertEqual(self.app.var_sp_l.get(), '280')
        # Update on the selected one of the series: it changes and only it
        self.set_params(l=300)
        self.app._sp_extract(False)
        self.assertIn('009-0201 300,4 ref 360,100', self.names())
        self.assertIn('009-0201 254,4 no ref', self.names())
        self.assertIn('009-0201 280,4 no ref', self.names())          # the first, untouched

    def test_non_extracted_trace_falls_back_to_the_last(self):
        # with a recorded channel selected (not extracted from the spectra) Update does not touch it
        self.app.var_include_ch.set(True)
        self.app.clear_all()
        self.app.process_file(SAMPLE)
        self.app.refresh_view()
        names = self.names()
        channel = next(n for n in names if n.endswith('dad1A'))
        before = self.app.chromatograms[channel]['df']['mAU'].to_numpy().copy()
        self.select_row(names.index(channel))
        self.set_params(l=300)
        self.app._sp_extract(False)
        np.testing.assert_array_equal(self.app.chromatograms[channel]['df']['mAU'], before)
        self.assertTrue(any(n.endswith('300,4 no ref') for n in self.names()))

    def test_session_keeps_the_parameters(self):
        self.set_params(l=254, b=8, ref=(400, 20))
        self.app._sp_extract(True)
        data = pickle.loads(pickle.dumps(self.app._session_data()))
        self.app.clear_all()
        self.app._apply_session(data)
        self.select_row(1)
        self.assertEqual((self.app.var_sp_l.get(), self.app.var_sp_b.get()), ('254', '8'))
        self.assertTrue(self.app.var_sp_use_ref.get())


if __name__ == "__main__":
    unittest.main()
