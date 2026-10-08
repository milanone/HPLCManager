"""Test of the spectra window: ticked list of the picked spectra.

Runs on the application with a hidden Tk window and the sample
'example data/009-0201.D' (gitignored); without it, it is skipped.

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


@unittest.skipUnless(os.path.isdir(SAMPLE), "sample 009-0201.D not present")
class TestSpectraWindow(unittest.TestCase):
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
        for t in (4.354, 8.8, 5.347):
            self.app._sp_spectrum_at(t)
        self.w = self.app._spec_win

    def tearDown(self):
        for w in self.root.winfo_children():
            w.destroy()

    def markers(self):
        return [a for a in self.app.ax.lines if a.get_gid() == 'marker']

    def visible_lines(self):
        return [l for l in self.w.ax.lines if l.get_visible()]

    def test_initial_state(self):
        self.assertEqual(len(self.app._chosen_spectra), 3)
        self.assertEqual(len(self.w.f_rows.winfo_children()), 3)
        self.assertEqual(len(self.markers()), 3)
        self.assertEqual(len(self.visible_lines()), 3)
        # three distinct, stable colors
        self.assertEqual(len({k['color'] for k in self.app._chosen_spectra}), 3)

    def test_untick_hides_line_marker_and_export(self):
        k = self.app._chosen_spectra[1]
        color = k['color']
        k['var'].set(False)
        self.app._sp_toggle(k, False)
        self.assertEqual(len(self.visible_lines()), 2)
        self.assertFalse(k['line'].get_visible())
        self.assertEqual(len(self.markers()), 2)
        self.assertNotIn(color, [m.get_color() for m in self.markers()])
        self.assertEqual(len(self.w.ax.get_legend().get_texts()), 2)
        self.assertEqual(self.app.spectra_table().shape[1], 2)     # the export skips the unticked one
        self.assertEqual(len(self.app._chosen_spectra), 3)           # but it is still in the list
        # it can be ticked again and comes back as it was, with the same color
        self.app._sp_toggle(k, True)
        self.assertEqual(len(self.visible_lines()), 3)
        self.assertEqual(len(self.markers()), 3)
        self.assertEqual(k['color'], color)
        self.assertEqual(self.app.spectra_table().shape[1], 3)

    def test_remove(self):
        k = self.app._chosen_spectra[2]
        self.app._sp_remove(k)
        self.assertEqual(len(self.app._chosen_spectra), 2)
        self.assertEqual(len(self.w.f_rows.winfo_children()), 2)
        self.assertEqual(len(self.w.ax.lines), 2)
        self.assertEqual(len(self.markers()), 2)
        self.assertEqual(self.app.spectra_table().shape[1], 2)

    def test_clear_all_and_new_click(self):
        self.app._sp_clear()
        self.assertEqual(self.app._chosen_spectra, [])
        self.assertEqual(len(self.markers()), 0)
        self.assertEqual(len(self.w.ax.lines), 0)
        self.assertTrue(self.w.winfo_exists())                 # the window stays open
        self.app._sp_spectrum_at(4.354)                          # and one can start again
        self.assertEqual(len(self.app._chosen_spectra), 1)
        self.assertEqual(len(self.markers()), 1)

    def test_window_zoom_is_kept(self):
        self.w.ax.set_xlim(250, 300)
        self.w.ax.set_ylim(0, 500)
        self.app._sp_toggle(self.app._chosen_spectra[0], False)
        self.assertEqual(self.w.ax.get_xlim(), (250, 300))
        self.assertEqual(self.w.ax.get_ylim(), (0, 500))

    def test_chromatogram_zoom_is_kept(self):
        self.app.ax.set_xlim(3, 6)
        self.app.canvas.draw()
        self.app._sp_toggle(self.app._chosen_spectra[0], False)
        self.app._sp_remove(self.app._chosen_spectra[0])
        self.assertEqual(self.app.ax.get_xlim(), (3, 6))

    def test_window_closing(self):
        self.w.destroy()
        self.root.update()
        self.assertEqual(self.app._chosen_spectra, [])
        self.assertIsNone(self.app._spec_win)
        self.assertEqual(len(self.markers()), 0)

    def test_app_clear_all_closes_the_window(self):
        self.app.clear_all()
        self.root.update()
        self.assertEqual(self.app._chosen_spectra, [])
        self.assertIsNone(self.app._spec_win)

    def test_export_only_the_selected(self):
        self.app._sp_toggle(self.app._chosen_spectra[0], False)
        with mock.patch.object(hplc.messagebox, 'showinfo') as info:
            with mock.patch.object(hplc.filedialog, 'asksaveasfilename', return_value=''):
                self.app.export_spectra()
        self.assertFalse(info.called)          # there were selected spectra: no "no spectrum" warning
        for k in self.app._chosen_spectra:
            self.app._sp_toggle(k, False)
        with mock.patch.object(hplc.messagebox, 'showinfo') as info:
            self.app.export_spectra()
        self.assertTrue(info.called)           # all unticked: warning


if __name__ == "__main__":
    unittest.main()
