"""How HPLCManager works with and without PlotStyleKit (the shared Origin style and figure editor).

The program carries the bootstrap of PlotStyleKit's integration module: when it cannot find PlotStyleKit it
asks for the folder, and without it every figure action has a fallback. The dialogs are replaced by answers;
the application runs on a hidden Tk window with the sample 'example data/009-0201.D' (gitignored: without it
these tests are skipped).

    py -m unittest discover -s tests -v
"""
import importlib.util
import json
import os
import pickle
import tempfile
import tkinter as tk
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("hplc_psk", os.path.join(HERE, "..", "HPLCManager.pyw"))
hplc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hplc)

SAMPLE = os.path.join(HERE, "..", "example data", "009-0201.D")
KIT = os.path.normpath(os.path.join(HERE, "..", "..", "PlotStyleKit"))


@unittest.skipUnless(os.path.isdir(SAMPLE) and os.path.isfile(os.path.join(KIT, "plotstylekit.py")),
                     "needs the sample data and the PlotStyleKit folder next to this project")
class TestWithAndWithoutPlotStyleKit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.root = tk.Tk()
        except tk.TclError as e:
            raise unittest.SkipTest("Tk not available: %s" % e)
        cls.root.withdraw()
        cls.app = hplc.HPLCManager(cls.root)
        cls.app.var_include_ch.set(True)
        cls.app.process_file(SAMPLE)
        cls.app.refresh_view()

    @classmethod
    def tearDownClass(cls):
        cls.root.destroy()

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.saved = (hplc.psk, hplc.PSK_PROBLEM, hplc.PSK_CONFIG)
        hplc.PSK_CONFIG = os.path.join(self.tmp.name, "psk.json")          # never touch the real ~/.plotstylekit.json
        self.addCleanup(self._restore)

    def _restore(self):
        hplc.psk, hplc.PSK_PROBLEM, hplc.PSK_CONFIG = self.saved
        for w in self.root.winfo_children():
            if isinstance(w, tk.Toplevel):             # the editor windows, not the application itself
                w.destroy()

    def path(self, name):
        return os.path.join(self.tmp.name, name)

    def test_the_kit_is_in_use_and_styles_the_plot(self):
        self.assertIsNotNone(hplc.psk)
        self.app.refresh_view()
        self.assertEqual(self.app.ax.xaxis.get_ticks_position(), 'bottom')        # Origin: ticks on bottom/left only
        self.assertFalse(any(g.get_visible() for g in self.app.ax.xaxis.get_gridlines()))

    def test_saved_figure_has_the_publication_size_and_keeps_the_legend_corner(self):
        with mock.patch.object(hplc.filedialog, 'asksaveasfilename', return_value=self.path('f.fig.pickle')), \
                mock.patch.object(hplc.messagebox, 'showinfo'):
            self.app.save_figure_pickle()
        with open(self.path('f.fig.pickle'), 'rb') as f:
            fig = pickle.load(f)
        self.assertAlmostEqual(fig.get_size_inches()[0] * 2.54, 9.0, places=3)
        self.assertEqual(fig.axes[0].get_legend()._loc, 1)                         # 'upper right', as fixed on screen

    def test_without_the_kit_the_figure_is_saved_as_it_is(self):
        hplc.psk = None
        live = tuple(self.app.fig.get_size_inches())
        with mock.patch.object(hplc.filedialog, 'asksaveasfilename', return_value=self.path('plain.fig.pickle')), \
                mock.patch.object(hplc.messagebox, 'showinfo') as info, \
                mock.patch.object(hplc.messagebox, 'showerror') as err:
            self.app.save_figure_pickle()
        self.assertFalse(err.called, err.call_args)
        self.assertIn('not restyled', info.call_args[0][1])
        with open(self.path('plain.fig.pickle'), 'rb') as f:
            fig = pickle.load(f)
        self.assertEqual(tuple(fig.get_size_inches()), live)
        self.assertFalse(any(a.get_gid() == 'marker' for a in fig.axes[0].lines))

    def test_without_the_kit_the_image_is_written_too(self):
        hplc.psk = None
        with mock.patch.object(hplc.filedialog, 'asksaveasfilename', return_value=self.path('plain.png')), \
                mock.patch.object(hplc.messagebox, 'showinfo'), mock.patch.object(hplc.messagebox, 'showerror') as err:
            self.app.save_figure_image()
        self.assertFalse(err.called, err.call_args)
        self.assertGreater(os.path.getsize(self.path('plain.png')), 5000)

    def test_without_the_kit_edit_figure_offers_to_locate_it(self):
        hplc.psk = None
        with mock.patch.object(hplc.messagebox, 'askyesno', return_value=False) as ask:
            self.app.open_figure_editor()
        self.assertTrue(ask.called)
        self.assertIsNone(hplc.psk)
        self.assertFalse(any(isinstance(w, tk.Toplevel) for w in self.root.winfo_children()))

    def test_cancelling_the_warning_is_remembered(self):
        hplc.psk = None
        with mock.patch.object(hplc.messagebox, 'askokcancel', return_value=False) as ask:
            self.app._offer_plotstylekit("PlotStyleKit was not found")
        self.assertEqual(ask.call_count, 1)
        self.assertIsNone(hplc.psk)
        with open(hplc.PSK_CONFIG, encoding='utf-8') as f:
            self.assertTrue(json.load(f)['declined'])

    def test_choosing_the_folder_starts_the_kit_and_remembers_it(self):
        hplc.psk = None
        with mock.patch.object(hplc.messagebox, 'askokcancel', return_value=True), \
                mock.patch.object(hplc.filedialog, 'askdirectory', return_value=KIT), \
                mock.patch.object(hplc.messagebox, 'showinfo'):
            self.app._offer_plotstylekit("PlotStyleKit was not found")
        self.assertIsNotNone(hplc.psk)
        with open(hplc.PSK_CONFIG, encoding='utf-8') as f:
            cfg = json.load(f)
        self.assertEqual(os.path.normpath(cfg['path']), KIT)
        self.assertFalse(cfg['declined'])

    def test_a_wrong_folder_is_refused_with_the_reason(self):
        hplc.psk = None
        empty = self.path('empty')
        os.makedirs(empty)
        texts = []
        answers = iter([True, False])
        with mock.patch.object(hplc.messagebox, 'askokcancel',
                               side_effect=lambda *a, **k: (texts.append(a[1]) or next(answers))), \
                mock.patch.object(hplc.filedialog, 'askdirectory', return_value=empty):
            self.app._offer_plotstylekit("PlotStyleKit was not found")
        self.assertEqual(len(texts), 2)
        self.assertIn("is not in", texts[1])
        self.assertIsNone(hplc.psk)


if __name__ == "__main__":
    unittest.main()
