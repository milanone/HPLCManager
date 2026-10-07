"""Test della finestra degli spettri: elenco spuntabile degli spettri scelti.

Gira sull'applicazione con una finestra Tk nascosta e il campione
'campioni di esempio/009-0201.D' (gitignorato); senza, viene saltato.

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

CAMPIONE = os.path.join(HERE, "..", "campioni di esempio", "009-0201.D")


@unittest.skipUnless(os.path.isdir(CAMPIONE), "campione 009-0201.D non presente")
class TestFinestraSpettri(unittest.TestCase):
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
        self.app.processa_file(CAMPIONE)
        self.app.aggiorna_vista()
        for t in (4.354, 8.8, 5.347):
            self.app._sp_spettro_a(t)
        self.w = self.app._spec_win

    def tearDown(self):
        for w in self.root.winfo_children():
            w.destroy()

    def marcatori(self):
        return [a for a in self.app.ax.lines if a.get_gid() == 'marcatore']

    def linee_visibili(self):
        return [l for l in self.w.ax.lines if l.get_visible()]

    def test_stato_iniziale(self):
        self.assertEqual(len(self.app._spettri_scelti), 3)
        self.assertEqual(len(self.w.f_righe.winfo_children()), 3)
        self.assertEqual(len(self.marcatori()), 3)
        self.assertEqual(len(self.linee_visibili()), 3)
        # tre colori distinti, stabili
        self.assertEqual(len({k['colore'] for k in self.app._spettri_scelti}), 3)

    def test_deseleziona_nasconde_linea_marcatore_e_esportazione(self):
        k = self.app._spettri_scelti[1]
        colore = k['colore']
        k['var'].set(False)
        self.app._sp_spunta(k, False)
        self.assertEqual(len(self.linee_visibili()), 2)
        self.assertFalse(k['linea'].get_visible())
        self.assertEqual(len(self.marcatori()), 2)
        self.assertNotIn(colore, [m.get_color() for m in self.marcatori()])
        self.assertEqual(len(self.w.ax.get_legend().get_texts()), 2)
        self.assertEqual(self.app.tabella_spettri().shape[1], 2)     # l'export salta quello deselezionato
        self.assertEqual(len(self.app._spettri_scelti), 3)           # ma e' ancora in elenco
        # lo si puo' riselezionare e ritorna com'era, con lo stesso colore
        self.app._sp_spunta(k, True)
        self.assertEqual(len(self.linee_visibili()), 3)
        self.assertEqual(len(self.marcatori()), 3)
        self.assertEqual(k['colore'], colore)
        self.assertEqual(self.app.tabella_spettri().shape[1], 3)

    def test_rimuovi(self):
        k = self.app._spettri_scelti[2]
        self.app._sp_rimuovi(k)
        self.assertEqual(len(self.app._spettri_scelti), 2)
        self.assertEqual(len(self.w.f_righe.winfo_children()), 2)
        self.assertEqual(len(self.w.ax.lines), 2)
        self.assertEqual(len(self.marcatori()), 2)
        self.assertEqual(self.app.tabella_spettri().shape[1], 2)

    def test_clear_all_e_nuovo_click(self):
        self.app._sp_svuota()
        self.assertEqual(self.app._spettri_scelti, [])
        self.assertEqual(len(self.marcatori()), 0)
        self.assertEqual(len(self.w.ax.lines), 0)
        self.assertTrue(self.w.winfo_exists())                 # la finestra resta aperta
        self.app._sp_spettro_a(4.354)                          # e si puo' ricominciare
        self.assertEqual(len(self.app._spettri_scelti), 1)
        self.assertEqual(len(self.marcatori()), 1)

    def test_zoom_della_finestra_si_conserva(self):
        self.w.ax.set_xlim(250, 300)
        self.w.ax.set_ylim(0, 500)
        self.app._sp_spunta(self.app._spettri_scelti[0], False)
        self.assertEqual(self.w.ax.get_xlim(), (250, 300))
        self.assertEqual(self.w.ax.get_ylim(), (0, 500))

    def test_zoom_del_cromatogramma_si_conserva(self):
        self.app.ax.set_xlim(3, 6)
        self.app.canvas.draw()
        self.app._sp_spunta(self.app._spettri_scelti[0], False)
        self.app._sp_rimuovi(self.app._spettri_scelti[0])
        self.assertEqual(self.app.ax.get_xlim(), (3, 6))

    def test_chiusura_finestra(self):
        self.w.destroy()
        self.root.update()
        self.assertEqual(self.app._spettri_scelti, [])
        self.assertIsNone(self.app._spec_win)
        self.assertEqual(len(self.marcatori()), 0)

    def test_clear_all_dell_app_chiude_la_finestra(self):
        self.app.clear_all()
        self.root.update()
        self.assertEqual(self.app._spettri_scelti, [])
        self.assertIsNone(self.app._spec_win)

    def test_esporta_solo_i_selezionati(self):
        self.app._sp_spunta(self.app._spettri_scelti[0], False)
        with mock.patch.object(hplc.messagebox, 'showinfo') as info:
            with mock.patch.object(hplc.filedialog, 'asksaveasfilename', return_value=''):
                self.app.esporta_spettri()
        self.assertFalse(info.called)          # c'erano spettri selezionati: nessun avviso "nessuno spettro"
        for k in self.app._spettri_scelti:
            self.app._sp_spunta(k, False)
        with mock.patch.object(hplc.messagebox, 'showinfo') as info:
            self.app.esporta_spettri()
        self.assertTrue(info.called)           # tutti deselezionati: avviso


if __name__ == "__main__":
    unittest.main()
