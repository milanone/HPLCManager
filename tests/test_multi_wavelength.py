"""Test dell'estrazione di piu' lunghezze d'onda in un colpo (dialogo 'Several wavelengths').

`estrai_serie` si prova senza finestra su spettri sintetici; il dialogo vero (numero di righe,
lunghezze d'onda, riferimento per riga, errori) gira con una finestra Tk nascosta sul campione
'campioni di esempio/009-0201.D' (gitignorato), se presente.

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


def app_sintetica():
    app = hplc.HPLCManager.__new__(hplc.HPLCManager)
    wl = 200 + 2 * np.arange(201, dtype=float)
    t = np.arange(10) / 60.0
    S = np.outer(np.arange(10) + 1.0, np.ones(201))      # spettro piatto, scala col tempo
    S[:, wl > 400] = 0.0                                    # sopra 400 nm vale 0
    app.cromatogrammi = {}
    app.spettri = {'D': {'t': t, 'wl': wl, 'S': S, 'info': {'Sample': 'X', 'Spectra': 's', 'Version': '31'}}}
    return app


class TestEstraiSerie(unittest.TestCase):
    def test_piu_tracce(self):
        app = app_sintetica()
        nomi = app.estrai_serie('D', [(250, 4, None, None), (300, 8, 500, 20), (350, 4, 500, 20)])
        self.assertEqual(nomi, ['D 250,4 no ref', 'D 300,8 ref 500,20', 'D 350,4 ref 500,20'])
        self.assertEqual(list(app.cromatogrammi), nomi)
        # spettro piatto = (i+1) fino a 400 nm, 0 oltre: con riferimento a 500 nm il segnale non cambia
        for n in nomi:
            np.testing.assert_allclose(app.cromatogrammi[n]['df']['mAU'], np.arange(10) + 1.0)
        self.assertEqual(app.cromatogrammi[nomi[0]]['info']['Sample'], 'X')
        self.assertNotIn('Spectra', app.cromatogrammi[nomi[0]]['info'])

    def test_tutto_o_niente(self):
        app = app_sintetica()
        with self.assertRaisesRegex(ValueError, "Riga 2"):
            app.estrai_serie('D', [(250, 4, None, None), (900, 4, None, None)])   # fuori dallo spettro
        self.assertEqual(app.cromatogrammi, {})
        with self.assertRaisesRegex(ValueError, "Riga 1"):
            app.estrai_serie('D', [(250, 0, None, None)])                         # banda nulla

    def test_stessa_lunghezza_due_volte(self):
        app = app_sintetica()
        nomi = app.estrai_serie('D', [(250, 4, None, None), (250, 4, None, None)])
        self.assertEqual(len(set(nomi)), 2)       # i nomi restano distinti


@unittest.skipUnless(os.path.isdir(CAMPIONE), "campione 009-0201.D non presente")
class TestDialogo(unittest.TestCase):
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

    def tearDown(self):
        for w in self.root.winfo_children():
            w.destroy()

    def riempi(self, top, valori):
        top._var_n.set(len(valori))
        top._imposta_n()
        self.assertEqual(len(top._righe), len(valori))
        for r, (l, b, rif) in zip(top._righe, valori):
            r['l'].set(l)
            r['b'].set(b)
            r['usa'].set(rif is not None)
            if rif is not None:
                r['rl'].set(rif[0])
                r['rb'].set(rif[1])

    def test_numero_e_valori_scelti_dall_utente(self):
        prima = len(self.app.cromatogrammi)
        top = self.app.apri_estrazione_multipla()
        self.riempi(top, [('254', '4', None), ('280', '4', ('360', '100')), ('310,5', '8', None), ('350', '4', ('400', '20'))])
        with mock.patch.object(hplc.messagebox, 'showerror') as err:
            top._estrai()
        self.assertFalse(err.called, err.call_args)
        nomi = list(self.app.cromatogrammi)
        self.assertEqual(len(nomi) - prima, 4)
        self.assertEqual(nomi[prima:], ['009-0201 254,4 no ref', '009-0201 280,4 ref 360,100',
                                        '009-0201 310.5,8 no ref', '009-0201 350,4 ref 400,20'])
        # la traccia 280/4 con riferimento 360/100 coincide con il canale registrato dad1A.ch
        t, y, _ = self.app.leggi_ch(os.path.join(CAMPIONE, 'dad1A.ch'))
        mia = self.app.cromatogrammi['009-0201 280,4 ref 360,100']['df']['mAU'].to_numpy()
        self.assertLess(np.std(mia - y[7:7 + len(mia)]), 0.3)
        self.assertTrue(self.app._dirty)

    def test_errore_non_aggiunge_nulla(self):
        prima = list(self.app.cromatogrammi)
        top = self.app.apri_estrazione_multipla()
        self.riempi(top, [('254', '4', None), ('', '4', None)])      # seconda riga vuota
        with mock.patch.object(hplc.messagebox, 'showerror') as err:
            top._estrai()
        self.assertTrue(err.called)
        self.assertIn("Riga 2", err.call_args[0][1])
        self.assertEqual(list(self.app.cromatogrammi), prima)
        self.assertTrue(top.winfo_exists())                         # il dialogo resta aperto per correggere

    def test_riapre_con_gli_ultimi_valori(self):
        top = self.app.apri_estrazione_multipla()
        self.riempi(top, [('254', '4', None), ('280', '6', ('360', '100'))])
        with mock.patch.object(hplc.messagebox, 'showerror'):
            top._estrai()
        top2 = self.app.apri_estrazione_multipla()
        self.assertEqual(len(top2._righe), 2)
        self.assertEqual([r['l'].get() for r in top2._righe], ['254', '280'])
        self.assertTrue(top2._righe[1]['usa'].get())
        self.assertEqual(top2._righe[1]['b'].get(), '6')

    def test_senza_spettri(self):
        self.app.spettri.clear()
        with mock.patch.object(hplc.messagebox, 'showinfo') as info:
            self.assertIsNone(self.app.apri_estrazione_multipla())
        self.assertTrue(info.called)


if __name__ == "__main__":
    unittest.main()
