"""Test: la traccia estratta selezionata nella lista e' quella che 'Update trace' modifica.

Riproduce la sequenza d'uso: estrai una traccia, aggiungine un'altra con 'Add as new', cambia
lunghezza d'onda o banda e premi Update: deve cambiare la traccia attiva, non sempre la prima.
La selezione nella lista e' quella vera (evento <<ListboxSelect>>).

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

CAMPIONE = os.path.join(HERE, "..", "campioni di esempio", "009-0201.D")


@unittest.skipUnless(os.path.isdir(CAMPIONE), "campione 009-0201.D non presente")
class TestTracciaAttiva(unittest.TestCase):
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
        self.root.update()

    def tearDown(self):
        for w in self.root.winfo_children():
            w.destroy()

    def nomi(self):
        return list(self.app.cromatogrammi)

    def seleziona(self, i):
        """Selezione vera nella lista: come un click sulla riga i."""
        self.app.lista.selection_clear(0, tk.END)
        self.app.lista.selection_set(i)
        self.app.lista.event_generate('<<ListboxSelect>>')
        self.root.update()

    def imposta(self, l=None, b=None, rif=None):
        if l is not None:
            self.app.var_sp_l.set(str(l))
        if b is not None:
            self.app.var_sp_b.set(str(b))
        if rif is not None:
            self.app.var_sp_uso_rif.set(rif is not False)
            if rif is not False:
                self.app.var_sp_rl.set(str(rif[0]))
                self.app.var_sp_rb.set(str(rif[1]))
            self.app._sp_stato_rif()

    def test_sequenza_dell_utente(self):
        # prima traccia: 280/4, senza riferimento (quella creata all'apertura)
        self.assertEqual(self.nomi(), ['009-0201 280,4 no ref'])
        # seconda: 254/4 con Add as new
        self.imposta(l=254)
        self.app._sp_estrai(True)
        self.assertEqual(self.nomi(), ['009-0201 280,4 no ref', '009-0201 254,4 no ref'])
        # la nuova e' la selezionata: Update ne cambia la lunghezza d'onda, la prima non si tocca
        self.assertEqual(self.app._selezionati(), ['009-0201 254,4 no ref'])
        prima = self.app.cromatogrammi['009-0201 280,4 no ref']['df']['mAU'].to_numpy().copy()
        self.imposta(l=300)
        self.app._sp_estrai(False)
        self.assertEqual(self.nomi(), ['009-0201 280,4 no ref', '009-0201 300,4 no ref'])
        np.testing.assert_array_equal(self.app.cromatogrammi['009-0201 280,4 no ref']['df']['mAU'], prima)
        # resta selezionata quella modificata: un secondo Update (banda) agisce ancora su di lei
        self.assertEqual(self.app._selezionati(), ['009-0201 300,4 no ref'])
        self.imposta(b=8)
        self.app._sp_estrai(False)
        self.assertEqual(self.nomi(), ['009-0201 280,4 no ref', '009-0201 300,8 no ref'])

    def test_selezionare_una_traccia_ricarica_i_parametri(self):
        self.imposta(l=254, b=8, rif=(400, 20))
        self.app._sp_estrai(True)                         # 254/8 con riferimento 400/20
        self.imposta(l=320, b=4, rif=False)
        self.app._sp_estrai(True)                         # 320/4 senza riferimento
        self.assertEqual(len(self.nomi()), 3)
        self.seleziona(1)                                 # click sulla seconda traccia
        self.assertEqual((self.app.var_sp_l.get(), self.app.var_sp_b.get()), ('254', '8'))
        self.assertTrue(self.app.var_sp_uso_rif.get())
        self.assertEqual((self.app.var_sp_rl.get(), self.app.var_sp_rb.get()), ('400', '20'))
        self.assertEqual(str(self.app.e_rl.cget('state')), 'normal')
        self.seleziona(2)                                 # terza: senza riferimento
        self.assertEqual(self.app.var_sp_l.get(), '320')
        self.assertFalse(self.app.var_sp_uso_rif.get())
        self.assertEqual(str(self.app.e_rl.cget('state')), 'disabled')
        self.seleziona(0)                                 # prima
        self.assertEqual((self.app.var_sp_l.get(), self.app.var_sp_b.get()), ('280', '4'))

    def test_modifica_la_traccia_cliccata_non_la_prima(self):
        self.imposta(l=254)
        self.app._sp_estrai(True)
        self.imposta(l=320)
        self.app._sp_estrai(True)                         # tre tracce: 280, 254, 320
        self.seleziona(1)                                 # torno sulla seconda (254)
        self.imposta(l=260)
        self.app._sp_estrai(False)
        self.assertEqual(self.nomi(), ['009-0201 280,4 no ref', '009-0201 260,4 no ref',
                                       '009-0201 320,4 no ref'])
        self.assertEqual(self.app._selezionati(), ['009-0201 260,4 no ref'])

    def test_serie_di_lunghezze_d_onda(self):
        nomi = self.app.estrai_serie('009-0201', [(254, 4, None, None), (280, 4, 360, 100)])
        self.assertEqual(self.app.cromatogrammi[nomi[1]]['estr'],
                         {'ds': '009-0201', 'l': 280, 'b': 4, 'rif': 360, 'rb': 100})
        self.app.aggiorna_vista()
        self.seleziona(self.nomi().index(nomi[1]))
        self.assertTrue(self.app.var_sp_uso_rif.get())
        self.assertEqual(self.app.var_sp_l.get(), '280')
        # Update sulla selezionata della serie: cambia lei e solo lei
        self.imposta(l=300)
        self.app._sp_estrai(False)
        self.assertIn('009-0201 300,4 ref 360,100', self.nomi())
        self.assertIn('009-0201 254,4 no ref', self.nomi())
        self.assertIn('009-0201 280,4 no ref', self.nomi())          # la prima, intatta

    def test_traccia_non_estratta_cade_sull_ultima(self):
        # con un canale registrato selezionato (non estratto dagli spettri) Update non lo tocca
        self.app.var_includi_ch.set(True)
        self.app.clear_all()
        self.app.processa_file(CAMPIONE)
        self.app.aggiorna_vista()
        nomi = self.nomi()
        canale = next(n for n in nomi if n.endswith('dad1A'))
        prima = self.app.cromatogrammi[canale]['df']['mAU'].to_numpy().copy()
        self.seleziona(nomi.index(canale))
        self.imposta(l=300)
        self.app._sp_estrai(False)
        np.testing.assert_array_equal(self.app.cromatogrammi[canale]['df']['mAU'], prima)
        self.assertTrue(any(n.endswith('300,4 no ref') for n in self.nomi()))

    def test_sessione_conserva_i_parametri(self):
        self.imposta(l=254, b=8, rif=(400, 20))
        self.app._sp_estrai(True)
        dati = pickle.loads(pickle.dumps(self.app._dati_sessione()))
        self.app.clear_all()
        self.app._applica_sessione(dati)
        self.seleziona(1)
        self.assertEqual((self.app.var_sp_l.get(), self.app.var_sp_b.get()), ('254', '8'))
        self.assertTrue(self.app.var_sp_uso_rif.get())


if __name__ == "__main__":
    unittest.main()
