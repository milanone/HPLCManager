"""Test di esportazione e sessione di HPLCManager.

Le tabelle (picchi, spettri) si provano senza finestra. Il salvataggio vero (CSV, sessione, figura)
gira sull'applicazione con una finestra Tk nascosta, sostituendo solo le finestre di dialogo
(file e messaggi); servono quindi uno schermo e il campione 'campioni di esempio/009-0201.D'
(gitignorato): senza, quei test vengono saltati.

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

CAMPIONE = os.path.join(HERE, "..", "campioni di esempio", "009-0201.D")


def traccia(picchi=(), picchi_cs=()):
    df = pd.DataFrame({'mAU': np.arange(5.0)}, index=pd.Index(np.arange(5.0), name='Time (min)'))
    return {'df': df, 'info': {}, 'picchi': list(picchi), 'picchi_cs': list(picchi_cs), 'nascosto': False}


class TestTabelle(unittest.TestCase):
    def setUp(self):
        self.app = hplc.HPLCManager.__new__(hplc.HPLCManager)

    def test_tabella_picchi(self):
        mio = {'n': 1, 'rt': 4.35, 'altezza': 100.0, 'area': 2000.0, 'pct': 80.0, 'inizio': 4.2, 'fine': 4.5}
        cs = {'n': 9, 'rt': 4.354, 'tipo': 'VV', 'larghezza': 0.18, 'area': 23183.0, 'altezza': 2003.0,
              'area_pct': 41.9}
        self.app.cromatogrammi = {'A': traccia([mio], [cs]), 'B': traccia()}
        df = self.app.tabella_picchi()
        self.assertEqual(list(df['Source']), ['HPLCManager', 'ChemStation'])
        self.assertEqual(list(df['Trace']), ['A', 'A'])
        self.assertEqual(df.loc[0, 'Area (mAU*s)'], 2000.0)
        self.assertEqual(df.loc[1, 'Type'], 'VV')
        self.assertEqual(df.loc[1, 'Width (min)'], 0.18)
        self.assertTrue(np.isnan(df.loc[0, 'Width (min)']))

    def test_tabella_picchi_vuota(self):
        self.app.cromatogrammi = {'A': traccia()}
        self.assertTrue(self.app.tabella_picchi().empty)

    def test_tabella_spettri(self):
        wl = np.array([200.0, 202.0, 204.0])
        self.app._spettri_scelti = [
            {'dataset': 'D1', 't': 4.357, 'wl': wl, 'y': np.array([1.0, 2.0, 3.0])},
            {'dataset': 'D1', 't': 8.8, 'wl': wl, 'y': np.array([4.0, 5.0, 6.0])}]
        df = self.app.tabella_spettri()
        self.assertEqual(df.index.name, 'Wavelength (nm)')
        self.assertEqual(list(df.columns), ['D1 t=4.357 min', 'D1 t=8.800 min'])
        np.testing.assert_allclose(df['D1 t=8.800 min'], [4.0, 5.0, 6.0])

    def test_tabella_spettri_vuota(self):
        self.app._spettri_scelti = []
        self.assertTrue(self.app.tabella_spettri().empty)


@unittest.skipUnless(os.path.isdir(CAMPIONE), "campione 009-0201.D non presente")
class TestConApplicazione(unittest.TestCase):
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
        self.app.var_includi_ch.set(True)                 # anche i canali con i picchi ChemStation
        self.app.processa_file(CAMPIONE)
        self.app.aggiorna_vista()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def tearDown(self):
        for w in self.root.winfo_children():
            w.destroy()

    def path(self, nome):
        return os.path.join(self.tmp.name, nome)

    def salva_con(self, funzione, nome):
        """Esegue `funzione` con la finestra 'salva con nome' e i messaggi simulati."""
        p = self.path(nome)
        with mock.patch.object(hplc.filedialog, 'asksaveasfilename', return_value=p), \
                mock.patch.object(hplc.messagebox, 'showinfo') as info, \
                mock.patch.object(hplc.messagebox, 'showerror') as err, \
                mock.patch.object(hplc.messagebox, 'showwarning') as warn:
            funzione()
        self.assertFalse(err.called, err.call_args)
        self.assertFalse(warn.called, warn.call_args)
        return p, info

    def test_export_picchi(self):
        n = [k for k in self.app.cromatogrammi if k.endswith('dad1A')][0]
        self.app.lista.selection_clear(0, tk.END)       # un click sostituisce la selezione
        self.app.lista.selection_set(list(self.app.cromatogrammi).index(n))
        self.app.apri_picchi()
        self.app.var_prom.set('5')
        self.app._pk_applica()
        p, _ = self.salva_con(self.app.esporta_picchi, 'picchi.csv')
        df = pd.read_csv(p, sep=';', encoding='latin-1')
        self.assertEqual(set(df['Source']), {'HPLCManager', 'ChemStation'})
        cs = df[(df['Source'] == 'ChemStation') & (df['Trace'] == n)]
        self.assertEqual(len(cs), 60)                     # picchi del REPORT01.CSV
        self.assertAlmostEqual(cs['Area (mAU*s)'].max(), 23183.0, delta=5)

    def test_export_picchi_senza_picchi(self):
        self.app.cromatogrammi = {'X': traccia()}
        with mock.patch.object(hplc.messagebox, 'showinfo') as info:
            self.app.esporta_picchi()
        self.assertTrue(info.called)

    def test_export_spettri(self):
        self.app._sp_dataset = lambda: self.app.spettri[next(iter(self.app.spettri))]
        self.app.var_sp_set = tk.StringVar(value=next(iter(self.app.spettri)))
        self.app._sp_spettro_a(4.354)
        self.app._sp_spettro_a(8.8)
        p, _ = self.salva_con(self.app.esporta_spettri, 'spettri.csv')
        df = pd.read_csv(p, sep=';', encoding='latin-1', index_col=0)
        self.assertEqual(df.shape, (201, 2))
        self.assertEqual(df.index[0], 200.0)
        # lo spettro a 4.35 min ha il massimo a ~272 nm
        self.assertAlmostEqual(df.iloc[:, 0].idxmax(), 272, delta=6)

    def test_export_dad_completo(self):
        p, _ = self.salva_con(self.app.esporta_dad_completo, 'dad.csv')
        df = pd.read_csv(p, sep=';', encoding='latin-1', index_col=0, nrows=3)
        self.assertEqual(df.shape, (3, 201))
        self.assertEqual(list(df.columns[:2]), ['200', '202'])

    def test_sessione_andata_e_ritorno(self):
        self.app.var_sp_l.set('254')
        self.app._sp_estrai(True)
        prima = list(self.app.cromatogrammi)
        picco = self.app.cromatogrammi[prima[0]]['picchi_cs'][:1]
        p, _ = self.salva_con(self.app.salva_sessione, 's.hplcsession')
        with gzip.open(p, 'rb') as f:
            self.assertEqual(pickle.load(f)['version'], 1)
        self.app.clear_all()
        self.assertEqual(self.app.cromatogrammi, {})
        with mock.patch.object(hplc.filedialog, 'askopenfilename', return_value=p):
            self.app.apri_sessione()
        self.assertEqual(list(self.app.cromatogrammi), prima)
        self.assertEqual(self.app.cromatogrammi[prima[0]]['picchi_cs'][:1], picco)
        self.assertIn(next(iter(self.app.spettri)), self.app.spettri)
        self.assertEqual(self.app.spettri[next(iter(self.app.spettri))]['S'].shape, (13493, 201))
        self.assertTrue(self.app.var_includi_ch.get())
        self.assertFalse(self.app._dirty)

    def test_apri_sessione_non_valida(self):
        p = self.path('x.hplcsession')
        with gzip.open(p, 'wb') as f:
            pickle.dump({'altro': 1}, f)
        with mock.patch.object(hplc.filedialog, 'askopenfilename', return_value=p), \
                mock.patch.object(hplc.messagebox, 'showerror') as err:
            self.app.apri_sessione()
        self.assertTrue(err.called)
        self.assertGreater(len(self.app.cromatogrammi), 0)   # lo stato corrente non viene toccato

    def test_figura(self):
        self.app._sp_dataset = lambda: self.app.spettri[next(iter(self.app.spettri))]
        self.app.var_sp_set = tk.StringVar(value=next(iter(self.app.spettri)))
        self.app._sp_spettro_a(4.354)                      # aggiunge una linea-marcatore
        self.assertTrue(any(a.get_gid() == 'marcatore' for a in self.app.ax.lines))
        p, _ = self.salva_con(self.app.salva_figura_immagine, 'f.png')
        self.assertGreater(os.path.getsize(p), 5000)
        p, _ = self.salva_con(self.app.salva_figura_pickle, 'f.fig.pickle')
        with open(p, 'rb') as f:
            fig = pickle.load(f)
        self.assertEqual(len(fig.axes), 1)
        # i marcatori degli spettri e il cursore non finiscono nella figura salvata...
        self.assertFalse(any(a.get_gid() == 'marcatore' for a in fig.axes[0].lines))
        # ...e restano invece sulla vista viva
        self.assertTrue(any(a.get_gid() == 'marcatore' for a in self.app.ax.lines))

    def test_editor_figura(self):
        if not os.path.isfile(os.path.join(HERE, "..", "..", "PlotStyleKit", "plot_editor.pyw")):
            self.skipTest("PlotStyleKit non presente")
        with mock.patch.object(hplc.messagebox, 'showerror') as err:
            self.app.apri_editor_figura()
        self.assertFalse(err.called, err.call_args)
        self.assertTrue(any(isinstance(w, tk.Toplevel) for w in self.root.winfo_children()))


if __name__ == "__main__":
    unittest.main()
