"""Test dei lettori dei registri ChemStation: LCDIAG.REG (curve dello strumento) e ACQRES.REG
(moduli e colonna), con la finestra delle curve.

Le curve si provano su un file sintetico (stessa struttura) e sul campione reale
'campioni di esempio/009-0201.D' (gitignorato), confrontando con i valori scritti da ChemStation in
Report00.CSV / RUN.LOG; senza campione i test relativi vengono saltati.

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

CAMPIONE = os.path.join(HERE, "..", "campioni di esempio", "009-0201.D")


def scrivi_diag(path, segnali):
    """Scrive un LCDIAG.REG sintetico: per ogni segnale (titolo, unita', valori, dt_min, fattore)
    prefisso con N e dt, 'min\\0unita'\\0', array uint32 e testata col titolo."""
    b = bytearray(b"\x02" + b"32\x00" + b"\x0dREGISTER FILE" + bytes(40))
    for titolo, unita, valori, dt, fattore in segnali:
        pre = bytearray(400)                         # prefisso di 400 byte: N e dt a offset fissi
        base = len(b) + 400                          # posizione di 'min' nel file
        b += pre
        b[base - 168:base - 164] = struct.pack("<I", len(valori))
        b[base - 116:base - 108] = struct.pack("<d", dt)
        b += b"min\x00" + unita.encode("latin-1") + b"\x00"
        b += np.round(np.asarray(valori) * fattore).astype("<u4").tobytes()
        b += b"ObjClass" + bytes(30) + b"Title" + bytes(20) + titolo.encode() + b"\x00arial\x00" + bytes(60)
    with open(path, "wb") as f:
        f.write(b)


class TestDiagnosticaSintetica(unittest.TestCase):
    def setUp(self):
        self.app = hplc.HPLCManager.__new__(hplc.HPLCManager)
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_scale_e_tempi(self):
        p = os.path.join(self.tmp.name, "LCDIAG.REG")
        scrivi_diag(p, [("PMP1, Pressure", "bar", [94.42, 95.0, 100.55], 0.005, 100),
                        ("PMP1, Flow", "ml/min", [1.0, 1.0, 0.95], 0.005, 1000),
                        ("PMP1, Solvent B", "%", [0.0, 12.5, 99.9], 0.005, 10),
                        ("THM1, Temperature (Left)", "\xb0C", [29.91, 30.0], 1 / 60, 100)])
        seg = self.app.leggi_diagnostica(p)
        self.assertEqual(list(seg), ["PMP1, Pressure", "PMP1, Flow", "PMP1, Solvent B", "THM1, Temperature (Left)"])
        np.testing.assert_allclose(seg["PMP1, Pressure"]["y"], [94.42, 95.0, 100.55])
        np.testing.assert_allclose(seg["PMP1, Flow"]["y"], [1.0, 1.0, 0.95])
        np.testing.assert_allclose(seg["PMP1, Solvent B"]["y"], [0.0, 12.5, 99.9])
        np.testing.assert_allclose(seg["THM1, Temperature (Left)"]["y"], [29.91, 30.0])
        np.testing.assert_allclose(seg["PMP1, Pressure"]["t"], [0, 0.005, 0.01])
        self.assertAlmostEqual(seg["THM1, Temperature (Left)"]["t"][1], 1 / 60)
        self.assertEqual(seg["PMP1, Pressure"]["unit"], "bar")

    def test_file_non_registro(self):
        p = os.path.join(self.tmp.name, "x.REG")
        with open(p, "wb") as f:
            f.write(b"niente" * 100)
        with self.assertRaises(ValueError):
            self.app.leggi_diagnostica(p)
        with self.assertRaises(ValueError):
            self.app.leggi_acqres(p)

    def test_tabella_strumento(self):
        t_veloce = np.arange(0, 3) * 0.005
        self.app.strumento = {'D': {'segnali': {
            'PMP1, Pressure': {'t': t_veloce, 'y': np.array([1.0, 2.0, 3.0]), 'unit': 'bar'},
            'THM1, Temperature (Left)': {'t': np.array([0.0, 0.01]), 'y': np.array([30.0, 31.0]), 'unit': '\xb0C'}}}}
        df = self.app.tabella_strumento('D')
        self.assertEqual(df.index.name, 'Time (min)')
        np.testing.assert_allclose(df['PMP1, Pressure [bar]'], [1, 2, 3])
        # la temperatura (piu' lenta) e' interpolata sulla scala dei tempi del segnale veloce
        np.testing.assert_allclose(df['THM1, Temperature (Left) [\xb0C]'], [30.0, 30.5, 31.0])


@unittest.skipUnless(os.path.isdir(CAMPIONE), "campione 009-0201.D non presente")
class TestCampioneReale(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = hplc.HPLCManager.__new__(hplc.HPLCManager)
        cls.seg = cls.app.leggi_diagnostica(os.path.join(CAMPIONE, "LCDIAG.REG"))
        cls.acq = cls.app.leggi_acqres(os.path.join(CAMPIONE, "ACQRES.REG"))

    def test_segnali_presenti(self):
        self.assertEqual(list(self.seg), ["PMP1, Pressure", "PMP1, Flow", "PMP1, Solvent A", "PMP1, Solvent B",
                                          "PMP1, Solvent C", "PMP1, Solvent D", "THM1, Temperature (Left)"])
        for t in ("PMP1, Pressure", "PMP1, Flow", "PMP1, Solvent A"):
            self.assertEqual(len(self.seg[t]['y']), 18000)
            self.assertAlmostEqual(self.seg[t]['t'][1], 0.005)
        self.assertEqual(len(self.seg["THM1, Temperature (Left)"]['y']), 5400)
        self.assertAlmostEqual(self.seg["PMP1, Pressure"]['t'][-1], 89.995)    # 90 min di corsa

    def test_coerente_con_report00_e_run_log(self):
        m = self.app.leggi_metadati_cartella(CAMPIONE)
        p = self.seg["PMP1, Pressure"]
        # Report00.CSV: Start Pressure 96.449997 bar (a pochi secondi dall'iniezione), Stop 104.720001 bar
        self.assertAlmostEqual(np.interp(4 / 60, p['t'], p['y']), 96.45, delta=1.0)
        self.assertAlmostEqual(p['y'][-1], 104.72, delta=0.5)
        self.assertTrue(m['Start Pressure'].startswith('96.449997'))
        # flusso 1 ml/min dall'inizio alla fine
        np.testing.assert_allclose(self.seg["PMP1, Flow"]['y'], 1.0)
        # RUN.LOG: temperatura colonna 29.9-30.0 C; la curva sta nello stesso intorno
        temp = self.seg["THM1, Temperature (Left)"]['y']
        self.assertGreaterEqual(temp.min(), 29.8)
        self.assertLessEqual(temp.max(), 30.2)
        self.assertEqual(self.seg["THM1, Temperature (Left)"]['unit'], '\xb0C')

    def test_gradiente_somma_100(self):
        tot = sum(self.seg["PMP1, Solvent %s" % c]['y'] for c in "ABCD")
        self.assertGreater(tot.min(), 99.0)
        self.assertLess(tot.max(), 101.0)
        self.assertEqual(self.seg["PMP1, Solvent A"]['y'][0], 100.0)       # parte da 100% A
        self.assertEqual(self.seg["PMP1, Solvent B"]['y'][0], 0.0)
        self.assertGreater(self.seg["PMP1, Solvent B"]['y'].max(), 60.0)   # gradiente verso B

    def test_moduli_e_colonna(self):
        nomi = [m['name'] for m in self.acq['moduli']]
        self.assertEqual(nomi, ['1100 Quaternary Pump', '1100 Autosampler', '1100 Diode Array Detector',
                                '1100 Column Thermostat'])
        det = self.acq['moduli'][2]
        self.assertEqual((det['part'], det['firmware'], det['build']), ('G1315B', 'A.05.09', '007'))
        self.assertRegex(det['serial'], r'^DE\d{8}$')          # numero di serie: non riportato nel repo
        self.assertEqual(self.acq['moduli'][0]['part'], 'G1311A')
        c = self.acq['colonna']
        self.assertEqual(c['description'], 'ODS Hypersil')
        self.assertEqual((c['length_mm'], c['diameter_mm'], c['particle_um']), (100.0, 2.1, 5.0))

    def test_metadati(self):
        m = self.app.leggi_metadati_cartella(CAMPIONE)
        self.assertEqual(m['Column'], 'ODS Hypersil')
        self.assertEqual(m['Column size'], '100 x 2.1 mm, 5 um')
        self.assertRegex(m['1100 Diode Array Detector'], r'^G1315B, S/N DE\d{8}, FW A\.05\.09$')


@unittest.skipUnless(os.path.isdir(CAMPIONE), "campione 009-0201.D non presente")
class TestApplicazione(unittest.TestCase):
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

    def test_caricamento(self):
        self.assertIn('009-0201', self.app.strumento)
        d = self.app.strumento['009-0201']
        self.assertEqual(len(d['segnali']), 7)
        self.assertEqual(len(d['moduli']), 4)
        # i metadati della traccia estratta mostrano la colonna
        n = list(self.app.cromatogrammi)[0]
        self.assertEqual(self.app.cromatogrammi[n]['info']['Column'], 'ODS Hypersil')

    def test_finestra_curve(self):
        w = self.app.apri_curve_strumento()
        self.root.update()
        etichette = [ax.get_ylabel() for ax in w.fig.axes]
        self.assertEqual(etichette, ['Pressure (bar)', 'Solvent composition (%)', 'Flow (ml/min)',
                                     'Temperature (\xb0C)'])
        # i solventi C e D restano a zero: nel pannello della composizione ci sono solo A e B
        self.assertEqual([l.get_label() for l in w.fig.axes[1].lines], ['Solvent A', 'Solvent B'])
        self.assertEqual(len(w.fig.axes[0].lines[0].get_xdata()), 18000)

    def test_esporta_curve(self):
        import pandas as pd
        p = os.path.join(tempfile.mkdtemp(), 'curve.csv')
        with mock.patch.object(hplc.filedialog, 'asksaveasfilename', return_value=p), \
                mock.patch.object(hplc.messagebox, 'showinfo'), \
                mock.patch.object(hplc.messagebox, 'showerror') as err:
            self.app.esporta_curve_strumento()
        self.assertFalse(err.called)
        df = pd.read_csv(p, sep=';', encoding='latin-1', index_col=0)
        self.assertEqual(df.shape, (18000, 7))
        self.assertAlmostEqual(df['PMP1, Pressure [bar]'].iloc[0], 94.42)
        self.assertAlmostEqual(df['PMP1, Flow [ml/min]'].iloc[100], 1.0)

    def test_sessione_porta_le_curve(self):
        import pickle
        dati = pickle.loads(pickle.dumps(self.app._dati_sessione()))   # come passando da file
        self.assertIn('009-0201', dati['strumento'])
        self.app.clear_all()
        self.assertEqual(self.app.strumento, {})
        self.app._applica_sessione(dati)
        self.assertIn('009-0201', self.app.strumento)

    def test_senza_curve(self):
        self.app.strumento.clear()
        with mock.patch.object(hplc.messagebox, 'showinfo') as info:
            self.assertIsNone(self.app.apri_curve_strumento())
        self.assertTrue(info.called)


if __name__ == "__main__":
    unittest.main()
