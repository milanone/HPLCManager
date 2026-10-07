"""Test del lettore ChemStation (.ch versione 30 e REPORTnn.CSV) di HPLCManager.

Due livelli: un file .ch sintetico costruito con lo stesso schema a delta/assoluti (non servono
dati veri ne' Tk) e, se presente, la cartella reale in 'campioni di esempio/009-0201.D'
(gitignorata) confrontata con le tabelle picchi scritte da ChemStation.

    py -m unittest discover -s tests -v
"""
import importlib.util
import os
import struct
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "hplc", os.path.join(HERE, "..", "HPLCManager.pyw"))
hplc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hplc)

CAMPIONE = os.path.join(HERE, "..", "campioni di esempio", "009-0201.D")


def scrivi_ch(path, valori, t0=-2560, dt=400, campione=b"Campione test"):
    """Scrive un .ch v30 minimo: intestazione a zeri + dati a segmenti da 25 voci."""
    b = bytearray(0x400)
    b[0:3] = b"\x0230"
    b[0x18:0x19 + len(campione)] = bytes([len(campione)]) + campione
    b[0x11A:0x122] = struct.pack(">ii", t0, t0 + dt * (len(valori) - 1))
    voci, prev = [], 0
    for v in valori:
        d = v - prev
        if -0x7FFF <= d <= 0x7FFF:
            voci.append(struct.pack(">h", d))
        else:
            voci.append(b"\x80\x00" + struct.pack(">i", v))
        prev = v
    for i in range(0, len(voci), 25):
        blocco = voci[i:i + 25]
        b += bytes([0x10, len(blocco)]) + b"".join(blocco)
    with open(path, "wb") as f:
        f.write(b)


class TestLeggiChSintetico(unittest.TestCase):
    def setUp(self):
        self.app = hplc.HPLCManager.__new__(hplc.HPLCManager)

    def test_delta_e_assoluti(self):
        # il salto da 10 a 200000 non sta in un int16: deve usare l'escape 0x8000
        valori = [0, 5, 10, 200000, 200010, 199990, 3, 0] * 7
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.ch")
            scrivi_ch(p, valori)
            t, y, info = self.app.leggi_ch(p)
        np.testing.assert_allclose(y, np.array(valori) / 2000.0)
        self.assertEqual(len(t), len(valori))
        self.assertAlmostEqual(t[0], -2560 / 60000.0)
        self.assertAlmostEqual(t[1] - t[0], 400 / 60000.0)
        self.assertEqual(info["Sample"], "Campione test")

    def test_versione_non_supportata(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.ch")
            with open(p, "wb") as f:
                f.write(b"\x02130" + bytes(0x500))
            with self.assertRaises(ValueError):
                self.app.leggi_ch(p)


@unittest.skipUnless(os.path.isdir(CAMPIONE), "campione 009-0201.D non presente")
class TestConReportChemStation(unittest.TestCase):
    """Ogni picco di REPORTnn.CSV deve cadere dove il segnale .ch ha il suo apice."""

    CANALI = {"dad1A.ch": "REPORT01.CSV", "dad1B.ch": "REPORT02.CSV", "dad1C.ch": "REPORT03.CSV",
              "dad1D.ch": "REPORT04.CSV", "dad1E.ch": "REPORT05.CSV"}

    def setUp(self):
        self.app = hplc.HPLCManager.__new__(hplc.HPLCManager)

    def test_intestazione(self):
        t, y, info = self.app.leggi_ch(os.path.join(CAMPIONE, "dad1A.ch"))
        self.assertEqual(len(y), 13501)
        self.assertAlmostEqual(t[1] - t[0], 400 / 60000.0)
        self.assertAlmostEqual(t[-1], 5397440 / 60000.0)   # circa 89.96 min
        self.assertEqual(info["Sample"], "Estr acq gambe")
        self.assertEqual(info["Module"], "G1315B")
        self.assertEqual(info["Method"], "POLIFENB.M")

    def test_picchi_alti_del_report(self):
        # Solo i picchi alti e prima di 30 min: li' la baseline di ChemStation sta pochi mAU
        # sotto il segnale. I picchi piccoli (baseline di integrazione sotto il minimo locale) e
        # quelli oltre 80 min (deriva negativa in coda al gradiente) non sono confrontabili
        # senza replicare l'integratore.
        controllati = 0
        for ch, rep in self.CANALI.items():
            t, y, _ = self.app.leggi_ch(os.path.join(CAMPIONE, ch))
            for p in self.app.leggi_report_csv(os.path.join(CAMPIONE, rep)):
                if p["altezza"] < 50 or p["rt"] > 30:
                    continue
                with self.subTest(canale=ch, picco=p["n"], rt=p["rt"]):
                    controllati += 1
                    m = np.abs(t - p["rt"]) <= 0.05
                    i = np.flatnonzero(m)[np.argmax(y[m])]
                    self.assertAlmostEqual(t[i], p["rt"], delta=0.01)   # apice al RT del report
                    rapporto = y[i] / p["altezza"]
                    self.assertTrue(1.0 <= rapporto <= 1.15, rapporto)  # segnale = altezza + baseline
        self.assertGreater(controllati, 20)

    def test_picco_principale(self):
        t, y, _ = self.app.leggi_ch(os.path.join(CAMPIONE, "dad1A.ch"))
        p = max(self.app.leggi_report_csv(os.path.join(CAMPIONE, "REPORT01.CSV")),
                key=lambda q: q["altezza"])
        self.assertAlmostEqual(p["rt"], 4.354, places=3)
        self.assertAlmostEqual(y.max(), 2100.5, delta=1.0)
        self.assertAlmostEqual(t[np.argmax(y)], p["rt"], delta=0.01)


class TestMetadatiCartella(unittest.TestCase):
    def setUp(self):
        self.app = hplc.HPLCManager.__new__(hplc.HPLCManager)

    def test_sintetico(self):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "Report00.CSV"), "w", encoding="latin-1", newline="") as f:
                f.write('"Sample Name","Campione X",""\r\n"Injection Date","30-Sep-26, 11:08:12",""\r\n'
                        '"Acq. Operator","30-Sep-26, 11:08:12",""\r\n"Location",9,""\r\n'
                        '"Inj Volume",20,"\xb5l"\r\n"Start Flow",1,"ml/min"\r\n'
                        '"Solvent 1","PMP1, Solvent A","H2O"\r\n"Column 1","Peak Number",""\r\n'
                        '"Signal 1","DAD1 A, Sig=280,4 Ref=360,100",""\r\n')
            with open(os.path.join(d, "RUN.LOG"), "w", encoding="latin-1") as f:
                f.write("1100 THM   1 Column temperature = 29.9 \xb0C   11:08:16 09/30/26\r\n"
                        "1100 THM   1 Column temperature = 30.1 \xb0C   12:38:16 09/30/26\r\n")
            m = self.app.leggi_metadati_cartella(d)
        self.assertEqual(m["Sample Name"], "Campione X")
        self.assertEqual(m["Location"], "9")
        self.assertEqual(m["Inj Volume"], "20 \xb5l")
        self.assertEqual(m["Start Flow"], "1 ml/min")
        self.assertEqual(m["Solvent 1"], "PMP1, Solvent A, H2O")
        self.assertEqual(m["Signal 1"], "DAD1 A, Sig=280,4 Ref=360,100")
        self.assertNotIn("Column 1", m)
        self.assertNotIn("Acq. Operator", m)       # contiene la data dell'iniezione, non l'operatore
        self.assertEqual(m["Column temperature"], "29.9-30.1 \xb0C")

    def test_cartella_vuota(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(self.app.leggi_metadati_cartella(d), {})

    @unittest.skipUnless(os.path.isdir(CAMPIONE), "campione 009-0201.D non presente")
    def test_campione_reale(self):
        m = self.app.leggi_metadati_cartella(CAMPIONE)
        self.assertEqual(m["Sample Name"], "Estr acq gambe")
        self.assertEqual(m["Location"], "9")
        self.assertEqual(m["Inj Volume"], "20 \xb5l")
        self.assertEqual(m["Start Pressure"], "96.449997 bar")
        self.assertEqual(m["Solvent 2"], "PMP1, Solvent B, ACN")
        self.assertEqual(m["Signal 5"], "DAD1 E, Sig=520,4 Ref=off")
        self.assertTrue(m["Data File"].endswith("009-0201.D"))
        self.assertEqual(m["Column temperature"], "29.9-30.0 \xb0C")


class TestLeggiReportCsv(unittest.TestCase):
    def test_riga(self):
        app = hplc.HPLCManager.__new__(hplc.HPLCManager)
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "REPORT01.CSV")
            with open(p, "w", encoding="latin-1", newline="") as f:
                f.write('1,2.78,"VV  ",0.085,100.7,16.9,0.18\r\n')
            r = app.leggi_report_csv(p)
        self.assertEqual(r, [{"n": 1, "rt": 2.78, "tipo": "VV", "larghezza": 0.085,
                              "area": 100.7, "altezza": 16.9, "area_pct": 0.18}])


if __name__ == "__main__":
    unittest.main()
