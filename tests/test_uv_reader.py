"""Test del lettore spettri DAD (.uv versione 31) e dell'estrazione per lunghezza d'onda.

Due livelli: un .uv sintetico (non servono dati veri ne' Tk) e, se presente, la cartella reale
'campioni di esempio/009-0201.D' (gitignorata): il cromatogramma estratto dal .uv a 280/295/320/350 nm
(banda 4 nm, riferimento 360/100 nm) deve riprodurre i segnali dad1A..D.ch calcolati da ChemStation.

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


def scrivi_uv(path, matrice, t0=240, dt=400, w0=4000, w1=12000, passo=40):
    """Scrive un .uv v31 minimo: intestazione a zeri + un record little-endian per spettro."""
    b = bytearray(0x200)
    b[0:3] = b"\x0231"
    for riga, ms in zip(matrice, (t0 + dt * np.arange(len(matrice)))):
        voci, prev = [], 0
        for v in riga:
            d = int(v) - prev
            if -0x7FFF <= d <= 0x7FFF:
                voci.append(struct.pack("<h", d))
            else:
                voci.append(b"\x00\x80" + struct.pack("<i", int(v)))
            prev = int(v)
        dati = b"".join(voci)
        lung = 22 + len(dati)
        b += struct.pack("<HHI", 0x43, lung, int(ms)) + struct.pack("<HHH", w0, w1, passo) + bytes(8) + dati
    b += bytes(100)    # coda dopo l'ultimo record
    with open(path, "wb") as f:
        f.write(b)


class TestLeggiUvSintetico(unittest.TestCase):
    def setUp(self):
        self.app = hplc.HPLCManager.__new__(hplc.HPLCManager)

    def test_delta_e_assoluti(self):
        rng = np.random.default_rng(1)
        m = rng.integers(-400, 400, size=(5, 201))
        m[2, 50:] += 600000          # salto che richiede il valore assoluto (0x8000 + int32)
        m[3] = np.arange(201) * 3000  # rampa: tutti i delta stanno in int16
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.uv")
            scrivi_uv(p, m)
            t, wl, S, info = self.app.leggi_uv(p)
        np.testing.assert_allclose(S, m / 2000.0)
        np.testing.assert_allclose(wl, 200 + 2 * np.arange(201))
        np.testing.assert_allclose(t, (240 + 400 * np.arange(5)) / 60000.0)
        self.assertEqual(info["Spectra"], "5 spectra, 200-600 nm, step 2.0 nm")

    def test_versione_non_supportata(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.uv")
            with open(p, "wb") as f:
                f.write(b"\x02131" + bytes(0x600))
            with self.assertRaises(ValueError):
                self.app.leggi_uv(p)

    def test_estrai_lambda(self):
        wl = 200 + 2 * np.arange(201, dtype=float)
        S = np.tile(wl, (3, 1)) * 0.0
        S[:, :] = 1.0                   # spettro piatto = 1 mAU
        S[1, :] += 10.0                 # il secondo spettro e' 11 mAU
        app = self.app
        y = app.estrai_lambda(wl, S, 280, 4)
        np.testing.assert_allclose(y, [1.0, 11.0, 1.0])
        y = app.estrai_lambda(wl, S, 280, 4, 360, 100)    # banda meno riferimento: tutto piatto -> 0
        np.testing.assert_allclose(y, 0.0, atol=1e-12)
        with self.assertRaises(ValueError):
            app.estrai_lambda(wl, S, 900, 4)               # fuori dallo spettro acquisito


@unittest.skipUnless(os.path.isdir(CAMPIONE), "campione 009-0201.D non presente")
class TestUvConChReale(unittest.TestCase):
    CANALI = {"dad1A.ch": 280, "dad1B.ch": 295, "dad1C.ch": 320, "dad1D.ch": 350}

    @classmethod
    def setUpClass(cls):
        cls.app = hplc.HPLCManager.__new__(hplc.HPLCManager)
        cls.t, cls.wl, cls.S, cls.info = cls.app.leggi_uv(os.path.join(CAMPIONE, "dad1.uv"))

    def test_struttura(self):
        self.assertEqual(self.S.shape, (13493, 201))
        self.assertEqual(self.wl[0], 200.0)
        self.assertEqual(self.wl[-1], 600.0)
        self.assertAlmostEqual(self.t[0], 240 / 60000.0)
        self.assertAlmostEqual(self.t[1] - self.t[0], 400 / 60000.0)

    def test_estrazione_come_chemstation(self):
        # il .ch parte a -2560 ms, il .uv a 240 ms: stessa scala dei tempi, 7 campioni di scarto
        for ch, lam in self.CANALI.items():
            with self.subTest(canale=ch):
                tc, yc, _ = self.app.leggi_ch(os.path.join(CAMPIONE, ch))
                y = self.app.estrai_lambda(self.wl, self.S, lam, 4, 360, 100)
                idx = np.round((self.t - tc[0]) * 60000.0 / 400.0).astype(int)
                np.testing.assert_allclose(tc[idx], self.t, atol=1e-6)
                scarto = y - yc[idx]
                # 4 canali, picco principale fino a 2100 mAU: scarto RMS sotto 0.3 mAU
                self.assertLess(scarto.std(), 0.3)
                self.assertLess(np.abs(scarto).max(), 8.0)


if __name__ == "__main__":
    unittest.main()
