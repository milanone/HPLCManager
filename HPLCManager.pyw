"""HPLCManager: visualizzazione e analisi di cromatogrammi HPLC (Agilent ChemStation)."""
import csv
import gzip
import os
import pickle
import re
import struct
import sys
import tkinter as tk
import traceback
from tkinter import ttk, filedialog, messagebox

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

try:
    from scipy.signal import find_peaks, savgol_filter
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False

try:
    from tkinterdnd2 import DND_FILES, TkinterDnD
    HAS_DND = True
except Exception:
    HAS_DND = False


def _carica_origin_style():
    import importlib.util as ilu
    here = os.path.dirname(os.path.abspath(__file__))
    for path in (os.path.join(here, 'origin_style.py'),
                 os.path.join(here, '..', 'PlotStyleKit', 'origin_style.py')):
        if os.path.isfile(path):
            spec = ilu.spec_from_file_location('origin_style', path)
            mod = ilu.module_from_spec(spec)
            spec.loader.exec_module(mod)
            return mod
    return None

try:
    origin_style = _carica_origin_style()
    if origin_style is not None:
        origin_style.applica_rcparams()
except Exception:
    origin_style = None


class HPLCManager:
    """Gestore di cromatogrammi HPLC.

    Modello dati: self.cromatogrammi = {nome: {'df': DataFrame indicizzato per tempo (min) con
    una colonna 'mAU', 'info': {...}, 'picchi': [...], 'picchi_cs': [...], 'nascosto': bool}},
    sulla falsariga di LabSpectrumManager. 'picchi' sono i picchi calcolati qui, 'picchi_cs'
    quelli scritti da ChemStation nel REPORTnn.CSV della stessa cartella .D (se c'e').
    """

    # Fattore di scala dei .ch di ChemStation: valore decodificato / 2000 = mAU
    # (dedotto dal confronto con Report.TXT, non letto dall'intestazione).
    CH_SCALA = 2000.0
    # Come i .ch, i valori del .uv sono in 1/2000 di mAU (verificato estraendo 280/295/320/350 nm
    # dal .uv e confrontando con dad1A..D.ch). L'estrazione di una banda usa l'asse delle lunghezze
    # d'onda spostato di UV_OFFSET_NM: valore empirico (minimo dello scarto rispetto ai .ch fra 0 e
    # 1 nm, ~0.45), non letto dall'intestazione.
    UV_SCALA = 2000.0
    UV_OFFSET_NM = 0.45

    def __init__(self, root):
        self.root = root
        self.root.title("HPLC Manager")
        self.root.geometry("1500x800")
        self.cromatogrammi = {}
        self._dirty = False
        self._cursor_artists = []
        self._tool_frame = None
        self.var_picchi_cs = tk.BooleanVar(value=True)
        self.spettri = {}          # {nome cartella .D: {'t', 'wl', 'S', 'info'}} dai file .uv
        self._spec_win = None
        self.strumento = {}        # {nome cartella .D: {'segnali': {titolo: {t, y, unit}}, 'moduli': [...], 'colonna': {...}}}
        self._strum_win = None
        self.var_click_spec = tk.BooleanVar(value=True)
        self._marcatori_artisti = []   # artisti matplotlib dei marcatori sul cromatogramma (per rimuoverli senza ridisegnare)
        self._sp_colore_n = 0          # contatore per assegnare un colore stabile a ogni spettro scelto
        self._sp_multi_ultimo = []     # valori dell'ultimo dialogo 'Several wavelengths' (si riaprono uguali)
        self._spettri_scelti = []      # spettri presi con click/Show: {'dataset','t','wl','y','colore','visibile','linea'}
        self._pe_module = None         # modulo plot_editor (PlotStyleKit), caricato al primo uso
        self.current_dir = os.getcwd()
        self._reset_vista = True       # prossimo ridisegno: assi automatici (nuovi dati) invece di conservare lo zoom
        self.var_includi_ch = tk.BooleanVar(value=False)   # carica anche i canali registrati (.ch)
        self._sp_live = {}         # {nome dataset: nome della traccia estratta 'viva'}

        self._costruisci_menu()

        self.paned = tk.PanedWindow(root, orient=tk.HORIZONTAL, sashrelief=tk.RAISED, sashwidth=4)
        self.paned.pack(fill=tk.BOTH, expand=True)

        # --- pannello sinistro: tabella dati ---
        self.f_left = tk.Frame(self.paned, width=400)
        tk.Label(self.f_left, text="Data Table", font=('Segoe UI', 10, 'bold')).pack(anchor='w', padx=4)
        self.tabella = ttk.Treeview(self.f_left, show='headings', selectmode='extended')
        sy = ttk.Scrollbar(self.f_left, orient='vertical', command=self.tabella.yview)
        sx = ttk.Scrollbar(self.f_left, orient='horizontal', command=self.tabella.xview)
        self.tabella.configure(yscrollcommand=sy.set, xscrollcommand=sx.set)
        sy.pack(side=tk.RIGHT, fill=tk.Y)
        sx.pack(side=tk.BOTTOM, fill=tk.X)
        self.tabella.pack(fill=tk.BOTH, expand=True)
        self.tabella.bind('<Control-c>', self._copia_tabella_dati)
        self.tabella.bind('<Control-a>', self._seleziona_tutto_tabella)
        self.paned.add(self.f_left, width=400, minsize=200, stretch='never')

        # --- pannello centrale: grafico ---
        self.f_plot = tk.Frame(self.paned, width=700)
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.f_plot)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.f_plot)
        self.toolbar.update()
        self.toolbar.pack(side=tk.TOP, fill=tk.X)
        self.plot_widget = self.canvas.get_tk_widget()
        self.plot_widget.pack(fill=tk.BOTH, expand=True)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('axes_leave_event', self.on_mouse_leave)
        self.canvas.mpl_connect('resize_event', self._ritaglia_margini)
        self.canvas.mpl_connect('button_press_event', self._premi_spettro)
        self.canvas.mpl_connect('button_release_event', self._click_spettro)
        self.paned.add(self.f_plot, width=700, minsize=300, stretch='always')

        # --- pannello destro: scrollabile, tutto dentro self.f_right_inner ---
        self.f_right = tk.Frame(self.paned, width=400)
        self.right_canvas = tk.Canvas(self.f_right, highlightthickness=0)
        rs = ttk.Scrollbar(self.f_right, orient='vertical', command=self.right_canvas.yview)
        self.right_canvas.configure(yscrollcommand=rs.set)
        rs.pack(side=tk.RIGHT, fill=tk.Y)
        self.right_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.f_right_inner = tk.Frame(self.right_canvas)
        self._right_win = self.right_canvas.create_window((0, 0), window=self.f_right_inner, anchor='nw')
        self.f_right_inner.bind('<Configure>', lambda e: self.right_canvas.configure(
            scrollregion=self.right_canvas.bbox('all')))
        self.right_canvas.bind('<Configure>', lambda e: self.right_canvas.itemconfigure(
            self._right_win, width=e.width))
        self.right_canvas.bind('<Enter>', lambda e: self.right_canvas.bind_all(
            '<MouseWheel>', lambda ev: self.right_canvas.yview_scroll(int(-ev.delta / 120), 'units')))
        self.right_canvas.bind('<Leave>', lambda e: self.right_canvas.unbind_all('<MouseWheel>'))
        self.paned.add(self.f_right, width=400, minsize=300, stretch='never')
        self._costruisci_pannello_destro()

        self.root.protocol("WM_DELETE_WINDOW", self.on_exit)
        if HAS_DND:
            try:
                self.root.drop_target_register(DND_FILES)
                self.root.dnd_bind('<<Drop>>', self.handle_drop)
            except Exception:
                pass
        self.aggiorna_vista()

    # ------------------------------------------------------------------ menu e pannelli
    def _costruisci_menu(self):
        self.menu_bar = tk.Menu(self.root)
        m = tk.Menu(self.menu_bar, tearoff=0)
        m.add_command(label="Open ChemStation folder (.D)...", command=self.carica_cartella_dialog)
        m.add_command(label="Open signal file (.ch)...", command=self.carica_file_dialog)
        m.add_checkbutton(label="Also load recorded channels (.ch) from .D folders",
                          variable=self.var_includi_ch)
        m.add_separator()
        m.add_command(label="Open Session...", command=self.apri_sessione)
        m.add_command(label="Save Session...", command=self.salva_sessione)
        m.add_separator()
        m.add_command(label="Export Traces (CSV)...", command=self.esporta_csv)
        m.add_command(label="Export Peaks (CSV)...", command=self.esporta_picchi)
        m.add_command(label="Export Spectra (CSV)...", command=self.esporta_spettri)
        m.add_command(label="Export Full DAD Data (CSV)...", command=self.esporta_dad_completo)
        m.add_command(label="Export Instrument Curves (CSV)...", command=self.esporta_curve_strumento)
        m.add_separator()
        m.add_command(label="Save Figure Image (PNG, PDF, SVG)...", command=self.salva_figura_immagine)
        m.add_command(label="Save Figure (pickle)...", command=self.salva_figura_pickle)
        m.add_command(label="Edit Figure...", command=self.apri_editor_figura)
        m.add_separator()
        m.add_command(label="Exit", command=self.on_exit)
        self.menu_bar.add_cascade(label="File", menu=m)

        m = tk.Menu(self.menu_bar, tearoff=0)
        m.add_command(label="Remove Selected", command=self.remove_selected)
        m.add_command(label="Clear All", command=self.clear_all)
        self.menu_bar.add_cascade(label="Edit", menu=m)

        m = tk.Menu(self.menu_bar, tearoff=0)
        m.add_command(label="Hide Selected", command=self.hide_selected)
        m.add_command(label="Show Selected", command=self.show_selected)
        m.add_separator()
        m.add_checkbutton(label="Show ChemStation peaks", variable=self.var_picchi_cs,
                          command=self.aggiorna_vista)
        self.menu_bar.add_cascade(label="View", menu=m)

        m = tk.Menu(self.menu_bar, tearoff=0)
        m.add_command(label="Peaks...", command=self.apri_picchi)
        m.add_command(label="Smoothing...", command=self.apri_smoothing)
        m.add_command(label="Trim...", command=self.apri_trim)
        m.add_command(label="Normalize", command=self.normalizza)
        m.add_separator()
        m.add_command(label="Spectra (DAD)...", command=self.apri_spettri)
        m.add_command(label="Several wavelengths...", command=self.apri_estrazione_multipla)
        m.add_command(label="Instrument curves (pressure, gradient...)", command=self.apri_curve_strumento)
        self.menu_bar.add_cascade(label="Tools", menu=m)
        self.root.config(menu=self.menu_bar)

    def _costruisci_pannello_destro(self):
        p = self.f_right_inner
        tk.Label(p, text="Chromatograms", font=('Segoe UI', 10, 'bold')).pack(anchor='w', padx=4, pady=(4, 0))
        fl = tk.Frame(p)
        fl.pack(fill=tk.X, padx=4)
        self.lista = tk.Listbox(fl, selectmode=tk.EXTENDED, height=6, exportselection=False)
        self.lista.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.lista.bind('<<ListboxSelect>>', self._su_selezione)
        self.lista.bind('<Button-3>', self._menu_lista)
        fb = tk.Frame(p)
        fb.pack(fill=tk.X, padx=4, pady=2)
        for txt, cmd in (("Remove", self.remove_selected), ("Clear", self.clear_all),
                         ("Hide", self.hide_selected), ("Show", self.show_selected)):
            tk.Button(fb, text=txt, command=cmd).pack(side=tk.LEFT, padx=1)
        tk.Button(p, text="Open ChemStation folder (.D)...",
                  command=self.carica_cartella_dialog).pack(fill=tk.X, padx=4, pady=2)

        tk.Label(p, text="Metadata", font=('Segoe UI', 10, 'bold')).pack(anchor='w', padx=4, pady=(6, 0))
        self.meta = tk.Text(p, height=12, wrap='none', state='disabled', font=('Consolas', 9))
        self.meta.pack(fill=tk.X, padx=4)
        sm = ttk.Scrollbar(p, orient='horizontal', command=self.meta.xview)
        self.meta.configure(xscrollcommand=sm.set)
        sm.pack(fill=tk.X, padx=4)
        self.meta.bind('<MouseWheel>', lambda e: (self.meta.yview_scroll(int(-e.delta / 120), 'units'), 'break')[1])

        tk.Label(p, text="Tools", font=('Segoe UI', 10, 'bold')).pack(anchor='w', padx=4, pady=(6, 0))
        fb = tk.Frame(p)
        fb.pack(fill=tk.X, padx=4)
        for txt, cmd in (("Peaks", self.apri_picchi), ("Smoothing", self.apri_smoothing),
                         ("Trim", self.apri_trim), ("Normalize", self.normalizza),
                         ("Spectra", self.apri_spettri)):
            tk.Button(fb, text=txt, command=cmd).pack(side=tk.LEFT, padx=1)
        tk.Button(p, text="Instrument curves (pressure, gradient...)",
                  command=self.apri_curve_strumento).pack(fill=tk.X, padx=4, pady=(3, 0))
        self.f_tool_host = tk.Frame(p)
        self.f_tool_host.pack(fill=tk.X, padx=4, pady=4)

    def _menu_lista(self, event):
        i = self.lista.nearest(event.y)
        if i >= 0 and i not in self.lista.curselection():
            self.lista.selection_clear(0, tk.END)
            self.lista.selection_set(i)
        menu = tk.Menu(self.root, tearoff=0)
        menu.add_command(label="Hide", command=self.hide_selected)
        menu.add_command(label="Show", command=self.show_selected)
        menu.add_command(label="Remove", command=self.remove_selected)
        menu.tk_popup(event.x_root, event.y_root)

    def _apri_pannello_tool(self, titolo):
        """Svuota l'area dei tool nel pannello destro e ritorna il frame in cui costruirne uno nuovo."""
        for w in self.f_tool_host.winfo_children():
            w.destroy()
        f = tk.LabelFrame(self.f_tool_host, text=titolo)
        f.pack(fill=tk.X)
        self._tool_frame = f
        return f

    # ------------------------------------------------------------------ lettura file
    def leggi_ch(self, path):
        """Legge un segnale ChemStation .ch versione 30 (big-endian).

        Ritorna (tempo_min, segnale_mAU, info). I dati partono da 0x400 in segmenti
        '0x10 n' seguiti da n voci: int16 = delta dal valore precedente, oppure
        0x8000 + int32 = valore assoluto (che azzera l'accumulo).
        """
        with open(path, 'rb') as f:
            b = f.read()
        if b[0:3] != b'\x0230':
            raise ValueError("File .ch non supportato (versione diversa da 30): "
                             + os.path.basename(path))

        def pascal(o):
            return b[o + 1:o + 1 + b[o]].decode('latin-1')

        # stringhe Pascal a offset fissi (versione, tipo file, nome campione, data, modulo, metodo)
        info = {'Version': pascal(0), 'Type': pascal(4), 'Sample': pascal(0x18)}
        info['Date'] = pascal(0xB2)
        info['Module'] = pascal(0xD0)
        info['Method'] = pascal(0xE4)
        # descrizione del segnale ("DAD A, Sig=280,4 Ref=360,100"): stringa Pascal prima di 0x400
        i = b.find(b'Sig=', 0, 0x400)
        if i > 0:
            s = i
            while s > 1 and b[s - 1] >= 32:
                s -= 1
            info['Signal'] = b[s:s + b[s - 1]].decode('latin-1')

        t0, t1 = struct.unpack('>ii', b[0x11A:0x122])
        dati = []
        v = 0
        pos = 0x400
        while pos + 2 <= len(b) and b[pos] == 0x10:
            n = b[pos + 1]
            pos += 2
            for _ in range(n):
                if pos + 2 > len(b):
                    break
                if b[pos:pos + 2] == b'\x80\x00':
                    v = struct.unpack('>i', b[pos + 2:pos + 6])[0]
                    pos += 6
                else:
                    v += struct.unpack('>h', b[pos:pos + 2])[0]
                    pos += 2
                dati.append(v)
        y = np.array(dati, dtype=float) / self.CH_SCALA
        dt = (t1 - t0) / (len(y) - 1) if len(y) > 1 else 0.0
        t = (t0 + np.arange(len(y)) * dt) / 60000.0
        return t, y, info

    def leggi_report_csv(self, path):
        """Legge un REPORTnn.CSV di ChemStation (una tabella picchi per segnale, senza intestazione).

        Ritorna una lista di dict: n, rt, tipo, larghezza, area, altezza, area_pct.
        """
        picchi = []
        with open(path, newline='', encoding='latin-1') as f:
            for r in csv.reader(f):
                if len(r) < 7:
                    continue
                picchi.append({'n': int(r[0]), 'rt': float(r[1]), 'tipo': r[2].strip(),
                               'larghezza': float(r[3]), 'area': float(r[4]),
                               'altezza': float(r[5]), 'area_pct': float(r[6])})
        return picchi

    def leggi_uv(self, path):
        """Legge un file spettri DAD .uv di ChemStation (versione 31).

        Ritorna (tempo_min, lunghezze_onda_nm, S, info) con S di forma (n_tempi, n_lunghezze) in mAU.
        Dopo l'intestazione (record a partire da 0x200) ogni spettro e' un record little-endian:
        uint16 tipo (0x43), uint16 lunghezza del record, uint32 tempo (ms); poi uint16 lambda
        iniziale, finale e passo (in 1/20 nm), 8 byte non decodificati; poi n valori codificati come
        nei .ch ma little-endian: int16 = delta dal precedente, 0x8000 + int32 = valore assoluto.
        Dopo l'ultimo record c'e' una coda (indice?) non letta.
        """
        with open(path, 'rb') as f:
            b = f.read()
        if b[0:3] != b'\x0231':
            raise ValueError("File .uv non supportato (versione diversa da 31): "
                             + os.path.basename(path))
        def pascal(o):
            return b[o + 1:o + 1 + b[o]].decode('latin-1')

        # stesse stringhe Pascal a offset fissi dei .ch (nome campione, data, modulo, metodo)
        info = {'Version': pascal(0), 'Type': pascal(4), 'Sample': pascal(0x18),
                'Date': pascal(0xB2), 'Module': pascal(0xD0), 'Method': pascal(0xE4)}
        tempi, righe = [], []
        s, wl = 0x200, None
        while s + 22 <= len(b):
            tipo, lung, ms = struct.unpack('<HHI', b[s:s + 8])
            if tipo != 0x43 or lung < 22 or s + lung > len(b):
                break
            w0, w1, st = struct.unpack('<HHH', b[s + 8:s + 14])
            if st == 0:
                break
            n = (w1 - w0) // st + 1
            if wl is None:
                wl = (w0 + st * np.arange(n)) / 20.0
            elif n != len(wl):
                break
            dati = b[s + 22:s + lung]
            if len(dati) == 2 * n:                       # nessun valore assoluto nel record
                v = np.cumsum(np.frombuffer(dati, dtype='<i2').astype(np.int64))
            else:
                v, val, o = np.empty(n, dtype=np.int64), 0, 0
                for k in range(n):
                    if dati[o:o + 2] == b'\x00\x80':
                        val = struct.unpack('<i', dati[o + 2:o + 6])[0]
                        o += 6
                    else:
                        val += struct.unpack('<h', dati[o:o + 2])[0]
                        o += 2
                    v[k] = val
            tempi.append(ms)
            righe.append(v)
            s += lung
        if not righe:
            raise ValueError("Nessuno spettro leggibile in " + os.path.basename(path))
        S = np.array(righe, dtype=float) / self.UV_SCALA
        info['Spectra'] = '%d spectra, %.0f-%.0f nm, step %.1f nm' % (len(righe), wl[0], wl[-1], wl[1] - wl[0])
        return np.array(tempi, dtype=float) / 60000.0, wl, S, info

    def estrai_lambda(self, wl, S, centro, banda, rif=None, banda_rif=None):
        """Cromatogramma estratto a `centro` (nm) con larghezza `banda`, meno la banda di riferimento.

        Come il segnale 'Sig=280,4 Ref=360,100' di ChemStation: media pesata della banda
        [centro - banda/2, centro + banda/2], con peso proporzionale alla sovrapposizione fra la
        banda e l'intervallo spettrale di ciascun punto (passo wl[1]-wl[0]); stessa cosa per il
        riferimento. Ritorna il segnale in mAU, un valore per ogni spettro.
        """
        passo = float(wl[1] - wl[0])
        x = wl + self.UV_OFFSET_NM

        def pesi(c, bw):
            lo, hi = c - bw / 2.0, c + bw / 2.0
            w = np.clip(np.minimum(x + passo / 2, hi) - np.maximum(x - passo / 2, lo), 0, None)
            if w.sum() == 0:
                raise ValueError("Banda %.1f +/- %.1f nm fuori dallo spettro acquisito (%.0f-%.0f nm)"
                                 % (c, bw / 2.0, wl[0], wl[-1]))
            return w / w.sum()

        y = S @ pesi(centro, banda)
        if rif is not None and banda_rif:
            y = y - S @ pesi(rif, banda_rif)
        return y

    # Unita' e fattori di scala dei segnali di LCDIAG.REG: il valore memorizzato (uint32) e' il valore
    # fisico moltiplicato per questo fattore (verificato su pressione, flusso, solventi, temperatura).
    UNITA_DIAG = {'bar': 100.0, 'ml/min': 1000.0, '%': 10.0, '\xb0C': 100.0}

    def leggi_diagnostica(self, path):
        """Legge LCDIAG.REG di ChemStation: i profili registrati dallo strumento durante la corsa
        (pressione, flusso, composizione dei solventi, temperatura colonna).

        Ritorna {titolo: {'t': minuti (array), 'y': valori (array), 'unit': 'bar'|'ml/min'|'%'|'\xb0C'}}
        in ordine di apparizione. Struttura (little-endian): per ogni segnale un array di uint32
        preceduto dalle stringhe 'min\\0<unita'>\\0'; il numero di punti sta a 168 byte prima di 'min'
        (uint32), l'intervallo di campionamento in minuti a 116 byte prima (double). Il titolo
        ('PMP1, Pressure') e' nella testata che SEGUE l'array. Il tempo parte da 0 (inizio corsa):
        e' un'ipotesi, coerente con la pressione a 4 s e la pressione finale del report.
        """
        with open(path, 'rb') as f:
            b = f.read()
        if b[3:4] != b'\x00' or b'REGISTER FILE' not in b[:24]:
            raise ValueError("File .REG non riconosciuto: " + os.path.basename(path))
        segnali = {}
        for k, m in enumerate(re.finditer(rb'min\x00(bar|ml/min|%|\xb0C)\x00', b)):
            u, e, unita = m.start(), m.end(), m.group(1).decode('latin-1')
            if u < 168:
                continue
            n = struct.unpack('<I', b[u - 168:u - 164])[0]
            dt = struct.unpack('<d', b[u - 116:u - 108])[0]
            if not (0 < n <= 10000000) or e + 4 * n > len(b) or not (0 < dt < 10):
                continue
            y = np.frombuffer(b[e:e + 4 * n], dtype='<u4').astype(float) / self.UNITA_DIAG[unita]
            testata = b[e + 4 * n:e + 4 * n + 1500]
            t = re.search(rb'((?:[A-Z]{3}\d), [^\x00]+)\x00arial', testata)
            titolo = t.group(1).decode('latin-1') if t else '%s signal %d' % (unita, k + 1)
            if titolo in segnali:
                titolo = '%s (%d)' % (titolo, k + 1)
            segnali[titolo] = {'t': np.arange(n) * dt, 'y': y, 'unit': unita}
        if not segnali:
            raise ValueError("Nessun segnale riconosciuto in " + os.path.basename(path))
        return segnali

    def leggi_acqres(self, path):
        """Legge ACQRES.REG di ChemStation: moduli dello strumento e dati della colonna.

        Ritorna {'moduli': [{'name', 'part', 'serial', 'firmware', 'build'}], 'colonna': {...}}.
        Registro con tabelle a colonne: le stringhe sono uint16 (lunghezza con NUL) + testo + NUL.
        Moduli: 5 colonne di stringhe consecutive (serialNumber, FWrevision, buildNumber, Name,
        PartNumber), una riga per modulo. Colonna: descrizione = la stringa che non e' versione del
        software, nome strumento o percorso; lunghezza, diametro e granulometria sono i primi tre double
        'tondi' della tabella 'Acquisition Results' (ColLength, ColDiameter, ParticleSize: assegnazione
        dedotta dall'ordine dei campi e da valori tipici, 100 x 2.1 mm, 5 um, non verificata su altri file).
        """
        with open(path, 'rb') as f:
            b = f.read()
        if b[3:4] != b'\x00' or b'REGISTER FILE' not in b[:24]:
            raise ValueError("File .REG non riconosciuto: " + os.path.basename(path))
        testo_ok = re.compile(r'^[\x20-\x7e\xa0-\xff]+$')

        def stringa(o):
            if o < 0 or o + 2 > len(b):
                return None
            n = struct.unpack('<H', b[o:o + 2])[0]
            if n < 2 or n > 200 or o + 2 + n > len(b) or b[o + 1 + n] != 0:
                return None
            t = b[o + 2:o + 1 + n].decode('latin-1')
            return t if testo_ok.match(t) and any(c.isalnum() for c in t) else None

        def catena(o):
            out = []
            while True:
                t = stringa(o)
                if t is None:
                    return out
                out.append((o, t))
                o += 2 + len(t) + 1
        migliore = []
        o = 0
        while o < len(b) - 2:
            c = catena(o)
            if len(c) > len(migliore):
                migliore = c
            o += max(1, sum(2 + len(t) + 1 for _, t in c) if c else 1)
        moduli = []
        colonne = [b'serialNumber', b'FWrevision', b'buildNumber', b'Name', b'PartNumber']
        posizioni = [b.find(c) for c in colonne]
        if len(migliore) >= 5 and len(migliore) % 5 == 0 and all(x >= 0 for x in posizioni) \
                and posizioni == sorted(posizioni):
            r = len(migliore) // 5
            val = [t for _, t in migliore]
            for i in range(r):
                moduli.append({'serial': val[i], 'firmware': val[r + i], 'build': val[2 * r + i],
                               'name': val[3 * r + i], 'part': val[4 * r + i]})
        fine_catena = (migliore[-1][0] + 3 + len(migliore[-1][1])) if migliore else 0
        inizio_catena = migliore[0][0] if migliore else 0
        altre = []
        o = 0
        while o < len(b) - 2:
            t = stringa(o)
            if t is not None and len(t) >= 3 and not (inizio_catena <= o < fine_catena):
                altre.append(t)
                o += 2 + len(t) + 1
            else:
                o += 1
        colonna = {}
        for t in altre:
            if t.startswith('Rev.') or '\\' in t or t.startswith('Instrument') or t.endswith('Results'):
                continue
            if re.fullmatch(r'[A-Za-z][A-Za-z0-9 ./-]{3,}', t) and not re.fullmatch(r'[A-Za-z]+[A-Z][a-z]+[A-Za-z]*', t):
                colonna['description'] = t
                break
        ini = b.find(b'WriteInfo')
        fine = inizio_catena if inizio_catena else len(b)
        valori = []
        o = ini + 9 if ini >= 0 else 0
        while o + 8 <= fine:
            v = struct.unpack('<d', b[o:o + 8])[0]
            if v == v and 0.4 <= abs(v) <= 1000 and abs(v * 10 - round(v * 10)) < 1e-9:
                valori.append(v)
                o += 8
            else:
                o += 1
        if len(valori) >= 2:
            colonna['length_mm'], colonna['diameter_mm'] = valori[0], valori[1]
        if len(valori) >= 3 and valori[2] <= 20:
            colonna['particle_um'] = valori[2]
        return {'moduli': moduli, 'colonna': colonna}

    def leggi_metadati_cartella(self, cartella):
        """Metadati dell'analisi dalla cartella .D: Report00.CSV (campione, iniezione, sequenza,
        pressione, flusso, solventi, segnali) e RUN.LOG (temperatura colonna).

        Ritorna un dict {etichetta: valore} in ordine di apparizione; vuoto se i file mancano.
        """
        info = {}
        rep00 = os.path.join(cartella, 'Report00.CSV')
        if os.path.isfile(rep00):
            saltare = ('Report Title', 'Sorted By', 'Multiplier', 'Dilution', 'Number of Columns',
                       'Method Info', 'Results Created by')
            with open(rep00, newline='', encoding='latin-1') as f:
                for r in csv.reader(f):
                    if len(r) < 2 or not r[0].strip() or r[0].startswith('Column ') or r[0] in saltare:
                        continue
                    valore = ''.join(x for x in r[1:]).strip()
                    if r[0] in ('Data File', 'Analysis Method', 'Sequence File'):
                        valore = ''.join(r[1:3])               # percorso + nome file
                    elif len(r) > 2 and r[2].strip() and r[0].startswith(('Inj Vol', 'Actual', 'Start', 'Stop')):
                        valore = '%s %s' % (r[1].strip(), r[2].strip())
                    elif r[0].startswith('Solvent'):
                        # le linee senza solvente (terzo campo vuoto) non si mostrano
                        valore = ', '.join(x.strip() for x in r[1:] if x.strip()) if len(r) > 2 and r[2].strip() else ''
                    elif r[0].startswith('Signal '):
                        valore = ', '.join(x.strip() for x in r[1:] if x.strip())
                    if valore:
                        info[r[0]] = valore
            # in questo export 'Acq. Operator' contiene la data dell'iniezione, non l'operatore
            if info.get('Acq. Operator') == info.get('Injection Date'):
                del info['Acq. Operator']
        log = os.path.join(cartella, 'RUN.LOG')
        if os.path.isfile(log):
            with open(log, encoding='latin-1') as f:
                temp = [float(x) for x in re.findall(r'Column temperature\s*=\s*([\d.]+)', f.read())]
            if temp:
                info['Column temperature'] = ('%.1f' % temp[0] if min(temp) == max(temp)
                                              else '%.1f-%.1f' % (min(temp), max(temp))) + ' \u00b0C'
        acq = os.path.join(cartella, 'ACQRES.REG')
        if os.path.isfile(acq):
            try:
                r = self.leggi_acqres(acq)
            except Exception:
                r = None
            if r:
                c = r['colonna']
                if c.get('description'):
                    info['Column'] = c['description']
                if 'length_mm' in c:
                    info['Column size'] = '%g x %g mm%s' % (c['length_mm'], c['diameter_mm'], (
                        ', %g um' % c['particle_um']) if 'particle_um' in c else '')
                for m in r['moduli']:
                    info[m['name']] = '%s, S/N %s, FW %s' % (m['part'], m['serial'], m['firmware'])
        return info

    def _unisci_metadati(self, info, cartella):
        """Aggiunge a `info` i metadati della cartella .D, senza ripetere quelli gia' letti
        dall'intestazione binaria (stesso valore sotto un'altra etichetta)."""
        extra = self.leggi_metadati_cartella(cartella)
        for nuovo, vecchio in (('Sample Name', 'Sample'), ('Injection Date', 'Date'),
                               ('Acq. Method', 'Method')):
            if extra.get(nuovo) == info.get(vecchio):
                del extra[nuovo]
        info.update(extra)

    def _carica_strumento(self, cartella, base):
        """Carica i profili dello strumento (LCDIAG.REG) e i dati di moduli e colonna (ACQRES.REG)."""
        dati = {'segnali': {}, 'moduli': [], 'colonna': {}}
        diag = os.path.join(cartella, 'LCDIAG.REG')
        if os.path.isfile(diag):
            try:
                dati['segnali'] = self.leggi_diagnostica(diag)
            except Exception as e:
                traceback.print_exc()
                messagebox.showwarning("Strumento", "Impossibile leggere LCDIAG.REG:\n%s" % e)
        acq = os.path.join(cartella, 'ACQRES.REG')
        if os.path.isfile(acq):
            try:
                r = self.leggi_acqres(acq)
                dati['moduli'], dati['colonna'] = r['moduli'], r['colonna']
            except Exception:
                traceback.print_exc()
        if dati['segnali'] or dati['moduli']:
            self.strumento[base] = dati

    def _carica_uv(self, cartella, base):
        """Carica gli spettri DAD (.uv) di una cartella .D, se ci sono."""
        uv = [f for f in os.listdir(cartella) if f.lower().endswith('.uv')]
        if not uv:
            return
        try:
            t, wl, S, info = self.leggi_uv(os.path.join(cartella, sorted(uv)[0]))
        except Exception as e:
            messagebox.showwarning("Spettri", "Impossibile leggere gli spettri DAD:\n%s" % e)
            return
        self._unisci_metadati(info, cartella)
        self.spettri[base] = {'t': t, 'wl': wl, 'S': S, 'info': info}

    def _picchi_chemstation(self, cartella, canale):
        """Picchi di ChemStation per il canale (es. 'dad1A.ch') nella cartella .D, oppure []."""
        lettera = os.path.splitext(canale)[0][-1].upper()
        numero = None
        rep = os.path.join(cartella, 'Report.TXT')
        if os.path.isfile(rep):
            with open(rep, encoding='latin-1') as f:
                for m in re.finditer(r'Signal\s+(\d+)\s*:\s*\w+?\d*\s+([A-Z])\b', f.read()):
                    if m.group(2) == lettera:
                        numero = int(m.group(1))
                        break
        if numero is None and 'A' <= lettera <= 'Z':
            numero = ord(lettera) - ord('A') + 1
        path = os.path.join(cartella, 'REPORT%02d.CSV' % (numero or 0))
        if numero and os.path.isfile(path):
            try:
                return self.leggi_report_csv(path)
            except Exception:
                return []
        return []

    def _nome_libero(self, nome):
        if nome not in self.cromatogrammi:
            return nome
        k = 2
        while '%s (%d)' % (nome, k) in self.cromatogrammi:
            k += 1
        return '%s (%d)' % (nome, k)

    def _aggiungi(self, nome, t, y, info, picchi_cs=None, estr=None):
        """Aggiunge una traccia. `estr` (solo per le tracce estratte dagli spettri DAD) sono i parametri
        con cui e' stata calcolata: {'ds', 'l', 'b', 'rif', 'rb'}; selezionandola nella lista tornano
        nel pannello, e 'Update trace' la modifica."""
        nome = self._nome_libero(nome)
        df = pd.DataFrame({'mAU': y}, index=pd.Index(t, name='Time (min)'))
        self.cromatogrammi[nome] = {'df': df, 'info': info, 'picchi': [],
                                    'picchi_cs': picchi_cs or [], 'nascosto': False}
        if estr is not None:
            self.cromatogrammi[nome]['estr'] = estr
        return nome

    def processa_file(self, path):
        """Carica un .ch singolo, una cartella .D (tutti i suoi segnali) o un CSV (tempo;segnale)."""
        self._reset_vista = True       # dati nuovi: il prossimo ridisegno mostra tutto il range
        if os.path.isdir(path):
            canali = sorted(f for f in os.listdir(path) if f.lower().endswith('.ch'))
            base = os.path.splitext(os.path.basename(os.path.normpath(path)))[0]
            self._carica_uv(path, base)
            self._carica_strumento(path, base)
            if base in self.spettri:
                # con gli spettri DAD i segnali registrati (.ch) sono solo alcune lunghezze d'onda
                # scelte all'acquisizione: si caricano su richiesta, di default si estrae da zero
                if self.var_includi_ch.get():
                    for c in canali:
                        self._carica_ch(os.path.join(path, c), '%s %s' % (base, os.path.splitext(c)[0]), path, c)
                self._estrazione_iniziale(base, path, canali)
            elif not canali:
                messagebox.showwarning("Apri", "Nessun file .ch ne' .uv nella cartella.")
                return
            else:
                for c in canali:
                    self._carica_ch(os.path.join(path, c), '%s %s' % (base, os.path.splitext(c)[0]), path, c)
        elif path.lower().endswith('.ch'):
            d = os.path.dirname(path)
            base = os.path.splitext(os.path.basename(d))[0] if d.lower().endswith('.d') else ''
            nome = os.path.splitext(os.path.basename(path))[0]
            self._carica_ch(path, (base + ' ' + nome).strip(), d, os.path.basename(path))
        elif path.lower().endswith(('.csv', '.txt')):
            self._carica_csv(path)
        else:
            messagebox.showwarning("Apri", "Formato non supportato: " + os.path.basename(path))

    def _carica_ch(self, path, nome, cartella, canale):
        try:
            t, y, info = self.leggi_ch(path)
        except Exception as e:
            messagebox.showerror("Errore", "Impossibile leggere %s:\n%s" % (os.path.basename(path), e))
            return
        info['Source'] = path
        if cartella.lower().endswith('.d'):
            self._unisci_metadati(info, cartella)
        picchi_cs = self._picchi_chemstation(cartella, canale) if cartella.lower().endswith('.d') else []
        self._aggiungi(nome, t, y, info, picchi_cs)

    def _carica_csv(self, path):
        try:
            df = pd.read_csv(path, sep=None, engine='python', encoding='latin-1', header=None,
                             usecols=[0, 1]).apply(pd.to_numeric, errors='coerce').dropna()
        except Exception as e:
            messagebox.showerror("Errore", "Impossibile leggere %s:\n%s" % (os.path.basename(path), e))
            return
        self._aggiungi(os.path.splitext(os.path.basename(path))[0], df.iloc[:, 0].to_numpy(),
                       df.iloc[:, 1].to_numpy(), {'Source': path, 'Type': 'CSV (time, signal)'})

    def carica_cartella_dialog(self):
        d = filedialog.askdirectory(title="Seleziona la cartella .D di ChemStation")
        if d:
            self.processa_file(d)
            self.aggiorna_vista()

    def carica_file_dialog(self):
        for p in filedialog.askopenfilenames(title="Apri segnale",
                                             filetypes=[("Segnali", "*.ch *.csv *.txt"), ("Tutti", "*.*")]):
            self.processa_file(p)
        self.aggiorna_vista()

    def handle_drop(self, event):
        for p in self.root.tk.splitlist(event.data):
            self.processa_file(p)
        self.aggiorna_vista()

    # ------------------------------------------------------------------ vista
    def _selezionati(self):
        nomi = list(self.cromatogrammi)
        return [nomi[i] for i in self.lista.curselection() if i < len(nomi)]

    def aggiorna_vista(self):
        """Aggiorna lista, tabella e grafico dopo ogni modifica ai dati."""
        sel = set(self._selezionati())
        self.lista.delete(0, tk.END)
        for i, (n, c) in enumerate(self.cromatogrammi.items()):
            self.lista.insert(tk.END, n)
            if c['nascosto']:
                self.lista.itemconfig(i, fg='grey')
            if n in sel:
                self.lista.selection_set(i)
        if len(self.cromatogrammi) == 1:
            self.lista.selection_set(0)
        self._popola_tabella_dati()
        self._ridisegna()
        self._aggiorna_metadati()

    def _ridisegna(self):
        self._cursor_artists = []
        # ax.clear() rimette gli assi in modalita' automatica e fa perdere lo zoom dell'utente a ogni
        # Hide/Show, Update trace, picchi, ecc. Se l'asse non e' automatico (zoom, pan o set_xlim)
        # i suoi limiti si conservano; con dati nuovi (self._reset_vista) si riparte dall'intero range.
        xl, yl = self.ax.get_xlim(), self.ax.get_ylim()
        tieni_x = not self.ax.get_autoscalex_on() and not self._reset_vista
        tieni_y = not self.ax.get_autoscaley_on() and not self._reset_vista
        self._reset_vista = False
        self.ax.clear()
        for n, c in self.cromatogrammi.items():
            if c['nascosto']:
                continue
            df = c['df']
            self.ax.plot(df.index, df['mAU'], lw=1, label=n)
            col = self.ax.lines[-1].get_color()
            for p in c['picchi']:
                m = (df.index >= p['inizio']) & (df.index <= p['fine'])
                x = df.index[m]
                base = np.interp(x, [p['inizio'], p['fine']], [p['y_inizio'], p['y_fine']])
                self.ax.fill_between(x, base, df['mAU'][m], color=col, alpha=0.3, lw=0)
                self.ax.annotate('%.2f' % p['rt'], (p['rt'], p['y_apice']), xytext=(0, 4),
                                 textcoords='offset points', ha='center', fontsize=8, color=col)
            if self.var_picchi_cs.get():
                for p in c['picchi_cs']:
                    y = np.interp(p['rt'], df.index, df['mAU'])
                    self.ax.plot([p['rt']], [y], marker='v', ms=4, color='red', ls='none')
        self._marcatori_artisti = []
        for t, col in self._marcatori_visibili():
            self._disegna_marcatore(t, col)
        self.ax.set_xlabel("Time (min)")
        self.ax.set_ylabel("Signal (mAU)")
        if any(not l.get_label().startswith('_') for l in self.ax.lines):   # almeno una traccia con etichetta
            # posizione fissa: con 'best' matplotlib la ricalcola a ogni redraw tenendo conto anche
            # della linea del cursore, e la legenda scapperebbe quando il cursore ci passa sopra
            self.ax.legend(loc='upper right')
        if origin_style is not None:
            origin_style.applica_stile_origin(self.ax, self.fig, set_size=False)
            self.ax.xaxis.label.set_fontsize(14)
            self.ax.yaxis.label.set_fontsize(14)
            self.ax.tick_params(axis='both', labelsize=12)
            leg = self.ax.get_legend()
            if leg is not None:
                for t in leg.get_texts():
                    t.set_fontsize(11)
            self._ritaglia_margini()
        else:
            self.ax.grid(True, linestyle=':', alpha=0.6)
        if tieni_x:
            self.ax.set_xlim(xl)
        if tieni_y:
            self.ax.set_ylim(yl)
        self.canvas.draw()

    def _ritaglia_margini(self, event=None):
        """Reset dei margini + tight_layout: con font ingranditi il titolo dell'asse Y
        finirebbe fuori dal canvas, e i margini sono frazioni (cambiano col ridimensionamento)."""
        if not self.cromatogrammi or origin_style is None:
            return
        try:
            self.fig.subplots_adjust(left=0.125, right=0.9, bottom=0.11, top=0.88,
                                     wspace=0.2, hspace=0.2)
            self.fig.tight_layout()
        except Exception:
            pass

    def _popola_tabella_dati(self):
        self.tabella.delete(*self.tabella.get_children())
        if not self.cromatogrammi:
            self.tabella['columns'] = ()
            return
        tab = pd.concat([c['df']['mAU'].rename(n) for n, c in self.cromatogrammi.items()], axis=1)
        tab.index = tab.index.round(4)
        tab = tab.sort_index().round(4)
        cols = ['Time (min)'] + list(tab.columns)
        self.tabella['columns'] = cols
        for c in cols:
            self.tabella.heading(c, text=c)
            self.tabella.column(c, width=90 if c == 'Time (min)' else 110, stretch=False)
        for t, riga in zip(tab.index, tab.itertuples(index=False)):
            self.tabella.insert('', tk.END, values=[t] + ['' if pd.isna(v) else v for v in riga])

    def _copia_tabella_dati(self, event=None):
        sel = self.tabella.selection()
        if not sel:
            return 'break'
        righe = ['\t'.join(self.tabella['columns'])]
        righe += ['\t'.join(str(v) for v in self.tabella.item(i, 'values')) for i in sel]
        self.root.clipboard_clear()
        self.root.clipboard_append('\n'.join(righe))
        return 'break'

    def _seleziona_tutto_tabella(self, event=None):
        self.tabella.selection_set(self.tabella.get_children())
        return 'break'

    def _aggiorna_metadati(self):
        self.meta.config(state='normal')
        self.meta.delete('1.0', tk.END)
        for n in self._selezionati():
            c = self.cromatogrammi[n]
            self.meta.insert(tk.END, '== %s\n' % n)
            # le informazioni piu' usate in cima (campione, metodo, colonna, iniezione), poi il resto
            prime = ('Sample', 'Date', 'Method', 'Column', 'Column size', 'Inj Volume', 'Location', 'Signal')
            ordine = [k for k in prime if k in c['info']] + [k for k in c['info'] if k not in prime]
            for k in ordine:
                self.meta.insert(tk.END, '%-19s %s\n' % (k, c['info'][k]))
            df = c['df']
            self.meta.insert(tk.END, '%-19s %d (%.2f - %.2f min)\n' % (
                'Points', len(df), df.index[0], df.index[-1]))
            if c['picchi_cs']:
                self.meta.insert(tk.END, '%-19s %d\n' % ('ChemStation peaks', len(c['picchi_cs'])))
        self.meta.config(state='disabled')
        # il nome del campione (o dei campioni) selezionato compare anche nel titolo della finestra
        campioni = []
        for n in self._selezionati():
            c = self.cromatogrammi[n]['info'].get('Sample')
            if c and c not in campioni:
                campioni.append(c)
        self.root.title("HPLC Manager" + (" - " + ", ".join(campioni) if campioni else ""))

    # ------------------------------------------------------------------ cursore
    def on_mouse_move(self, event):
        self._pulisci_cursore()
        if event.inaxes is not self.ax or event.xdata is None:
            return
        vis = [(n, c['df']) for n, c in self.cromatogrammi.items() if not c['nascosto']]
        if not vis:
            return
        x = event.xdata
        righe = ['t = %.3f min' % x]
        for n, df in vis:
            righe.append('%s: %.2f mAU' % (n, np.interp(x, df.index, df['mAU'])))
        # il riquadro dei valori sta in alto a destra sotto la legenda, ancorato al suo bordo
        # destro: la posizione non dipende dal cursore, cambia solo il contenuto
        y_alto = 0.98
        leg = self.ax.get_legend()
        if leg is not None:
            bb = leg.get_window_extent(self.canvas.get_renderer())
            y_alto = self.ax.transAxes.inverted().transform((bb.x0, bb.y0))[1] - 0.015
        self._cursor_artists.append(self.ax.axvline(x, color='grey', lw=0.8, ls=':', zorder=1))
        self._cursor_artists.append(self.ax.text(
            0.98, y_alto, '\n'.join(righe), transform=self.ax.transAxes, va='top', ha='right',
            fontsize=9, zorder=6, bbox=dict(boxstyle='round', fc='white', ec='grey', alpha=0.85)))
        self.canvas.draw_idle()

    def on_mouse_leave(self, event):
        self._pulisci_cursore()
        self.canvas.draw_idle()

    def _pulisci_cursore(self):
        for a in self._cursor_artists:
            try:
                a.remove()
            except Exception:
                pass
        self._cursor_artists = []

    # ------------------------------------------------------------------ gestione lista
    def remove_selected(self):
        for n in self._selezionati():
            del self.cromatogrammi[n]
        self.aggiorna_vista()

    def clear_all(self):
        self._reset_vista = True
        self.cromatogrammi.clear()
        self.spettri.clear()
        self.strumento.clear()
        self._strum_chiudi()
        self._sp_chiudi_finestra()
        self._sp_live.clear()
        self._dirty = False
        for w in self.f_tool_host.winfo_children():
            w.destroy()
        self.aggiorna_vista()

    def hide_selected(self):
        for n in self._selezionati():
            self.cromatogrammi[n]['nascosto'] = True
        self.aggiorna_vista()

    def show_selected(self):
        for n in self._selezionati():
            self.cromatogrammi[n]['nascosto'] = False
        self.aggiorna_vista()

    def esporta_csv(self):
        if not self.cromatogrammi:
            return
        p = filedialog.asksaveasfilename(defaultextension='.csv', filetypes=[("CSV", "*.csv")])
        if not p:
            return
        tab = pd.concat([c['df']['mAU'].rename(n) for n, c in self.cromatogrammi.items()], axis=1)
        tab.sort_index().to_csv(p, sep=';', na_rep='', encoding='latin-1')
        self._dirty = False

    # ------------------------------------------------------------------ tool: picchi
    def _una_selezione(self, titolo):
        sel = self._selezionati()
        if len(sel) != 1:
            messagebox.showinfo(titolo, "Seleziona un solo cromatogramma nella lista.")
            return None
        return sel[0]

    def apri_picchi(self):
        f = self._apri_pannello_tool("Peaks")
        if not HAS_SCIPY:
            tk.Label(f, text="scipy non installato: rilevamento picchi non disponibile.").pack()
            return
        self.var_prom = tk.StringVar(value='5')
        self.var_tmin = tk.StringVar(value='')
        self.var_tmax = tk.StringVar(value='')
        for r, (txt, var) in enumerate((("Min prominence (mAU):", self.var_prom),
                                        ("From (min):", self.var_tmin), ("To (min):", self.var_tmax))):
            tk.Label(f, text=txt).grid(row=r, column=0, sticky='w')
            tk.Entry(f, textvariable=var, width=8).grid(row=r, column=1, padx=4)
        tk.Button(f, text="Detect & integrate", command=self._pk_applica).grid(row=3, column=0, columnspan=2, pady=3)
        tk.Button(f, text="Clear peaks", command=self._pk_pulisci).grid(row=3, column=2)
        cols = ('n', 'rt', 'h', 'area', 'pct')
        self.tab_picchi = ttk.Treeview(f, columns=cols, show='headings', height=12)
        for c, t, w in (('n', '#', 30), ('rt', 'RT (min)', 70), ('h', 'Height', 65),
                        ('area', 'Area (mAU*s)', 90), ('pct', 'Area %', 60)):
            self.tab_picchi.heading(c, text=t)
            self.tab_picchi.column(c, width=w, anchor='e')
        self.tab_picchi.grid(row=4, column=0, columnspan=3, sticky='ew', pady=3)
        tk.Button(f, text="Copy table", command=self._pk_copia).grid(row=5, column=0, sticky='w')
        tk.Button(f, text="Export CSV", command=self.esporta_picchi).grid(row=5, column=1, columnspan=2, sticky='w')
        tk.Label(f, text="Baseline: retta fra le basi dei picchi.\nTriangoli rossi: picchi ChemStation.",
                 justify='left', fg='grey30').grid(row=6, column=0, columnspan=3, sticky='w')
        sel = self._selezionati()
        if len(sel) == 1:
            self._pk_riempi(sel[0])

    def _pk_applica(self):
        n = self._una_selezione("Peaks")
        if n is None:
            return
        try:
            prom = float(self.var_prom.get().replace(',', '.'))
            tmin = float(self.var_tmin.get().replace(',', '.')) if self.var_tmin.get().strip() else None
            tmax = float(self.var_tmax.get().replace(',', '.')) if self.var_tmax.get().strip() else None
        except ValueError:
            messagebox.showerror("Peaks", "Valori numerici non validi.")
            return
        df = self.cromatogrammi[n]['df']
        t, y = df.index.to_numpy(), df['mAU'].to_numpy()
        idx, _ = find_peaks(y, prominence=prom)
        # basi del picco = minimi locali piu' vicini all'apice (su una copia leggermente
        # smussata, per non fermarsi sul rumore). Le 'left_bases' di scipy possono arrivare
        # fino a inizio cromatogramma per un picco alto e non vanno bene per l'integrazione.
        ys = savgol_filter(y, 7, 2) if len(y) > 7 else y
        picchi = []
        for i in idx:
            if (tmin is not None and t[i] < tmin) or (tmax is not None and t[i] > tmax):
                continue
            j = max(i - 3, 0) + int(np.argmax(ys[max(i - 3, 0):i + 4]))   # apice della curva smussata
            l = j
            while l > 0 and ys[l - 1] <= ys[l]:
                l -= 1
            r = j
            while r < len(y) - 1 and ys[r + 1] <= ys[r]:
                r += 1
            tt, yy = t[l:r + 1], y[l:r + 1]
            base = np.linspace(yy[0], yy[-1], len(yy))
            area = float(np.trapezoid(yy - base, tt)) * 60.0   # mAU*min -> mAU*s
            # RT al vertice della parabola per i 3 punti attorno all'apice
            rt = t[i]
            if 0 < i < len(t) - 1:
                a, b_, c = y[i - 1], y[i], y[i + 1]
                den = a - 2 * b_ + c
                if den != 0:
                    rt = t[i] + 0.5 * (a - c) / den * (t[i + 1] - t[i])
            picchi.append({'rt': float(rt), 'altezza': float(y[i] - base[i - l]), 'area': area,
                           'inizio': float(tt[0]), 'fine': float(tt[-1]), 'y_inizio': float(yy[0]),
                           'y_fine': float(yy[-1]), 'y_apice': float(y[i])})
        tot = sum(p['area'] for p in picchi) or 1.0
        for k, p in enumerate(picchi, 1):
            p['n'] = k
            p['pct'] = 100.0 * p['area'] / tot
        self.cromatogrammi[n]['picchi'] = picchi
        self._pk_riempi(n)
        self._ridisegna()

    def _pk_riempi(self, n):
        self.tab_picchi.delete(*self.tab_picchi.get_children())
        for p in self.cromatogrammi[n]['picchi']:
            self.tab_picchi.insert('', tk.END, values=(
                p['n'], '%.3f' % p['rt'], '%.2f' % p['altezza'], '%.2f' % p['area'], '%.2f' % p['pct']))

    def _pk_pulisci(self):
        n = self._una_selezione("Peaks")
        if n is not None:
            self.cromatogrammi[n]['picchi'] = []
            self._pk_riempi(n)
            self._ridisegna()

    def _pk_copia(self):
        righe = ['#\tRT (min)\tHeight\tArea (mAU*s)\tArea %']
        righe += ['\t'.join(str(v) for v in self.tab_picchi.item(i, 'values'))
                  for i in self.tab_picchi.get_children()]
        self.root.clipboard_clear()
        self.root.clipboard_append('\n'.join(righe))

    # ------------------------------------------------------------------ tool: smoothing, trim, normalize
    def apri_smoothing(self):
        f = self._apri_pannello_tool("Smoothing (Savitzky-Golay)")
        self.var_sg_w = tk.StringVar(value='11')
        self.var_sg_p = tk.StringVar(value='3')
        tk.Label(f, text="Window (points, odd):").grid(row=0, column=0, sticky='w')
        tk.Entry(f, textvariable=self.var_sg_w, width=6).grid(row=0, column=1)
        tk.Label(f, text="Polynomial order:").grid(row=1, column=0, sticky='w')
        tk.Entry(f, textvariable=self.var_sg_p, width=6).grid(row=1, column=1)
        tk.Button(f, text="Apply (new trace)", command=self._sg_applica).grid(row=2, column=0, columnspan=2, pady=3)

    def _sg_applica(self):
        n = self._una_selezione("Smoothing")
        if n is None:
            return
        try:
            w, p = int(self.var_sg_w.get()), int(self.var_sg_p.get())
        except ValueError:
            messagebox.showerror("Smoothing", "Valori non validi.")
            return
        if w % 2 == 0 or w <= p:
            messagebox.showerror("Smoothing", "La finestra deve essere dispari e maggiore dell'ordine.")
            return
        df = self.cromatogrammi[n]['df']
        if HAS_SCIPY:
            y = savgol_filter(df['mAU'].to_numpy(), w, p)
        else:
            y = np.convolve(df['mAU'].to_numpy(), np.ones(w) / w, mode='same')
        info = dict(self.cromatogrammi[n]['info'], Processing='Smoothing w=%d p=%d' % (w, p))
        self._aggiungi('%s_sg' % n, df.index.to_numpy(), y, info)
        self._dirty = True
        self.aggiorna_vista()

    def apri_trim(self):
        f = self._apri_pannello_tool("Trim")
        self.var_tr_a = tk.StringVar()
        self.var_tr_b = tk.StringVar()
        sel = self._selezionati()
        if len(sel) == 1:
            ix = self.cromatogrammi[sel[0]]['df'].index
            self.var_tr_a.set('%.2f' % ix[0])
            self.var_tr_b.set('%.2f' % ix[-1])
        tk.Label(f, text="From (min):").grid(row=0, column=0, sticky='w')
        tk.Entry(f, textvariable=self.var_tr_a, width=8).grid(row=0, column=1)
        tk.Label(f, text="To (min):").grid(row=1, column=0, sticky='w')
        tk.Entry(f, textvariable=self.var_tr_b, width=8).grid(row=1, column=1)
        tk.Button(f, text="Apply (new trace)", command=self._tr_applica).grid(row=2, column=0, columnspan=2, pady=3)

    def _tr_applica(self):
        n = self._una_selezione("Trim")
        if n is None:
            return
        try:
            a, b = float(self.var_tr_a.get().replace(',', '.')), float(self.var_tr_b.get().replace(',', '.'))
        except ValueError:
            messagebox.showerror("Trim", "Valori non validi.")
            return
        df = self.cromatogrammi[n]['df']
        sub = df[(df.index >= min(a, b)) & (df.index <= max(a, b))]
        if len(sub) < 2:
            messagebox.showerror("Trim", "Intervallo vuoto.")
            return
        info = dict(self.cromatogrammi[n]['info'], Processing='Trim %.2f-%.2f min' % (a, b))
        self._aggiungi('%s_trim' % n, sub.index.to_numpy(), sub['mAU'].to_numpy(), info)
        self._dirty = True
        self.aggiorna_vista()

    def normalizza(self):
        sel = self._selezionati()
        if not sel:
            messagebox.showinfo("Normalize", "Seleziona almeno un cromatogramma.")
            return
        for n in sel:
            df = self.cromatogrammi[n]['df']
            m = float(np.abs(df['mAU']).max()) or 1.0
            info = dict(self.cromatogrammi[n]['info'], Processing='Normalized to max = 100')
            self._aggiungi('%s_norm' % n, df.index.to_numpy(), df['mAU'].to_numpy() * 100.0 / m, info)
        self._dirty = True
        self.aggiorna_vista()

    # ------------------------------------------------------------------ tool: spettri DAD
    def _estrazione_iniziale(self, base, cartella, canali):
        """Dopo l'apertura di una .D con spettri: estrae una traccia (lunghezza d'onda del primo
        segnale registrato, banda 4 nm, riferimento spento) e apre il pannello di controllo."""
        lam = 254.0
        for c in canali[:1]:
            try:
                sig = self.leggi_ch(os.path.join(cartella, c))[2].get('Signal', '')
                m = re.search(r'Sig=\s*([\d.]+)', sig)
                if m:
                    lam = float(m.group(1))
            except Exception:
                pass
        self.var_sp_set = tk.StringVar(value=base)
        self._sp_valori_iniziali(lam)
        self.apri_spettri()
        self._sp_estrai()

    def _sp_valori_iniziali(self, lam):
        self.var_sp_l = tk.StringVar(value='%g' % lam)
        self.var_sp_b = tk.StringVar(value='4')
        self.var_sp_rl = tk.StringVar(value='360')
        self.var_sp_rb = tk.StringVar(value='100')
        self.var_sp_uso_rif = tk.BooleanVar(value=False)     # riferimento spento di default
        self.var_sp_t = tk.StringVar(value='4.35')

    def apri_spettri(self):
        f = self._apri_pannello_tool("Detector signal (DAD)")
        if not self.spettri:
            tk.Label(f, text="Nessuno spettro caricato: apri una cartella .D\ncon il file .uv.",
                     justify='left').pack(anchor='w')
            return
        nomi = list(self.spettri)
        if not hasattr(self, 'var_sp_l'):
            self.var_sp_set = tk.StringVar(value=nomi[0])
            self._sp_valori_iniziali(254.0)
        if self.var_sp_set.get() not in self.spettri:
            self.var_sp_set.set(nomi[0])
        cb = ttk.Combobox(f, textvariable=self.var_sp_set, values=nomi, state='readonly', width=18)
        cb.grid(row=0, column=0, columnspan=3, sticky='w')
        d = self.spettri[self.var_sp_set.get()]
        tk.Label(f, text="%s" % d['info']['Spectra'], fg='grey30').grid(row=1, column=0, columnspan=3, sticky='w')
        tk.Label(f, text="Wavelength (nm)").grid(row=2, column=0, sticky='w')
        e_l = tk.Entry(f, textvariable=self.var_sp_l, width=7)
        e_l.grid(row=2, column=1)
        tk.Label(f, text="Bandwidth (nm)").grid(row=3, column=0, sticky='w')
        e_b = tk.Entry(f, textvariable=self.var_sp_b, width=7)
        e_b.grid(row=3, column=1)
        self.chk_rif = tk.Checkbutton(f, text="Reference (nm / bw)", variable=self.var_sp_uso_rif,
                                      command=self._sp_stato_rif)
        self.chk_rif.grid(row=4, column=0, sticky='w')
        self.e_rl = tk.Entry(f, textvariable=self.var_sp_rl, width=7)
        self.e_rl.grid(row=4, column=1)
        self.e_rb = tk.Entry(f, textvariable=self.var_sp_rb, width=7)
        self.e_rb.grid(row=4, column=2)
        self._sp_stato_rif()
        tk.Button(f, text="Update trace", command=lambda: self._sp_estrai(False)).grid(
            row=5, column=0, columnspan=2, pady=3, sticky='ew')
        tk.Button(f, text="Add as new", command=lambda: self._sp_estrai(True)).grid(
            row=5, column=2, pady=3)
        for e in (e_l, e_b, self.e_rl, self.e_rb):
            e.bind('<Return>', lambda ev: self._sp_estrai(False))
        tk.Button(f, text="Several wavelengths...", command=self.apri_estrazione_multipla).grid(
            row=6, column=0, columnspan=3, pady=1, sticky='ew')
        tk.Label(f, text="Spectrum at (min)").grid(row=7, column=0, sticky='w')
        tk.Entry(f, textvariable=self.var_sp_t, width=7).grid(row=7, column=1)
        tk.Button(f, text="Show", command=self._sp_mostra).grid(row=7, column=2)
        tk.Checkbutton(f, text="Click on plot shows the spectrum\n(Shift+click works even if unticked)", justify="left",
                       variable=self.var_click_spec).grid(row=8, column=0, columnspan=3, sticky='w')
        tk.Button(f, text="Time-wavelength map", command=self._sp_mappa).grid(
            row=9, column=0, columnspan=3, pady=3)
        tk.Button(f, text="Export picked spectra (CSV)", command=self.esporta_spettri).grid(
            row=10, column=0, columnspan=3, pady=1)

    def _sp_stato_rif(self):
        stato = 'normal' if self.var_sp_uso_rif.get() else 'disabled'
        self.e_rl.config(state=stato)
        self.e_rb.config(state=stato)

    def _sp_dataset(self):
        return self.spettri.get(self.var_sp_set.get()) if hasattr(self, 'var_sp_set') else None

    def _sp_estrai(self, nuova=False):
        """Estrae il cromatogramma alla lunghezza d'onda/banda/riferimento del pannello.

        Di default modifica la traccia estratta selezionata nella lista (la sostituisce nella stessa
        posizione); se non ce n'e', l'ultima estratta; con nuova=True ne aggiunge una separata, che
        diventa la selezionata e quindi quella che l'Update successivo modifica."""
        d = self._sp_dataset()
        if d is None:
            return
        try:
            l, b = float(self.var_sp_l.get().replace(',', '.')), float(self.var_sp_b.get().replace(',', '.'))
            if b <= 0:
                raise ValueError
            rif = rb = None
            if self.var_sp_uso_rif.get():
                rif, rb = float(self.var_sp_rl.get().replace(',', '.')), float(self.var_sp_rb.get().replace(',', '.'))
                if rb <= 0:
                    raise ValueError
        except ValueError:
            messagebox.showerror("Spectra", "Valori numerici non validi.")
            return
        try:
            y = self.estrai_lambda(d['wl'], d['S'], l, b, rif, rb)
        except ValueError as e:
            messagebox.showerror("Spectra", str(e))
            return
        ds = self.var_sp_set.get()
        etichetta = self._etichetta_estrazione(l, b, rif, rb)
        info = self._info_estrazione(d, ds, etichetta)
        estr = {'ds': ds, 'l': l, 'b': b, 'rif': rif, 'rb': rb}
        vecchio = None if nuova else self._sp_bersaglio(ds)
        if vecchio is None:
            nome = self._aggiungi('%s %s' % (ds, etichetta), d['t'], y, info, estr=estr)
        else:
            # sostituisce la traccia selezionata mantenendo la sua posizione nella lista
            nome = '%s %s' % (ds, etichetta)
            if nome != vecchio:
                nome = self._nome_libero(nome)
            df = pd.DataFrame({'mAU': y}, index=pd.Index(d['t'], name='Time (min)'))
            nuovo = dict(self.cromatogrammi[vecchio], df=df, info=info, picchi=[], picchi_cs=[],
                         nascosto=False, estr=estr)
            self.cromatogrammi = {(nome if k == vecchio else k): (nuovo if k == vecchio else v)
                                  for k, v in self.cromatogrammi.items()}
        self._sp_live[ds] = nome          # l'ultima toccata: serve se non c'e' una selezione utile
        if nuova:
            self._dirty = True            # una traccia 'viva' si ricalcola sempre: non serve esportarla
        self.aggiorna_vista()
        self._seleziona_nome(nome)        # la traccia appena estratta/modificata diventa quella attiva

    def _sp_bersaglio(self, ds):
        """Traccia che 'Update trace' deve modificare: quella selezionata nella lista se e' una traccia
        estratta dal dataset `ds`; altrimenti l'ultima estratta/modificata; None se non ce n'e'."""
        sel = self._selezionati()
        if len(sel) == 1:
            e = self.cromatogrammi[sel[0]].get('estr')
            if e and e['ds'] == ds:
                return sel[0]
        ultimo = self._sp_live.get(ds)
        return ultimo if ultimo in self.cromatogrammi else None

    def _seleziona_nome(self, nome):
        """Seleziona nella lista la traccia `nome` (e aggiorna i metadati)."""
        nomi = list(self.cromatogrammi)
        if nome not in nomi:
            return
        i = nomi.index(nome)
        self.lista.selection_clear(0, tk.END)
        self.lista.selection_set(i)
        self.lista.see(i)
        self._aggiorna_metadati()

    def _su_selezione(self, event=None):
        """Cambio di selezione nella lista: metadati, e per una traccia estratta i suoi parametri
        (lunghezza d'onda, banda, riferimento) tornano nei campi del pannello, pronti da modificare."""
        self._aggiorna_metadati()
        sel = self._selezionati()
        if len(sel) != 1 or not hasattr(self, 'var_sp_l'):
            return
        e = self.cromatogrammi[sel[0]].get('estr')
        if not e or e['ds'] not in self.spettri:
            return
        self.var_sp_set.set(e['ds'])
        self.var_sp_l.set('%g' % e['l'])
        self.var_sp_b.set('%g' % e['b'])
        self.var_sp_uso_rif.set(e['rif'] is not None)
        if e['rif'] is not None:
            self.var_sp_rl.set('%g' % e['rif'])
            self.var_sp_rb.set('%g' % e['rb'])
        try:
            self._sp_stato_rif()
        except tk.TclError:           # pannello non (piu') costruito
            pass

    def _disegna_marcatore(self, t, col):
        """Linea tratteggiata + etichetta del tempo sul cromatogramma (senza ridisegnare il resto,
        cosi' lo zoom dell'utente non si perde)."""
        self._marcatori_artisti.append(self.ax.axvline(t, color=col, lw=0.9, ls='--', zorder=1, gid='marcatore'))
        self._marcatori_artisti.append(self.ax.annotate(
            '%.2f' % t, (t, 1.0), xycoords=('data', 'axes fraction'), xytext=(2, -3),
            textcoords='offset points', rotation=90, va='top', ha='left', fontsize=8, color=col, gid='marcatore'))

    @staticmethod
    def _etichetta_estrazione(l, b, rif, rb):
        return '%g,%g' % (l, b) + (' ref %g,%g' % (rif, rb) if rif is not None else ' no ref')

    @staticmethod
    def _info_estrazione(d, ds, etichetta):
        info = {k: v for k, v in d['info'].items() if k not in ('Version', 'Type', 'Spectra')}
        info.update({'Source': ds + ' (dad .uv)', 'Signal': 'Extracted Sig=' + etichetta})
        return info

    def estrai_serie(self, ds, specifiche):
        """Estrae piu' tracce dal dataset `ds` in un colpo e le aggiunge alla lista.

        `specifiche`: lista di (lunghezza d'onda, banda, riferimento o None, banda del riferimento).
        Tutto o niente: se una riga non e' valida (banda non positiva, fuori dallo spettro acquisito)
        solleva ValueError con il numero della riga e non aggiunge nulla. Ritorna i nomi delle tracce."""
        d = self.spettri[ds]
        calcolate = []
        for i, (l, b, rif, rb) in enumerate(specifiche, 1):
            if b <= 0 or (rif is not None and rb <= 0):
                raise ValueError("Riga %d: la banda deve essere positiva." % i)
            try:
                y = self.estrai_lambda(d['wl'], d['S'], l, b, rif, rb)
            except ValueError as e:
                raise ValueError("Riga %d: %s" % (i, e))
            calcolate.append((self._etichetta_estrazione(l, b, rif, rb), y,
                              {'ds': ds, 'l': l, 'b': b, 'rif': rif, 'rb': rb}))
        nomi = []
        for etichetta, y, estr in calcolate:
            nomi.append(self._aggiungi('%s %s' % (ds, etichetta), d['t'], y,
                                       self._info_estrazione(d, ds, etichetta), estr=estr))
        return nomi

    def apri_estrazione_multipla(self):
        """Dialogo: quante tracce vuoi, e per ognuna lunghezza d'onda, banda ed eventuale
        riferimento (ognuna col suo, o nessuno). Le tracce create restano nella lista come
        'Add as new'; la traccia 'viva' del pannello non cambia."""
        if not self.spettri:
            messagebox.showinfo("Several wavelengths", "Nessuno spettro caricato: apri una cartella .D con il .uv.")
            return
        ds = self.var_sp_set.get() if hasattr(self, 'var_sp_set') and self.var_sp_set.get() in self.spettri \
            else next(iter(self.spettri))
        top = tk.Toplevel(self.root)
        top.title("Extract several wavelengths")
        top.transient(self.root)
        tk.Label(top, text="Dataset: %s" % ds, fg='grey30').grid(row=0, column=0, columnspan=6, sticky='w', padx=8, pady=(8, 0))
        tk.Label(top, text="Number of traces:").grid(row=1, column=0, columnspan=2, sticky='w', padx=8, pady=6)
        var_n = tk.IntVar(value=len(self._sp_multi_ultimo) or 3)
        for c, t in enumerate(("#", "Wavelength (nm)", "Bandwidth (nm)", "Reference", "Ref. wavelength", "Ref. bandwidth")):
            tk.Label(top, text=t, font=('Segoe UI', 9, 'bold')).grid(row=2, column=c, padx=4)
        righe = []
        base_b = self.var_sp_b.get() if hasattr(self, 'var_sp_b') else '4'
        base_rl = self.var_sp_rl.get() if hasattr(self, 'var_sp_rl') else '360'
        base_rb = self.var_sp_rb.get() if hasattr(self, 'var_sp_rb') else '100'
        base_usa = self.var_sp_uso_rif.get() if hasattr(self, 'var_sp_uso_rif') else False

        def aggiorna_stato(r):
            st = 'normal' if r['usa'].get() else 'disabled'
            r['e_rl'].config(state=st)
            r['e_rb'].config(state=st)

        def aggiungi_riga():
            i = len(righe)
            prec = self._sp_multi_ultimo[i] if i < len(self._sp_multi_ultimo) else None
            r = {'l': tk.StringVar(value=prec['l'] if prec else ''),
                 'b': tk.StringVar(value=prec['b'] if prec else base_b),
                 'usa': tk.BooleanVar(value=prec['usa'] if prec else base_usa),
                 'rl': tk.StringVar(value=prec['rl'] if prec else base_rl),
                 'rb': tk.StringVar(value=prec['rb'] if prec else base_rb)}
            r['w'] = [tk.Label(top, text=str(i + 1)),
                      tk.Entry(top, textvariable=r['l'], width=10),
                      tk.Entry(top, textvariable=r['b'], width=10),
                      tk.Checkbutton(top, variable=r['usa'], command=lambda r=r: aggiorna_stato(r))]
            r['e_rl'] = tk.Entry(top, textvariable=r['rl'], width=10)
            r['e_rb'] = tk.Entry(top, textvariable=r['rb'], width=10)
            r['w'] += [r['e_rl'], r['e_rb']]
            for c, w in enumerate(r['w']):
                w.grid(row=3 + i, column=c, padx=4, pady=1)
            aggiorna_stato(r)
            righe.append(r)

        def imposta_n(*_):
            try:
                n = max(1, min(30, int(var_n.get())))
            except (tk.TclError, ValueError):
                return
            while len(righe) < n:
                aggiungi_riga()
            while len(righe) > n:
                for w in righe.pop()['w']:
                    w.destroy()
            top.update_idletasks()

        spin = tk.Spinbox(top, from_=1, to=30, width=4, textvariable=var_n, command=imposta_n)
        spin.grid(row=1, column=2, sticky='w')
        spin.bind('<Return>', imposta_n)
        spin.bind('<FocusOut>', imposta_n)
        frame_btn = tk.Frame(top)
        frame_btn.grid(row=100, column=0, columnspan=6, pady=8)

        def leggi():
            imposta_n()
            specifiche, salvate = [], []
            for i, r in enumerate(righe, 1):
                try:
                    l, b = float(r['l'].get().replace(',', '.')), float(r['b'].get().replace(',', '.'))
                    rif = rb = None
                    if r['usa'].get():
                        rif, rb = float(r['rl'].get().replace(',', '.')), float(r['rb'].get().replace(',', '.'))
                except ValueError:
                    raise ValueError("Riga %d: valori numerici non validi." % i)
                specifiche.append((l, b, rif, rb))
                salvate.append({k: (r[k].get()) for k in ('l', 'b', 'usa', 'rl', 'rb')})
            return specifiche, salvate

        def estrai():
            try:
                specifiche, salvate = leggi()
                self.estrai_serie(ds, specifiche)
            except ValueError as e:
                messagebox.showerror("Several wavelengths", str(e), parent=top)
                return
            self._sp_multi_ultimo = salvate
            self._dirty = True
            top.destroy()
            self.aggiorna_vista()

        tk.Button(frame_btn, text="Extract", width=12, command=estrai).pack(side=tk.LEFT, padx=6)
        tk.Button(frame_btn, text="Cancel", width=12, command=top.destroy).pack(side=tk.LEFT, padx=6)
        top.bind('<Escape>', lambda e: top.destroy())
        imposta_n()
        top._estrai = estrai          # per i test
        top._righe = righe
        top._var_n = var_n
        top._imposta_n = imposta_n
        return top

    def _marcatori_visibili(self):
        return [(k['t'], k['colore']) for k in self._spettri_scelti if k.get('visibile', True)]

    def _aggiorna_marcatori(self):
        """Ridisegna le linee tratteggiate dei soli spettri selezionati, senza toccare lo zoom."""
        for a in self._marcatori_artisti:
            try:
                a.remove()
            except Exception:
                pass
        self._marcatori_artisti = []
        for t, col in self._marcatori_visibili():
            self._disegna_marcatore(t, col)
        self.canvas.draw_idle()

    def _sp_chiudi_finestra(self):
        """Chiude la finestra degli spettri (se aperta) e azzera gli spettri scelti."""
        w = self._spec_win
        if w is not None:
            try:
                w.destroy()
            except Exception:
                pass
        self._spettri_scelti = []
        self._spec_win = None

    def _sp_finestra(self):
        """Finestra (Toplevel) con il grafico degli spettri e, a destra, l'elenco di quelli scelti
        (spuntati = visibili); creata al primo uso."""
        if self._spec_win is not None and self._spec_win.winfo_exists():
            return self._spec_win
        w = tk.Toplevel(self.root)
        w.title("DAD spectra")
        w.geometry("920x500")
        # elenco a destra (larghezza fissa, scrollabile), grafico a sinistra
        w.f_dx = tk.Frame(w, width=200)
        w.f_dx.pack(side=tk.RIGHT, fill=tk.Y)
        w.f_dx.pack_propagate(False)
        tk.Label(w.f_dx, text="Picked spectra", font=('Segoe UI', 10, 'bold')).pack(anchor='w', padx=6, pady=(6, 0))
        tk.Label(w.f_dx, text="untick = hide, \u2715 = remove", fg='grey30').pack(anchor='w', padx=6)
        tk.Button(w.f_dx, text="Clear all", command=self._sp_svuota).pack(side=tk.BOTTOM, pady=6)
        cv = tk.Canvas(w.f_dx, highlightthickness=0)
        sb = ttk.Scrollbar(w.f_dx, orient='vertical', command=cv.yview)
        cv.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        cv.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(4, 0))
        w.f_righe = tk.Frame(cv)
        cv.create_window((0, 0), window=w.f_righe, anchor='nw')
        w.f_righe.bind('<Configure>', lambda e: cv.configure(scrollregion=cv.bbox('all')))
        w.f_sx = tk.Frame(w)
        w.f_sx.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        w.fig = Figure()
        w.ax = w.fig.add_subplot(111)
        w.canvas = FigureCanvasTkAgg(w.fig, master=w.f_sx)
        NavigationToolbar2Tk(w.canvas, w.f_sx).update()
        w.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        w.ax.set_xlabel("Wavelength (nm)")
        w.ax.set_ylabel("Absorbance (mAU)")
        w.bind('<Destroy>', lambda e: self._sp_finestra_chiusa(e, w))
        self._spec_win = w
        return w

    def _sp_finestra_chiusa(self, event, w):
        """Alla chiusura della finestra degli spettri spariscono anche le linee sul cromatogramma."""
        if event.widget is not w:
            return
        self._spettri_scelti = []
        self._spec_win = None
        try:
            self._aggiorna_marcatori()
        except Exception:
            pass

    def _sp_aggiorna_legenda(self):
        w = self._spec_win
        if w is None:
            return
        visibili = [k['linea'] for k in self._spettri_scelti if k.get('visibile', True)]
        if visibili:
            w.ax.legend(handles=visibili, loc='upper right')
        elif w.ax.get_legend() is not None:
            w.ax.get_legend().remove()
        w.ax.relim(visible_only=True)
        w.ax.autoscale_view()       # rispetta lo zoom: agisce solo sugli assi ancora automatici
        w.canvas.draw_idle()

    def _sp_riga_elenco(self, w, k):
        """Aggiunge all'elenco della finestra la riga dello spettro `k` (spunta + rimuovi)."""
        riga = tk.Frame(w.f_righe)
        riga.pack(fill=tk.X, anchor='w')
        var = tk.BooleanVar(value=k.get('visibile', True))
        tk.Checkbutton(riga, text='%.3f min' % k['t'], variable=var, fg=k['colore'], anchor='w',
                       font=('Segoe UI', 9, 'bold'), command=lambda: self._sp_spunta(k, var.get())
                       ).pack(side=tk.LEFT)
        tk.Button(riga, text='\u2715', width=2, relief=tk.FLAT,
                  command=lambda: self._sp_rimuovi(k)).pack(side=tk.RIGHT, padx=2)
        k['riga'] = riga
        k['var'] = var

    def _sp_spunta(self, k, visibile):
        k['visibile'] = bool(visibile)
        k['linea'].set_visible(bool(visibile))
        self._sp_aggiorna_legenda()
        self._aggiorna_marcatori()

    def _sp_rimuovi(self, k):
        if k not in self._spettri_scelti:
            return
        self._spettri_scelti.remove(k)
        try:
            k['linea'].remove()
            k['riga'].destroy()
        except Exception:
            pass
        self._sp_aggiorna_legenda()
        self._aggiorna_marcatori()

    def _sp_svuota(self):
        for k in list(self._spettri_scelti):
            self._sp_rimuovi(k)

    def _sp_spettro_a(self, t_min):
        """Aggiunge alla finestra degli spettri lo spettro piu' vicino a t_min."""
        d = self._sp_dataset()
        if d is None:
            return
        k = int(np.argmin(np.abs(d['t'] - t_min)))
        w = self._sp_finestra()
        colori = matplotlib.rcParams['axes.prop_cycle'].by_key()['color']
        colore = colori[self._sp_colore_n % len(colori)]
        self._sp_colore_n += 1
        t = float(d['t'][k])
        linea, = w.ax.plot(d['wl'], d['S'][k], lw=1, color=colore, label='%.3f min' % t)
        scelto = {'dataset': self.var_sp_set.get(), 't': t, 'wl': d['wl'].copy(), 'y': d['S'][k].copy(),
                  'colore': colore, 'visibile': True, 'linea': linea}
        self._spettri_scelti.append(scelto)
        self._sp_riga_elenco(w, scelto)
        self._sp_aggiorna_legenda()
        self._aggiorna_marcatori()       # linea tratteggiata sul cromatogramma, dello stesso colore
        w.lift()

    def _sp_mostra(self):
        try:
            self._sp_spettro_a(float(self.var_sp_t.get().replace(',', '.')))
        except ValueError:
            messagebox.showerror("Spectra", "Tempo non valido.")

    def _premi_spettro(self, event):
        # si ricorda dove e' stato premuto il tasto: serve a distinguere un click da un trascinamento
        self._press_pos = (event.x, event.y) if event.button == 1 else None

    def _click_spettro(self, event):
        """Click sul cromatogramma = spettro a quel tempo.

        Si decide al rilascio del tasto: e' un click se il mouse si e' mosso meno di 5 pixel. Cosi'
        funziona anche con la lente (zoom) o il pan della barra attivi: quegli strumenti lavorano
        solo trascinando, e un trascinamento (zoom su un rettangolo, spostamento) non e' un click.
        Maiusc+click vale anche con la casella spenta."""
        press, self._press_pos = getattr(self, '_press_pos', None), None
        if (press is None or event.inaxes is not self.ax or event.button != 1
                or event.xdata is None or not self.spettri):
            return
        if np.hypot(event.x - press[0], event.y - press[1]) >= 5:
            return
        maiusc = event.key is not None and 'shift' in str(event.key)
        if not maiusc and not self.var_click_spec.get():
            return
        self.var_sp_t.set('%.3f' % event.xdata)
        self._sp_spettro_a(event.xdata)

    def _sp_mappa(self):
        d = self._sp_dataset()
        if d is None:
            return
        w = tk.Toplevel(self.root)
        w.title("Time-wavelength map")
        w.geometry("800x520")
        fig = Figure()
        ax = fig.add_subplot(111)
        lim = float(np.percentile(np.abs(d['S']), 99))
        im = ax.imshow(d['S'].T, aspect='auto', origin='lower', cmap='viridis', vmin=0, vmax=lim,
                       extent=[d['t'][0], d['t'][-1], d['wl'][0], d['wl'][-1]])
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Wavelength (nm)")
        fig.colorbar(im, ax=ax, label="mAU")
        fig.tight_layout()
        cv = FigureCanvasTkAgg(fig, master=w)
        NavigationToolbar2Tk(cv, w).update()
        cv.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        cv.draw()

    # ------------------------------------------------------------------ esportazione e sessione
    def tabella_picchi(self):
        """DataFrame con i picchi di tutte le tracce: quelli calcolati qui (Source = HPLCManager) e
        quelli scritti da ChemStation nel REPORTnn.CSV (Source = ChemStation), uno per riga."""
        righe = []
        for n, c in self.cromatogrammi.items():
            for p in c['picchi']:
                righe.append({'Trace': n, 'Source': 'HPLCManager', 'Peak': p['n'], 'RT (min)': p['rt'],
                              'Height (mAU)': p['altezza'], 'Area (mAU*s)': p['area'],
                              'Area %': p['pct'], 'Start (min)': p['inizio'], 'End (min)': p['fine'],
                              'Width (min)': np.nan, 'Type': ''})
            for p in c['picchi_cs']:
                righe.append({'Trace': n, 'Source': 'ChemStation', 'Peak': p['n'], 'RT (min)': p['rt'],
                              'Height (mAU)': p['altezza'], 'Area (mAU*s)': p['area'],
                              'Area %': p['area_pct'], 'Start (min)': np.nan, 'End (min)': np.nan,
                              'Width (min)': p['larghezza'], 'Type': p['tipo']})
        return pd.DataFrame(righe, columns=['Trace', 'Source', 'Peak', 'RT (min)', 'Height (mAU)',
                                            'Area (mAU*s)', 'Area %', 'Start (min)', 'End (min)',
                                            'Width (min)', 'Type'])

    def tabella_spettri(self):
        """DataFrame degli spettri scelti e selezionati (spuntati) nella finestra degli spettri:
        lunghezza d'onda + una colonna per spettro, intestata con dataset e tempo."""
        scelti = [k for k in self._spettri_scelti if k.get('visibile', True)]
        if not scelti:
            return pd.DataFrame()
        colonne = []
        for k in scelti:
            nome = '%s t=%.3f min' % (k['dataset'], k['t'])
            colonne.append(pd.Series(k['y'], index=pd.Index(k['wl'], name='Wavelength (nm)'), name=nome))
        return pd.concat(colonne, axis=1).sort_index()

    def _scegli_file_salvataggio(self, titolo, estensione, tipi):
        return filedialog.asksaveasfilename(title=titolo, initialdir=self.current_dir,
                                            defaultextension=estensione, filetypes=tipi)

    def _scrivi_csv(self, df, path, titolo, indice=True):
        try:
            df.to_csv(path, sep=';', na_rep='', encoding='latin-1', index=indice)
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror(titolo, "Salvataggio fallito:\n%s" % e)
            return False
        self.current_dir = os.path.dirname(path)
        return True

    def esporta_picchi(self):
        df = self.tabella_picchi()
        if df.empty:
            messagebox.showinfo("Export Peaks", "Nessun picco: usa Tools > Peaks (Detect & integrate) "
                                "oppure apri una cartella .D con i canali registrati (picchi ChemStation).")
            return
        path = self._scegli_file_salvataggio("Export Peaks", '.csv', [("CSV", "*.csv")])
        if path and self._scrivi_csv(df, path, "Export Peaks", indice=False):
            messagebox.showinfo("Export Peaks", "%d picchi salvati:\n%s" % (len(df), path))

    def esporta_spettri(self):
        df = self.tabella_spettri()
        if df.empty:
            messagebox.showinfo("Export Spectra", "Nessuno spettro scelto: clicca sul cromatogramma "
                                "(Tools > Spectra) oppure usa 'Spectrum at (min)' e Show.")
            return
        path = self._scegli_file_salvataggio("Export Spectra", '.csv', [("CSV", "*.csv")])
        if path and self._scrivi_csv(df, path, "Export Spectra"):
            messagebox.showinfo("Export Spectra", "%d spettri salvati:\n%s" % (df.shape[1], path))

    def esporta_dad_completo(self):
        """Tutto il cubo DAD del dataset selezionato: una riga per tempo, una colonna per lunghezza
        d'onda (con 13493 spettri x 201 lunghezze d'onda il file e' di circa 25 MB)."""
        d = self._sp_dataset()
        if d is None:
            messagebox.showinfo("Export Full DAD Data", "Nessun dato DAD caricato (apri una cartella .D con il .uv).")
            return
        path = self._scegli_file_salvataggio("Export Full DAD Data", '.csv', [("CSV", "*.csv")])
        if not path:
            return
        df = pd.DataFrame(d['S'], index=pd.Index(d['t'], name='Time (min)'),
                          columns=['%g' % w for w in d['wl']])
        self.root.config(cursor='watch')
        self.root.update_idletasks()
        try:
            ok = self._scrivi_csv(df, path, "Export Full DAD Data")
        finally:
            self.root.config(cursor='')
        if ok:
            messagebox.showinfo("Export Full DAD Data", "%d spettri x %d lunghezze d'onda salvati:\n%s"
                                % (df.shape[0], df.shape[1], path))

    def _dati_sessione(self):
        return {'version': 1, 'cromatogrammi': self.cromatogrammi, 'spettri': self.spettri,
                'sp_live': dict(self._sp_live), 'includi_ch': bool(self.var_includi_ch.get()),
                'strumento': self.strumento}

    def _applica_sessione(self, data):
        self.cromatogrammi = data['cromatogrammi']
        self.spettri = data.get('spettri', {})
        self.strumento = data.get('strumento', {})
        self._strum_chiudi()
        self._sp_live = dict(data.get('sp_live', {}))
        self.var_includi_ch.set(bool(data.get('includi_ch', False)))
        self._sp_chiudi_finestra()
        self._reset_vista = True
        for w in self.f_tool_host.winfo_children():
            w.destroy()
        for attr in ('var_sp_set', 'var_sp_l'):          # il pannello Spectra si ricostruisce da zero
            if hasattr(self, attr):
                delattr(self, attr)
        self.aggiorna_vista()
        if self.spettri:
            self.apri_spettri()

    def salva_sessione(self):
        """Salva tracce (comprese quelle derivate), picchi e spettri DAD in un unico file, per
        riprendere il lavoro piu' avanti senza rifare i passaggi."""
        if not self.cromatogrammi:
            messagebox.showwarning("Save Session", "Nessun cromatogramma caricato.")
            return
        path = self._scegli_file_salvataggio("Save Session", '.hplcsession',
                                             [("HPLC Manager Session", "*.hplcsession"), ("All Files", "*.*")])
        if not path:
            return
        try:
            with gzip.open(path, 'wb') as f:
                pickle.dump(self._dati_sessione(), f)
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror("Save Session", "Salvataggio fallito:\n%s" % e)
            return
        self.current_dir = os.path.dirname(path)
        self._dirty = False
        messagebox.showinfo("Save Session", "Sessione salvata:\n%s" % path)

    def apri_sessione(self):
        """Carica una sessione salvata, sostituendo tutto quanto e' caricato."""
        if self._dirty and not messagebox.askyesno(
                "Open Session", "Ci sono tracce derivate non esportate ne' salvate. Scartarle e aprire la sessione?"):
            return
        path = filedialog.askopenfilename(initialdir=self.current_dir, title="Open Session",
                                          filetypes=[("HPLC Manager Session", "*.hplcsession"), ("All Files", "*.*")])
        if not path:
            return
        try:
            with gzip.open(path, 'rb') as f:
                data = pickle.load(f)
            if not isinstance(data, dict) or 'cromatogrammi' not in data:
                raise ValueError("non e' una sessione di HPLC Manager")
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror("Open Session", "Apertura sessione fallita:\n%s" % e)
            return
        self.current_dir = os.path.dirname(path)
        self._applica_sessione(data)
        self._dirty = False

    def _copia_figura(self):
        """Copia della figura corrente per salvataggio/editor: senza il cursore e senza le linee
        degli spettri scelti, ristilizzata Origin 'single' (4:3) invece della dimensione del pannello."""
        self._pulisci_cursore()
        fig = pickle.loads(pickle.dumps(self.fig))
        for ax in fig.axes:
            for a in list(ax.lines) + list(ax.texts):
                if a.get_gid() == 'marcatore':
                    a.remove()
        if origin_style is not None and fig.axes:
            origin_style.applica_stile_origin(fig.axes[0], fig, set_size=True, preset='single')
        return fig

    def salva_figura_immagine(self):
        if not self.cromatogrammi:
            messagebox.showwarning("Save Figure Image", "Nessun cromatogramma caricato.")
            return
        path = self._scegli_file_salvataggio("Save Figure Image", '.png', [
            ("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg"), ("All Files", "*.*")])
        if not path:
            return
        try:
            self._copia_figura().savefig(path, dpi=300)
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror("Save Figure Image", "Salvataggio fallito:\n%s" % e)
            return
        self.current_dir = os.path.dirname(path)
        messagebox.showinfo("Save Figure Image", "Figura salvata:\n%s" % path)

    def salva_figura_pickle(self):
        """Salva la figura corrente come pickle matplotlib: riapribile con plot_editor.pyw come
        oggetto Figure/Axes vivo (non un raster)."""
        if not self.cromatogrammi:
            messagebox.showwarning("Save Figure", "Nessun cromatogramma caricato.")
            return
        path = self._scegli_file_salvataggio("Save Figure", '.fig.pickle', [
            ("Matplotlib Figure (pickle)", "*.pickle *.pkl"), ("All Files", "*.*")])
        if not path:
            return
        try:
            with open(path, 'wb') as f:
                pickle.dump(self._copia_figura(), f)
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror("Save Figure", "Salvataggio fallito:\n%s" % e)
            return
        self.current_dir = os.path.dirname(path)
        messagebox.showinfo("Save Figure", "Figura salvata:\n%s" % path)

    def _carica_plot_editor(self):
        """Importa (una sola volta) plot_editor: prima una copia locale, poi il repo fratello PlotStyleKit."""
        if self._pe_module is None:
            import importlib.util
            here = os.path.dirname(os.path.abspath(__file__))
            candidati = [os.path.join(here, 'plot_editor.pyw'),
                         os.path.join(here, '..', 'PlotStyleKit', 'plot_editor.pyw')]
            path = next((c for c in candidati if os.path.isfile(c)), None)
            if path is None:
                raise FileNotFoundError(
                    "plot_editor.pyw non trovato (repo PlotStyleKit mancante accanto a questo progetto)")
            spec = importlib.util.spec_from_file_location('plot_editor', path)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            self._pe_module = mod
        return self._pe_module

    def apri_editor_figura(self):
        """Apre il Plot Editor di PlotStyleKit sulla figura corrente (una copia indipendente)."""
        if not self.cromatogrammi:
            messagebox.showwarning("Edit Figure", "Nessun cromatogramma caricato.")
            return
        try:
            pe = self._carica_plot_editor()
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror("Edit Figure", "plot_editor.pyw non disponibile:\n%s" % e)
            return
        try:
            fig = self._copia_figura()
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror("Edit Figure", "Impossibile duplicare la figura:\n%s" % e)
            return
        top = tk.Toplevel(self.root)
        top.geometry("1300x820")
        editor = pe.PlotEditor(top)
        editor.carica_figura(fig, title="figura corrente")

    # ------------------------------------------------------------------ curve dello strumento
    GRUPPI_STRUM = (('bar', 'Pressure (bar)'), ('%', 'Solvent composition (%)'),
                    ('ml/min', 'Flow (ml/min)'), ('\xb0C', 'Temperature (\xb0C)'))

    def _strum_chiudi(self):
        w = self._strum_win
        self._strum_win = None
        if w is not None:
            try:
                w.destroy()
            except Exception:
                pass

    def _strum_dataset(self):
        """Dataset (cartella .D) di cui mostrare le curve: quello scelto nella finestra, o il primo."""
        v = getattr(self, 'var_strum_set', None)
        if v is not None and v.get() in self.strumento and self.strumento[v.get()]['segnali']:
            return v.get()
        for k, d in self.strumento.items():
            if d['segnali']:
                return k
        return None

    def tabella_strumento(self, ds):
        """DataFrame delle curve del dataset `ds` sulla scala dei tempi del segnale piu' veloce
        (pressione, flusso e solventi a 0.3 s); i segnali piu' lenti (temperatura, 1 s) sono
        interpolati linearmente su quella scala."""
        seg = self.strumento[ds]['segnali']
        base = min(seg.values(), key=lambda d: d['t'][1] - d['t'][0] if len(d['t']) > 1 else 1e9)
        out = pd.DataFrame(index=pd.Index(np.round(base['t'], 6), name='Time (min)'))
        for titolo, d in seg.items():
            out['%s [%s]' % (titolo, d['unit'])] = np.interp(base['t'], d['t'], d['y'])
        return out

    def apri_curve_strumento(self):
        ds = self._strum_dataset()
        if ds is None:
            messagebox.showinfo("Instrument curves", "Nessuna curva dello strumento: apri una cartella .D "
                                "che contenga LCDIAG.REG.")
            return
        self._strum_chiudi()
        w = tk.Toplevel(self.root)
        w.title("Instrument curves - %s" % ds)
        w.geometry("900x720")
        top = tk.Frame(w)
        top.pack(fill=tk.X, padx=6, pady=4)
        self.var_strum_set = tk.StringVar(value=ds)
        nomi = [k for k, d in self.strumento.items() if d['segnali']]
        if len(nomi) > 1:
            cb = ttk.Combobox(top, textvariable=self.var_strum_set, values=nomi, state='readonly', width=24)
            cb.pack(side=tk.LEFT)
            cb.bind('<<ComboboxSelected>>', lambda e: self._strum_disegna())
        tk.Button(top, text="Export CSV", command=self.esporta_curve_strumento).pack(side=tk.RIGHT)
        w.fig = Figure()
        w.canvas = FigureCanvasTkAgg(w.fig, master=w)
        NavigationToolbar2Tk(w.canvas, w).update()
        w.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self._strum_win = w
        w.bind('<Destroy>', lambda e: setattr(self, '_strum_win', None) if e.widget is w else None)
        self._strum_disegna()
        return w

    def _strum_disegna(self):
        w = self._strum_win
        ds = self._strum_dataset()
        if w is None or ds is None:
            return
        w.title("Instrument curves - %s" % ds)
        seg = self.strumento[ds]['segnali']
        gruppi = []
        for unita, titolo in self.GRUPPI_STRUM:
            membri = [(t, d) for t, d in seg.items() if d['unit'] == unita]
            if unita == '%':      # solo i solventi che non restano a zero
                membri = [(t, d) for t, d in membri if np.any(d['y'] > 0)]
            if membri:
                gruppi.append((titolo, membri))
        w.fig.clear()
        assi = w.fig.subplots(len(gruppi), 1, sharex=True) if gruppi else []
        assi = np.atleast_1d(assi)
        for ax, (titolo, membri) in zip(assi, gruppi):
            for t, d in membri:
                etichetta = t.split(', ', 1)[-1]
                ax.plot(d['t'], d['y'], lw=1, label=etichetta)
            ax.set_ylabel(titolo, fontsize=9)
            ax.tick_params(labelsize=8)
            if len(membri) > 1:
                ax.legend(fontsize=8, loc='upper right')
        if len(assi):
            assi[-1].set_xlabel("Time (min)")
        w.fig.tight_layout()
        w.canvas.draw()

    def esporta_curve_strumento(self):
        ds = self._strum_dataset()
        if ds is None:
            messagebox.showinfo("Export Instrument Curves", "Nessuna curva dello strumento caricata.")
            return
        path = self._scegli_file_salvataggio("Export Instrument Curves", '.csv', [("CSV", "*.csv")])
        if not path:
            return
        df = self.tabella_strumento(ds)
        if self._scrivi_csv(df, path, "Export Instrument Curves"):
            messagebox.showinfo("Export Instrument Curves", "%d punti, %d curve salvati:\n%s"
                                % (len(df), df.shape[1], path))

    def on_exit(self):
        if self._dirty and not messagebox.askyesno(
                "Esci", "Ci sono tracce derivate non esportate. Uscire?"):
            return
        self.root.destroy()


def main():
    root = TkinterDnD.Tk() if HAS_DND else tk.Tk()
    app = HPLCManager(root)
    for a in sys.argv[1:]:
        app.processa_file(a)
    if len(sys.argv) > 1:
        app.aggiorna_vista()
    root.mainloop()


if __name__ == '__main__':
    main()
