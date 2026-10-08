"""HPLCManager: viewer and analysis tool for HPLC chromatograms (Agilent ChemStation)."""
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


def _load_origin_style():
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
    origin_style = _load_origin_style()
    if origin_style is not None:
        origin_style.applica_rcparams()
except Exception:
    origin_style = None


class HPLCManager:
    """Manager for HPLC chromatograms.

    Data model: self.chromatograms = {name: {'df': DataFrame indexed by time (min) with
    a single 'mAU' column, 'info': {...}, 'peaks': [...], 'cs_peaks': [...], 'hidden': bool}},
    modeled on LabSpectrumManager. 'peaks' are the peaks computed here, 'cs_peaks'
    those written by ChemStation in the REPORTnn.CSV of the same .D folder (if present).
    """

    # Scale factor of the ChemStation .ch files: decoded value / 2000 = mAU
    # (inferred by comparison with Report.TXT, not read from the header).
    CH_SCALE = 2000.0
    # Like the .ch files, the .uv values are in 1/2000 mAU (verified by extracting 280/295/320/350 nm
    # from the .uv and comparing with dad1A..D.ch). Band extraction uses the wavelength axis
    # shifted by UV_OFFSET_NM: an empirical value (the deviation from the .ch files is minimal
    # between 0 and 1 nm, ~0.45), not read from the header.
    UV_SCALE = 2000.0
    UV_OFFSET_NM = 0.45

    def __init__(self, root):
        self.root = root
        self.root.title("HPLC Manager")
        self.root.geometry("1500x800")
        self.chromatograms = {}
        self._dirty = False
        self._cursor_artists = []
        self._tool_frame = None
        self.var_cs_peaks = tk.BooleanVar(value=True)
        self.spectra = {}          # {.D folder name: {'t', 'wl', 'S', 'info'}} from the .uv files
        self._spec_win = None
        self.instrument = {}        # {.D folder name: {'signals': {title: {t, y, unit}}, 'modules': [...], 'column': {...}}}
        self._instr_win = None
        self.var_click_spec = tk.BooleanVar(value=True)
        self._marker_artists = []   # matplotlib artists of the markers on the chromatogram (to remove them without redrawing)
        self._sp_color_n = 0          # counter that assigns a stable color to each picked spectrum
        self._sp_multi_last = []     # values of the last 'Several wavelengths' dialog (reopened unchanged)
        self._chosen_spectra = []      # spectra picked with click/Show: {'dataset','t','wl','y','color','visible','line'}
        self._pe_module = None         # plot_editor module (PlotStyleKit), loaded on first use
        self.current_dir = os.getcwd()
        self._reset_view = True       # next redraw: automatic axes (new data) instead of keeping the zoom
        self.var_include_ch = tk.BooleanVar(value=False)   # also load the recorded channels (.ch)
        self._sp_live = {}         # {dataset name: name of the 'live' extracted trace}

        self._build_menu()

        self.paned = tk.PanedWindow(root, orient=tk.HORIZONTAL, sashrelief=tk.RAISED, sashwidth=4)
        self.paned.pack(fill=tk.BOTH, expand=True)

        # --- left panel: data table ---
        self.f_left = tk.Frame(self.paned, width=400)
        tk.Label(self.f_left, text="Data Table", font=('Segoe UI', 10, 'bold')).pack(anchor='w', padx=4)
        self.table = ttk.Treeview(self.f_left, show='headings', selectmode='extended')
        sy = ttk.Scrollbar(self.f_left, orient='vertical', command=self.table.yview)
        sx = ttk.Scrollbar(self.f_left, orient='horizontal', command=self.table.xview)
        self.table.configure(yscrollcommand=sy.set, xscrollcommand=sx.set)
        sy.pack(side=tk.RIGHT, fill=tk.Y)
        sx.pack(side=tk.BOTTOM, fill=tk.X)
        self.table.pack(fill=tk.BOTH, expand=True)
        self.table.bind('<Control-c>', self._copy_data_table)
        self.table.bind('<Control-a>', self._select_all_table)
        self.paned.add(self.f_left, width=400, minsize=200, stretch='never')

        # --- central panel: plot ---
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
        self.canvas.mpl_connect('resize_event', self._trim_margins)
        self.canvas.mpl_connect('button_press_event', self._press_spectrum)
        self.canvas.mpl_connect('button_release_event', self._click_spectrum)
        self.paned.add(self.f_plot, width=700, minsize=300, stretch='always')

        # --- right panel: scrollable, everything inside self.f_right_inner ---
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
        self._build_right_panel()

        self.root.protocol("WM_DELETE_WINDOW", self.on_exit)
        if HAS_DND:
            try:
                self.root.drop_target_register(DND_FILES)
                self.root.dnd_bind('<<Drop>>', self.handle_drop)
            except Exception:
                pass
        self.refresh_view()

    # ------------------------------------------------------------------ menus and panels
    def _build_menu(self):
        self.menu_bar = tk.Menu(self.root)
        m = tk.Menu(self.menu_bar, tearoff=0)
        m.add_command(label="Open ChemStation folder (.D)...", command=self.load_folder_dialog)
        m.add_command(label="Open signal file (.ch)...", command=self.load_file_dialog)
        m.add_checkbutton(label="Also load recorded channels (.ch) from .D folders",
                          variable=self.var_include_ch)
        m.add_separator()
        m.add_command(label="Open Session...", command=self.open_session)
        m.add_command(label="Save Session...", command=self.save_session)
        m.add_separator()
        m.add_command(label="Export Traces (CSV)...", command=self.export_csv)
        m.add_command(label="Export Peaks (CSV)...", command=self.export_peaks)
        m.add_command(label="Export Spectra (CSV)...", command=self.export_spectra)
        m.add_command(label="Export Full DAD Data (CSV)...", command=self.export_full_dad)
        m.add_command(label="Export Instrument Curves (CSV)...", command=self.export_instrument_curves)
        m.add_separator()
        m.add_command(label="Save Figure Image (PNG, PDF, SVG)...", command=self.save_figure_image)
        m.add_command(label="Save Figure (pickle)...", command=self.save_figure_pickle)
        m.add_command(label="Edit Figure...", command=self.open_figure_editor)
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
        m.add_checkbutton(label="Show ChemStation peaks", variable=self.var_cs_peaks,
                          command=self.refresh_view)
        self.menu_bar.add_cascade(label="View", menu=m)

        m = tk.Menu(self.menu_bar, tearoff=0)
        m.add_command(label="Peaks...", command=self.open_peaks)
        m.add_command(label="Smoothing...", command=self.open_smoothing)
        m.add_command(label="Trim...", command=self.open_trim)
        m.add_command(label="Normalize", command=self.normalize)
        m.add_separator()
        m.add_command(label="Spectra (DAD)...", command=self.open_spectra)
        m.add_command(label="Several wavelengths...", command=self.open_multi_extraction)
        m.add_command(label="Instrument curves (pressure, gradient...)", command=self.open_instrument_curves)
        self.menu_bar.add_cascade(label="Tools", menu=m)
        self.root.config(menu=self.menu_bar)

    def _build_right_panel(self):
        p = self.f_right_inner
        tk.Label(p, text="Chromatograms", font=('Segoe UI', 10, 'bold')).pack(anchor='w', padx=4, pady=(4, 0))
        fl = tk.Frame(p)
        fl.pack(fill=tk.X, padx=4)
        self.listbox = tk.Listbox(fl, selectmode=tk.EXTENDED, height=6, exportselection=False)
        self.listbox.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.listbox.bind('<<ListboxSelect>>', self._on_selection)
        self.listbox.bind('<Button-3>', self._list_menu)
        fb = tk.Frame(p)
        fb.pack(fill=tk.X, padx=4, pady=2)
        for txt, cmd in (("Remove", self.remove_selected), ("Clear", self.clear_all),
                         ("Hide", self.hide_selected), ("Show", self.show_selected)):
            tk.Button(fb, text=txt, command=cmd).pack(side=tk.LEFT, padx=1)
        tk.Button(p, text="Open ChemStation folder (.D)...",
                  command=self.load_folder_dialog).pack(fill=tk.X, padx=4, pady=2)

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
        for txt, cmd in (("Peaks", self.open_peaks), ("Smoothing", self.open_smoothing),
                         ("Trim", self.open_trim), ("Normalize", self.normalize),
                         ("Spectra", self.open_spectra)):
            tk.Button(fb, text=txt, command=cmd).pack(side=tk.LEFT, padx=1)
        tk.Button(p, text="Instrument curves (pressure, gradient...)",
                  command=self.open_instrument_curves).pack(fill=tk.X, padx=4, pady=(3, 0))
        self.f_tool_host = tk.Frame(p)
        self.f_tool_host.pack(fill=tk.X, padx=4, pady=4)

    def _list_menu(self, event):
        i = self.listbox.nearest(event.y)
        if i >= 0 and i not in self.listbox.curselection():
            self.listbox.selection_clear(0, tk.END)
            self.listbox.selection_set(i)
        menu = tk.Menu(self.root, tearoff=0)
        menu.add_command(label="Hide", command=self.hide_selected)
        menu.add_command(label="Show", command=self.show_selected)
        menu.add_command(label="Remove", command=self.remove_selected)
        menu.tk_popup(event.x_root, event.y_root)

    def _open_tool_panel(self, title):
        """Clear the tool area in the right panel and return the frame in which to build a new one."""
        for w in self.f_tool_host.winfo_children():
            w.destroy()
        f = tk.LabelFrame(self.f_tool_host, text=title)
        f.pack(fill=tk.X)
        self._tool_frame = f
        return f

    # ------------------------------------------------------------------ file reading
    def read_ch(self, path):
        """Read a ChemStation .ch signal, version 30 (big-endian).

        Returns (time_min, signal_mAU, info). The data start at 0x400, in segments
        '0x10 n' followed by n entries: int16 = delta from the previous value, or
        0x8000 + int32 = absolute value (which resets the accumulation).
        """
        with open(path, 'rb') as f:
            b = f.read()
        if b[0:3] != b'\x0230':
            raise ValueError(".ch file not supported (version other than 30): "
                             + os.path.basename(path))

        def pascal(o):
            return b[o + 1:o + 1 + b[o]].decode('latin-1')

        # Pascal strings at fixed offsets (version, file type, sample name, date, module, method)
        info = {'Version': pascal(0), 'Type': pascal(4), 'Sample': pascal(0x18)}
        info['Date'] = pascal(0xB2)
        info['Module'] = pascal(0xD0)
        info['Method'] = pascal(0xE4)
        # signal description ("DAD A, Sig=280,4 Ref=360,100"): Pascal string before 0x400
        i = b.find(b'Sig=', 0, 0x400)
        if i > 0:
            s = i
            while s > 1 and b[s - 1] >= 32:
                s -= 1
            info['Signal'] = b[s:s + b[s - 1]].decode('latin-1')

        t0, t1 = struct.unpack('>ii', b[0x11A:0x122])
        data = []
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
                data.append(v)
        y = np.array(data, dtype=float) / self.CH_SCALE
        dt = (t1 - t0) / (len(y) - 1) if len(y) > 1 else 0.0
        t = (t0 + np.arange(len(y)) * dt) / 60000.0
        return t, y, info

    def read_report_csv(self, path):
        """Read a ChemStation REPORTnn.CSV (one peak table per signal, no header).

        Returns a list of dicts: n, rt, kind, width, area, height, area_pct.
        """
        peaks = []
        with open(path, newline='', encoding='latin-1') as f:
            for r in csv.reader(f):
                if len(r) < 7:
                    continue
                peaks.append({'n': int(r[0]), 'rt': float(r[1]), 'kind': r[2].strip(),
                               'width': float(r[3]), 'area': float(r[4]),
                               'height': float(r[5]), 'area_pct': float(r[6])})
        return peaks

    def read_uv(self, path):
        """Read a ChemStation DAD spectra file (.uv), version 31.

        Returns (time_min, wavelengths_nm, S, info) with S of shape (n_times, n_wavelengths) in mAU.
        After the header (records starting at 0x200) each spectrum is a little-endian record:
        uint16 type (0x43), uint16 record length, uint32 time (ms); then uint16 start
        wavelength, end wavelength and step (in 1/20 nm), 8 undecoded bytes; then n values encoded as
        in the .ch files but little-endian: int16 = delta from the previous one, 0x8000 + int32 = absolute value.
        After the last record there is a tail (index?) that is not read.
        """
        with open(path, 'rb') as f:
            b = f.read()
        if b[0:3] != b'\x0231':
            raise ValueError(".uv file not supported (version other than 31): "
                             + os.path.basename(path))
        def pascal(o):
            return b[o + 1:o + 1 + b[o]].decode('latin-1')

        # same Pascal strings at fixed offsets as the .ch files (sample name, date, module, method)
        info = {'Version': pascal(0), 'Type': pascal(4), 'Sample': pascal(0x18),
                'Date': pascal(0xB2), 'Module': pascal(0xD0), 'Method': pascal(0xE4)}
        times, rows = [], []
        s, wl = 0x200, None
        while s + 22 <= len(b):
            kind, length, ms = struct.unpack('<HHI', b[s:s + 8])
            if kind != 0x43 or length < 22 or s + length > len(b):
                break
            w0, w1, st = struct.unpack('<HHH', b[s + 8:s + 14])
            if st == 0:
                break
            n = (w1 - w0) // st + 1
            if wl is None:
                wl = (w0 + st * np.arange(n)) / 20.0
            elif n != len(wl):
                break
            data = b[s + 22:s + length]
            if len(data) == 2 * n:                       # no absolute value in the record
                v = np.cumsum(np.frombuffer(data, dtype='<i2').astype(np.int64))
            else:
                v, val, o = np.empty(n, dtype=np.int64), 0, 0
                for k in range(n):
                    if data[o:o + 2] == b'\x00\x80':
                        val = struct.unpack('<i', data[o + 2:o + 6])[0]
                        o += 6
                    else:
                        val += struct.unpack('<h', data[o:o + 2])[0]
                        o += 2
                    v[k] = val
            times.append(ms)
            rows.append(v)
            s += length
        if not rows:
            raise ValueError("No readable spectrum in " + os.path.basename(path))
        S = np.array(rows, dtype=float) / self.UV_SCALE
        info['Spectra'] = '%d spectra, %.0f-%.0f nm, step %.1f nm' % (len(rows), wl[0], wl[-1], wl[1] - wl[0])
        return np.array(times, dtype=float) / 60000.0, wl, S, info

    def extract_wavelength(self, wl, S, center, bandwidth, ref=None, ref_bandwidth=None):
        """Chromatogram extracted at `center` (nm) with width `bandwidth`, minus the reference band.

        Like the ChemStation signal 'Sig=280,4 Ref=360,100': weighted mean of the band
        [center - bandwidth/2, center + bandwidth/2], with weight proportional to the overlap between
        the band and the spectral interval of each point (step wl[1]-wl[0]); the same for the
        reference. Returns the signal in mAU, one value per spectrum.
        """
        step = float(wl[1] - wl[0])
        x = wl + self.UV_OFFSET_NM

        def weights(c, bw):
            lo, hi = c - bw / 2.0, c + bw / 2.0
            w = np.clip(np.minimum(x + step / 2, hi) - np.maximum(x - step / 2, lo), 0, None)
            if w.sum() == 0:
                raise ValueError("Band %.1f +/- %.1f nm outside the acquired spectrum (%.0f-%.0f nm)"
                                 % (c, bw / 2.0, wl[0], wl[-1]))
            return w / w.sum()

        y = S @ weights(center, bandwidth)
        if ref is not None and ref_bandwidth:
            y = y - S @ weights(ref, ref_bandwidth)
        return y

    # Units and scale factors of the LCDIAG.REG signals: the stored value (uint32) is the
    # physical value multiplied by this factor (verified on pressure, flow, solvents, temperature).
    DIAG_UNITS = {'bar': 100.0, 'ml/min': 1000.0, '%': 10.0, '\xb0C': 100.0}

    def read_diagnostics(self, path):
        """Read LCDIAG.REG of ChemStation: the profiles recorded by the instrument during the run
        (pressure, flow, solvent composition, column temperature).

        Returns {title: {'t': minutes (array), 'y': values (array), 'unit': 'bar'|'ml/min'|'%'|'\xb0C'}}
        in order of appearance. Structure (little-endian): for each signal an array of uint32
        preceded by the strings 'min\\0<unit>\\0'; the number of points is 168 bytes before 'min'
        (uint32), the sampling interval in minutes 116 bytes before (double). The title
        ('PMP1, Pressure') is in the header that FOLLOWS the array. Time starts at 0 (start of run):
        this is an assumption, consistent with the pressure at 4 s and the final pressure in the report.
        """
        with open(path, 'rb') as f:
            b = f.read()
        if b[3:4] != b'\x00' or b'REGISTER FILE' not in b[:24]:
            raise ValueError(".REG file not recognized: " + os.path.basename(path))
        signals = {}
        for k, m in enumerate(re.finditer(rb'min\x00(bar|ml/min|%|\xb0C)\x00', b)):
            u, e, unit = m.start(), m.end(), m.group(1).decode('latin-1')
            if u < 168:
                continue
            n = struct.unpack('<I', b[u - 168:u - 164])[0]
            dt = struct.unpack('<d', b[u - 116:u - 108])[0]
            if not (0 < n <= 10000000) or e + 4 * n > len(b) or not (0 < dt < 10):
                continue
            y = np.frombuffer(b[e:e + 4 * n], dtype='<u4').astype(float) / self.DIAG_UNITS[unit]
            header_block = b[e + 4 * n:e + 4 * n + 1500]
            t = re.search(rb'((?:[A-Z]{3}\d), [^\x00]+)\x00arial', header_block)
            title = t.group(1).decode('latin-1') if t else '%s signal %d' % (unit, k + 1)
            if title in signals:
                title = '%s (%d)' % (title, k + 1)
            signals[title] = {'t': np.arange(n) * dt, 'y': y, 'unit': unit}
        if not signals:
            raise ValueError("No signal recognized in " + os.path.basename(path))
        return signals

    def read_acqres(self, path):
        """Read ChemStation ACQRES.REG: instrument modules and column data.

        Returns {'modules': [{'name', 'part', 'serial', 'firmware', 'build'}], 'column': {...}}.
        Register with column-based tables: strings are uint16 (length including NUL) + text + NUL.
        Modules: 5 columns of consecutive strings (serialNumber, FWrevision, buildNumber, Name,
        PartNumber), one row per module. Column: description = the string that is not a software
        version, instrument name or path; length, diameter and particle size are the first three 'round'
        doubles of the 'Acquisition Results' table (ColLength, ColDiameter, ParticleSize: assignment
        inferred from the field order and typical values, 100 x 2.1 mm, 5 um, not verified on other files).
        """
        with open(path, 'rb') as f:
            b = f.read()
        if b[3:4] != b'\x00' or b'REGISTER FILE' not in b[:24]:
            raise ValueError(".REG file not recognized: " + os.path.basename(path))
        text_ok = re.compile(r'^[\x20-\x7e\xa0-\xff]+$')

        def string(o):
            if o < 0 or o + 2 > len(b):
                return None
            n = struct.unpack('<H', b[o:o + 2])[0]
            if n < 2 or n > 200 or o + 2 + n > len(b) or b[o + 1 + n] != 0:
                return None
            t = b[o + 2:o + 1 + n].decode('latin-1')
            return t if text_ok.match(t) and any(c.isalnum() for c in t) else None

        def chain(o):
            out = []
            while True:
                t = string(o)
                if t is None:
                    return out
                out.append((o, t))
                o += 2 + len(t) + 1
        best = []
        o = 0
        while o < len(b) - 2:
            c = chain(o)
            if len(c) > len(best):
                best = c
            o += max(1, sum(2 + len(t) + 1 for _, t in c) if c else 1)
        modules = []
        columns = [b'serialNumber', b'FWrevision', b'buildNumber', b'Name', b'PartNumber']
        positions = [b.find(c) for c in columns]
        if len(best) >= 5 and len(best) % 5 == 0 and all(x >= 0 for x in positions) \
                and positions == sorted(positions):
            r = len(best) // 5
            val = [t for _, t in best]
            for i in range(r):
                modules.append({'serial': val[i], 'firmware': val[r + i], 'build': val[2 * r + i],
                               'name': val[3 * r + i], 'part': val[4 * r + i]})
        chain_end = (best[-1][0] + 3 + len(best[-1][1])) if best else 0
        chain_start = best[0][0] if best else 0
        others = []
        o = 0
        while o < len(b) - 2:
            t = string(o)
            if t is not None and len(t) >= 3 and not (chain_start <= o < chain_end):
                others.append(t)
                o += 2 + len(t) + 1
            else:
                o += 1
        column = {}
        for t in others:
            if t.startswith('Rev.') or '\\' in t or t.startswith('Instrument') or t.endswith('Results'):
                continue
            if re.fullmatch(r'[A-Za-z][A-Za-z0-9 ./-]{3,}', t) and not re.fullmatch(r'[A-Za-z]+[A-Z][a-z]+[A-Za-z]*', t):
                column['description'] = t
                break
        ini_pos = b.find(b'WriteInfo')
        end = chain_start if chain_start else len(b)
        values = []
        o = ini_pos + 9 if ini_pos >= 0 else 0
        while o + 8 <= end:
            v = struct.unpack('<d', b[o:o + 8])[0]
            if v == v and 0.4 <= abs(v) <= 1000 and abs(v * 10 - round(v * 10)) < 1e-9:
                values.append(v)
                o += 8
            else:
                o += 1
        if len(values) >= 2:
            column['length_mm'], column['diameter_mm'] = values[0], values[1]
        if len(values) >= 3 and values[2] <= 20:
            column['particle_um'] = values[2]
        return {'modules': modules, 'column': column}

    def read_folder_metadata(self, folder):
        """Analysis metadata from the .D folder: Report00.CSV (sample, injection, sequence,
        pressure, flow, solvents, signals) and RUN.LOG (column temperature).

        Returns a dict {label: value} in order of appearance; empty if the files are missing.
        """
        info = {}
        rep00 = os.path.join(folder, 'Report00.CSV')
        if os.path.isfile(rep00):
            skip = ('Report Title', 'Sorted By', 'Multiplier', 'Dilution', 'Number of Columns',
                    'Method Info', 'Results Created by')
            with open(rep00, newline='', encoding='latin-1') as f:
                for r in csv.reader(f):
                    if len(r) < 2 or not r[0].strip() or r[0].startswith('Column ') or r[0] in skip:
                        continue
                    value = ''.join(x for x in r[1:]).strip()
                    if r[0] in ('Data File', 'Analysis Method', 'Sequence File'):
                        value = ''.join(r[1:3])               # path + file name
                    elif len(r) > 2 and r[2].strip() and r[0].startswith(('Inj Vol', 'Actual', 'Start', 'Stop')):
                        value = '%s %s' % (r[1].strip(), r[2].strip())
                    elif r[0].startswith('Solvent'):
                        # lines without solvent (empty third field) are not shown
                        value = ', '.join(x.strip() for x in r[1:] if x.strip()) if len(r) > 2 and r[2].strip() else ''
                    elif r[0].startswith('Signal '):
                        value = ', '.join(x.strip() for x in r[1:] if x.strip())
                    if value:
                        info[r[0]] = value
            # in this export 'Acq. Operator' holds the injection date, not the operator
            if info.get('Acq. Operator') == info.get('Injection Date'):
                del info['Acq. Operator']
        log = os.path.join(folder, 'RUN.LOG')
        if os.path.isfile(log):
            with open(log, encoding='latin-1') as f:
                temp = [float(x) for x in re.findall(r'Column temperature\s*=\s*([\d.]+)', f.read())]
            if temp:
                info['Column temperature'] = ('%.1f' % temp[0] if min(temp) == max(temp)
                                              else '%.1f-%.1f' % (min(temp), max(temp))) + ' \u00b0C'
        acq = os.path.join(folder, 'ACQRES.REG')
        if os.path.isfile(acq):
            try:
                r = self.read_acqres(acq)
            except Exception:
                r = None
            if r:
                c = r['column']
                if c.get('description'):
                    info['Column'] = c['description']
                if 'length_mm' in c:
                    info['Column size'] = '%g x %g mm%s' % (c['length_mm'], c['diameter_mm'], (
                        ', %g um' % c['particle_um']) if 'particle_um' in c else '')
                for m in r['modules']:
                    info[m['name']] = '%s, S/N %s, FW %s' % (m['part'], m['serial'], m['firmware'])
        return info

    def _merge_metadata(self, info, folder):
        """Add to `info` the metadata of the .D folder, without repeating those already read
        from the binary header (same value under a different label)."""
        extra = self.read_folder_metadata(folder)
        for alt_label, old in (('Sample Name', 'Sample'), ('Injection Date', 'Date'),
                               ('Acq. Method', 'Method')):
            if extra.get(alt_label) == info.get(old):
                del extra[alt_label]
        info.update(extra)

    def _load_instrument(self, folder, base):
        """Load the instrument profiles (LCDIAG.REG) and the module and column data (ACQRES.REG)."""
        data = {'signals': {}, 'modules': [], 'column': {}}
        diag = os.path.join(folder, 'LCDIAG.REG')
        if os.path.isfile(diag):
            try:
                data['signals'] = self.read_diagnostics(diag)
            except Exception as e:
                traceback.print_exc()
                messagebox.showwarning("Instrument", "Cannot read LCDIAG.REG:\n%s" % e)
        acq = os.path.join(folder, 'ACQRES.REG')
        if os.path.isfile(acq):
            try:
                r = self.read_acqres(acq)
                data['modules'], data['column'] = r['modules'], r['column']
            except Exception:
                traceback.print_exc()
        if data['signals'] or data['modules']:
            self.instrument[base] = data

    def _load_uv(self, folder, base):
        """Load the DAD spectra (.uv) of a .D folder, if present."""
        uv = [f for f in os.listdir(folder) if f.lower().endswith('.uv')]
        if not uv:
            return
        try:
            t, wl, S, info = self.read_uv(os.path.join(folder, sorted(uv)[0]))
        except Exception as e:
            messagebox.showwarning("Spectra", "Cannot read the DAD spectra:\n%s" % e)
            return
        self._merge_metadata(info, folder)
        self.spectra[base] = {'t': t, 'wl': wl, 'S': S, 'info': info}

    def _chemstation_peaks(self, folder, channel):
        """ChemStation peaks for the channel (e.g. 'dad1A.ch') in the .D folder, or []."""
        letter = os.path.splitext(channel)[0][-1].upper()
        number = None
        rep = os.path.join(folder, 'Report.TXT')
        if os.path.isfile(rep):
            with open(rep, encoding='latin-1') as f:
                for m in re.finditer(r'Signal\s+(\d+)\s*:\s*\w+?\d*\s+([A-Z])\b', f.read()):
                    if m.group(2) == letter:
                        number = int(m.group(1))
                        break
        if number is None and 'A' <= letter <= 'Z':
            number = ord(letter) - ord('A') + 1
        path = os.path.join(folder, 'REPORT%02d.CSV' % (number or 0))
        if number and os.path.isfile(path):
            try:
                return self.read_report_csv(path)
            except Exception:
                return []
        return []

    def _free_name(self, name):
        if name not in self.chromatograms:
            return name
        k = 2
        while '%s (%d)' % (name, k) in self.chromatograms:
            k += 1
        return '%s (%d)' % (name, k)

    def _add_trace(self, name, t, y, info, cs_peaks=None, estr=None):
        """Add a trace. `estr` (only for traces extracted from the DAD spectra) holds the parameters
        it was computed with: {'ds', 'l', 'b', 'ref', 'rb'}; selecting it in the list brings them back
        into the panel, and 'Update trace' modifies it."""
        name = self._free_name(name)
        df = pd.DataFrame({'mAU': y}, index=pd.Index(t, name='Time (min)'))
        self.chromatograms[name] = {'df': df, 'info': info, 'peaks': [],
                                    'cs_peaks': cs_peaks or [], 'hidden': False}
        if estr is not None:
            self.chromatograms[name]['estr'] = estr
        return name

    def process_file(self, path):
        """Load a single .ch, a .D folder (all its signals) or a CSV (time;signal)."""
        self._reset_view = True       # new data: the next redraw shows the whole range
        if os.path.isdir(path):
            channels = sorted(f for f in os.listdir(path) if f.lower().endswith('.ch'))
            base = os.path.splitext(os.path.basename(os.path.normpath(path)))[0]
            self._load_uv(path, base)
            self._load_instrument(path, base)
            if base in self.spectra:
                # with the DAD spectra, the recorded signals (.ch) are just a few wavelengths
                # chosen at acquisition: they are loaded on request, by default we extract from scratch
                if self.var_include_ch.get():
                    for c in channels:
                        self._load_ch(os.path.join(path, c), '%s %s' % (base, os.path.splitext(c)[0]), path, c)
                self._initial_extraction(base, path, channels)
            elif not channels:
                messagebox.showwarning("Open", "No .ch or .uv file in the folder.")
                return
            else:
                for c in channels:
                    self._load_ch(os.path.join(path, c), '%s %s' % (base, os.path.splitext(c)[0]), path, c)
        elif path.lower().endswith('.ch'):
            d = os.path.dirname(path)
            base = os.path.splitext(os.path.basename(d))[0] if d.lower().endswith('.d') else ''
            name = os.path.splitext(os.path.basename(path))[0]
            self._load_ch(path, (base + ' ' + name).strip(), d, os.path.basename(path))
        elif path.lower().endswith(('.csv', '.txt')):
            self._load_csv(path)
        else:
            messagebox.showwarning("Open", "Unsupported format: " + os.path.basename(path))

    def _load_ch(self, path, name, folder, channel):
        try:
            t, y, info = self.read_ch(path)
        except Exception as e:
            messagebox.showerror("Error", "Cannot read %s:\n%s" % (os.path.basename(path), e))
            return
        info['Source'] = path
        if folder.lower().endswith('.d'):
            self._merge_metadata(info, folder)
        cs_peaks = self._chemstation_peaks(folder, channel) if folder.lower().endswith('.d') else []
        self._add_trace(name, t, y, info, cs_peaks)

    def _load_csv(self, path):
        try:
            df = pd.read_csv(path, sep=None, engine='python', encoding='latin-1', header=None,
                             usecols=[0, 1]).apply(pd.to_numeric, errors='coerce').dropna()
        except Exception as e:
            messagebox.showerror("Error", "Cannot read %s:\n%s" % (os.path.basename(path), e))
            return
        self._add_trace(os.path.splitext(os.path.basename(path))[0], df.iloc[:, 0].to_numpy(),
                       df.iloc[:, 1].to_numpy(), {'Source': path, 'Type': 'CSV (time, signal)'})

    def load_folder_dialog(self):
        d = filedialog.askdirectory(title="Select the ChemStation .D folder")
        if d:
            self.process_file(d)
            self.refresh_view()

    def load_file_dialog(self):
        for p in filedialog.askopenfilenames(title="Open signal",
                                             filetypes=[("Signals", "*.ch *.csv *.txt"), ("All Files", "*.*")]):
            self.process_file(p)
        self.refresh_view()

    def handle_drop(self, event):
        for p in self.root.tk.splitlist(event.data):
            self.process_file(p)
        self.refresh_view()

    # ------------------------------------------------------------------ view
    def _selected(self):
        names = list(self.chromatograms)
        return [names[i] for i in self.listbox.curselection() if i < len(names)]

    def refresh_view(self):
        """Update list, table and plot after every change to the data."""
        sel = set(self._selected())
        self.listbox.delete(0, tk.END)
        for i, (n, c) in enumerate(self.chromatograms.items()):
            self.listbox.insert(tk.END, n)
            if c['hidden']:
                self.listbox.itemconfig(i, fg='grey')
            if n in sel:
                self.listbox.selection_set(i)
        if len(self.chromatograms) == 1:
            self.listbox.selection_set(0)
        self._populate_data_table()
        self._redraw()
        self._update_metadata()

    def _redraw(self):
        self._cursor_artists = []
        # ax.clear() puts the axes back in automatic mode and loses the user's zoom at every
        # Hide/Show, Update trace, peaks, etc. If the axis is not automatic (zoom, pan or set_xlim)
        # its limits are kept; with new data (self._reset_view) we start again from the whole range.
        xl, yl = self.ax.get_xlim(), self.ax.get_ylim()
        keep_x = not self.ax.get_autoscalex_on() and not self._reset_view
        keep_y = not self.ax.get_autoscaley_on() and not self._reset_view
        self._reset_view = False
        self.ax.clear()
        for n, c in self.chromatograms.items():
            if c['hidden']:
                continue
            df = c['df']
            self.ax.plot(df.index, df['mAU'], lw=1, label=n)
            col = self.ax.lines[-1].get_color()
            for p in c['peaks']:
                m = (df.index >= p['start']) & (df.index <= p['end'])
                x = df.index[m]
                base = np.interp(x, [p['start'], p['end']], [p['y_start'], p['y_end']])
                self.ax.fill_between(x, base, df['mAU'][m], color=col, alpha=0.3, lw=0)
                self.ax.annotate('%.2f' % p['rt'], (p['rt'], p['y_apex']), xytext=(0, 4),
                                 textcoords='offset points', ha='center', fontsize=8, color=col)
            if self.var_cs_peaks.get():
                for p in c['cs_peaks']:
                    y = np.interp(p['rt'], df.index, df['mAU'])
                    self.ax.plot([p['rt']], [y], marker='v', ms=4, color='red', ls='none')
        self._marker_artists = []
        for t, col in self._visible_markers():
            self._draw_marker(t, col)
        self.ax.set_xlabel("Time (min)")
        self.ax.set_ylabel("Signal (mAU)")
        if any(not l.get_label().startswith('_') for l in self.ax.lines):   # at least one labeled trace
            # fixed position: with 'best' matplotlib recomputes it at every redraw, taking the
            # cursor line into account too, and the legend would jump when the cursor passes over it
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
            self._trim_margins()
        else:
            self.ax.grid(True, linestyle=':', alpha=0.6)
        if keep_x:
            self.ax.set_xlim(xl)
        if keep_y:
            self.ax.set_ylim(yl)
        self.canvas.draw()

    def _trim_margins(self, event=None):
        """Reset the margins + tight_layout: with enlarged fonts the Y axis title
        would end up outside the canvas, and the margins are fractions (they change on resize)."""
        if not self.chromatograms or origin_style is None:
            return
        try:
            self.fig.subplots_adjust(left=0.125, right=0.9, bottom=0.11, top=0.88,
                                     wspace=0.2, hspace=0.2)
            self.fig.tight_layout()
        except Exception:
            pass

    def _populate_data_table(self):
        self.table.delete(*self.table.get_children())
        if not self.chromatograms:
            self.table['columns'] = ()
            return
        tab = pd.concat([c['df']['mAU'].rename(n) for n, c in self.chromatograms.items()], axis=1)
        tab.index = tab.index.round(4)
        tab = tab.sort_index().round(4)
        cols = ['Time (min)'] + list(tab.columns)
        self.table['columns'] = cols
        for c in cols:
            self.table.heading(c, text=c)
            self.table.column(c, width=90 if c == 'Time (min)' else 110, stretch=False)
        for t, row in zip(tab.index, tab.itertuples(index=False)):
            self.table.insert('', tk.END, values=[t] + ['' if pd.isna(v) else v for v in row])

    def _copy_data_table(self, event=None):
        sel = self.table.selection()
        if not sel:
            return 'break'
        rows = ['\t'.join(self.table['columns'])]
        rows += ['\t'.join(str(v) for v in self.table.item(i, 'values')) for i in sel]
        self.root.clipboard_clear()
        self.root.clipboard_append('\n'.join(rows))
        return 'break'

    def _select_all_table(self, event=None):
        self.table.selection_set(self.table.get_children())
        return 'break'

    def _update_metadata(self):
        self.meta.config(state='normal')
        self.meta.delete('1.0', tk.END)
        for n in self._selected():
            c = self.chromatograms[n]
            self.meta.insert(tk.END, '== %s\n' % n)
            # the most used information first (sample, method, column, injection), then the rest
            first_keys = ('Sample', 'Date', 'Method', 'Column', 'Column size', 'Inj Volume', 'Location', 'Signal')
            order = [k for k in first_keys if k in c['info']] + [k for k in c['info'] if k not in first_keys]
            for k in order:
                self.meta.insert(tk.END, '%-19s %s\n' % (k, c['info'][k]))
            df = c['df']
            self.meta.insert(tk.END, '%-19s %d (%.2f - %.2f min)\n' % (
                'Points', len(df), df.index[0], df.index[-1]))
            if c['cs_peaks']:
                self.meta.insert(tk.END, '%-19s %d\n' % ('ChemStation peaks', len(c['cs_peaks'])))
        self.meta.config(state='disabled')
        # the selected sample name (or names) also appears in the window title
        samples = []
        for n in self._selected():
            c = self.chromatograms[n]['info'].get('Sample')
            if c and c not in samples:
                samples.append(c)
        self.root.title("HPLC Manager" + (" - " + ", ".join(samples) if samples else ""))

    # ------------------------------------------------------------------ cursor
    def on_mouse_move(self, event):
        self._clear_cursor()
        if event.inaxes is not self.ax or event.xdata is None:
            return
        vis = [(n, c['df']) for n, c in self.chromatograms.items() if not c['hidden']]
        if not vis:
            return
        x = event.xdata
        rows = ['t = %.3f min' % x]
        for n, df in vis:
            rows.append('%s: %.2f mAU' % (n, np.interp(x, df.index, df['mAU'])))
        # the values box sits at the top right under the legend, anchored to its right
        # edge: its position does not depend on the cursor, only the content changes
        y_top = 0.98
        leg = self.ax.get_legend()
        if leg is not None:
            bb = leg.get_window_extent(self.canvas.get_renderer())
            y_top = self.ax.transAxes.inverted().transform((bb.x0, bb.y0))[1] - 0.015
        self._cursor_artists.append(self.ax.axvline(x, color='grey', lw=0.8, ls=':', zorder=1))
        self._cursor_artists.append(self.ax.text(
            0.98, y_top, '\n'.join(rows), transform=self.ax.transAxes, va='top', ha='right',
            fontsize=9, zorder=6, bbox=dict(boxstyle='round', fc='white', ec='grey', alpha=0.85)))
        self.canvas.draw_idle()

    def on_mouse_leave(self, event):
        self._clear_cursor()
        self.canvas.draw_idle()

    def _clear_cursor(self):
        for a in self._cursor_artists:
            try:
                a.remove()
            except Exception:
                pass
        self._cursor_artists = []

    # ------------------------------------------------------------------ list handling
    def remove_selected(self):
        for n in self._selected():
            del self.chromatograms[n]
        self.refresh_view()

    def clear_all(self):
        self._reset_view = True
        self.chromatograms.clear()
        self.spectra.clear()
        self.instrument.clear()
        self._instr_close()
        self._sp_close_window()
        self._sp_live.clear()
        self._dirty = False
        for w in self.f_tool_host.winfo_children():
            w.destroy()
        self.refresh_view()

    def hide_selected(self):
        for n in self._selected():
            self.chromatograms[n]['hidden'] = True
        self.refresh_view()

    def show_selected(self):
        for n in self._selected():
            self.chromatograms[n]['hidden'] = False
        self.refresh_view()

    def export_csv(self):
        if not self.chromatograms:
            return
        p = filedialog.asksaveasfilename(defaultextension='.csv', filetypes=[("CSV", "*.csv")])
        if not p:
            return
        tab = pd.concat([c['df']['mAU'].rename(n) for n, c in self.chromatograms.items()], axis=1)
        tab.sort_index().to_csv(p, sep=';', na_rep='', encoding='latin-1')
        self._dirty = False

    # ------------------------------------------------------------------ tool: peaks
    def _single_selection(self, title):
        sel = self._selected()
        if len(sel) != 1:
            messagebox.showinfo(title, "Select a single chromatogram in the list.")
            return None
        return sel[0]

    def open_peaks(self):
        f = self._open_tool_panel("Peaks")
        if not HAS_SCIPY:
            tk.Label(f, text="scipy not installed: peak detection not available.").pack()
            return
        self.var_prom = tk.StringVar(value='5')
        self.var_tmin = tk.StringVar(value='')
        self.var_tmax = tk.StringVar(value='')
        for r, (txt, var) in enumerate((("Min prominence (mAU):", self.var_prom),
                                        ("From (min):", self.var_tmin), ("To (min):", self.var_tmax))):
            tk.Label(f, text=txt).grid(row=r, column=0, sticky='w')
            tk.Entry(f, textvariable=var, width=8).grid(row=r, column=1, padx=4)
        tk.Button(f, text="Detect & integrate", command=self._pk_apply).grid(row=3, column=0, columnspan=2, pady=3)
        tk.Button(f, text="Clear peaks", command=self._pk_clear).grid(row=3, column=2)
        cols = ('n', 'rt', 'h', 'area', 'pct')
        self.tab_peaks = ttk.Treeview(f, columns=cols, show='headings', height=12)
        for c, t, w in (('n', '#', 30), ('rt', 'RT (min)', 70), ('h', 'Height', 65),
                        ('area', 'Area (mAU*s)', 90), ('pct', 'Area %', 60)):
            self.tab_peaks.heading(c, text=t)
            self.tab_peaks.column(c, width=w, anchor='e')
        self.tab_peaks.grid(row=4, column=0, columnspan=3, sticky='ew', pady=3)
        tk.Button(f, text="Copy table", command=self._pk_copy).grid(row=5, column=0, sticky='w')
        tk.Button(f, text="Export CSV", command=self.export_peaks).grid(row=5, column=1, columnspan=2, sticky='w')
        tk.Label(f, text="Baseline: straight line between the peak bases.\nRed triangles: ChemStation peaks.",
                 justify='left', fg='grey30').grid(row=6, column=0, columnspan=3, sticky='w')
        sel = self._selected()
        if len(sel) == 1:
            self._pk_fill(sel[0])

    def _pk_apply(self):
        n = self._single_selection("Peaks")
        if n is None:
            return
        try:
            prom = float(self.var_prom.get().replace(',', '.'))
            tmin = float(self.var_tmin.get().replace(',', '.')) if self.var_tmin.get().strip() else None
            tmax = float(self.var_tmax.get().replace(',', '.')) if self.var_tmax.get().strip() else None
        except ValueError:
            messagebox.showerror("Peaks", "Invalid numeric values.")
            return
        df = self.chromatograms[n]['df']
        t, y = df.index.to_numpy(), df['mAU'].to_numpy()
        idx, _ = find_peaks(y, prominence=prom)
        # peak bases = local minima nearest to the apex (on a slightly smoothed
        # copy, so as not to stop on noise). scipy's 'left_bases' can reach
        # back to the start of the chromatogram for a tall peak and are not suitable for integration.
        ys = savgol_filter(y, 7, 2) if len(y) > 7 else y
        peaks = []
        for i in idx:
            if (tmin is not None and t[i] < tmin) or (tmax is not None and t[i] > tmax):
                continue
            j = max(i - 3, 0) + int(np.argmax(ys[max(i - 3, 0):i + 4]))   # apex of the smoothed curve
            l = j
            while l > 0 and ys[l - 1] <= ys[l]:
                l -= 1
            r = j
            while r < len(y) - 1 and ys[r + 1] <= ys[r]:
                r += 1
            tt, yy = t[l:r + 1], y[l:r + 1]
            base = np.linspace(yy[0], yy[-1], len(yy))
            area = float(np.trapezoid(yy - base, tt)) * 60.0   # mAU*min -> mAU*s
            # RT at the vertex of the parabola through the 3 points around the apex
            rt = t[i]
            if 0 < i < len(t) - 1:
                a, b_, c = y[i - 1], y[i], y[i + 1]
                den = a - 2 * b_ + c
                if den != 0:
                    rt = t[i] + 0.5 * (a - c) / den * (t[i + 1] - t[i])
            peaks.append({'rt': float(rt), 'height': float(y[i] - base[i - l]), 'area': area,
                           'start': float(tt[0]), 'end': float(tt[-1]), 'y_start': float(yy[0]),
                           'y_end': float(yy[-1]), 'y_apex': float(y[i])})
        tot = sum(p['area'] for p in peaks) or 1.0
        for k, p in enumerate(peaks, 1):
            p['n'] = k
            p['pct'] = 100.0 * p['area'] / tot
        self.chromatograms[n]['peaks'] = peaks
        self._pk_fill(n)
        self._redraw()

    def _pk_fill(self, n):
        self.tab_peaks.delete(*self.tab_peaks.get_children())
        for p in self.chromatograms[n]['peaks']:
            self.tab_peaks.insert('', tk.END, values=(
                p['n'], '%.3f' % p['rt'], '%.2f' % p['height'], '%.2f' % p['area'], '%.2f' % p['pct']))

    def _pk_clear(self):
        n = self._single_selection("Peaks")
        if n is not None:
            self.chromatograms[n]['peaks'] = []
            self._pk_fill(n)
            self._redraw()

    def _pk_copy(self):
        rows = ['#\tRT (min)\tHeight\tArea (mAU*s)\tArea %']
        rows += ['\t'.join(str(v) for v in self.tab_peaks.item(i, 'values'))
                  for i in self.tab_peaks.get_children()]
        self.root.clipboard_clear()
        self.root.clipboard_append('\n'.join(rows))

    # ------------------------------------------------------------------ tool: smoothing, trim, normalize
    def open_smoothing(self):
        f = self._open_tool_panel("Smoothing (Savitzky-Golay)")
        self.var_sg_w = tk.StringVar(value='11')
        self.var_sg_p = tk.StringVar(value='3')
        tk.Label(f, text="Window (points, odd):").grid(row=0, column=0, sticky='w')
        tk.Entry(f, textvariable=self.var_sg_w, width=6).grid(row=0, column=1)
        tk.Label(f, text="Polynomial order:").grid(row=1, column=0, sticky='w')
        tk.Entry(f, textvariable=self.var_sg_p, width=6).grid(row=1, column=1)
        tk.Button(f, text="Apply (new trace)", command=self._sg_apply).grid(row=2, column=0, columnspan=2, pady=3)

    def _sg_apply(self):
        n = self._single_selection("Smoothing")
        if n is None:
            return
        try:
            w, p = int(self.var_sg_w.get()), int(self.var_sg_p.get())
        except ValueError:
            messagebox.showerror("Smoothing", "Invalid values.")
            return
        if w % 2 == 0 or w <= p:
            messagebox.showerror("Smoothing", "The window must be odd and larger than the order.")
            return
        df = self.chromatograms[n]['df']
        if HAS_SCIPY:
            y = savgol_filter(df['mAU'].to_numpy(), w, p)
        else:
            y = np.convolve(df['mAU'].to_numpy(), np.ones(w) / w, mode='same')
        info = dict(self.chromatograms[n]['info'], Processing='Smoothing w=%d p=%d' % (w, p))
        self._add_trace('%s_sg' % n, df.index.to_numpy(), y, info)
        self._dirty = True
        self.refresh_view()

    def open_trim(self):
        f = self._open_tool_panel("Trim")
        self.var_tr_a = tk.StringVar()
        self.var_tr_b = tk.StringVar()
        sel = self._selected()
        if len(sel) == 1:
            ix = self.chromatograms[sel[0]]['df'].index
            self.var_tr_a.set('%.2f' % ix[0])
            self.var_tr_b.set('%.2f' % ix[-1])
        tk.Label(f, text="From (min):").grid(row=0, column=0, sticky='w')
        tk.Entry(f, textvariable=self.var_tr_a, width=8).grid(row=0, column=1)
        tk.Label(f, text="To (min):").grid(row=1, column=0, sticky='w')
        tk.Entry(f, textvariable=self.var_tr_b, width=8).grid(row=1, column=1)
        tk.Button(f, text="Apply (new trace)", command=self._tr_apply).grid(row=2, column=0, columnspan=2, pady=3)

    def _tr_apply(self):
        n = self._single_selection("Trim")
        if n is None:
            return
        try:
            a, b = float(self.var_tr_a.get().replace(',', '.')), float(self.var_tr_b.get().replace(',', '.'))
        except ValueError:
            messagebox.showerror("Trim", "Invalid values.")
            return
        df = self.chromatograms[n]['df']
        sub = df[(df.index >= min(a, b)) & (df.index <= max(a, b))]
        if len(sub) < 2:
            messagebox.showerror("Trim", "Empty range.")
            return
        info = dict(self.chromatograms[n]['info'], Processing='Trim %.2f-%.2f min' % (a, b))
        self._add_trace('%s_trim' % n, sub.index.to_numpy(), sub['mAU'].to_numpy(), info)
        self._dirty = True
        self.refresh_view()

    def normalize(self):
        sel = self._selected()
        if not sel:
            messagebox.showinfo("Normalize", "Select at least one chromatogram.")
            return
        for n in sel:
            df = self.chromatograms[n]['df']
            m = float(np.abs(df['mAU']).max()) or 1.0
            info = dict(self.chromatograms[n]['info'], Processing='Normalized to max = 100')
            self._add_trace('%s_norm' % n, df.index.to_numpy(), df['mAU'].to_numpy() * 100.0 / m, info)
        self._dirty = True
        self.refresh_view()

    # ------------------------------------------------------------------ tool: DAD spectra
    def _initial_extraction(self, base, folder, channels):
        """After opening a .D with spectra: extract one trace (wavelength of the first
        recorded signal, 4 nm band, reference off) and open the control panel."""
        lam = 254.0
        for c in channels[:1]:
            try:
                sig = self.read_ch(os.path.join(folder, c))[2].get('Signal', '')
                m = re.search(r'Sig=\s*([\d.]+)', sig)
                if m:
                    lam = float(m.group(1))
            except Exception:
                pass
        self.var_sp_set = tk.StringVar(value=base)
        self._sp_initial_values(lam)
        self.open_spectra()
        self._sp_extract()

    def _sp_initial_values(self, lam):
        self.var_sp_l = tk.StringVar(value='%g' % lam)
        self.var_sp_b = tk.StringVar(value='4')
        self.var_sp_rl = tk.StringVar(value='360')
        self.var_sp_rb = tk.StringVar(value='100')
        self.var_sp_use_ref = tk.BooleanVar(value=False)     # reference off by default
        self.var_sp_t = tk.StringVar(value='4.35')

    def open_spectra(self):
        f = self._open_tool_panel("Detector signal (DAD)")
        if not self.spectra:
            tk.Label(f, text="No spectra loaded: open a .D folder\nwith the .uv file.",
                     justify='left').pack(anchor='w')
            return
        names = list(self.spectra)
        if not hasattr(self, 'var_sp_l'):
            self.var_sp_set = tk.StringVar(value=names[0])
            self._sp_initial_values(254.0)
        if self.var_sp_set.get() not in self.spectra:
            self.var_sp_set.set(names[0])
        cb = ttk.Combobox(f, textvariable=self.var_sp_set, values=names, state='readonly', width=18)
        cb.grid(row=0, column=0, columnspan=3, sticky='w')
        d = self.spectra[self.var_sp_set.get()]
        tk.Label(f, text="%s" % d['info']['Spectra'], fg='grey30').grid(row=1, column=0, columnspan=3, sticky='w')
        tk.Label(f, text="Wavelength (nm)").grid(row=2, column=0, sticky='w')
        e_l = tk.Entry(f, textvariable=self.var_sp_l, width=7)
        e_l.grid(row=2, column=1)
        tk.Label(f, text="Bandwidth (nm)").grid(row=3, column=0, sticky='w')
        e_b = tk.Entry(f, textvariable=self.var_sp_b, width=7)
        e_b.grid(row=3, column=1)
        self.chk_ref = tk.Checkbutton(f, text="Reference (nm / bw)", variable=self.var_sp_use_ref,
                                      command=self._sp_ref_state)
        self.chk_ref.grid(row=4, column=0, sticky='w')
        self.e_rl = tk.Entry(f, textvariable=self.var_sp_rl, width=7)
        self.e_rl.grid(row=4, column=1)
        self.e_rb = tk.Entry(f, textvariable=self.var_sp_rb, width=7)
        self.e_rb.grid(row=4, column=2)
        self._sp_ref_state()
        tk.Button(f, text="Update trace", command=lambda: self._sp_extract(False)).grid(
            row=5, column=0, columnspan=2, pady=3, sticky='ew')
        tk.Button(f, text="Add as new", command=lambda: self._sp_extract(True)).grid(
            row=5, column=2, pady=3)
        for e in (e_l, e_b, self.e_rl, self.e_rb):
            e.bind('<Return>', lambda ev: self._sp_extract(False))
        tk.Button(f, text="Several wavelengths...", command=self.open_multi_extraction).grid(
            row=6, column=0, columnspan=3, pady=1, sticky='ew')
        tk.Label(f, text="Spectrum at (min)").grid(row=7, column=0, sticky='w')
        tk.Entry(f, textvariable=self.var_sp_t, width=7).grid(row=7, column=1)
        tk.Button(f, text="Show", command=self._sp_show).grid(row=7, column=2)
        tk.Checkbutton(f, text="Click on plot shows the spectrum\n(Shift+click works even if unticked)", justify="left",
                       variable=self.var_click_spec).grid(row=8, column=0, columnspan=3, sticky='w')
        tk.Button(f, text="Time-wavelength map", command=self._sp_map).grid(
            row=9, column=0, columnspan=3, pady=3)
        tk.Button(f, text="Export picked spectra (CSV)", command=self.export_spectra).grid(
            row=10, column=0, columnspan=3, pady=1)

    def _sp_ref_state(self):
        state = 'normal' if self.var_sp_use_ref.get() else 'disabled'
        self.e_rl.config(state=state)
        self.e_rb.config(state=state)

    def _sp_dataset(self):
        return self.spectra.get(self.var_sp_set.get()) if hasattr(self, 'var_sp_set') else None

    def _sp_extract(self, new=False):
        """Extract the chromatogram at the wavelength/bandwidth/reference of the panel.

        By default it modifies the extracted trace selected in the list (replacing it in the same
        position); if there is none, the last one extracted; with new=True it adds a separate one, which
        becomes the selected one and therefore the one the next Update modifies."""
        d = self._sp_dataset()
        if d is None:
            return
        try:
            l, b = float(self.var_sp_l.get().replace(',', '.')), float(self.var_sp_b.get().replace(',', '.'))
            if b <= 0:
                raise ValueError
            ref = rb = None
            if self.var_sp_use_ref.get():
                ref, rb = float(self.var_sp_rl.get().replace(',', '.')), float(self.var_sp_rb.get().replace(',', '.'))
                if rb <= 0:
                    raise ValueError
        except ValueError:
            messagebox.showerror("Spectra", "Invalid numeric values.")
            return
        try:
            y = self.extract_wavelength(d['wl'], d['S'], l, b, ref, rb)
        except ValueError as e:
            messagebox.showerror("Spectra", str(e))
            return
        ds = self.var_sp_set.get()
        label = self._extraction_label(l, b, ref, rb)
        info = self._extraction_info(d, ds, label)
        estr = {'ds': ds, 'l': l, 'b': b, 'ref': ref, 'rb': rb}
        old = None if new else self._sp_target(ds)
        if old is None:
            name = self._add_trace('%s %s' % (ds, label), d['t'], y, info, estr=estr)
        else:
            # replace the selected trace, keeping its position in the list
            name = '%s %s' % (ds, label)
            if name != old:
                name = self._free_name(name)
            df = pd.DataFrame({'mAU': y}, index=pd.Index(d['t'], name='Time (min)'))
            alt_label = dict(self.chromatograms[old], df=df, info=info, peaks=[], cs_peaks=[],
                         hidden=False, estr=estr)
            self.chromatograms = {(name if k == old else k): (alt_label if k == old else v)
                                  for k, v in self.chromatograms.items()}
        self._sp_live[ds] = name          # the last one touched: used if there is no useful selection
        if new:
            self._dirty = True            # a 'live' trace is always recomputed: no need to export it
        self.refresh_view()
        self._select_name(name)        # the trace just extracted/modified becomes the active one

    def _sp_target(self, ds):
        """Trace that 'Update trace' must modify: the one selected in the list if it is a trace
        extracted from dataset `ds`; otherwise the last extracted/modified one; None if there is none."""
        sel = self._selected()
        if len(sel) == 1:
            e = self.chromatograms[sel[0]].get('estr')
            if e and e['ds'] == ds:
                return sel[0]
        last = self._sp_live.get(ds)
        return last if last in self.chromatograms else None

    def _select_name(self, name):
        """Select trace `name` in the list (and update the metadata)."""
        names = list(self.chromatograms)
        if name not in names:
            return
        i = names.index(name)
        self.listbox.selection_clear(0, tk.END)
        self.listbox.selection_set(i)
        self.listbox.see(i)
        self._update_metadata()

    def _on_selection(self, event=None):
        """Selection change in the list: metadata, and for an extracted trace its parameters
        (wavelength, bandwidth, reference) go back into the panel fields, ready to be modified."""
        self._update_metadata()
        sel = self._selected()
        if len(sel) != 1 or not hasattr(self, 'var_sp_l'):
            return
        e = self.chromatograms[sel[0]].get('estr')
        if not e or e['ds'] not in self.spectra:
            return
        self.var_sp_set.set(e['ds'])
        self.var_sp_l.set('%g' % e['l'])
        self.var_sp_b.set('%g' % e['b'])
        self.var_sp_use_ref.set(e['ref'] is not None)
        if e['ref'] is not None:
            self.var_sp_rl.set('%g' % e['ref'])
            self.var_sp_rb.set('%g' % e['rb'])
        try:
            self._sp_ref_state()
        except tk.TclError:           # panel not (or no longer) built
            pass

    def _draw_marker(self, t, col):
        """Dashed line + time label on the chromatogram (without redrawing the rest,
        so that the user's zoom is not lost)."""
        self._marker_artists.append(self.ax.axvline(t, color=col, lw=0.9, ls='--', zorder=1, gid='marker'))
        self._marker_artists.append(self.ax.annotate(
            '%.2f' % t, (t, 1.0), xycoords=('data', 'axes fraction'), xytext=(2, -3),
            textcoords='offset points', rotation=90, va='top', ha='left', fontsize=8, color=col, gid='marker'))

    @staticmethod
    def _extraction_label(l, b, ref, rb):
        return '%g,%g' % (l, b) + (' ref %g,%g' % (ref, rb) if ref is not None else ' no ref')

    @staticmethod
    def _extraction_info(d, ds, label):
        info = {k: v for k, v in d['info'].items() if k not in ('Version', 'Type', 'Spectra')}
        info.update({'Source': ds + ' (dad .uv)', 'Signal': 'Extracted Sig=' + label})
        return info

    def extract_series(self, ds, specs):
        """Extract several traces from dataset `ds` in one go and add them to the list.

        `specs`: list of (wavelength, bandwidth, reference or None, reference bandwidth).
        All or nothing: if a row is not valid (non-positive bandwidth, outside the acquired spectrum)
        it raises ValueError with the row number and adds nothing. Returns the trace names."""
        d = self.spectra[ds]
        computed = []
        for i, (l, b, ref, rb) in enumerate(specs, 1):
            if b <= 0 or (ref is not None and rb <= 0):
                raise ValueError("Row %d: the bandwidth must be positive." % i)
            try:
                y = self.extract_wavelength(d['wl'], d['S'], l, b, ref, rb)
            except ValueError as e:
                raise ValueError("Row %d: %s" % (i, e))
            computed.append((self._extraction_label(l, b, ref, rb), y,
                              {'ds': ds, 'l': l, 'b': b, 'ref': ref, 'rb': rb}))
        names = []
        for label, y, estr in computed:
            names.append(self._add_trace('%s %s' % (ds, label), d['t'], y,
                                       self._extraction_info(d, ds, label), estr=estr))
        return names

    def open_multi_extraction(self):
        """Dialog: how many traces you want, and for each one wavelength, bandwidth and optional
        reference (each with its own, or none). The created traces stay in the list as
        'Add as new'; the panel's 'live' trace does not change."""
        if not self.spectra:
            messagebox.showinfo("Several wavelengths", "No spectra loaded: open a .D folder with the .uv.")
            return
        ds = self.var_sp_set.get() if hasattr(self, 'var_sp_set') and self.var_sp_set.get() in self.spectra \
            else next(iter(self.spectra))
        top = tk.Toplevel(self.root)
        top.title("Extract several wavelengths")
        top.transient(self.root)
        tk.Label(top, text="Dataset: %s" % ds, fg='grey30').grid(row=0, column=0, columnspan=6, sticky='w', padx=8, pady=(8, 0))
        tk.Label(top, text="Number of traces:").grid(row=1, column=0, columnspan=2, sticky='w', padx=8, pady=6)
        var_n = tk.IntVar(value=len(self._sp_multi_last) or 3)
        for c, t in enumerate(("#", "Wavelength (nm)", "Bandwidth (nm)", "Reference", "Ref. wavelength", "Ref. bandwidth")):
            tk.Label(top, text=t, font=('Segoe UI', 9, 'bold')).grid(row=2, column=c, padx=4)
        rows = []
        base_b = self.var_sp_b.get() if hasattr(self, 'var_sp_b') else '4'
        base_rl = self.var_sp_rl.get() if hasattr(self, 'var_sp_rl') else '360'
        base_rb = self.var_sp_rb.get() if hasattr(self, 'var_sp_rb') else '100'
        base_use = self.var_sp_use_ref.get() if hasattr(self, 'var_sp_use_ref') else False

        def update_state(r):
            st = 'normal' if r['use'].get() else 'disabled'
            r['e_rl'].config(state=st)
            r['e_rb'].config(state=st)

        def add_row():
            i = len(rows)
            previous = self._sp_multi_last[i] if i < len(self._sp_multi_last) else None
            r = {'l': tk.StringVar(value=previous['l'] if previous else ''),
                 'b': tk.StringVar(value=previous['b'] if previous else base_b),
                 'use': tk.BooleanVar(value=previous['use'] if previous else base_use),
                 'rl': tk.StringVar(value=previous['rl'] if previous else base_rl),
                 'rb': tk.StringVar(value=previous['rb'] if previous else base_rb)}
            r['w'] = [tk.Label(top, text=str(i + 1)),
                      tk.Entry(top, textvariable=r['l'], width=10),
                      tk.Entry(top, textvariable=r['b'], width=10),
                      tk.Checkbutton(top, variable=r['use'], command=lambda r=r: update_state(r))]
            r['e_rl'] = tk.Entry(top, textvariable=r['rl'], width=10)
            r['e_rb'] = tk.Entry(top, textvariable=r['rb'], width=10)
            r['w'] += [r['e_rl'], r['e_rb']]
            for c, w in enumerate(r['w']):
                w.grid(row=3 + i, column=c, padx=4, pady=1)
            update_state(r)
            rows.append(r)

        def set_n(*_):
            try:
                n = max(1, min(30, int(var_n.get())))
            except (tk.TclError, ValueError):
                return
            while len(rows) < n:
                add_row()
            while len(rows) > n:
                for w in rows.pop()['w']:
                    w.destroy()
            top.update_idletasks()

        spin = tk.Spinbox(top, from_=1, to=30, width=4, textvariable=var_n, command=set_n)
        spin.grid(row=1, column=2, sticky='w')
        spin.bind('<Return>', set_n)
        spin.bind('<FocusOut>', set_n)
        frame_btn = tk.Frame(top)
        frame_btn.grid(row=100, column=0, columnspan=6, pady=8)

        def read_it():
            set_n()
            specs, saved = [], []
            for i, r in enumerate(rows, 1):
                try:
                    l, b = float(r['l'].get().replace(',', '.')), float(r['b'].get().replace(',', '.'))
                    ref = rb = None
                    if r['use'].get():
                        ref, rb = float(r['rl'].get().replace(',', '.')), float(r['rb'].get().replace(',', '.'))
                except ValueError:
                    raise ValueError("Row %d: invalid numeric values." % i)
                specs.append((l, b, ref, rb))
                saved.append({k: (r[k].get()) for k in ('l', 'b', 'use', 'rl', 'rb')})
            return specs, saved

        def extract():
            try:
                specs, saved = read_it()
                self.extract_series(ds, specs)
            except ValueError as e:
                messagebox.showerror("Several wavelengths", str(e), parent=top)
                return
            self._sp_multi_last = saved
            self._dirty = True
            top.destroy()
            self.refresh_view()

        tk.Button(frame_btn, text="Extract", width=12, command=extract).pack(side=tk.LEFT, padx=6)
        tk.Button(frame_btn, text="Cancel", width=12, command=top.destroy).pack(side=tk.LEFT, padx=6)
        top.bind('<Escape>', lambda e: top.destroy())
        set_n()
        top._extract = extract          # for the tests
        top._rows = rows
        top._var_n = var_n
        top._set_n = set_n
        return top

    def _visible_markers(self):
        return [(k['t'], k['color']) for k in self._chosen_spectra if k.get('visible', True)]

    def _update_markers(self):
        """Redraw the dashed lines of the selected spectra only, without touching the zoom."""
        for a in self._marker_artists:
            try:
                a.remove()
            except Exception:
                pass
        self._marker_artists = []
        for t, col in self._visible_markers():
            self._draw_marker(t, col)
        self.canvas.draw_idle()

    def _sp_close_window(self):
        """Close the spectra window (if open) and clear the picked spectra."""
        w = self._spec_win
        if w is not None:
            try:
                w.destroy()
            except Exception:
                pass
        self._chosen_spectra = []
        self._spec_win = None

    def _sp_window(self):
        """Window (Toplevel) with the spectra plot and, on the right, the list of the picked ones
        (ticked = visible); created on first use."""
        if self._spec_win is not None and self._spec_win.winfo_exists():
            return self._spec_win
        w = tk.Toplevel(self.root)
        w.title("DAD spectra")
        w.geometry("920x500")
        # list on the right (fixed width, scrollable), plot on the left
        w.f_picked = tk.Frame(w, width=200)
        w.f_picked.pack(side=tk.RIGHT, fill=tk.Y)
        w.f_picked.pack_propagate(False)
        tk.Label(w.f_picked, text="Picked spectra", font=('Segoe UI', 10, 'bold')).pack(anchor='w', padx=6, pady=(6, 0))
        tk.Label(w.f_picked, text="untick = hide, \u2715 = remove", fg='grey30').pack(anchor='w', padx=6)
        tk.Button(w.f_picked, text="Clear all", command=self._sp_clear).pack(side=tk.BOTTOM, pady=6)
        cv = tk.Canvas(w.f_picked, highlightthickness=0)
        sb = ttk.Scrollbar(w.f_picked, orient='vertical', command=cv.yview)
        cv.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        cv.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(4, 0))
        w.f_rows = tk.Frame(cv)
        cv.create_window((0, 0), window=w.f_rows, anchor='nw')
        w.f_rows.bind('<Configure>', lambda e: cv.configure(scrollregion=cv.bbox('all')))
        w.f_spec_plot = tk.Frame(w)
        w.f_spec_plot.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        w.fig = Figure()
        w.ax = w.fig.add_subplot(111)
        w.canvas = FigureCanvasTkAgg(w.fig, master=w.f_spec_plot)
        NavigationToolbar2Tk(w.canvas, w.f_spec_plot).update()
        w.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        w.ax.set_xlabel("Wavelength (nm)")
        w.ax.set_ylabel("Absorbance (mAU)")
        w.bind('<Destroy>', lambda e: self._sp_window_closed(e, w))
        self._spec_win = w
        return w

    def _sp_window_closed(self, event, w):
        """When the spectra window closes, the lines on the chromatogram disappear too."""
        if event.widget is not w:
            return
        self._chosen_spectra = []
        self._spec_win = None
        try:
            self._update_markers()
        except Exception:
            pass

    def _sp_update_legend(self):
        w = self._spec_win
        if w is None:
            return
        visible_ones = [k['line'] for k in self._chosen_spectra if k.get('visible', True)]
        if visible_ones:
            w.ax.legend(handles=visible_ones, loc='upper right')
        elif w.ax.get_legend() is not None:
            w.ax.get_legend().remove()
        w.ax.relim(visible_only=True)
        w.ax.autoscale_view()       # respects the zoom: acts only on the axes that are still automatic
        w.canvas.draw_idle()

    def _sp_list_row(self, w, k):
        """Add to the window's list the row for spectrum `k` (tick + remove)."""
        row = tk.Frame(w.f_rows)
        row.pack(fill=tk.X, anchor='w')
        var = tk.BooleanVar(value=k.get('visible', True))
        tk.Checkbutton(row, text='%.3f min' % k['t'], variable=var, fg=k['color'], anchor='w',
                       font=('Segoe UI', 9, 'bold'), command=lambda: self._sp_toggle(k, var.get())
                       ).pack(side=tk.LEFT)
        tk.Button(row, text='\u2715', width=2, relief=tk.FLAT,
                  command=lambda: self._sp_remove(k)).pack(side=tk.RIGHT, padx=2)
        k['row'] = row
        k['var'] = var

    def _sp_toggle(self, k, visible):
        k['visible'] = bool(visible)
        k['line'].set_visible(bool(visible))
        self._sp_update_legend()
        self._update_markers()

    def _sp_remove(self, k):
        if k not in self._chosen_spectra:
            return
        self._chosen_spectra.remove(k)
        try:
            k['line'].remove()
            k['row'].destroy()
        except Exception:
            pass
        self._sp_update_legend()
        self._update_markers()

    def _sp_clear(self):
        for k in list(self._chosen_spectra):
            self._sp_remove(k)

    def _sp_spectrum_at(self, t_min):
        """Add to the spectra window the spectrum nearest to t_min."""
        d = self._sp_dataset()
        if d is None:
            return
        k = int(np.argmin(np.abs(d['t'] - t_min)))
        w = self._sp_window()
        colors = matplotlib.rcParams['axes.prop_cycle'].by_key()['color']
        color = colors[self._sp_color_n % len(colors)]
        self._sp_color_n += 1
        t = float(d['t'][k])
        line, = w.ax.plot(d['wl'], d['S'][k], lw=1, color=color, label='%.3f min' % t)
        chosen = {'dataset': self.var_sp_set.get(), 't': t, 'wl': d['wl'].copy(), 'y': d['S'][k].copy(),
                  'color': color, 'visible': True, 'line': line}
        self._chosen_spectra.append(chosen)
        self._sp_list_row(w, chosen)
        self._sp_update_legend()
        self._update_markers()       # dashed line on the chromatogram, same color
        w.lift()

    def _sp_show(self):
        try:
            self._sp_spectrum_at(float(self.var_sp_t.get().replace(',', '.')))
        except ValueError:
            messagebox.showerror("Spectra", "Invalid time.")

    def _press_spectrum(self, event):
        # remember where the button was pressed: needed to tell a click from a drag
        self._press_pos = (event.x, event.y) if event.button == 1 else None

    def _click_spectrum(self, event):
        """Click on the chromatogram = spectrum at that time.

        The decision is taken on button release: it is a click if the mouse moved less than 5 pixels. This
        way it also works with the zoom (magnifier) or pan tool active: those tools work only
        by dragging, and a drag (zoom on a rectangle, panning) is not a click.
        Shift+click counts even with the checkbox off."""
        press, self._press_pos = getattr(self, '_press_pos', None), None
        if (press is None or event.inaxes is not self.ax or event.button != 1
                or event.xdata is None or not self.spectra):
            return
        if np.hypot(event.x - press[0], event.y - press[1]) >= 5:
            return
        shift = event.key is not None and 'shift' in str(event.key)
        if not shift and not self.var_click_spec.get():
            return
        self.var_sp_t.set('%.3f' % event.xdata)
        self._sp_spectrum_at(event.xdata)

    def _sp_map(self):
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

    # ------------------------------------------------------------------ export and session
    def peaks_table(self):
        """DataFrame with the peaks of all traces: those computed here (Source = HPLCManager) and
        those written by ChemStation in the REPORTnn.CSV (Source = ChemStation), one per row."""
        rows = []
        for n, c in self.chromatograms.items():
            for p in c['peaks']:
                rows.append({'Trace': n, 'Source': 'HPLCManager', 'Peak': p['n'], 'RT (min)': p['rt'],
                              'Height (mAU)': p['height'], 'Area (mAU*s)': p['area'],
                              'Area %': p['pct'], 'Start (min)': p['start'], 'End (min)': p['end'],
                              'Width (min)': np.nan, 'Type': ''})
            for p in c['cs_peaks']:
                rows.append({'Trace': n, 'Source': 'ChemStation', 'Peak': p['n'], 'RT (min)': p['rt'],
                              'Height (mAU)': p['height'], 'Area (mAU*s)': p['area'],
                              'Area %': p['area_pct'], 'Start (min)': np.nan, 'End (min)': np.nan,
                              'Width (min)': p['width'], 'Type': p['kind']})
        return pd.DataFrame(rows, columns=['Trace', 'Source', 'Peak', 'RT (min)', 'Height (mAU)',
                                            'Area (mAU*s)', 'Area %', 'Start (min)', 'End (min)',
                                            'Width (min)', 'Type'])

    def spectra_table(self):
        """DataFrame of the picked and selected (ticked) spectra in the spectra window:
        wavelength + one column per spectrum, headed with dataset and time."""
        chosen = [k for k in self._chosen_spectra if k.get('visible', True)]
        if not chosen:
            return pd.DataFrame()
        columns = []
        for k in chosen:
            name = '%s t=%.3f min' % (k['dataset'], k['t'])
            columns.append(pd.Series(k['y'], index=pd.Index(k['wl'], name='Wavelength (nm)'), name=name))
        return pd.concat(columns, axis=1).sort_index()

    def _choose_save_file(self, title, extension, types):
        return filedialog.asksaveasfilename(title=title, initialdir=self.current_dir,
                                            defaultextension=extension, filetypes=types)

    def _write_csv(self, df, path, title, index_pos=True):
        try:
            df.to_csv(path, sep=';', na_rep='', encoding='latin-1', index=index_pos)
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror(title, "Save failed:\n%s" % e)
            return False
        self.current_dir = os.path.dirname(path)
        return True

    def export_peaks(self):
        df = self.peaks_table()
        if df.empty:
            messagebox.showinfo("Export Peaks", "No peaks: use Tools > Peaks (Detect & integrate) "
                                "or open a .D folder with the recorded channels (ChemStation peaks).")
            return
        path = self._choose_save_file("Export Peaks", '.csv', [("CSV", "*.csv")])
        if path and self._write_csv(df, path, "Export Peaks", index_pos=False):
            messagebox.showinfo("Export Peaks", "%d peaks saved:\n%s" % (len(df), path))

    def export_spectra(self):
        df = self.spectra_table()
        if df.empty:
            messagebox.showinfo("Export Spectra", "No spectrum picked: click on the chromatogram "
                                "(Tools > Spectra) or use 'Spectrum at (min)' and Show.")
            return
        path = self._choose_save_file("Export Spectra", '.csv', [("CSV", "*.csv")])
        if path and self._write_csv(df, path, "Export Spectra"):
            messagebox.showinfo("Export Spectra", "%d spectra saved:\n%s" % (df.shape[1], path))

    def export_full_dad(self):
        """The whole DAD cube of the selected dataset: one row per time, one column per wavelength
        (with 13493 spectra x 201 wavelengths the file is about 25 MB)."""
        d = self._sp_dataset()
        if d is None:
            messagebox.showinfo("Export Full DAD Data", "No DAD data loaded (open a .D folder with the .uv).")
            return
        path = self._choose_save_file("Export Full DAD Data", '.csv', [("CSV", "*.csv")])
        if not path:
            return
        df = pd.DataFrame(d['S'], index=pd.Index(d['t'], name='Time (min)'),
                          columns=['%g' % w for w in d['wl']])
        self.root.config(cursor='watch')
        self.root.update_idletasks()
        try:
            ok = self._write_csv(df, path, "Export Full DAD Data")
        finally:
            self.root.config(cursor='')
        if ok:
            messagebox.showinfo("Export Full DAD Data", "%d spectra x %d wavelengths saved:\n%s"
                                % (df.shape[0], df.shape[1], path))

    def _session_data(self):
        return {'version': 1, 'chromatograms': self.chromatograms, 'spectra': self.spectra,
                'sp_live': dict(self._sp_live), 'include_ch': bool(self.var_include_ch.get()),
                'instrument': self.instrument}

    def _apply_session(self, data):
        self.chromatograms = data['chromatograms']
        self.spectra = data.get('spectra', {})
        self.instrument = data.get('instrument', {})
        self._instr_close()
        self._sp_live = dict(data.get('sp_live', {}))
        self.var_include_ch.set(bool(data.get('include_ch', False)))
        self._sp_close_window()
        self._reset_view = True
        for w in self.f_tool_host.winfo_children():
            w.destroy()
        for attr in ('var_sp_set', 'var_sp_l'):          # the Spectra panel is rebuilt from scratch
            if hasattr(self, attr):
                delattr(self, attr)
        self.refresh_view()
        if self.spectra:
            self.open_spectra()

    def save_session(self):
        """Save traces (including derived ones), peaks and DAD spectra in a single file, to
        resume the work later without redoing the steps."""
        if not self.chromatograms:
            messagebox.showwarning("Save Session", "No chromatogram loaded.")
            return
        path = self._choose_save_file("Save Session", '.hplcsession',
                                             [("HPLC Manager Session", "*.hplcsession"), ("All Files", "*.*")])
        if not path:
            return
        try:
            with gzip.open(path, 'wb') as f:
                pickle.dump(self._session_data(), f)
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror("Save Session", "Save failed:\n%s" % e)
            return
        self.current_dir = os.path.dirname(path)
        self._dirty = False
        messagebox.showinfo("Save Session", "Session saved:\n%s" % path)

    def open_session(self):
        """Load a saved session, replacing everything currently loaded."""
        if self._dirty and not messagebox.askyesno(
                "Open Session", "There are derived traces that were neither exported nor saved. Discard them and open the session?"):
            return
        path = filedialog.askopenfilename(initialdir=self.current_dir, title="Open Session",
                                          filetypes=[("HPLC Manager Session", "*.hplcsession"), ("All Files", "*.*")])
        if not path:
            return
        try:
            with gzip.open(path, 'rb') as f:
                data = pickle.load(f)
            if not isinstance(data, dict) or 'chromatograms' not in data:
                raise ValueError("not an HPLC Manager session")
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror("Open Session", "Opening the session failed:\n%s" % e)
            return
        self.current_dir = os.path.dirname(path)
        self._apply_session(data)
        self._dirty = False

    def _copy_figure(self):
        """Copy of the current figure for saving/editor: without the cursor and without the lines
        of the picked spectra, restyled Origin 'single' (4:3) instead of the panel size."""
        self._clear_cursor()
        fig = pickle.loads(pickle.dumps(self.fig))
        for ax in fig.axes:
            for a in list(ax.lines) + list(ax.texts):
                if a.get_gid() == 'marker':
                    a.remove()
        if origin_style is not None and fig.axes:
            origin_style.applica_stile_origin(fig.axes[0], fig, set_size=True, preset='single')
        return fig

    def save_figure_image(self):
        if not self.chromatograms:
            messagebox.showwarning("Save Figure Image", "No chromatogram loaded.")
            return
        path = self._choose_save_file("Save Figure Image", '.png', [
            ("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg"), ("All Files", "*.*")])
        if not path:
            return
        try:
            self._copy_figure().savefig(path, dpi=300)
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror("Save Figure Image", "Save failed:\n%s" % e)
            return
        self.current_dir = os.path.dirname(path)
        messagebox.showinfo("Save Figure Image", "Figure saved:\n%s" % path)

    def save_figure_pickle(self):
        """Save the current figure as a matplotlib pickle: reopenable with plot_editor.pyw as a
        live Figure/Axes object (not a raster)."""
        if not self.chromatograms:
            messagebox.showwarning("Save Figure", "No chromatogram loaded.")
            return
        path = self._choose_save_file("Save Figure", '.fig.pickle', [
            ("Matplotlib Figure (pickle)", "*.pickle *.pkl"), ("All Files", "*.*")])
        if not path:
            return
        try:
            with open(path, 'wb') as f:
                pickle.dump(self._copy_figure(), f)
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror("Save Figure", "Save failed:\n%s" % e)
            return
        self.current_dir = os.path.dirname(path)
        messagebox.showinfo("Save Figure", "Figure saved:\n%s" % path)

    def _load_plot_editor(self):
        """Import plot_editor (only once): first a local copy, then the sibling repo PlotStyleKit."""
        if self._pe_module is None:
            import importlib.util
            here = os.path.dirname(os.path.abspath(__file__))
            candidates = [os.path.join(here, 'plot_editor.pyw'),
                         os.path.join(here, '..', 'PlotStyleKit', 'plot_editor.pyw')]
            path = next((c for c in candidates if os.path.isfile(c)), None)
            if path is None:
                raise FileNotFoundError(
                    "plot_editor.pyw not found (PlotStyleKit repo missing next to this project)")
            spec = importlib.util.spec_from_file_location('plot_editor', path)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            self._pe_module = mod
        return self._pe_module

    def open_figure_editor(self):
        """Open the PlotStyleKit Plot Editor on the current figure (an independent copy)."""
        if not self.chromatograms:
            messagebox.showwarning("Edit Figure", "No chromatogram loaded.")
            return
        try:
            pe = self._load_plot_editor()
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror("Edit Figure", "plot_editor.pyw not available:\n%s" % e)
            return
        try:
            fig = self._copy_figure()
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror("Edit Figure", "Cannot duplicate the figure:\n%s" % e)
            return
        top = tk.Toplevel(self.root)
        top.geometry("1300x820")
        editor = pe.PlotEditor(top)
        editor.carica_figura(fig, title="current figure")

    # ------------------------------------------------------------------ instrument curves
    INSTRUMENT_GROUPS = (('bar', 'Pressure (bar)'), ('%', 'Solvent composition (%)'),
                    ('ml/min', 'Flow (ml/min)'), ('\xb0C', 'Temperature (\xb0C)'))

    def _instr_close(self):
        w = self._instr_win
        self._instr_win = None
        if w is not None:
            try:
                w.destroy()
            except Exception:
                pass

    def _instr_dataset(self):
        """Dataset (.D folder) whose curves to show: the one chosen in the window, or the first."""
        v = getattr(self, 'var_instr_set', None)
        if v is not None and v.get() in self.instrument and self.instrument[v.get()]['signals']:
            return v.get()
        for k, d in self.instrument.items():
            if d['signals']:
                return k
        return None

    def instrument_table(self, ds):
        """DataFrame of the curves of dataset `ds` on the time scale of the fastest signal
        (pressure, flow and solvents at 0.3 s); the slower signals (temperature, 1 s) are
        linearly interpolated onto that scale."""
        seg = self.instrument[ds]['signals']
        base = min(seg.values(), key=lambda d: d['t'][1] - d['t'][0] if len(d['t']) > 1 else 1e9)
        out = pd.DataFrame(index=pd.Index(np.round(base['t'], 6), name='Time (min)'))
        for title, d in seg.items():
            out['%s [%s]' % (title, d['unit'])] = np.interp(base['t'], d['t'], d['y'])
        return out

    def open_instrument_curves(self):
        ds = self._instr_dataset()
        if ds is None:
            messagebox.showinfo("Instrument curves", "No instrument curves: open a .D folder "
                                "that contains LCDIAG.REG.")
            return
        self._instr_close()
        w = tk.Toplevel(self.root)
        w.title("Instrument curves - %s" % ds)
        w.geometry("900x720")
        top = tk.Frame(w)
        top.pack(fill=tk.X, padx=6, pady=4)
        self.var_instr_set = tk.StringVar(value=ds)
        names = [k for k, d in self.instrument.items() if d['signals']]
        if len(names) > 1:
            cb = ttk.Combobox(top, textvariable=self.var_instr_set, values=names, state='readonly', width=24)
            cb.pack(side=tk.LEFT)
            cb.bind('<<ComboboxSelected>>', lambda e: self._instr_draw())
        tk.Button(top, text="Export CSV", command=self.export_instrument_curves).pack(side=tk.RIGHT)
        w.fig = Figure()
        w.canvas = FigureCanvasTkAgg(w.fig, master=w)
        NavigationToolbar2Tk(w.canvas, w).update()
        w.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self._instr_win = w
        w.bind('<Destroy>', lambda e: setattr(self, '_instr_win', None) if e.widget is w else None)
        self._instr_draw()
        return w

    def _instr_draw(self):
        w = self._instr_win
        ds = self._instr_dataset()
        if w is None or ds is None:
            return
        w.title("Instrument curves - %s" % ds)
        seg = self.instrument[ds]['signals']
        groups = []
        for unit, title in self.INSTRUMENT_GROUPS:
            members = [(t, d) for t, d in seg.items() if d['unit'] == unit]
            if unit == '%':      # only the solvents that do not stay at zero
                members = [(t, d) for t, d in members if np.any(d['y'] > 0)]
            if members:
                groups.append((title, members))
        w.fig.clear()
        axes_list = w.fig.subplots(len(groups), 1, sharex=True) if groups else []
        axes_list = np.atleast_1d(axes_list)
        for ax, (title, members) in zip(axes_list, groups):
            for t, d in members:
                label = t.split(', ', 1)[-1]
                ax.plot(d['t'], d['y'], lw=1, label=label)
            ax.set_ylabel(title, fontsize=9)
            ax.tick_params(labelsize=8)
            if len(members) > 1:
                ax.legend(fontsize=8, loc='upper right')
        if len(axes_list):
            axes_list[-1].set_xlabel("Time (min)")
        w.fig.tight_layout()
        w.canvas.draw()

    def export_instrument_curves(self):
        ds = self._instr_dataset()
        if ds is None:
            messagebox.showinfo("Export Instrument Curves", "No instrument curves loaded.")
            return
        path = self._choose_save_file("Export Instrument Curves", '.csv', [("CSV", "*.csv")])
        if not path:
            return
        df = self.instrument_table(ds)
        if self._write_csv(df, path, "Export Instrument Curves"):
            messagebox.showinfo("Export Instrument Curves", "%d points, %d curves saved:\n%s"
                                % (len(df), df.shape[1], path))

    def on_exit(self):
        if self._dirty and not messagebox.askyesno(
                "Exit", "There are derived traces that were not exported. Quit?"):
            return
        self.root.destroy()


def main():
    root = TkinterDnD.Tk() if HAS_DND else tk.Tk()
    app = HPLCManager(root)
    for a in sys.argv[1:]:
        app.process_file(a)
    if len(sys.argv) > 1:
        app.refresh_view()
    root.mainloop()


if __name__ == '__main__':
    main()
