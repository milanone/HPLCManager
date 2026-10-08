# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.
General rules for all of Francesco's projects (stack, code conventions, repository rules, git/GitHub, the
Drive mirror, how he works) are in `..\CLAUDE.md`; this file has the project-specific details.

## Purpose and priorities

Viewer for Agilent ChemStation HPLC-DAD data (Rev. A.10.02, 1100 series, DAD G1315B). The central use case:
**open a raw `.D` folder, choose wavelength, bandwidth and reference (on/off), see the extracted trace, click
the chromatogram to see the spectrum at that time**. ChemStation only exports the signals chosen before the
acquisition even though the DAD spectra are all in the folder; OpenChrom (tried before) failed Francesco on
exactly this. So: the basic operations must work first time, the extracted trace is the default view (the
recorded `.ch` channels are optional), and every action edits what is selected in the list.

## Running / tests

```bash
HPLCManager.bat [file_or_.D_folder]
pythonw HPLCManager.pyw [file_or_.D_folder]
py -m unittest discover -s tests -v     # real-sample tests are skipped if the sample folder is missing
```

Sample data: `example data/009-0201.D` (gitignored: real data, never commit). Tk tests need a screen.
Screenshots for the README are taken from the real program (`ImageGrab` on the window, window topmost) and
checked before publishing.

## Architecture

Single file `HPLCManager.pyw`, one class `HPLCManager` (Tkinter + matplotlib), 3-pane `tk.PanedWindow` like
LabSpectrumManager (each pane given an explicit `width` and `stretch`: left table 400 px fixed, centre plot
elastic, right panel 400 px fixed): data table (left), plot with toolbar and cursor readout (centre), scrollable
right panel (`self.f_right_inner`: chromatogram list, metadata, tool buttons, tool panel in `self.f_tool_host`).

**Data model.** `self.chromatograms = {name: {'df': DataFrame (index = time in min, column 'mAU'), 'info': {...},
'peaks': [...], 'cs_peaks': [...], 'hidden': bool, 'estr': {...}}}`. `peaks` = peaks computed here,
`cs_peaks` = peaks from the ChemStation `REPORTnn.CSV` of the same folder, `estr` (only for traces extracted from
the spectra) = `{ds, l, b, ref, rb}`, the parameters it was computed with. `self.spectra` = `{folder name: {'t',
'wl', 'S', 'info'}}` (DAD cube), `self.instrument` = `{folder name: {'signals', 'modules', 'column'}}`.

**Loading** (`process_file`): a `.D` folder loads `.uv` (`_load_uv`), the instrument registers
(`_load_instrument`), and runs `_initial_extraction()` (first recorded signal's wavelength, 4 nm, reference off);
its `.ch` files load only if File -> "Also load recorded channels" is on (or if there is no `.uv`). A single `.ch`
or a two-column CSV also loads. `_reset_view` makes the next redraw autoscale (new data); otherwise
`_redraw()` keeps the user's zoom (it restores the axis limits if the axes are not on autoscale).

**Extraction** (the heart): `extract_wavelength(wl, S, center, bandwidth, ref, ref_bandwidth)` = overlap-weighted mean of the
band minus the same for the reference, with the empirical `UV_OFFSET_NM`. Detector signal panel (`open_spectra`):
`_sp_extract(False)` (*Update trace*, also Enter in a field) replaces **the selected extracted trace** in place
(`_sp_target()`: the selected trace if it has `estr` and belongs to the dataset; else the last touched one,
`_sp_live[ds]`); `_sp_extract(True)` (*Add as new*) adds one and selects it. Selecting a trace in the list
(`_on_selection()`) loads its wavelength/bandwidth/reference into the panel. `open_multi_extraction()` is the
"Several wavelengths" dialog (N rows, each with wavelength, bandwidth and an optional reference) over
`extract_series()` (all or nothing).

**Spectra.** Click on the chromatogram = spectrum at that time: decided on button *release* (moved < 5 px =
click, so it also works with the toolbar's zoom/pan active; Shift+click works even with the checkbox off). Picked
spectra live in `self._chosen_spectra` (`dataset, t, wl, y, color, visible, line`, plus the Tk `row`/`var`);
the spectra window (`_sp_window`) lists them with a checkbox (`_sp_toggle`) and a remove cross
(`_sp_remove`). The dashed pick lines on the chromatogram (`gid='marker'`) come from
`_visible_markers()` / `_update_markers()`, without a full redraw (zoom kept) and are removed from saved
figures (`_copy_figure`). Time-wavelength map: `_sp_map`.

**Plot.** The trace legend is fixed at `upper right` (with `best` matplotlib recomputes it on every redraw and it
jumps when the cursor line passes over it); the cursor readout sits right under it, anchored to the same edge.

**Peaks** (`open_peaks`, `_pk_apply`): scipy `find_peaks`, peak bases = nearest local minima on a
Savitzky-Golay copy (scipy's `left_bases` can reach the start of the chromatogram), linear baseline, area in
mAU*s. Close to ChemStation for isolated peaks (main peak 101%), 50-90% for overlapping/small ones: its
integrator is not replicated. Also Smoothing, Trim, Normalize (each derived trace is a new entry and sets
`self._dirty`).

**Instrument data.** `read_diagnostics()` (LCDIAG.REG) and `read_acqres()` (ACQRES.REG) read the `.REG`
registers; `read_folder_metadata()` adds Report00.CSV, RUN.LOG and the column/modules to the metadata (shown
first: Sample, Date, Method, Column...). `open_instrument_curves()` is the stacked-curves window (grouped by unit);
`instrument_table()` interpolates the slower signals onto the pump time base for the CSV export.

**Export / session** (File menu): `peaks_table()` / `spectra_table()` (only the ticked spectra) build the
CSVs; `export_full_dad`; `export_instrument_curves`; `save_session` / `open_session` pickle
`_session_data()` into a gzip `.hplcsession`; `_copy_figure()` makes the Origin-styled copy used by
`save_figure_image`, `save_figure_pickle` and `open_figure_editor` (PlotStyleKit `plot_editor.pyw`, loaded by
path). CSVs use `;` and latin-1.

## ChemStation file formats (notes in Italian; all verified on one run, `009-0201.D`)

Software: **Agilent ChemStation Rev. A.10.02 [1757]** (file `.ch` versione 30, `.uv` versione 31), "Instrument 1"
con **DAD (G1315B)**, moduli serie 1100, solventi H2O+H3PO4 / ACN. Campione: 30/09/2026, metodo `POLIFENB.M`,
20 uL, 90 min, 5 segnali DAD A-E. Una cartella `.D` = una iniezione (sul PC di acquisizione
`C:\HPCHEM\1\DATA\<giorno>\<nome>.D`). I numeri di serie degli strumenti non vanno scritti nel repo pubblico.

### Contenuto della cartella `.D`
| File | Contenuto |
|---|---|
| `dad1A.ch` ... `dad1E.ch` | segnali cromatografici binari, uno per canale DAD (A = 280,4 ref 360,100; B = 295,4 ref 360,100; ...) |
| `dad1.uv` | spettri DAD completi (tempo x lunghezza d'onda), ~6 MB, versione 31 |
| `Report.TXT` | report di integrazione leggibile (Area Percent), **latin-1**, CRLF |
| `Report00.CSV` | metadati: campione, metodo, pressione, flusso, solventi; formato `"campo","valore","unita'"` |
| `REPORT01..05.CSV` | tabelle picchi, **una per segnale**, senza intestazione |
| `RUN.LOG` | log dello strumento, latin-1 |
| `ACQRES.REG`, `LCDIAG.REG` | registri binari (moduli e colonna; curve dello strumento) |
| `SAMPLE.MAC` | macro |

`REPORTnn.CSV`, colonne: `#, RT [min], Tipo, Larghezza [min], Area [mAU*s], Altezza [mAU], Area %`
(es. `1,2.7838,"VV  ",0.0854,100.785,16.941,0.1824`). Il tipo e' il codice ChemStation (BB, BV, VV, VB, PB, BP,
VP, PV...). Valori a precisione float piena: sono la **verita' di riferimento** per validare l'integratore.
`REPORT04.CSV` e' vuoto (segnale senza picchi); `REPORT05` ha un solo picco (segnale a 520 nm).

### Formato `.ch` versione 30
Tutto **big-endian**.
- `0x000`: stringhe Pascal (lunghezza + testo) a offset fissi: versione (`"30"`), tipo file (`LC DATA FILE`, 0x04),
  nome campione (0x18), data (0xB2, `30-Sep-26, 11:08:12`), modulo (0xD0, `G1315B`), metodo (0xE4,
  `POLIFENB.M`); la descrizione del segnale (`DAD A, Sig=280,4 Ref=360,100`) e' prima di 0x400.
- `0x11A` int32 = tempo iniziale (ms, qui `-2560`); `0x11E` int32 = tempo finale (ms, qui `5397440`). Passo 400 ms
  (2,5 Hz): `13501` punti = `(fine - inizio) / 400 + 1`. Tempo del punto i: `(inizio + i * 400) / 60000` min.
- Dati da `0x400` a fine file, in **segmenti**: `0x10 n` (marcatore + numero di voci, qui sempre 25) seguito da `n`
  voci, ciascuna: `int16` con segno = **delta** dal valore precedente; oppure `0x8000` + `int32` = valore
  **assoluto** (resetta l'accumulo; 6 byte in tutto). L'accumulo parte da 0 e prosegue tra i segmenti; l'ultimo
  segmento puo' essere piu' corto.
- Valore decodificato / **2000** = mAU (massimo a 4,357 min: 4201053 -> 2100,5 mAU; il report dice altezza
  2003,05 mAU, cioe' al netto della baseline di integrazione, circa -97 mAU li').
- Confronto con i REPORT (`tests/test_chemstation_reader.py`): per i picchi >= 50 mAU entro 30 min l'apice cade a
  meno di 0,01 min dall'RT del report e `segnale / altezza` sta tra 1,00 e 1,15. Per i picchi piccoli la baseline di
  ChemStation sta **sotto** il minimo locale del segnale e oltre 80 min e' negativa (-10/-40 mAU): per confrontare
  altezze e aree servirebbe replicare l'integratore. Molti picchi sono spalle (apice non coincidente con l'RT).
- Aperto: il fattore 1/2000 e' dedotto, non letto da un campo dell'intestazione. Altre versioni `.ch`
  (130/179/181) hanno layout diverso: confermare su altri campioni.

### Formato `dad1.uv` versione 31
Il cubo DAD completo: uno spettro per ogni istante, qui **13493 spettri x 201 lunghezze d'onda (200-600 nm, passo
2 nm)**, ogni 400 ms. I `dad1A..E.ch` sono solo profili estratti durante l'acquisizione (`Sig=280,4 Ref=360,100` =
media 278-282 nm meno media 310-410 nm). **Little-endian** (al contrario dei `.ch`).
- `0x000`: stringa Pascal con la versione (`"31"`) e le stesse stringhe a offset fissi dei `.ch`; i record partono da
  `0x200`.
- Record (uno per spettro, lunghezza variabile 424-452 byte): `uint16` tipo (`0x43`), `uint16` lunghezza **totale**
  del record, `uint32` tempo (ms; 240, 640, 1040...), `uint16` lambda iniziale, finale e passo in 1/20 nm (4000,
  12000, 40), 8 byte non decodificati (contengono `0x50/0x51`, `0x04`, `0x190`=400 ms), poi i valori.
- Valori: come nei `.ch`, little-endian: `int16` = delta dal valore precedente (parte da 0 a ogni record),
  `0x8000` + `int32` = valore assoluto. Tutti i 13493 record danno esattamente 201 valori e consumano esattamente i
  byte dichiarati. Diviso **2000** = mAU. `read_uv()` ci mette ~1 s (percorso veloce numpy per i record senza
  valori assoluti).
- Dopo l'ultimo record c'e' una coda di ~135 kB (forse un indice, ~10 byte per spettro): non letta.
- **Tempo**: `.uv` e `.ch` usano lo stesso asse in ms; il `.ch` parte a -2560 ms, il `.uv` a 240 ms, quindi
  `indice_ch = indice_uv + 7`.
- **Estrazione di una banda**: media pesata, peso = sovrapposizione fra la banda `centro +/- banda/2` e l'intervallo
  di 2 nm di ciascun punto spettrale, meno lo stesso per il riferimento. Scarto RMS rispetto a `dad1A..D.ch`: 0,17 /
  0,03 / 0,01 / 0,02 mAU (A = 280 nm, picco principale a 2100 mAU; scarto massimo ~7 mAU = 0,3%).
- L'asse delle lunghezze d'onda va spostato di **+0,45 nm** (`UV_OFFSET_NM`): valore **empirico** (minimo dello
  scarto fra 0 e 1 nm), non letto dall'intestazione. Senza, lo scarto RMS sale di circa 6 volte (A 1,1; B 1,8 mAU).
  Origine non chiarita (centro del pixel? calibrazione della lampada?): da riverificare con altri campioni.

### Metadati: cosa e dove
- intestazione binaria di `.ch`/`.uv` (vedi sopra);
- `Report00.CSV`: campione, file dati, strumento, metodo e modifiche, data iniezione, riga sequenza, vial, volume
  iniettato, file di sequenza, pressione e flusso (inizio/fine), solventi A-D, elenco segnali. La riga
  `Acq. Operator` contiene la data dell'iniezione, non l'operatore (si scarta);
- `RUN.LOG`: temperatura colonna (min-max sui messaggi del termostato).
- **Non letti**: il resto di `RUN.LOG` (pressioni nel tempo, eventi), `Report.TXT` (ripete i REPORTnn.CSV), gli 8 byte
  non decodificati di ogni record `.uv`, la coda del `.uv`, il file del metodo `.M` (non e' nella cartella `.D`).

### Registri `.REG`
Stesso contenitore per `ACQRES.REG` e `LCDIAG.REG`: stringa Pascal `"32"`, `REGISTER FILE`, versione `A.00.01`, poi
tabelle con testate (`ObjClass`, `Title`, `DataType`...) e dati. Little-endian.

**`LCDIAG.REG`** = i profili che lo strumento ha registrato durante la corsa (nel campione 7 segnali):
- ogni segnale e' un **array di `uint32`** preceduto dalle stringhe `min\0<unita'>\0` (`bar`, `ml/min`, `%`, `°C`);
  il valore fisico e' `uint32 / fattore` con fattore **100** (bar), **1000** (ml/min), **10** (%), **100** (°C);
- numero di punti = `uint32` a **168 byte prima** di `min`; intervallo di campionamento = `double` (in **minuti**) a
  **116 byte prima** (0,005 min = 0,3 s per la pompa, 1/60 min = 1 s per il termostato);
- il **titolo** (`PMP1, Pressure`...) sta nella testata che **segue** l'array, prima di `arial`;
- segnali del campione: pressione (18000 punti), flusso, solventi A-D (composizione %), temperatura colonna sinistra
  (5400 punti). Il tempo parte da 0 (ipotesi: coerente con la pressione a 4 s, 96,05 bar, contro 96,45 bar in
  `Report00.CSV`, e con la pressione finale, 104,5 bar contro 104,72). A+B+C+D = 100 ±0,5%.
- Altri segnali di moduli diversi (autocampionatore, DAD: temperature lampada, ore di accensione...) sono nella parte
  iniziale (tabelle `Start/Stop Conditions`, `VisBurnTime`, `UVOnTime`...): **non letti**.

**`ACQRES.REG`** = moduli e colonna: stringhe come `uint16` lunghezza (con NUL) + testo + NUL.
- **Moduli** (nel campione 4): cinque colonne di stringhe consecutive, in questo ordine: serialNumber, FWrevision,
  buildNumber, Name, PartNumber, una riga per modulo (G1311A pompa quaternaria FW A.05.06; G1313A autocampionatore
  FW A.05.06; G1315B DAD FW A.05.09; G1316A termostato FW A.05.09).
- **Colonna**: descrizione `ODS Hypersil` (la stringa che non e' versione software, strumento o percorso); lunghezza
  100 mm, diametro 2,1 mm, granulometria 5 um sono i primi tre double "tondi" dei dati numerici. L'assegnazione ai
  campi (ColLength, ColDiameter, ParticleSize) e' **dedotta** (ordine dei campi e valori tipici di una Hypersil ODS
  100 x 2.1 mm, 5 um), non verificata su altri file; un quarto double (68) non e' identificato e non viene mostrato.

## Open points / ideas

- Verify everything on more `.D` folders (other methods, columns, spectral ranges, ChemStation versions): so far
  one run only.
- Compare area/height computed here with `REPORTnn.CSV` (tolerance to decide); replicate ChemStation's integrator or
  keep our own? Calibration with external standards, system suitability, peak purity (spectra on the peak),
  spectral library, gradient (% B) overlay on the chromatogram on a secondary axis, batch extraction over a series of
  `.D` folders: all asked about or considered, none done. Ask Francesco what his analyses are (here: polyphenols,
  90 min gradient) before building them.

## Editing conventions (project)
- Follow `..\CLAUDE.md`: everything in English (identifiers, comments, docstrings, UI and `messagebox` text,
  tests), surgical edits. The code was translated from Italian to English in one pass (identifiers, comments,
  messages, tests, the `example data/` folder). Only this file's format notes (sections below) stay in Italian.
  External names from the sibling `PlotStyleKit` repo (`applica_rcparams`, `applica_stile_origin`,
  `carica_figura`) are still Italian there and are called as they are.
- Sessions saved before the translation (`.hplcsession` with Italian keys such as `cromatogrammi`) cannot be opened.
- Test GUI features with real events (list selection via `<<ListboxSelect>>`, mouse events on the canvas), and look
  at the result in a screenshot; the Tk tests must clean up their windows.
