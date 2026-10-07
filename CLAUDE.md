# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Knowledge base

General working style lives in the root folder: `..\knowledge\modo-di-lavorare.md` (also
loaded via `..\CLAUDE.md`). Project-specific notes live here:
- `knowledge/chemstation.md` — Agilent ChemStation data formats (`.ch` v30 and `.uv` v31 decoded and verified on `009-0201.D`).

When the user states a preference or a decision, add it to the relevant `knowledge/` file.
The closest reference project is `..\LabSpectrumManager` (same layout, same conventions).

## Running the Application

```bash
HPLCManager.bat [optional_chromatogram_file]
pythonw HPLCManager.pyw [optional_chromatogram_file]
```

## Architecture

Single file `HPLCManager.pyw`, one class `HPLCManager` (Tkinter + matplotlib), 3-pane layout like
LabSpectrumManager: data table (left), plot with toolbar and cursor readout (centre), scrollable
right panel (`self.f_right_inner`: chromatogram list, metadata, tool panels in `self.f_tool_host`).

Data model: `self.cromatogrammi = {name: {'df': DataFrame (index = time in min, column 'mAU'),
'info': {...}, 'picchi': [...], 'picchi_cs': [...], 'nascosto': bool}}`. `picchi` = peaks computed
here, `picchi_cs` = peaks from the ChemStation `REPORTnn.CSV` of the same `.D` folder.

Readers: `leggi_ch()` (`.ch` v30), `leggi_uv()` (DAD spectra `.uv` v31, little-endian), `leggi_report_csv()`, `processa_file()` (folder `.D` / `.ch` / CSV). Opening a `.D` folder also loads its `.uv` into `self.spettri` and (unless File -> "Also load recorded channels" is on) skips its `.ch` files, doing `_estrazione_iniziale()` instead: one 'live' trace per dataset (`self._sp_live`), replaced in place by `_sp_estrai(False)` (Update) or kept apart by `_sp_estrai(True)` (Add as new); reference band off by default (`{folder name: {'t','wl','S','info'}}`); `estrai_lambda()` builds a chromatogram at any wavelength/bandwidth (+ optional reference band) like ChemStation's `Sig=280,4 Ref=360,100`, with the empirical `UV_OFFSET_NM`.
Tools: Spectra (`apri_spettri`: extract chromatogram, spectrum at time / click on plot, time-wavelength map), Peaks (`apri_picchi`, `_pk_applica`: scipy `find_peaks`, bases = nearest local minima on a
Savitzky-Golay copy, linear baseline, area in mAU*s), Smoothing, Trim, Normalize; each derived trace
is a new entry and sets `self._dirty`. Peak areas are close to ChemStation's for isolated peaks
(main peak 101%) but 50-90% for overlapping/small ones: ChemStation's integrator is not replicated.

Export/session (File menu): `tabella_picchi()` / `tabella_spettri()` build the DataFrames written by
`esporta_picchi` / `esporta_spettri`; `esporta_dad_completo` writes the whole DAD matrix; `salva_sessione` /
`apri_sessione` pickle `_dati_sessione()` into a gzip `.hplcsession`; `_copia_figura()` makes the Origin-styled
copy (without cursor and the `gid='marcatore'` pick lines) used by `salva_figura_immagine`,
`salva_figura_pickle` and `apri_editor_figura` (PlotStyleKit `plot_editor.pyw`, loaded by path).
Exports are CSV with `;` separator and latin-1 encoding, like LabSpectrumManager.

`apri_estrazione_multipla()` is the 'Several wavelengths' dialog (N rows: wavelength, bandwidth, optional reference per row) over `estrai_serie()` (all-or-nothing, adds the traces like 'Add as new'; `_sp_multi_ultimo` restores the last values).

Picked spectra live in `self._spettri_scelti` (`dataset, t, wl, y, colore, visibile, linea`, plus the Tk `riga`/`var`); the spectra window (`_sp_finestra`) lists them with a checkbox (hide/show: `_sp_spunta`) and a remove cross (`_sp_rimuovi`); the dashed pick lines on the chromatogram come from `_marcatori_visibili()`/`_aggiorna_marcatori()` (no full redraw, zoom kept), and `tabella_spettri()` exports only the ticked ones.

Tests: `py -m unittest discover -s tests -v` (the real-sample tests are skipped if
`campioni di esempio/009-0201.D` is missing).

## Editing conventions
- Edits must be surgical and non-destructive; never refactor or rename existing methods unless asked.
- Preserve existing comments and docstrings; when in doubt, ask.
- Code comments and identifiers are Italian; button/label/panel text is English; `messagebox`
  body text is Italian.
