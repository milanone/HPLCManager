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

**Extraction** (the heart): `extract_wavelength(wl, S, center, bandwidth, ref, ref_bandwidth)` = overlap-weighted
mean of the band minus the same for the reference, with the empirical `UV_OFFSET_NM`. Detector signal panel
(`open_spectra`): `_sp_extract(False)` (*Update trace*, also Enter in a field) replaces **the selected extracted
trace** in place (`_sp_target()`: the selected trace if it has `estr` and belongs to the dataset; else the last
touched one, `_sp_live[ds]`); `_sp_extract(True)` (*Add as new*) adds one and selects it. Selecting a trace in the
list (`_on_selection()`) loads its wavelength/bandwidth/reference into the panel. `open_multi_extraction()` is the
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

## ChemStation file formats (all verified on one run, `009-0201.D`)

Software: **Agilent ChemStation Rev. A.10.02 [1757]** (`.ch` version 30, `.uv` version 31), "Instrument 1" with a
**DAD (G1315B)**, 1100 series modules, solvents H2O+H3PO4 / ACN. Sample: 30/09/2026, method `POLIFENB.M`, 20 uL,
90 min, 5 DAD signals A-E. One `.D` folder = one injection (on the acquisition PC
`C:\HPCHEM\1\DATA\<day>\<name>.D`). Instrument serial numbers must not be written in the public repo.

### Contents of the `.D` folder
| File | Contents |
|---|---|
| `dad1A.ch` ... `dad1E.ch` | binary chromatographic signals, one per DAD channel (A = 280,4 ref 360,100; B = 295,4 ref 360,100; ...) |
| `dad1.uv` | full DAD spectra (time x wavelength), ~6 MB, version 31 |
| `Report.TXT` | readable integration report (Area Percent), **latin-1**, CRLF |
| `Report00.CSV` | metadata: sample, method, pressure, flow, solvents; format `"field","value","unit"` |
| `REPORT01..05.CSV` | peak tables, **one per signal**, no header |
| `RUN.LOG` | instrument log, latin-1 |
| `ACQRES.REG`, `LCDIAG.REG` | binary registers (modules and column; instrument curves) |
| `SAMPLE.MAC` | macro |

`REPORTnn.CSV` columns: `#, RT [min], Type, Width [min], Area [mAU*s], Height [mAU], Area %`
(e.g. `1,2.7838,"VV  ",0.0854,100.785,16.941,0.1824`). The type is the ChemStation code (BB, BV, VV, VB, PB, BP,
VP, PV...). Values are at full float precision: they are the **reference truth** for validating the integrator.
`REPORT04.CSV` is empty (signal without peaks); `REPORT05` has a single peak (signal at 520 nm).

### `.ch` format, version 30
Everything **big-endian**.
- `0x000`: Pascal strings (length + text) at fixed offsets: version (`"30"`), file type (`LC DATA FILE`, 0x04),
  sample name (0x18), date (0xB2, `30-Sep-26, 11:08:12`), module (0xD0, `G1315B`), method (0xE4, `POLIFENB.M`);
  the signal description (`DAD A, Sig=280,4 Ref=360,100`) is before 0x400.
- `0x11A` int32 = start time (ms, here `-2560`); `0x11E` int32 = end time (ms, here `5397440`). Step 400 ms
  (2.5 Hz): `13501` points = `(end - start) / 400 + 1`. Time of point i: `(start + i * 400) / 60000` min.
- Data from `0x400` to end of file, in **segments**: `0x10 n` (marker + number of entries, always 25 here)
  followed by `n` entries, each either a signed `int16` = **delta** from the previous value, or `0x8000` + `int32`
  = **absolute** value (resets the accumulation; 6 bytes in all). The accumulation starts at 0 and continues
  across segments; the last segment may be shorter.
- Decoded value / **2000** = mAU (maximum at 4.357 min: 4201053 -> 2100.5 mAU; the report says height
  2003.05 mAU, i.e. net of the integration baseline, about -97 mAU there).
- Comparison with the REPORTs (`tests/test_chemstation_reader.py`): for peaks >= 50 mAU within 30 min the apex falls
  less than 0.01 min from the report RT and `signal / height` is between 1.00 and 1.15. For small peaks the
  ChemStation baseline lies **below** the local minimum of the signal, and beyond 80 min it is negative (-10/-40
  mAU): comparing heights and areas would require replicating the integrator. Many peaks are shoulders (apex not
  coinciding with the RT).
- Open: the 1/2000 factor is inferred, not read from a header field. Other `.ch` versions (130/179/181) have a
  different layout: confirm on other samples.

### `dad1.uv` format, version 31
The full DAD cube: one spectrum per instant, here **13493 spectra x 201 wavelengths (200-600 nm, step 2 nm)**,
every 400 ms. The `dad1A..E.ch` are just profiles extracted during the acquisition (`Sig=280,4 Ref=360,100` = mean
278-282 nm minus mean 310-410 nm). **Little-endian** (unlike the `.ch`).
- `0x000`: Pascal string with the version (`"31"`) and the same fixed-offset strings as the `.ch`; the records start
  at `0x200`.
- Record (one per spectrum, variable length 424-452 bytes): `uint16` type (`0x43`), `uint16` **total** record
  length, `uint32` time (ms; 240, 640, 1040...), `uint16` start wavelength, end wavelength and step in 1/20 nm
  (4000, 12000, 40), 8 undecoded bytes (they contain `0x50/0x51`, `0x04`, `0x190`=400 ms), then the values.
- Values: as in the `.ch`, little-endian: `int16` = delta from the previous value (starts at 0 in every record),
  `0x8000` + `int32` = absolute value. All 13493 records give exactly 201 values and consume exactly the declared
  bytes. Divided by **2000** = mAU. `read_uv()` takes ~1 s (numpy fast path for records without absolute values).
- After the last record there is a ~135 kB tail (maybe an index, ~10 bytes per spectrum): not read.
- **Time**: `.uv` and `.ch` use the same axis in ms; the `.ch` starts at -2560 ms, the `.uv` at 240 ms, so
  `ch_index = uv_index + 7`.
- **Band extraction**: weighted mean, weight = overlap between the band `center +/- bandwidth/2` and the 2 nm
  interval of each spectral point, minus the same for the reference. RMS deviation from `dad1A..D.ch`: 0.17 /
  0.03 / 0.01 / 0.02 mAU (A = 280 nm, main peak at 2100 mAU; maximum deviation ~7 mAU = 0.3%).
- The wavelength axis must be shifted by **+0.45 nm** (`UV_OFFSET_NM`): an **empirical** value (minimum of the
  deviation between 0 and 1 nm), not read from the header. Without it the RMS deviation grows about 6 times (A 1.1;
  B 1.8 mAU). Origin unclear (pixel centre? lamp calibration?): to be rechecked with other samples.

### Metadata: what and where
- binary header of `.ch`/`.uv` (see above);
- `Report00.CSV`: sample, data file, instrument, method and changes, injection date, sequence line, vial, injected
  volume, sequence file, pressure and flow (start/end), solvents A-D, signal list. The `Acq. Operator` line holds the
  injection date, not the operator (it is discarded);
- `RUN.LOG`: column temperature (min-max over the thermostat messages).
- **Not read**: the rest of `RUN.LOG` (pressures over time, events), `Report.TXT` (repeats the REPORTnn.CSV), the 8
  undecoded bytes of each `.uv` record, the `.uv` tail, the `.M` method file (it is not in the `.D` folder).

### `.REG` registers
Same container for `ACQRES.REG` and `LCDIAG.REG`: Pascal string `"32"`, `REGISTER FILE`, version `A.00.01`, then
tables with headers (`ObjClass`, `Title`, `DataType`...) and data. Little-endian.

**`LCDIAG.REG`** = the profiles the instrument recorded during the run (7 signals in the sample):
- each signal is an **array of `uint32`** preceded by the strings `min\0<unit>\0` (`bar`, `ml/min`, `%`, `°C`); the
  physical value is `uint32 / factor` with factor **100** (bar), **1000** (ml/min), **10** (%), **100** (°C);
- number of points = `uint32` **168 bytes before** `min`; sampling interval = `double` (in **minutes**) **116 bytes
  before** (0.005 min = 0.3 s for the pump, 1/60 min = 1 s for the thermostat);
- the **title** (`PMP1, Pressure`...) is in the header that **follows** the array, before `arial`;
- signals in the sample: pressure (18000 points), flow, solvents A-D (composition %), left column temperature
  (5400 points). Time starts at 0 (assumption: consistent with the pressure at 4 s, 96.05 bar, against 96.45 bar in
  `Report00.CSV`, and with the final pressure, 104.5 bar against 104.72). A+B+C+D = 100 ±0.5%.
- Other signals from different modules (autosampler, DAD: lamp temperatures, burn hours...) are in the initial part
  (tables `Start/Stop Conditions`, `VisBurnTime`, `UVOnTime`...): **not read**.

**`ACQRES.REG`** = modules and column: strings as `uint16` length (with NUL) + text + NUL.
- **Modules** (4 in the sample): five columns of consecutive strings, in this order: serialNumber, FWrevision,
  buildNumber, Name, PartNumber, one row per module (G1311A quaternary pump FW A.05.06; G1313A autosampler
  FW A.05.06; G1315B DAD FW A.05.09; G1316A thermostat FW A.05.09).
- **Column**: description `ODS Hypersil` (the string that is not a software version, instrument or path); length
  100 mm, diameter 2.1 mm, particle size 5 um are the first three "round" doubles of the numeric data. The
  assignment to the fields (ColLength, ColDiameter, ParticleSize) is **inferred** (field order and typical values
  of a Hypersil ODS 100 x 2.1 mm, 5 um), not verified on other files; a fourth double (68) is not identified and is
  not shown.

## Open points / ideas

- Verify everything on more `.D` folders (other methods, columns, spectral ranges, ChemStation versions): so far
  one run only.
- Compare area/height computed here with `REPORTnn.CSV` (tolerance to decide); replicate ChemStation's integrator or
  keep our own? Calibration with external standards, system suitability, peak purity (spectra on the peak),
  spectral library, gradient (% B) overlay on the chromatogram on a secondary axis, batch extraction over a series of
  `.D` folders: all asked about or considered, none done. Ask Francesco what his analyses are (here: polyphenols,
  90 min gradient) before building them.

## Editing conventions (project)
- Follow `..\CLAUDE.md`: everything in English (identifiers, comments, docstrings, UI and `messagebox` text, tests,
  commits, README and this file), surgical edits. The code was translated from Italian to English in one pass
  (identifiers, comments, messages, tests, the `example data/` folder).
- External names from the sibling `PlotStyleKit` repo (`applica_rcparams`, `applica_stile_origin`,
  `carica_figura`) are still Italian there and are called as they are.
- Sessions saved before the translation (`.hplcsession` with Italian keys such as `cromatogrammi`) cannot be opened.
- Test GUI features with real events (list selection via `<<ListboxSelect>>`, mouse events on the canvas), and look
  at the result in a screenshot; the Tk tests must clean up their windows.
