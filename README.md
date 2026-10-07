# HPLCManager

Viewer and analysis tool for HPLC chromatograms exported by Agilent ChemStation: load, overlay and
compare chromatograms, peak detection and integration, quantification.

Open a ChemStation `.D` folder (File menu, button or drag-and-drop) and all its signals load together;
the red triangles mark the peaks ChemStation reported. Tools: Peaks (detection + integration, table),
Smoothing (Savitzky-Golay), Trim, Normalize; data table, cursor readout, CSV export.
Opening a `.D` folder loads its DAD spectra (`.uv`) and shows one extracted trace; the recorded
channels (`.ch`) are loaded only if File -> "Also load recorded channels" is ticked. The *Detector
signal* panel sets wavelength, bandwidth and an optional reference band (checkbox: off = no reference);
Enter or *Update trace* recomputes the trace, *Add as new* keeps a copy for comparison. *Several wavelengths...* opens a dialog where you choose how many traces you want and, for each,
wavelength, bandwidth and an optional reference (each row has its own); all are extracted in one go
(all or nothing if a row is invalid). The same panel
shows the spectrum at a given time or at a click on the plot (a click, also while the toolbar's zoom/pan is active (dragging still zooms/pans); Shift+click works even with the checkbox off; a dashed line marks each picked time; in the spectra window a list on the right lets you untick a spectrum to hide it - its dashed line and the CSV export follow - or remove it with the cross), and draws a time-wavelength map.

**Instrument data.** The metadata panel also shows the column (description, size, particle size) and the
instrument modules (part number, serial number, firmware), read from `ACQRES.REG`. *Instrument curves*
(button in the right panel, or Tools menu) plots what the instrument recorded during the run, from
`LCDIAG.REG`: pressure, solvent composition (gradient), flow and column temperature, on a common time axis;
*Export Instrument Curves* saves them as CSV (slower signals interpolated onto the pump time base).

**Saving and exporting** (File menu): *Save/Open Session* (`.hplcsession`: all traces, derived ones
included, peaks and DAD spectra); *Export Traces* (merged CSV); *Export Peaks* (HPLCManager and
ChemStation peaks, one per row, with a Source column); *Export Spectra* (the spectra picked on the plot);
*Export Full DAD Data* (time x wavelength matrix); *Save Figure Image* (PNG/PDF/SVG, 300 dpi);
*Save Figure (pickle)* and *Edit Figure...* with the shared PlotStyleKit editor. Figures and images
leave out the cursor and the dashed lines of the picked spectra.

Supported: ChemStation `.ch` version 30 and DAD `.uv` version 31 (decoded natively, see `knowledge/chemstation.md`) and
two-column CSV. Peak areas are close to ChemStation's for isolated peaks, lower for overlapping ones.

## Running

```
pythonw HPLCManager.pyw [chromatogram_file]
```

## Dependencies

`numpy`, `pandas`, `matplotlib`, `tkinter` (standard library); `scipy` and `tkinterdnd2` optional.
[PlotStyleKit](https://github.com/milanone/PlotStyleKit) is an optional sibling repo (`../PlotStyleKit`).
