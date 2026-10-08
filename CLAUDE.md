# CLAUDE.md

Questo file guida Claude Code (claude.ai/code) quando lavora su questo repository. Tutti i `CLAUDE.md` sono in
italiano. Le regole comuni a tutti i progetti di Francesco (stack, convenzioni di codice, regole del repository,
git/GitHub, mirror su Drive, come lavora) stanno in `..\CLAUDE.md`; qui ci sono i dettagli del progetto.

## Scopo e priorita'

Visualizzatore di dati HPLC-DAD Agilent ChemStation (Rev. A.10.02, serie 1100, DAD G1315B). Caso d'uso centrale:
**aprire una cartella `.D` grezza, scegliere lunghezza d'onda, banda e riferimento (acceso/spento), vedere la traccia
estratta, cliccare sul cromatogramma per vedere lo spettro a quel tempo**. ChemStation esporta solo i segnali scelti
prima dell'acquisizione, anche se gli spettri DAD completi sono tutti nella cartella; OpenChrom (provato prima) ha
fallito con Francesco proprio su questo. Quindi: le operazioni base devono riuscire al primo tentativo, la traccia
estratta e' la vista predefinita (i canali `.ch` registrati sono opzionali) e ogni azione modifica cio' che e'
selezionato nella lista.

## Avvio e test

```bash
HPLCManager.bat [file_o_cartella_.D]
pythonw HPLCManager.pyw [file_o_cartella_.D]
py -m unittest discover -s tests -v     # i test sul campione reale si saltano se manca la cartella
```

Dati di esempio: `example data/009-0201.D` (gitignorata: dati reali, mai da committare). I test Tk richiedono uno
schermo. Gli screenshot del README si fanno dal programma vero (`ImageGrab` sulla finestra, finestra in primo piano)
e si controllano prima di pubblicare.

## Architettura

File unico `HPLCManager.pyw`, una classe `HPLCManager` (Tkinter + matplotlib), `tk.PanedWindow` a 3 pannelli come
LabSpectrumManager (a ogni pannello si danno `width` e `stretch` espliciti: tabella a sinistra 400 px fissi, grafico
al centro elastico, pannello destro 400 px fissi): tabella dati (sinistra), grafico con barra strumenti e lettura
del cursore (centro), pannello destro scrollabile (`self.f_right_inner`: lista cromatogrammi, metadati, bottoni
degli strumenti, pannello dello strumento in `self.f_tool_host`).

**Modello dati.** `self.chromatograms = {name: {'df': DataFrame (indice = tempo in min, colonna 'mAU'), 'info': {...},
'peaks': [...], 'cs_peaks': [...], 'hidden': bool, 'estr': {...}}}`. `peaks` = picchi calcolati qui, `cs_peaks` =
picchi del `REPORTnn.CSV` di ChemStation nella stessa cartella, `estr` (solo per le tracce estratte dagli spettri) =
`{ds, l, b, ref, rb}`, i parametri con cui e' stata calcolata. `self.spectra` = `{nome cartella: {'t', 'wl', 'S',
'info'}}` (cubo DAD), `self.instrument` = `{nome cartella: {'signals', 'modules', 'column'}}`.

**Caricamento** (`process_file`): una cartella `.D` carica il `.uv` (`_load_uv`), i registri dello strumento
(`_load_instrument`) ed esegue `_initial_extraction()` (lunghezza d'onda del primo segnale registrato, 4 nm,
riferimento spento); i suoi `.ch` si caricano solo se File -> "Also load recorded channels" e' attivo (o se manca il
`.uv`). Si carica anche un singolo `.ch` o un CSV a due colonne. `_reset_view` fa riscalare automaticamente il
prossimo ridisegno (dati nuovi); altrimenti `_redraw()` conserva lo zoom dell'utente (ripristina i limiti degli assi
se non sono in scala automatica).

**Estrazione** (il cuore): `extract_wavelength(wl, S, center, bandwidth, ref, ref_bandwidth)` = media pesata per
sovrapposizione della banda meno la stessa per il riferimento, con l'empirico `UV_OFFSET_NM`. Pannello del segnale del
rivelatore (`open_spectra`): `_sp_extract(False)` (*Update trace*, anche Invio in un campo) sostituisce **la traccia
estratta selezionata** al suo posto (`_sp_target()`: la traccia selezionata se ha `estr` e appartiene al dataset;
altrimenti l'ultima toccata, `_sp_live[ds]`); `_sp_extract(True)` (*Add as new*) ne aggiunge una e la seleziona.
Selezionare una traccia nella lista (`_on_selection()`) carica la sua lunghezza d'onda/banda/riferimento nel
pannello. `open_multi_extraction()` e' il dialogo "Several wavelengths" (N righe, ognuna con lunghezza d'onda, banda e
riferimento facoltativo) su `extract_series()` (tutto o niente).

**Spettri.** Click sul cromatogramma = spettro a quel tempo: si decide al *rilascio* del tasto (mosso < 5 px = click,
quindi funziona anche con zoom/pan della barra attivi; Maiusc+click funziona anche con la casella spenta). Gli
spettri scelti stanno in `self._chosen_spectra` (`dataset, t, wl, y, color, visible, line`, piu' `row`/`var` di Tk);
la finestra degli spettri (`_sp_window`) li elenca con una casella (`_sp_toggle`) e una crocetta di rimozione
(`_sp_remove`). Le linee tratteggiate sul cromatogramma (`gid='marker'`) vengono da `_visible_markers()` /
`_update_markers()`, senza ridisegno completo (zoom conservato), e si tolgono dalle figure salvate (`_copy_figure`).
Mappa tempo-lunghezza d'onda: `_sp_map`.

**Grafico.** La legenda delle tracce e' fissa in `upper right` (con `best` matplotlib la ricalcola a ogni ridisegno e
salta quando la linea del cursore ci passa sopra); la lettura del cursore sta subito sotto, ancorata allo stesso bordo.

**Picchi** (`open_peaks`, `_pk_apply`): `find_peaks` di scipy, basi dei picchi = minimi locali piu' vicini su una
copia Savitzky-Golay (le `left_bases` di scipy possono arrivare all'inizio del cromatogramma), baseline lineare, area
in mAU*s. Vicino a ChemStation per i picchi isolati (picco principale 101%), 50-90% per quelli sovrapposti o piccoli:
il suo integratore non e' replicato. Ci sono anche Smoothing, Trim, Normalize (ogni traccia derivata e' una nuova
voce e imposta `self._dirty`).

**Dati dello strumento.** `read_diagnostics()` (LCDIAG.REG) e `read_acqres()` (ACQRES.REG) leggono i registri `.REG`;
`read_folder_metadata()` aggiunge ai metadati Report00.CSV, RUN.LOG e colonna/moduli (mostrati per primi: Sample,
Date, Method, Column...). `open_instrument_curves()` e' la finestra delle curve impilate (raggruppate per unita');
`instrument_table()` interpola i segnali piu' lenti sulla base dei tempi della pompa per l'esportazione CSV.

**Esportazione / sessione** (menu File): `peaks_table()` / `spectra_table()` (solo gli spettri spuntati) costruiscono
i CSV; `export_full_dad`; `export_instrument_curves`; `save_session` / `open_session` salvano con pickle
`_session_data()` in un `.hplcsession` gzip; `_copy_figure()` crea la copia in stile Origin usata da
`save_figure_image`, `save_figure_pickle` e `open_figure_editor` (`plot_editor.pyw` di PlotStyleKit, caricato per
percorso). I CSV usano `;` e latin-1.

## Formati dei file ChemStation (tutto verificato su una sola corsa, `009-0201.D`)

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

## Punti aperti / idee

- Verificare tutto su piu' cartelle `.D` (altri metodi, colonne, intervalli spettrali, versioni di ChemStation): finora
  una sola corsa.
- Confrontare area/altezza calcolate qui con `REPORTnn.CSV` (tolleranza da decidere); replicare l'integratore di
  ChemStation o tenere il nostro? Calibrazione con standard esterni, idoneita' del sistema, purezza dei picchi
  (spettri sul picco), libreria spettrale, sovrapposizione del gradiente (% B) sul cromatogramma su un asse
  secondario, estrazione in batch su una serie di cartelle `.D`: tutto chiesto o considerato, niente fatto. Chiedere a
  Francesco quali sono le sue analisi (qui: polifenoli, gradiente di 90 min) prima di costruirle.

## Convenzioni di modifica (progetto)
- Seguire `..\CLAUDE.md`: codice, commenti, docstring, testi dell'interfaccia e dei `messagebox`, test, commit e
  README in inglese; i `CLAUDE.md` in italiano; modifiche chirurgiche. Il codice di questo progetto e' stato tradotto
  dall'italiano all'inglese in un unico passaggio (identificatori, commenti, messaggi, test, cartella `example data/`).
- I nomi esterni del repo fratello `PlotStyleKit` (`applica_rcparams`, `applica_stile_origin`, `carica_figura`) sono
  ancora in italiano li' e si chiamano cosi' come sono.
- Le sessioni salvate prima della traduzione (`.hplcsession` con chiavi italiane come `cromatogrammi`) non si
  possono aprire.
- Provare le funzioni con eventi veri (selezione nella lista con `<<ListboxSelect>>`, eventi del mouse sul canvas) e
  guardare il risultato in uno screenshot; i test Tk devono ripulire le proprie finestre.
