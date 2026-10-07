# Agilent ChemStation: formati dati

Software: **Agilent ChemStation Rev. A.10.02 [1757]** (file `.ch` versione 30, `.uv` versione 31),
"Instrument 1" con **DAD (G1315B)**, moduli serie 1100, solventi H2O+H3PO4 / ACN.
Campione verificato: `campioni di esempio/009-0201.D` (30/09/2026, metodo `POLIFENB.M`,
"Estr acq gambe", 20 uL, 90 min, 5 segnali DAD A-E). Sul disco di acquisizione:
`C:\HPCHEM\1\DATA\MD260929\009-0201.D` (una cartella `.D` = una iniezione).

## Contenuto della cartella `.D`
| File | Contenuto |
|---|---|
| `dad1A.ch` ... `dad1E.ch` | segnali cromatografici binari, uno per canale DAD (A = 280,4 ref 360,100; B = 295,4 ref 360,100; ...) |
| `dad1.uv` | spettri DAD completi (tempo x lunghezza d'onda), ~6 MB, versione 31 |
| `Report.TXT` | report di integrazione leggibile (Area Percent), **latin-1**, CRLF |
| `Report00.CSV` | metadati: campione, metodo, pressione, flusso, solventi, colonna; formato `"campo","valore","unita'"` |
| `REPORT01..05.CSV` | tabelle picchi, **una per segnale**, senza intestazione |
| `RUN.LOG` | log dello strumento, latin-1 |
| `ACQRES.REG`, `LCDIAG.REG` | registri binari (diagnostica, configurazione moduli) |
| `SAMPLE.MAC` | macro |

`REPORTnn.CSV`, colonne: `#, RT [min], Tipo, Larghezza [min], Area [mAU*s], Altezza [mAU], Area %`
(es. `1,2.7838,"VV  ",0.0854,100.785,16.941,0.1824`). Il tipo e' il codice ChemStation
(BB, BV, VV, VB, PB, BP, VP, PV...). Valori a precisione float piena: sono la **verita' di
riferimento** per validare l'integratore di HPLCManager. `REPORT04.CSV` e' vuoto (segnale senza
picchi); `REPORT05` ha un solo picco (segnale a 520 nm).

## Formato `.ch` versione 30 (decodificato e verificato)
Tutto **big-endian**.
- `0x000`: stringhe Pascal (lunghezza + testo): versione (`"30"`), tipo file (`LC DATA FILE`),
  nome campione, data (`30-Sep-26, 11:08:12`), modulo (`G1315B`), metodo (`POLIFENB.M`),
  software e versione, unita'.
- `0x11A` int32 = tempo iniziale (ms, qui `-2560`); `0x11E` int32 = tempo finale (ms, qui
  `5397440`). Passo 400 ms (2,5 Hz): `13501` punti = `(fine - inizio) / 400 + 1`.
- Dati da `0x400` a fine file, in **segmenti**: `0x10 n` (marcatore + numero di voci, qui sempre 25)
  seguito da `n` voci, ciascuna:
  - `int16` con segno = **delta** rispetto al valore precedente;
  - oppure `0x8000` + `int32` = valore **assoluto** (resetta l'accumulo; 6 byte in tutto).
  L'accumulo parte da 0 e prosegue tra i segmenti; l'ultimo segmento puo' essere piu' corto.
- Valore decodificato / **2000** = mAU. Verifica: massimo a 4,357 min, 4201053 -> 2100,5 mAU; il
  report dice altezza 2003,05 mAU, cioe' al netto della baseline di integrazione (circa -97 mAU li').
- Tempo del punto i: `(inizio + i * 400) / 60000` min.
- Confronto con i REPORT (test `tests/test_chemstation_reader.py`): per i picchi >= 50 mAU entro 30 min
  l'apice del segnale cade a meno di 0,01 min dall'RT del report e `segnale / altezza` sta tra 1,00 e
  1,15. Per i picchi piccoli la baseline di ChemStation sta **sotto** il minimo locale del segnale e
  oltre 80 min e' negativa (-10/-40 mAU): per confrontare altezze e aree servira' replicare
  l'integratore. Molti picchi sono spalle (apice non coincidente con l'RT).
- Aperto: il fattore 1/2000 e' dedotto, non letto da un campo dell'intestazione. Altre versioni
  `.ch` (130/179/181) hanno layout diverso: confermare su altri campioni.

## Formato `dad1.uv` versione 31 (decodificato e verificato)
Il cubo DAD completo: uno spettro per ogni istante, qui **13493 spettri x 201 lunghezze d'onda
(200-600 nm, passo 2 nm)**, ogni 400 ms. I `dad1A..E.ch` sono solo profili da esso estratti
durante l'acquisizione (`Sig=280,4 Ref=360,100` = media 278-282 nm meno media 310-410 nm).
**Little-endian** (al contrario dei `.ch`).
- `0x000`: stringa Pascal con la versione (`"31"`), poi intestazione; i record partono da `0x200`.
- Record (uno per spettro, lunghezza variabile 424-452 byte): `uint16` tipo (`0x43`), `uint16`
  lunghezza **totale** del record, `uint32` tempo (ms; 240, 640, 1040...), `uint16` lambda iniziale,
  finale e passo in 1/20 nm (4000, 12000, 40), 8 byte non decodificati (contengono `0x50/0x51`,
  `0x04`, `0x190`=400 ms), poi i valori.
- Valori: stessa idea dei `.ch`, little-endian: `int16` = delta dal valore precedente (parte da 0 a
  ogni record), `0x8000` + `int32` = valore assoluto. Tutti i 13493 record danno esattamente 201 valori
  e consumano esattamente i byte dichiarati. Diviso **2000** = mAU.
- Dopo l'ultimo record c'e' una coda di ~135 kB (forse un indice, ~10 byte per spettro): non letta.
- **Tempo**: il `.uv` e il `.ch` usano lo stesso asse in ms; il `.ch` parte a -2560 ms, il `.uv` a
  240 ms, quindi `indice_ch = indice_uv + 7`.
- **Estrazione di una banda** (`estrai_lambda`): media pesata, peso = sovrapposizione fra la banda
  `centro +/- banda/2` e l'intervallo di 2 nm di ciascun punto spettrale, meno lo stesso per il
  riferimento. Verifica sul campione: scarto RMS rispetto a `dad1A..D.ch` di 0,17 / 0,03 / 0,01 /
  0,02 mAU (A = 280 nm, picco principale a 2100 mAU; scarto massimo ~7 mAU = 0,3%).
- Per ottenere quel risultato l'asse delle lunghezze d'onda va spostato di **+0,45 nm**
  (`UV_OFFSET_NM`): valore **empirico**, minimo dello scarto fra 0 e 1 nm, non letto dall'intestazione.
  Senza spostamento lo scarto RMS sale di circa 6 volte (A 1,1; B 1,8 mAU). Origine non chiarita
  (centro del pixel? calibrazione della lampada?): da riverificare con altri campioni.
- `leggi_uv()` ci mette ~1 s per l'intero file (percorso veloce numpy per i record senza valori assoluti).

## Metadati (cosa e dove)
Letti da HPLCManager (`leggi_metadati_cartella`, `leggi_ch`, `leggi_uv`):
- intestazione binaria di `.ch`/`.uv` (stringhe Pascal a offset fissi 0x18 campione, 0xB2 data, 0xD0 modulo,
  0xE4 metodo; `.ch`: anche la descrizione del segnale "DAD A, Sig=280,4 Ref=360,100");
- `Report00.CSV`: campione, file dati, strumento, metodo e modifiche, data iniezione, riga sequenza, vial,
  volume iniettato, file di sequenza, pressione e flusso (inizio/fine), solventi A-D, elenco segnali. La riga
  `Acq. Operator` contiene la data dell'iniezione, non l'operatore (si scarta);
- `RUN.LOG`: temperatura colonna (min-max sui messaggi del termostato).

**Non letti**: `ACQRES.REG` (registro binario: moduli 1100 con numeri di serie e firmware, colonna
"ODS Hypersil", ecc.), `LCDIAG.REG` (465 kB di diagnostica), il resto di `RUN.LOG` (pressioni nel tempo,
eventi), `Report.TXT` (ripete i REPORTnn.CSV), gli 8 byte non decodificati di ogni record `.uv`, la coda
del `.uv`, il file del metodo `.M` (non e' nella cartella `.D`).

## Da fare / da chiedere
- Confrontare area e altezza calcolate da HPLCManager con `REPORTnn.CSV` (tolleranza da decidere).
- Verificare gli spettri DAD con altri campioni (altri intervalli, passo 1 nm, altre versioni `.uv`).
- Purezza di picco (confronto degli spettri sul picco), libreria di spettri.
- Integrazione: replicare l'integratore ChemStation o farne uno proprio? Con quali parametri
  (soglia, pendenza, larghezza minima)?
- Chiedere a Francesco: tipo di analisi (qui polifenoli, gradiente 90 min), quantificazione con
  standard esterno, system suitability.

## Decisioni
(nessuna ancora)
