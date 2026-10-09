# Status of the candidate validation data (after the first fetch attempt, 2026-10-09)

Fetched by opencode following `docs/opencode_fetch_prompt.md` and uploaded for analysis. Checked against the files received, **not** against opencode's summaries (some of which contain errors, noted below).

| id | received | verdict |
|---|---|---|
| `azcona2017` | open-access PDF (35 pages), opencode summary and test matrix | **used**: `docs/azcona_validation.md`. Its test matrix lists case 3 as "harmonic" (the paper calls 1.58 s the Snap Condition for both configurations) and d = 19.870 (Table 1: 19.872); the paper's numbers were taken from the PDF. |
| `lopezolocco2022` | summary only; **no paper** (MDPI blocked the download; the raw file is a 0-byte placeholder) | not usable yet. The paper is open access (CC BY 4.0): it needs a manual download of the PDF. The summary's "findings from citing papers" is second-hand and was not used. |
| `nlr_mooring_1_100` | `Dastaset Description.pdf`, `Metadata.xlsx` (sheets: Channel List, Tank Testing Matrix, Spring Constants Test Matrix, Configuration Properties, Load Cell Layout), saved MHKDR page, CSV export of the configuration sheet; **not** the 4.75 GB data archive | metadata only. The opencode summary's table of line properties has shifted columns for the semi-taut configurations (missing `-` entries): read the CSV/xlsx, not the summary. See below for a targeted download. |
| `tud_oc5_taut_leg` | metadata page and small scripts; data download failed (HTTP 503) | no data. Retry later or download manually from the 4TU dataset page. |
| `tud_metsch_schreier` | metadata page and small scripts; data download failed (HTTP 503) | no data. The "10.9 s / 5.01 N/m" figures remain **unverified**. |
| `oc5_phase2` | model-definition PDF (2.5 MB), opencode summary | platform and mooring definition only (no measured time series); useful for building the OC5 configuration, not yet used. Its licence was not confirmed. |

## Targeted download of the NLR (MHKDR 677) data
`Data.zip` is 4.75 GB. Its structure (Dataset Description): one folder per mooring/load-cell combination (e.g. `CatenaryWithFloat2_FH`) with `ForceDisp`, `FreeDecay`, `WaveRuns`,
each holding `EDASS` (1000 Hz load cell + wave height) and `Qualisys` (240 Hz 6-DOF) data as raw (`.tdms`/`.tsv`) and processed (`.parquet`) files, plus `SpringConstants` and `Calibration Waves`.
A useful minimum: the processed `.parquet` files of `FreeDecay` and `WaveRuns` for `CatenaryWithFloat2_FH` and one nylon/polyester configuration, and `SpringConstants`. Wave cases: 132 mm / 1.6 s,
35 mm / 1.1 s, 17.5 mm / 0.9 s. Line properties are in the Configuration Properties sheet (e.g. total line length 436 cm, chain 420.5 cm for the catenary with float, rope diameter 0.078 cm).
