# Prompt for opencode: fetch the candidate validation data

Paste everything below the line into opencode, run from the root of a local clone of `markonika/markonika` (branch `claude/peaceful-clarke-vofpwc`), on a machine with normal internet access.

---

You are helping build validation data for a proprietary C++ mooring-line solver in `mooring/` (read `mooring/README.md` and `mooring/docs/validation_data_candidates.md` first). Your job is ONLY to fetch, store, inventory and summarise public experimental data. Do not change any source code, tests or configs in `mooring/src`, `mooring/include`, `mooring/tests`, `mooring/examples`.

## Rules
1. Everything lives in `mooring/data/external/<id>/` (ids below). Raw downloads are never modified. Add `mooring/data/external/raw/` and `*.zip` to `.gitignore`; commit only manifests, summaries and small derived CSV files.
2. For every file record in `mooring/data/external/MANIFEST.json`: source URL, final URL after redirects, retrieval date (UTC), size in bytes, SHA-256, HTTP status, and the licence text or licence URL shown on the page. Save a copy of the landing page as HTML next to the data.
3. The URLs in `validation_data_candidates.md` came from web-search summaries and were never opened. Verify each one. If a link is dead, wrong, needs a login/registration/request email, or the licence forbids redistribution: do NOT circumvent anything. Record the problem in `mooring/docs/validation_data_status.md` and continue with the next item.
4. Never invent numbers. Anything you cannot read from the source is `null` with a note. Values read off a figure by eye are marked `"source": "digitised"` with an uncertainty estimate; values from a table or a data file are marked `"source": "table"` / `"file"`.
5. Respect licences (keep attribution texts); do not scrape beyond the listed items; be polite with rate limits (one request at a time, small delay).
6. Do not run anything that tunes or fits the solver. No git push. At the end only `git add` the allowed files and commit locally; tell me what is large and left untracked.

## Items, in priority order
**A. `azcona2017`**: Azcona, Munduate, González, Nygaard (2017), Ocean Engineering 129, 415-427, DOI 10.1016/j.oceaneng.2016.10.051. Try the open copy https://oa.upm.es/38826/1/INVE_MEM_2015_213336.pdf. Save the PDF. Extract into `SUMMARY.md` and `test_matrix.csv`: chain properties (length, mass per length in air, wet weight, diameter or link dimensions, axial stiffness EA), basin depth, positions of the anchor and the suspension point, excitation (every amplitude and period), drag/added-mass coefficients used, and every measured tension value that is printed in a table (max, min, mean per case). List separately which results exist only as figures. If any supplementary data or a repository link is mentioned, follow it.

**B. `lopezolocco2022`**: Lopez-Olocco et al. (2022), J. Mar. Sci. Eng. 10(5), 676, DOI 10.3390/jmse10050676, https://www.mdpi.com/2077-1312/10/5/676. Save the PDF and the Data Availability Statement text verbatim. Extract: scale and line properties, clump-weight mass, size and position(s) along the line, excitation amplitudes and periods, measured fairlead tensions and dissipated energies from any table, and the instrumentation. Note where data are only in figures.

**C. `nlr_mooring_1_100`**: MHKDR submission 677, https://mhkdr.openei.org/submissions/677 ("Evaluation of 1/100th Scale Mooring Systems for Wave Energy Converters"). Download every resource file (the dataset description PDF, e.g. https://mhkdr.openei.org/files/677/Dastaset%20Description.pdf, `Metadata.xlsx`, and the data archives). Write `SUMMARY.md`: licence, test matrix (free decay, forced displacement, wave cases, spring-constant characterisation), the four mooring designs, channel names/units/sample rates for each file type, body mass/inertia/CG, fairlead and anchor coordinates, line properties from the sheet "Configuration Properties" (export every sheet to CSV in `derived/`), and whether any record contains a slack-snap event (peak tension vs mean). Do not unpack archives larger than 2 GB without asking.

**D. `tud_oc5_taut_leg`**: 4TU dataset https://data.4tu.nl/datasets/9f32fe66-ec63-4cad-bf6a-7e9fe62d68e8 (DOI 10.4121/9f32fe66-ec63-4cad-bf6a-7e9fe62d68e8). Same procedure: licence, file list, channels, sample rates, test matrix (regular, irregular, multi-sine), mooring configurations (linear springs and nonlinear nylon: record their force-extension data if given), platform mass/inertia/CG, fairlead/anchor coordinates, and whether free-decay records exist.

**E. `tud_metsch_schreier`**: 4TU dataset https://data.4tu.nl/datasets/a4ad6557-14a4-4c26-b699-a59314472f12 (DOI 10.4121/a4ad6557-14a4-4c26-b699-a59314472f12.v1). Same procedure. The page says processed data are shared on request by email: do NOT send any email; just record the instruction and the address in `validation_data_status.md`. Look for the thesis PDF and find the sections that give the surge natural period and mooring stiffness; quote them with page numbers (the brief for this project claimed "surge natural period 10.9 s, stiffness 5.01 N/m", unverified: confirm or contradict).

**F. `oc5_phase2`**: https://a2e.energy.gov/api/datasets/oc5/oc5.phase2/files/oc5.phase2.model.definition-semisubmersible-floating-system-phase2-oc5-ver15.pdf. Save the model-definition PDF. Extract platform mass/inertia/CG, fairlead/anchor coordinates, line properties and the list of available measured channels. Record whether the MARIN tension time series are downloadable and under what terms (registration is a stop condition, see rule 3).

**G. Optional, only if A-F are done:** the Zenodo record 10.5281/zenodo.3377120 (MARMOK-A-5, OPERA), and the Texas A&M inserted-springs paper (Applied Ocean Research / Ocean Engineering, 2001, "Dynamic analysis of mooring lines with inserted springs") for its digitised time series. Inventory only.

## Deliverables
* `mooring/data/external/<id>/` with `raw/`, `SUMMARY.md`, and `derived/*.csv` (tables you extracted), plus the saved landing page.
* `mooring/data/external/MANIFEST.json` covering every file.
* `mooring/docs/validation_data_status.md`: one section per item with: reachable (yes/no, why), licence, what is in it, what the solver would need to use it (line properties, coordinates, excitation, measured quantities), what is missing, and a verdict: ready / needs author request / figures only / not usable.
* A short final message listing: what was downloaded (sizes), what failed and why, anything that contradicts `validation_data_candidates.md`, and a one-line recommendation for which dataset to compare first.
* A zip `mooring/data/external_raw.zip` of `raw/` folders that I can upload for analysis (untracked).
