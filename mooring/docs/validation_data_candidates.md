# Candidate experimental data for further validation

**Update 2026-10-09:** a first fetch was done; see `docs/validation_data_status.md` (Azcona et al. used in `docs/azcona_validation.md`; the TU Delft data and the Lopez-Olocco paper could not be obtained; the NLR data archive was not downloaded).

Compiled from web search results on 2026-10-08. **Nothing here was opened or downloaded**: every host (4TU, MHKDR, Zenodo, MDPI, UPM, TU Delft, DOE A2e)
was unreachable from the development environment (network policy). Licences, file formats, channel lists, sample rates and whether files are downloadable are
therefore **unconfirmed** unless a row says otherwise. Statements below are what the search results said, not checked facts. Links were returned by the search and not tested.

## Ranked by fit with the current code

| # | Data | What it would test | Status / what is known |
|---|---|---|---|
| 1 | Azcona et al. 2017, submerged chain, suspension point driven by horizontal harmonic motion, in a basin. [UPM open-access PDF](https://oa.upm.es/38826/1/INVE_MEM_2015_213336.pdf), [journal](https://www.sciencedirect.com/science/article/abs/pii/S0029801816304942) | Closest to the Chalmers validation and to snap loads: tension at the suspension point and line motion, including cases where the line loses and regains tension. The paper reports good prediction with a slight over-prediction of the maximum tension. | Free PDF reported. No separate dataset found: numbers would come from figures/tables. |
| 2 | NREL/NLR 1/100-scale WEC moorings, [MHKDR submission 677](https://mhkdr.openei.org/submissions/677) | Four designs (catenary with float, taut nylon, semitaut nylon, semitaut polyester); free decay, forced displacement, waves, line spring-constant characterisation. Tests platform coupling, the floater element, decay checks. | Reported public. Channels: wave height, line tension, 6-DOF of two bodies. Licence not confirmed (MHKDR is generally CC BY 4.0). Snap loads not confirmed. Line specs reportedly in `Metadata.xlsx`, sheet "Configuration Properties". |
| 3 | Lopez-Olocco et al. 2022, 1:30 chain with/without clumped weights, CEHIPAR. [paper](https://www.mdpi.com/2077-1312/10/5/676), [abstract](https://doaj.org/article/5e715e689c2040208a04b0c8eb7037c6) | Clump-weight point element under forced fairlead motion at several amplitudes and periods; fairlead tension, dissipated energy. | Open-access paper. Data-availability statement not seen. |
| 4 | TU Delft 1:96 OC5 semisubmersible, taut-leg mooring. [4TU dataset](https://data.4tu.nl/datasets/9f32fe66-ec63-4cad-bf6a-7e9fe62d68e8) (doi 10.4121/9f32fe66-ec63-4cad-bf6a-7e9fe62d68e8) | Platform coupling with linear springs and with **nonlinear nylon** lines; 6-DOF motion and three tensions in regular, irregular and multi-sine waves. | Reported public. Free-decay tests: not known. Needs hull excitation coefficients (see below). |
| 5 | Metsch and Schreier MSc data, same semisubmersible at 1:96. [4TU dataset](https://data.4tu.nl/datasets/a4ad6557-14a4-4c26-b699-a59314472f12) (doi 10.4121/a4ad6557-14a4-4c26-b699-a59314472f12.v1) | Low-frequency surge: mooring matching OC5 stiffness and the surge natural frequency; monochromatic and bichromatic waves; 6-DOF, three line forces, four wave gauges. | Reported CC BY 4.0, CSV; processed data shared on request by email. |
| 6 | OC5 Phase II, MARIN 1:50 DeepCwind. [model definition](https://a2e.energy.gov/api/datasets/oc5/oc5.phase2/files/oc5.phase2.model.definition-semisubmersible-floating-system-phase2-oc5-ver15.pdf), [lumped-mass validation paper](https://www.researchgate.net/publication/279520507_Validation_of_a_lumped-mass_mooring_line_model_with_DeepCwind_semisubmersible_model_test_data) | Lumped-mass model driven by prescribed fairlead motion vs fairlead tensions (the paper reports very good agreement). | Whether the tension time series are downloadable: not confirmed. |

## Weaker or unconfirmed
* [2023 slack-line tank study](https://www.sciencedirect.com/science/article/pii/S0029801823001877): fairlead and anchor tension at 100 Hz (as described in the search summary; the finite-volume snap-load paper at that link is a related numerical study).
* Shallow-water [mooring damping and snap study](https://asmedigitalcollection.asme.org/offshoremechanics/article-abstract/141/5/051603/476892/An-Experimental-Study-of-Mooring-Line-Damping-and?redirectedFrom=fulltext): wavemaker-driven fairlead, line geometry by image processing.
* 1:60 [thesis](https://memorial.scholaris.ca/bitstreams/be7419ac-6c62-4449-b107-011d622f7b19/download): forced-oscillation studless chain, slack-snap events in the buoy-assisted case.
* Coastlab 2024 [wave-tank database](https://proceedings.open.tudelft.nl/coastlab24/article/view/803): includes snap loads in mooring lines; hosting not confirmed.
* [MARMOK-A-5 OPERA data (Zenodo 10.5281/zenodo.3377120)](https://scienceportal.tecnalia.com/en/datasets/h2020-opera-project-mooring-system-experimental-data-from-marmok-/): field data (polyester and elastomeric tethers on an OWC device); hard to use for validation.
* Texas A&M forced-oscillation study of [lines with inserted springs](https://www.sciencedirect.com/science/article/abs/pii/S0141118701000232): reported to publish digitised time series; hosting not confirmed.
* Chalmers group: no open dataset found beyond the data already used ([Paredes et al.](https://research.chalmers.se/en/publication/236872), [taut-moored WEC study](https://research.chalmers.se/en/publication/223579)).

## Corrections and cautions
* The Metsch data are **not** a wave-energy-converter test: they are the OC5 floating wind semisubmersible at 1:96. The brief's "Metsch, 1:96, surge period 10.9 s, stiffness 5.01 N/m
  (thesis sec. 3.3)" figures were taken from the original brief and have **not been verified** against the thesis; treat them as unverified.
* Items 4-6 are floating-platform tests in waves. The hull excitation coefficients for those platforms (and added mass/damping) are not supplied by what was found; a BEM run
  (e.g. Capytaine, read as data only per the licensing rules) would be needed.

## Code features each data set would require
* Nonlinear / viscoelastic tension law for nylon and polyester (#2, #4): the user-defined tension-strain hook of the brief is **not implemented yet**; it would be built first.
* Non-uniform line properties along the length (inserted springs, mixed chain/rope lines): the cable currently has one EA, mass and diameter for its whole length.

## Suggested order and network access
Start with #1, #3, #2: closest to what the code already does, and #2/#3 are the cleanest test of the point-element and platform work.
Hosts to allow in the environment's network settings: `oa.upm.es`, `www.mdpi.com`, `mhkdr.openei.org`, `data.4tu.nl`, `a2e.energy.gov`, `zenodo.org`, `proceedings.open.tudelft.nl`.
