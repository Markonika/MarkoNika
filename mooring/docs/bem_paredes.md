# BEM coefficients for the Paredes buoy hull

**These numbers are computed, not measured.** The NLR/TU Delft data hosts are not reachable from the build environment, so no real
data was obtained. Instead a linear potential-flow (BEM) solution for the hull alone was computed with Capytaine 3.0.0
(`scripts/bem_paredes_buoy.py`; result in `docs/data_paredes_bem.json`, 18 periods, 0.7-5.0 s). Nothing was tuned.

## Set-up
* Geometry reconstructed from the thesis: C33 = 2043.494 N/m gives D = 0.5150 m; draft 0.176 m; CG 0.1002 m below the still-water line;
  water depth 0.9 m; rho = 1000. The hull is a flat-bottomed cylinder (the exact thesis hull shape was not available), clipped at z = 0.
* Output is in the code's convention: moments about the CG, F(t) = w a sin(wt + delta) for eta = a sin(wt) at x = 0.
  Excitation = Froude-Krylov + diffraction (Capytaine's diffraction result alone is NOT the excitation; this was a bug in the first run).
* **Not done:** T = 1.5 s and CON2 RAO runs, viscous damping, comparison with the true hull shape.

## Comparison with thesis Table 3.4 (hull alone, T = 1.3 s)
| quantity | BEM | Table 3.4 |
|---|---|---|
| A11 / B11 | 26.04 / 35.68 | 24.55 / 31.73 |
| A33 / B33 | 28.18 / 39.54 | 26.23 / 37.94 |
| w1 / delta1 | 1104 / 1.468 | 1058 / 1.475 |
| w3 / delta3 | 863 / 0.267 | 818 / 0.258 |
| A55 / A15 / B15 (about CG) | 0.312 / 1.33 / 2.36 | 0.260 / 1.17 / 1.03 |
| w5 / delta5 (about CG) | 80.5 / 1.468 | 34.3 / -1.667 |

Surge and heave agree to about 5 % (A, B within 6-15 %). Pitch differs strongly **because of the reference point**: moving the BEM
moment to the still-water line (M_wl = M_cg + z_cg F1) gives |w5| = 30.2 / 25.6 / 22.0 and delta5 = -1.674 / -1.646 / -1.627 at
T = 1.3 / 1.4 / 1.5 s, against Table 3.4's 34.3 / 29.7 / 25.9 and -1.667 / -1.641 / -1.623. The phase matches to 0.01 rad, the magnitude is
12-15 % low. This strongly suggests the thesis pitch coefficients refer to the waterline, whereas this code's pitch DOF is about the CG.
The shifted A15 and B15 match in magnitude (8-17 %) but with the **opposite sign** to the table; this is unresolved.

## Moored wave RAOs with the BEM coefficients (H = 0.08 m)
`paredes_report examples/paredes waves <cfg> - docs/data_paredes_bem.json`. Format: Table 3.4 run -> BEM run / measured (read by eye, see `waves_validation.md`).

| config, T | surge | heave | pitch/(ka) |
|---|---|---|---|
| CON1, 1.3 | 0.69 -> 0.87 / 0.25-0.33 | 1.29 -> 1.45 / 1.45-1.50 | 1.63 -> 1.96 / 2.1-2.3 |
| CAT, 1.3 | 0.75 -> 0.85 / 0.25-0.33 | 1.22 -> 1.32 / 1.45-1.50 | 1.38 -> 1.89 / 2.1-2.55 |
| CON1, 1.4 | 0.77 -> 0.90 / 0.50 | 1.14 -> 1.26 / 1.40-1.45 | 1.26 -> 1.63 / 1.75-1.90 |
| CAT, 1.4 | 0.82 -> 0.89 / 0.46-0.50 | 1.09 -> 1.19 / 1.33 | 1.02 -> 1.66 / 1.85-2.05 |

## Findings
* **Heave and pitch move toward the measurements** (heave by 8-13 %, pitch by 20-60 %): CON1 heave at 1.3 s is now inside the measured band, and
  the CAT pitch deficit shrinks from 35-50 % to 10-20 %. The Table 3.4 coefficients (waterline-referenced pitch used about the CG) are therefore
  a plausible partial cause of the earlier pitch disagreement.
* **Surge gets worse** (over-predicted by a factor 2.6-3.5 at 1.3 s). The BEM does not explain the surge disagreement; the measured surge
  minimum is probably a viscous/mooring effect that inviscid constant-coefficient hydrodynamics cannot capture.
* Caveats: the BEM hull is a reconstruction, and the by-eye measured values carry their own uncertainty, so these are
  indications, not a validation. Real validation data (MHKDR 677 / OC5, 4TU) still requires the hosts listed in the session notes to be allowed.

## Mesh convergence
Same hull, `mesh_resolution_factor` 0.6 / 1.0 / 1.5 / 2.0 = 364 / ~1000 / 2232 / 3936 wetted panels (displaced volume 0.03636 / - / 0.03662 / 0.03664 m^3
against 0.03666 exact). Values at T = 1.3 s (similar at 1.4 and 1.5 s):

| quantity | 0.6 | 1.0 | 1.5 | 2.0 | change 1.0 -> 2.0 |
|---|---|---|---|---|---|
| A11 | 26.73 | 26.04 | 25.72 | 25.56 | -1.9 % |
| A33 | 28.10 | 28.18 | 28.26 | 28.30 | +0.5 % |
| B11 | 37.06 | 35.68 | 35.03 | 34.68 | -2.8 % |
| B33 | 37.79 | 39.54 | 40.30 | 40.69 | +2.9 % |
| A15 / B15 | 1.281 / 2.291 | 1.331 / 2.364 | 1.350 / 2.387 | 1.360 / 2.399 | +2.2 % / +1.5 % |
| w1 / w3 / w5 | 1115 / 858 / 79.1 | 1105 / 863 / 80.5 | 1099 / 865 / 80.7 | 1097 / 865 / 80.7 | -0.7 % / +0.2 % / +0.2 % |
| delta1 / delta3 | 1.464 / 0.264 | 1.468 / 0.267 | 1.470 / 0.268 | 1.470 / 0.268 | < 0.003 rad |

* Excitation is converged to < 1 %; added mass and damping still drift by up to 3 % between 1.0 and 2.0 (B33 and B11 move in opposite directions), and
  A, B are still changing slowly at 2.0, so the finest mesh is not fully converged in B (probably another ~1 %).
* **Effect on the moored RAOs** (CAT, mesh 1.0 -> 2.0, `docs/data_paredes_bem_fine.json`): T = 1.3 s surge 0.852 -> 0.857, heave 1.319 -> 1.321, pitch 1.892 -> 1.892;
  T = 1.4 s 0.893 -> 0.897, 1.185 -> 1.189, 1.660 -> 1.652. All changes <= 0.6 %, so the discretisation is not a cause of the conclusions above
  (the 10-20 % pitch and heave shortfalls and the surge excess); the hull-shape reconstruction and viscous effects remain unquantified.
  CON1/CON2 were not re-run on the fine mesh.
