# BEM coefficients for the Paredes buoy hull

**These numbers are computed, not measured.** The NLR/TU Delft data hosts are not reachable from the build environment, so no real
data was obtained. Instead a linear potential-flow (BEM) solution for the hull alone was computed with Capytaine 3.0.0
(`scripts/bem_paredes_buoy.py`; result in `docs/data_paredes_bem.json`, 18 periods, 0.7-5.0 s). Nothing was tuned.

## Set-up
* Geometry reconstructed from the thesis: C33 = 2043.494 N/m gives D = 0.5150 m; draft 0.176 m; CG 0.1002 m below the still-water line;
  water depth 0.9 m; rho = 1000. The hull is a flat-bottomed cylinder (the exact thesis hull shape was not available), clipped at z = 0.
* Output is in the code's convention: moments about the CG, F(t) = w a sin(wt + delta) for eta = a sin(wt) at x = 0.
  Excitation = Froude-Krylov + diffraction (Capytaine's diffraction result alone is NOT the excitation; this was a bug in the first run).
* **Not done:** viscous damping, comparison with the true hull shape.

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

## CON2 and T = 1.5 s (added after the first write-up)
Same procedure; Table 3.4 run -> BEM run / measured. No measured RAO exists for T = 1.5 s (the thesis figures contain 1.3 and 1.4 s only), so those rows are model-to-model only.

| config, T | surge | heave | pitch/(ka) |
|---|---|---|---|
| CON2, 1.3 | 0.69 -> 0.92 / 0.40-0.45 | 1.29 -> 1.48 / 1.45-1.50 | 2.20 -> 2.35 / 1.7-1.95 |
| CON2, 1.4 | 0.78 -> 0.93 / 0.62-0.68 | 1.17 -> 1.30 / 1.30-1.48 | 1.64 -> 1.89 / 1.55-1.75 |
| CAT, 1.5 | 0.88 -> 0.94 / - | 1.03 -> 1.12 / - | 0.83 -> 1.54 / - |
| CON1, 1.5 | 0.84 -> 0.94 / - | 1.07 -> 1.17 / - | 1.06 -> 1.46 / - |
| CON2, 1.5 | 0.81 -> 0.93 / - | 1.11 -> 1.22 / - | 1.38 -> 1.63 / - |

* **CON2 is the counter-example to the pitch story.** Its Table 3.4 pitch (2.20 / 1.64) was already 13-17 % above the measured band, and the BEM
  coefficients raise it further (2.35 / 1.89, 20-38 % above at 1.3 s and 8-22 % above at 1.4 s). So the pitch result is a mixed improvement:
  CON1 and CAT move toward the data, CON2 moves away from it. Across the six measured cases the BEM pitch is now within or near the band for
  CON1 at 1.3 s (1.96 vs 2.1-2.3, -7 %) and CON1/CAT at 1.4 s are still 5-20 % low, while CON2 is high. The waterline-reference explanation therefore
  does not by itself make the model consistent; the mooring pitch stiffness/compliance differences between configurations are a likelier remaining factor.
* CON2 heave with the BEM coefficients is inside the measured band at both periods, and the surge excess persists (about a factor 2 at 1.3 s).
* At T = 1.5 s the BEM raises pitch strongly (+0.7 for CAT/CON1) because Table 3.4 combined with the CG-referenced stiffness gives a very low pitch there;
  there is no measurement to check against.

## Why is the CON2 pitch too high? (diagnostic study, T = 1.3 s, BEM coefficients)
CON2 is the only configuration whose pitch RAO is above the measured band (2.35 vs 1.7-1.95). All runs below use the same BEM table; the variants are
**diagnostics, not calibrations** (the base configuration files are unchanged).

| run | pitch/(ka) at 1.3 s | reading |
|---|---|---|
| CON2 as configured | 2.348 | |
| closed form, hull + quasi-static mooring stiffness `K` from `mooring_modes` | 2.047 | what the mooring does *statically* |
| same closed form, CON1 / CAT | 1.941 / 2.046 (time domain 1.957 / 1.892) | CON1 reproduced to 1 %; CAT is 8 % below its static value |
| CON2 with 46 instead of 23 segments | 2.340 | not a discretisation effect |
| CON2 with clump drag Cd = 0 | 2.371 | not clump drag |
| CON2 with the clump removed | 1.957 (= CON1 exactly) | the clump alone carries the whole excess |
| CON2 with zero drag on lines, floater and clump | 4.097 | drag is what limits the response |

Facts established:
* The hull is on a steep resonance flank: BEM pitch natural period 1.11 s against waves of 1.3-1.4 s, radiation damping ratio 1.3 %. Free-hull pitch RAO is
  2.14 (1.3 s) and moves to 2.97 for I55 = 1.0 instead of 0.87 kg m^2 (the inertia implied by the measured 1.170 s pitch period), i.e. the absolute pitch level is
  extremely sensitive to the poorly known pitch inertia and to the (absent) viscous pitch damping. This affects all configurations, not CON2 in particular.
* Quasi-static mooring stiffness is nearly identical in the three configurations (surge 22-27 N/m, surge-pitch coupling 3.7-5.0 N, pitch 38.9-39.5 N m/rad) and
  predicts pitch within 1 % for CON1 only. CON2's time-domain pitch is **15 % above** its quasi-static prediction, so the excess is not quasi-static (later shown to be amplitude dependent, see below).
* Line modes (`mooring_modes --lines`, undamped, rope law included): CON1 leg lowest mode **0.868 s**; CON2 leg modes **1.331 s and 1.500 s** (clump on the rope), then 0.691 s.
  The T = 1.3 s wave is 2.3 % from the 1.331 s mode; T = 1.4 s lies between the two modes.

~~Interpretation: clump-rope resonance at 1.33 s~~ **(retracted, see the period sweep and amplitude test below)**: the first write-up suggested that the 1.331 s line mode
is excited almost at resonance by the 1.3 s wave. The sweep below does not support this.

### Period sweep through the 1.331 s line mode (CON2, H = 0.08 m; CON1 as control)
BEM rows at the extra periods were computed with the same mesh (`scripts/bem_paredes_buoy.py out.json 1.0 <periods>`); `static-K` is the closed-form 2-DOF response with the
quasi-static linearised stiffness of `mooring_modes` (CON2: Kxx 26.63, Kxp 4.11, Kpp 38.94; CON1: 24.72, 3.72, 39.51).

| T [s] | CON2 time domain | CON2 static-K | ratio | CON1 time domain | CON1 static-K | ratio |
|---|---|---|---|---|---|---|
| 1.20 | 3.927 | 3.021 | 1.300 | | | |
| 1.25 | 2.850 | 2.387 | 1.194 | 2.259 | 2.240 | 1.009 |
| 1.28 | 2.511 | 2.162 | 1.161 | | | |
| 1.30 | 2.348 | 2.047 | 1.147 | 1.957 | 1.941 | 1.008 |
| 1.31 | 2.281 | 1.998 | 1.142 | | | |
| 1.33 | 2.167 | 1.911 | 1.134 | 1.834 | 1.819 | 1.008 |
| 1.35 | 2.072 | 1.838 | 1.127 | | | |
| 1.37 | 1.992 | 1.775 | 1.122 | 1.703 | 1.698 | 1.003 |
| 1.40 | 1.889 | 1.697 | 1.113 | 1.630 | 1.627 | 1.002 |
| 1.45 | 1.745 | 1.598 | 1.092 | | | |
| 1.50 | 1.628 | 1.524 | 1.068 | 1.462 | 1.469 | 0.995 |

* **There is no peak near 1.33 s.** The CON2 pitch falls monotonically with period, and the ratio to the quasi-static value increases smoothly toward short periods
  (1.07 at 1.5 s, 1.13 at 1.33 s, 1.30 at 1.2 s). A resonance of the 1.331 s line mode would give a local maximum of the ratio there; it is absent.
* The CON1 control is quasi-static to within 1 % at every period, so the hull-plus-stiffness model and the BEM table are not the source of the CON2 excess.
* A single frequency-independent reduction of the CON2 pitch stiffness of about **1.5 N m/rad (4 %; range 1.2-1.7 over the sweep)** reproduces all eleven time-domain values
  (equivalently +0.06...+0.08 kg m^2 of inertia). The excess therefore looks like a small softening of the effective pitch stiffness amplified by the steep resonance flank
  (hull pitch period 1.11 s), not a dynamic line effect at the wave frequency. (Diagnostic fit only; it is not used anywhere.)

### Amplitude test (T = 1.3 s)
| H [m] | CON2 pitch/(ka) | CON2 heave | CON1 pitch/(ka) |
|---|---|---|---|
| 0.02 | 2.128 | 1.393 | |
| 0.04 | 2.285 | 1.445 | 1.947 |
| 0.08 | 2.348 | 1.450 | 1.957 |
| 0.16 | 2.326 | 1.490 | |

The static-K value is 2.047. CON1 is independent of wave height (as a linear system would be); **CON2 is not**: its pitch/(ka) grows from 4 % above the quasi-static value
at H = 0.02 m to 15 % at 0.08 m and then saturates (heave varies too, 1.39 -> 1.49). The CON2 excess is therefore an **amplitude-dependent (nonlinear) effect of the clump leg**, not a
linear resonance. Candidate mechanisms, none tested: loss of tension (the line law is bilinear, T = max(0, ...)) or large changes of the leg geometry with the clump swinging at 13-23 deg of pitch;
drag nonlinearity of the floater/clump (the no-drag run, 4.10, shows the drag controls the level but the clump drag alone does not change it); and the small-angle rigid-body model itself, which is being used at
pitch amplitudes of about 0.2 rad. The tension time series below tests the first of these.


### Tension time series of the legs (T = 1.3 s; `paredes_report ... <cfg> - <bem.json> 1.30 <H> <out.csv>`)
`LumpedMassCable::minSegmentTension` (new diagnostic) gives, for every output sample (0.01 s), the smallest segment tension, its segment and the number of slack
segments (strain <= 0, i.e. zero tension in the bilinear law). Statistics over the last 15 wave cycles of the 30 s run:

| run | fairlead tension range [N] (3 legs) | smallest segment tension [N] | slack segments | pitch/(ka) |
|---|---|---|---|---|
| CON2, H = 0.08 m | 9.78 .. 12.84 | 5.62 / 5.98 / 5.98 | none (0 % of the time) | 2.348 |
| CON2, H = 0.02 m | 10.52 .. 11.53 | 6.22 / 6.28 / 6.28 | none | 2.128 |
| CON1, H = 0.08 m | 2.63 .. 2.94 | 2.62 / 2.63 / 2.62 | none | 1.957 |

(the static values are 10.9-11.0 N fairlead and about 6.3-6.8 N minimum for CON2, and 2.74 N for CON1).

* **Loss of tension is ruled out** as the cause: no segment of any leg slackens in the analysis window; the CON2 minimum is 5.6 N, i.e. about 90 % of its static value.
  The bilinear tension law is never active there.
* **Start-up artefact:** slack segments do appear in the first samples (CON2 up to t = 0.13 s, CON1 up to 0.04 s, up to all 23 segments). They come from the harness
  imposing the steady-state pose and velocity of the free body at t = 0 on lines still in their static shape; they are over by 0.13 s and the fit window starts at 10.5 s.
* **The response is stationary** in the fit window: the CON2 pitch amplitude is 0.229-0.233 rad cycle by cycle after the 4th cycle (the 15-cycle fit is not biased by a slow transient);
  CON1 fluctuates by about +-3 % cycle to cycle (0.185-0.202 rad), so its RAO carries a corresponding scatter.
* **The line force is nearly linear in the wave height:** the amplitude of the leg-1 fairlead tension is 0.36 N at H = 0.02 m and 1.29 N at 0.08 m (a factor 3.6 for a factor 4,
  i.e. about 10 % sub-linear), while the pitch RAO changes by 10 % between those heights. The nonlinearity of the CON2 pitch is thus small in the line force itself.
* Pitch amplitude at H = 0.08 m is 0.23 rad (13 deg) for CON2 and 0.19 rad for CON1; the rigid body of this code is a small-angle model.

What this leaves: the amplitude dependence appears only in CON2, whose pretension is four times that of CON1 (11 N against 2.7 N per leg). The pattern is consistent with a geometric
(large-rotation) effect of the pretension on the fairlead moment, which a small-angle body formulation does not represent, but this is **not tested** (it would need a large-angle
body formulation or a CON2 variant with the pretension scaled). The assumed rope/floater/clump drag (not given in the thesis) remains the only damping and the pitch level is
very sensitive to the pitch inertia and damping of the hull (see above). Status: the CON2 pitch excess is a nonlinear effect specific to the heavily pretensioned configuration;
mechanism not identified; loss of tension and clump drag are excluded.

### Scaled-pretension variants of CON2 (T = 1.3 s, BEM coefficients)
The CON2 pretension comes from the clump (submerged weight 9.97 N, mass 1.115 kg at 1.6 m); CON1 is the same rope without it. The clump was scaled by `s`
(**diagnostic variants, not calibrations**; the base files are unchanged): **A** scales submerged weight and mass together, **B** scales the submerged weight only (mass fixed), which separates pretension
from clump inertia. For each variant the quasi-static stiffness was recomputed with `mooring_modes` (the static stiffness does not depend on the clump mass:
A and B have identical `K`) and the closed-form 2-DOF pitch RAO of that stiffness is the reference.

| variant | leg-1 fairlead tension [N] | Kxx / Kxp / Kpp | H [m] | pitch/(ka) time domain | quasi-static | ratio | min segment tension [N] |
|---|---|---|---|---|---|---|---|
| CON1 (s = 0) | 2.63 .. 2.93 | 24.7 / 3.72 / 39.5 | 0.08 | 1.957 | 1.941 | 1.008 | 2.62 |
| A s = 0.25 | 4.13 .. 4.29 | 18.1 / 2.77 / 39.3 | 0.08 | 2.014 | 1.951 | 1.032 | 3.55 |
| A s = 0.25 | 4.02 .. 4.29 | same | 0.02 | 2.078 | 1.951 | 1.065 | 3.39 |
| B s = 0.25 | 3.51 .. 4.79 | same | 0.08 | 1.990 | 1.951 | 1.020 | 3.40 |
| A s = 0.5 | 5.94 .. 6.87 | 19.4 / 3.01 / 39.3 | 0.08 | 2.050 | 1.961 | 1.045 | 4.47 |
| A s = 0.5 | 6.17 .. 6.45 | same | 0.02 | 2.043 | 1.961 | 1.042 | 4.61 |
| B s = 0.5 | 5.55 .. 7.17 | same | 0.08 | 2.040 | 1.961 | 1.040 | 4.42 |
| A s = 0.75 | 7.89 .. 9.59 | 22.6 / 3.51 / 39.1 | 0.08 | 2.165 | 1.996 | 1.085 | 5.28 |
| CON2 (s = 1) | 9.78 .. 12.84 | 26.6 / 4.11 / 38.9 | 0.08 | 2.348 | 2.047 | 1.147 | 5.62 |
| CON2 (s = 1) | 10.52 .. 11.53 | same | 0.02 | 2.128 | 2.047 | 1.040 | 6.22 |
| A s = 1.5 | 13.48 .. 18.78 | 36.2 / 5.49 / 38.4 | 0.08 | 2.963 | 2.185 | 1.356 | 7.10 |

No segment is slack in the analysis window in any variant.

* **The excess over the quasi-static response is controlled by the pretension.** At H = 0.08 m the ratio time domain / quasi-static is 1.008, 1.032, 1.045, 1.085, 1.147 and 1.356 for
  s = 0, 0.25, 0.5, 0.75, 1, 1.5 (mean fairlead tension about 2.7, 4.2, 6.4, 8.7, 11, 16 N). For s >= 0.5 the excess (ratio - 1) grows roughly as the square of the pretension
  (exponent about 2.1-2.4 between neighbouring points; a rough reading of six points, not a fitted law); s = 0.25 lies above that trend.
* **The clump inertia is not the cause:** the B variants (mass fixed) give 1.990 and 2.040 against 2.014 and 2.050 for A (within 1-2 %).
* **The amplitude dependence only appears at full pretension:** for s = 0.5 the RAO is the same at H = 0.02 and 0.08 m (2.043 vs 2.050); for s = 0.25 it is 2.078 vs 2.014 (3 % in the opposite
  direction, which shows that amplitude-independent scatter of this size exists); for s = 1 it is 2.128 vs 2.348 (10 %). Together with the earlier wave-height test, the
  excess is amplitude- and pretension-dependent in a way that fits a nonlinearity that scales with (pretension x amplitude), the signature of a geometric stiffness effect
  at finite rotation. A residual of about 4-6 % above the quasi-static value remains at the smallest amplitudes and is not explained.

Status: the CON2 pitch excess is a **pretension-controlled nonlinear effect**; inertia, drag of the clump, discretisation, line resonance and loss of tension have each been excluded.
~~The geometric large-rotation explanation~~ **(refuted by the quasi-static pitch curve below)**.

### Quasi-static pitch curve: the geometric (large-rotation) explanation is refuted
`paredes_report examples/paredes pitchcurve <cfg>` holds the pitch at offsets up to +-0.25 rad from the equilibrium angle, relaxes the other DOFs (lines in static equilibrium) and prints the
line pitch moment. In this code the **fairlead positions already use the exact rotation** (`rotationFromVector`), so the geometric stiffness of the pretension at finite rotation is part of the
line force; only the body's own dynamics (linear A, B, C, omega = d theta/dt) are small-angle.

| theta [rad] | CON2 line moment [N m] | CON2 secant stiffness [N m/rad] | CON1 line moment [N m] | CON1 secant stiffness [N m/rad] |
|---|---|---|---|---|
| +0.05 | +0.0070 | 0.537 | -0.0588 | 1.188 |
| +0.10 | -0.0194 | 0.533 | -0.1180 | 1.187 |
| +0.20 | -0.0718 | 0.528 | -0.2369 | 1.187 |
| +0.25 | -0.0979 | 0.527 | -0.2968 | 1.190 |
| -0.05 | +0.0612 | 0.547 | +0.0603 (offset 0.0006) | 1.193 |
| -0.10 | +0.0893 | 0.554 | +0.1204 | 1.198 |
| -0.25 | +0.1794 | 0.582 | +0.3053 | 1.219 |

(CON2 moment at theta = 0: +0.0338 N m, CON1: +0.0006 N m; tangent stiffness 0.541 and 1.19.) The secant stiffness varies by at most 2-8 % over +-0.25 rad (0.1 N m/rad or less), and the
**total** line pitch moment at 0.25 rad is 0.1-0.3 N m. The CON2 time-domain excess is equivalent to a pitch-stiffness reduction of about 1.5 N m/rad (0.35 N m at 0.23 rad), i.e. several times the
entire line moment. A large-angle rigid-body formulation cannot therefore be expected to remove it: **the geometric/large-rotation hypothesis is refuted, and the earlier statements that it was the most
likely explanation (here, in the paper and in the roadmap) were wrong.** (Not covered by this test: simultaneous surge, heave and pitch of the actual motion, and the hull's own nonlinear hydrostatics, which affect
every configuration equally and so cannot create a CON2-specific excess.)

What remains: the excess is a **dynamic** line-body interaction (it is absent quasi-statically), grows with pretension (about its square) and with wave height, saturating at the highest height tested, and
is not due to slack, clump inertia, clump drag, resonance of the 1.331 s line mode, or discretisation. The candidates left are the dynamics of the lines themselves (inertia, added mass and
tension-dependent transverse stiffness at 1.3 s, which a quasi-static stiffness omits) and the quadratic drag of lines, floater and clump (without any drag the pitch is 4.1). The next
discriminating tests are the linear dynamic stiffness of the legs at 1.3 s from the linearised line model (drag-free) and time-domain runs with drag, or inertia, switched off one at a time.
