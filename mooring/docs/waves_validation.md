# Milestone 9 (optional module): wave kinematics and regular-wave RAOs

## What was added
* `WaveField` (`include/mooring/waves.hpp`): any number of Airy components (regular or a superposition), finite-depth dispersion, Wheeler or no
  stretching, optional cosine ramp. Elevation of a component at (x, y, t): `A sin(w t - k s + phi)`, s along the propagation direction (so the
  elevation at the origin is A sin(w t), the convention of Eq. 3.60). Velocity and acceleration at any point, zero above the instantaneous surface.
  `asWaterField()` feeds `Environment::water`, so the Morison drag/added mass of the lines **and** of the point elements use the wave kinematics.
* `CoupledSystem::setExcitation`: time-dependent body load. With `body.wave_force {w, delta}` and `waves` in the JSON the body is excited by
  `f_i = w_i A sin(w t + delta_i)` (Eq. 3.60; w in N per m of amplitude, N m per m for pitch), at the mean position only, exactly as in the thesis.
* JSON: `waves {depth_m, surface_z_m, stretching, ramp_time_s, components[{height_m|amplitude_m, period_s, phase_rad, direction_deg}]}`.
* `paredesWaveRAO` (+ `paredes_report ... waves`): moored RAO at T = 1.30 / 1.40 / 1.50 s with the Table 3.4 coefficients of the thesis.

## Verification (closed forms; `tests/test_waves.cpp`, all in the default suite except the moored runs)
| Check | Result |
|---|---|
| Dispersion relation (3 depths x 3 periods), deep- and shallow-water limits | residual < 1e-12; limits to 1e-9 / 2e-3 |
| Crest velocity: Wheeler = A w coth(kh) (SWL value), linear = A w cosh(k(A+h))/sinh(kh); bed A w/sinh(kh) | exact to 1e-8 |
| Wheeler kinematic surface condition: w at z = eta equals d(eta)/dt | 1e-6 |
| Linear theory: a = du/dt at a fixed point; continuity du/dx + dw/dz = 0 | 1e-6 |
| Wheeler acceleration (stretched linear value) vs exact du/dt | 3.6 % of A w^2 at SWL + A/2 (a known approximation of Wheeler stretching) |
| Superposition of two components in different directions; zero kinematics above the surface; ramp | exact |
| Free buoy in regular waves (Table 3.4 coefficients incl. surge-pitch coupling), simulated vs closed-form RAO | surge/heave/pitch agree to < 1e-3 (T = 1.30: 0.6603 / 1.3007 / 1.8992; T = 1.40: 0.7436 / 1.1493 / 1.4540) |
| Same, started from the closed-form steady state | max deviation from the analytic time history 8.6e-6 of the amplitude over 20 s |
Note on the free-body test: the coupled surge-pitch mode is so lightly damped that a run from rest is still 0.26 % off at 600 s and 3e-5 at 1200 s,
independent of dt (checked for dt = 2 ... 0.25 ms). The moored runs therefore start from the free-body steady state (a closed form), so that only
the (small) mooring-induced transient remains; for CAT at 1.3 s this gives the same RAOs as a 60 s run from rest with a ramp (0.749 / 1.223 / 1.376 vs 0.748 / 1.221 / 1.366).

## Validation against the thesis: moored buoy RAO, T = 1.3 and 1.4 s (`ctest -R waves_validation`)
Model: configurations of milestone 8, Table 3.4 coefficients (A, B incl. surge-pitch coupling about the CG, w, delta), H = 0.08 m, Wheeler stretching on the lines,
30 s from the free steady state, first-order amplitude by least-squares fit over the last 15 wave periods (the thesis also fitted the first-order component).
Surge and heave are normalised by the wave amplitude a, pitch by k a (as in Figs 5.22-5.24). Measured values are read **by eye** from the plots as the range between
the H = 0.04 m and 0.08 m panels (uncertainty about +-0.07 surge, +-0.1 heave, +-0.25 pitch; the thesis stresses the RAO decreases with wave height).

| config, T | surge model / measured | heave model / measured | pitch/(ka) model / measured |
|---|---|---|---|
| CON1, 1.3 | 0.69 / 0.25-0.33 | 1.29 / 1.45-1.50 | 1.63 / 2.1-2.3 |
| CON2, 1.3 | 0.69 / 0.40-0.45 | 1.29 / 1.45-1.50 | 2.20 / 1.7-1.95 |
| CAT, 1.3 | 0.75 / 0.25-0.33 | 1.22 / 1.45-1.50 | 1.38 / 2.1-2.55 |
| CON1, 1.4 | 0.77 / 0.50 | 1.14 / 1.40-1.45 | 1.26 / 1.75-1.90 |
| CON2, 1.4 | 0.78 / 0.62-0.68 | 1.17 / 1.30-1.48 | 1.64 / 1.55-1.75 |
| CAT, 1.4 | 0.82 / 0.46-0.50 | 1.09 / 1.33 | 1.02 / 1.85-2.05 |
Free hull (closed form, no mooring): 0.66 / 1.30 / 1.90 at 1.3 s and 0.74 / 1.15 / 1.45 at 1.4 s.

A priori acceptance criteria (fixed before the runs): heave within 25 %, pitch within 30 %, surge "of the right order". Outcome:
* **Heave: pass in 6/6** but systematically 10-21 % low (it is the free-hull value; the mooring barely matters in heave).
* **Pitch: pass in 3/6** (CON2 at both periods, CON1 at 1.3 s marginally); **fails** for CAT (35-50 % low) and CON1 at 1.4 s (-34 %).
* **Surge: no quantitative agreement**; over-predicted by 14 % (CON2, 1.4 s) up to a factor 3 (CAT, 1.3 s). The measured surge RAO has a sharp minimum
  near T = 1.3 s that the model (and even the free hull) does not reproduce.
The test asserts only the documented envelope of these known disagreements (it does not claim agreement).

Hypotheses, none tested (no BEM data or raw measurements available): (i) the Table 3.4 coefficients belong to the Chapter 3 buoy set-up (same hull and depth
but possibly different draft/ballast); (ii) the reference point of the surge-pitch coupling and of the pitch coefficients (Table 3.4 gives no reference point; the thesis
refers motions to the centre of flotation); (iii) viscous damping and second-order effects, which the thesis notes lower the measured RAO with wave height; (iv) the
mooring in the experiment (CAT chain on the bed, floater/clump motion) is more compliant in pitch than the model; (v) the by-eye readings. The sign convention of the coupling
terms is supported by the pitch results: with the tabulated signs the free-hull pitch RAO is 1.90 / 1.45 (measured 1.75-2.5 / 1.6-2.0), with the coupling sign flipped 0.47 / 0.34 and without
coupling 1.24 / 0.92.

## Limitations of the wave implementation (also in README)
* Linear (first-order) waves only; no second-order drift (the thesis attributes the experimental surge offset to it), no wave-induced mean loads on the body.
* Submergence of the line nodes (weight/buoyancy blend and hydrodynamic forces) uses the **still-water** level by default (optional `waves.instantaneous_surface`, tested below: negligible effect here); Wheeler stretching is applied
  to the kinematics of the submerged part only. Floater and clump forces share that approximation.
* Wheeler accelerations are the stretched linear values (about 4 % of A w^2 off the exact derivative near the surface).
* The body excitation is a constant-coefficient force at the mean position (no dependence on the instantaneous position, no frequency-domain memory); irregular waves
  would need w(omega) from BEM output (not implemented; the component list supports irregular kinematics for the lines).
* Only T = 1.30, 1.40, 1.50 s have body coefficients (Table 3.4). The measured RAO figures contain 1.3 and 1.4 s (and not 1.5 s).

## Sensitivity: instantaneous free-surface submergence (hypothesis tested and not supported)
Hypothesis (listed above as a limitation): using the still-water level for the submergence of line nodes and point elements causes some of the RAO disagreements. `waves.instantaneous_surface = true` (assumptions item 44) makes the blend follow the instantaneous surface. Same runs as the table above (steady start, H = 0.08 m), surge / heave / pitch RAO, still-water level -> instantaneous surface:

| config, T | surge | heave | pitch/(ka) |
|---|---|---|---|
| CON1, 1.3 | 0.691 -> 0.691 | 1.287 -> 1.287 | 1.631 -> 1.630 |
| CON2, 1.3 | 0.693 -> 0.693 | 1.287 -> 1.287 | 2.201 -> 2.201 |
| CAT, 1.3 | 0.749 -> 0.750 | 1.223 -> 1.226 | 1.376 -> 1.392 |
| CON1, 1.4 | 0.771 -> 0.771 | 1.141 -> 1.141 | 1.262 -> 1.261 |
| CON2, 1.4 | 0.775 -> 0.775 | 1.169 -> 1.169 | 1.639 -> 1.637 |
| CAT, 1.4 | 0.815 -> 0.815 | 1.092 -> 1.094 | 1.020 -> 1.032 |

The largest change is 1.2 % (CAT pitch), far below the 10-50 % disagreements, so the still-water-level approximation is not their cause (the lines are nearly fully submerged in these set-ups; the effect could matter for lines crossing the surface). The default stays off.
