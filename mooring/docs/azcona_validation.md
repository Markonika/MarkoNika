# Validation against Azcona et al. (2017): submerged chain, forced fairlead motion

Reference: Azcona, Munduate, González, Nygaard (2017), *Experimental validation of a dynamic mooring lines code with tension and motion measurements of a submerged chain*,
Ocean Engineering 129:415-427 (open copy: UPM repository). Tests at the Ecole Centrale de Nantes wave tank (50 x 30 x 5 m). Data used here: **tables of the paper and values read by
eye from its Figs 5 and 6** (`examples/azcona/measured.csv`). No raw time series were available.
Config files: `examples/azcona/conf1.json`, `conf2.json`; run: `scripts/azcona_compare.py OUTDIR`; tests: `tests/test_azcona.cpp` (`ctest -R azcona_validation`).

## Set-up taken from the paper (no values tuned)
Chain: length 21 m (incl. the fairlead load cell), 69 g/m, steel 7850 kg/m^3, equivalent hydrodynamic diameter 3.4 mm, EA = 3.4e5 N, Cm = 1.0, Cdn = 1.4, Cdt = 0.67 (Table 3).
Basin depth 5 m, anchor on the floor at the origin, fairlead at the free surface (z = 5 m) with mean horizontal distance d = 19.364 m (configuration 1) and 19.872 m (configuration 2; the text says
19.870, Table 1 says 19.872), horizontal sinusoidal motion in the plane of the catenary at T = 1.58, 3.16, 4.74 s. Seabed (Table 4): stiffness 20 N/m^2 per unit length (our Ks*D1 = 20,
i.e. Ks = 5882 Pa/m), damping 0.1 Ns/m^2 (our zeta = 0.0426 from 2 zeta sqrt(Ks D1 m_l)), friction 0.5. 60 segments, RK4, CFL 0.5, 12 cycles with a 2-cycle ramp; the maximum is the mean of the per-step
maxima of the last 4 cycles. Static start: touchdown catenary + relaxation.

### Interpretations and assumptions (not stated unambiguously in the paper)
* **Amplitude.** Table 1 lists "0.25 m". The axes of Figs 5 and 6 span only 19.24-19.50 m and 19.745-20.0 m, i.e. +-0.125 m about the mean, so 0.25 m is read as the peak-to-peak stroke and
  **0.125 m is used**. Supporting evidence from the model itself: with 0.25 m as the amplitude the configuration 2 / 1.58 s maximum would be 268 N against about 44.5 N measured.
* The paper's "structural damping 0.1 %" has no usable definition here: internal damping c_int = 0 in the baseline (sensitivity below).
* Seabed friction: the paper multiplies the friction coefficient by the vertical seabed reaction; this code uses the submerged weight with a velocity ramp (v_lim = 0.01 m/s, not given by the paper). The normal
  friction coefficient of 0.5 is not implemented (tangential only).
* Cm = 1 is used as the normal added-mass coefficient (added mass = Cm rho A (1+eps)), as in the paper's wording "added mass coefficient 1".
* The fairlead sits exactly at the still-water level, so its node is half submerged in the blend of this code.

## Results
**Static fairlead tension** (Table 5; paper's own computation in brackets):
| | measured | model | difference |
|---|---|---|---|
| Configuration 1 | 8.13 N | 8.199 N | +0.85 % (paper's code: 8.10 N) |
| Configuration 2 | 14.48 N | 15.020 N | +3.7 % (paper's code: 14.70 N) |
The static tension depends strongly on the anchor distance (7 N/m and 22 N/m): the configuration 2 difference corresponds to 2.4 cm. The paper itself adjusted the anchor position (5.305 m from the wall) to fit the shape;
the nominal d is used here. Quasi-static end tensions over the stroke: model 7.35-9.26 N and 12.6-18.5 N vs about 7.7-9.3 N and 13.5-18 N read from the grey lines of the figures.

**Dynamic maximum fairlead tension** (60 segments; measured values read by eye, +-0.5 N; the criterion fixed in advance was +-8 %):
| case | T [s] | max model | max measured | diff | min model | min measured |
|---|---|---|---|---|---|---|
| conf1 | 1.58 | 15.12 | 14.5 | +4.3 % | 3.36 | 3.0 |
| conf1 | 3.16 | 10.23 | 9.8 | +4.4 % | 6.47 | 6.1 |
| conf1 | 4.74 | 9.35 | 9.2 | +1.6 % | 7.29 | 7.3 |
| conf2 | 1.58 (snap, total loss of tension) | **57.37** | 44.5 | **+28.9 %** | 0.00 | ~1.5 |
| conf2 | 3.16 | 25.07 | 23.5 | +6.7 % | 7.14 | 7.5 |
| conf2 | 4.74 | 20.35 | 19.5 | +4.4 % | 10.85 | 11.0 |
Five of six cases meet the +-8 % criterion; the paper reports +4.5 % for its own code in both 1.58 s cases (here +4.3 % for configuration 1). The model's minima agree to about 0.4 N where the line stays taut and
reproduce the total loss of tension in configuration 2 (the measured plateau is about 1.5 N in the figure, the text says 0).

### The snap case (configuration 2, 1.58 s): the +29 % is not a converged number
| variation (one at a time) | max [N] |
|---|---|
| baseline: 60 segments, RK4, CFL 0.5 | 57.4 |
| 30 / 120 / 240 segments | 53.8 / 50.3 / 48.5 |
| RK4 CFL 0.25 (60 seg.) / CFL 0.1 (60 seg.) | 59.8 / 56.3 |
| velocity-Verlet, 60 / 120 segments | 66.8 / 54.6 |
| RK4 CFL 0.25, 120 segments | 54.4 |
| d = 19.870 (text value) / d = 19.848 (reproduces the measured static tension 14.48 N) | 54.7 / 46.5 |
| friction 1.0 / 0.2 | 49.9 / 59.4 |
| internal damping c_int = 5 N s | 48.3 |
| Cdn 2.4 (wire diameter instead of 1.4) / Cdt 1.15 / Cm 0 | 98.0 / 49.1 / 54.0 |
| soil stiffness x10 / v_lim 0.001 | 58.9 / 59.0 |
The discretisation alone spans about 48.5-66.8 N (the best-converged runs give 48-55 N, i.e. +9 to +24 % of the measured value), and physical inputs the paper does not pin down (anchor position by 2.4 cm,
friction, internal damping) move the result by -10 to -20 %. The case is therefore consistent with the measurement within the model's own uncertainty, but this study cannot claim agreement better than about +10-25 %.
Configuration 1 at 1.58 s converges cleanly (14.99, 15.12, 15.13, 15.15 N for 30, 60, 120, 240 segments, identical for RK4 and Verlet).
Nothing in the table was adopted as a calibration. The Cdn sensitivity shows how strongly the snap peak depends on the (guideline) drag coefficient, as the paper also notes.

## Not covered
* The paper's marker trajectories (Figs 10-25, eight positions along the chain) were not compared (figures only, no tabulated data).
* The computed tension loops (tension vs fairlead position) were compared through maxima, minima and quasi-static ends, not point by point.
