# Milestone 4: Chalmers/SSPA 33 m chain validation

Reference: Bergdahl, Eskilsson & Palm (2016), JMSE 4(1):5; set-up as in Paredes (2016) sec. 3.5.2, Table 3.1.
Inputs: `examples/chalmers/chalmers_config.json` (from `chalmers_model_inputs.json`), measured maxima from
`chalmers_table7_max_tension.csv` (digitised, +/-5 % reading error).

## Model as run
N = 66 segments (l0 = 0.5 m), RK4, CFL 0.5 (the seabed contact limit dt = 1.113e-4 s governs), no internal damping
(c_int = 0), Cm = 3.8, Cdt = 0.5, Cdn = 2.5, D0 = D1 = 2.2 mm, Ks = 3 GPa/m, zeta_s = 1, mu = 0.3, v_lim = 0.01 m/s,
anchor on the seabed at x = 0, top end at rest (32.554, 0, 3.3) m, water surface z = 3 m, seabed z = 0.
Initial state: touchdown catenary, then static relaxation. Motion: circle of radius r in the vertical x-z plane
about the rest position, amplitude ramped by a cosine over 2 cycles, 10 cycles in total; the reported value is the
mean of the per-cycle maxima of cycles 5-9 (maxima taken from every time step, no filtering).

**Rotation direction.** The inputs do not state it. Figure 3.5 of the thesis shows the motor with the arrow on the
left of the circle pointing up, i.e. the top of the circle moves toward the right (clockwise in x-z, `direction = -1`).
**Disclosure:** the first full grid was run counter-clockwise (an arbitrary default) and over-predicted by ~10 %. A
spot run with the sense reversed came out much closer to the data, and only then was Figure 3.5 consulted to see
whether the thesis states the sense. The figure supports clockwise, so it is used, but the check was motivated by the
mismatch and the choice cannot be called blind. No coefficient (Cd, Cm, mu, EA, ...) was changed. The opposite
direction is reported below as a sensitivity because it matters a great deal.

## Results

| Check | Result |
|---|---|
| Static top tension (published 22.68 N) | N = 16: 22.300, 33: 22.549, 66: 22.653, 132: 22.685, 264: 22.709 N (66: -0.12 %, 132: +0.02 %) |
| Regression r^2, simulated vs measured, 30 cases (published MooDy DG: 0.98) | **0.985** (slope 0.952, intercept +2.2 N) |
| r^2 against the 1:1 line (no fit) | 0.984 |
| RMSE / mean bias | 1.41 N / +0.12 N |
| Cases within the +/-5 % reading error | 24 of 30; all 30 within 6.7 % |

Figure: `docs/figures/chalmers_scatter.png`; per-case values: `docs/data/chalmers_grid_results_cw.csv`.

Pattern in the residuals (not tuned away): small-radius cases (r = 0.075, 0.1 m) are under-predicted at short
periods (-4 to -6 % at T = 1.25-1.5 s) and over-predicted at long periods (+4 to +7 % at T = 2.5-3.5 s). Both are
the largest deviations of the grid and are at the edge of the reading error.

### Sensitivity (diagnostic, not used for the headline numbers)
| Variation | Effect on the 30-case grid or on spot cases |
|---|---|
| Opposite rotation direction (`motion.direction = +1`) | r^2 (regression) 0.973, r^2 vs 1:1 **0.806**, mean bias **+4.4 N (+9.8 %)**; e.g. r = 0.2, T = 3.5: 59.0 N vs 52.0 N |
| Cdn 2.5 -> 1.2 (clockwise, spot cases) | r = 0.2, T = 3.5: 51.96 -> 51.4 N; r = 0.1, T = 2.5: 37.8 -> 36.1 N (grid-wide effect not run) |
| Friction mu 0.3 -> 0 / 0.1 (counter-clockwise run) | 59.0 -> 63.2 N / 58.8 N at r = 0.2, T = 3.5 |
| Cm 3.8 -> 0 (counter-clockwise run) | <= 0.4 N in three spot cases (consistent with Bergdahl: Cm = 0 also fits) |
| c_int = 50 N s (counter-clockwise run) | 59.0 -> 58.8 N (T = 3.5), 74.7 -> 70.5 N (T = 1.25) |
| Static tolerances: span +/-5 mm, vertical +/-3 mm | static top tension +/-0.33 N, +/-0.05 N |
| EA 10 000 -> 5 000 / 20 000 N (static, top offset 0.2 m) | 41.1 -> 30.7 / 52.8 N: the result is very sensitive to EA |

The spot-case sensitivities of Cm, friction and c_int were measured with the (wrong-direction) counter-clockwise
set-up before the direction was identified, so they show relative size only.

## Time series (r = 0.2 m)
`docs/figures/chalmers_timeseries.png` (raw, unfiltered, 1 ms output). T = 3.5 s: smooth cycle, 0-52 N, with a
low-tension ripple of a few N while the line is slack and a visible kink as the anchor segment re-tensions. T = 1.25 s:
broad, noisy loading with sharp rises up to ~72 N and long intervals at (near) zero tension. This agrees qualitatively
with the description in the thesis (smooth vs. quick motion with sharp peak loads, ripple near zero tension).
**No measured time series was supplied (only the maxima), so no quantitative shape comparison was possible.**

## Convergence (clockwise; `docs/data/chalmers_convergence.json`, `docs/figures/chalmers_convergence.png`)
Mean maximum top tension [N]:

| | N = 16 | 33 | 66 | 132 | 264 |
|---|---|---|---|---|---|
| r 0.2, T 3.5 | 52.11 | 52.02 | 51.96 | 51.99 | 52.00 |
| r 0.1, T 2.0 | 40.95 | 40.88 | 40.90 | 40.87 | 40.84 |
| r 0.2, T 1.25 (snap) | 70.83 | 70.56 | 72.29 | 71.57 | 72.13 |

| dt [s] (N = 66) | 1.1e-4 | 5.5e-5 | 2.75e-5 | 1.4e-5 |
|---|---|---|---|---|
| r 0.2, T 3.5 | 51.964 | 51.964 | 51.962 | 51.967 |
| r 0.1, T 2.0 | 40.900 | 40.898 | 40.899 | 40.900 |
| r 0.2, T 1.25 (snap) | 71.61 | 71.23 | 71.35 | 71.25 |

* Smooth cases stop changing at N >= 33 (l0 <= 1 m): variation <= 0.15 % up to N = 264; dt has no visible effect
  (<= 0.02 %) below the soil-limited 1.1e-4 s.
* The snap case does **not** converge monotonically: N >= 33 scatters between 70.6 and 72.3 N (+/-1.2 %) and dt
  between 71.2 and 71.6 N (+/-0.3 %), because the maximum is a one-or-few-sample spike. Everything stays inside the
  5 % reading error, but a single run of this case cannot be trusted better than about +/-1.5 %.
* RK4 and velocity-Verlet agree (51.966 vs 51.964 N; 40.904 vs 40.900 N; snap case 72.49 vs 72.29 N at default dt).
* The static top tension converges slowly with N (-1.7 % at N = 16 to +0.1 % at N = 264, limit about 22.72 N).

## Things that disagree with the reference or are uncertain
* The simulated static tension converges to ~22.72 N, 0.2 % above the published 22.68 N (within the geometry tolerance).
* The headline r^2 depends on the rotation sense being clockwise. With the other sense the model over-predicts by 10 %.
* Nominal chain areas A1 (only used with Cm != 0) are assumed pi/4 D0^2; D0 = D1 = 2.2 mm is the link thickness, not
  an equivalent hydrodynamic diameter. Results are sensitive to Cdn (see above), so Cd and D0 are partly interchangeable.
* The measured data are digitised chart values (+/-5 %), so the r^2 comparison cannot discriminate small model errors.
* The initial phase of the circle and the ramp (2 cycles) are our choices; the maxima converge by cycle 1-2 after the ramp.
