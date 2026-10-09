# Numerical code for stabilising floating offshore wind turbines — design

Status: Layer 1 core implemented and tested; Layers 2–4 are design only.
Everything below the "verified" heading is checked by `tests/`; everything else is intent.

## 1. Purpose and scope

Compare stabilisation methods for floating wind turbine platforms, starting from the
gyrostabiliser (passive → PID → MPC → sliding mode / learning) and later the ATMD, under one
plant model so controllers are compared on equal terms.

Gaps this code is meant to address, taken from the literature work so far (each is "none found"
in my searches, not proven absent):

1. A wave-memory (state-space radiation) model inside the controller's predictor.
2. Interaction between the stabiliser and the blade-pitch loop / negative aerodynamic damping.
3. Sliding-mode or learning-based precession control with explicit actuator limits.
4. Energy accounting from a model rather than a heuristic penalty.
5. Several random seeds and a common baseline for controller comparison.

Closest prior work to differentiate from: Wang et al., Ocean Eng. 343 (2026) 123190 (MPC with
ATMD + gyro, linearised table, Morison waves, no controlled-case validation).

## 2. Layered plan

| Layer | Plant | Controllers | Purpose |
|---|---|---|---|
| 1 (done) | Platform pitch + gyro array, 4 states | Passive, saturated rate feedback | Test integrator, gyro model, wave generator |
| 2 | Pitch + roll + surge, state-space radiation memory, ATMD option | Add PID, MPC | Wave memory, 3D gyro coupling, QP controller |
| 3 | Add flexible tower modes, quasi-steady thrust, blade-pitch PI | Add sliding mode, learning | Aerodynamic damping and blade-pitch interaction |
| 4 | 13-DOF Kane model, BEM aerodynamics | All | Cross-check against OpenFAST; design optimisation (NSGA-II) |

Rule: a layer is not extended until the one below passes its verification tests.

## 3. Models

### Layer 1 (implemented): `plant.py`
State `[phi, phi_dot, alpha, alpha_dot]`, n identical gyros lumped into one equivalent unit:

    I_F phi''  + C_F phi' + K_F phi = M_wave - n Izz w_g alpha' cos(alpha)
    n Ixx alpha'' + n Cp alpha'     = n Izz w_g phi' cos(alpha) + u

Form follows Wang C. et al. (Ocean Eng. 329, 121147, Eq. 5). Precession angle hard stop at ±70°
applied after each RK step. Parameters are placeholders shaped like OC4-DeepCwind (Izz 1.5e5,
I_F 8.35e9, T_n 25 s); `K_F` is set only to give the 25 s period.

### Layer 2 (design): wave memory
Cummins equation with the radiation convolution replaced by a state-space fit:

    (M + A_inf) x'' + C x' + K x + ∫ h(t-τ) x'(τ) dτ = F_exc
    →  z' = A_r z + B_r x',   F_rad = C_r z          (fit order 4–8 per DOF pair)

Fit method: frequency-domain least squares / vector fitting to the BEM added mass and damping
(Capytaine or NEMOH). Check: reproduce the BEM frequency response and the convolution result in
the time domain. Reduce with balanced truncation if the MPC state grows.

### Layer 3 (design): aerodynamics and the blade-pitch loop
Quasi-steady thrust `T(V_rel, θ)` from a table, linearised per operating point, with the PI
blade-pitch loop closed. Checks: negative-damping onset when pitch gains rise (Nyquist margin as
in Yu et al. 2024, Fig. 2), and the effect of adding the stabiliser on that margin.

### Layer 4 (design): full model
13-DOF Kane model as in the 2025/2026 Wang papers, BEM aerodynamics, mooring by quasi-static
lookup first. Cross-check against OpenFAST with the stabiliser *off* and then *on*, since the
published validations covered the off case only.

## 4. Numerics

| Item | Choice | Reason |
|---|---|---|
| Time stepping, now | Fixed-step RK4, dt = 0.01 s | 4th order confirmed by test; fast enough |
| Stiff additions (lumped-mass mooring, actuator lags) | Implicit (BDF/Radau) or generalised-α | RK4 stability limit |
| Controller | Own sample time Ts = 0.05 s, zero-order hold | Matches MPC literature |
| Non-smooth elements | Hard stop outside RK stages; tanh/boundary layer for sliding mode | Avoid numerical chattering |
| Waves | JONSWAP, jittered frequencies, random phases, multiple seeds | Avoids record repetition; statistics |
| MPC (Layer 2) | Linearise (automatic differentiation), exact discretisation, QP via OSQP, warm start; RTI nonlinear MPC if cos(alpha) matters | Real-time capable |
| State estimation | Kalman / EKF on measured pitch, rate, precession | Wave moment unmeasured |
| Optimisation (Layer 4) | NSGA-II, reduced model in the loop | As in Yu et al. 2024 |

## 5. Verified in Layer 1 (`tests/test_plant.py`, 6 tests passing)

- Energy conserved to 1e-6 relative (n = 1, no damping, no torque): the gyroscopic coupling is
  skew-symmetric, so the model has no spurious dissipation.
- Free-decay period with the gyro off matches the target within 0.1%.
- RK4 error ratio between dt = 0.02 and 0.01 lies in 12–20 (≈16 for 4th order).
- Precession hard stop is respected.
- JONSWAP reproduces the target Hs.
- A passive gyro reduces pitch in irregular waves.

### Plausibility sweep (`examples/sweep_passive.py`)
Wave moment is calibrated so gyro-off pitch std is 0.5°. Result (3 seeds):
reduction grows with spin speed with diminishing returns (6000 rpm 86%, 12000 rpm 97%), is poor
at very low precession damping (Cp 4e3–4e4: 11–18%), and falls at high damping (1.4e6: 77%).
Qualitatively this matches Wang C. et al.; numbers do **not** reproduce theirs (e.g. their
4e4–2e5 range is near-optimal and 1.4e6 gives 25%; here 4e4 is poor). The wave-moment transfer
function and platform parameters are placeholders, so do not read these values as results.

## 6. Known limitations of Layer 1

- Single DOF (pitch). No roll, surge, tower, wind, or mooring.
- Wave moment is a flat-times-roll-off placeholder, calibrated by one scale factor.
- Radiation is a constant damping ratio, no memory yet.
- No energy/power model yet (heuristic penalties in the literature are not used).
- `SaturatedRateFeedback` is a baseline of my own, not a published controller.

## 7. Open decisions (not yet answered)

1. Platform: semi-submersible first (assumed) or spar as well?
2. Goal: controller comparison on reduced models (assumed) or load prediction?
3. Language: Python/NumPy/SciPy assumed; CasADi and OSQP when MPC arrives (not installed yet).
4. Where this lives: currently `floating-stabilisation/` inside the profile repository. A
   separate repository would be cleaner if the project grows.

## 8. Next steps

1. Fit a state-space radiation model to BEM data (needs Capytaine; install and test offline).
2. Extend to pitch + roll with the counter-rotating pair's yaw cancellation checked explicitly.
3. PID precession controller, then MPC with constraints on alpha, alpha_dot, u.
4. Add the 15 MW spar parameters once the spar study is read in full.

## 9. Run

    cd floating-stabilisation
    pip install numpy scipy matplotlib pytest
    PYTHONPATH=src python3 -m pytest -q
    PYTHONPATH=src python3 examples/sweep_passive.py
