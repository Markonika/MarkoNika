# Frequency-dependent added mass and radiation damping (Cummins equation)

## Model
With frequency-dependent hydrodynamic data A(omega), B(omega) (6x6, from a BEM/diffraction solution) the body equation is the Cummins equation

    (M + A_inf) xi'' + int_0^t K(t - tau) xi'(tau) dtau + B_extra xi' + Dq |xi'| xi' + C xi = F(t),

    K(t) = (2/pi) int_0^inf B(omega) cos(omega t) domega,                  (retardation function)
    B(omega) = int_0^inf K(t) cos(omega t) dt,   A(omega) = A_inf - (1/omega) int_0^inf K(t) sin(omega t) dt   (Ogilvie relations).

* `body.hydro_table` = {`rows`: [{`omega_rad_s` or `period_s`, `A` / `A_diag`, `B` / `B_diag`}, ...], `A_inf` / `A_inf_diag` (optional), `retardation_length_s` (30), `retardation_step_s` (0.01), `high_frequency_tail` (true)}. The body's own constant `A` (and `B`, `Dq`) stay as additional terms: A_inf from the table is added to `body.A` in the mass matrix.
* **Retardation function** (`RadiationModel`): K_ij(t) for every pair with a non-zero B_ij, on a uniform grid of step `retardation_step_s` up to `retardation_length_s`, by the exact integral of the *piecewise-linear* B(omega) (closed form per segment, no quadrature in omega), with B(0) = 0 and above the highest tabulated frequency B_N (omega_N/omega)^2 (tail integral through the sine integral Si; switch off with `high_frequency_tail`). Linear interpolation in t; zero beyond the length.
* **A_inf**: from the table, or evaluated from the Ogilvie relation at the highest tabulated frequency, A_inf = A(omega_N) + (1/omega_N) int K sin(omega_N t) dt (exact integral of the piecewise-linear K).
* **Time stepping** (`RadiationMemory`, explicit coupling): the memory force I(t_{n+1}) = int K(t_{n+1} - tau) v(tau) dtau is evaluated with the trapezoid rule over the stored velocity history (non-uniform steps allowed; entries older than the kernel length are dropped). The zero-lag part (dt/2) K(0) v_{n+1} is treated implicitly together with the damping in the velocity-Verlet update, the rest of I enters the right-hand side. Cost: history length / dt per non-zero channel pair per step (e.g. 600 points for a 3 s kernel at dt = 5 ms).
* **Linearisation**: `bodyModesRadiation` (used by `mooring_modes` when a table is present) solves omega_k^2 (Mrb + A(omega_k)) x = K x per mode by fixed-point iteration on omega_k, with A(omega) interpolated from the table.

## Verification (`tests/test_radiation.cpp`)
Test data with an exactly known kernel: K(t) = k0 exp(-a t), i.e. B(omega) = k0 a/(a^2 + omega^2) and A(omega) = A_inf - k0/(a^2 + omega^2) (k0 = 200, a = 5, A_inf = 8), tabulated to 300 rad/s with spacing 0.1.
| Check | Result |
|---|---|
| K(t) from the table vs k0 exp(-a t), t = 0 ... 2 s | within 1e-4 k0 (limited by the table); K(0) = 200 |
| A_inf from the table (Ogilvie at 300 rad/s) vs 8 | 8.00000 |
| A, B from the kernel (Ogilvie) at 1, 3, 6, 12 rad/s vs the tabulated functions | within 4e-3 A_inf / B(0) |
| Time-domain heave (M = 10, C = 900) forced at 3 rad/s vs the frequency-domain steady state with A(omega), B(omega) | 3.5e-4, 1.1e-4, 5.4e-5 of the amplitude for dt = 40, 20, 10 ms (second order down to a floor of ~4e-5 set by the tabulated kernel) |
| Same near resonance (7 rad/s): amplitude with A(omega), B(omega) vs constant A_inf without damping | 0.281 m vs 2.78 m (the frequency dependence is not a small correction) |
| Free body in a 30-component JONSWAP sea through the JSON runner (table rows at 0.5 rad/s spacing above 0.5 rad/s, A_inf given, force table): heave vs the sum of the component steady states with the kernel's own A(omega), B(omega) | 1.1e-5 m for a 3.1e-2 m amplitude (4e-4) |
| Natural heave frequency: constant A_inf vs A(omega) (fixed point) vs a direct root-find on the same table | 7.071 vs 7.600 rad/s; root 7.60008 rad/s, equal to 6 digits |
| Implicit coupling together with the memory force | rejected (logic error) |

## Limitations
* The data must come from a diffraction/BEM solution; there is none for the Paredes buoy (its thesis gives A and B at one frequency per configuration), so this feature is verified against an analytic kernel, not against measurements.
* Explicit coupling only. The velocity history before `begin()` is taken as zero: starting from a steady-state velocity (as in `paredesWaveRAO`) gives a start-up transient of the memory force.
* The retardation length must cover the decay of K (the table-derived kernel is truncated, not windowed); a table that stops where B is still large gives Gibbs-type errors in K at short times (use the high-frequency tail and extend the table).
* Coupled surge-pitch terms etc. are handled through the full 6x6 pairs; the code assumes a body-fixed (not time-varying) set of coefficients and linear motions.
* Forward-speed and second-order terms are not included.
