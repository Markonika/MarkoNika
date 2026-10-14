# Irregular waves and depth-varying currents

## Irregular waves
* **Spectra**: JONSWAP (`type: "jonswap"`, peak enhancement `gamma`, default 3.3) and Pierson-Moskowitz (`type: "pm"`, gamma = 1) in terms of Hs and Tp, S(w) = norm w^-5 exp(-1.25 (wp/w)^4) gamma^exp(-(w - wp)^2/(2 sigma^2 wp^2)), sigma = 0.07 (w <= wp) / 0.09 (w > wp). `norm` is fixed numerically so that the zeroth moment is exactly Hs^2/16 (`JonswapSpectrum`). The gamma = 1 form equals the closed form 5/16 Hs^2 wp^4 w^-5 exp(-5/4 (wp/w)^4) to 1e-4 (the truncation of the numerical integral).
* **Discretisation** (`WaveField::addSpectrum`): `n_components` (default 50) bins of equal width over [omega_min_factor, omega_max_factor] x wp (default 0.5-4), one Airy component per bin at a uniformly random frequency inside the bin (`seed`, deterministic) with amplitude sqrt(2 S dw) and a random phase. The random frequency inside the bin avoids the artificial repetition of equally spaced lines; the amplitudes are rescaled so that sum(a^2/2) = Hs^2/16 holds exactly.
* **Kinematics**: the existing linear superposition (Wheeler stretching optional, finite depth). The evaluation was optimised (one sincos per component, one exponential for cosh/sinh, precomputed direction cosines): 2.7 us per node evaluation with 100 components against 10.5 us before (0.1 us with one component). A 60-segment line with 50 components costs about 3 s of wall time per simulated second with explicit RK4 (dt 1.2e-4 s); use fewer components (30 is usually enough for the line kinematics) or the implicit scheme for long records.
* **Body excitation**: the first-order load of the body now accepts a **frequency-dependent table**, `body.wave_force_table` = [{"period_s", "w": [6], "delta": [6]}, ...] (w in N per m of amplitude; N m for moments; delta the phase in rad), interpolated linearly in omega at every component: F_i(t) = ramp(t) sum_k w_i(omega_k) a_k sin(omega_k t + phi_k + delta_i(omega_k)). A component outside the table range is an error (no silent extrapolation); the load coefficients are for waves along +x (direction 0). `body.wave_force` (one regular wave) is unchanged, now including the component's phase, and is rejected together with several components. The body added mass A and damping B stay constant (as in the Paredes model, taken at one frequency): a broad-band spectrum needs A(omega), B(omega) with retardation functions, which is not implemented (see limitations).
* Waves: `waves.spectrum` {type, hs_m, tp_s, gamma, n_components, seed, omega_min_factor, omega_max_factor, direction_deg}; `waves.components` may be given in addition or instead.

## Depth-varying currents
`environment.current_profile`:
* `{"type": "table", "points": [{"depth_m": 0, "velocity_m_s": [u, v, w]}, ...]}`: depth below the still-water level, strictly increasing, linear interpolation, surface value above the surface, last value below the deepest point (give an explicit zero point at the seabed for a boundary layer);
* `{"type": "power_law", "surface_velocity_m_s": [...], "exponent": 0.143, "water_depth_m": h}`: U(d) = Us ((h - d)/h)^alpha, zero at and below the seabed.
`environment.current_m_s` (uniform) is unchanged. The current is superposed on the wave orbital velocity (`withCurrent`), for both the lines and the point elements; the accelerations are those of the waves only. It is available in the single-line runner and in the platform runner (with or without waves; with waves the lines' environment is rebuilt from the wave field plus the current).

## Verification (`tests/test_irregular.cpp`)
| Check | Result |
|---|---|
| PM density vs the closed form (5 frequencies); m0 of PM/JONSWAP (gamma 1, 3.3, 7) vs Hs^2/16 | 1e-4; 1e-5 |
| Discretised Hs; same seed reproducible, other seed different; one component per bin | exact (1e-12) |
| Variance of a 20 000 s elevation record (80 components, Hs 3 m) vs Hs^2/16 | 0.5618 m^2 vs 0.5625 m^2 (Hs 2.998 m) |
| Velocity / acceleration of the spectrum field (unstretched) vs the sum of the single-component fields | 1e-12 |
| Current profile: table interpolation, clamping, power law, composition with waves; JSON helper | 1e-14 |
| Vertical straight line in a sheared current: Morison drag at every node and the total vs the analytic sum | exact (0 difference), total 146.22 N |
| Free body in a 30-component JONSWAP sea with a 4-row force table: heave vs the sum of the steady-state responses of the components | max difference 2.0e-6 m for an amplitude of 2.2e-2 m (1e-4) |
| Component outside the force table / single-wave input with several components | rejected |
| Existing regular-wave tests and the Paredes wave validation after the kinematics optimisation | unchanged |

## Limitations
* Linear waves: no second-order wave loads (drift, slowly varying) on the body, no wave-current interaction (no Doppler shift of the waves, the current is simply added), no wave spreading (one direction per spectrum; several spectra can be built from `components`), no wave-induced current.
* The body load coefficients must be supplied as a table; A and B are constant. The single set of Paredes coefficients (1.3, 1.4, 1.5 s) covers only a narrow-band sea.
* The current acts on the lines and point elements only, not on the body (no current drag on the hull).
* The random-frequency discretisation does not repeat within any simulation time, but a record of finite length has the sampling variability of a random process (variance within a few percent for 80 components over 20 000 s; much more for a 3-hour storm statistic of extremes, not studied here).
