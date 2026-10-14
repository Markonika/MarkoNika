# Fully implicit integration (generalised-alpha / Newmark)

`numerics.scheme = "genalpha"` (or `"newmark"`, rho_inf = 1) integrates the whole line implicitly (`Scheme::GenAlpha`, `LumpedMassCable::genAlphaStep`). The explicit RK4 / Verlet schemes remain the default.

## Method
Chung & Hulbert (1993) generalised-alpha with the Newmark updates, applied to the acceleration form of the nodal equations, A(r, v, t) = M(r)^-1 F(r, v, t) (added mass, drag, seabed, point elements, rope law included):

    (1 - alpha_m) a1 + alpha_m a0 = A( (1 - alpha_f) r1 + alpha_f r0, (1 - alpha_f) v1 + alpha_f v0, t + (1 - alpha_f) dt )

with alpha_m, alpha_f, gamma, beta from rho_inf (second order, unconditionally stable for the linear problem; rho_inf = 1 is the non-dissipative average-acceleration Newmark; rho_inf < 1 damps the highest frequencies). The unknowns are the interior node positions, found by Newton iteration:
* Jacobian J = (1 - alpha_m)/(beta dt^2) I - (1 - alpha_f)(dA/dr + gamma/(beta dt) dA/dv): block tridiagonal (a node depends on itself and its two neighbours). dA/dr and dA/dv are assembled by finite differences with three node colours (18 evaluations of A, 12 in planar mode) and reused for `jacobian_interval` steps (default 10); the blocks do not depend on dt. They are rebuilt when a damped Newton step signals a poor Jacobian or convergence was slow, and once more within a step before it is declared failed. Block-Thomas solve with 3x3 blocks.
* Backtracking line search on the residual norm (the bilinear tension makes the plain Newton step cycle across the slack/taut kink). A step that cannot find a descent direction is accepted only if the residual has dropped by 1e-3 (noise floor of the finite-difference Jacobian), otherwise it counts as failed.
* Convergence: the accepted correction < `newton_tol` (default 1e-8) x the smallest segment length.
* Step control: the step is `dt_s` or `implicit_dt_factor` (default 10) x the axial-wave CFL step. A failed step halves a sub-step factor (down to 1/256, where the last iterate is accepted and counted in `newton_failures`); every three successes it grows by 1.5 back to 1. So the step shrinks only around slack/snap events (`step_halvings` in the output).
* The internal damping c_int is treated implicitly through the velocity; the c_int and soil step limits disappear. Maxwell internal strains are held during the step and advanced exactly afterwards (as in the explicit schemes).

## Verification (`tests/test_genalpha.cpp`)
| Check | Result |
|---|---|
| Vibrating string, Newmark vs RK4 reference, dt = 0.04 / 0.02 / 0.01 | error ratios 4.04, 4.01 (second order) |
| rho_inf = 0.8, dt = 0.02 / 0.01 | ratio 4.0 (second order) |
| Damped longitudinal mode (beta = 2e-4), dt = 16 x the explicit limit | energy decay within 2 % of the exact exp(-beta omega^2 t) (error from the period error omega dt = 0.4) |
| Strongly damped mode (beta = 5e-3, explicit limit 40 x below the CFL step) | stable, finite, monotonically dissipative |
| Lopez-Olocco WO_CW, A = 0.225 m, T = 2.8 s (seabed, Morison, c_int) | 24.302 N vs 24.320 N (RK4), -0.07 %; the same with the Jacobian rebuilt every step |

## Efficiency and accuracy on the validation cases
Lopez-Olocco WO_CW, A = 0.225 m, T = 2.8 s (explicit RK4 with implicit damping: dt = 1.0e-4 s, 24.320 N, about 10 s; explicit damping: 1.4e-5 s, about 50 s). Wall time includes about 0.7 s of static relaxation. `implicit_dt_factor` f:

| f | dt [s] | step halvings | max tension [N] | difference | wall [s] |
|---|---|---|---|---|---|
| 10 | 1.0e-3 | 0 | 24.2998 | -0.08 % | 2.4 |
| 20 | 2.0e-3 | 0 | 24.2366 | -0.34 % | 1.7 |
| 30 | 3.0e-3 | 641 | 24.4196 | +0.41 % | 2.7 |
| 40 | 4.0e-3 | 418 | 24.3053 | -0.06 % | 1.8 |
| 60 | 6.0e-3 | 1728 | 24.9415 | +2.6 % | 3.5 |

The error is not monotone in the step once step halvings occur; f = 10-20 is the range checked. About 4-5 times faster than RK4 with implicit damping and about 20 times faster than the original explicit damping at errors below 0.4 %.

**The snap case does not benefit.** Azcona configuration 2, T = 1.58 s (60 segments, maximum fairlead tension, measured about 44.5 N), explicit RK4 gave 57.4 N:

| scheme | dt [s] | step halvings | max tension [N] | wall [s] |
|---|---|---|---|---|
| GenAlpha rho 0.8, f = 1 | 7.9e-5 | 0 | 60.7 | 15 |
| f = 2 | 1.6e-4 | 27 | 53.7 | 10 |
| f = 5 | 3.9e-4 | 1723 | 51.2 | 9 |
| f = 10 | 7.9e-4 | 2509 | 50.4 | 10 |
| Newmark (rho 1), f = 2 | 1.6e-4 | 11199 | 56.3 | 35 |

The implicit scheme is not faster here (the sub-stepping drops to the explicit step around the snap) and the maximum still scatters between 50 and 61 N, as with the explicit schemes and the discretisation study (48-67 N). So the time integrator is not the cause of the snap-case spread; the sensitivity lies in the model and its discretisation (slack/snap, seabed contact, anchor position).

## Limitations
* The Jacobian is finite-difference: one rebuild costs 12-18 evaluations of the acceleration function, so the scheme only pays when steps are much larger than the explicit limit and the response is smooth (taut/slowly slack lines).
* Not a step-size error controller: the step is fixed apart from the Newton-failure halving; no local error estimate.
* Large steps skip fast physics (axial waves, contact chatter); high-frequency tension content (snap peaks) is smoothed with rho_inf < 1.
* `implicit_damping` (Strang splitting for the explicit schemes, assumptions item 45) is separate and not used with GenAlpha.
