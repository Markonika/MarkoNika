# Implicit body integration with iterated coupling

`numerics.coupling = "implicit"` (platform runner; `CoupledSystem::setImplicitCoupling`) replaces the explicit partitioned scheme (velocity-Verlet for the body, lines advanced to the new time with the fairlead moving linearly from the old to the new pose) by an implicit one:

* The body is integrated by Newmark-beta (`newmark_gamma` = 0.5, average acceleration, no numerical damping; gamma > 0.5 damps high frequencies) with the unknown pose xi_{n+1}: R(xi) = M a + B v + Dq |v| v - F_ext(xi, lines, t_{n+1}) = 0.
* Newton iterations with the approximate Jacobian M/(beta dt^2) + gamma/(beta dt) (B + 2 Dq |v|) + K, where K is the static linearised stiffness incl. hydrostatics (`CoupledSystem::stiffness`, built at the first step and rebuilt when an iteration needs more than half of `coupling_max_iter`, once more on failure).
* Every iteration restores all lines to the state at the start of the step (`CableModel::saveState` / `restoreState`) and re-advances them with the trial fairlead path, so the line forces are the genuine dynamic ones; the lines are advanced on `numerics.threads` threads. Convergence: Newton correction < `coupling_tol` (1e-9) x (1 + |xi|). The report gives line re-advances, stiffness builds and unconverged steps.
* The line integrator is chosen independently (explicit RK4 or the implicit generalised-alpha with adaptive steps).

## Why (and when)
The explicit scheme is stable only for dt_body < 2/omega_s (omega_s = sqrt(K_line/M)): stiff lines on light bodies (tendons, taut synthetic legs) force steps far below the physical period. The Paredes buoy cases are not in this regime (omega_s about 33 rad/s, dt_body 2 ms is chosen for accuracy), so **this option does not speed up the existing validation cases**: it costs several line re-advances per body step. It pays only when the body step can be raised by more than the iteration count and the line integration scales with the body step (implicit lines).

## Verification (`tests/test_coupling.cpp`)
* Stiff line (EA = 2e6 N, 10 m, body 2 kg: omega_s = 316 rad/s, explicit limit 6.3 ms), slow sinusoidal heave force: implicit coupling at dt_body = 50 ms (8 x the explicit limit) gives the heave response -2.478e-5 m against -2.503e-5 m from the explicit scheme at 2 ms (1 %; 30 steps, 114 line advances, 3.8 per step) while the explicit scheme at 50 ms diverges.
* Pendulum swing on a moderately stiff line (EA = 2e4 N), x(3 s) against the explicit scheme at 5 ms: differences 1.2 % of the amplitude at dt_body = 0.2 s and 0.1 s, 0.06 % at 0.025 s. The accuracy at large steps is limited by the line's own lowest transverse mode (period about 0.14 s), which the coarse body step does not resolve; it is not the Newmark phase error (second order, 0.3 % in frequency at dt = 0.2). With 40 instead of 10 segments the error at dt = 0.05 s is 0.1 % of the amplitude.

## Limitations
* Newton with an approximate (static, occasionally refreshed) stiffness converges linearly: 4-7 line advances per step in the tests.
* The line force is evaluated once per iteration from a fresh advance, so a line with many sub-steps makes each iteration expensive; use implicit line integration with a large step to benefit.
* Only models that implement `saveState` / `restoreState` (the lumped-mass cable) can be coupled implicitly; not a monolithic scheme (the line states are not part of the Newton unknowns).
