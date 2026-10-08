# Mooring line dynamics solver (lumped-mass, C++17)

Copyright Marko Nika. All rights reserved. Proprietary and confidential. No licence is granted;
the GPL `LICENSE` at the repository root belongs to the profile repository, not to this directory.

## Build and test
    cmake -S . -B build && cmake --build build -j && ./build/mooring_tests   # about 1 min; full Chalmers grid: ctest -R chalmers_grid

Warnings: `-Wall -Wextra -Wpedantic`. Dependencies: see `THIRD_PARTY.md`.

## Status
Milestones 1-8 done (static solver, dynamics, hydrodynamics + seabed, Chalmers validation, planar/3D regression, point elements, 6-DOF platform with multi-line coupling, Paredes free-buoy and moored-buoy validation).

| Test (milestone 1) | Result |
|---|---|
| Elastic catenary self-consistency (Eqs. 3.54-3.56 form) | pass |
| Unequal end heights | pass |
| Static LM cable vs catenary, EA = 200 kN, N = 10..160 | rel. L2 3.6e-4 -> 1.5e-6, observed order 1.99 |
| Same, EA = K*L (diagnostic) | rel. L2 2.6e-4 -> 1.1e-6, order 1.99 |
| Vibrating string, t = 2 s, dt = 1e-4, N = 5..80 | rel. L2 8.3e-2 -> 2.9e-4, order 2.01 (RK4 and Verlet) |
| Energy, 200 s undamped kicked catenary, N = 40, CFL 0.25 | max drift 5.4e-5 (RK4), 1.0e-4 (Verlet) of excitation energy; at CFL 0.5: 5.1e-4 / 4.2e-4 |
| Added-mass node solve, Morison drag, seabed force law (unit tests) | exact to round-off |
| String in water with Cm rho A1 (1+eps) = m_l | omega 9.8828 vs 9.8827 predicted; wet/dry ratio 0.7071 |
| Quadratic drag decay, 1/a = 1/a0 + kappa t | 0.1-0.7 % for t >= 1 s; 3.2 % at t = 0.5 s (first-order theory, bound 3.6 %) |
| Partly submerged catenary | support z sum = sum of node weights (1e-6) |
| Chain with seabed touchdown (a = H/w = 5 m) | H 21.387 vs 21.381 N, V 17.765 vs 17.760 N (N = 60); 37 nodes on bed |
| Dynamic drop onto seabed at true Ks | stable, settles to H = 21.373 vs 21.381 N, no penetration > 1e-7 m |
| Chalmers static top tension (22.68 N published) | 22.653 N (N = 66), 22.685 N (N = 132) |
| Chalmers 30-case grid vs Table 7 (clockwise motion) | regression r2 = 0.985 (published 0.98), RMSE 1.4 N, bias +0.12 N, 24/30 within 5 %, max 6.7 % |
| Same grid, opposite rotation sense (sensitivity) | r2 vs 1:1 = 0.81, bias +4.4 N (+9.8 %) - direction matters, see docs/chalmers_validation.md |
| Planar mode vs 3D, Chalmers (3 cases) | identical (difference 0.0); 3D run stays in plane (< 1e-12 m) |
| Rotation of the set-up by 37 deg about z | smooth cases 1e-4..2.4e-4, snap case 0.1-0.5 %: all inside the measured round-off noise floor of the metric (1e-12 change in dt gives 3e-5..1.8e-4) |
| Lateral current +/-0.3 m/s | tension identical, deflection mirrored (max \|y\| 1.160 m both); no-current run stays planar |
| 1 mm out-of-plane perturbation, 6 cycles | never exceeds 1 mm, 0.3-0.5 mm at the end; maxima change by <= 0.8 % |
| Circularly polarised string (3D analytic) | radial error 1.6e-8, phase error 1e-7 rad after 0.92 periods |
| Clump 300 N at mid-span vs exact elastic half-catenaries | rel. L2 6.0e-5 -> 9.6e-7 (N = 20..160), order 1.99; support loads sum to weight (< 3e-9) |
| Floater (2500 N buoyancy, arch) vs exact | rel. L2 2.0e-5 -> 3.1e-7, order 1.99 |
| Point drag in a 1.5 m/s current | support load 10.602875 N = Morison drag (exact) |
| Point added mass, 2-segment oscillator | period 3.31136 s = exact; Cm=1/Cm=0 ratio 1.10554 = exact |
| Rotation kinematics, free heave decay with A, B, C | orthonormal R; decay error 2.4e-5 / 6.0e-6 / 1.5e-6 for dt = 2 / 1 / 0.5 ms (order 2) |
| Body on a vertical line: heave frequency | 9.99230 vs exact end-mass spring 9.99167 rad/s (6e-5); ratio 1-50 differ by <= 1e-3 |
| Body on a line: pendulum period | 6.37436 vs 6.37495 s (-9e-5); independent of sub-step ratio 1-100 |
| JSON platform runner, same pendulum | period error -1.1e-4, sub-step ratio 51 |
| 3 legs at 120 deg, asymmetric: static equilibrium | Newton in 4 iterations, residual 9e-10 N / 3e-11 N m; independent force/moment balance to 1e-6 |
| Coupling stability | stable at 0.4 x, unstable at 2 x the estimate 2/omega_s (two body masses) |
| Paredes free buoy: heave / pitch damped period (potential-theory A, B) | 1.0997 s vs 1.112 +- 0.006 (-1.1 %); 1.0979 s vs 1.170 +- 0.005 (**-6.2 %, known disagreement**) |
| Paredes moored statics, CON1 / CON2 / CAT (no tuning) | leg tension 2.75 / 10.9-11.0 / 2.93 N vs 2.8-3.1 / 10.6-11.0 / 3.0-3.1 N; draft change +0.6 / +12.9 / +3.4 mm vs +1 / +13 / +4 mm |
| Paredes surge secant stiffness vs thesis Fig. 5.21 curves | CON2 within 3 %, CAT 7-10 % low (U-shape reproduced), **CON1 15-19 % low** (all below the 41 N/m design value, as measured) |
| Paredes damped surge periods CON1 / CON2 / CAT | 9.286 / 9.207 / 9.349 s vs 8.561 / 9.22 / 9.14 s (+8.5 % / -0.1 % / +2.3 %) |
| Support reactions, static line | vertical sum = weight (1e-6), horizontal sum = 0 |

Details: `docs/chalmers_validation.md` (results, convergence, caveats), `docs/config.md` (configuration fields).
Run the Chalmers case: `build/mooring_run examples/chalmers/chalmers_config.json`; full grid:
`ctest --test-dir build -R chalmers_grid` or `scripts/chalmers_grid.py OUTDIR` (plots need matplotlib; analysis only).

Paredes benchmark: `build/paredes_report examples/paredes` (details and caveats: `docs/paredes_validation.md`).
Platform: `build/mooring_platform examples/platform/three_leg_example.json` (config fields in `docs/config.md`).

## Equation-to-code map
| Equation (Paredes 2016) | Code |
|---|---|
| 3.25 bilinear tension (+ internal damping) | `LumpedMassCable::computeForces` |
| 3.26 submerged weight | `CableParams::submergedWeight` |
| 3.27 added mass | `addedMassSolve` (`hydro.hpp`), `acceleration` |
| 3.28-3.29 drag | `morisonDrag` |
| 3.32-3.36 seabed | `seabedForce` |
| 3.59 body equation (constant A, B, C) | `RigidBody6DOF`, `CoupledSystem::step` |
| A.4-A.13 static equilibrium of CON1/CON2 | solved by the full line statics + `CoupledSystem::solveEquilibrium` |
| 3.57 / 3.58 string wave | `tests/test_dynamics.cpp` (exact solution) |
| 3.54-3.56 elastic catenary | `ElasticCatenary` (`src/catenary.cpp`) |

## Limitations (so far)
Perfectly flexible cable (no bending/torsion/VIV); small-angle rigid-body dynamics (no quaternion/Euler nonlinearity), constant A and B, still water only (no wave kinematics until milestone 9), flat seabed, no soil dynamics beyond the spring-damper; momentum conservation not tested (ends are held).
See `docs/assumptions.md`.
