# Mooring line dynamics solver (lumped-mass, C++17)

Copyright Marko Nika. All rights reserved. Proprietary and confidential. No licence is granted;
the GPL `LICENSE` at the repository root belongs to the profile repository, not to this directory.

## Build and test
    cmake -S . -B build && cmake --build build -j && ./build/mooring_tests

Warnings: `-Wall -Wextra -Wpedantic`. Dependencies: see `THIRD_PARTY.md`.

## Status
Milestones 1 (static solver + catenary), 2 (dynamic solver) and 3 (hydrodynamics + seabed) done.

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
| Support reactions, static line | vertical sum = weight (1e-6), horizontal sum = 0 |

## Equation-to-code map
| Equation (Paredes 2016) | Code |
|---|---|
| 3.25 bilinear tension (+ internal damping) | `LumpedMassCable::computeForces` |
| 3.26 submerged weight | `CableParams::submergedWeight` |
| 3.27 added mass | `addedMassSolve` (`hydro.hpp`), `acceleration` |
| 3.28-3.29 drag | `morisonDrag` |
| 3.32-3.36 seabed | `seabedForce` |
| 3.57 / 3.58 string wave | `tests/test_dynamics.cpp` (exact solution) |
| 3.54-3.56 elastic catenary | `ElasticCatenary` (`src/catenary.cpp`) |

## Limitations (so far)
Perfectly flexible cable (no bending/torsion/VIV); still water only (no wave kinematics until milestone 9), flat seabed, no soil dynamics beyond the spring-damper; momentum conservation not tested (ends are held).
See `docs/assumptions.md`.
