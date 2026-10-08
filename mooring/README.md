# Mooring line dynamics solver (lumped-mass, C++17)

Copyright Marko Nika. All rights reserved. Proprietary and confidential. No licence is granted;
the GPL `LICENSE` at the repository root belongs to the profile repository, not to this directory.

## Build and test
    cmake -S . -B build && cmake --build build -j && ./build/mooring_tests

Warnings: `-Wall -Wextra -Wpedantic`. Dependencies: see `THIRD_PARTY.md`.

## Status
Milestones 1 (static solver + catenary) and 2 (dynamic solver) done.

| Test (milestone 1) | Result |
|---|---|
| Elastic catenary self-consistency (Eqs. 3.54-3.56 form) | pass |
| Unequal end heights | pass |
| Static LM cable vs catenary, EA = 200 kN, N = 10..160 | rel. L2 3.6e-4 -> 1.5e-6, observed order 1.99 |
| Same, EA = K*L (diagnostic) | rel. L2 2.6e-4 -> 1.1e-6, order 1.99 |
| Vibrating string, t = 2 s, dt = 1e-4, N = 5..80 | rel. L2 8.3e-2 -> 2.9e-4, order 2.01 (RK4 and Verlet) |
| Energy, 200 s undamped kicked catenary, N = 40, CFL 0.25 | max drift 5.4e-5 (RK4), 1.0e-4 (Verlet) of excitation energy; at CFL 0.5: 5.1e-4 / 4.2e-4 |
| Support reactions, static line | vertical sum = weight (1e-6), horizontal sum = 0 |

## Equation-to-code map
| Equation (Paredes 2016) | Code |
|---|---|
| 3.25 bilinear tension (+ internal damping) | `LumpedMassCable::computeForces` |
| 3.57 / 3.58 string wave | `tests/test_dynamics.cpp` (exact solution) |
| 3.54-3.56 elastic catenary | `ElasticCatenary` (`src/catenary.cpp`) |

## Limitations (so far)
Perfectly flexible cable (no bending/torsion/VIV); no hydrodynamics or seabed yet (milestone 3); momentum conservation not tested (ends are held).
See `docs/assumptions.md`.
