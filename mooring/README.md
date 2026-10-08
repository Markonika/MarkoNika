# Mooring line dynamics solver (lumped-mass, C++17)

Copyright Marko Nika. All rights reserved. Proprietary and confidential. No licence is granted;
the GPL `LICENSE` at the repository root belongs to the profile repository, not to this directory.

## Build and test
    cmake -S . -B build && cmake --build build -j && ./build/mooring_tests

Warnings: `-Wall -Wextra -Wpedantic`. Dependencies: see `THIRD_PARTY.md`.

## Status
Milestone 1 (static solver + catenary test): done.

| Test (milestone 1) | Result |
|---|---|
| Elastic catenary self-consistency (Eqs. 3.54-3.56 form) | pass |
| Unequal end heights | pass |
| Static LM cable vs catenary, EA = 200 kN, N = 10..160 | rel. L2 3.6e-4 -> 1.5e-6, observed order 1.99 |
| Same, EA = K*L (diagnostic) | rel. L2 2.6e-4 -> 1.1e-6, order 1.99 |

## Equation-to-code map
| Equation (Paredes 2016) | Code |
|---|---|
| 3.25 bilinear tension | `LumpedMassCable::computeForces` |
| 3.54-3.56 elastic catenary | `ElasticCatenary` (`src/catenary.cpp`) |

## Limitations (so far)
Perfectly flexible cable (no bending/torsion/VIV); no dynamics, hydrodynamics or seabed yet.
See `docs/assumptions.md`.
