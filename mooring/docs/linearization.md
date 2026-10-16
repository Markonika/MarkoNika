# Linearisation about the static state (`mooring_modes`, `symmetricEigen`, `CoupledSystem::stiffness`, `LumpedMassCable::naturalModes`)

## What is computed
* **Line natural modes** (`LumpedMassCable::naturalModes`): K phi = omega^2 M phi with K = -M dA/dr from the finite-difference Jacobian blocks of the implicit solver (A = M^-1 F, block tridiagonal), M the nodal inertia (mass + added mass normal to the line, point elements), anchor and fairlead fixed. K is symmetrised (the asymmetry is reported; 1e-7 .. 1e-3 in the cases below); the eigenproblem is solved by a cyclic Jacobi method after scaling by M^-1/2. The stiffness is the *dynamic* (instantaneous) one: the Maxwell branches of the rope law are included. Modal damping ratio zeta_k = phi_k^T C phi_k / (2 omega_k) with C = -M dA/dv (internal damping; quadratic drag has no linearisation at rest, so the hydrodynamic damping is not included).
* **Body stiffness** (`CoupledSystem::stiffness`, `lineStiffness`): K = -dG/dxi by central differences (step 1e-4 m/rad) of the static line solutions (the lines are solved statically at every perturbed pose; the evaluations run on `numerics.threads` threads); this is the same Jacobian that the equilibrium Newton uses, so it is consistent with the model's definition of the line load, including the geometric stiffness T/L and the rotation coupling.
* **Body natural frequencies** (`bodyModes`): (C + K_lines) x = omega^2 (M + A) x, undamped, by Cholesky + Jacobi.
* App: `build/mooring_modes config.json [key=value ...] [--lines N]` for a single-line or a platform config.

## Verification (`tests/test_linearize.cpp`)
| Check | Result |
|---|---|
| Jacobi vs the 1-D Laplacian eigenvalues (n = 30), orthonormality, residual | 1e-12 |
| Taut chain (N = 40, planar): transverse modes 1-3 vs (2/l0) sqrt(T/(m_l(1+eps))) sin(k pi/2N) | 9.92707 / 19.83884 / 29.72001 rad/s, equal to 5 digits |
| Same chain: first longitudinal mode vs (2/l0) sqrt(EA/m_l) sin(pi/2N); modal damping vs beta omega/2 | 314.0785 (exact to 6 digits); 0.01570 vs 0.01570 |
| Two-section string, lowest frequency vs the exact two-medium eigenvalue equation | 8.9498 rad/s, -0.0007 % |
| Body hanging on a line: K(heave) vs EA/L0; K(surge) vs T/L; both natural frequencies found | 200.0000 vs 200.0000; 1.948 vs 1.943 (the 0.24 % is the line's own weight) |

## Paredes buoy (no tuning; configs `examples/paredes/*.json`, constant A at the measured frequency of each configuration)
Linearised natural periods at the static equilibrium, against the damped periods of the time-domain surge decay (release from 0.1 m) and the measured ones (Table 5.11):

| | linear surge period (undamped) | decay test, model (damped) | measured |
|---|---|---|---|
| CON1 | 9.340 s | 9.286 s | 8.561 +- 0.006 s |
| CON2 | 9.156 s | 9.207 s | 9.22 +- 0.08 s |
| CAT | 10.018 s | 9.349 s | 9.14 +- 0.01 s |

Heave 1.097 / 1.095 / 1.096 s and pitch 1.074 / 1.094 / 1.082 s for CON1 / CON2 / CAT (undamped; the free-buoy decay gave heave 1.0997 s damped, pitch 1.0979 s). CON1 and CON2 are consistent with the decay periods (the undamped period is slightly below the damped one for CON2, +0.6 %). For CAT the linear period at zero offset is 7 % above the decay period: the catenary stiffness is amplitude dependent (U-shaped, docs/paredes_validation.md: 33.9 N/m at -0.1 m, 20.4 N/m at +0.1 m), so the 0.1 m decay averages a stiffer line than the stiffness at rest. The mooring line modes are 1.5 s (leg 1, CON2) and above; the lowest line mode periods (0.7-1.5 s) are close to the wave periods of the Paredes RAO tests (1.3-1.5 s), which is worth noting for the wave-run comparison.

Lopez-Olocco chain at its static state (60 segments, planar not set, 177 dof): lowest modes 1.25, 2.46, 3.67, 4.37, 4.89, 6.11 rad/s (periods 5.03, 2.56, 1.71, 1.44, 1.28, 1.03 s), i.e. the excitation periods 2.8-5.5 s of the tests lie above the first two modes.

## Fairlead dynamic stiffness (`LumpedMassCable::fairleadDynamicStiffness`, `CoupledSystem::dynamicLineStiffness`, `mooring_modes --dynstiff w1,w2,...`)
The line seen from its fairlead at angular frequency omega, about the static shape, anchor fixed, interior nodes eliminated:
`f = f0 - D(omega) u`, `D = -(G_tt - G_ti (G_ii + omega^2 M)^-1 G_it)`, `G = dF/dr + i omega dF/dv` (central differences of the nodal forces of `computeForces`: stretch-rate/internal damping
included; quadratic Morison drag has no linearisation at rest in still water; wave forces on the line are not included), `M` the nodal inertia of `naturalModes` (mass + added mass). For a platform the legs are assembled
into a complex 6 x 6 body stiffness `sum_k J_k^T D_k J_k + K_geo` (`J = [I, -[R a]x]`, `K_geo` the geometric stiffness of the static fairlead force through the rotating lever arm).

Verification (`tests/test_linearize.cpp`, all in the default suite; plus checks on the Paredes legs):
| Check | Result |
|---|---|
| taut beaded chain (N = 60), transverse D at 0.3 / 0.6 / 0.9 / 0.97 / 1.03 times the first natural frequency vs the exact lumped result `k (1 - sin((N-1) th)/sin(N th))`, `2 - 2 cos th = m w^2/k` | 6.9143 / -5.8235 / -86.3203 / -322.3061 / 341.2961, equal to the exact result to all digits shown (test tolerance 1e-3) |
| same chain vs the continuous string `T k cot(k L)` | within the O(1/N) lumping error (the fairlead node's half mass belongs to the body: 6.914 vs 6.841 at 0.3 w1) |
| static limit | transverse T/L to 1e-4, axial EA/L to 2e-3 |
| axial bead chain (link stiffness EA/l0) at 0.3 and 0.7 of the axial reference frequency | within 2e-3 |
| pole at the first natural frequency | sign change across 0.98 / 1.02 w1, magnitude > 10 T/L |
| Paredes CON2 platform, omega -> 0 | equals the finite-difference `lineStiffness` to 4 digits (26.6267 / 4.1106 / 17.9203 / 1.1760 ...), geometric term included |
| Paredes CON2 leg, omega at the `naturalModes` frequencies 4.1896 and 4.7207 rad/s | stiffness changes sign across each and reaches 6e4 N/m (surge, 4.1896) and 2e4 N/m (heave, 4.7207) there: poles at the natural frequencies |

## Limitations
* Linearisation at one static state; no amplitude dependence, no hydrodynamic (drag) damping (the dynamic stiffness above is undamped apart from the internal damping, which is why it has sharp poles), no frequency-dependent added mass. The body model has constant A and B anyway.
* Jacobi eigen-decomposition is O(n^3) per sweep: fine for lines up to about 240 segments (n about 700); not meant for much finer models.
* Finite-difference stiffness: a line near the slack/taut kink gives a stiffness that depends on the difference step; use it at taut, well-defined states.
