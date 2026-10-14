# Bending stiffness

## Model
Discrete curvature energy on the line, for unit segment tangents t_j = (r_{j+1} - r_j)/|r_{j+1} - r_j|:

    U_b = sum over interior nodes i of  EI_i |t_i - t_{i-1}|^2 / (2 lbar_i)

with lbar_i the mean *unstretched* length of the two segments and EI_i the mean of the two segment values. For small angles |t_i - t_{i-1}| = theta_i, so U_b = sum of 1/2 EI kappa^2 lbar with kappa = theta/lbar. The nodal forces are minus the gradient (`LumpedMassCable::addBendingForces`): with g = EI (t_i - t_{i-1})/lbar and P = I - t t^T, dU/dd_i = P_i g/|d_i| and dU/dd_{i-1} = -P_{i-1} g/|d_{i-1}|. Each node therefore depends on the two nodes on either side (**five-node stencil**).

* Ends: by default the end segments are pinned (no bending term at the anchor/fairlead nodes). A **clamped** end (built-in end, bend stiffener, hang-off angle) is selected by giving its direction c; it adds EI (1 - (c.t)^2)/l0 for the end segment, i.e. the support node sees the mirror image of the first segment (curvature 2 phi / l0 over a half cell).
* Bending acts at any tension, including slack lines, which is its purpose. The force on the fairlead (`endForce`) includes the shear from bending; `endSegmentTension` remains the axial tension.
* Non-uniform lines: EI per section.
* Config: `line.bending_stiffness_Nm2` (also in `line.sections[]` and `line_defaults`), `line.bending_anchor_direction`, `line.bending_fairlead_direction` ([x, y, z]).

## Numerics
* Explicit schemes: the highest bending frequency of the five-point stencil is omega = 4 sqrt(EI/(m_l l0^4)), so dt <= cfl * 2/omega = cfl l0^2 sqrt(m_l/EI)/2 is added to the step limit (`stableDt`).
* Static relaxation: the fictitious mass includes the bending stiffness bound 16 EI/l0^3.
* Implicit scheme: the finite-difference Jacobian is assembled with five colours and the Newton system is solved by a banded LU (half bandwidth 11 in dof units, no pivoting); without bending the tridiagonal block solver is used, so lines without bending are unchanged bit for bit. The linearised modes (`naturalModes`) use the same blocks.

## Verification (`tests/test_bending.cpp`)
| Check | Result |
|---|---|
| Force = -gradient of the energy (random 3-D shape, both ends clamped) | max error 4.5e-7 N for components up to 27 N |
| Pinned beam under uniform load (N = 40): mid-span deflection vs 5 q L^4/(384 EI) | +0.01 % |
| Clamped-clamped beam: vs q L^4/(384 EI) | +0.49 % |
| Taut pinned beam, natural frequencies 1-3 vs sqrt((T k^2 + EI k^4)/rho), k = j pi/L (N = 100, EI = 50, T = 100 N) | 10.1710 / 21.7251 / 35.7751 rad/s, -0.002 / -0.011 / -0.033 % |
| Generalised-alpha (banded Jacobian) vs RK4 on a vibrating beam | 1e-7 of a 5e-4 m amplitude; the explicit step is the bending limit |

Static deflection tests need a stiff axial spring and a small load: with a soft EA the beam is axially soft (the end shear has an axial component and the tension is non-uniform; at EA = 0.5 N the pinned beam deflected 5 % more than beam theory), and with a large load the membrane tension stiffens it. That is physics of the fixed-end beam-column, not a discretisation error (the deflection did not change with N).

## Sensitivity: does bending explain the Azcona snap-case excess?
Azcona configuration 2, T = 1.58 s (60 segments, adaptive implicit integration, tolerance 1e-3, measured maximum about 44.5 N, converged value without bending 48.6 N). The chain of the experiment has no defined bending stiffness (a solid steel rod of 3.4 mm diameter would have EI = 1.4 N m^2, a studless chain much less), so EI is varied as an exploration, not a calibration:

| EI [N m^2] | 0 | 1e-4 | 1e-3 | 1e-2 | 1e-1 | 1 |
|---|---|---|---|---|---|---|
| max tension [N] | 48.6 | 48.4 | 49.7 | 58.5 | 63.0 | 81.4 |

Static top tension is unchanged up to EI = 1e-1 (15.019 -> 15.016 N). Bending has no effect for EI <= 1e-3 and *raises* the snap peak above that, so it cannot explain the +12 % excess over the measurement; the candidates remain the anchor position, friction, drag coefficients and internal damping (docs/azcona_validation.md).

## Limitations
* Small-strain Euler-Bernoulli bending of a straight reference configuration: no torsion, no moment-curvature nonlinearity, no coupling with tension (stiffening is captured through the geometry only), no bending damping.
* The bending term uses the unstretched lbar; the O(strain) difference to the stretched length is neglected.
* Pinned/clamped end conditions only (no rotational spring at the supports); the clamp direction is fixed in space.
* The banded LU has no pivoting: a Newton failure (singular pivot) triggers the usual step reduction.
