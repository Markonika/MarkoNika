# Assumptions and open items

1. **K = EA.** Paredes (2016) Eqs. 3.54-3.55 contain `m_l g s / K`, which equals `H s / EA`, so the
   "stiffness K = 200 kN/m" of the sec. 3.5.1 static test is the axial stiffness EA = 200 kN.
   The test additionally runs EA = K*L (stiff) as a diagnostic.
2. Static test is in air: weight per length w = m_l * g (dry).
3. Static relaxation uses fictitious nodal mass (massFactor * EA / l0, pseudo-dt = 1); only the
   converged equilibrium is physical. massFactor = 2 was unstable, 32 is the default.
4. Reported L2 error is relative, over nodes at Lagrangian coordinate s_i = i*l0.
5. Still needed: fairlead coordinates and line data for the Paredes buoy (thesis Figs 5.4-5.5,
   Tables 5.3, 5.5-5.10) - to be read from the supplied thesis PDF at milestone 7.
6. Chalmers inputs received and stored in `examples/chalmers/`.
7. **Vibrating string.** Paredes Eq. 3.57 gives c^2 = tau/(m_l(1+eps)) in Lagrangian s. For the stated
   case (L0 = 0.5 m stretched to 1 m, eps = 1, tau = 1 N, m_l = 1 kg/m, so EA = 1 N) c = 0.7071 m/s and
   the fundamental period is 1.414 s, not 2 s; the thesis' "2 s" is the simulated duration. Eq. 3.58 as
   printed (omega = c pi / (L(1+eps)) with L the stretched length) does not reproduce this, so the
   consistent Lagrangian form omega = pi c / L0 is used. Amplitude 1e-4 m keeps nonlinearity negligible.
8. **Internal damping / slack regularisation.** T = EA eps + c_int d(eps)/dt for eps > 0, clipped at T >= 0
   (a cable never pushes). c_int = 0 by default. Slack segments (eps <= 0, T = 0) and damping-induced
   clipping are counted in `DynStats` (`slackSegmentEvals`, `clippedTensionEvals`) and must be reported by
   every run. Explicit damping adds a stability limit, included in `stableDt()`:
   dt <= cfl * m_l l0^2 / (2 c_int).
9. **Time step.** dt = cfl * l0 / c with c = sqrt(EA/m_l), default cfl = 0.5 (both schemes stable).
   Velocity-Verlet treats c_int damping with the half-step velocity (lagged); it is only first-order
   accurate in the damping term.
10. **Force on body** (`forceOnBody`, `endForce`) is the net force on the end node (tension + half-node
    weight; hydrodynamics from milestone 3), not the end-segment tension alone: the latter misses
    w*l0/2 vertically (2.5 % for N = 40). The end node's own inertia is not passed to the body.
    Coupling is explicit/partitioned with linear interpolation of the fairlead motion; stability limits
    to be established at milestone 7.
11. **Momentum check (test 3).** The ends are always held (fixed or prescribed), so linear momentum is not
    conserved; only energy is checked. A free-end momentum test is deferred to the platform milestone.

12. **Hydrodynamics (milestone 3).** Forces follow Eqs. 3.27-3.31 per node. Water velocity in Eq. 3.30/3.31 is
    read as v_rel = v_w - dr/dt (the typeset "a_w" is a misprint). The tangential drag term is the
    sign-preserving 0.5 Cdt rho D (v.t)|v.t| t (the typeset absolute value is lost). The (1+eps) factor is
    realised by using the actual (stretched) tributary length of each node. Added mass uses
    c_i = Cm rho_w A1 L_stretched with A1 = pi/4 D0^2 unless given; solved exactly per node as
    M^-1 = t t^T/m + (I - t t^T)/(m + c), the water-acceleration term entering the right-hand side.
    Node tangent = central difference of the neighbours (one-sided at the ends).
13. **Surface.** Nodes blend linearly between dry weight m_l g (above) and submerged weight gamma_l (Eq. 3.26,
    `CableParams::submergedWeight`) over one segment length around `surfaceZ`; hydrodynamic forces are
    scaled by the same submerged fraction. Wave-surface kinematics (Wheeler) come in milestone 9.
14. **Seabed.** Flat, normal +z. Ks is the thesis' "GPa/m" quantity: normal force per length = Ks D1 dH
    (Pa/m * m * m = N/m); damping 2 zeta sqrt(Ks D1 m_l) while penetrating downward only (Eq. 3.32). Friction is
    -gamma_l mu min(|v_st|/v_lim, 1) t_s (Eq. 3.34: scaled by submerged weight, not by the normal force) and acts
    while in contact. Soil forces use the unstretched tributary length (no (1+eps), as in the equations).
15. **Time step with a seabed.** The contact oscillator has omega^2 = Ks D1/m_l independent of l0 (9.0e3 rad/s for the
    Chalmers chain), so `stableDt()` adds dt <= cfl*2/(omega (zeta + sqrt(zeta^2-1))) = 1.1e-4 s there. This
    limit is real; the drop test runs at dt = 2.8e-5 s because c_int = 50 N s adds a tighter one.
16. **Static relaxation with a seabed** caps the soil stiffness at EA/(l0^2 D1) (node contact stiffness = EA/l0) to keep
    the fictitious-mass scheme stable. Equilibrium penetration is then w l0^2/EA (7e-8 m in the tests) instead
    of w/(Ks D1); both are negligible. Dynamic runs use the true Ks.
17. **Chalmers chain diameters.** Table 3.1 gives D0 = D1 = link thickness = 0.0022 m; the JSON agrees. The nominal
    area A1 is not tabulated: pi/4 D0^2 is assumed (only matters for Cm != 0; Cm = 0 also fits the data).
18. Energy accounting (`energy()`) is valid only without hydro/seabed/variable weight (it includes no dissipation
    and assumes uniform w).
19. **Chalmers rotation direction** is not in the supplied inputs; clockwise in x-z (`direction = -1`) was read from
    thesis Fig. 3.5 (arrow on the motor circle) after a counter-clockwise first run over-predicted (not a blind choice). Results depend strongly on it (see docs/chalmers_validation.md).
20. **Per-cycle maxima** are taken from every time step; "mean maximum" averages cycles 5-9 of a 10-cycle run after a
    2-cycle cosine ramp of the radius. The thesis ran 15 cycles; maxima were stationary to <0.1 % after the ramp.
21. Top/anchor tension in outputs is the raw end-segment tension (EA eps + c_int d(eps)/dt, clipped at 0), not the
    support force (which also contains the end node's half-weight and hydrodynamic force).
22. **Planar mode** (`CableParams::planar`, `numerics.planar`): y-forces are set to zero and y is flattened to the anchor's y in the
    constructor and `setInitialState`; a top-motion function must stay in the plane. It is a consistency/cost option: 3D runs
    with in-plane data give identical results (bit-for-bit in the tests), so there is no separate 2D code path.
23. **Rotation invariance** holds to 1e-4 .. 2.4e-4 (smooth cases) and 0.1-0.5 % (snap case) in the maximum tension. **This is the noise floor of the
    metric, not a rotation error:** changing the time step by 1 part in 1e12 or the circle phase by 1e-9 deg moves the same maxima by 3e-5 .. 1.8e-4
    (smooth, 3 averaged cycles). The slack/contact events amplify round-off (sensitive dependence), so any maximum from a short window
    carries ~1e-4 relative noise (more for snap cases). An earlier note here blamed the static-relaxation tolerance; a direct test refuted
    that (tightening the tolerance reduced the static difference from 1.6e-6 to 2.8e-8 N but not the dynamic one). `maxOutOfPlane` is sampled every 64
    steps, so it can miss short excursions.
24. The out-of-plane stability check is empirical (three Chalmers cases, 1 mm perturbation, 6 cycles). It shows no growth
    there; it is **not** a proof that no parametric out-of-plane instability exists elsewhere in the parameter space.
25. **Point elements** (`PointElement`, `LumpedMassCable::addPointElement`): weight m g always; buoyancy rho_w g V, drag
    0.5 rho_w Cd A |v_rel| v_rel and added mass Cm rho_w V (a_w - a) only on the submerged fraction of the node and only with
    `Environment::hydro`. Added mass is isotropic (point body) and added to the node inertia; the Froude-Krylov term is
    omitted, consistently with Eq. 3.27 for the cable. `floater()`/`clump()` derive V from the thesis' table values
    (buoyancy force / submerged weight and mass; Tables 5.7-5.8). Elements have no seabed contact of their own: a clump on
    the bed rests via the stiff contact of its node (penetration ~ W/(Ks D1 l0)). `energy()` ignores point elements.
26. Point-element drag acts in still water as well as in a current via `Environment::water`; there is no wave-surface
    interaction before milestone 9.
27. **Platform model (milestone 7).** `RigidBody6DOF`/`CoupledSystem`: (M + A) xi'' + B xi' + C xi + Dq|xi'|xi' = F0 + line loads,
    with constant (frequency-independent) 6x6 A, B, C, as in Paredes (2016) Eq. 3.59. xi = CG translation and a rotation vector
    measured from the free-floating reference pose; the rotational dynamics use omega ~ d(theta)/dt and a constant inertia tensor
    (small-angle model). Fairlead positions use the exact exponential map R(theta) (orthonormal), fairlead velocity is
    u' + theta' x (R a). **Not implemented yet:** full nonlinear rigid-body dynamics (quaternion + Euler equations), and any
    frequency dependence of A and B. The small-angle model is adequate up to roughly 15-20 deg, which is also the limit of the linear
    hydrostatic stiffness C. BEM coefficients are to be supplied as constants (file reader for Capytaine/NEMOH output: milestone 8/9).
28. **Mean loads.** Weight and buoyancy are not modelled separately: the reference pose is the free-floating equilibrium of the hull alone
    (weight = buoyancy there), C holds the hydrostatic restoring stiffness, and `F0` carries any other constant load. The lines
    then shift the equilibrium (pretension), found by `solveEquilibrium`. For a body hanging in air use F0 = -m g and C = 0.
29. **Coupling scheme** (explicit, partitioned, conventional serial staggered): per body step dt, (i) body kick-drift with the old
    acceleration (velocity-Verlet half step), (ii) every line is advanced to t+dt inside `forceOnBody` with its own CFL sub-steps and
    the fairlead moving linearly from its old to its new position, (iii) the body velocity is completed with the new line force
    (linear damping B implicit). The line force is therefore evaluated at the new pose but the pose used the old force: second-order
    for smooth problems (free heave decay error ratio 4.0 per dt halving) with no iteration. The sub-step ratio reported is
    ceil(dt_body / shortest line step).
30. **Stability limit of the coupling.** Measured for eight combinations of EA, N and M (heave of a body on a vertical line): the largest stable
    dt_body is about 0.5-0.9 x 2/omega_s with omega_s^2 = k_line / (M + A + m_line/3), k_line = EA/L (static stiffness of the taut line at the
    fairlead, summed over lines in the load direction); coarse lines (N = 5) sit at the lower end. It is the body step, not the line's CFL step,
    that is limited, so the sub-step ratio can be several hundred without problem. Recommendation: dt_body <= 0.25 * 2/omega_s. For the
    Paredes buoy (M + A ~ 62 kg, k ~ tens of N/m heave from the lines but hydrostatic 2040 N/m) the hydrostatic stiffness dominates, so the
    limit is far above the 1e-3 s one would choose for accuracy. Tests: stable at 0.4 x and unstable at 2 x the estimate.
    The explicit scheme can also be destabilised by a light body on heavy lines (added-mass-like effect: line inertia m_line
    comparable to M + A); this was not explored beyond M = 0.2 kg with a 0.01 kg line.
31. **Equilibrium solve** is a Newton iteration with a finite-difference Jacobian (step 1e-4 m / rad) on F0 - C xi + sum(line loads(xi)),
    solving every line statically at each evaluation (warm start from the previous shape). It needs taut lines at the start: a slack line
    is a one-sided spring with a zero Jacobian (found the hard way in a test). Static line solves use forceTol 1e-7 of the force scale
    max(w l0, m_l g l0, 1e-7 EA), i.e. a few 1e-10 N for the test lines; 1e-9 of the scale is below the round-off floor and never converges.
32. **Static relaxation tolerance scale** changed (also for the single-line runner): it was w*l0 and is now max(w l0, m_l g l0, 1e-7 EA) so
    that weightless taut lines converge; results for weighted lines are unchanged.
33. **Paredes benchmark inputs** (milestone 8): see docs/paredes_validation.md for the list of assumptions that are not in the thesis (line hydrodynamic
    coefficients and diameters, floater/clump drag, yaw inertia, CAT soil, anchor height, horizontal anchor distance convention). The geometry is
    read from Figs 5.4-5.5 and Annex A of the thesis; the labelled rope sections add up consistently (0.700 + 0.900 + 0.685 = 0.700 + 1.585 = 2.285 m).
34. **Constrained equilibrium** (`solveEquilibrium(..., fixed)`): a held DOF carries a reaction and is excluded from the residual; used for the
    surge restoring force (surge imposed, other 5 DOF free) and for the heave/pitch stiffness of the mooring (+-0.01 m, +-1 deg).
35. **Secant stiffness definition**: K(x) = -(F_x(x) - F_x(0))/x with the net line force at the free equilibrium subtracted (it is zero when there is no mean load).
    The thesis' Table 5.12 values carry no stated displacement; comparisons use the curves of its Fig. 5.21 (digitised by eye, about +-1.5 N/m).
36. **Waves** (milestone 9): see docs/waves_validation.md. Linear Airy waves; Wheeler stretching of the kinematics (stretched linear accelerations);
    submergence of nodes and point elements from the still-water level; no second-order forces; body excitation `f = w A sin(omega t + delta)` at the mean position.
37. **Approx tolerance pitfall (found in milestone 9):** doctest's default `Approx` adds an absolute scale of 1.0, so a 'relative' epsilon on quantities of order 0.1
    or smaller is effectively absolute. Milestone 9 tests use `.scale(0)` for such comparisons. All `Approx` uses of milestones 1-8 were re-checked: the low-valued ones
    (soil friction forces 0.035-0.05 N, the ratio 0.7071, periods of 3-6 s) have effective tolerances of about 1e-5 absolute or 0.12-0.27 % relative instead of the nominal
    0.1-0.2 %; the rest compare quantities of order 1 to 1e3 and are unaffected. No quoted result is based on such a tolerance: every number in the README and the docs comes from the
    values printed by the tests, not from a pass/fail threshold.
38. **Steady-state start for moored wave runs:** pose and velocity are initialised from the closed-form response of the free body (6-DOF model without lines)
    because the surge-pitch mode needs ~1000 s to settle from rest. The remaining transient is the mooring-induced one; checked on CAT (60 s from rest vs 30 s from the steady state: same RAOs within 1 %).
39. **Azcona et al. validation** (docs/azcona_validation.md): the amplitude of Table 1 ("0.25 m") is read as the peak-to-peak stroke (+-0.125 m), from the axes of Figs 5/6; seabed
    stiffness 20 N/m^2 per length is Ks*D1 = 20; damping 0.1 Ns/m^2 gives zeta = 0.0426; "structural damping 0.1 %" is not implemented (c_int = 0); v_lim = 0.01 m/s and the friction law (weight-based, tangential only) are
    this code's, not the paper's. Dynamic tensions of the paper exist only as figures: the comparison values were read by eye (+-0.5 N).
40. New runner option `motion.type = "harmonic"` (sinusoidal top-end motion along `direction_vec`, cosine ramp).

41. **Anisotropic point elements.** Point elements may carry separate tangential/normal drag (`Cd_t`, `area_t_m2`, `Cd_n`, `area_n_m2`) and added mass (`added_mass_t_kg`, `added_mass_n_kg`); the tangent is taken from the neighbouring nodes. If these keys are absent the earlier isotropic behaviour is unchanged.
42. **Lopez-Olocco clump.** The clump tangential drag coefficient 1.17 is assumed (not in the source); static CW2 pretension rise is +6.9 % vs the paper's ≈ +10 % (docs/lopezolocco_validation.md).
