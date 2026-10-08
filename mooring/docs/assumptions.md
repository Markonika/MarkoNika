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
