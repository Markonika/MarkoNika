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

