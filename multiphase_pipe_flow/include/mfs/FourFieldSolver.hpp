#pragma once

#include "mfs/Closures.hpp"
#include "mfs/FlowState.hpp"
#include "mfs/FluidProperties.hpp"
#include "mfs/FluxLimiter.hpp"
#include "mfs/PipeGeometry.hpp"

#include <functional>
#include <vector>

namespace mfs {

// Inlet / outlet boundary conditions. The inlet is a fixed-mass-flow
// (Dirichlet velocity + phase distribution) boundary; the outlet is a
// fixed-pressure (Dirichlet pressure, floating velocity) boundary. This is
// the standard, well-posed pairing for a subsonic 1D two-phase pipe flow
// solved with a segregated pressure-correction scheme.
struct BoundaryConditions {
    double inletSuperficialLiquid = 0.5;   // Usl [m/s]
    double inletSuperficialGas = 3.0;      // Usg [m/s]
    double inletLiquidHoldup = 0.25;       // eL at the inlet
    double inletDropletFractionOfLiquid = 0.0; // ed / eL at the inlet
    double inletBubbleFractionOfGas = 0.0;     // eb / eG at the inlet
    double outletPressure = 1.0e5;         // [Pa]

    // Small-amplitude seed disturbance superposed on the inlet holdup to
    // trigger the natural growth of interfacial waves when starting from a
    // perfectly flat stratified interface (Section 4 of the paper notes
    // this is needed so that transition to slug/large-wave flow, when the
    // conditions call for it, is not suppressed by the exact symmetry of a
    // perfectly uniform initial/inlet condition).
    bool seedDisturbance = true;
    double disturbanceAmplitude = 0.03; // relative amplitude on eL
    double disturbanceFrequency = 0.5;  // [Hz]
};

struct SolverOptions {
    FrictionCorrelation wallCorrelation = FrictionCorrelation::TaitelDukler1976;
    FrictionCorrelation interfacialCorrelation = FrictionCorrelation::AndreussiPersen1987;

    // Data-driven correction to the interfacial friction closure's
    // enhancement above the baseline gas-wall friction factor, off by
    // default so existing behaviour is unchanged unless opted into.
    // Motivated directly by this codebase's own 349-case quantitative
    // validation (see VALIDATION.md): predicted holdup tracks measured
    // holdup well in the stratified-like regime (MAE 0.09 on Kokal 1987's
    // 168 cases, 0.11 on Newton 1997's 55) but degrades sharply on cases
    // closer to the slug/annular transition (MAE 0.30 across 21
    // independent Mendeley-archived campaigns, 126 cases) -- exactly the
    // regime where interfacialCorrelation's AndreussiPersen1987 branch
    // applies its enhancement above the baseline gas-wall friction factor
    // (active whenever the Kelvin-Helmholtz parameter F exceeds the
    // correlation's own F0=0.36 threshold). Whatever correlation is
    // selected, let fgw be the gas-wall friction factor and fi the
    // correlation's returned (possibly enhanced) interfacial friction
    // factor; the correction actually applied is
    //   fi_corrected = fgw + closureCorrectionScale * (fi - fgw)
    // -- i.e. it rescales ONLY the enhancement a correlation adds above
    // the no-enhancement baseline, automatically a no-op (fi_corrected ==
    // fi) wherever a case is far enough from the transition that fi==fgw
    // already, and reduces to the UNCORRECTED closure exactly when
    // closureCorrectionScale == 1.0 (its default value) even with the
    // option enabled -- the option's own on/off switch, not this
    // parameter's value, is what determines whether the correction is
    // live, matching this codebase's house convention of gating new
    // behaviour behind an explicit enable flag rather than a
    // defaults-to-neutral parameter alone.
    //
    // closureCorrectionScale was fit by a line search over the held-out
    // validation data itself (train on Kokal 1987 + Newton 1997, 223
    // cases; evaluate on the fully independent 126-case Mendeley set
    // never used for fitting) and found NOT to generalize -- a clean
    // negative result, reported as such rather than tuned to look like
    // a win; see VALIDATION.md for the full train/test methodology and
    // numbers. The default stays at 1.0 (the unmodified closure)
    // because no tested value actually improved the held-out set.
    bool enableDataDrivenClosureCorrection = false;
    double closureCorrectionScale = 1.0;

    // A second, independent knob on the SAME closure, tried after the
    // magnitude rescale above came back negative: rather than rescaling
    // the enhancement's SIZE, shift the Kelvin-Helmholtz onset threshold
    // (F0=0.36) itself that decides WHETHER the enhancement applies at
    // all -- F > F0 + closureF0Shift triggers it instead of F > F0 (see
    // mfs::InterfacialFrictionInputs::f0Shift, and Closures.cpp's
    // AndreussiPersen1987 branch). Also gated behind
    // enableDataDrivenClosureCorrection (both knobs are corrections to
    // the same AndreussiPersen1987 closure, so share the one on/off
    // switch); zero by default, a true no-op. See VALIDATION.md for the
    // same train/test methodology applied to this parameter, and
    // whether, unlike the magnitude rescale, it actually generalized.
    double closureF0Shift = 0.0;

    // A third, independent knob, tried after BOTH knobs above on the
    // interfacial friction closure came back negative: rather than the
    // interfacial correlation, this one scales the droplet DEPOSITION
    // velocity (the rate liquid droplets suspended in the gas core
    // return to the continuous liquid film -- see
    // DropletDepositionInputs/dropletDepositionRate() in Closures.hpp,
    // called with a hardcoded depositionVelocity=0.1 in
    // FourFieldSolver.cpp's computeClosures()). Motivated by the
    // validation driver's own holdup metric, which sums el+ed (the
    // CONTINUOUS liquid fraction plus the DISPERSED droplet fraction):
    // if deposition returns droplets to the film faster than real
    // droplets actually settle, suspended droplet holdup (ed) would be
    // systematically underestimated, a plausible contributor (among
    // others) to the solver's documented under-prediction of holdup in
    // the slug-dominated regime (see VALIDATION.md). Also gated behind
    // enableDataDrivenClosureCorrection; 1.0 by default, a true no-op
    // (depositionVelocity_used = 0.1 * closureDepositionVelocityScale).
    double closureDepositionVelocityScale = 1.0;

    double courantTarget = 0.5;    // target Courant number, Eq. (23) requires < 1
    double minTimeStep = 1.0e-6;
    double maxTimeStep = 5.0e-2;
    double wallRoughness = 0.0;    // [m], for rough-wall correlations

    // Under-relaxation applied to the pressure correction (and the
    // resulting velocity correction), standard SIMPLE-family practice.
    // Without it, cells where a phase's volume fraction is momentarily
    // near zero (e.g. the thin gas layer under a forming liquid slug) give
    // the linearized pressure-correction system a very weak local
    // compressibility coefficient there, which can overshoot to
    // nonphysical (even negative) pressure in a single step.
    double pressureRelaxation = 0.2;
    double minPressure = 1.0e3; // [Pa] hard floor, applied after each correction
    double maxPressureChangeFraction = 0.05; // hard cap: |dP| <= this * P, per step
    // Absolute ceiling, as a multiple of the outlet pressure: a blunt but
    // effective safety net against runaway divergence compounding over
    // many steps in a numerically stiff cell (see FourFieldSolver.cpp,
    // solvePressureCorrection()), independent of the per-step relative cap.
    double maxPressureFactor = 20.0;

    // Defensive velocity limiter. The dispersed-phase algebraic momentum
    // balance (Eq. 10) is well-posed in the entrainment-dominated limit
    // (see the derivation in FourFieldSolver.cpp, updateDispersedMomentum)
    // but can still become locally stiff -- e.g. a cell carrying advected
    // droplet/bubble fraction with negligible *local* entrainment and a
    // degenerate drag coefficient -- and momentarily overshoot to an
    // unphysical velocity before the next step's drag term pulls it back.
    // Clamping to a generous multiple of any physically expected pipe-flow
    // velocity (an order of magnitude above the ~90 m/s peaks reported in
    // the source paper's own large-wave/slug-formation transients) keeps
    // such transients bounded without affecting normal operating ranges.
    double maxVelocity = 100.0; // [m/s]

    // Floor applied to e1/e2 ONLY where they appear as a denominator in the
    // layer momentum equations (Eq. 6-7), as opposed to their role
    // elsewhere as the actual tracked volume fraction (which uses the
    // finer Constants::small_e floor). Momentum-equation source terms
    // (wall/interfacial shear, entrainment momentum flux) do not shrink
    // proportionally as a layer's cross-section is pinched down -- they
    // are set by the OTHER phase's conditions -- so dividing by an
    // already-tiny e1/e2 can amplify them well past what the velocity
    // limiter alone comfortably absorbs. A dedicated, coarser floor here
    // keeps that amplification bounded (to 1/momentumFractionFloor)
    // without perturbing the volume-fraction fields themselves.
    double momentumFractionFloor = 0.05;

    // Higher-order / flux-limited (MUSCL/TVD, Sweby 1984 form) face-value
    // reconstruction for the field continuity equations' advective fluxes
    // (updateContinuity(), the el/ed/eg/eb transport), in place of the
    // plain first-order donor-cell upwinding used when this is None (the
    // default, so existing behaviour and the VALIDATION.md numbers are
    // unchanged unless this is opted into). First-order upwind is highly
    // numerically diffusive on the wave/front-sharpening physics this model
    // targets -- exactly the mechanism flagged in VALIDATION.md as a
    // candidate explanation for why some F>F0 (Kelvin-Helmholtz unstable)
    // conditions fail to grow into the slug/roll-wave regime in this
    // explicit time-marching scheme. A limited high-resolution
    // reconstruction is formally second-order in smooth regions while
    // remaining Total-Variation-Diminishing (no new overshoot/oscillation)
    // at a front, via the limiter itself falling back to first order there.
    // The pressure-correction system (solvePressureCorrection()) is
    // deliberately left on plain upwind regardless of this setting -- a
    // standard "deferred correction" split: the correction step only needs
    // a stable, well-conditioned linearization to drive the pressure
    // iteration, not the final transport accuracy, and its Gf coefficients
    // were tuned/validated against the robustness work in this file.
    // Van Leer is a reasonable general-purpose default once enabled
    // (less compressive/more diffusive than Superbee, so less prone to
    // artificially steepening an already-sharp front, but still markedly
    // less diffusive than first-order upwind).
    FluxLimiterType advectionLimiter = FluxLimiterType::None;

    // Surface-tension-motivated short-wave regularization of the field
    // continuity equation, off by default (existing behaviour unchanged
    // unless opted into). The paper's own Eq. 6-7 (and hence this
    // implementation) carry no interfacial-curvature term, which leaves
    // the continuum model short-wave ILL-POSED: linearizing the inviscid
    // layer equations around a uniform stratified base state gives a
    // growth rate that increases without bound as wavelength shrinks (see
    // VALIDATION.md, "Recommended follow-up" item 1, Update 5, and the
    // README section this option is documented in) -- the same short-wave
    // ill-posedness documented for two-fluid models lacking such a term
    // (Stewart & Wendroff 1984; Ramshaw & Trapp 1978).
    //
    // The physically literal fix -- an interfacial pressure jump
    // (Young-Laplace, P2 = P1 + sigma*d^2(h1)/dz^2) entering the gas
    // momentum equation as -sigma/rho2 * d^3(h1)/dz^3 -- was tried and
    // abandoned: a third derivative is a DISPERSIVE term, and (as verified
    // both by direct von Neumann analysis and by empirical blow-up in
    // testing) no one-sided/upwind-biased discretization of it is
    // unconditionally stable under this solver's explicit forward-Euler
    // time-stepping the way upwind differencing of an ADVECTIVE
    // (first-derivative) term is -- a real, non-obvious pitfall of that
    // approach, not just a matter of shrinking dt further. See
    // VALIDATION.md item 1, Update 6 for the full account.
    //
    // What's implemented instead is a biharmonic ("hyperdiffusion") proxy
    // added directly to the liquid-holdup continuity equation:
    // d(eL)/dt += -nu4 * d^4(eL)/dz^4, with the local coefficient
    // nu4 = hyperdiffusionCoefficient * |ul| * dz^3 (dz = local cell
    // width, ul = local liquid velocity). This is NOT a literal capillary-
    // wave model -- it doesn't use fluid.sigma at all -- but a fourth
    // derivative has real, non-positive eigenvalues everywhere (confirmed
    // by von Neumann analysis: the standard centred 4-point stencil's
    // symbol is 16*sin^4(theta/2)/dz^4, always >= 0, so as a damping term
    // -nu4*D4 it is genuinely diffusive, not dispersive), giving a
    // standard, unconditionally-safe-for-a-finite-dt explicit stability
    // bound with none of the third derivative's pitfalls. Verified (via
    // the same linearized dispersion relation used throughout this
    // investigation, evaluated numerically, not just asserted) to
    // reproduce the qualitatively right behaviour: strong, monotonically
    // increasing suppression of growth at short wavelength, negligible
    // effect at long wavelength -- capping this model's otherwise-
    // unbounded short-wave growth rate the same way real surface tension
    // would, without claiming to model surface tension itself. See
    // FourFieldSolver.cpp, updateContinuity() and fourthDerivativeAt(),
    // and stableTimeStep() for the (standard, diffusion-type) stability
    // cap this also engages once enabled.
    bool enableSurfaceTension = false;
    double hyperdiffusionCoefficient = 5.0e-4; // dimensionless; see above

    // Physics-based turbulent-viscosity regularization of the RELATIVE
    // velocity between the two layers (u2 - u1), off by default. This is a
    // separate mechanism from enableSurfaceTension above and independently
    // toggleable: the two are complementary, not alternatives.
    // enableSurfaceTension's biharmonic term targets this model's LINEAR
    // short-wave ill-posedness; this option targets its NONLINEAR
    // behaviour once a wave has already grown into a sharp front -- the
    // "stuck case" investigated in the KH case study (VALIDATION.md item
    // 1) never gets there, so the two options address different stages of
    // the same underlying problem.
    //
    // Physical basis: Lopez-de-Bertodano & Clausse, "Nonlinear Stability
    // in the Two-Fluid Model of Two-Phase Flow" (Physics of Fluids, 2026;
    // arXiv:2509.04679), derive a turbulent eddy viscosity
    //   nu_t = l_m * |u_r|                                    (their Eq. 14)
    // (l_m a constant "mixing length", u_r the relative velocity between
    // phases) and add it, alongside the kinematic viscosity, as a
    // diffusive force on the relative-momentum equation:
    //   F_visc,W = d/dx(nu * dW/dx),  nu = nu_k + nu_t         (their Eq. 11-13)
    // In their simplest (no-inertial-coupling) reduction, the relative
    // momentum variable W reduces exactly to u_r, so this is literally a
    // diffusion of the slip velocity between phases. They show
    // analytically (their Eq. 90-96) that this term's viscous force scales
    // as (mesh spacing)^-2 as a forming wave narrows, which beats the
    // (mesh spacing)^-1 scaling of the advective growth term -- unlike
    // surface tension, whose contribution is independent of the narrowing
    // width and therefore fails to arrest the nonlinear "shock-spike"
    // blowup their surface-tension-only model still exhibits (their
    // Section V.A/V.C). With both surface tension (linear) and this
    // turbulent viscosity (nonlinear) together, their model produces
    // KH-instability-triggered SLUG FLOW emerging with no imposed slug
    // structure (their Section V.D, Fig. 19) -- precisely the qualitative
    // behaviour ("stuck" stratified flow that never develops into
    // slug/annular) this codebase's own validation has repeatedly found
    // missing.
    //
    // Translation to this solver's four-field, dimensional, distinct-
    // density formulation (their model is a simplified, implicitly
    // non-dimensionalized two-fluid square-channel formulation, not
    // literally this codebase's own equations, so this is an honest
    // adaptation of their closure, not a line-by-line reproduction): the
    // diffusive acceleration on the relative layer velocity is
    //   a_diff = nu_t * d^2(u2-u1)/dz^2,   nu_t = l_m * |u2-u1|
    // with l_m = turbulentMixingLengthFraction * D (a diameter-based
    // length scale, consistent with how every other closure in this
    // codebase non-dimensionalizes length -- their own toy test used a
    // fixed l_m=1mm in an idealized unit-density square channel, which
    // does not transfer to this solver's dimensional air-water pipe
    // scales). This acceleration is then split between the two layers by
    // reduced-mass weighting so the shared coupling conserves total
    // momentum exactly while reducing to a_diff on the relative velocity:
    //   du1/dt += -a_diff * m2/(m1+m2),  du2/dt += +a_diff * m1/(m1+m2)
    // with m1 = rho1*e1, m2 = rho2*e2 (per-unit-volume layer "masses" at
    // the face). See FourFieldSolver.cpp, updateLayerMomentum() and
    // secondDerivativeAtFace(), and stableTimeStep() for the standard
    // explicit-diffusion stability cap (dt <= dz^2/(2*nu_t), verified by
    // von Neumann analysis of the same 3-point stencil used here) this
    // also engages once enabled.
    bool enableTurbulentViscosity = false;
    double turbulentMixingLengthFraction = 0.1; // l_m = this * D; dimensionless

    // Implicit (IMEX) treatment of wall and interfacial friction in the
    // layer momentum equations, off by default so existing behaviour is
    // unchanged unless opted into. These are the classic stiff source
    // terms in this model class: their coefficients (friction factors x
    // density x |velocity|) can be large -- especially interfacial drag
    // near a flow-regime transition -- and under pure explicit forward-
    // Euler that stiffness forces dt down independent of, and often well
    // below, the material (advective) Courant limit that
    // SolverOptions::courantTarget alone accounts for. Standard practice
    // in this class of two-fluid code (e.g. RELAP5/TRAC/CATHARE-family
    // system codes) is to treat friction implicitly while leaving
    // advection explicit or separately semi-implicit -- the IMEX split
    // this option performs.
    //
    // Implementation: at each face, the wall-friction and interfacial-
    // friction accelerations are linearised by freezing their nonlinear
    // parts (friction factor, and the relative/absolute velocity
    // magnitude that makes wall and interfacial shear stress quadratic in
    // velocity) at the OLD time level, extracting an effective LINEAR
    // drag rate k = (old acceleration)/(old velocity or relative
    // velocity). This turns the two layer-momentum equations, restricted
    // to just their friction terms, into an exactly-solvable local 2x2
    // linear system in (u1_new, u2_new) at that face (no spatial coupling
    // between faces, so no banded/tridiagonal solve is needed -- unlike
    // enableTurbulentViscosity's genuinely spatial diffusion term, which
    // remains explicit with its own stability cap regardless of this
    // option). All other terms (advection, pressure gradient, gravity,
    // mass-transfer sources, and the turbulent-viscosity term if enabled)
    // stay explicit as before. See FourFieldSolver.cpp,
    // updateLayerMomentum() for the derivation and the exact 2x2 solve.
    //
    // Because backward-Euler friction is unconditionally stable in the
    // drag-coefficient magnitude, enabling this removes any dt
    // restriction that magnitude alone would otherwise impose; the
    // material Courant limit on advection (SolverOptions::courantTarget)
    // is untouched by this option and remains the governing constraint.
    // A full SETS-type treatment that also relaxes the material Courant
    // limit itself (by additionally semi-implicating the mass/pressure
    // propagation, not just the friction source terms) is a larger,
    // separate undertaking not attempted here.
    bool enableImplicitFriction = false;

    // Exponential time differencing (ETD) treatment of wall and
    // interfacial friction, off by default. An alternative to
    // enableImplicitFriction, not a replacement: both linearize the
    // SAME friction terms the SAME way (freeze the friction factor and
    // the |velocity| that makes each shear stress quadratic at the OLD
    // state, giving an effective linear drag-rate matrix
    // A = [[-(kw1+ki1), ki1], [ki2, -(kw2+ki2)]] and constant forcing
    // b = (E1, E2), exactly as enableImplicitFriction derives it -- see
    // that option's own documentation for the derivation of kw1, kw2,
    // ki1, ki2), but integrate the resulting LINEAR ODE
    // dU/dt = A*U + b (U = (u1, u2)) EXACTLY over the step via the
    // matrix exponential, instead of backward-Euler's first-order
    // approximation of it. For any constant-coefficient linear ODE this
    // is the best possible local accuracy a method built on this exact
    // linearization can achieve (zero truncation error in the
    // linearized system itself; all remaining error comes from the
    // linearization -- freezing the friction factor/|velocity| at the
    // old state -- shared identically with enableImplicitFriction), at
    // the same O(1)-per-face cost: a 2x2 matrix function evaluated in
    // closed form, no linear solve needed at all.
    //
    // Closed form: exp(A*dt) = alpha0*I + alpha1*A, and the forcing
    // integral integral_0^dt exp(A*s) ds = beta0*I + beta1*A, with
    // alpha0, alpha1, beta0, beta1 built from A's eigenvalues via the
    // standard divided-difference (confluent Vandermonde) representation
    // of a 2x2 matrix function -- see FourFieldSolver.cpp,
    // expAndPhiCoeffs2x2(). A's eigenvalues are GUARANTEED real and
    // non-positive for this specific matrix: writing s1=kw1+ki1,
    // s2=kw2+ki2 (both >=0), the trace is -(s1+s2) <= 0 and the
    // discriminant (s1-s2)^2 + 4*ki1*ki2 is a sum of squares, hence
    // always >= 0 -- so no complex-eigenvalue branch is needed, and the
    // repeated-eigenvalue case (a measure-zero but reachable edge case,
    // e.g. ki1 == ki2 == 0 with kw1 == kw2) is handled by its own
    // confluent limit formula, not skipped.
    bool enableETDFriction = false;

    // Multi-stage IMEX Runge-Kutta treatment of wall and interfacial
    // friction, off by default. A third alternative alongside
    // enableImplicitFriction and enableETDFriction, all three sharing
    // the SAME linearization (same A, b as described above) and
    // differing only in how that linear ODE is integrated over the
    // step: enableImplicitFriction uses single-stage backward Euler
    // (1st order in the linearized system), enableETDFriction the exact
    // matrix exponential (zero error in the linearized system), and this
    // option a 2-stage, 2nd-order, L-stable diagonally-implicit
    // Runge-Kutta (SDIRK) scheme for the implicit part, paired with the
    // frozen explicit forcing below -- an instance of the IMEX-RK
    // framework surveyed from Pareschi & Russo (2005, "Implicit-explicit
    // Runge-Kutta schemes and applications to hyperbolic systems with
    // relaxation", J. Sci. Comput. 25), the scheme that survey explicitly
    // names as "a more rigorous completion" of the single-stage split
    // already in place. Butcher tableau: c=(gamma,1),
    // A_im=[[gamma,0],[1-gamma,gamma]], b=(1-gamma,gamma), with
    // gamma = 1 - 1/sqrt(2):
    //   (I - dt*gamma*A) U1 = U0
    //   (I - dt*gamma*A) U2 = U0 + dt*b + dt*(1-gamma)*A*U1
    //   U_new = U0 + dt*b + dt*[(1-gamma)*A*U1 + gamma*A*U2]
    // L-stability was confirmed directly, not assumed: substituting a
    // scalar test equation dy/dt=lambda*y gives the stability function
    // R(z) = (1 - gamma^2*z) / (1 - gamma*z)^2 for z=dt*lambda, which
    // ->0 as z->-infinity, as required. An EARLIER version of this
    // scheme used a21=(1-2*gamma) and b=(1/2,1/2) instead of the
    // (1-gamma, gamma) above -- a plausible-looking but WRONG tableau
    // (its own R(infinity) works out to -1, not L-stable at all) caught
    // not by re-deriving the stability function first, but because the
    // direct stiff dt-sweep test (see README.md) showed it blowing up at
    // exactly the dt where enableImplicitFriction's plain backward Euler
    // stays stable -- the opposite of what a strictly more accurate
    // integrator of the same linear system should ever do, which is what
    // prompted re-deriving R(z) and finding the actual error.
    //
    // Both stage solves share the same 2x2 matrix (I - dt*gamma*A), only
    // the right-hand side differs, so this costs two small linear solves
    // per face (same closed-form Cramer's-rule style as
    // enableImplicitFriction's single solve) plus one explicit
    // combination -- still O(1) per face, no iteration.
    //
    // SCOPE, stated explicitly because it is the one simplification that
    // distinguishes this from a genuine whole-equation multi-stage
    // scheme: only the LOCAL, per-face friction ODE is advanced through
    // the two stages. The non-stiff forcing b = (E1, E2)
    // (advection, pressure gradient, gravity, mass-transfer and
    // turbulent-viscosity sources) is frozen at the start-of-step state
    // throughout BOTH stages -- it is never re-evaluated at an
    // intermediate stage value, which would require re-coupling with the
    // pressure-correction and advection steps at each stage, a much
    // larger undertaking outside this option's scope. This is the same
    // kind of explicit scoping restriction enableImplicitFriction and
    // the JFNK solver's own lagged source terms already use, applied
    // here for the same reason: it keeps the new code self-contained and
    // independently testable without touching the rest of the time-
    // advance pipeline.
    bool enableIMEXRKFriction = false;

    // Well-balanced discretization of the interface-slope gravity term
    // (g*cos(theta)*dh1/dz + g*sin(theta) in Eq. 6-7), off by default so
    // existing behaviour is unchanged unless opted into. Relevant ONLY to
    // terrain-following runs where theta varies along the pipe
    // (setInclinationProfile() with a non-constant profile) -- for any
    // constant-theta run (the overwhelming majority of this codebase's own
    // validation sweeps) this option changes essentially nothing, since the
    // defect it fixes vanishes identically when theta doesn't vary.
    //
    // The defect: at a true static (zero-velocity) equilibrium in a
    // terrain-following pipe, setting the layer-momentum equations' RHS to
    // zero gives dh1/dz = -tan(theta(z)) -- the liquid height within the
    // pipe cross-section must vary along z to keep the physical free
    // surface level as the pipe tilts. The default discretization evaluates
    // this term from a face-averaged theta (thetaF = 0.5*(theta[cL]+
    // theta[cR])) applied to a plain centred difference of the two
    // neighbouring cells' h1 values -- consistent (convergent as the mesh
    // refines) but NOT exact at any fixed finite resolution wherever theta
    // differs between adjacent cells: initializing the solver at the true
    // equilibrium and stepping forward with zero inflow produces a small,
    // spurious velocity drift confirmed (by direct test, not just derived)
    // to shrink at close to first order as the mesh is refined -- i.e.
    // ordinary truncation error, not a structural non-convergence, but
    // still a real, avoidable discrepancy specifically in the
    // terrain-following capability Section 5 of the paper demonstrates.
    //
    // The fix (a discrete analogue of hydrostatic reconstruction, e.g.
    // Audusse et al. 2004 for shallow water, adapted to this term's
    // algebraic structure rather than copied from it): build a per-cell
    // REFERENCE profile h1Eq(i) that exactly satisfies the equilibrium
    // relation using each cell's OWN theta (not a face average), by a
    // simple O(N) running sum from an arbitrary anchor (only differences of
    // h1Eq ever enter the scheme, so the anchor itself is irrelevant).
    // Replacing the bare h1 gradient with the gradient of the DEVIATION
    // eta = h1 - h1Eq, and folding the separate sin(theta) term into that
    // same replacement, is an exact algebraic identity in continuous form
    // (g*cos(theta)*d(eta)/dz == g*cos(theta)*dh1/dz + g*sin(theta) when
    // dh1Eq/dz=-tan(theta) exactly) -- and, built this way, is ALSO exact at
    // the discrete level: if the state is initialised with h1 == h1Eq
    // pointwise (the true per-cell equilibrium), eta is identically zero at
    // every cell center, so its gradient -- and hence the entire term -- is
    // exactly zero at every face, for any mesh coarseness and any theta
    // profile. See FourFieldSolver.cpp, computeH1Equilibrium() and
    // updateLayerMomentum().
    bool enableWellBalancedGravity = false;

    // Interfacial pressure-jump regularization (Stewart & Wendroff, 1984;
    // Toumi & Kumbaro, 1996; adapted here to this solver's two-LAYER,
    // not general multi-fluid, acceleration-form momentum equations), off
    // by default so existing behaviour is unchanged unless opted into.
    // Found via literature search specifically for mechanisms that give a
    // continuous, area-averaged two-fluid model an actual slug/
    // intermittency onset -- see README.md, "Interfacial pressure-jump
    // regularization" for the full account of the search and the result.
    //
    // NOT the same mechanism as enableSurfaceTension's abandoned Young-
    // Laplace attempt (see that option's own documentation), even though
    // both are called "an interfacial pressure jump" in parts of the
    // literature -- an easy mix-up this comment exists partly to head
    // off. The Young-Laplace term is P2-P1 = sigma*d^2(h1)/dz^2, a SPATIAL
    // CURVATURE term requiring a third-derivative stencil in the momentum
    // equation, which was tried and abandoned because no one-sided
    // discretization of a dispersive third derivative is unconditionally
    // stable under this solver's explicit time-stepping. This option's
    // term is a purely LOCAL, ALGEBRAIC closure (no new derivative order
    // at all -- it only modifies what pressure value the EXISTING
    // first-derivative pressure-gradient term in each layer's own
    // momentum equation reads), derived from how the two layers' local
    // relative motion past the interface gives them physically different
    // local pressures (a Bernoulli/streamline-curvature effect), not from
    // surface tension.
    //
    // Derivation: let Dp = Ci * rho1*rho2/(rho1*e2+rho2*e1) * (u1-u2)^2
    // (u1, u2, e1, e2 all per-cell; rho1, rho2 the per-cell layer mixture
    // densities already computed by rho1Of/rho2Of). Define per-layer
    // pressures P1 = P - e2*Dp, P2 = P + e1*Dp (P the existing shared
    // mixture pressure this solver already solves for via
    // solvePressureCorrection(), Eq. 22 -- UNCHANGED by this option: since
    // e1+e2=1 identically, e1*P1+e2*P2 = e1*P-e1*e2*Dp+e2*P+e1*e2*Dp = P,
    // so the mixture/total pressure equation's own derivation is
    // unaffected and does not need to change). Each layer's own momentum
    // equation then reads its OWN layer's pressure gradient
    // (dP1/dz/rho1, dP2/dz/rho2) instead of the shared dP/dz/rho1,
    // dP/dz/rho2 both layers use when this is off -- see
    // FourFieldSolver.cpp, updateLayerMomentum(). When Ci=0 (or this
    // option is off), P1==P2==P exactly and every downstream use is
    // bit-for-bit identical to the pre-existing shared-pressure code.
    //
    // Verified by re-deriving this solver's own linearized inviscid KH
    // dispersion relation (the same analysis used throughout VALIDATION.md
    // item 1 and the slug-capturing investigation) WITH this term
    // included, via sympy, and confirming two things before writing any
    // solver code: (1) at Ci=0 the extended derivation reproduces the
    // pre-existing (no-term) dispersion relation EXACTLY (checked
    // numerically, ratio 1.0 to machine precision, across several
    // wavenumbers and base states) -- the re-derivation itself is sound,
    // not just the Ci=0 special case being trivially right; (2) because
    // e1+e2=1 identically, P2-P1 collapses to exactly Dp (no further
    // algebra needed), and Dp's own linearization contributes to the
    // dispersion relation's omega^1 and omega^0 coefficients but NOT its
    // omega^2 coefficient -- meaning this term does NOT change the
    // existing short-wave (k -> infinity) UNBOUNDED growth-rate scaling
    // (still grows linearly in k, same role as items 1/38-47's surface
    // tension and biharmonic terms, left untouched by this option) but
    // DOES shift the discriminant's SIGN uniformly across all
    // wavelengths, i.e. this is an ONSET-threshold mechanism (does
    // instability exist at this base state at all), complementary to
    // rather than competing with those short-wave-capping terms.
    //
    // Evaluated at this codebase's own three worst-under-predicted
    // Mendeley cases (Brito 2012, Baba et al. 2017, Ekinci 2015 -- same
    // cases as the slug-capturing investigation) across the literature's
    // own recommended range (Evje & Flatten assume 1 < Ci <= 2 for
    // guaranteed hyperbolicity in the general multi-fluid case): the
    // linear KH onset holdup drops substantially as Ci increases from 0
    // (where it reproduces the ~0.70-0.73 onset the slug_kh_check.cpp
    // diagnostic already found) down to roughly 0.40-0.48 at Ci=1, but
    // stays well ABOVE the solver's own predicted thin-film holdup for
    // these cases (0.02-0.05) even at Ci=2 -- i.e. in its literature-
    // motivated range, this term does NOT by itself linearly destabilize
    // the thin-film branch the solver actually converges to; whatever
    // effect it has on the documented holdup under-prediction gap is
    // necessarily a NONLINEAR one (interacting with the fine-resolution/
    // AMR/bistable-branch dynamics the slug-capturing investigation
    // already found), to be checked empirically, not assumed from this
    // linear result alone. See VALIDATION.md for the full numeric
    // account and whatever empirical result follows.
    bool enableInterfacialPressureJump = false;
    double interfacialPressureJumpCoefficient = 1.5; // Ci, dimensionless; see above

    // Adaptive (non-uniform) mesh refinement, off by default so existing
    // behaviour at a fixed uniform resolution is unchanged unless opted
    // into. When enabled, every `adaptEveryNSteps` steps each cell's
    // Kelvin-Helmholtz stability parameter F (the same quantity computed
    // for the Andreussi & Persen interfacial friction closure, see
    // Closures.hpp) is used as the refinement indicator: cells where the
    // interface is going unstable (F approaching or past the closure's own
    // F0 threshold) are split; cells deep in a stable, quiescent state are
    // merged with a like neighbour. This follows the same indicator used
    // by Gourma, Jia & Thompson (2013), "Two-Fluid Model for 1D Gas-Liquid
    // Slug Flows: Realizable Mean Slug Characteristics", Multiphase
    // Science and Technology 25(1), for AMR on this class of 1D two-fluid
    // model. See FourFieldSolver.cpp, adaptMesh(), for the refine/coarsen
    // mechanics (a plain non-uniform 1D grid, no hanging nodes).
    struct AdaptiveMeshOptions {
        bool enabled = false;
        int adaptEveryNSteps = 20;
        double refineThreshold = 0.36;    // F0 (Andreussi & Persen 1987)
        double coarsenThreshold = 0.15;   // hysteresis gap below F0, avoids refine/coarsen chatter
        // Default to 1.0: don't refine FINER than the initial mesh (which
        // the user presumably already chose to be adequate, e.g. ~1
        // diameter per the paper's own convergence finding) -- AMR's job
        // by default is purely to coarsen quiescent regions, so it can
        // only reduce cost relative to a fixed run at the initial
        // resolution, never increase it by over-refining. Lower this
        // explicitly to refine beyond the initial resolution in hot zones.
        double minCellWidthFraction = 1.0; // as a fraction of the INITIAL uniform spacing
        double maxCellWidthFraction = 4.0;  // as a fraction of the INITIAL uniform spacing
        double maxCellCountFactor = 4.0;    // cap on N, as a multiple of the initial N
    } amr;

    // Moving-mesh (r-adaptive) node tracking: the literature's other AMR
    // strategy for this class of model (see e.g. Nydal & Banerjee 1996; De
    // Leebeeck 2010, "A roll wave and slug tracking scheme for gas-liquid
    // pipe flow"), and architecturally distinct from `amr` above. Where
    // h-refinement keeps a fixed background resolution and adds/removes
    // cells within a budget, this keeps the cell COUNT fixed and instead
    // continuously relocates all N node positions toward wherever a
    // monitor function is largest, via the classical equidistribution
    // principle (de Boor, 1974): the new mesh is chosen so that the
    // integral of the monitor function is equal over every cell. Because
    // resolution is limited only by the total point budget N, not by any
    // minimum-cell-width floor, this can concentrate resolution on a sharp
    // front far more tightly than h-refinement (whose default
    // `minCellWidthFraction` explicitly forbids refining past the initial
    // spacing) -- at the cost of remapping the entire field onto the new
    // mesh every relocation, and of losing the "coarsen everywhere else"
    // computational saving h-refinement gets from actually reducing cell
    // count. See FourFieldSolver.cpp, computeMonitorFunction() and
    // relocateMesh(), and mutually exclusive with `amr` by convention
    // (both change the mesh; combining them is untested).
    struct MovingMeshOptions {
        bool enabled = false;
        int relocateEveryNSteps = 1; // true tracking needs to keep pace with front motion
        double relaxation = 0.5;     // under-relaxation of node motion per relocation, in [0,1]

        // Monitor function M(z) = 1 + holdupGradientWeight*|d(eL)/dz|*L
        //                            + khIndicatorWeight*min(F/F0, khIndicatorCap)
        // The first term is the classical "arc-length" monitor (Huang &
        // Russell) that concentrates points where the liquid holdup itself
        // is changing sharply -- i.e. directly on a wave/slug front once
        // one exists. The second, smaller term adds anticipatory pull
        // toward cells the interfacial closure already flags as unstable
        // (the same F used by `amr`), attracting points to where a front
        // is ABOUT to form, before its holdup gradient is yet sharp.
        double holdupGradientWeight = 8.0;
        double khIndicatorWeight = 0.5;
        double khIndicatorCap = 5.0;      // caps the F/F0 contribution, as a multiple of F0
        double monitorCap = 50.0;         // hard cap on M itself, bounding max achievable clustering
        int monitorSmoothingPasses = 2;   // 3-point smoothing of M before equidistributing (standard practice, avoids noisy/oscillatory mesh motion)
    } movingMesh;
};

// Core solver implementing the simplified four-field model of Bonizzi,
// Andreussi & Banerjee (2009):
//  - Eq. (6),(7): layer-1 / layer-2 mixture momentum, WITHOUT the slip flux
//    term (shown in the paper's Appendix A to be an order of magnitude or
//    more smaller than the other terms retained).
//  - Eq. (10): dispersed-phase (droplet, bubble) momentum, reduced to an
//    algebraic balance because the inertial terms are negligible (as shown
//    a posteriori in the paper).
//  - Eq. (11)-(16): the four field continuity equations plus the total
//    liquid/gas continuity equations used for mass bookkeeping.
//  - Eq. (22): the pressure equation obtained from the combined total
//    liquid/gas continuity equations, solved with a standard segregated
//    pressure-correction (SIMPLE-family) scheme, per Ferziger & Peric
//    (1999) as cited in the paper.
//
// Generalisation to arbitrary inclination: the paper's own momentum
// equations already carry a general inclination angle theta (the sin(theta)
// gravity-along-pipe term and the cos(theta) hydrostatic "level-gradient"
// term). This implementation keeps theta as a per-cell field, so it applies
// unchanged from horizontal (theta=0) through any incline up to vertical
// (theta=+-90 deg), and can vary along the pipe to represent terrain
// (exactly the motivating scenario the paper's introduction describes).
// At theta=90 deg the cos(theta) "level-gradient" term vanishes identically
// -- there is no longer a preferred low side of the pipe to drive
// gravity-segregated stratification -- which is the correct physical limit.
// The closures in Table 1-3, however, were derived and validated for
// near-horizontal stratified/slug/bubbly/annular flows; applying them
// unchanged at steep or vertical inclination is an extrapolation beyond the
// paper's validated range (the paper's own conclusion flags this
// extension -- vertical pipelines -- as future work). See README.md.
class FourFieldSolver {
public:
    FourFieldSolver(double diameter, double length, int nCells,
                     FluidProperties fluid, SolverOptions options = {});

    // Inclination profile theta(z) [rad from horizontal]; positive = uphill
    // in the flow direction. Default: horizontal everywhere.
    void setInclinationProfile(const std::function<double(double z)>& thetaOfZ);
    void setInclinationConstant(double thetaRad);

    void setBoundaryConditions(const BoundaryConditions& bc);

    // Initialise a uniform stratified field at the given total liquid
    // holdup eL0, all liquid continuous / all gas continuous (el=eL0,
    // eg=1-eL0, ed=eb=0), with velocities from the inlet superficial
    // velocities.
    void initializeStratified(double eL0);

    // Advance one explicit (+ implicit pressure-correction) time step of
    // size dt. Returns the actual dt used (may be reduced internally for
    // robustness, though normally equals the input).
    double step(double dt);

    // Eq. (23): dt such that max(|u|)*dt/dz = courantTarget.
    double stableTimeStep() const;

    // Jacobian-free Newton-Krylov (JFNK) fully implicit alternative to
    // step()'s segregated explicit/semi-implicit pressure-velocity update,
    // addressing the pressure-velocity coupling limitation flagged (but
    // left unaddressed) by this codebase's own earlier investigation (see
    // VALIDATION.md, "Recommended follow-up" item 4 and README.md,
    // "Pressure-velocity coupling investigation"). Backward-Euler, solved
    // simultaneously and nonlinearly for the layer velocities (u1, u2) and
    // pressure (P) via Newton's method, each Newton step's linear solve
    // performed by restarted GMRES using only matrix-free, finite-difference
    // directional derivatives of the residual -- no Jacobian is ever formed
    // or stored. See FourFieldSolverImplicit.cpp for the full derivation,
    // residual definition, and scoping notes (volume fractions, dispersed-
    // phase velocities, and mass-transfer/entrainment source terms remain
    // explicit/lagged from the start-of-step state, exactly as in step();
    // only the layer-velocity/pressure subsystem -- advection, the pressure
    // gradient, gravity, and wall/interfacial friction, all evaluated
    // nonlinearly at the trial state, no hand-linearization -- is implicit).
    //
    // Unlike step(), this does NOT silently fall back or clamp through a
    // non-converged solve: ImplicitStepResult::converged reports whether
    // Newton reached options_.jfnk.newtonTol within
    // options_.jfnk.maxNewtonIters, and the solver's state is left
    // UNCHANGED (as if this call never happened) if it did not -- the
    // caller is expected to retry at a smaller dt, standard
    // pseudo-transient-continuation practice for this class of method
    // (Knoll & Keyes, 2004, "Jacobian-free Newton-Krylov methods: a survey
    // of approaches and applications", J. Comput. Phys. 193).
    struct JFNKOptions {
        int maxNewtonIters = 30;
        // ||R||_2, relative to the initial residual norm. 1e-6 (not a
        // tighter value like 1e-8) is a deliberate choice, confirmed by
        // direct test: a first-order forward-difference Jacobian-free
        // directional derivative has its own truncation-error floor, and
        // demanding many more orders of magnitude of reduction than that
        // floor allows just burns Newton iterations without further
        // progress -- see stepImplicitPressureVelocityAdaptive()'s
        // documentation for the full account of how this was diagnosed.
        double newtonTol = 1.0e-6;
        int gmresRestart = 30;
        int gmresMaxIters = 60;
        double gmresTol = 1.0e-6;   // relative linear-residual tolerance per Newton step
        double fdEpsilon = 1.0e-7;  // relative finite-difference step for J*v
    } jfnk;

    struct ImplicitStepResult {
        bool converged = false;
        int newtonIterations = 0;
        int totalGmresIterations = 0;
        double initialResidualNorm = 0.0;
        double finalResidualNorm = 0.0;
    };
    ImplicitStepResult stepImplicitPressureVelocity(double dt);

    // Pseudo-transient-continuation wrapper around
    // stepImplicitPressureVelocity(): if a sub-step's Newton solve fails
    // to converge, halve it and retry, accumulating converged sub-steps
    // until the full requested dt is covered (or the sub-step width falls
    // below options_.minTimeStep, at which point it gives up and leaves
    // the solver exactly at the end of the last successfully converged
    // sub-step -- never at a non-converged, potentially garbage state).
    //
    // This is not a cosmetic add-on: it is the standard, expected way to
    // use a JFNK solver at all, and is necessary here specifically, not
    // just generically "good practice" -- confirmed by direct test, not
    // assumed. Newton's convergence basin shrinks as the step gets more
    // nonlinear (larger dt): on a mild test case, a direct attempt at
    // dt = 1e-4 * stableTimeStep() converges to machine precision in 2
    // Newton iterations, while the SAME physical case at
    // dt = 0.3 * stableTimeStep() stalls after a large initial reduction
    // and never reaches even a modest tolerance in 30 iterations -- not
    // because of a bug (the residual evaluator and the matrix-free GMRES
    // linear solve were each independently verified correct), but because
    // the finite-difference-Jacobian Newton step's basin of convergence is
    // genuinely smaller than the full step at that stiffness. Repeated
    // smaller sub-steps, each well inside that basin, is the textbook
    // remedy (see e.g. Knoll & Keyes 2004's own discussion of continuation
    // strategies for exactly this failure mode), not a workaround for a
    // defect.
    struct ImplicitAdaptiveResult {
        bool converged = false;       // true iff the FULL requested dt was covered
        double dtCovered = 0.0;       // actual simulated time advanced (== dt if converged)
        int subSteps = 0;
        int totalNewtonIterations = 0;
        int totalGmresIterations = 0;
    };
    ImplicitAdaptiveResult stepImplicitPressureVelocityAdaptive(double dt);

    // Anderson-accelerated Picard (fixed-point) solve of the EXACT SAME
    // backward-Euler (u1, u2, P) residual system stepImplicitPressureVelocity()
    // solves by Newton-GMRES -- same scoping, same snapshot/lagging
    // conventions, same scaled unknown/residual (see
    // FourFieldSolverImplicit.cpp for the shared formulation). The
    // difference is purely the solution method: instead of a Newton step
    // with a matrix-free GMRES linear solve, each iteration takes the
    // trivial Richardson/Picard step xHat_{k+1} = xHat_k - Fhat(xHat_k)
    // (one scaled residual evaluation, no linear solve at all), optionally
    // extrapolated by Anderson mixing (Walker & Ni, 2011, "Anderson
    // acceleration for fixed-point iterations") over a short window of
    // past iterates -- the standard, much cheaper alternative to a fully
    // coupled Newton-Krylov solve for accelerating a segregated fixed-point
    // loop, surveyed alongside JFNK itself. This exists specifically to
    // answer the question JFNK's own measured cost raises: how much of its
    // stability benefit can a far cheaper fixed-point accelerator recover,
    // without ever forming a Jacobian-vector product or running GMRES.
    //
    // Like stepImplicitPressureVelocity(), a non-converged call leaves
    // state_ completely unchanged.
    struct SegregatedOptions {
        int maxIters = 200;
        double tol = 1.0e-6;       // ||Fhat||, relative to the initial residual norm
        int andersonDepth = 0;     // 0 = plain Picard/Richardson; >0 = Anderson(m)
        double andersonBeta = 1.0; // mixing damping, (0,1]; 1.0 = undamped
        double andersonReg = 1.0e-10; // Tikhonov reg., AS A FRACTION of the least-squares Gram matrix's own average diagonal magnitude (relative, not absolute -- see andersonLeastSquares())
    } segregated;

    struct SegregatedStepResult {
        bool converged = false;
        int iterations = 0;
        double initialResidualNorm = 0.0;
        double finalResidualNorm = 0.0;
    };
    SegregatedStepResult stepSegregatedAccelerated(double dt);

    double time() const { return time_; }
    const FlowState& state() const { return state_; }
    FlowState& mutableState() { return state_; }

    const PipeGeometry& geometry() const { return geometry_; }
    const FluidProperties& fluid() const { return fluid_; }

    // Total liquid / gas mass currently stored in the domain [kg per unit
    // pipe cross-section... actually kg, integrating volume fraction *
    // density * cell volume (A*dz)].
    double totalLiquidMass() const;
    double totalGasMass() const;
    double totalDropletMass() const;
    double totalBubbleMass() const;

    // Diagnostics for Eq. (24)-style mass-conservation bookkeeping.
    struct MassFluxes { double liquidIn, liquidOut, gasIn, gasOut; };
    MassFluxes boundaryMassFluxes() const;

    // Per-cell Kelvin-Helmholtz stability parameter F, using the CURRENT
    // state -- the same quantity used as the AMR refinement indicator when
    // SolverOptions::amr.enabled is set. Exposed regardless of whether AMR
    // is on, since it's a useful diagnostic (e.g. to plot alongside the
    // mesh) on its own. See Closures.hpp, kelvinHelmholtzParameterF().
    std::vector<double> refinementIndicatorProfile() const;

private:
    PipeGeometry geometry_;
    FluidProperties fluid_;
    SolverOptions options_;
    BoundaryConditions bc_;
    FlowState state_;
    double time_ = 0.0;
    double initialDz_ = 0.0; // reference spacing (L / initial N), for AMR min/max width bounds
    int initialN_ = 0;
    long stepCount_ = 0;

    // Per-cell cached geometry from the last computeGeometry() call.
    std::vector<PipeGeometry::StratifiedGeometry> geom_;

    // Per-cell reference liquid-height equilibrium profile (see
    // SolverOptions::enableWellBalancedGravity), rebuilt from the current
    // theta field every step inside computeGeometry() -- cheap (O(N)) and
    // consistent with this file's existing no-stale-caching style. Only
    // read when enableWellBalancedGravity is set.
    std::vector<double> h1Eq_;

    // Face-centred shear stresses (size N+1) and cell-centred entrainment /
    // disengagement / deposition rates (size N), refreshed each step by
    // computeClosures() from the previous step's velocity field (lagged
    // coefficients, standard for an explicit / semi-implicit scheme).
    std::vector<double> tauW1_, tauW2_, tauI_;
    std::vector<double> Ue_, Ud_, phiE_, phiDe_;
    std::vector<double> rhoGasCell_;

    // Time-dependent inlet volume fractions (includes the optional seed
    // disturbance on the liquid holdup).
    double inletEL() const;
    double inletEd() const;
    double inletEl() const;
    double inletEG() const;
    double inletEb() const;
    double inletEg() const;

    void computeGeometry();
    void computeH1Equilibrium();
    void computeClosures();
    void applyInletBoundary();
    void applyOutletBoundary();
    void backSubstitutePhaseVelocities();
    void updateLayerMomentum(double dt);
    void updateDispersedMomentum(double dt);
    void solvePressureCorrection(double dt);
    void updateContinuity(double dt);
    void clampVolumeFractions();

    // Adaptive mesh: computeIndicator() gives the per-cell F used by both
    // refinementIndicatorProfile() and adaptMesh() itself; adaptMesh()
    // rebuilds every state vector in lock-step onto a new face array,
    // splitting cells flagged for refinement and merging adjacent pairs
    // flagged for coarsening. Returns true if the mesh changed.
    std::vector<double> computeIndicator() const;
    bool adaptMesh();

    // Moving mesh: computeMonitorFunction() gives the per-cell M used to
    // equidistribute node positions; relocateMesh() computes the new
    // (relaxed) node positions and conservatively remaps every state
    // vector onto them, keeping N fixed. See FourFieldSolver.cpp for the
    // equidistribution/remap mechanics.
    std::vector<double> computeMonitorFunction() const;
    bool relocateMesh();

    // JFNK (see stepImplicitPressureVelocity() above and
    // FourFieldSolverImplicit.cpp for the full derivation and the residual
    // definition/unknown layout). Kept in a separate translation unit from
    // the rest of this file so the existing, heavily regression-tested
    // explicit/semi-implicit pipeline above is untouched by this addition --
    // implemented entirely inside stepImplicitPressureVelocity() itself
    // (private free functions/lambdas in that .cpp, not additional class
    // methods), reading the same private state (geom_, rhoGasCell_, etc.)
    // but never writing it except at the very end of a CONVERGED call.
};

} // namespace mfs
