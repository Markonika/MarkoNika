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
};

} // namespace mfs
