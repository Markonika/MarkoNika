#pragma once

#include "mfs/Closures.hpp"
#include "mfs/FlowState.hpp"
#include "mfs/FluidProperties.hpp"
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

private:
    PipeGeometry geometry_;
    FluidProperties fluid_;
    SolverOptions options_;
    BoundaryConditions bc_;
    FlowState state_;
    double time_ = 0.0;

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
};

} // namespace mfs
