#include "mfs/FourFieldSolver.hpp"
#include "mfs/Constants.hpp"
#include "mfs/TridiagonalSolver.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace mfs {

using constants::gravity;
using constants::pi;
using constants::small_e;
using constants::tiny;

FourFieldSolver::FourFieldSolver(double diameter, double length, int nCells,
                                  FluidProperties fluid, SolverOptions options)
    : geometry_(diameter), fluid_(fluid), options_(options) {
    state_.resize(nCells, length);
    initialDz_ = length / nCells;
    initialN_ = nCells;
    geom_.resize(nCells);
    tauW1_.assign(nCells + 1, 0.0);
    tauW2_.assign(nCells + 1, 0.0);
    tauI_.assign(nCells + 1, 0.0);
    Ue_.assign(nCells, 0.0);
    Ud_.assign(nCells, 0.0);
    phiE_.assign(nCells, 0.0);
    phiDe_.assign(nCells, 0.0);
    rhoGasCell_.assign(nCells, fluid_.rhoGas(bc_.outletPressure));
}

void FourFieldSolver::setInclinationProfile(const std::function<double(double)>& thetaOfZ) {
    for (int i = 0; i < state_.N; ++i) state_.theta[i] = thetaOfZ(state_.cellCenter(i));
}

void FourFieldSolver::setInclinationConstant(double thetaRad) {
    std::fill(state_.theta.begin(), state_.theta.end(), thetaRad);
}

void FourFieldSolver::setBoundaryConditions(const BoundaryConditions& bc) { bc_ = bc; }

double FourFieldSolver::inletEL() const {
    double eL = bc_.inletLiquidHoldup;
    if (bc_.seedDisturbance) {
        eL *= (1.0 + bc_.disturbanceAmplitude * std::sin(2.0 * pi * bc_.disturbanceFrequency * time_));
    }
    return std::clamp(eL, small_e, 1.0 - small_e);
}
double FourFieldSolver::inletEd() const { return inletEL() * bc_.inletDropletFractionOfLiquid; }
double FourFieldSolver::inletEl() const { return inletEL() - inletEd(); }
double FourFieldSolver::inletEG() const { return 1.0 - inletEL(); }
double FourFieldSolver::inletEb() const { return inletEG() * bc_.inletBubbleFractionOfGas; }
double FourFieldSolver::inletEg() const { return inletEG() - inletEb(); }

void FourFieldSolver::initializeStratified(double eL0) {
    const double eL = std::clamp(eL0, small_e, 1.0 - small_e);
    for (int i = 0; i < state_.N; ++i) {
        state_.el[i] = eL;
        state_.ed[i] = 0.0;
        state_.eg[i] = 1.0 - eL;
        state_.eb[i] = 0.0;
        state_.P[i] = bc_.outletPressure;
    }
    const double ul0 = bc_.inletSuperficialLiquid / eL;
    const double ug0 = bc_.inletSuperficialGas / (1.0 - eL);
    std::fill(state_.u1.begin(), state_.u1.end(), ul0);
    std::fill(state_.u2.begin(), state_.u2.end(), ug0);
    std::fill(state_.ud.begin(), state_.ud.end(), ul0);
    std::fill(state_.ub.begin(), state_.ub.end(), ug0);
    std::fill(state_.ul.begin(), state_.ul.end(), ul0);
    std::fill(state_.ug.begin(), state_.ug.end(), ug0);
    std::fill(rhoGasCell_.begin(), rhoGasCell_.end(), fluid_.rhoGas(bc_.outletPressure));
    time_ = 0.0;
    // geom_ is otherwise only refreshed inside step() (via computeGeometry()
    // at the top of the pipeline); populate it here too so a stableTimeStep()
    // call made before the first step() -- the conventional usage, e.g.
    // `dt = solver.stableTimeStep(); solver.step(dt);` -- sees valid
    // geometry if it needs it (the surface-tension stability cap does; see
    // SolverOptions::enableSurfaceTension).
    computeGeometry();
}

void FourFieldSolver::computeGeometry() {
    for (int i = 0; i < state_.N; ++i) {
        geom_[i] = geometry_.fromAreaFraction(state_.e1(i));
        rhoGasCell_[i] = fluid_.rhoGas(state_.P[i]);
    }
    if (options_.enableWellBalancedGravity) computeH1Equilibrium();
}

// See SolverOptions::enableWellBalancedGravity for the full derivation.
// Builds the per-cell-centre reference profile h1Eq satisfying
// dh1Eq/dz = -tan(theta(z)) exactly, using EACH cell's own theta (not a
// face-averaged value) -- a simple running sum from faceZ[0], anchored
// arbitrarily at 0 (only differences of h1Eq ever enter the scheme, so the
// anchor value itself is never read).
void FourFieldSolver::computeH1Equilibrium() {
    const int N = state_.N;
    h1Eq_.resize(N);
    double faceVal = 0.0;
    for (int i = 0; i < N; ++i) {
        const double tanTheta = std::tan(state_.theta[i]);
        h1Eq_[i] = faceVal - tanTheta * (state_.cellCenter(i) - state_.faceZ[i]);
        faceVal -= tanTheta * state_.cellWidth(i);
    }
}

double FourFieldSolver::stableTimeStep() const {
    double umax = tiny;
    for (double v : state_.u1) umax = std::max(umax, std::fabs(v));
    for (double v : state_.u2) umax = std::max(umax, std::fabs(v));
    for (double v : state_.ud) umax = std::max(umax, std::fabs(v));
    for (double v : state_.ub) umax = std::max(umax, std::fabs(v));
    double dt = options_.courantTarget * state_.minCellWidth() / umax;

    if (options_.enableSurfaceTension) {
        // Stability bound for the biharmonic hyperdiffusion term
        // (SolverOptions::enableSurfaceTension / hyperdiffusionCoefficient,
        // see FourFieldSolver.cpp's updateContinuity()). Unlike the
        // dispersive third-derivative term this replaced, a fourth
        // derivative's discrete symbol is real and non-negative everywhere
        // (von Neumann analysis of the standard centred stencil:
        // 16*sin^4(theta/2)/dz^4), so -nu4*D4 is genuinely diffusive and
        // explicit forward-Euler stability is the standard, textbook
        // bound dt <= dz^4/(8*nu4), evaluated per cell with its own local
        // nu4 = hyperdiffusionCoefficient*|ul|*dz^3, taking the minimum
        // (most restrictive) over the domain. safety keeps a margin below
        // that exact bound for the coupling with the rest of this
        // (nonlinear) system, which the bound itself doesn't account for.
        const double safety = 0.5;
        double dtCap = options_.maxTimeStep;
        for (int i = 0; i < state_.N; ++i) {
            const double dz = std::max(state_.cellWidth(i), tiny);
            const double ulLocal = 0.5 * (state_.ul[i] + state_.ul[i + 1]);
            const double nu4 = options_.hyperdiffusionCoefficient * std::fabs(ulLocal) * dz * dz * dz;
            const double bound = dz * dz * dz * dz / (8.0 * std::max(nu4, tiny));
            dtCap = std::min(dtCap, safety * bound);
        }
        dt = std::min(dt, dtCap);
    }

    if (options_.enableTurbulentViscosity) {
        // Stability bound for the turbulent-viscosity diffusion of the
        // relative layer velocity (SolverOptions::enableTurbulentViscosity,
        // updateLayerMomentum()). This is a genuine (not hyper-) diffusion
        // term, so its discrete symbol (von Neumann analysis of the same
        // centred 3-point stencil used by secondDerivativeAtFace(): a
        // real, non-positive -4*sin^2(theta/2)/dz^2) gives the standard
        // textbook explicit-diffusion bound dt <= dz^2/(2*nu_t), evaluated
        // per face with its own local nu_t = mixingLength*|u2-u1|, taking
        // the minimum (most restrictive) over the domain. The reduced-mass
        // split in updateLayerMomentum() only ever reduces each layer's
        // own effective diffusivity below nu_t, so using the unscaled nu_t
        // here is a conservative (safe) bound.
        const double safety = 0.5;
        const double mixingLength = options_.turbulentMixingLengthFraction * geometry_.diameter();
        double dtCap = options_.maxTimeStep;
        for (int f = 1; f < state_.N; ++f) {
            const double dz = std::max(std::min(state_.cellWidth(f - 1), state_.cellWidth(f)), tiny);
            const double urLocal = state_.u2[f] - state_.u1[f];
            const double nuT = mixingLength * std::fabs(urLocal);
            const double bound = dz * dz / (2.0 * std::max(nuT, tiny));
            dtCap = std::min(dtCap, safety * bound);
        }
        dt = std::min(dt, dtCap);
    }

    return std::clamp(dt, options_.minTimeStep, options_.maxTimeStep);
}

void FourFieldSolver::computeClosures() {
    const double D = geometry_.diameter();
    const double A = geometry_.area();
    const double rhoL = fluid_.rhoLiquid;
    const int N = state_.N;

    for (int f = 0; f <= N; ++f) {
        const int cL = (f == 0) ? 0 : f - 1;
        const int cR = (f == N) ? N - 1 : f;

        const double D1f = 0.5 * (geom_[cL].D1 + geom_[cR].D1);
        const double D2f = 0.5 * (geom_[cL].D2 + geom_[cR].D2);
        const double h1f = 0.5 * (geom_[cL].h1 + geom_[cR].h1);
        const double A1f = 0.5 * (geom_[cL].A1 + geom_[cR].A1);
        const double A2f = 0.5 * (geom_[cL].A2 + geom_[cR].A2);
        const double Sif = 0.5 * (geom_[cL].Si + geom_[cR].Si);
        const double rhoGf = 0.5 * (rhoGasCell_[cL] + rhoGasCell_[cR]);
        const double thetaF = 0.5 * (state_.theta[cL] + state_.theta[cR]);
        const double e2f = A2f / A;

        const double ulO = state_.ul[f];
        const double ugO = state_.ug[f];

        const double flw = wallFrictionFactor({rhoL, ulO, D1f, fluid_.muLiquid, options_.wallRoughness},
                                               options_.wallCorrelation);
        const double fgw = wallFrictionFactor({rhoGf, ugO, D2f, fluid_.muGas, options_.wallRoughness},
                                               options_.wallCorrelation);

        InterfacialFrictionInputs ifi{};
        ifi.densityGas = rhoGf;
        ifi.densityLiquid = rhoL;
        ifi.viscosityGas = fluid_.muGas;
        ifi.velocityGas = ugO;
        ifi.velocityLiquid = ulO;
        ifi.hydraulicDiameter1 = D1f;
        ifi.hydraulicDiameter2 = D2f;
        ifi.liquidHeight = h1f;
        ifi.pipeDiameter = D;
        ifi.area1 = A1f;
        ifi.area2 = A2f;
        ifi.areaFraction2 = e2f;
        ifi.dA1dh1 = Sif;
        ifi.inclination = thetaF;
        ifi.localPressure = state_.P[cL];
        ifi.gasWallFriction = fgw;
        // See SolverOptions::closureF0Shift: shifts the AndreussiPersen1987
        // branch's own F0 onset threshold, zero (true no-op) by default.
        // Set on `ifi` BEFORE the call (unlike closureCorrectionScale
        // below, which rescales the already-computed `fi`) because this
        // changes which regime `fi` itself falls into, not just its size.
        ifi.f0Shift = options_.enableDataDrivenClosureCorrection ? options_.closureF0Shift : 0.0;
        const double fi = interfacialFrictionFactor(ifi, options_.interfacialCorrelation);
        // See SolverOptions::enableDataDrivenClosureCorrection: rescales
        // ONLY the correlation's enhancement above the baseline gas-wall
        // friction factor (fi-fgw), a no-op wherever that enhancement is
        // already zero. Computed as a separate local rather than folded
        // into `fi` itself so the option-off path below is a single
        // untouched read of `fi`, not an arithmetic identity relying on
        // closureCorrectionScale's default happening to be neutral.
        const double fiUsed = options_.enableDataDrivenClosureCorrection
                                   ? fgw + options_.closureCorrectionScale * (fi - fgw)
                                   : fi;

        tauW1_[f] = wallShearStress(flw, rhoL, ulO);
        tauW2_[f] = wallShearStress(fgw, rhoGf, ugO);
        tauI_[f] = interfacialShearStress(fiUsed, rhoGf, ugO, ulO);
    }

    for (int i = 0; i < N; ++i) {
        const double ulC = 0.5 * (state_.ul[i] + state_.ul[i + 1]);
        const double ugC = 0.5 * (state_.ug[i] + state_.ug[i + 1]);
        const double udC = 0.5 * (state_.ud[i] + state_.ud[i + 1]);
        const double ubC = 0.5 * (state_.ub[i] + state_.ub[i + 1]);
        const double umix = state_.el[i] * ulC + state_.ed[i] * udC + state_.eg[i] * ugC + state_.eb[i] * ubC;

        WaveCelerityInputs wci{umix, D, state_.theta[i]};
        const double uwave = waveCelerity(wci);

        BubbleEntrainmentInputs bei{rhoGasCell_[i], geom_[i].Si, uwave, ulC, A};
        phiE_[i] = bubbleEntrainmentRate(bei, D);

        BubbleDisengagementInputs bdi{rhoGasCell_[i], rhoL, fluid_.sigma, geom_[i].Si, state_.el[i], 0.28};
        phiDe_[i] = bubbleDisengagementRate(bdi);

        DropletEntrainmentInputs dei{D, rhoGasCell_[i], rhoL, ugC, fluid_.muLiquid, 7.7e-8};
        Ue_[i] = dropletEntrainmentRate(dei);

        // See SolverOptions::closureDepositionVelocityScale: true no-op
        // (0.1*1.0 == 0.1 exactly) unless enabled.
        const double depositionVelocity = 0.1 *
            (options_.enableDataDrivenClosureCorrection ? options_.closureDepositionVelocityScale : 1.0);
        DropletDepositionInputs ddi{D, rhoL, state_.ed[i], state_.eg[i], depositionVelocity};
        Ud_[i] = dropletDepositionRate(ddi);
    }
}

namespace {
    // MUSCL/TVD-limited face value (Sweby 1984 "high resolution" flux
    // form): donor + 0.5*psi(r)*(accept-donor), where donor/accept are the
    // upwind/downwind cell values straddling the face and farUpwindVal is
    // the next cell further upwind of donor (or the inlet boundary value
    // standing in for it at the first interior face), used to form
    // r = (donor-farUpwind)/(accept-donor). Falls back to first-order
    // upwind (returns donorVal unchanged) when: the limiter is None, no
    // far-upwind value is available (one cell shy of the far boundary,
    // where the 3-point stencil doesn't exist), or the local gradient
    // (accept-donor) is ~0 (r would be ill-conditioned there, but the
    // reconstruction reduces to first order anyway since a limiter's whole
    // point is moot on a flat profile).
    inline double limitedFaceValue(double donorVal, double acceptVal, double farUpwindVal,
                                    bool haveFarUpwind, FluxLimiterType limiter) {
        if (!haveFarUpwind || limiter == FluxLimiterType::None) return donorVal;
        const double denom = acceptVal - donorVal;
        if (std::fabs(denom) < 1.0e-12) return donorVal;
        const double r = (donorVal - farUpwindVal) / denom;
        const double psi = fluxLimiterPsi(limiter, r);
        return donorVal + 0.5 * psi * denom;
    }

    inline double rho1Of(const FlowState& s, const std::vector<double>& rhoGasCell, double rhoL, int i) {
        const double e1 = s.e1(i);
        return (s.el[i] * rhoL + s.eb[i] * rhoGasCell[i]) / std::max(e1, small_e);
    }
    inline double rho2Of(const FlowState& s, const std::vector<double>& rhoGasCell, double rhoL, int i) {
        const double e2 = s.e2(i);
        return (s.eg[i] * rhoGasCell[i] + s.ed[i] * rhoL) / std::max(e2, small_e);
    }

    // For the linear ODE dU/dt = A*U + b with constant 2x2 matrix
    // A = [[A11,A12],[A21,A22]], returns (alpha0, alpha1, beta0, beta1)
    // such that the EXACT solution over a step dt is
    //   U(dt) = (alpha0*I + alpha1*A)*U0 + (beta0*I + beta1*A)*b.
    // exp(A*dt) = alpha0*I + alpha1*A and the forcing integral
    // integral_0^dt exp(A*s) ds = beta0*I + beta1*A are both standard 2x2
    // matrix FUNCTIONS of A, built from A's eigenvalues via the confluent
    // divided-difference (Cayley-Hamilton) representation f(A) = beta0*I
    // + beta1*A: for distinct eigenvalues lam1 != lam2, beta1 =
    // (f(lam1)-f(lam2))/(lam1-lam2), beta0 = f(lam1) - lam1*beta1; for a
    // repeated eigenvalue, the confluent limit beta1 = f'(lam),
    // beta0 = f(lam) - lam*f'(lam). Used by SolverOptions::enableETDFriction
    // (f = exp(.*dt) for alpha, f = phi1-type (exp(.*dt)-1)/. for beta)
    // and, via the same eigenvalues, by enableIMEXRKFriction's two stage
    // solves. Assumes (and both callers' own documentation proves, for
    // the specific physical drag matrices this is applied to) that A's
    // eigenvalues are always real -- the discriminant is clamped to >= 0
    // defensively rather than handling a complex case that should not
    // arise.
    inline void expAndPhiCoeffs2x2(double A11, double A12, double A21, double A22, double dt,
                                    double& alpha0, double& alpha1, double& beta0, double& beta1) {
        const double T = A11 + A22;
        const double D = A11 * A22 - A12 * A21;
        const double disc = std::max(T * T - 4.0 * D, 0.0);
        const double sq = std::sqrt(disc);
        const double lam1 = 0.5 * (T + sq);
        const double lam2 = 0.5 * (T - sq);

        // phi1-type function f(z) = (e^z - 1)/lam for z = lam*dt, with a
        // direct (not 0/0) Taylor-series formula near z=0 -- the true
        // removable singularity there, reached whenever a face has zero
        // net drag along one of A's eigen-directions (e.g. no wall
        // friction at all), not just a theoretical edge case.
        auto fPhi = [dt](double lam) {
            const double z = lam * dt;
            if (std::fabs(z) < 1.0e-4) return dt * (1.0 + z * (0.5 + z * (1.0 / 6.0 + z * (1.0 / 24.0))));
            return (std::exp(z) - 1.0) / lam;
        };
        // d(fPhi)/dlam, needed only in the repeated-eigenvalue branch
        // below; same near-zero series safeguard, this time against the
        // exact formula's own z->0 cancellation.
        auto fPhiDeriv = [dt](double lam) {
            const double z = lam * dt;
            if (std::fabs(z) < 1.0e-4) return dt * dt * (0.5 + z * (1.0 / 3.0 + z * (1.0 / 8.0)));
            return (dt * std::exp(z) * lam - (std::exp(z) - 1.0)) / (lam * lam);
        };

        if (std::fabs(lam1 - lam2) > 1.0e-9 * std::max(1.0, std::fabs(lam1) + std::fabs(lam2))) {
            const double e1 = std::exp(lam1 * dt), e2 = std::exp(lam2 * dt);
            alpha1 = (e1 - e2) / (lam1 - lam2);
            alpha0 = e1 - lam1 * alpha1;
            const double p1 = fPhi(lam1), p2 = fPhi(lam2);
            beta1 = (p1 - p2) / (lam1 - lam2);
            beta0 = p1 - lam1 * beta1;
        } else {
            const double lam = 0.5 * (lam1 + lam2);
            const double e = std::exp(lam * dt);
            alpha0 = e * (1.0 - lam * dt);
            alpha1 = dt * e;
            beta1 = fPhiDeriv(lam);
            beta0 = fPhi(lam) - lam * beta1;
        }
    }

    // Flux-form (Burgers-type) MUSCL/TVD reconstruction for a face-centred
    // velocity field's OWN self-advection term u*du/dz, used by
    // updateLayerMomentum() in place of plain non-conservative upwind
    // differencing when a limiter is selected. u*du/dz = d(u^2/2)/dz for
    // smooth (differentiable) u -- the same identity that makes Burgers'
    // equation's non-conservative and conservative forms equivalent -- so
    // this recasts the term as a genuine flux divergence, letting it reuse
    // limitedFaceValue() exactly as continuity's fluxes do, rather than
    // needing a separately-derived non-conservative-form limiter.
    //
    // uFace (size N+1) lives at the mesh's FACES; the natural "flux point"
    // for a field staggered that way is the CELL CENTRE between two
    // consecutive face samples (uFace[i], uFace[i+1] straddle cell i), so
    // this builds one reconstructed value uHat(i) per CELL from the local
    // 3-point upwind-biased stencil (donor/accept faces of cell i, plus
    // one more face further upwind for the limiter ratio -- unavailable
    // only at the one cell adjacent to whichever end is upwind, where it
    // falls back to first order, same convention as updateContinuity's
    // upwind()). F(i) = 0.5*uHat(i)^2 is then differenced over the
    // CENTRE-to-centre distance in updateLayerMomentum(), matching the
    // spacing already used there for dPdz/dh1dz.
    //
    // This is a genuinely different BASE discretization from the plain
    // non-conservative first-order scheme (it uses the local flow
    // direction and a flux difference rather than u1c times a one-sided
    // difference of u1 itself), not purely "the same first-order scheme
    // plus a limiter" -- which is why it is only substituted in when a
    // limiter is actually selected, leaving the default (None) path
    // untouched and bit-for-bit unchanged from before this was added.
    std::vector<double> buildAdvectiveFlux(const std::vector<double>& uFace, int N, FluxLimiterType limiter) {
        std::vector<double> F(N);
        for (int i = 0; i < N; ++i) {
            const double flowDir = uFace[i] + uFace[i + 1];
            double donor, accept, farVal;
            bool haveFar;
            if (flowDir >= 0.0) {
                donor = uFace[i];
                accept = uFace[i + 1];
                haveFar = (i - 1 >= 0);
                farVal = haveFar ? uFace[i - 1] : 0.0;
            } else {
                donor = uFace[i + 1];
                accept = uFace[i];
                haveFar = (i + 2 <= N);
                farVal = haveFar ? uFace[i + 2] : 0.0;
            }
            const double uHat = limitedFaceValue(donor, accept, farVal, haveFar, limiter);
            F[i] = 0.5 * uHat * uHat;
        }
        return F;
    }

    // Non-uniform-mesh second derivative of a cell-centred field at cell i
    // (3-point, standard formula: reduces to the familiar
    // (v[i+1]-2v[i]+v[i-1])/dz^2 on a uniform mesh). Used twice in
    // succession by fourthDerivativeEL() below to build a biharmonic
    // ("hyperdiffusion") operator -- see SolverOptions::enableSurfaceTension
    // for why a plain second derivative alone isn't what's wanted here (it
    // would damp long, physically meaningful wavelengths too; the repeated
    // application targets short wavelengths much more selectively, exactly
    // the k^4-vs-k^2 scaling real surface tension provides -- see
    // VALIDATION.md item 1, Update 5/6 for the derivation).
    double secondDerivativeAt(const FlowState& s, const std::vector<double>& v, int i) {
        const double dzL = s.centerDistance(i - 1, i);
        const double dzR = s.centerDistance(i, i + 1);
        return 2.0 * ((v[i + 1] - v[i]) / (dzR * (dzL + dzR)) - (v[i] - v[i - 1]) / (dzL * (dzL + dzR)));
    }

    // Biharmonic ("hyperdiffusion") operator on a cell-centred field,
    // D4 = D2[D2[v]] (apply secondDerivativeAt() twice). Needs v at
    // i-2..i+2, so valid only for i in [2, N-3]; returns 0 outside that
    // range, consistent with this file's other boundary fallbacks (see
    // e.g. updateContinuity()'s upwind() lambda).
    double fourthDerivativeAt(const FlowState& s, const std::vector<double>& v, int i) {
        const int N = s.N;
        if (i < 2 || i > N - 3) return 0.0;
        const double d2Lo = secondDerivativeAt(s, v, i - 1);
        const double d2Mid = secondDerivativeAt(s, v, i);
        const double d2Hi = secondDerivativeAt(s, v, i + 1);
        const double dzL = s.centerDistance(i - 1, i);
        const double dzR = s.centerDistance(i, i + 1);
        return 2.0 * ((d2Hi - d2Mid) / (dzR * (dzL + dzR)) - (d2Mid - d2Lo) / (dzL * (dzL + dzR)));
    }

    // Same 3-point non-uniform second-derivative formula as
    // secondDerivativeAt() above, but for a FACE-indexed field (u1, u2 and
    // combinations thereof are stored one value per face, index 0..N, not
    // one per cell). The spacing between consecutive faces f-1,f and f,f+1
    // is exactly the width of the cell between them (cellWidth(f-1) and
    // cellWidth(f) respectively) rather than a center-to-center distance.
    // Used by SolverOptions::enableTurbulentViscosity's diffusive term on
    // the relative layer velocity (u2-u1); see updateLayerMomentum().
    // Valid for f in [1, N-1]; returns 0 outside that range.
    double secondDerivativeAtFace(const FlowState& s, const std::vector<double>& v, int f) {
        const int N = s.N;
        if (f < 1 || f > N - 1) return 0.0;
        const double dzL = s.cellWidth(f - 1);
        const double dzR = s.cellWidth(f);
        return 2.0 * ((v[f + 1] - v[f]) / (dzR * (dzL + dzR)) - (v[f] - v[f - 1]) / (dzL * (dzL + dzR)));
    }
}

void FourFieldSolver::updateLayerMomentum(double dt) {
    const double A = geometry_.area();
    const double rhoL = fluid_.rhoLiquid;
    const int N = state_.N;

    std::vector<double> u1New = state_.u1;
    std::vector<double> u2New = state_.u2;

    // See buildAdvectiveFlux()'s comment: only switches to the flux-form
    // TVD scheme when a limiter is actually selected, otherwise the loop
    // below keeps using the original plain upwind differencing unchanged.
    const bool useFluxForm = options_.advectionLimiter != FluxLimiterType::None;
    std::vector<double> F1flux, F2flux;
    if (useFluxForm) {
        F1flux = buildAdvectiveFlux(state_.u1, N, options_.advectionLimiter);
        F2flux = buildAdvectiveFlux(state_.u2, N, options_.advectionLimiter);
    }

    // See SolverOptions::enableTurbulentViscosity: relative-velocity field
    // (u2-u1) at the OLD time level, built once before the face loop so
    // secondDerivativeAtFace() can read neighbouring faces regardless of
    // loop order.
    std::vector<double> urOld;
    const double mixingLength = options_.turbulentMixingLengthFraction * geometry_.diameter();
    if (options_.enableTurbulentViscosity) {
        urOld.resize(N + 1);
        for (int i = 0; i <= N; ++i) urOld[i] = state_.u2[i] - state_.u1[i];
    }

    for (int f = 1; f < N; ++f) {
        const int cL = f - 1, cR = f;
        const double e1L = state_.e1(cL), e1R = state_.e1(cR);
        const double e2L = state_.e2(cL), e2R = state_.e2(cR);
        const double e1f = 0.5 * (e1L + e1R);
        const double e2f = 0.5 * (e2L + e2R);
        // See SolverOptions::momentumFractionFloor: used only where e1f/e2f
        // divide a momentum-equation source term, not for the terms
        // themselves or for the tracked volume fractions.
        const double e1fSafe = std::max(e1f, options_.momentumFractionFloor);
        const double e2fSafe = std::max(e2f, options_.momentumFractionFloor);
        const double rho1f = 0.5 * (rho1Of(state_, rhoGasCell_, rhoL, cL) + rho1Of(state_, rhoGasCell_, rhoL, cR));
        const double rho2f = 0.5 * (rho2Of(state_, rhoGasCell_, rhoL, cL) + rho2Of(state_, rhoGasCell_, rhoL, cR));

        // Cell-centred gradients (P, h1) use the distance between the two
        // neighbouring cell CENTRES, which equals a shared dz only on a
        // uniform mesh.
        const double centerDz = state_.centerDistance(cL, cR);
        const double dPdz = (state_.P[cR] - state_.P[cL]) / centerDz;
        const double thetaF = 0.5 * (state_.theta[cL] + state_.theta[cR]);
        const double dh1dz = (geom_[cR].h1 - geom_[cL].h1) / centerDz;

        // See SolverOptions::enableWellBalancedGravity. gravCosTerm and
        // gravSinTerm are extracted as named sub-expressions ONLY -- when
        // the option is off, each is computed via the exact same operation
        // order as the original inline "gravity * std::cos(thetaF) *
        // dh1dz" / "gravity * std::sin(thetaF)" expressions below, so
        // "- gravCosTerm - gravSinTerm" in the sums further down is
        // bit-for-bit identical to the pre-existing inline arithmetic (see
        // the enableImplicitFriction comment above on why term order is
        // preserved deliberately, not just "equivalent").
        double gravCosTerm, gravSinTerm;
        if (!options_.enableWellBalancedGravity) {
            gravCosTerm = gravity * std::cos(thetaF) * dh1dz;
            gravSinTerm = gravity * std::sin(thetaF);
        } else {
            const double detadz = ((geom_[cR].h1 - h1Eq_[cR]) - (geom_[cL].h1 - h1Eq_[cL])) / centerDz;
            gravCosTerm = gravity * std::cos(thetaF) * detadz;
            gravSinTerm = 0.0; // folded into detadz already -- see derivation above
        }
        const double Swp1f = 0.5 * (geom_[cL].Swp1 + geom_[cR].Swp1);
        const double Swp2f = 0.5 * (geom_[cL].Swp2 + geom_[cR].Swp2);
        const double Sif = 0.5 * (geom_[cL].Si + geom_[cR].Si);

        const double UeF = 0.5 * (Ue_[cL] + Ue_[cR]);
        const double UdF = 0.5 * (Ud_[cL] + Ud_[cR]);
        const double phiEF = 0.5 * (phiE_[cL] + phiE_[cR]);
        const double phiDeF = 0.5 * (phiDe_[cL] + phiDe_[cR]);

        const double ulF = state_.ul[f], ugF = state_.ug[f], udF = state_.ud[f], ubF = state_.ub[f];

        // See SolverOptions::enableTurbulentViscosity for the derivation:
        // diffusive acceleration on the relative velocity a_diff =
        // nu_t*d^2(u2-u1)/dz^2, split between the two layers by
        // reduced-mass weighting so total momentum is conserved exactly
        // while d(u2-u1)/dt still receives the full a_diff.
        double turb1 = 0.0, turb2 = 0.0;
        if (options_.enableTurbulentViscosity) {
            const double nuT = mixingLength * std::fabs(urOld[f]);
            const double aDiff = nuT * secondDerivativeAtFace(state_, urOld, f);
            const double m1 = rho1f * e1fSafe;
            const double m2 = rho2f * e2fSafe;
            const double mSum = m1 + m2;
            turb1 = -aDiff * m2 / mSum;
            turb2 = +aDiff * m1 / mSum;
        }

        // --- Layer 1 momentum (Eq. 6, slip flux term dropped per Appendix A) ---
        // Upwind advection differences two FACE values (f against its
        // upwind neighbour f-1 or f+1); the correct spacing is the width of
        // the CELL between those two faces (cL=f-1 when going backward,
        // cR=f when going forward), not the center-to-center distance above.
        // With a limiter selected, this switches to the flux-form
        // d(u1^2/2)/dz built by buildAdvectiveFlux() instead -- see that
        // function's comment.
        double advectiveTerm1;
        if (useFluxForm) {
            advectiveTerm1 = (F1flux[cR] - F1flux[cL]) / centerDz;
        } else {
            const double u1c = state_.u1[f];
            const double du1dzAdv = (u1c >= 0.0) ? (state_.u1[f] - state_.u1[f - 1]) / state_.cellWidth(cL)
                                                  : (state_.u1[f + 1] - state_.u1[f]) / state_.cellWidth(cR);
            advectiveTerm1 = u1c * du1dzAdv;
        }
        const double massSrc1 = (-UeF * ulF + UdF * udF + phiEF * ugF - phiDeF * ubF) / (e1fSafe * rho1f);
        const double massSrc2 = (UeF * ulF - UdF * udF - phiEF * ugF + phiDeF * ubF) / (e2fSafe * rho2f);

        // --- Layer 2 advection (layer 1's is computed just above) ---
        double advectiveTerm2;
        if (useFluxForm) {
            advectiveTerm2 = (F2flux[cR] - F2flux[cL]) / centerDz;
        } else {
            const double u2c = state_.u2[f];
            const double du2dzAdv = (u2c >= 0.0) ? (state_.u2[f] - state_.u2[f - 1]) / state_.cellWidth(cL)
                                                  : (state_.u2[f + 1] - state_.u2[f]) / state_.cellWidth(cR);
            advectiveTerm2 = u2c * du2dzAdv;
        }

        if (options_.enableIMEXRKFriction || options_.enableETDFriction) {
            // Shared with enableImplicitFriction's own derivation (see that
            // option's documentation): E1, E2 are everything EXCEPT wall
            // and interfacial friction; freezing the friction factor and
            // the |velocity| that makes each shear stress quadratic at
            // the OLD state turns the two friction terms into a LINEAR
            // drag-rate matrix A = [[-(kw1+ki1), ki1], [ki2, -(kw2+ki2)]]
            // acting on U=(u1,u2), with constant forcing b=(E1,E2) -- the
            // SAME linear ODE dU/dt = A*U + b that enableImplicitFriction
            // solves by backward Euler, but integrated differently below.
            // Computed fresh here (not shared with the branches above/
            // below) deliberately -- same self-containment rationale as
            // rho1Of/rho2Of's own duplication across translation units in
            // this codebase: it keeps the two already bit-exact-verified
            // branches below completely untouched by this addition.
            const double E1 = -advectiveTerm1 - dPdz / rho1f - gravCosTerm - gravSinTerm + massSrc1 + turb1;
            const double E2 = -advectiveTerm2 - dPdz / rho2f - gravCosTerm - gravSinTerm + massSrc2 + turb2;
            const double epsVel = 1.0e-6;
            const double u1Old = state_.u1[f], u2Old = state_.u2[f];
            const double urMag = std::max(std::fabs(u2Old - u1Old), epsVel);
            const double kw1 = std::fabs(tauW1_[f] * Swp1f / (A * e1fSafe * rho1f)) / std::max(std::fabs(u1Old), epsVel);
            const double kw2 = std::fabs(tauW2_[f] * Swp2f / (A * e2fSafe * rho2f)) / std::max(std::fabs(u2Old), epsVel);
            const double ki1 = std::fabs(tauI_[f] * Sif / (A * e1fSafe * rho1f)) / urMag;
            const double ki2 = std::fabs(tauI_[f] * Sif / (A * e2fSafe * rho2f)) / urMag;
            const double lA11 = -(kw1 + ki1), lA12 = ki1, lA21 = ki2, lA22 = -(kw2 + ki2);

            if (options_.enableIMEXRKFriction) {
                // Pareschi & Russo (2005) IMEX-SSP2(2,2,2): 2-stage,
                // 2nd-order, L-stable (gamma = 1 - 1/sqrt(2)) -- see
                // SolverOptions::enableIMEXRKFriction for the full
                // derivation and explicit scoping note (E1, E2 frozen
                // across both stages, not re-evaluated at U1). Butcher
                // tableau: c=(gamma,1), A_im=[[gamma,0],[1-gamma,gamma]],
                // b=(1-gamma,gamma) -- the standard 2-stage, 2nd-order,
                // L-stable SDIRK, confirmed L-stable here (not assumed)
                // by deriving its scalar stability function
                // R(z)=(1-gamma^2*z)/(1-gamma*z)^2 and checking R(z)->0
                // as z=dt*lambda->-infinity. An EARLIER version of this
                // code used a21=(1-2*gamma) and b=(1/2,1/2) instead --
                // algebraically plausible-looking but WRONG: its own
                // stability function works out to R(infinity)=-1, not
                // L-stable at all, and it was caught exactly the way this
                // project catches this class of error -- not by
                // re-deriving it on paper first, but because the direct
                // stiff dt-sweep test below showed it blowing up at the
                // same dt where plain backward-Euler (enableImplicitFriction)
                // remains stable, the opposite of what a MORE accurate
                // integrator of the same stiff term should ever do.
                const double gamma = 1.0 - 1.0 / std::sqrt(2.0);
                // Stage matrix M = I - dt*gamma*A, shared by both stages.
                const double m11 = 1.0 - dt * gamma * lA11, m12 = -dt * gamma * lA12;
                const double m21 = -dt * gamma * lA21, m22 = 1.0 - dt * gamma * lA22;
                const double detM = m11 * m22 - m12 * m21;
                // Stage 1: M*U1 = U0.
                const double u1_1 = (u1Old * m22 - m12 * u2Old) / detM;
                const double u2_1 = (m11 * u2Old - m21 * u1Old) / detM;
                // Stage 2: M*U2 = U0 + dt*b + dt*(1-gamma)*A*U1.
                const double rhs1 = u1Old + dt * E1 + dt * (1.0 - gamma) * (lA11 * u1_1 + lA12 * u2_1);
                const double rhs2 = u2Old + dt * E2 + dt * (1.0 - gamma) * (lA21 * u1_1 + lA22 * u2_1);
                const double u1_2 = (rhs1 * m22 - m12 * rhs2) / detM;
                const double u2_2 = (m11 * rhs2 - m21 * rhs1) / detM;
                // Final combination: U_new = U0 + dt*b + dt*[(1-gamma)*A*U1 + gamma*A*U2].
                const double Au1 = lA11 * u1_1 + lA12 * u2_1;
                const double Au2 = lA21 * u1_1 + lA22 * u2_1;
                const double Av1 = lA11 * u1_2 + lA12 * u2_2;
                const double Av2 = lA21 * u1_2 + lA22 * u2_2;
                u1New[f] = u1Old + dt * E1 + dt * ((1.0 - gamma) * Au1 + gamma * Av1);
                u2New[f] = u2Old + dt * E2 + dt * ((1.0 - gamma) * Au2 + gamma * Av2);
            } else {
                // Exponential time differencing: integrate dU/dt=A*U+b
                // EXACTLY (zero error beyond the linearization itself) --
                // see SolverOptions::enableETDFriction.
                double alpha0, alpha1, beta0, beta1;
                expAndPhiCoeffs2x2(lA11, lA12, lA21, lA22, dt, alpha0, alpha1, beta0, beta1);
                const double AU1 = lA11 * u1Old + lA12 * u2Old;
                const double AU2 = lA21 * u1Old + lA22 * u2Old;
                const double Ab1 = lA11 * E1 + lA12 * E2;
                const double Ab2 = lA21 * E1 + lA22 * E2;
                u1New[f] = alpha0 * u1Old + alpha1 * AU1 + beta0 * E1 + beta1 * Ab1;
                u2New[f] = alpha0 * u2Old + alpha1 * AU2 + beta0 * E2 + beta1 * Ab2;
            }
        } else if (!options_.enableImplicitFriction) {
            // Exact same left-to-right term order as before this option
            // existed (verified bit-for-bit reproducible against the
            // pre-existing benchmark) -- floating-point addition is not
            // strictly associative, so this is deliberately NOT rebuilt
            // from the E1/E2 decomposition below, which reorders terms
            // and was found to (harmlessly, but needlessly) perturb exact
            // step counts on chaotic-sensitive long runs.
            const double du1dt = -advectiveTerm1
                                  - dPdz / rho1f
                                  - gravCosTerm
                                  - gravSinTerm
                                  - tauW1_[f] * Swp1f / (A * e1fSafe * rho1f)
                                  + tauI_[f] * Sif / (A * e1fSafe * rho1f)
                                  + massSrc1
                                  + turb1;
            const double du2dt = -advectiveTerm2
                                  - dPdz / rho2f
                                  - gravCosTerm
                                  - gravSinTerm
                                  - tauW2_[f] * Swp2f / (A * e2fSafe * rho2f)
                                  - tauI_[f] * Sif / (A * e2fSafe * rho2f)
                                  + massSrc2
                                  + turb2;
            u1New[f] = state_.u1[f] + dt * du1dt;
            u2New[f] = state_.u2[f] + dt * du2dt;
        } else {
            // E1, E2: everything EXCEPT wall and interfacial friction,
            // used only by the implicit solve below (a separate code
            // path with no prior bit-exact baseline, so term order here
            // is free to differ from the explicit branch above).
            const double E1 = -advectiveTerm1
                               - dPdz / rho1f
                               - gravCosTerm
                               - gravSinTerm
                               + massSrc1
                               + turb1;
            const double E2 = -advectiveTerm2
                               - dPdz / rho2f
                               - gravCosTerm
                               - gravSinTerm
                               + massSrc2
                               + turb2;
            // Local 2x2 implicit (backward-Euler) solve for wall and
            // interfacial friction, all other terms (E1, E2) explicit.
            // Linearise each nonlinear friction force by extracting an
            // effective linear drag rate k = (old acceleration)/(old
            // velocity or relative velocity) -- i.e. freeze the friction
            // factor and the |velocity| that made the original closure
            // quadratic, at the OLD state, and solve for the new
            // velocities exactly given those frozen coefficients. See
            // SolverOptions::enableImplicitFriction for the full
            // derivation.
            // wallShearStress()/interfacialShearStress() return a SIGNED
            // force (proportional to |u|*u or |u_r|*u_r, not |u| alone),
            // matching the sign of the velocity/relative-velocity it
            // opposes. The linear coefficient k must therefore be built
            // from |force|/|velocity| (always >= 0), NOT signed-force
            // divided by |velocity| -- the latter would silently flip
            // sign whenever the old velocity was negative, turning
            // friction into anti-friction. Caught by inspecting
            // wallShearStress()'s actual implementation before relying on
            // this, not assumed.
            const double epsVel = 1.0e-6;
            const double u1Old = state_.u1[f], u2Old = state_.u2[f];
            const double urMag = std::max(std::fabs(u2Old - u1Old), epsVel);
            const double kw1 = std::fabs(tauW1_[f] * Swp1f / (A * e1fSafe * rho1f)) / std::max(std::fabs(u1Old), epsVel);
            const double kw2 = std::fabs(tauW2_[f] * Swp2f / (A * e2fSafe * rho2f)) / std::max(std::fabs(u2Old), epsVel);
            const double ki1 = std::fabs(tauI_[f] * Sif / (A * e1fSafe * rho1f)) / urMag;
            const double ki2 = std::fabs(tauI_[f] * Sif / (A * e2fSafe * rho2f)) / urMag;

            const double a11 = 1.0 + dt * kw1 + dt * ki1;
            const double a12 = -dt * ki1;
            const double a21 = -dt * ki2;
            const double a22 = 1.0 + dt * kw2 + dt * ki2;
            const double b1 = u1Old + dt * E1;
            const double b2 = u2Old + dt * E2;
            const double det = a11 * a22 - a12 * a21; // always >= 1 here: a11,a22>=1, a12*a21<=0
            u1New[f] = (b1 * a22 - a12 * b2) / det;
            u2New[f] = (a11 * b2 - a21 * b1) / det;
        }
    }

    const double vmax = options_.maxVelocity;
    for (double& v : u1New) v = std::clamp(v, -vmax, vmax);
    for (double& v : u2New) v = std::clamp(v, -vmax, vmax);
    state_.u1 = u1New;
    state_.u2 = u2New;
    applyInletBoundary();
    applyOutletBoundary();
}

void FourFieldSolver::backSubstitutePhaseVelocities() {
    // Eq. (3): u1 = (cb*rhoG*ub + (1-cb)*rhoL*ul)/rho1  =>  ul = (rho1*u1 - cb*rhoG*ub) / ((1-cb)*rhoL)
    // Eq. (4): u2 = (cd*rhoL*ud + (1-cd)*rhoG*ug)/rho2  =>  ug = (rho2*u2 - cd*rhoL*ud) / ((1-cd)*rhoG)
    const double rhoL = fluid_.rhoLiquid;
    const int N = state_.N;
    for (int f = 0; f <= N; ++f) {
        const int cL = (f == 0) ? 0 : f - 1;
        const int cR = (f == N) ? N - 1 : f;
        const double el = 0.5 * (state_.el[cL] + state_.el[cR]);
        const double eb = 0.5 * (state_.eb[cL] + state_.eb[cR]);
        const double eg = 0.5 * (state_.eg[cL] + state_.eg[cR]);
        const double ed = 0.5 * (state_.ed[cL] + state_.ed[cR]);
        const double e1 = std::max(el + eb, small_e);
        const double e2 = std::max(eg + ed, small_e);
        const double rhoG = 0.5 * (rhoGasCell_[cL] + rhoGasCell_[cR]);
        const double rho1 = (el * rhoL + eb * rhoG) / e1;
        const double rho2 = (eg * rhoG + ed * rhoL) / e2;
        const double cb = std::clamp(eb / e1, 0.0, 1.0);
        const double cd = std::clamp(ed / e2, 0.0, 1.0);

        // Floored at the coarser momentum-fraction scale (not small_e):
        // (1-cb)/(1-cd) appear here as a genuine denominator amplifying
        // u1/u2 into ul/ug, so an overly fine floor lets a cell that is
        // almost entirely dispersed phase (cb or cd near 1) blow the
        // back-substituted continuous-phase velocity up well past what the
        // velocity limiter alone comfortably absorbs.
        const double floor = options_.momentumFractionFloor;
        const double denomL = std::max((1.0 - cb) * rhoL, floor * rhoL);
        state_.ul[f] = (rho1 * state_.u1[f] - cb * rhoG * state_.ub[f]) / denomL;

        const double denomG = std::max((1.0 - cd) * rhoG, floor * rhoG);
        state_.ug[f] = (rho2 * state_.u2[f] - cd * rhoL * state_.ud[f]) / denomG;
    }
}

void FourFieldSolver::updateDispersedMomentum(double /*dt*/) {
    // Eq. (10) is reduced to an algebraic balance (inertial terms
    // negligible, per the paper's own a posteriori check) so no explicit
    // time-step dependence appears here; dt is kept in the signature for
    // interface symmetry with the other update*() steps.
    const double D = geometry_.diameter();
    const double rhoL = fluid_.rhoLiquid;
    const int N = state_.N;

    std::vector<double> udNew = state_.ud;
    std::vector<double> ubNew = state_.ub;

    for (int f = 1; f < N; ++f) {
        const int cL = f - 1, cR = f;
        const double edF = 0.5 * (state_.ed[cL] + state_.ed[cR]);
        const double elF = 0.5 * (state_.el[cL] + state_.el[cR]);
        const double ebF = 0.5 * (state_.eb[cL] + state_.eb[cR]);
        const double rhoGf = 0.5 * (rhoGasCell_[cL] + rhoGasCell_[cR]);
        const double dPdz = (state_.P[cR] - state_.P[cL]) / state_.centerDistance(cL, cR);
        const double thetaF = 0.5 * (state_.theta[cL] + state_.theta[cR]);
        const double XeDrop = 0.5 * (Ue_[cL] + Ue_[cR]);
        const double XeBub = 0.5 * (phiE_[cL] + phiE_[cR]);
        // Note: the disengagement/deposition rates (Ud_, phiDe_) do NOT
        // appear below. Writing the algebraic balance for the VELOCITY
        // (rather than the momentum em*qm*um) via the quotient rule and
        // substituting the field's own continuity equation shows that a
        // term removing mass at the population's own mean velocity leaves
        // that mean velocity unchanged -- only entrainment (which adds
        // mass at a *different* source velocity, uk) and drag change um.
        // Concretely, starting from Eq. (10) and Eq. (12)/(14):
        //   d(um)/dt = [-em*dP/dz - em*qm*g*sin(theta) + Xe*(uk-um) + Fdrag(uc-um)] / (em*qm)
        // so the well-posed algebraic (steady-relaxation) closure is
        //   um = (-em*dP/dz - em*qm*g*sin(theta) + Xe*uk + Cdrag*uc) / (Xe + Cdrag)
        // This also fixes what would otherwise be a spurious 0/0 blow-up
        // whenever the dispersed fraction is (locally) zero but active
        // entrainment is feeding it: with Xe>0 the denominator stays
        // bounded away from zero and um relaxes to the source velocity uk,
        // exactly the correct physical limit for a nascent droplet/bubble
        // population.

        const double ul = state_.ul[f], ug = state_.ug[f];
        const double udOld = state_.ud[f], ubOld = state_.ub[f];

        // --- Droplets: continuous medium is gas ---
        {
            const double ustarG = std::sqrt(std::fabs(tauW2_[f]) / std::max(rhoGf, tiny));
            const double dd = dropletDiameter({std::max(ustarG, 1.0e-3), rhoGf, rhoL, fluid_.sigma, D});
            const double Red = reynoldsNumber(rhoGf, ug - udOld, dd, fluid_.muGas);
            const double CDd = dropletDragCoefficient(Red);
            const double Cdrag = 0.75 * (CDd / std::max(dd, tiny)) * rhoGf * edF * std::fabs(ug - udOld);
            const double S = -edF * dPdz - edF * rhoL * gravity * std::sin(thetaF);
            const double denom = XeDrop + Cdrag + tiny;
            udNew[f] = (S + XeDrop * ul + Cdrag * ug) / denom;
        }

        // --- Bubbles: continuous medium is liquid ---
        {
            const double flw = tauW1_[f] / std::max(0.5 * rhoL * ul * ul, tiny); // recover f_lw from stored stress
            const double u1c = state_.u1[f];
            const double db = bubbleDiameter({std::fabs(flw) > tiny ? std::fabs(flw) : 0.02, rhoL, u1c, elF, 1.05},
                                              fluid_.sigma);
            const double Reb = reynoldsNumber(rhoL, u1c - ubOld, db, fluid_.muLiquid);
            const double Eo = eotvosNumber(rhoL, rhoGf, db, fluid_.sigma);
            const double CDb = bubbleDragCoefficient(Reb, Eo);
            const double Cdrag = 0.75 * (CDb / std::max(db, tiny)) * rhoL * ebF * std::fabs(ul - ubOld);
            const double S = -ebF * dPdz - ebF * rhoGf * gravity * std::sin(thetaF);
            const double denom = XeBub + Cdrag + tiny;
            ubNew[f] = (S + XeBub * ug + Cdrag * ul) / denom;
        }
    }

    const double vmax = options_.maxVelocity;
    for (double& v : udNew) v = std::clamp(v, -vmax, vmax);
    for (double& v : ubNew) v = std::clamp(v, -vmax, vmax);
    state_.ud = udNew;
    state_.ub = ubNew;
    applyInletBoundary();
    applyOutletBoundary();
}

void FourFieldSolver::solvePressureCorrection(double dt) {
    const int N = state_.N;
    const double rhoL = fluid_.rhoLiquid;

    std::vector<double> a(N, 0.0), b(N, 0.0), c(N, 0.0), R(N, 0.0);

    auto rho1At = [&](int i) { return rho1Of(state_, rhoGasCell_, rhoL, i); };
    auto rho2At = [&](int i) { return rho2Of(state_, rhoGasCell_, rhoL, i); };

    // Interior faces: coupling between cell f-1 and cell f.
    for (int f = 1; f < N; ++f) {
        const int cL = f - 1, cR = f;
        const double rho1f = 0.5 * (rho1At(cL) + rho1At(cR));
        const double rho2f = 0.5 * (rho2At(cL) + rho2At(cR));

        const double elF = (state_.ul[f] >= 0.0) ? state_.el[cL] : state_.el[cR];
        const double edF = (state_.ud[f] >= 0.0) ? state_.ed[cL] : state_.ed[cR];
        const double egF = (state_.ug[f] >= 0.0) ? state_.eg[cL] : state_.eg[cR];
        const double ebF = (state_.ub[f] >= 0.0) ? state_.eb[cL] : state_.eb[cR];

        const double qLiquid = elF * state_.ul[f] + edF * state_.ud[f];
        const double qGas = egF * state_.ug[f] + ebF * state_.ub[f];
        const double qTotal = qLiquid + qGas;

        // Each cell's row normalises the flux divergence by ITS OWN width
        // (a proper finite-volume balance); the pressure-correction
        // sensitivity Gf normalises by the distance between the two cell
        // CENTRES instead (it comes from a gradient, see
        // updateLayerMomentum's dPdz for the same distinction). The two
        // coincide only on a uniform mesh.
        const double widthL = state_.cellWidth(cL);
        const double widthR = state_.cellWidth(cR);
        const double centerDz = state_.centerDistance(cL, cR);

        R[cL] += qTotal / widthL;
        R[cR] -= qTotal / widthR;

        const double Gf = (dt * (elF + ebF) / rho1f + dt * (edF + egF) / rho2f) / centerDz;
        b[cL] += Gf / widthL;
        c[cL] += -Gf / widthL;
        a[cR] += -Gf / widthR;
        b[cR] += Gf / widthR;
    }

    // Inlet face (fixed flow, no pressure sensitivity): contributes to R[0] only.
    {
        const double qLiquidIn = inletEl() * state_.ul[0] + inletEd() * state_.ud[0];
        const double qGasIn = inletEg() * state_.ug[0] + inletEb() * state_.ub[0];
        R[0] -= (qLiquidIn + qGasIn) / state_.cellWidth(0);
    }

    // Outlet face (fixed pressure, P'=0 there): contributes to R[N-1] and its own diagonal.
    // The "ghost" point beyond the outlet is taken one full last-cell-width
    // past the last cell centre (matching the velocity-correction loop
    // below and preserving the original uniform-mesh convention exactly
    // when the mesh happens to be uniform), not a half-width to the actual
    // outlet face -- a standard, simple boundary treatment, not something
    // this refactor changes.
    {
        const int c = N - 1;
        const double rho1f = rho1At(c);
        const double rho2f = rho2At(c);
        const double elF = state_.el[c], edF = state_.ed[c], egF = state_.eg[c], ebF = state_.eb[c];
        const double qLiquidOut = elF * state_.ul[N] + edF * state_.ud[N];
        const double qGasOut = egF * state_.ug[N] + ebF * state_.ub[N];
        const double widthC = state_.cellWidth(c);
        R[c] += (qLiquidOut + qGasOut) / widthC;

        const double Gf = (dt * (elF + ebF) / rho1f + dt * (edF + egF) / rho2f) / widthC;
        b[c] += Gf / widthC;
    }

    // Local gas-compressibility accumulation term.
    for (int i = 0; i < N; ++i) {
        const double eG = state_.eG(i);
        const double coeff = eG * fluid_.drhoGasdP() / (std::max(rhoGasCell_[i], tiny) * dt);
        b[i] += coeff;
    }

    std::vector<double> rhs(N);
    for (int i = 0; i < N; ++i) rhs[i] = -R[i];

    std::vector<double> Pprime = solveTridiagonal(a, b, c, rhs);
    // Under-relax the correction (see SolverOptions::pressureRelaxation):
    // a cell whose gas volume fraction is momentarily near zero (a thin
    // gas layer being pinched off under a forming slug) gives its
    // accumulation coefficient almost no weight, which can otherwise let
    // a single step's correction overshoot to a nonphysical pressure. On
    // top of relaxation, hard-clamp the fractional change per step: in the
    // same near-singular cells the linear system can still be marginally
    // conditioned enough to produce a large *relaxed* correction, and
    // without a magnitude cap that can compound step over step into
    // unbounded growth (in either sign) well before it is ever damped out.
    for (int i = 0; i < N; ++i) {
        Pprime[i] *= options_.pressureRelaxation;
        const double limit = options_.maxPressureChangeFraction * state_.P[i];
        Pprime[i] = std::clamp(Pprime[i], -limit, limit);
    }

    const double pMax = options_.maxPressureFactor * bc_.outletPressure;
    for (int i = 0; i < N; ++i)
        state_.P[i] = std::clamp(state_.P[i] + Pprime[i], options_.minPressure, pMax);
    // Hard-pin the last cell to the prescribed outlet pressure. The
    // correction above was built assuming a zero correction at a "ghost"
    // node just beyond the outlet face (a standard fixed-pressure outlet
    // treatment); explicitly resetting the last cell here keeps the
    // pressure level from drifting over many steps and is a common,
    // pragmatic simplification for this class of segregated solver (the
    // resulting small one-cell inconsistency is absorbed the same way the
    // paper's own method tolerates the <1% mass error of Eq. 24).
    state_.P[N - 1] = bc_.outletPressure;

    for (int f = 1; f < N; ++f) {
        const int cL = f - 1, cR = f;
        const double rho1f = 0.5 * (rho1At(cL) + rho1At(cR));
        const double rho2f = 0.5 * (rho2At(cL) + rho2At(cR));
        const double dP = (Pprime[cR] - Pprime[cL]) / state_.centerDistance(cL, cR);
        const double du1 = -dt / rho1f * dP;
        const double du2 = -dt / rho2f * dP;
        state_.u1[f] += du1;
        state_.u2[f] += du2;
        state_.ub[f] += du1; // bubbles grouped with layer 1
        state_.ud[f] += du2; // droplets grouped with layer 2
    }
    {
        const int c = N - 1;
        const double rho1f = rho1At(c);
        const double rho2f = rho2At(c);
        const double dP = (0.0 - Pprime[c]) / state_.cellWidth(c); // same ghost convention as above
        const double du1 = -dt / rho1f * dP;
        const double du2 = -dt / rho2f * dP;
        state_.u1[N] += du1;
        state_.u2[N] += du2;
        state_.ub[N] += du1;
        state_.ud[N] += du2;
    }

    for (int i = 0; i < N; ++i) rhoGasCell_[i] = fluid_.rhoGas(state_.P[i]);
}

void FourFieldSolver::updateContinuity(double dt) {
    const int N = state_.N;
    const double rhoL = fluid_.rhoLiquid;

    // Interior faces get a MUSCL/TVD-limited reconstruction (see
    // limitedFaceValue() above) when options_.advectionLimiter != None;
    // the inlet Dirichlet value stands in as the far-upwind ghost value at
    // the first interior face (f==1, flow forward), extending the
    // high-resolution stencil all the way to the first cell. The very last
    // interior face under backflow (f==N-1, vel<0) has no far-upwind
    // neighbour within the domain and falls back to first order there.
    auto upwind = [&](const std::vector<double>& cellVals, int f, double vel, double inletVal) -> double {
        if (f == 0) return inletVal;
        if (f == N) return cellVals[N - 1];
        const int cL = f - 1, cR = f;
        if (vel >= 0.0) {
            const bool haveFar = true; // inlet value substitutes for cL-1 when cL==0
            const double farVal = (cL - 1 >= 0) ? cellVals[cL - 1] : inletVal;
            return limitedFaceValue(cellVals[cL], cellVals[cR], farVal, haveFar, options_.advectionLimiter);
        } else {
            const bool haveFar = (cR + 1 <= N - 1);
            const double farVal = haveFar ? cellVals[cR + 1] : 0.0;
            return limitedFaceValue(cellVals[cR], cellVals[cL], farVal, haveFar, options_.advectionLimiter);
        }
    };

    std::vector<double> edNew(N), ebNew(N), eLNew(N), eGNew(N);

    // Surface-tension-motivated regularization (SolverOptions::
    // enableSurfaceTension): a biharmonic ("hyperdiffusion") term added to
    // the liquid-holdup continuity equation, -nu4*d^4(eL)/dz^4, with a
    // local coefficient nu4 = hyperdiffusionCoefficient*|ul|*dz^3. Built
    // from the START-OF-STEP eL field (captured once, before this loop
    // mutates anything), consistent with this scheme's other lagged
    // coefficients. See fourthDerivativeAt()'s comment and
    // SolverOptions::enableSurfaceTension for the full derivation/history.
    std::vector<double> eLField;
    if (options_.enableSurfaceTension) {
        eLField.resize(N);
        for (int i = 0; i < N; ++i) eLField[i] = state_.eL(i);
    }

    for (int i = 0; i < N; ++i) {
        const int fL = i, fR = i + 1;
        const double dz = state_.cellWidth(i); // this cell's own width: a proper FV divergence

        const double edL = upwind(state_.ed, fL, state_.ud[fL], inletEd());
        const double edR = upwind(state_.ed, fR, state_.ud[fR], inletEd());
        const double fluxEdR = edR * state_.ud[fR];
        const double fluxEdL = edL * state_.ud[fL];
        edNew[i] = state_.ed[i] + dt * (-(fluxEdR - fluxEdL) / dz + (Ue_[i] - Ud_[i]) / rhoL);

        const double ebL = upwind(state_.eb, fL, state_.ub[fL], inletEb());
        const double ebR = upwind(state_.eb, fR, state_.ub[fR], inletEb());
        const double fluxEbR = ebR * state_.ub[fR];
        const double fluxEbL = ebL * state_.ub[fL];
        ebNew[i] = state_.eb[i] + dt * (-(fluxEbR - fluxEbL) / dz + (phiE_[i] - phiDe_[i]) / std::max(rhoGasCell_[i], tiny));

        const double elL = upwind(state_.el, fL, state_.ul[fL], inletEl());
        const double elR = upwind(state_.el, fR, state_.ul[fR], inletEl());
        const double udAtL = upwind(state_.ed, fL, state_.ud[fL], inletEd()); // reuse ed upwind for the ed*ud term
        const double udAtR = upwind(state_.ed, fR, state_.ud[fR], inletEd());
        const double fluxLR = elR * state_.ul[fR] + udAtR * state_.ud[fR];
        const double fluxLL = elL * state_.ul[fL] + udAtL * state_.ud[fL];
        eLNew[i] = state_.eL(i) + dt * (-(fluxLR - fluxLL) / dz);

        const double egL = upwind(state_.eg, fL, state_.ug[fL], inletEg());
        const double egR = upwind(state_.eg, fR, state_.ug[fR], inletEg());
        const double ubAtL = upwind(state_.eb, fL, state_.ub[fL], inletEb());
        const double ubAtR = upwind(state_.eb, fR, state_.ub[fR], inletEb());
        const double fluxGR = egR * state_.ug[fR] + ubAtR * state_.ub[fR];
        const double fluxGL = egL * state_.ug[fL] + ubAtL * state_.ub[fL];
        eGNew[i] = state_.eG(i) + dt * (-(fluxGR - fluxGL) / dz);

        if (options_.enableSurfaceTension) {
            const double ulLocal = 0.5 * (state_.ul[fL] + state_.ul[fR]);
            const double nu4 = options_.hyperdiffusionCoefficient * std::fabs(ulLocal) * dz * dz * dz;
            const double biharmonic = -nu4 * fourthDerivativeAt(state_, eLField, i);
            // Equal and opposite on eG: this is a redistribution of the
            // interface position, not a mass source, and keeps
            // eLNew+eGNew invariant before the renormalisation below even
            // has to act (the renormalisation would otherwise absorb any
            // imbalance asymmetrically between the two fields).
            eLNew[i] += dt * biharmonic;
            eGNew[i] -= dt * biharmonic;
        }
    }

    for (int i = 0; i < N; ++i) {
        state_.ed[i] = std::max(edNew[i], 0.0);
        state_.eb[i] = std::max(ebNew[i], 0.0);
        const double eL = std::clamp(eLNew[i], small_e, 1.0 - small_e);
        const double eG = std::clamp(eGNew[i], small_e, 1.0 - small_e);
        // Renormalise eL, eG to sum to 1 exactly (kinematic constraint); the
        // pre-renormalisation drift is the model's own mass-conservation
        // error and can be inspected via boundaryMassFluxes() /
        // totalLiquidMass() etc. before this step, per Eq. (24).
        const double sum = eL + eG;
        state_.el[i] = std::clamp(eL / sum, state_.ed[i], 1.0);
        state_.eg[i] = std::clamp(eG / sum, state_.eb[i], 1.0);
    }
    clampVolumeFractions();
}

void FourFieldSolver::clampVolumeFractions() {
    for (int i = 0; i < state_.N; ++i) {
        double eL = state_.el[i] + state_.ed[i];
        double eG = state_.eg[i] + state_.eb[i];
        // Eq. (1): eL + eG = 1 identically; re-derive el, eg (Eq. 21).
        eL = std::clamp(eL, small_e, 1.0 - small_e);
        eG = 1.0 - eL;
        state_.ed[i] = std::clamp(state_.ed[i], 0.0, eL);
        state_.el[i] = eL - state_.ed[i];
        state_.eb[i] = std::clamp(state_.eb[i], 0.0, eG);
        state_.eg[i] = eG - state_.eb[i];
    }
}

void FourFieldSolver::applyInletBoundary() {
    const double ul0 = bc_.inletSuperficialLiquid / std::max(inletEL(), small_e);
    const double ug0 = bc_.inletSuperficialGas / std::max(inletEG(), small_e);
    state_.u1[0] = ul0;
    state_.u2[0] = ug0;
    state_.ud[0] = ul0;
    state_.ub[0] = ug0;
    state_.ul[0] = ul0;
    state_.ug[0] = ug0;
}

void FourFieldSolver::applyOutletBoundary() {
    const int N = state_.N;
    state_.u1[N] = state_.u1[N - 1];
    state_.u2[N] = state_.u2[N - 1];
    state_.ud[N] = state_.ud[N - 1];
    state_.ub[N] = state_.ub[N - 1];
    state_.ul[N] = state_.ul[N - 1];
    state_.ug[N] = state_.ug[N - 1];
    // Outlet pressure itself is enforced in solvePressureCorrection().
}

namespace {
    // Lightweight opt-in diagnostic: set the MFS_DEBUG environment variable
    // to trace which field and which pipeline stage first produces a
    // NaN/Inf, e.g. when experimenting with custom closures. Zero overhead
    // when unset (checked once, cached).
    bool hasBad(const std::vector<double>& v, int& idx) {
        for (std::size_t i = 0; i < v.size(); ++i)
            if (std::isnan(v[i]) || std::isinf(v[i])) { idx = static_cast<int>(i); return true; }
        idx = -1;
        return false;
    }
    void debugCheck(const char* stage, const FlowState& s) {
        static const bool on = std::getenv("MFS_DEBUG") != nullptr;
        static long stepCounter = 0;
        if (!on) return;
        if (std::string(stage) == "layer momentum") ++stepCounter;
        int idx;
        const std::pair<const char*, const std::vector<double>*> fields[] = {
            {"el", &s.el}, {"eb", &s.eb}, {"P", &s.P},
            {"u1", &s.u1}, {"u2", &s.u2}, {"ud", &s.ud}, {"ub", &s.ub},
        };
        for (const auto& [name, vec] : fields) {
            if (hasBad(*vec, idx)) {
                std::fprintf(stderr, "[MFS_DEBUG] step=%ld stage='%s': NaN/Inf in %s at index %d\n",
                             stepCounter, stage, name, idx);
            }
        }
    }
}

double FourFieldSolver::step(double dt) {
    computeGeometry();
    computeClosures();

    updateLayerMomentum(dt);
    debugCheck("layer momentum", state_);
    backSubstitutePhaseVelocities();
    updateDispersedMomentum(dt);
    debugCheck("dispersed momentum", state_);
    backSubstitutePhaseVelocities();

    solvePressureCorrection(dt);
    debugCheck("pressure correction", state_);
    backSubstitutePhaseVelocities();

    updateContinuity(dt);
    debugCheck("continuity", state_);

    time_ += dt;
    ++stepCount_;
    if (options_.amr.enabled && stepCount_ % std::max(1, options_.amr.adaptEveryNSteps) == 0) {
        adaptMesh();
    }
    if (options_.movingMesh.enabled &&
        stepCount_ % std::max(1, options_.movingMesh.relocateEveryNSteps) == 0) {
        relocateMesh();
    }
    return dt;
}

double FourFieldSolver::totalLiquidMass() const {
    double m = 0.0;
    const double A = geometry_.area();
    for (int i = 0; i < state_.N; ++i) m += state_.eL(i) * fluid_.rhoLiquid * A * state_.cellWidth(i);
    return m;
}
double FourFieldSolver::totalGasMass() const {
    double m = 0.0;
    const double A = geometry_.area();
    for (int i = 0; i < state_.N; ++i) m += state_.eG(i) * rhoGasCell_[i] * A * state_.cellWidth(i);
    return m;
}
double FourFieldSolver::totalDropletMass() const {
    double m = 0.0;
    const double A = geometry_.area();
    for (int i = 0; i < state_.N; ++i) m += state_.ed[i] * fluid_.rhoLiquid * A * state_.cellWidth(i);
    return m;
}
double FourFieldSolver::totalBubbleMass() const {
    double m = 0.0;
    const double A = geometry_.area();
    for (int i = 0; i < state_.N; ++i) m += state_.eb[i] * rhoGasCell_[i] * A * state_.cellWidth(i);
    return m;
}

FourFieldSolver::MassFluxes FourFieldSolver::boundaryMassFluxes() const {
    const double A = geometry_.area();
    const int N = state_.N;
    MassFluxes mf{};
    mf.liquidIn = (inletEl() * state_.ul[0] + inletEd() * state_.ud[0]) * fluid_.rhoLiquid * A;
    mf.liquidOut = (state_.el[N - 1] * state_.ul[N] + state_.ed[N - 1] * state_.ud[N]) * fluid_.rhoLiquid * A;
    mf.gasIn = (inletEg() * state_.ug[0] + inletEb() * state_.ub[0]) * rhoGasCell_[0] * A;
    mf.gasOut = (state_.eg[N - 1] * state_.ug[N] + state_.eb[N - 1] * state_.ub[N]) * rhoGasCell_[N - 1] * A;
    return mf;
}

std::vector<double> FourFieldSolver::computeIndicator() const {
    const int N = state_.N;
    const double D = geometry_.diameter();
    const double rhoL = fluid_.rhoLiquid;
    std::vector<double> F(N, 0.0);
    for (int i = 0; i < N; ++i) {
        const auto g = geometry_.fromAreaFraction(state_.e1(i));
        const double rhoG = fluid_.rhoGas(state_.P[i]);
        const double ulC = 0.5 * (state_.ul[i] + state_.ul[i + 1]);
        const double ugC = 0.5 * (state_.ug[i] + state_.ug[i + 1]);

        InterfacialFrictionInputs ifi{};
        ifi.densityGas = rhoG;
        ifi.densityLiquid = rhoL;
        ifi.velocityGas = ugC;
        ifi.velocityLiquid = ulC;
        ifi.hydraulicDiameter1 = g.D1;
        ifi.hydraulicDiameter2 = g.D2;
        ifi.liquidHeight = g.h1;
        ifi.pipeDiameter = D;
        ifi.area1 = g.A1;
        ifi.area2 = g.A2;
        ifi.areaFraction2 = g.A2 / geometry_.area();
        ifi.dA1dh1 = g.Si;
        ifi.inclination = state_.theta[i];
        F[i] = kelvinHelmholtzParameterF(ifi);
    }
    return F;
}

std::vector<double> FourFieldSolver::refinementIndicatorProfile() const { return computeIndicator(); }

bool FourFieldSolver::adaptMesh() {
    const int N = state_.N;
    if (N < 2) return false;

    const std::vector<double> F = computeIndicator();

    const double minWidth = initialDz_ * options_.amr.minCellWidthFraction;
    const double maxWidth = initialDz_ * options_.amr.maxCellWidthFraction;
    const int maxCells = std::max(initialN_, static_cast<int>(initialN_ * options_.amr.maxCellCountFactor));

    std::vector<char> refine(N, 0), coarsen(N, 0);
    int cellBudget = maxCells - N; // extra cells still allowed; each split adds exactly one
    for (int i = 0; i < N; ++i) {
        const double w = state_.cellWidth(i);
        if (F[i] > options_.amr.refineThreshold && w > 2.0 * minWidth && cellBudget > 0) {
            refine[i] = 1;
            --cellBudget;
        } else if (F[i] < options_.amr.coarsenThreshold && w < maxWidth) {
            coarsen[i] = 1;
        }
    }

    // Rebuild the mesh in a single left-to-right pass. No hanging nodes: a
    // "group" is either one refined cell (split into two half-width
    // cells), a pair of adjacent coarsen-flagged cells (merged into one,
    // width-weighted average of their state -- exactly mass-conservative
    // for the volume fractions since mass = e*rho*A*width sums correctly),
    // or a single unchanged cell. Splitting duplicates the parent's
    // cell-centred state into both children, which is likewise exactly
    // conservative (two half-width cells at the parent's own value
    // integrate to the same total as the parent). Face velocities at a new
    // interior split face are linearly interpolated between the cell's own
    // two original bounding faces.
    std::vector<double> nFaceZ, nEl, nEd, nEg, nEb, nP, nTheta;
    std::vector<double> nU1, nU2, nUd, nUb, nUl, nUg;
    nFaceZ.reserve(N + 2);
    nEl.reserve(N); nEd.reserve(N); nEg.reserve(N); nEb.reserve(N); nP.reserve(N); nTheta.reserve(N);
    nU1.reserve(N + 2); nU2.reserve(N + 2); nUd.reserve(N + 2); nUb.reserve(N + 2);
    nUl.reserve(N + 2); nUg.reserve(N + 2);

    auto emitFace = [&](int origFaceIdx) {
        nFaceZ.push_back(state_.faceZ[origFaceIdx]);
        nU1.push_back(state_.u1[origFaceIdx]);
        nU2.push_back(state_.u2[origFaceIdx]);
        nUd.push_back(state_.ud[origFaceIdx]);
        nUb.push_back(state_.ub[origFaceIdx]);
        nUl.push_back(state_.ul[origFaceIdx]);
        nUg.push_back(state_.ug[origFaceIdx]);
    };
    auto emitCellCopy = [&](int i) {
        nEl.push_back(state_.el[i]); nEd.push_back(state_.ed[i]);
        nEg.push_back(state_.eg[i]); nEb.push_back(state_.eb[i]);
        nP.push_back(state_.P[i]); nTheta.push_back(state_.theta[i]);
    };

    emitFace(0);
    bool changed = false;
    int i = 0;
    while (i < N) {
        if (coarsen[i] && i + 1 < N && coarsen[i + 1]) {
            const double w0 = state_.cellWidth(i), w1 = state_.cellWidth(i + 1);
            const double wSum = w0 + w1;
            auto blend = [&](const std::vector<double>& v) { return (v[i] * w0 + v[i + 1] * w1) / wSum; };
            nEl.push_back(blend(state_.el));
            nEd.push_back(blend(state_.ed));
            nEg.push_back(blend(state_.eg));
            nEb.push_back(blend(state_.eb));
            nP.push_back(blend(state_.P));
            nTheta.push_back(blend(state_.theta));
            emitFace(i + 2); // drop the shared interior face i+1
            changed = true;
            i += 2;
        } else if (refine[i]) {
            emitCellCopy(i);
            emitCellCopy(i);
            const double zMid = 0.5 * (state_.faceZ[i] + state_.faceZ[i + 1]);
            nFaceZ.push_back(zMid);
            nU1.push_back(0.5 * (state_.u1[i] + state_.u1[i + 1]));
            nU2.push_back(0.5 * (state_.u2[i] + state_.u2[i + 1]));
            nUd.push_back(0.5 * (state_.ud[i] + state_.ud[i + 1]));
            nUb.push_back(0.5 * (state_.ub[i] + state_.ub[i + 1]));
            nUl.push_back(0.5 * (state_.ul[i] + state_.ul[i + 1]));
            nUg.push_back(0.5 * (state_.ug[i] + state_.ug[i + 1]));
            emitFace(i + 1);
            changed = true;
            i += 1;
        } else {
            emitCellCopy(i);
            emitFace(i + 1);
            i += 1;
        }
    }

    if (!changed) return false;

    const int newN = static_cast<int>(nEl.size());
    state_.N = newN;
    state_.faceZ = std::move(nFaceZ);
    state_.el = std::move(nEl);
    state_.ed = std::move(nEd);
    state_.eg = std::move(nEg);
    state_.eb = std::move(nEb);
    state_.P = std::move(nP);
    state_.theta = std::move(nTheta);
    state_.u1 = std::move(nU1);
    state_.u2 = std::move(nU2);
    state_.ud = std::move(nUd);
    state_.ub = std::move(nUb);
    state_.ul = std::move(nUl);
    state_.ug = std::move(nUg);

    // Cache arrays are fully recomputed from state_ every step (see
    // computeGeometry()/computeClosures()); only their SIZE matters here.
    geom_.resize(newN);
    tauW1_.assign(newN + 1, 0.0);
    tauW2_.assign(newN + 1, 0.0);
    tauI_.assign(newN + 1, 0.0);
    Ue_.assign(newN, 0.0);
    Ud_.assign(newN, 0.0);
    phiE_.assign(newN, 0.0);
    phiDe_.assign(newN, 0.0);
    rhoGasCell_.assign(newN, fluid_.rhoGas(bc_.outletPressure));

    return true;
}

namespace {
    // Exact conservative interpolation of a piecewise-constant
    // cell-centred field from one 1D partition onto another: for each new
    // cell, integrate the old field over its overlap with each old cell it
    // spans, divided by the new cell's width. Standard two-pointer sweep
    // over two sorted face arrays sharing the same endpoints.
    std::vector<double> conservativeRemapCellField(const std::vector<double>& oldVals,
                                                     const std::vector<double>& oldFaceZ,
                                                     const std::vector<double>& newFaceZ) {
        const int nOld = static_cast<int>(oldVals.size());
        const int nNew = static_cast<int>(newFaceZ.size()) - 1;
        std::vector<double> out(nNew, 0.0);
        int oi = 0;
        for (int ni = 0; ni < nNew; ++ni) {
            const double lo = newFaceZ[ni], hi = newFaceZ[ni + 1];
            double acc = 0.0;
            while (oi < nOld && oldFaceZ[oi] < hi) {
                const double oLo = std::max(lo, oldFaceZ[oi]);
                const double oHi = std::min(hi, oldFaceZ[oi + 1]);
                if (oHi > oLo) acc += oldVals[oi] * (oHi - oLo);
                if (oldFaceZ[oi + 1] > hi) break; // this old cell still overlaps the NEXT new cell too
                ++oi;
            }
            out[ni] = acc / std::max(hi - lo, tiny);
        }
        return out;
    }

    // Linear interpolation of a face-centred (point-sampled) field from
    // the old face positions onto the new ones. Not "conservative" in the
    // integral sense above -- there is no natural integral of a point
    // sample -- linear interpolation is the standard, appropriate
    // treatment for this kind of data (matching how new interior faces
    // were handled when splitting a cell in adaptMesh()).
    std::vector<double> linearRemapFaceField(const std::vector<double>& oldVals,
                                              const std::vector<double>& oldFaceZ,
                                              const std::vector<double>& newFaceZ) {
        const int nOldFaces = static_cast<int>(oldFaceZ.size());
        std::vector<double> out(newFaceZ.size());
        int k = 0;
        for (std::size_t j = 0; j < newFaceZ.size(); ++j) {
            const double z = newFaceZ[j];
            while (k < nOldFaces - 2 && oldFaceZ[k + 1] < z) ++k;
            const double z0 = oldFaceZ[k], z1 = oldFaceZ[k + 1];
            const double t = (z1 > z0) ? (z - z0) / (z1 - z0) : 0.0;
            out[j] = oldVals[k] + std::clamp(t, 0.0, 1.0) * (oldVals[k + 1] - oldVals[k]);
        }
        return out;
    }
}

std::vector<double> FourFieldSolver::computeMonitorFunction() const {
    const int N = state_.N;
    std::vector<double> M(N, 1.0);
    if (N < 2) return M;

    const std::vector<double> F = computeIndicator();
    const double L = state_.L;
    const auto& opt = options_.movingMesh;
    // Reused purely as a fixed reference scale (Andreussi & Persen's own
    // F0), independent of whether h-refinement AMR is itself enabled.
    const double F0 = options_.amr.refineThreshold;

    for (int i = 0; i < N; ++i) {
        const int iL = std::max(i - 1, 0);
        const int iR = std::min(i + 1, N - 1);
        const double dz = state_.centerDistance(iL, iR);
        const double gradEl = (dz > 0.0) ? std::fabs(state_.eL(iR) - state_.eL(iL)) / dz : 0.0;
        const double fTerm = std::min(F[i] / F0, opt.khIndicatorCap);
        M[i] = 1.0 + opt.holdupGradientWeight * gradEl * L + opt.khIndicatorWeight * fTerm;
    }

    for (int pass = 0; pass < opt.monitorSmoothingPasses; ++pass) {
        std::vector<double> Msmooth(N);
        for (int i = 0; i < N; ++i) {
            const int iL = std::max(i - 1, 0);
            const int iR = std::min(i + 1, N - 1);
            Msmooth[i] = 0.25 * M[iL] + 0.5 * M[i] + 0.25 * M[iR];
        }
        M = std::move(Msmooth);
    }

    for (double& v : M) v = std::min(v, opt.monitorCap);
    return M;
}

bool FourFieldSolver::relocateMesh() {
    const int N = state_.N;
    if (N < 3) return false;

    const std::vector<double> M = computeMonitorFunction();

    // Cumulative integral of M (piecewise-constant per cell) over the OLD
    // mesh, i.e. Theta at each OLD face.
    std::vector<double> theta(N + 1, 0.0);
    for (int i = 0; i < N; ++i) theta[i + 1] = theta[i] + M[i] * state_.cellWidth(i);
    const double thetaTotal = theta[N];
    if (thetaTotal <= 0.0) return false;

    // Equidistribution (de Boor, 1974): the target new face j is where the
    // cumulative monitor integral reaches (j/N) of its total. Theta is
    // piecewise LINEAR in z (since M is piecewise constant), so inverting
    // it within the bracketing old cell is a direct linear solve.
    std::vector<double> targetZ(N + 1);
    targetZ[0] = state_.faceZ[0];
    targetZ[N] = state_.faceZ[N];
    int oldCell = 0;
    for (int j = 1; j < N; ++j) {
        const double target = (static_cast<double>(j) / N) * thetaTotal;
        while (oldCell < N - 1 && theta[oldCell + 1] < target) ++oldCell;
        const double slope = M[oldCell];
        const double frac = (slope > 0.0) ? (target - theta[oldCell]) / (slope * state_.cellWidth(oldCell)) : 0.0;
        targetZ[j] = state_.faceZ[oldCell] + std::clamp(frac, 0.0, 1.0) * state_.cellWidth(oldCell);
    }

    // Under-relax toward the target (a full jump to the equidistributed
    // mesh every relocation would let the grid itself oscillate violently
    // near a sharp, still-forming front), then repair monotonicity: a
    // heavily relaxed point could in principle end up on the wrong side of
    // a neighbour if consecutive targets differ a lot, so nudge onto a
    // strictly increasing sequence with a small minimum gap afterward.
    std::vector<double> newZ(N + 1);
    newZ[0] = state_.faceZ[0];
    newZ[N] = state_.faceZ[N];
    const double relax = std::clamp(options_.movingMesh.relaxation, 0.0, 1.0);
    for (int j = 1; j < N; ++j) newZ[j] = state_.faceZ[j] + relax * (targetZ[j] - state_.faceZ[j]);

    const double minGap = 1.0e-6 * state_.L;
    for (int j = 1; j < N; ++j) newZ[j] = std::max(newZ[j], newZ[j - 1] + minGap);
    for (int j = N - 1; j >= 1; --j) newZ[j] = std::min(newZ[j], newZ[j + 1] - minGap);
    for (int j = 1; j <= N; ++j)
        if (newZ[j] <= newZ[j - 1]) return false; // degenerate; leave the mesh untouched this step

    const std::vector<double> oldFaceZ = state_.faceZ;
    state_.el = conservativeRemapCellField(state_.el, oldFaceZ, newZ);
    state_.ed = conservativeRemapCellField(state_.ed, oldFaceZ, newZ);
    state_.eg = conservativeRemapCellField(state_.eg, oldFaceZ, newZ);
    state_.eb = conservativeRemapCellField(state_.eb, oldFaceZ, newZ);
    state_.P = conservativeRemapCellField(state_.P, oldFaceZ, newZ);
    state_.theta = conservativeRemapCellField(state_.theta, oldFaceZ, newZ);

    state_.u1 = linearRemapFaceField(state_.u1, oldFaceZ, newZ);
    state_.u2 = linearRemapFaceField(state_.u2, oldFaceZ, newZ);
    state_.ud = linearRemapFaceField(state_.ud, oldFaceZ, newZ);
    state_.ub = linearRemapFaceField(state_.ub, oldFaceZ, newZ);
    state_.ul = linearRemapFaceField(state_.ul, oldFaceZ, newZ);
    state_.ug = linearRemapFaceField(state_.ug, oldFaceZ, newZ);

    state_.faceZ = newZ;
    // N is unchanged by relocation, so the cache array sizes need no
    // adjustment; their contents are recomputed fresh every step regardless.
    return true;
}

} // namespace mfs
