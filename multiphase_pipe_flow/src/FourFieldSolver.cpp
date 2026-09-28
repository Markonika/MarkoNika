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
    state_.resize(nCells);
    state_.L = length;
    state_.dz = length / nCells;
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
}

void FourFieldSolver::computeGeometry() {
    for (int i = 0; i < state_.N; ++i) {
        geom_[i] = geometry_.fromAreaFraction(state_.e1(i));
        rhoGasCell_[i] = fluid_.rhoGas(state_.P[i]);
    }
}

double FourFieldSolver::stableTimeStep() const {
    double umax = tiny;
    for (double v : state_.u1) umax = std::max(umax, std::fabs(v));
    for (double v : state_.u2) umax = std::max(umax, std::fabs(v));
    for (double v : state_.ud) umax = std::max(umax, std::fabs(v));
    for (double v : state_.ub) umax = std::max(umax, std::fabs(v));
    double dt = options_.courantTarget * state_.dz / umax;
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
        const double fi = interfacialFrictionFactor(ifi, options_.interfacialCorrelation);

        tauW1_[f] = wallShearStress(flw, rhoL, ulO);
        tauW2_[f] = wallShearStress(fgw, rhoGf, ugO);
        tauI_[f] = interfacialShearStress(fi, rhoGf, ugO, ulO);
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

        DropletDepositionInputs ddi{D, rhoL, state_.ed[i], state_.eg[i], 0.1};
        Ud_[i] = dropletDepositionRate(ddi);
    }
}

namespace {
    inline double rho1Of(const FlowState& s, const std::vector<double>& rhoGasCell, double rhoL, int i) {
        const double e1 = s.e1(i);
        return (s.el[i] * rhoL + s.eb[i] * rhoGasCell[i]) / std::max(e1, small_e);
    }
    inline double rho2Of(const FlowState& s, const std::vector<double>& rhoGasCell, double rhoL, int i) {
        const double e2 = s.e2(i);
        return (s.eg[i] * rhoGasCell[i] + s.ed[i] * rhoL) / std::max(e2, small_e);
    }
}

void FourFieldSolver::updateLayerMomentum(double dt) {
    const double A = geometry_.area();
    const double rhoL = fluid_.rhoLiquid;
    const double dz = state_.dz;
    const int N = state_.N;

    std::vector<double> u1New = state_.u1;
    std::vector<double> u2New = state_.u2;

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

        const double dPdz = (state_.P[cR] - state_.P[cL]) / dz;
        const double thetaF = 0.5 * (state_.theta[cL] + state_.theta[cR]);
        const double dh1dz = (geom_[cR].h1 - geom_[cL].h1) / dz;
        const double Swp1f = 0.5 * (geom_[cL].Swp1 + geom_[cR].Swp1);
        const double Swp2f = 0.5 * (geom_[cL].Swp2 + geom_[cR].Swp2);
        const double Sif = 0.5 * (geom_[cL].Si + geom_[cR].Si);

        const double UeF = 0.5 * (Ue_[cL] + Ue_[cR]);
        const double UdF = 0.5 * (Ud_[cL] + Ud_[cR]);
        const double phiEF = 0.5 * (phiE_[cL] + phiE_[cR]);
        const double phiDeF = 0.5 * (phiDe_[cL] + phiDe_[cR]);

        const double ulF = state_.ul[f], ugF = state_.ug[f], udF = state_.ud[f], ubF = state_.ub[f];

        // --- Layer 1 momentum (Eq. 6, slip flux term dropped per Appendix A) ---
        const double u1c = state_.u1[f];
        const double du1dzAdv = (u1c >= 0.0) ? (state_.u1[f] - state_.u1[f - 1]) / dz
                                              : (state_.u1[f + 1] - state_.u1[f]) / dz;
        const double massSrc1 = (-UeF * ulF + UdF * udF + phiEF * ugF - phiDeF * ubF) / (e1fSafe * rho1f);
        const double du1dt = -u1c * du1dzAdv
                              - dPdz / rho1f
                              - gravity * std::cos(thetaF) * dh1dz
                              - gravity * std::sin(thetaF)
                              - tauW1_[f] * Swp1f / (A * e1fSafe * rho1f)
                              + tauI_[f] * Sif / (A * e1fSafe * rho1f)
                              + massSrc1;
        u1New[f] = state_.u1[f] + dt * du1dt;

        // --- Layer 2 momentum (Eq. 7, slip flux term dropped) ---
        const double u2c = state_.u2[f];
        const double du2dzAdv = (u2c >= 0.0) ? (state_.u2[f] - state_.u2[f - 1]) / dz
                                              : (state_.u2[f + 1] - state_.u2[f]) / dz;
        const double massSrc2 = (UeF * ulF - UdF * udF - phiEF * ugF + phiDeF * ubF) / (e2fSafe * rho2f);
        const double du2dt = -u2c * du2dzAdv
                              - dPdz / rho2f
                              - gravity * std::cos(thetaF) * dh1dz
                              - gravity * std::sin(thetaF)
                              - tauW2_[f] * Swp2f / (A * e2fSafe * rho2f)
                              - tauI_[f] * Sif / (A * e2fSafe * rho2f)
                              + massSrc2;
        u2New[f] = state_.u2[f] + dt * du2dt;
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
    const double dz = state_.dz;
    const int N = state_.N;

    std::vector<double> udNew = state_.ud;
    std::vector<double> ubNew = state_.ub;

    for (int f = 1; f < N; ++f) {
        const int cL = f - 1, cR = f;
        const double edF = 0.5 * (state_.ed[cL] + state_.ed[cR]);
        const double elF = 0.5 * (state_.el[cL] + state_.el[cR]);
        const double ebF = 0.5 * (state_.eb[cL] + state_.eb[cR]);
        const double rhoGf = 0.5 * (rhoGasCell_[cL] + rhoGasCell_[cR]);
        const double dPdz = (state_.P[cR] - state_.P[cL]) / dz;
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
    const double dz = state_.dz;
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

        R[cL] += qTotal / dz;
        R[cR] -= qTotal / dz;

        const double Gf = (dt * (elF + ebF) / rho1f + dt * (edF + egF) / rho2f) / dz;
        b[cL] += Gf / dz;
        c[cL] += -Gf / dz;
        a[cR] += -Gf / dz;
        b[cR] += Gf / dz;
    }

    // Inlet face (fixed flow, no pressure sensitivity): contributes to R[0] only.
    {
        const double qLiquidIn = inletEl() * state_.ul[0] + inletEd() * state_.ud[0];
        const double qGasIn = inletEg() * state_.ug[0] + inletEb() * state_.ub[0];
        R[0] -= (qLiquidIn + qGasIn) / dz;
    }

    // Outlet face (fixed pressure, P'=0 there): contributes to R[N-1] and its own diagonal.
    {
        const int c = N - 1;
        const double rho1f = rho1At(c);
        const double rho2f = rho2At(c);
        const double elF = state_.el[c], edF = state_.ed[c], egF = state_.eg[c], ebF = state_.eb[c];
        const double qLiquidOut = elF * state_.ul[N] + edF * state_.ud[N];
        const double qGasOut = egF * state_.ug[N] + ebF * state_.ub[N];
        R[c] += (qLiquidOut + qGasOut) / dz;

        const double Gf = (dt * (elF + ebF) / rho1f + dt * (edF + egF) / rho2f) / dz;
        b[c] += Gf / dz;
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
        const double dP = (Pprime[cR] - Pprime[cL]) / dz;
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
        const double dP = (0.0 - Pprime[c]) / dz;
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
    const double dz = state_.dz;
    const double rhoL = fluid_.rhoLiquid;

    auto upwind = [&](const std::vector<double>& cellVals, int f, double vel, double inletVal) -> double {
        if (f == 0) return inletVal;
        if (f == N) return cellVals[N - 1];
        return (vel >= 0.0) ? cellVals[f - 1] : cellVals[f];
    };

    std::vector<double> edNew(N), ebNew(N), eLNew(N), eGNew(N);

    for (int i = 0; i < N; ++i) {
        const int fL = i, fR = i + 1;

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
    return dt;
}

double FourFieldSolver::totalLiquidMass() const {
    double m = 0.0;
    const double A = geometry_.area();
    for (int i = 0; i < state_.N; ++i) m += state_.eL(i) * fluid_.rhoLiquid * A * state_.dz;
    return m;
}
double FourFieldSolver::totalGasMass() const {
    double m = 0.0;
    const double A = geometry_.area();
    for (int i = 0; i < state_.N; ++i) m += state_.eG(i) * rhoGasCell_[i] * A * state_.dz;
    return m;
}
double FourFieldSolver::totalDropletMass() const {
    double m = 0.0;
    const double A = geometry_.area();
    for (int i = 0; i < state_.N; ++i) m += state_.ed[i] * fluid_.rhoLiquid * A * state_.dz;
    return m;
}
double FourFieldSolver::totalBubbleMass() const {
    double m = 0.0;
    const double A = geometry_.area();
    for (int i = 0; i < state_.N; ++i) m += state_.eb[i] * rhoGasCell_[i] * A * state_.dz;
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

} // namespace mfs
