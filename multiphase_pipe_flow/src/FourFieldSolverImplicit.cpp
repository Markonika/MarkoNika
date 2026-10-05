// Jacobian-free Newton-Krylov (JFNK) fully implicit pressure-velocity
// solve. See FourFieldSolver.hpp, stepImplicitPressureVelocity() for the
// public-facing documentation and scoping summary. Kept as a separate
// translation unit from FourFieldSolver.cpp so none of the existing,
// heavily regression-tested explicit/semi-implicit pipeline is touched by
// this addition.
//
// SCOPE. Only (u1, u2, P) are Newton unknowns. Volume fractions
// (el, ed, eg, eb), dispersed-phase velocities (ud, ub), and the
// entrainment/deposition/disengagement mass-transfer rates are held fixed
// at their start-of-step values throughout the Newton iteration -- exactly
// the same lagging step()'s own segregated scheme already uses for these
// quantities (updateDispersedMomentum and updateContinuity both run AFTER
// the momentum/pressure update, using its result). What genuinely changes
// relative to step() is that advection, the pressure gradient, and BOTH
// wall and interfacial friction are evaluated NONLINEARLY at the trial
// (u1, u2) -- not hand-linearised the way enableImplicitFriction's local
// 2x2 solve does it -- and the gas-compressibility term in the pressure
// equation uses the true nonlinear accumulation (rho_gas(P_new) -
// rho_gas(P_old))/dt rather than solvePressureCorrection's linearised
// coefficient. This is the real pressure-velocity-friction coupling this
// codebase's own prior investigation (VALIDATION.md item 4) identified as
// unaddressed.
//
// RESIDUAL. Unknown vector x, size 3N-2:
//   x[0 .. N-2]     = u1 at interior faces 1..N-1
//   x[N-1 .. 2N-3]  = u2 at interior faces 1..N-1
//   x[2N-2 .. 3N-3] = P  at cells 0..N-1
// Face 0 (inlet) and face N (outlet) velocities are not unknowns: the
// inlet is a fixed Dirichlet velocity (applyInletBoundary()'s own
// formula), the outlet a zero-gradient extrapolation of the last interior
// face (applyOutletBoundary()'s own convention) -- both applied directly
// inside the residual evaluator. Cell N-1's pressure equation is replaced
// by the Dirichlet outlet-pressure residual (P[N-1] - outletPressure),
// matching how solvePressureCorrection already hard-pins that cell.
//
// Momentum residual at face f: R_u[f] = (u_trial[f] - u_old[f]) -
// dt * du/dt(trial), backward Euler on exactly the same physics as
// updateLayerMomentum's explicit branch (same advection, pressure
// gradient, gravity-term treatment including the opt-in well-balanced
// discretization, and mass-transfer/turbulent-viscosity source terms,
// all lagged per the scoping above), except friction, which is
// RECOMPUTED from the closures at the trial velocity rather than read
// from the lagged tauW1_/tauW2_/tauI_ arrays.
//
// Pressure residual at cell i<N-1: the same net-volumetric-flux-divergence
// construction solvePressureCorrection() already uses (qLiquid, qGas,
// qTotal, inlet/outlet boundary treatment), evaluated at the trial
// velocities and fixed (lagged) upwind volume fractions, plus the true
// nonlinear gas-accumulation term in place of solvePressureCorrection's
// linearised coefficient.
#include "mfs/Closures.hpp"
#include "mfs/Constants.hpp"
#include "mfs/FourFieldSolver.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <vector>

namespace mfs {

using constants::gravity;
using constants::small_e;
using constants::tiny;

namespace {

// Local duplicates of FourFieldSolver.cpp's own anonymous-namespace
// rho1Of/rho2Of (small enough that duplicating is simpler and safer than
// sharing a header just for these two one-line functions across the
// translation-unit boundary this file deliberately keeps).
inline double rho1Of(const FlowState& s, const std::vector<double>& rhoGasCell, double rhoL, int i) {
    const double e1 = s.e1(i);
    return (s.el[i] * rhoL + s.eb[i] * rhoGasCell[i]) / std::max(e1, small_e);
}
inline double rho2Of(const FlowState& s, const std::vector<double>& rhoGasCell, double rhoL, int i) {
    const double e2 = s.e2(i);
    return (s.eg[i] * rhoGasCell[i] + s.ed[i] * rhoL) / std::max(e2, small_e);
}

// Restarted GMRES (Saad & Schultz 1986) with Givens-rotation least-squares,
// no preconditioner. matVec(v, out) must compute out = J*v. Solves
// J*x = b approximately, starting from (and overwriting) x. Returns the
// number of matrix-vector products consumed. This is the ONLY linear
// solve used by the Newton loop below; its correctness is cross-checked
// indirectly and robustly regardless of any subtle implementation error
// here, since the outer Newton loop independently re-evaluates the true
// (nonlinear) residual after every step and simply iterates (or reports
// non-convergence) rather than trusting the linear solve blindly.
template <typename MatVec>
int gmresSolve(MatVec matVec, const std::vector<double>& b, std::vector<double>& x,
               int restart, int maxIters, double relTol) {
    const int n = static_cast<int>(b.size());
    double bNorm = 0.0;
    for (double v : b) bNorm += v * v;
    bNorm = std::sqrt(bNorm);
    if (bNorm < 1.0e-300) {
        std::fill(x.begin(), x.end(), 0.0);
        return 0;
    }

    int totalMatVecs = 0;
    while (totalMatVecs < maxIters) {
        std::vector<double> Ax(n);
        matVec(x, Ax);
        ++totalMatVecs;
        std::vector<double> r(n);
        for (int i = 0; i < n; ++i) r[i] = b[i] - Ax[i];
        double rNorm = 0.0;
        for (double v : r) rNorm += v * v;
        rNorm = std::sqrt(rNorm);
        if (rNorm <= relTol * bNorm) return totalMatVecs;

        const int m = std::min(restart, maxIters - totalMatVecs);
        if (m <= 0) break;

        std::vector<std::vector<double>> V(m + 1, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> H(m + 1, std::vector<double>(m, 0.0));
        std::vector<double> cs(m, 0.0), sn(m, 0.0), g(m + 1, 0.0);
        for (int i = 0; i < n; ++i) V[0][i] = r[i] / rNorm;
        g[0] = rNorm;

        int k = 0;
        for (; k < m; ++k) {
            std::vector<double> w(n);
            matVec(V[k], w);
            ++totalMatVecs;
            for (int i = 0; i <= k; ++i) {
                double hik = 0.0;
                for (int j = 0; j < n; ++j) hik += w[j] * V[i][j];
                H[i][k] = hik;
                for (int j = 0; j < n; ++j) w[j] -= hik * V[i][j];
            }
            double wNorm = 0.0;
            for (double v : w) wNorm += v * v;
            wNorm = std::sqrt(wNorm);
            H[k + 1][k] = wNorm;
            if (wNorm > 1.0e-300) {
                for (int j = 0; j < n; ++j) V[k + 1][j] = w[j] / wNorm;
            }
            for (int i = 0; i < k; ++i) {
                const double temp = cs[i] * H[i][k] + sn[i] * H[i + 1][k];
                H[i + 1][k] = -sn[i] * H[i][k] + cs[i] * H[i + 1][k];
                H[i][k] = temp;
            }
            const double denom = std::sqrt(H[k][k] * H[k][k] + H[k + 1][k] * H[k + 1][k]);
            if (denom > 1.0e-300) {
                cs[k] = H[k][k] / denom;
                sn[k] = H[k + 1][k] / denom;
            } else {
                cs[k] = 1.0;
                sn[k] = 0.0;
            }
            H[k][k] = cs[k] * H[k][k] + sn[k] * H[k + 1][k];
            H[k + 1][k] = 0.0;
            const double gTemp = cs[k] * g[k];
            g[k + 1] = -sn[k] * g[k];
            g[k] = gTemp;

            const double resid = std::fabs(g[k + 1]);
            ++k; // k now counts completed Arnoldi steps (1-indexed count)
            if (resid <= relTol * bNorm || totalMatVecs >= maxIters) break;
        }

        std::vector<double> y(k, 0.0);
        for (int i = k - 1; i >= 0; --i) {
            double sum = g[i];
            for (int j = i + 1; j < k; ++j) sum -= H[i][j] * y[j];
            y[i] = (std::fabs(H[i][i]) > 1.0e-300) ? sum / H[i][i] : 0.0;
        }
        for (int i = 0; i < k; ++i)
            for (int j = 0; j < n; ++j) x[j] += y[i] * V[i][j];
    }
    return totalMatVecs;
}

} // namespace

FourFieldSolver::ImplicitStepResult FourFieldSolver::stepImplicitPressureVelocity(double dt) {
    ImplicitStepResult result;
    const int N = state_.N;
    const int nUnk = 3 * N - 2;
    if (N < 2) return result;

    computeGeometry(); // geometry + h1Eq_ (if enabled) from the CURRENT (old) volume fractions
    computeClosures();  // tauW1_/tauW2_/tauI_/Ue_/Ud_/phiE_/phiDe_ at the OLD velocity field --
                         // only the mass-transfer rates (Ue_, Ud_, phiE_, phiDe_) are actually
                         // used below (lagged, per the scoping note above); tauW1_/tauW2_/tauI_
                         // are deliberately NOT used -- friction is recomputed at the trial
                         // velocity inside the residual evaluator instead.

    const double A = geometry_.area();
    const double D = geometry_.diameter();
    const double rhoL = fluid_.rhoLiquid;

    // Snapshot everything this residual needs from the OLD state. Captured
    // by value (not reference to state_) so later trial evaluations, which
    // never touch state_ itself, cannot see their own side effects.
    const std::vector<double> u1Old = state_.u1, u2Old = state_.u2, POld = state_.P;
    const std::vector<double> rhoGasOld = rhoGasCell_;
    const std::vector<double> elOld = state_.el, edOld = state_.ed, egOld = state_.eg, ebOld = state_.eb;
    const std::vector<double> ulOld = state_.ul, ugOld = state_.ug, udOld = state_.ud, ubOld = state_.ub;
    const auto geomOld = geom_;
    const double outletPressure = bc_.outletPressure;
    const double maxV = options_.maxVelocity;

    const double inletEl_ = inletEl(), inletEd_ = inletEd(), inletEg_ = inletEg(), inletEb_ = inletEb();
    const double inletEL_ = inletEL(), inletEG_ = inletEG();
    const double ul0 = bc_.inletSuperficialLiquid / std::max(inletEL_, small_e);
    const double ug0 = bc_.inletSuperficialGas / std::max(inletEG_, small_e);

    // Mixing-length diffusion term (SolverOptions::enableTurbulentViscosity),
    // if enabled: lagged, built once from the OLD (u2-u1) field, exactly
    // like massSrc1/massSrc2 below -- out of scope for this first JFNK
    // implementation to make implicit in its own right (see
    // FourFieldSolver.hpp's scoping note).
    std::vector<double> turb1(N + 1, 0.0), turb2(N + 1, 0.0);
    if (options_.enableTurbulentViscosity) {
        const double mixingLength = options_.turbulentMixingLengthFraction * D;
        std::vector<double> urOld(N + 1);
        for (int i = 0; i <= N; ++i) urOld[i] = u2Old[i] - u1Old[i];
        for (int f = 1; f < N; ++f) {
            const double dzL = state_.cellWidth(f - 1), dzR = state_.cellWidth(f);
            const double d2 = 2.0 * ((urOld[f + 1] - urOld[f]) / (dzR * (dzL + dzR)) -
                                      (urOld[f] - urOld[f - 1]) / (dzL * (dzL + dzR)));
            const double nuT = mixingLength * std::fabs(urOld[f]);
            const double aDiff = nuT * d2;
            const double e1fSafe = std::max(0.5 * (state_.e1(f - 1) + state_.e1(f)), options_.momentumFractionFloor);
            const double e2fSafe = std::max(0.5 * (state_.e2(f - 1) + state_.e2(f)), options_.momentumFractionFloor);
            const double rho1f = 0.5 * (rho1Of(state_, rhoGasCell_, rhoL, f - 1) + rho1Of(state_, rhoGasCell_, rhoL, f));
            const double rho2f = 0.5 * (rho2Of(state_, rhoGasCell_, rhoL, f - 1) + rho2Of(state_, rhoGasCell_, rhoL, f));
            const double m1 = rho1f * e1fSafe, m2 = rho2f * e2fSafe;
            turb1[f] = -aDiff * m2 / (m1 + m2);
            turb2[f] = +aDiff * m1 / (m1 + m2);
        }
    }

    // Lagged (old-state) dimensionless friction FACTORS, precomputed once
    // per face -- NOT recomputed at the trial velocity inside the Newton
    // iteration. These correlations (wallFrictionFactor,
    // interfacialFrictionFactor) switch branch at Reynolds-number and
    // Kelvin-Helmholtz-F thresholds (e.g. the Andreussi-Persen F0=0.36
    // enhanced-friction onset), making them piecewise, non-smooth
    // functions of velocity -- confirmed by direct test to stall Newton's
    // convergence outright when left inside the trial-state evaluation (a
    // tiny step along the steepest-descent direction was found to
    // increase the residual by 3x, and a step 100x smaller than that by
    // over 1000x, the signature of a sharp kink, not just a steep but
    // smooth landscape). Freezing the FACTOR at its old-state value while
    // still recomputing the resulting shear stress at the TRIAL velocity
    // (wallShearStress/interfacialShearStress are themselves smooth,
    // quadratic |u|*u functions with no branches) keeps the dominant,
    // stiffness-causing nonlinearity genuinely implicit -- the actual
    // target of this method -- while removing the correlation-switch
    // discontinuity from the Newton iteration entirely. This mirrors
    // standard practice in real system codes (RELAP/TRAC/CATHARE-class),
    // which likewise re-evaluate friction CORRELATIONS only once per time
    // step even under implicit velocity solves.
    std::vector<double> flwLagged(N + 1, 0.0), fgwLagged(N + 1, 0.0), fiLagged(N + 1, 0.0);
    {
        const double epsVel = 1.0e-3;
        for (int f = 1; f < N; ++f) {
            const double rhoGf = 0.5 * (rhoGasOld[f - 1] + rhoGasOld[f]);
            const double u1o = u1Old[f], u2o = u2Old[f], uro = u2Old[f] - u1Old[f];
            flwLagged[f] = (std::fabs(u1o) > epsVel) ? tauW1_[f] / (0.5 * rhoL * std::fabs(u1o) * u1o) : 0.0;
            fgwLagged[f] = (std::fabs(u2o) > epsVel) ? tauW2_[f] / (0.5 * rhoGf * std::fabs(u2o) * u2o) : 0.0;
            fiLagged[f] = (std::fabs(uro) > epsVel) ? tauI_[f] / (0.5 * rhoGf * std::fabs(uro) * uro) : 0.0;
        }
    }

    // --- Residual evaluator: x (trial) -> R. Pure function of its
    // arguments and the OLD-state snapshots captured above; never touches
    // state_/rhoGasCell_/geom_ (those still hold the OLD-state values
    // throughout, read-only, as snapshotted/cached before this point).
    auto evalResidual = [&](const std::vector<double>& x, std::vector<double>& R) {
        auto u1At = [&](int f) -> double {
            if (f == 0) return ul0;
            if (f == N) return x[N - 2]; // zero-gradient outlet: equals last interior face
            return x[f - 1];
        };
        auto u2At = [&](int f) -> double {
            if (f == 0) return ug0;
            if (f == N) return x[(N - 1) + (N - 2)];
            return x[(N - 1) + (f - 1)];
        };
        auto pAt = [&](int i) -> double { return x[2 * (N - 1) + i]; };

        R.assign(nUnk, 0.0);

        // --- Momentum residuals, interior faces f=1..N-1 ---
        for (int f = 1; f < N; ++f) {
            const int cL = f - 1, cR = f;
            const double e1L = state_.e1(cL), e1R = state_.e1(cR);
            const double e2L = state_.e2(cL), e2R = state_.e2(cR);
            const double e1fSafe = std::max(0.5 * (e1L + e1R), options_.momentumFractionFloor);
            const double e2fSafe = std::max(0.5 * (e2L + e2R), options_.momentumFractionFloor);
            const double rho1f = 0.5 * (rho1Of(state_, rhoGasCell_, rhoL, cL) + rho1Of(state_, rhoGasCell_, rhoL, cR));
            const double rho2f = 0.5 * (rho2Of(state_, rhoGasCell_, rhoL, cL) + rho2Of(state_, rhoGasCell_, rhoL, cR));
            const double rhoGf = 0.5 * (rhoGasOld[cL] + rhoGasOld[cR]);

            const double centerDz = state_.centerDistance(cL, cR);
            const double thetaF = 0.5 * (state_.theta[cL] + state_.theta[cR]);
            const double Swp1f = 0.5 * (geomOld[cL].Swp1 + geomOld[cR].Swp1);
            const double Swp2f = 0.5 * (geomOld[cL].Swp2 + geomOld[cR].Swp2);
            const double Sif = 0.5 * (geomOld[cL].Si + geomOld[cR].Si);

            // Pressure gradient and gravity term: trial P, lagged geometry.
            const double dPdz = (pAt(cR) - pAt(cL)) / centerDz;
            double gravCosTerm, gravSinTerm;
            if (!options_.enableWellBalancedGravity) {
                const double dh1dz = (geomOld[cR].h1 - geomOld[cL].h1) / centerDz;
                gravCosTerm = gravity * std::cos(thetaF) * dh1dz;
                gravSinTerm = gravity * std::sin(thetaF);
            } else {
                const double detadz = ((geomOld[cR].h1 - h1Eq_[cR]) - (geomOld[cL].h1 - h1Eq_[cL])) / centerDz;
                gravCosTerm = gravity * std::cos(thetaF) * detadz;
                gravSinTerm = 0.0;
            }

            // Advection: upwind on the TRIAL velocity field (same formula
            // as updateLayerMomentum's non-flux-form branch).
            const double u1f = u1At(f);
            const double du1dzAdv = (u1f >= 0.0) ? (u1At(f) - u1At(f - 1)) / state_.cellWidth(cL)
                                                  : (u1At(f + 1) - u1At(f)) / state_.cellWidth(cR);
            const double advectiveTerm1 = u1f * du1dzAdv;
            const double u2f = u2At(f);
            const double du2dzAdv = (u2f >= 0.0) ? (u2At(f) - u2At(f - 1)) / state_.cellWidth(cL)
                                                  : (u2At(f + 1) - u2At(f)) / state_.cellWidth(cR);
            const double advectiveTerm2 = u2f * du2dzAdv;

            // Friction: dimensionless factors LAGGED (precomputed above,
            // from the old state); the resulting shear stress is still
            // recomputed at the TRIAL velocity -- a genuinely nonlinear
            // (quadratic, smooth) implicit treatment of the dominant
            // |velocity|*velocity dependence, without the correlation's
            // own branch discontinuities inside the Newton iteration. See
            // the comment on flwLagged/fgwLagged/fiLagged above.
            const double tauW1 = wallShearStress(flwLagged[f], rhoL, u1f);
            const double tauW2 = wallShearStress(fgwLagged[f], rhoGf, u2f);
            const double tauI = interfacialShearStress(fiLagged[f], rhoGf, u2f, u1f);

            // Mass-transfer source terms: lagged (OLD entrainment/
            // deposition rates and OLD continuous/dispersed velocities),
            // same simplification as the turbulent-viscosity term above.
            const double UeF = 0.5 * (Ue_[cL] + Ue_[cR]);
            const double UdF = 0.5 * (Ud_[cL] + Ud_[cR]);
            const double phiEF = 0.5 * (phiE_[cL] + phiE_[cR]);
            const double phiDeF = 0.5 * (phiDe_[cL] + phiDe_[cR]);
            const double massSrc1 = (-UeF * ulOld[f] + UdF * udOld[f] + phiEF * ugOld[f] - phiDeF * ubOld[f]) / (e1fSafe * rho1f);
            const double massSrc2 = (UeF * ulOld[f] - UdF * udOld[f] - phiEF * ugOld[f] + phiDeF * ubOld[f]) / (e2fSafe * rho2f);

            const double du1dt = -advectiveTerm1 - dPdz / rho1f - gravCosTerm - gravSinTerm
                                  - tauW1 * Swp1f / (A * e1fSafe * rho1f)
                                  + tauI * Sif / (A * e1fSafe * rho1f)
                                  + massSrc1 + turb1[f];
            const double du2dt = -advectiveTerm2 - dPdz / rho2f - gravCosTerm - gravSinTerm
                                  - tauW2 * Swp2f / (A * e2fSafe * rho2f)
                                  - tauI * Sif / (A * e2fSafe * rho2f)
                                  + massSrc2 + turb2[f];

            R[f - 1] = (u1f - u1Old[f]) - dt * du1dt;
            R[(N - 1) + (f - 1)] = (u2f - u2Old[f]) - dt * du2dt;
        }

        // --- Pressure / mass-conservation residuals, cells 0..N-1 ---
        auto rho1At = [&](int i) { return rho1Of(state_, rhoGasCell_, rhoL, i); };
        auto rho2At = [&](int i) { return rho2Of(state_, rhoGasCell_, rhoL, i); };
        std::vector<double> Rp(N, 0.0);

        for (int f = 1; f < N; ++f) {
            const int cL = f - 1, cR = f;
            const double u1f = u1At(f), u2f = u2At(f);
            // Back-substitute continuous-phase velocities from the trial
            // (u1,u2) using the SAME formula as backSubstitutePhaseVelocities(),
            // with volume fractions/densities lagged at the OLD state.
            const double el = 0.5 * (elOld[cL] + elOld[cR]);
            const double eb = 0.5 * (ebOld[cL] + ebOld[cR]);
            const double eg = 0.5 * (egOld[cL] + egOld[cR]);
            const double ed = 0.5 * (edOld[cL] + edOld[cR]);
            const double e1 = std::max(el + eb, small_e);
            const double e2 = std::max(eg + ed, small_e);
            const double rhoG = 0.5 * (rhoGasOld[cL] + rhoGasOld[cR]);
            const double rho1 = (el * rhoL + eb * rhoG) / e1;
            const double rho2 = (eg * rhoG + ed * rhoL) / e2;
            const double cb = std::clamp(eb / e1, 0.0, 1.0);
            const double cd = std::clamp(ed / e2, 0.0, 1.0);
            const double floor = options_.momentumFractionFloor;
            const double denomL = std::max((1.0 - cb) * rhoL, floor * rhoL);
            const double ulf = (rho1 * u1f - cb * rhoG * ubOld[f]) / denomL;
            const double denomG = std::max((1.0 - cd) * rhoG, floor * rhoG);
            const double ugf = (rho2 * u2f - cd * rhoL * udOld[f]) / denomG;

            const double elF = (ulf >= 0.0) ? elOld[cL] : elOld[cR];
            const double edF = (udOld[f] >= 0.0) ? edOld[cL] : edOld[cR];
            const double egF = (ugf >= 0.0) ? egOld[cL] : egOld[cR];
            const double ebF = (ubOld[f] >= 0.0) ? ebOld[cL] : ebOld[cR];
            const double qTotal = elF * ulf + edF * udOld[f] + egF * ugf + ebF * ubOld[f];

            const double widthL = state_.cellWidth(cL), widthR = state_.cellWidth(cR);
            Rp[cL] += qTotal / widthL;
            Rp[cR] -= qTotal / widthR;
        }
        {
            const double qLiquidIn = inletEl_ * ul0 + inletEd_ * udOld[0];
            const double qGasIn = inletEg_ * ug0 + inletEb_ * ubOld[0];
            Rp[0] -= (qLiquidIn + qGasIn) / state_.cellWidth(0);
        }
        {
            const int c = N - 1;
            const double u1N = u1At(N), u2N = u2At(N);
            const double rhoG = rhoGasOld[c];
            const double rho1 = (elOld[c] * rhoL + ebOld[c] * rhoG) / std::max(elOld[c] + ebOld[c], small_e);
            const double rho2 = (egOld[c] * rhoG + edOld[c] * rhoL) / std::max(egOld[c] + edOld[c], small_e);
            const double cb = std::clamp(ebOld[c] / std::max(elOld[c] + ebOld[c], small_e), 0.0, 1.0);
            const double cd = std::clamp(edOld[c] / std::max(egOld[c] + edOld[c], small_e), 0.0, 1.0);
            const double floor = options_.momentumFractionFloor;
            const double ulN = (rho1 * u1N - cb * rhoG * ubOld[N]) / std::max((1.0 - cb) * rhoL, floor * rhoL);
            const double ugN = (rho2 * u2N - cd * rhoL * udOld[N]) / std::max((1.0 - cd) * rhoG, floor * rhoG);
            const double qLiquidOut = elOld[c] * ulN + edOld[c] * udOld[N];
            const double qGasOut = egOld[c] * ugN + ebOld[c] * ubOld[N];
            Rp[c] += (qLiquidOut + qGasOut) / state_.cellWidth(c);
        }

        for (int i = 0; i < N - 1; ++i) {
            const double eG = egOld[i] + ebOld[i];
            const double rhoGasNew = fluid_.rhoGas(pAt(i));
            const double accumulation = eG * (rhoGasNew - rhoGasOld[i]) / (std::max(rhoGasOld[i], tiny) * dt);
            R[2 * (N - 1) + i] = Rp[i] + accumulation;
        }
        // Outlet cell: Dirichlet pressure residual, not a flux balance,
        // matching solvePressureCorrection's own hard-pin of this cell.
        R[2 * (N - 1) + (N - 1)] = pAt(N - 1) - outletPressure;
        (void)rho1At; (void)rho2At; // kept for symmetry/readability; not otherwise needed
    };

    // --- Newton loop, operating on a DIAGONALLY SCALED unknown/residual,
    // not the raw physical (u1,u2,P) vector. ---
    //
    // The physical unknowns span wildly different magnitudes (velocities
    // O(0.1-10) m/s, pressures O(1e5-1e6) Pa), and the momentum and
    // pressure residual BLOCKS likewise have different natural units
    // (velocity, and inverse time respectively). Working directly in
    // physical units was tried first and FAILED outright -- confirmed by
    // direct test, not assumed: Newton's residual norm did not decrease
    // at all across 30 iterations on even a mild, non-stiff case, because
    // a single finite-difference step size applied uniformly across such
    // mismatched scales is simultaneously far too large for the velocity
    // components and (relatively) far too small to resolve the pressure
    // components' sensitivity, badly conditioning both the finite-
    // difference Jacobian estimate and GMRES's own orthogonalization
    // (confirmed innocent in isolation: see the standalone GMRES unit
    // test in this session's working notes, which solves several small
    // dense systems to near machine precision with the exact same GMRES
    // code). Diagonal scaling -- xHat = x / scale, with the residual
    // correspondingly divided by rscale -- is standard practice for
    // exactly this reason in any real multi-physics JFNK implementation.
    const double uScale = std::max(1.0e-3, [&] {
        double m = 0.0;
        for (double v : u1Old) m = std::max(m, std::fabs(v));
        for (double v : u2Old) m = std::max(m, std::fabs(v));
        return m;
    }());
    const double pScale = std::max(1.0, outletPressure);
    const double Lscale = state_.L / N;
    const double rScaleP = std::max(uScale / Lscale, 1.0 / std::max(dt, tiny));

    std::vector<double> scale(nUnk), rscale(nUnk);
    for (int i = 0; i < 2 * (N - 1); ++i) scale[i] = uScale;
    for (int i = 2 * (N - 1); i < nUnk; ++i) scale[i] = pScale;
    for (int i = 0; i < 2 * (N - 1); ++i) rscale[i] = uScale;
    for (int i = 2 * (N - 1); i < nUnk - 1; ++i) rscale[i] = rScaleP;
    rscale[nUnk - 1] = pScale; // outlet cell's Dirichlet-pressure residual: units of Pa, not 1/s

    auto evalScaled = [&](const std::vector<double>& xHat, std::vector<double>& Fhat) {
        std::vector<double> xPhys(nUnk);
        for (int i = 0; i < nUnk; ++i) xPhys[i] = xHat[i] * scale[i];
        std::vector<double> Rphys;
        evalResidual(xPhys, Rphys);
        Fhat.resize(nUnk);
        for (int i = 0; i < nUnk; ++i) Fhat[i] = Rphys[i] / rscale[i];
    };

    std::vector<double> xHat(nUnk);
    for (int f = 1; f < N; ++f) {
        xHat[f - 1] = u1Old[f] / uScale;
        xHat[(N - 1) + (f - 1)] = u2Old[f] / uScale;
    }
    for (int i = 0; i < N; ++i) xHat[2 * (N - 1) + i] = POld[i] / pScale;

    std::vector<double> F0;
    evalScaled(xHat, F0);
    double r0norm = 0.0;
    for (double v : F0) r0norm += v * v;
    r0norm = std::sqrt(r0norm);
    result.initialResidualNorm = r0norm;
    if (r0norm < 1.0e-300) {
        result.converged = true;
        result.finalResidualNorm = r0norm;
        return result; // already at a fixed point (e.g. dt effectively 0)
    }

    const double absTol = std::max(jfnk.newtonTol * r0norm, 1.0e-14);
    std::vector<double> Fcur = F0;
    double rNorm = r0norm;
    int newtonIter = 0;
    for (; newtonIter < jfnk.maxNewtonIters; ++newtonIter) {
        if (rNorm <= absTol) break;

        // Matrix-free Jacobian-vector product via forward finite
        // difference, now entirely in the well-scaled xHat/Fhat space, so
        // a single relative step size (jfnk.fdEpsilon) is meaningful
        // uniformly across all components.
        auto matVec = [&](const std::vector<double>& v, std::vector<double>& out) {
            double vNorm = 0.0;
            for (double c : v) vNorm += c * c;
            vNorm = std::sqrt(vNorm);
            if (vNorm < 1.0e-300) { std::fill(out.begin(), out.end(), 0.0); return; }
            const double h = jfnk.fdEpsilon * std::max(1.0, std::sqrt(std::inner_product(xHat.begin(), xHat.end(), xHat.begin(), 0.0)) / std::sqrt((double)nUnk)) / vNorm;
            std::vector<double> xPerturbed(nUnk);
            for (int i = 0; i < nUnk; ++i) xPerturbed[i] = xHat[i] + h * v[i];
            std::vector<double> Fp;
            evalScaled(xPerturbed, Fp);
            out.resize(nUnk);
            for (int i = 0; i < nUnk; ++i) out[i] = (Fp[i] - Fcur[i]) / h;
        };

        // A left Jacobi (diagonal) preconditioner was tried here first --
        // this is a coupled pressure-velocity ("saddle-point") system, a
        // class well known to be badly conditioned for plain Krylov
        // methods -- but confirmed by direct test to make no measurable
        // difference to the stall described below, while costing an
        // extra nUnk residual evaluations per Newton iteration; removed
        // rather than kept as dead weight. What DOES fix the stall is
        // addressed one level up, in
        // stepImplicitPressureVelocityAdaptive()'s pseudo-transient
        // continuation (sub-stepping) -- see its own documentation for
        // the full diagnosis (the direct, undamped step was shown by
        // further testing to stall not from ill-conditioning but from
        // Newton's convergence basin shrinking at larger dt, confirmed by
        // the SAME case converging to machine precision in 2 iterations
        // at a much smaller dt).
        std::vector<double> negF(nUnk);
        for (int i = 0; i < nUnk; ++i) negF[i] = -Fcur[i];
        std::vector<double> dxHat(nUnk, 0.0);
        const int gm = gmresSolve(matVec, negF, dxHat, jfnk.gmresRestart,
                                   jfnk.gmresMaxIters, jfnk.gmresTol);
        result.totalGmresIterations += gm;

        // Lightweight opt-in trace (set MFS_JFNK_DEBUG), same convention as
        // FourFieldSolver.cpp's own MFS_DEBUG -- zero overhead when unset.
        if (std::getenv("MFS_JFNK_DEBUG")) {
            double dxNorm = 0.0;
            for (double v : dxHat) dxNorm += v * v;
            dxNorm = std::sqrt(dxNorm);
            std::fprintf(stderr, "[jfnk] newtonIter=%d rNorm=%.6e gmresIters=%d ||dxHat||=%.6e\n",
                         newtonIter, rNorm, gm, dxNorm);
        }

        // Simple damping/globalization: halve the step if it fails to
        // reduce the residual, up to a handful of times, rather than
        // accepting a Newton step that makes things worse (plain,
        // standard practice for a first robust JFNK implementation --
        // more sophisticated line search is a possible future refinement).
        double lambda = 1.0;
        std::vector<double> xTrial(nUnk), Ftrial;
        double rTrialNorm = rNorm;
        for (int ls = 0; ls < 8; ++ls) {
            for (int i = 0; i < nUnk; ++i) xTrial[i] = xHat[i] + lambda * dxHat[i];
            evalScaled(xTrial, Ftrial);
            rTrialNorm = 0.0;
            for (double v : Ftrial) rTrialNorm += v * v;
            rTrialNorm = std::sqrt(rTrialNorm);
            if (rTrialNorm < rNorm || lambda < 1.0 / 64.0) break;
            lambda *= 0.5;
        }
        xHat = xTrial;
        Fcur = Ftrial;
        rNorm = rTrialNorm;
    }

    // Back to physical units for the convergence check / commit below.
    std::vector<double> x(nUnk);
    for (int i = 0; i < nUnk; ++i) x[i] = xHat[i] * scale[i];

    result.newtonIterations = newtonIter;
    result.finalResidualNorm = rNorm;
    result.converged = (rNorm <= absTol);
    if (!result.converged) return result; // state_ left untouched, as documented

    // Commit: unpack x into state_, clamp velocities (same policy as
    // step()'s own momentum update), then run the existing (unmodified)
    // dispersed-momentum and continuity stages to finish the step exactly
    // as step() itself would after its own momentum/pressure update.
    for (int f = 1; f < N; ++f) {
        state_.u1[f] = std::clamp(x[f - 1], -maxV, maxV);
        state_.u2[f] = std::clamp(x[(N - 1) + (f - 1)], -maxV, maxV);
    }
    state_.u1[0] = ul0; state_.u2[0] = ug0;
    state_.u1[N] = state_.u1[N - 1]; state_.u2[N] = state_.u2[N - 1];
    for (int i = 0; i < N; ++i) {
        const double pMax = options_.maxPressureFactor * outletPressure;
        state_.P[i] = std::clamp(x[2 * (N - 1) + i], options_.minPressure, pMax);
    }
    state_.P[N - 1] = outletPressure;
    for (int i = 0; i < N; ++i) rhoGasCell_[i] = fluid_.rhoGas(state_.P[i]);

    backSubstitutePhaseVelocities();
    applyInletBoundary();
    applyOutletBoundary();
    updateDispersedMomentum(dt);
    backSubstitutePhaseVelocities();
    updateContinuity(dt);

    time_ += dt;
    ++stepCount_;
    return result;
}

FourFieldSolver::ImplicitAdaptiveResult FourFieldSolver::stepImplicitPressureVelocityAdaptive(double dt) {
    ImplicitAdaptiveResult agg;
    double remaining = dt;
    double trial = dt;
    while (remaining > 1.0e-15) {
        trial = std::min(trial, remaining);
        if (trial < options_.minTimeStep) {
            return agg; // give up: leaves the solver at the end of the last converged sub-step
        }
        const auto r = stepImplicitPressureVelocity(trial);
        agg.totalNewtonIterations += r.newtonIterations;
        agg.totalGmresIterations += r.totalGmresIterations;
        if (r.converged) {
            ++agg.subSteps;
            remaining -= trial;
            agg.dtCovered += trial;
            trial *= 2.0; // grow back toward the target once Newton is converging easily
        } else {
            trial *= 0.5; // shrink and retry at this same point in time
        }
    }
    agg.converged = true;
    return agg;
}

} // namespace mfs
