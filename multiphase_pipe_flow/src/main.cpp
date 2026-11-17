// Command-line driver for the four-field near-horizontal-to-vertical
// gas-liquid pipe flow solver (mfs::FourFieldSolver), implementing the
// model of Bonizzi, Andreussi & Banerjee (2009), Int. J. Multiphase Flow
// 35, 34-46, generalised to arbitrary pipe inclination.
//
// Usage:
//   mfs_demo horizontal   Reproduces the qualitative Section 4/5.2 setup:
//                         horizontal stratified air/oil flow with a
//                         liquid-rate step increase that triggers slugging.
//   mfs_demo vertical     A vertical (theta=90 deg) bubbly-flow case,
//                         demonstrating the same equations applied outside
//                         their originally validated near-horizontal range.
//   mfs_demo terrain      A terrain-following pipe (inclination varies
//                         along z), the scenario the paper's introduction
//                         cites as its core motivation.
//   mfs_demo horizontal_amr
//                         The same horizontal case, at a fixed resolution
//                         and with adaptive mesh refinement
//                         (SolverOptions::amr) enabled under each of its
//                         two refine/coarsen indicators, timed side by
//                         side. The default, `KelvinHelmholtz`, refines
//                         using the same F parameter already computed for
//                         the interfacial friction closure (Closures.hpp)
//                         -- a PHYSICAL instability criterion, coarsening
//                         cells where the interface is quiescent and never
//                         refining finer than the baseline resolution by
//                         default. `WaveletDetail` instead refines/
//                         coarsens on a purely NUMERICAL representability
//                         criterion in the adaptive-multiresolution/
//                         wavelet-threshold tradition (Harten, 1994; Cohen,
//                         Kaber, Müller & Postel, 2003) -- see
//                         FourFieldSolver.hpp, SolverOptions::amr, and
//                         FourFieldSolver.cpp's computeWaveletDetail().
//   mfs_demo horizontal_movingmesh
//                         The same horizontal case, fixed grid vs. a moving
//                         (r-adaptive) mesh that keeps the cell count fixed
//                         and continuously relocates all N points toward
//                         wherever a monitor function is largest (see
//                         SolverOptions::movingMesh). Unlike AMR this is
//                         NOT about speed (it costs more per step, since
//                         the whole field is remapped every relocation) --
//                         the point is a much sharper captured front,
//                         since resolution isn't bounded by any minimum
//                         cell width.
//   mfs_demo horizontal_limiter
//                         The same horizontal case, on the SAME fixed
//                         N=300 grid throughout, comparing first-order
//                         upwind against MUSCL/TVD flux-limited advection
//                         (SolverOptions::advectionLimiter: Minmod, Van
//                         Leer, Superbee) in the field continuity
//                         equations. Reduces numerical diffusion of a
//                         forming front without touching the mesh itself.
//   mfs_demo mass_conservation
//                         A Kelvin-Helmholtz case study (Brito 2012, the
//                         same one used elsewhere in this project's
//                         regularization checks), tracking a mass-error
//                         ratio (Kim, Lee & Kim 2022's own methodology)
//                         against time, across three mesh resolutions,
//                         with SolverOptions::enableSurfaceTension off vs.
//                         on -- verifying that the biharmonic term's
//                         eL/eG "equal and opposite" construction adds no
//                         measurable mass leak beyond the baseline
//                         scheme's own. See VALIDATION.md, item 27.
//   mfs_demo wang2019_friction
//                         This solver's own self-consistent predicted
//                         gas-wall/liquid-wall friction factor (Re_G/f_G,
//                         Re_L/f_L, from SolverOptions::wallCorrelation's
//                         Taitel & Dukler 1976 closure evaluated at the
//                         solver's OWN predicted holdup, not an externally
//                         imposed one) over a (Vsl,Vsg) sweep at atmosphere
//                         and 2 MPa, for comparison against the digitized
//                         experimental data from Wang, Liu, Chang, Zhao &
//                         Guo (2019), Int. J. Heat Mass Transfer 143:118537.
//
// Each case writes a CSV time-series of cell-centred fields to
// "<case>_output.csv" (one block of rows per snapshot) and prints a
// flow-regime summary and mass-balance diagnostic to stdout.

#include "mfs/Constants.hpp"
#include "mfs/FlowRegimeClassifier.hpp"
#include "mfs/FourFieldSolver.hpp"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace {

void writeSnapshotHeader(std::ofstream& out) {
    out << "time,z,eL,eG,el,ed,eg,eb,P,u1,u2,ul,ug,ud,ub,theta_deg\n";
}

void writeSnapshot(std::ofstream& out, const mfs::FourFieldSolver& solver) {
    const auto& s = solver.state();
    for (int i = 0; i < s.N; ++i) {
        const double z = s.cellCenter(i);
        const double u1c = 0.5 * (s.u1[i] + s.u1[i + 1]);
        const double u2c = 0.5 * (s.u2[i] + s.u2[i + 1]);
        const double ulc = 0.5 * (s.ul[i] + s.ul[i + 1]);
        const double ugc = 0.5 * (s.ug[i] + s.ug[i + 1]);
        const double udc = 0.5 * (s.ud[i] + s.ud[i + 1]);
        const double ubc = 0.5 * (s.ub[i] + s.ub[i + 1]);
        out << solver.time() << ',' << z << ',' << s.eL(i) << ',' << s.eG(i) << ','
            << s.el[i] << ',' << s.ed[i] << ',' << s.eg[i] << ',' << s.eb[i] << ','
            << s.P[i] << ',' << u1c << ',' << u2c << ',' << ulc << ',' << ugc << ','
            << udc << ',' << ubc << ',' << s.theta[i] * 180.0 / mfs::constants::pi << '\n';
    }
}

void printRegimeSummary(const mfs::FourFieldSolver& solver) {
    const auto& s = solver.state();
    // Crude single-snapshot fluctuation estimate from the |d(eg)/dz|
    // between neighbouring cells, used only for this printed summary.
    std::vector<double> fluct(s.N, 0.0);
    for (int i = 0; i < s.N; ++i) {
        const double left = s.eg[std::max(i - 1, 0)];
        const double right = s.eg[std::min(i + 1, s.N - 1)];
        fluct[i] = std::fabs(right - left);
    }
    const auto regimes = mfs::classifyInstant(s.eg, fluct);
    int counts[4] = {0, 0, 0, 0};
    for (auto r : regimes) counts[static_cast<int>(r)]++;
    std::cout << "  regime cell counts -> stratified=" << counts[0] << " annular=" << counts[1]
              << " slug=" << counts[2] << " bubbly=" << counts[3] << "\n";
}

void printMassBalance(const mfs::FourFieldSolver& solver) {
    const auto mf = solver.boundaryMassFluxes();
    std::cout << "  liquid mass in pipe = " << solver.totalLiquidMass()
              << " kg, boundary flux in=" << mf.liquidIn << " out=" << mf.liquidOut << " kg/s\n";
    std::cout << "  gas    mass in pipe = " << solver.totalGasMass()
              << " kg, boundary flux in=" << mf.gasIn << " out=" << mf.gasOut << " kg/s\n";
}

mfs::FluidProperties airWater() {
    mfs::FluidProperties f;
    f.rhoLiquid = 998.0;
    f.muLiquid = 1.0e-3;
    f.sigma = 0.072;
    f.gasConstant = 287.0;
    f.temperature = 293.0;
    f.muGas = 1.8e-5;
    return f;
}

// useAmr: enables SolverOptions::amr at the (default) conservative setting
// -- coarsen quiescent regions, never refine finer than the baseline N.
// useMovingMesh: enables SolverOptions::movingMesh instead (mutually
// exclusive with useAmr by convention -- see FourFieldSolver.hpp).
// quiet: suppress the per-snapshot stdout diagnostics (used when this is
// called back-to-back for the comparison runs below).
// Returns stats for those comparisons, including the steepest captured
// |d(eL)/dz| at the end of the run -- the moving-mesh comparison's point
// is front sharpness (reduced numerical diffusion), not cell count or speed.
struct HorizontalRunStats { int finalN; int steps; double maxHoldupGradient; };

HorizontalRunStats runHorizontalCase(bool useAmr, bool useMovingMesh, const std::string& outputPath,
                                      bool quiet = false,
                                      mfs::FluxLimiterType limiter = mfs::FluxLimiterType::None,
                                      mfs::SolverOptions::AdaptiveMeshOptions::RefinementIndicator amrIndicator =
                                          mfs::SolverOptions::AdaptiveMeshOptions::RefinementIndicator::KelvinHelmholtz) {
    if (!quiet) {
        std::cout << "=== Horizontal stratified -> slug transient (cf. Section 4/5.2)"
                   << (useAmr ? ", with AMR ===\n" : useMovingMesh ? ", with moving mesh ===\n" : " ===\n");
    }
    const double D = 0.08, L = 30.0;
    const int N = 300; // ~1 diameter per cell, matching the paper's grid-independence finding (Fig. 2)

    mfs::SolverOptions opt;
    opt.amr.enabled = useAmr;
    opt.amr.indicator = amrIndicator;
    opt.movingMesh.enabled = useMovingMesh;
    opt.advectionLimiter = limiter;
    mfs::FourFieldSolver solver(D, L, N, airWater(), opt);
    solver.setInclinationConstant(0.0);

    mfs::BoundaryConditions bc;
    bc.inletSuperficialLiquid = 0.2;
    bc.inletSuperficialGas = 3.0;
    bc.inletLiquidHoldup = 0.15;
    bc.outletPressure = 1.0e5;
    bc.seedDisturbance = true;
    bc.disturbanceAmplitude = 0.04;
    bc.disturbanceFrequency = 0.7;
    solver.setBoundaryConditions(bc);
    solver.initializeStratified(0.15);

    std::ofstream out(outputPath);
    writeSnapshotHeader(out);

    const double tRampStart = 3.0;
    const double tEnd = 12.0;
    const double snapshotInterval = 1.0;
    double nextSnapshot = 0.0;
    int steps = 0;

    while (solver.time() < tEnd) {
        // Liquid-rate ramp, mirroring the paper's Section 4 transient
        // (0.2 -> 0.5 m/s), rescaled in time for a short demo run.
        if (solver.time() > tRampStart) {
            mfs::BoundaryConditions bcNow = bc;
            const double frac = std::min(1.0, (solver.time() - tRampStart) / 2.0);
            bcNow.inletSuperficialLiquid = 0.2 + frac * (0.5 - 0.2);
            solver.setBoundaryConditions(bcNow);
        }

        const double dt = solver.stableTimeStep();
        solver.step(dt);
        ++steps;

        if (solver.time() >= nextSnapshot) {
            writeSnapshot(out, solver);
            if (!quiet) {
                std::cout << "t=" << solver.time() << " s, N=" << solver.state().N << " cells\n";
                printRegimeSummary(solver);
                printMassBalance(solver);
            }
            nextSnapshot += snapshotInterval;
        }
    }
    double maxGrad = 0.0;
    const auto& s = solver.state();
    for (int i = 1; i < s.N; ++i) {
        maxGrad = std::max(maxGrad, std::fabs(s.eL(i) - s.eL(i - 1)) / s.centerDistance(i - 1, i));
    }
    if (!quiet) std::cout << "Wrote " << outputPath << "\n\n";
    return {s.N, steps, maxGrad};
}

void runHorizontalAmrComparison() {
    using Indicator = mfs::SolverOptions::AdaptiveMeshOptions::RefinementIndicator;
    std::cout << "=== AMR vs. fixed-grid comparison on the horizontal slug-formation case:\n"
                 "    the default Kelvin-Helmholtz indicator (a PHYSICAL instability\n"
                 "    criterion) against the WaveletDetail one (a NUMERICAL\n"
                 "    representability criterion, adaptive-multiresolution-style -- see\n"
                 "    SolverOptions::amr and FourFieldSolver.cpp's computeWaveletDetail()).\n"
                 "    ===\n";
    std::cout << "Running fixed grid (N=300)...\n";
    auto t0 = std::chrono::steady_clock::now();
    const auto fixedStats = runHorizontalCase(false, false, "horizontal_fixed_output.csv", /*quiet=*/true);
    auto t1 = std::chrono::steady_clock::now();

    std::cout << "Running with AMR enabled (Kelvin-Helmholtz indicator)...\n";
    const auto khStats = runHorizontalCase(true, false, "horizontal_amr_output.csv", /*quiet=*/true,
                                            mfs::FluxLimiterType::None, Indicator::KelvinHelmholtz);
    auto t2 = std::chrono::steady_clock::now();

    std::cout << "Running with AMR enabled (WaveletDetail indicator)...\n";
    const auto waveletStats = runHorizontalCase(true, false, "horizontal_amr_wavelet_output.csv", /*quiet=*/true,
                                                 mfs::FluxLimiterType::None, Indicator::WaveletDetail);
    auto t3 = std::chrono::steady_clock::now();

    const double fixedSec = std::chrono::duration<double>(t1 - t0).count();
    const double khSec = std::chrono::duration<double>(t2 - t1).count();
    const double waveletSec = std::chrono::duration<double>(t3 - t2).count();
    std::cout << "\nFixed grid     : " << fixedStats.steps << " steps, N=" << fixedStats.finalN
              << " cells throughout, " << fixedSec << " s, steepest |d(eL)/dz| = "
              << fixedStats.maxHoldupGradient << " /m\n";
    std::cout << "AMR (KH)       : " << khStats.steps << " steps, N=" << khStats.finalN
              << " cells at end, " << khSec << " s (" << (fixedSec / khSec)
              << "x), steepest |d(eL)/dz| = " << khStats.maxHoldupGradient << " /m ("
              << (khStats.maxHoldupGradient / fixedStats.maxHoldupGradient) << "x fixed)\n";
    std::cout << "AMR (Wavelet)  : " << waveletStats.steps << " steps, N=" << waveletStats.finalN
              << " cells at end, " << waveletSec << " s (" << (fixedSec / waveletSec)
              << "x), steepest |d(eL)/dz| = " << waveletStats.maxHoldupGradient << " /m ("
              << (waveletStats.maxHoldupGradient / fixedStats.maxHoldupGradient) << "x fixed, "
              << (waveletStats.maxHoldupGradient / khStats.maxHoldupGradient) << "x KH)\n";
    std::cout << "Wrote horizontal_fixed_output.csv, horizontal_amr_output.csv and\n"
                 "horizontal_amr_wavelet_output.csv (each AMR file's per-cell rows are at\n"
                 "whatever local resolution that snapshot's mesh had -- z spacing varies\n"
                 "row to row within a snapshot).\n\n";
}

void runHorizontalMovingMeshComparison() {
    std::cout << "=== Moving-mesh vs. fixed-grid comparison on the horizontal slug-formation\n"
                 "    case: this method's point is reduced numerical diffusion of a forming\n"
                 "    front (steeper captured gradient), NOT speed -- it relocates all N\n"
                 "    points and remaps the full field every step, so it costs MORE per step\n"
                 "    than the fixed grid, unlike AMR above. See README. ===\n";
    std::cout << "Running fixed grid (N=300)...\n";
    auto t0 = std::chrono::steady_clock::now();
    const auto fixedStats = runHorizontalCase(false, false, "horizontal_fixed2_output.csv", /*quiet=*/true);
    auto t1 = std::chrono::steady_clock::now();

    std::cout << "Running with moving mesh enabled (this is slower; relocates every step)...\n";
    const auto mmStats = runHorizontalCase(false, true, "horizontal_movingmesh_output.csv", /*quiet=*/true);
    auto t2 = std::chrono::steady_clock::now();

    const double fixedSec = std::chrono::duration<double>(t1 - t0).count();
    const double mmSec = std::chrono::duration<double>(t2 - t1).count();
    std::cout << "\nFixed grid   : " << fixedStats.steps << " steps, " << fixedSec
              << " s, steepest captured |d(eL)/dz| = " << fixedStats.maxHoldupGradient << " /m\n";
    std::cout << "Moving mesh  : " << mmStats.steps << " steps, " << mmSec << " s ("
              << (mmSec / fixedSec) << "x slower), steepest captured |d(eL)/dz| = "
              << mmStats.maxHoldupGradient << " /m ("
              << (mmStats.maxHoldupGradient / fixedStats.maxHoldupGradient) << "x sharper)\n";
    std::cout << "Wrote horizontal_fixed2_output.csv and horizontal_movingmesh_output.csv\n\n";
}

void runHorizontalLimiterComparison() {
    std::cout << "=== Flux-limiter (MUSCL/TVD) vs. first-order-upwind comparison on the\n"
                 "    horizontal slug-formation case: unlike the mesh-adaptivity demos above,\n"
                 "    this changes the RECONSTRUCTION used by both the field continuity\n"
                 "    equations (updateContinuity()) and the layer momentum equations' own\n"
                 "    self-advection term (updateLayerMomentum()) on the SAME fixed N=300\n"
                 "    grid -- the point is less numerical diffusion of a forming front, at a\n"
                 "    modest per-step cost. See README and\n"
                 "    SolverOptions::advectionLimiter. ===\n";

    struct Case { const char* name; mfs::FluxLimiterType limiter; const char* path; };
    const Case cases[] = {
        {"First-order upwind (None)", mfs::FluxLimiterType::None, "horizontal_limiter_none_output.csv"},
        {"Minmod", mfs::FluxLimiterType::Minmod, "horizontal_limiter_minmod_output.csv"},
        {"Van Leer", mfs::FluxLimiterType::VanLeer, "horizontal_limiter_vanleer_output.csv"},
        {"Superbee", mfs::FluxLimiterType::Superbee, "horizontal_limiter_superbee_output.csv"},
    };

    HorizontalRunStats noneStats{};
    double noneSec = 0.0;
    for (const auto& c : cases) {
        std::cout << "Running " << c.name << "...\n";
        auto t0 = std::chrono::steady_clock::now();
        const auto stats = runHorizontalCase(false, false, c.path, /*quiet=*/true, c.limiter);
        auto t1 = std::chrono::steady_clock::now();
        const double sec = std::chrono::duration<double>(t1 - t0).count();
        if (c.limiter == mfs::FluxLimiterType::None) { noneStats = stats; noneSec = sec; }
        std::cout << "  " << stats.steps << " steps, " << sec
                  << " s, steepest captured |d(eL)/dz| = " << stats.maxHoldupGradient << " /m";
        if (c.limiter != mfs::FluxLimiterType::None && noneSec > 0.0) {
            std::cout << " (" << (stats.maxHoldupGradient / noneStats.maxHoldupGradient)
                       << "x sharper than upwind, " << (sec / noneSec) << "x the cost)";
        }
        std::cout << "\n";
    }
    std::cout << "Wrote horizontal_limiter_{none,minmod,vanleer,superbee}_output.csv\n\n";
}

void runVerticalCase() {
    std::cout << "=== Vertical bubbly flow (theta = 90 deg; extrapolation beyond the paper's\n"
                 "    validated near-horizontal range -- see README) ===\n";
    const double D = 0.05, L = 10.0;
    const int N = 100;

    mfs::SolverOptions opt;
    mfs::FourFieldSolver solver(D, L, N, airWater(), opt);
    solver.setInclinationConstant(90.0 * mfs::constants::pi / 180.0);

    mfs::BoundaryConditions bc;
    bc.inletSuperficialLiquid = 0.5;
    bc.inletSuperficialGas = 0.2;
    bc.inletLiquidHoldup = 0.6;
    bc.inletBubbleFractionOfGas = 0.5; // a mix of continuous gas and bubbles
    bc.outletPressure = 3.0e5;
    bc.seedDisturbance = false;
    solver.setBoundaryConditions(bc);
    solver.initializeStratified(0.6);

    std::ofstream out("vertical_output.csv");
    writeSnapshotHeader(out);

    const double tEnd = 6.0;
    const double snapshotInterval = 1.0;
    double nextSnapshot = 0.0;
    while (solver.time() < tEnd) {
        const double dt = solver.stableTimeStep();
        solver.step(dt);
        if (solver.time() >= nextSnapshot) {
            writeSnapshot(out, solver);
            std::cout << "t=" << solver.time() << " s\n";
            printMassBalance(solver);
            nextSnapshot += snapshotInterval;
        }
    }
    std::cout << "Wrote vertical_output.csv\n\n";
}

void runTerrainCase() {
    std::cout << "=== Terrain-following inclined pipeline (the scenario motivating the\n"
                 "    paper's introduction: flow regimes evolving as inclination changes) ===\n";
    const double D = 0.1, L = 60.0;
    const int N = 400;

    mfs::SolverOptions opt;
    mfs::FourFieldSolver solver(D, L, N, airWater(), opt);
    // A shallow V-section: -5 deg for the first half, +5 deg for the second.
    solver.setInclinationProfile([&](double z) {
        const double deg = (z < L / 2.0) ? -5.0 : 5.0;
        return deg * mfs::constants::pi / 180.0;
    });

    mfs::BoundaryConditions bc;
    bc.inletSuperficialLiquid = 0.6;
    bc.inletSuperficialGas = 4.0;
    bc.inletLiquidHoldup = 0.2;
    bc.outletPressure = 2.0e5;
    bc.seedDisturbance = true;
    bc.disturbanceAmplitude = 0.03;
    bc.disturbanceFrequency = 0.5;
    solver.setBoundaryConditions(bc);
    solver.initializeStratified(0.2);

    std::ofstream out("terrain_output.csv");
    writeSnapshotHeader(out);

    const double tEnd = 10.0;
    const double snapshotInterval = 2.0;
    double nextSnapshot = 0.0;
    while (solver.time() < tEnd) {
        const double dt = solver.stableTimeStep();
        solver.step(dt);
        if (solver.time() >= nextSnapshot) {
            writeSnapshot(out, solver);
            std::cout << "t=" << solver.time() << " s\n";
            printRegimeSummary(solver);
            printMassBalance(solver);
            nextSnapshot += snapshotInterval;
        }
    }
    std::cout << "Wrote terrain_output.csv\n\n";
}

// Mass-conservation benchmark for SolverOptions::enableSurfaceTension,
// methodology directly from Kim, Lee & Kim (2022, Nucl. Eng. Technol.
// 54:4322-4328), Eq. (30): track a mass-error ratio
// e(t) = |m(t) - (m(0)+m_in-m_out)| / (m(0)+m_in-m_out) * 100 against
// time, across mesh resolutions, with the regularization on vs. off, on
// the same Brito (2012) Kelvin-Helmholtz case study this project's other
// regularization checks already use. The paper finds a prior artificial-
// viscosity well-posedness fix (Fullmer et al. 2014) leaks total mass
// because it diffuses each phase's mass equation with its own
// coefficient; this project's enableSurfaceTension already applies its
// biharmonic term to eL and the EXACT SAME quantity, negated, to eG (see
// FourFieldSolver.cpp, updateContinuity()), the four-field analogue of
// that paper's fix (a shared coefficient so the phases' diffusion terms
// cancel when summed, since the phase fractions sum to 1 identically).
// This benchmark verifies that claim empirically rather than by
// inspection: see VALIDATION.md, item 27, for the full account and the
// result already recorded from this exact methodology.
void runMassConservationBenchmark() {
    std::cout << "=== Mass-conservation benchmark for SolverOptions::enableSurfaceTension\n"
                 "    (Kim, Lee & Kim 2022 methodology): mass-error ratio vs. time, across\n"
                 "    mesh resolutions, regularization on vs. off. ===\n";

    const double D = 0.051, thetaDeg = 0.0;
    const double Vsl = 0.055365784383465656, Vsg = 2.6420445461549353;
    const double rhoL = 859.8021132828677, rhoG = 1.1694722577783099;
    const double muL = 0.054, muG = 1.02e-05, sigmaST = 0.0275;
    const double e10 = 0.5;
    const double L = std::clamp(60.0 * D, 0.6, 20.0);
    const double tEnd = 3.0;
    const double snapshotInterval = tEnd / 150.0;

    std::ofstream out("mass_conservation_output.csv");
    out << "N,surfaceTension,time,eLiq_pct,eGas_pct\n";

    for (int N : {60, 150, 300}) {
        for (bool surfaceTension : {false, true}) {
            mfs::FluidProperties fluid;
            fluid.rhoLiquid = rhoL; fluid.muLiquid = muL; fluid.sigma = sigmaST; fluid.muGas = muG;
            fluid.temperature = 293.0;
            fluid.gasConstant = 101325.0 / (rhoG * fluid.temperature);
            mfs::SolverOptions opt;
            opt.enableSurfaceTension = surfaceTension;
            mfs::FourFieldSolver solver(D, L, N, fluid, opt);
            solver.setInclinationConstant(thetaDeg * mfs::constants::pi / 180.0);

            mfs::BoundaryConditions bc;
            bc.inletSuperficialLiquid = Vsl;
            bc.inletSuperficialGas = Vsg;
            bc.inletLiquidHoldup = e10;
            bc.outletPressure = 101325.0;
            bc.seedDisturbance = true;
            bc.disturbanceAmplitude = 0.05;
            bc.disturbanceFrequency = 0.8;
            solver.setBoundaryConditions(bc);
            solver.initializeStratified(e10);

            const double mLiq0 = solver.totalLiquidMass();
            const double mGas0 = solver.totalGasMass();
            double liqInIntegral = 0.0, liqOutIntegral = 0.0, gasInIntegral = 0.0, gasOutIntegral = 0.0;
            double eLiqMax = 0.0, eGasMax = 0.0;
            double nextSnapshot = 0.0;

            while (solver.time() < tEnd) {
                const double dt = solver.stableTimeStep();
                const auto mfBefore = solver.boundaryMassFluxes();
                solver.step(dt);
                liqInIntegral += mfBefore.liquidIn * dt;
                liqOutIntegral += mfBefore.liquidOut * dt;
                gasInIntegral += mfBefore.gasIn * dt;
                gasOutIntegral += mfBefore.gasOut * dt;

                const double mLiqExpected = mLiq0 + liqInIntegral - liqOutIntegral;
                const double mGasExpected = mGas0 + gasInIntegral - gasOutIntegral;
                const double eLiq = std::fabs(solver.totalLiquidMass() - mLiqExpected) /
                                     std::max(std::fabs(mLiqExpected), 1e-30) * 100.0;
                const double eGas = std::fabs(solver.totalGasMass() - mGasExpected) /
                                     std::max(std::fabs(mGasExpected), 1e-30) * 100.0;
                eLiqMax = std::max(eLiqMax, eLiq);
                eGasMax = std::max(eGasMax, eGas);

                if (solver.time() >= nextSnapshot) {
                    out << N << ',' << (surfaceTension ? 1 : 0) << ',' << solver.time() << ','
                        << eLiq << ',' << eGas << '\n';
                    nextSnapshot += snapshotInterval;
                }
            }
            std::cout << "N=" << N << "  surfaceTension=" << surfaceTension
                      << "  eLiq_max=" << eLiqMax << "%  eGas_max=" << eGasMax << "%\n";
        }
    }
    std::cout << "Wrote mass_conservation_output.csv\n\n";
}

// Self-consistent gas-wall / liquid-wall friction-factor prediction, for
// comparison against Wang, Liu, Chang, Zhao & Guo (2019, Int. J. Heat Mass
// Transfer 143:118537) Figs. 12(a)/14(a) (gas-wall) and 17 (liquid-wall).
// That paper reports MEASURED friction factors from their own air-water
// wavy-stratified-flow experiments (50 mm ID pipe, 11.6 m test section) at
// atmosphere and 2 MPa, vs. the gas/liquid-phase Reynolds number -- not a
// superficial-velocity sweep, so this runs this project's own solver at
// matched (D, L, fluid properties) over a grid of (Vsl, Vsg) spanning their
// reported wave-stratified flow envelope (their Fig. 4), lets it reach a
// quasi-steady holdup the same way validate_quantitative.cpp does, and
// THEN derives the actual phase Reynolds numbers from the solver's own
// predicted holdup (not an externally imposed one) before evaluating this
// project's existing wallFrictionFactor() (SolverOptions::wallCorrelation,
// Taitel & Dukler 1976) at those self-consistent conditions -- i.e. this is
// genuinely what the solver itself predicts for f_G/f_L, not a standalone
// re-evaluation of the closure formula at arbitrary Re values.
void runWang2019FrictionComparison() {
    std::cout << "=== Wang et al. (2019, IJHMT 143:118537) gas-wall / liquid-wall\n"
                 "    friction-factor comparison: this solver's own self-consistent\n"
                 "    predicted Re_G/f_G, Re_L/f_L across a (Vsl,Vsg) sweep spanning\n"
                 "    their reported wave-stratified envelope, atmosphere and 2 MPa. ===\n";

    struct Condition {
        std::string label;
        double rhoL, muL, rhoG, muG, outletPressure;
        mfs::FrictionCorrelation gasWallLit; // the paper's own fitted correlation, this pressure
    };
    const std::vector<Condition> conditions = {
        {"atm",  998.0,   1.0e-3,  1.2,  1.8e-5,  101325.0, mfs::FrictionCorrelation::WangGuo2019GasWallAtm},
        {"2MPa", 1005.4,  7.24e-4, 23.4, 1.87e-5, 2.0e6,    mfs::FrictionCorrelation::WangGuo2019GasWall2MPa},
    };

    const double D = 0.050;
    const double L = 11.6;
    const int N = 60;

    std::ofstream out("wang2019_friction_output.csv");
    out << "pressure,Vsl,Vsg,eL,ReG,fG,fG_wangguo2019,ReL,fL,fL_wangguo2019,blew_up\n";

    for (const auto& cond : conditions) {
        const std::vector<double> VlsList = (cond.label == "atm")
            ? std::vector<double>{0.023, 0.034, 0.049, 0.07, 0.1}
            : std::vector<double>{0.01, 0.023, 0.034, 0.049, 0.07, 0.1};
        const std::vector<double> VgsList = (cond.label == "atm")
            ? std::vector<double>{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}
            : std::vector<double>{0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5};

        for (double Vsl : VlsList) {
            for (double Vsg : VgsList) {
                mfs::FluidProperties fluid;
                fluid.rhoLiquid = cond.rhoL;
                fluid.muLiquid = cond.muL;
                fluid.muGas = cond.muG;
                fluid.sigma = 0.072;
                fluid.temperature = 293.0;
                fluid.gasConstant = cond.outletPressure / (cond.rhoG * fluid.temperature);

                mfs::SolverOptions opt;
                mfs::FourFieldSolver solver(D, L, N, fluid, opt);
                solver.setInclinationConstant(0.0);

                mfs::BoundaryConditions bc;
                bc.inletSuperficialLiquid = Vsl;
                bc.inletSuperficialGas = Vsg;
                const double homogeneousGuess = Vsl / std::max(Vsl + Vsg, 1e-6);
                bc.inletLiquidHoldup = std::clamp(homogeneousGuess, 0.02, 0.98);
                bc.outletPressure = cond.outletPressure;
                bc.seedDisturbance = true;
                bc.disturbanceAmplitude = 0.03;
                bc.disturbanceFrequency = 0.6;
                solver.setBoundaryConditions(bc);
                solver.initializeStratified(bc.inletLiquidHoldup);

                const double umix = Vsl + Vsg;
                const double residenceTime = L / std::max(umix, 0.02);
                const double tEnd = std::clamp(10.0 * residenceTime, 1.0, 8.0);
                const double warmupEnd = 0.5 * tEnd;
                const int maxSteps = 300000;
                const int iStart = (2 * N) / 3;

                std::vector<double> holdupSamples;
                double nextSample = warmupEnd;
                const double sampleInterval = (tEnd - warmupEnd) / 30.0;
                bool blewUp = false;
                int step = 0;
                while (solver.time() < tEnd && step < maxSteps) {
                    const double dt = solver.stableTimeStep();
                    solver.step(dt);
                    ++step;
                    for (double v : solver.state().el) {
                        if (std::isnan(v) || std::isinf(v)) { blewUp = true; break; }
                    }
                    if (blewUp) break;
                    if (solver.time() >= nextSample) {
                        const auto& s = solver.state();
                        double hSum = 0.0;
                        int hN = 0;
                        for (int i = iStart; i < N; ++i) { hSum += s.el[i] + s.ed[i]; ++hN; }
                        if (hN > 0) holdupSamples.push_back(hSum / hN);
                        nextSample += sampleInterval;
                    }
                }

                if (blewUp || holdupSamples.empty()) {
                    out << cond.label << ',' << Vsl << ',' << Vsg << ",,,,,,,,1\n";
                    continue;
                }
                double hMean = 0.0;
                for (double v : holdupSamples) hMean += v;
                const double eL = hMean / holdupSamples.size();

                const mfs::PipeGeometry geom(D);
                const auto sg = geom.fromAreaFraction(eL);
                const double uG = Vsg / std::max(1.0 - eL, 1e-6);
                const double uL = Vsl / std::max(eL, 1e-6);
                const double ReG = mfs::reynoldsNumber(cond.rhoG, uG, sg.D2, cond.muG);
                const double ReL = mfs::reynoldsNumber(cond.rhoL, uL, sg.D1, cond.muL);
                const double fG = mfs::wallFrictionFactor({cond.rhoG, uG, sg.D2, cond.muG},
                                                            mfs::FrictionCorrelation::TaitelDukler1976);
                const double fL = mfs::wallFrictionFactor({cond.rhoL, uL, sg.D1, cond.muL},
                                                            mfs::FrictionCorrelation::TaitelDukler1976);
                // The paper's own fitted correlation (Eqs. 17/18/20),
                // evaluated at this same solver-predicted, self-consistent
                // (ReG, ReL) -- see Closures.hpp's WangGuo2019* comments
                // for why these are literature-specific, not general
                // physics. No atmosphere liquid-wall fit was published, so
                // fL_wangguo2019 is left blank for that condition.
                const double fG_lit = mfs::wallFrictionFactor({cond.rhoG, uG, sg.D2, cond.muG}, cond.gasWallLit);
                std::ostringstream fLlit;
                if (cond.label == "2MPa") {
                    const double v = mfs::wallFrictionFactor(
                        {cond.rhoL, uL, sg.D1, cond.muL}, mfs::FrictionCorrelation::WangGuo2019LiquidWall2MPa);
                    fLlit << v;
                }

                out << cond.label << ',' << Vsl << ',' << Vsg << ',' << eL << ','
                    << ReG << ',' << fG << ',' << fG_lit << ','
                    << ReL << ',' << fL << ',' << fLlit.str() << ",0\n";
            }
        }
    }
    std::cout << "Wrote wang2019_friction_output.csv\n\n";
}

} // namespace

int main(int argc, char** argv) {
    const std::string caseName = (argc > 1) ? argv[1] : "horizontal";

    if (caseName == "horizontal") {
        runHorizontalCase(/*useAmr=*/false, /*useMovingMesh=*/false, "horizontal_output.csv");
    } else if (caseName == "vertical") {
        runVerticalCase();
    } else if (caseName == "terrain") {
        runTerrainCase();
    } else if (caseName == "horizontal_amr") {
        runHorizontalAmrComparison();
    } else if (caseName == "horizontal_movingmesh") {
        runHorizontalMovingMeshComparison();
    } else if (caseName == "horizontal_limiter") {
        runHorizontalLimiterComparison();
    } else if (caseName == "mass_conservation") {
        runMassConservationBenchmark();
    } else if (caseName == "wang2019_friction") {
        runWang2019FrictionComparison();
    } else if (caseName == "all") {
        runHorizontalCase(/*useAmr=*/false, /*useMovingMesh=*/false, "horizontal_output.csv");
        runVerticalCase();
        runTerrainCase();
        runHorizontalAmrComparison();
        runHorizontalMovingMeshComparison();
        runHorizontalLimiterComparison();
        runMassConservationBenchmark();
        runWang2019FrictionComparison();
    } else {
        std::cerr << "Unknown case '" << caseName
                   << "'. Use one of: horizontal, vertical, terrain, horizontal_amr, "
                      "horizontal_movingmesh, horizontal_limiter, mass_conservation, "
                      "wang2019_friction, all\n";
        return 1;
    }
    return 0;
}
