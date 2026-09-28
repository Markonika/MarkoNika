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
//                         The same horizontal case, once at a fixed
//                         resolution and once with adaptive mesh
//                         refinement (SolverOptions::amr) enabled, timed
//                         side by side. AMR refines using the same
//                         Kelvin-Helmholtz F parameter already computed for
//                         the interfacial friction closure (Closures.hpp)
//                         as its indicator, coarsening cells where the
//                         interface is quiescent and never refining finer
//                         than the baseline resolution by default -- see
//                         FourFieldSolver.hpp, SolverOptions::amr.
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
                                      bool quiet = false) {
    if (!quiet) {
        std::cout << "=== Horizontal stratified -> slug transient (cf. Section 4/5.2)"
                   << (useAmr ? ", with AMR ===\n" : useMovingMesh ? ", with moving mesh ===\n" : " ===\n");
    }
    const double D = 0.08, L = 30.0;
    const int N = 300; // ~1 diameter per cell, matching the paper's grid-independence finding (Fig. 2)

    mfs::SolverOptions opt;
    opt.amr.enabled = useAmr;
    opt.movingMesh.enabled = useMovingMesh;
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
    std::cout << "=== AMR vs. fixed-grid comparison on the horizontal slug-formation case ===\n";
    std::cout << "Running fixed grid (N=300)...\n";
    auto t0 = std::chrono::steady_clock::now();
    const auto fixedStats = runHorizontalCase(false, false, "horizontal_fixed_output.csv", /*quiet=*/true);
    auto t1 = std::chrono::steady_clock::now();

    std::cout << "Running with AMR enabled...\n";
    const auto amrStats = runHorizontalCase(true, false, "horizontal_amr_output.csv", /*quiet=*/true);
    auto t2 = std::chrono::steady_clock::now();

    const double fixedSec = std::chrono::duration<double>(t1 - t0).count();
    const double amrSec = std::chrono::duration<double>(t2 - t1).count();
    std::cout << "\nFixed grid : " << fixedStats.steps << " steps, N=" << fixedStats.finalN
              << " cells throughout, " << fixedSec << " s\n";
    std::cout << "AMR        : " << amrStats.steps << " steps, N=" << amrStats.finalN
              << " cells at end, " << amrSec << " s (" << (fixedSec / amrSec) << "x)\n";
    std::cout << "Wrote horizontal_fixed_output.csv and horizontal_amr_output.csv\n"
                 "(the AMR file's per-cell rows are at whatever local resolution that\n"
                 "snapshot's mesh had -- z spacing varies row to row within a snapshot).\n\n";
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
    } else if (caseName == "all") {
        runHorizontalCase(/*useAmr=*/false, /*useMovingMesh=*/false, "horizontal_output.csv");
        runVerticalCase();
        runTerrainCase();
        runHorizontalAmrComparison();
        runHorizontalMovingMeshComparison();
    } else {
        std::cerr << "Unknown case '" << caseName
                   << "'. Use one of: horizontal, vertical, terrain, horizontal_amr, "
                      "horizontal_movingmesh, all\n";
        return 1;
    }
    return 0;
}
