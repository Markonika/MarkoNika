// Validation driver: compares the four-field solver's predicted flow
// regime (via mfs::classifyFromHistory, Section 5.1 criteria) against the
// experimentally observed flow pattern reported in Shoham's (1982) PhD
// thesis dataset ("Flow pattern transition and characterization in
// gas-liquid two-phase flow in inclined pipes", Tel Aviv University),
// spanning pipe inclination from -90 deg (vertical downward) to +90 deg
// (vertical upward). The dataset (air-water, D = 25.4 mm and 51 mm) was
// obtained from the public compilation:
//   BioAITeam, "Machine learning applications to predict two-phase flow
//   patterns", https://github.com/BioAITeam/Machine-learning-applications-to-predict-two-phase-flow-patterns
//   (Databases/ShohamDB.csv), itself citing Shoham (1982).
//
// This is NOT a reproduction of the original paper's own validation
// (which used specific published cases from Nydal et al. 1992, Andritsos
// et al. 1989, and Taitel & Dukler 1976 -- see README.md and
// VALIDATION.md). It is an independent, additional check of the same
// flow-regime-independent model against a much larger, angle-resolved
// experimental set, run at coarse resolution and reduced pipe length for
// tractability (see VALIDATION.md for the methodology and its limits).
//
// Usage: validate_shoham <input.csv> <output.csv>

#include "mfs/Constants.hpp"
#include "mfs/FlowRegimeClassifier.hpp"
#include "mfs/FourFieldSolver.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Case {
    double Vsl, Vsg, VisL, VisG, DenL, DenG, ST, Ang, ID;
    std::string FlowPattern;
};

std::vector<Case> readCases(const std::string& path) {
    std::ifstream in(path);
    std::string line;
    std::getline(in, line); // header
    std::vector<Case> cases;
    while (std::getline(in, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string tok;
        std::vector<std::string> f;
        while (std::getline(ss, tok, ',')) f.push_back(tok);
        if (f.size() < 10) continue;
        Case c;
        c.Vsl = std::stod(f[0]);
        c.Vsg = std::stod(f[1]);
        c.VisL = std::stod(f[2]);
        c.VisG = std::stod(f[3]);
        c.DenL = std::stod(f[4]);
        c.DenG = std::stod(f[5]);
        c.ST = std::stod(f[6]);
        c.Ang = std::stod(f[7]);
        c.ID = std::stod(f[8]);
        c.FlowPattern = f[9];
        cases.push_back(c);
    }
    return cases;
}

// Map both the experimental label and the model's classifier output onto a
// common 4-category scheme so they can be compared directly.
// Shoham labels: I=Intermittent(slug/plug), A=Annular, SW=Stratified Wavy,
// SS=Stratified Smooth, DB=Dispersed Bubble, B=Bubble.
std::string canonicalLabel(const std::string& shohamLabel) {
    if (shohamLabel == "I") return "slug";
    if (shohamLabel == "A") return "annular";
    if (shohamLabel == "SW" || shohamLabel == "SS") return "stratified";
    if (shohamLabel == "DB" || shohamLabel == "B") return "bubbly";
    return "unknown";
}

std::string canonicalLabel(mfs::FlowRegime r) {
    switch (r) {
        case mfs::FlowRegime::Stratified: return "stratified";
        case mfs::FlowRegime::Annular: return "annular";
        case mfs::FlowRegime::Slug: return "slug";
        case mfs::FlowRegime::Bubbly: return "bubbly";
    }
    return "unknown";
}

// Tunable via extra CLI args so the same driver can run the fast, large-N
// sweep and a slower, long-development-length spot check without a
// separate copy of this file. Defaults match the original fast sweep.
double g_lengthInDiameters = 60.0;
int g_cellsPerRun = 40;
double g_maxResidenceTimes = 8.0;
double g_maxSimTime = 6.0;

std::string runOneCase(const Case& c) {
    const double D = c.ID;
    const double L = std::clamp(g_lengthInDiameters * D, 0.6, 20.0);
    const int N = g_cellsPerRun;

    mfs::FluidProperties fluid;
    fluid.rhoLiquid = c.DenL;
    fluid.muLiquid = c.VisL;
    fluid.sigma = c.ST;
    fluid.muGas = c.VisG;
    fluid.temperature = 293.0;
    const double outletPressure = 101325.0;
    fluid.gasConstant = outletPressure / (c.DenG * fluid.temperature); // reproduces DenG at outletPressure

    mfs::SolverOptions opt;
    mfs::FourFieldSolver solver(D, L, N, fluid, opt);
    solver.setInclinationConstant(c.Ang * mfs::constants::pi / 180.0);

    mfs::BoundaryConditions bc;
    bc.inletSuperficialLiquid = c.Vsl;
    bc.inletSuperficialGas = c.Vsg;
    const double homogeneousGuess = c.Vsl / std::max(c.Vsl + c.Vsg, 1e-6);
    bc.inletLiquidHoldup = std::clamp(homogeneousGuess, 0.05, 0.95);
    bc.outletPressure = outletPressure;
    bc.seedDisturbance = true;
    bc.disturbanceAmplitude = 0.03;
    bc.disturbanceFrequency = 0.6;
    solver.setBoundaryConditions(bc);
    solver.initializeStratified(bc.inletLiquidHoldup);

    // Run for several residence times, hard-capped on both simulated time
    // and step count for tractability across a large sweep. Snapshots for
    // the regime classifier are taken at fixed simulated-time intervals
    // (not fixed step counts, since dt varies enormously across this
    // sweep) over the back half of the run, letting transients settle.
    const double umix = c.Vsl + c.Vsg;
    const double residenceTime = L / std::max(umix, 0.05);
    const double tEnd = std::clamp(g_maxResidenceTimes * residenceTime, 1.0, g_maxSimTime);
    const double warmupEnd = 0.5 * tEnd;
    const int maxSteps = 200000;
    const double sampleInterval = (tEnd - warmupEnd) / 40.0;

    std::vector<std::vector<double>> egHistory;
    double nextSample = warmupEnd;
    int step = 0;
    bool blewUp = false;
    while (solver.time() < tEnd && step < maxSteps) {
        const double dt = solver.stableTimeStep();
        solver.step(dt);
        ++step;
        if (solver.time() >= nextSample) {
            egHistory.push_back(solver.state().eg);
            nextSample += sampleInterval;
        }
        for (double v : solver.state().el) {
            if (std::isnan(v) || std::isinf(v)) { blewUp = true; break; }
        }
        if (blewUp) break;
    }
    if (blewUp) return "unknown";

    if (egHistory.size() < 3) return "unknown"; // blew up or too short to classify
    const auto regimes = mfs::classifyFromHistory(egHistory);

    // Classify using only the back third of the pipe (closest to fully
    // developed, per the paper's own emphasis on development length).
    const int n = static_cast<int>(regimes.size());
    const int start = (2 * n) / 3;
    std::map<mfs::FlowRegime, int> counts;
    for (int i = start; i < n; ++i) counts[regimes[i]]++;
    mfs::FlowRegime majority = mfs::FlowRegime::Stratified;
    int best = -1;
    for (auto& [r, cnt] : counts) if (cnt > best) { best = cnt; majority = r; }
    return canonicalLabel(majority);
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: validate_shoham <input.csv> <output.csv> "
                     "[lengthInDiameters] [cellsPerRun] [maxResidenceTimes] [maxSimTime]\n";
        return 1;
    }
    if (argc > 3) g_lengthInDiameters = std::stod(argv[3]);
    if (argc > 4) g_cellsPerRun = std::stoi(argv[4]);
    if (argc > 5) g_maxResidenceTimes = std::stod(argv[5]);
    if (argc > 6) g_maxSimTime = std::stod(argv[6]);

    const auto cases = readCases(argv[1]);
    std::ofstream out(argv[2]);
    out << "Ang,Vsl,Vsg,ID,DenL,DenG,Labeled,LabeledCanonical,Predicted,Agree\n";

    std::map<std::pair<std::string, std::string>, int> confusion;
    int total = 0, agree = 0;

    for (std::size_t i = 0; i < cases.size(); ++i) {
        const auto& c = cases[i];
        const std::string predicted = runOneCase(c);
        const std::string labeled = canonicalLabel(c.FlowPattern);
        const bool ok = (predicted == labeled);
        confusion[{labeled, predicted}]++;
        total++;
        if (ok) agree++;

        out << c.Ang << ',' << c.Vsl << ',' << c.Vsg << ',' << c.ID << ',' << c.DenL << ',' << c.DenG
            << ',' << c.FlowPattern << ',' << labeled << ',' << predicted << ',' << (ok ? 1 : 0) << '\n';
        out.flush();

        std::cout << "[" << (i + 1) << "/" << cases.size() << "] Ang=" << c.Ang
                  << " Vsl=" << c.Vsl << " Vsg=" << c.Vsg << " D=" << c.ID
                  << " labeled=" << labeled << " predicted=" << predicted
                  << (ok ? "  MATCH" : "  ---") << std::endl;
    }

    std::cout << "\n=== Overall accuracy: " << agree << "/" << total << " = "
              << (100.0 * agree / std::max(total, 1)) << "% ===\n";
    std::cout << "\nConfusion matrix (rows=labeled, cols=predicted):\n";
    for (const auto& [key, cnt] : confusion) {
        std::cout << "  labeled=" << key.first << " predicted=" << key.second << " : " << cnt << "\n";
    }
    return 0;
}
