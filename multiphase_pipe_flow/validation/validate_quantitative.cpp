// Quantitative validation driver: runs the four-field solver at matched
// operating conditions (Vsl, Vsg, D, inclination, fluid properties) drawn
// from real experimental campaigns and reports the solver's predicted
// mean liquid holdup and mean frictional pressure gradient, for later
// comparison against the experimental values already attached to each
// case in the input CSV.
//
// Input CSV columns (see scratchpad/quant_validation/unified_cases.csv):
//   case_id,source,D_m,angle_deg,Vsl,Vsg,rhoL,rhoG,muL,muG,sigma,
//   exp_holdup,exp_dPdz_Pa_m,exp_pattern
//
// Output CSV adds: pred_holdup,pred_dPdz_Pa_m,blew_up
//
// Usage: validate_quantitative <input.csv> <output.csv> [closureCorrectionScale] [closureF0Shift] [closureDepositionVelocityScale] [enableTurbulentViscosity(0|1)] [turbulentMixingLengthFraction] [enableInterfacialPressureJump(0|1)] [interfacialPressureJumpCoefficient] [enableInterfacialSlopeClosure(0|1)] [interfacialSlopeCoefficient] [N] [enableRoeScheme(0|1)]
// The optional third/fourth/fifth arguments, if any is given, set
// SolverOptions::enableDataDrivenClosureCorrection=true and
// closureCorrectionScale/closureF0Shift/closureDepositionVelocityScale
// to the given values (default 1.0/0.0/1.0, i.e. no-ops, for whichever
// are omitted) for every case -- used to line-search the fits in
// VALIDATION.md's data-driven closure correction section; omitted
// entirely, behavior is the original, unmodified closure. The optional
// sixth/seventh arguments set SolverOptions::enableTurbulentViscosity
// and turbulentMixingLengthFraction, used to test the turbulent-
// viscosity regularization against the same holdup validation. The
// optional eighth/ninth arguments, if the eighth is nonzero, set
// SolverOptions::enableInterfacialPressureJump=true and
// interfacialPressureJumpCoefficient (default 1.5 if omitted) for every
// case -- used to test the Stewart-Wendroff/Toumi-Kumbaro-style
// interfacial pressure-jump regularization against the same holdup
// validation. The optional tenth/eleventh arguments, if the tenth is
// nonzero, set SolverOptions::enableInterfacialSlopeClosure=true and
// interfacialSlopeCoefficient (default 1.0 if omitted) -- used to test
// the Buist & Sanderse-motivated interfacial-slope closure input. The
// optional twelfth argument overrides the mesh resolution N (default
// 60, the value every prior invocation of this driver implicitly used);
// added specifically to run the SAME 223/126-case train/test split at
// multiple resolutions for the multi-resolution curriculum fit of
// interfacialSlopeCoefficient (see VALIDATION.md), adapting the coarse-
// to-fine curriculum-learning idea of Garnier, Viquerat & Hachem,
// "Curriculum Learning for Mesh-based Simulations" (arXiv:2509.13138),
// to this driver's scalar line-search methodology rather than a trained
// network (see VALIDATION.md for the full account and the honest
// distinction between the two).

#include "mfs/Constants.hpp"
#include "mfs/FourFieldSolver.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct Case {
    std::string case_id, source, exp_pattern;
    double D, angle_deg, Vsl, Vsg, rhoL, rhoG, muL, muG, sigma;
    double exp_holdup, exp_dPdz; // exp_dPdz may be NaN (not reported)
};

std::vector<std::string> splitCsvLine(const std::string& line) {
    // Split on ',' preserving trailing empty fields (std::getline on a
    // stringstream silently drops a final empty token, which matters here
    // since exp_pattern -- the last column -- is empty for every non-Kokal
    // row).
    std::vector<std::string> out;
    std::size_t start = 0;
    while (true) {
        const std::size_t comma = line.find(',', start);
        if (comma == std::string::npos) {
            out.push_back(line.substr(start));
            break;
        }
        out.push_back(line.substr(start, comma - start));
        start = comma + 1;
    }
    return out;
}

double parseOrNan(const std::string& s) {
    if (s.empty()) return std::nan("");
    try { return std::stod(s); } catch (...) { return std::nan(""); }
}

std::vector<Case> readCases(const std::string& path) {
    std::ifstream in(path);
    std::string line;
    std::getline(in, line); // header
    std::vector<Case> cases;
    while (std::getline(in, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
        if (line.empty()) continue;
        auto f = splitCsvLine(line);
        if (f.size() < 14) continue;
        Case c;
        c.case_id = f[0]; c.source = f[1];
        c.D = std::stod(f[2]); c.angle_deg = std::stod(f[3]);
        c.Vsl = std::stod(f[4]); c.Vsg = std::stod(f[5]);
        c.rhoL = std::stod(f[6]); c.rhoG = std::stod(f[7]);
        c.muL = std::stod(f[8]); c.muG = std::stod(f[9]);
        c.sigma = std::stod(f[10]);
        c.exp_holdup = parseOrNan(f[11]);
        c.exp_dPdz = parseOrNan(f[12]);
        c.exp_pattern = f[13];
        cases.push_back(c);
    }
    return cases;
}

struct RunResult {
    double predHoldup = std::nan("");
    double predDPdz = std::nan("");
    bool blewUp = false;
};

RunResult runOneCase(const Case& c, bool haveClosureScale, double closureScale,
                      bool haveF0Shift, double f0Shift,
                      bool haveDepositionScale, double depositionScale,
                      bool enableTurbVisc, bool haveMixingLength, double mixingLengthFraction,
                      bool enableInterfacialPressureJump, bool haveCi, double Ci,
                      bool enableInterfacialSlopeClosure, bool haveSlopeCoeff, double slopeCoeff,
                      int N, bool enableRoeScheme) {
    const double D = c.D;
    const double L = std::clamp(60.0 * D, 0.6, 20.0);

    mfs::FluidProperties fluid;
    fluid.rhoLiquid = c.rhoL;
    fluid.muLiquid = c.muL;
    fluid.sigma = c.sigma;
    fluid.muGas = c.muG;
    fluid.temperature = 293.0;
    const double outletPressure = 101325.0;
    fluid.gasConstant = outletPressure / (c.rhoG * fluid.temperature);

    mfs::SolverOptions opt;
    if (haveClosureScale || haveF0Shift || haveDepositionScale) {
        opt.enableDataDrivenClosureCorrection = true;
        if (haveClosureScale) opt.closureCorrectionScale = closureScale;
        if (haveF0Shift) opt.closureF0Shift = f0Shift;
        if (haveDepositionScale) opt.closureDepositionVelocityScale = depositionScale;
    }
    if (enableTurbVisc) {
        opt.enableTurbulentViscosity = true;
        if (haveMixingLength) opt.turbulentMixingLengthFraction = mixingLengthFraction;
    }
    if (enableInterfacialPressureJump) {
        opt.enableInterfacialPressureJump = true;
        if (haveCi) opt.interfacialPressureJumpCoefficient = Ci;
    }
    if (enableInterfacialSlopeClosure) {
        opt.enableInterfacialSlopeClosure = true;
        if (haveSlopeCoeff) opt.interfacialSlopeCoefficient = slopeCoeff;
    }
    opt.enableRoeScheme = enableRoeScheme;
    mfs::FourFieldSolver solver(D, L, N, fluid, opt);
    solver.setInclinationConstant(c.angle_deg * mfs::constants::pi / 180.0);

    mfs::BoundaryConditions bc;
    bc.inletSuperficialLiquid = c.Vsl;
    bc.inletSuperficialGas = c.Vsg;
    const double homogeneousGuess = c.Vsl / std::max(c.Vsl + c.Vsg, 1e-6);
    bc.inletLiquidHoldup = std::clamp(homogeneousGuess, 0.02, 0.98);
    bc.outletPressure = outletPressure;
    bc.seedDisturbance = true;
    bc.disturbanceAmplitude = 0.03;
    bc.disturbanceFrequency = 0.6;
    solver.setBoundaryConditions(bc);
    solver.initializeStratified(bc.inletLiquidHoldup);

    const double umix = c.Vsl + c.Vsg;
    const double residenceTime = L / std::max(umix, 0.02);
    const double tEnd = std::clamp(10.0 * residenceTime, 1.0, 8.0);
    const double warmupEnd = 0.5 * tEnd;
    const int maxSteps = 300000;
    const double sampleInterval = (tEnd - warmupEnd) / 30.0;

    // back third of the domain, away from inlet/outlet boundary effects
    const int iStart = (2 * N) / 3;

    std::vector<double> holdupSamples;
    std::vector<double> dpdzSamples;
    double nextSample = warmupEnd;
    int step = 0;
    RunResult result;

    while (solver.time() < tEnd && step < maxSteps) {
        const double dt = solver.stableTimeStep();
        solver.step(dt);
        ++step;
        for (double v : solver.state().el) {
            if (std::isnan(v) || std::isinf(v)) { result.blewUp = true; break; }
        }
        if (result.blewUp) break;

        if (solver.time() >= nextSample) {
            const auto& s = solver.state();
            double hSum = 0.0;
            int hN = 0;
            for (int i = iStart; i < N; ++i) {
                hSum += s.el[i] + s.ed[i];
                ++hN;
            }
            if (hN > 0) holdupSamples.push_back(hSum / hN);

            // linear regression of P vs z over the back third
            double sz = 0, sp = 0, szz = 0, szp = 0;
            int nreg = 0;
            for (int i = iStart; i < N; ++i) {
                const double z = s.cellCenter(i);
                const double p = s.P[i];
                sz += z; sp += p; szz += z * z; szp += z * p;
                ++nreg;
            }
            if (nreg >= 2) {
                const double denom = nreg * szz - sz * sz;
                if (std::fabs(denom) > 1e-12) {
                    const double slope = (nreg * szp - sz * sp) / denom;
                    dpdzSamples.push_back(-slope); // report as a positive pressure DROP gradient
                }
            }
            nextSample += sampleInterval;
        }
    }

    if (result.blewUp || holdupSamples.empty()) {
        result.blewUp = true;
        return result;
    }
    double hMean = 0.0;
    for (double v : holdupSamples) hMean += v;
    result.predHoldup = hMean / holdupSamples.size();

    if (!dpdzSamples.empty()) {
        double dMean = 0.0;
        for (double v : dpdzSamples) dMean += v;
        result.predDPdz = dMean / dpdzSamples.size();
    }
    return result;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: validate_quantitative <input.csv> <output.csv> [closureCorrectionScale] [closureF0Shift] [closureDepositionVelocityScale] [enableTurbulentViscosity(0|1)] [turbulentMixingLengthFraction] [enableInterfacialPressureJump(0|1)] [interfacialPressureJumpCoefficient] [enableInterfacialSlopeClosure(0|1)] [interfacialSlopeCoefficient] [N] [enableRoeScheme(0|1)]\n";
        return 1;
    }
    const bool haveClosureScale = (argc >= 4);
    const double closureScale = haveClosureScale ? std::stod(argv[3]) : 1.0;
    const bool haveF0Shift = (argc >= 5);
    const double f0Shift = haveF0Shift ? std::stod(argv[4]) : 0.0;
    const bool haveDepositionScale = (argc >= 6);
    const double depositionScale = haveDepositionScale ? std::stod(argv[5]) : 1.0;
    const bool enableTurbVisc = (argc >= 7) && std::stoi(argv[6]) != 0;
    const bool haveMixingLength = (argc >= 8);
    const double mixingLengthFraction = haveMixingLength ? std::stod(argv[7]) : 0.1;
    const bool enableInterfacialPressureJump = (argc >= 9) && std::stoi(argv[8]) != 0;
    const bool haveCi = (argc >= 10);
    const double Ci = haveCi ? std::stod(argv[9]) : 1.5;
    const bool enableInterfacialSlopeClosure = (argc >= 11) && std::stoi(argv[10]) != 0;
    const bool haveSlopeCoeff = (argc >= 12);
    const double slopeCoeff = haveSlopeCoeff ? std::stod(argv[11]) : 1.0;
    // Mesh resolution N, default 60 (exactly the prior hardcoded value --
    // bit-exact with every earlier invocation that omits this argument).
    // Added to let a single driver be reused across mesh resolutions
    // instead of a one-off scratch program per resolution -- see
    // VALIDATION.md's multi-resolution curriculum fit for why this
    // mattered in practice.
    const int N = (argc >= 13) ? std::stoi(argv[12]) : 60;
    // See SolverOptions::enableRoeScheme; no extra coefficient needed
    // (a scheme choice, not a tunable).
    const bool enableRoeScheme = (argc >= 14) && std::stoi(argv[13]) != 0;
    const auto cases = readCases(argv[1]);

    // Each case is an independent FourFieldSolver instance (no shared
    // mutable state), so the outer loop is embarrassingly parallel. Results
    // are collected into a preallocated vector indexed by case so that the
    // CSV is still written out in the original, deterministic case order
    // regardless of which thread finished which case first or when.
    std::vector<RunResult> results(cases.size());
#ifdef _OPENMP
    std::cout << "[OpenMP] running " << cases.size() << " cases across "
              << omp_get_max_threads() << " threads\n";
#pragma omp parallel for schedule(dynamic)
#endif
    for (std::size_t i = 0; i < cases.size(); ++i) {
        results[i] = runOneCase(cases[i], haveClosureScale, closureScale, haveF0Shift, f0Shift,
                                 haveDepositionScale, depositionScale, enableTurbVisc,
                                 haveMixingLength, mixingLengthFraction,
                                 enableInterfacialPressureJump, haveCi, Ci,
                                 enableInterfacialSlopeClosure, haveSlopeCoeff, slopeCoeff, N,
                                 enableRoeScheme);
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            std::cout << "[" << (i + 1) << "/" << cases.size() << "] " << cases[i].case_id
                      << "  predHoldup=" << results[i].predHoldup
                      << "  expHoldup=" << cases[i].exp_holdup
                      << (results[i].blewUp ? "  BLEW UP" : "") << std::endl;
        }
    }

    std::ofstream out(argv[2]);
    out << "case_id,source,D_m,angle_deg,Vsl,Vsg,exp_holdup,exp_dPdz_Pa_m,exp_pattern,"
           "pred_holdup,pred_dPdz_Pa_m,blew_up\n";
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const auto& c = cases[i];
        const auto& r = results[i];
        out << c.case_id << ',' << c.source << ',' << c.D << ',' << c.angle_deg << ','
            << c.Vsl << ',' << c.Vsg << ','
            << (std::isnan(c.exp_holdup) ? "" : std::to_string(c.exp_holdup)) << ','
            << (std::isnan(c.exp_dPdz) ? "" : std::to_string(c.exp_dPdz)) << ','
            << c.exp_pattern << ','
            << (std::isnan(r.predHoldup) ? "" : std::to_string(r.predHoldup)) << ','
            << (std::isnan(r.predDPdz) ? "" : std::to_string(r.predDPdz)) << ','
            << (r.blewUp ? 1 : 0) << '\n';
    }
    return 0;
}
