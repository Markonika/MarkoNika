// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Usage: paredes_report [dir_with_configs] [statics|stiffness|free|decay|waves|all] [con1|con2|cat] [inst|-] [coefficient_json|-] [period_s] [wave_height_m] [tension_csv]
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include "mooring/paredes.hpp"

using namespace mooring;
using json = nlohmann::json;

static json load(const std::string& dir, const std::string& name) {
    std::ifstream in(dir + "/" + name + ".json");
    if (!in) { std::fprintf(stderr, "cannot open %s/%s.json\n", dir.c_str(), name.c_str()); std::exit(2); }
    return json::parse(in);
}

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? argv[1] : "examples/paredes", what = argc > 2 ? argv[2] : "all";
    const bool all = what == "all";
    if (all || what == "free") {
        const json fb = load(dir, "free_buoy");
        std::printf("== free buoy decay (potential-theory A, B; measured heave 1.112+-0.006 s, pitch 1.170+-0.005 s)\n");
        std::printf("   heave T_d = %.4f s,  pitch T_d = %.4f s\n", paredesDecayPeriod(fb, 2, 0.02, 20.0, 1e-3), paredesDecayPeriod(fb, 4, 0.05, 20.0, 1e-3));
    }
    json coef; const bool haveCoef = argc > 5 && std::string(argv[5]) != "-";
    if (haveCoef) { std::ifstream ci(argv[5]); if (!ci) { std::fprintf(stderr, "cannot open %s\n", argv[5]); return 2; } coef = json::parse(ci); }
    for (const char* n : {"con1", "con2", "cat"}) {
        if (argc > 3 && std::string(argv[3]) != n) continue;               // optional: run a single configuration
        const json cfg = load(dir, n);
        if (all || what == "statics") {
            const ParedesStatics s = paredesStatics(cfg);
            std::printf("== %s statics: converged %d, xi = [%.5f %.5f %.5f | %.5f %.5f %.5f], draft change %+.4f m\n   leg tension at fairlead (N):", n, s.eq.converged,
                        s.xi[0], s.xi[1], s.xi[2], s.xi[3], s.xi[4], s.xi[5], s.draftChange);
            for (double t : s.topTension) std::printf(" %.3f", t);
            std::printf("   |force on buoy|:"); for (double t : s.fairleadForce) std::printf(" %.3f", t);
            std::printf("\n");
        }
        if (all || what == "stiffness") {
            const ParedesStiffness k = paredesSurgeStiffness(cfg, {-0.1, -0.05, 0.05, 0.1, 0.2});
            std::printf("== %s surge: F0 = %.3f N;  x [m] / Fx [N] / K [N/m]:", n, k.F0);
            for (size_t i = 0; i < k.x.size(); ++i) std::printf("  (%.2f %.3f %.2f%s)", k.x[i], k.Fx[i], k.K[i], k.converged[i] ? "" : "!");
            std::printf("\n   heave Km = %.2f N/m, pitch Km = %.3f N m/rad\n", paredesHeaveStiffness(cfg), paredesPitchStiffness(cfg));
        }
        if (what == "waves") {
            std::vector<double> periods = {1.30, 1.40};
            if (argc > 6) periods = {std::atof(argv[6])};                                // optional: a single period (any period present in the coefficient file)
            for (double T : periods) {
                CoupledResult detail;
                const ParedesRAO r = paredesWaveRAO(cfg, T, argc > 7 ? std::atof(argv[7]) : 0.08, 30.0, 15, 2e-3, true, argc > 4 && std::string(argv[4]) == "inst", haveCoef ? &coef : nullptr, argc > 8 ? &detail : nullptr);
                std::printf("== %s regular waves H = %.3f m, T = %.2f s: RAO surge %.3f, heave %.3f, pitch/(ka) %.3f (mean surge offset %.4f m)\n", n, argc > 7 ? std::atof(argv[7]) : 0.08, T, r.surge, r.heave, r.pitch, r.meanSurge);
                if (argc > 8) {
                    // Line-tension time series (all samples) and a summary over the last 15 wave cycles: slack means a segment with strain <= 0.
                    std::ofstream o(argv[8]);
                    o << "time_s,surge,heave,pitch";
                    for (const auto& ln : detail.lineNames) o << "," << ln << "_top_N," << ln << "_min_segment_N," << ln << "_min_segment," << ln << "_slack_segments";
                    o << "\n";
                    const double t0 = detail.samples.back().t - 15.0 * T;
                    std::vector<double> topMin(detail.lineNames.size(), 1e300), topMax(detail.lineNames.size(), 0.0), segMin(detail.lineNames.size(), 1e300);
                    std::vector<int> maxSlack(detail.lineNames.size(), 0); std::vector<double> slackTime(detail.lineNames.size(), 0.0); double span = 0.0;
                    for (size_t i = 0; i < detail.samples.size(); ++i) {
                        const auto& sm = detail.samples[i];
                        o << sm.t << "," << sm.xi[0] << "," << sm.xi[2] << "," << sm.xi[4];
                        for (size_t k = 0; k < detail.lineNames.size(); ++k) o << "," << sm.topTension[k] << "," << sm.minTension[k] << "," << sm.minSegment[k] << "," << sm.numSlack[k];
                        o << "\n";
                        if (sm.t < t0) continue;
                        const double dt = i + 1 < detail.samples.size() ? detail.samples[i + 1].t - sm.t : 0.0; span += dt;
                        for (size_t k = 0; k < detail.lineNames.size(); ++k) {
                            topMin[k] = std::min(topMin[k], sm.topTension[k]); topMax[k] = std::max(topMax[k], sm.topTension[k]);
                            segMin[k] = std::min(segMin[k], sm.minTension[k]); maxSlack[k] = std::max(maxSlack[k], sm.numSlack[k]);
                            if (sm.numSlack[k] > 0) slackTime[k] += dt;
                        }
                    }
                    for (size_t k = 0; k < detail.lineNames.size(); ++k)
                        std::printf("   %s: fairlead tension %.3f .. %.3f N, smallest segment tension %.4f N, max slack segments %d, time with any slack segment %.1f %% of the last 15 cycles\n",
                                    detail.lineNames[k].c_str(), topMin[k], topMax[k], segMin[k], maxSlack[k], span > 0 ? 100.0 * slackTime[k] / span : 0.0);
                }
                std::fflush(stdout);
            }
        }
        if (what == "pitchcurve") {
            const std::vector<double> th = {0.0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, -0.02, -0.05, -0.1, -0.15, -0.2, -0.25};
            const std::vector<double> M = paredesPitchCurve(cfg, th);
            std::printf("== %s line pitch moment at held pitch (others relaxed): theta [rad] / M [N m] / secant stiffness -(M-M0)/theta [N m/rad]\n", n);
            for (size_t i = 0; i < th.size(); ++i) std::printf("   %+.2f  %+.4f  %s\n", th[i], M[i], i == 0 ? "-" : (std::to_string(-(M[i] - M[0]) / th[i])).c_str());
        }
        if (all || what == "decay") {
            std::printf("== %s surge decay: T_d = %.3f s\n", n, paredesDecayPeriod(cfg, 0, 0.1, 45.0, 2e-3));
        }
        std::fflush(stdout);
    }
    return 0;
}
