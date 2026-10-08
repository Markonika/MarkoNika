// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Usage: paredes_report [dir_with_configs] [statics|stiffness|free|decay|waves|all] [con1|con2|cat]
#include <cstdio>
#include <fstream>
#include <string>
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
            for (double T : {1.30, 1.40}) {
                const ParedesRAO r = paredesWaveRAO(cfg, T, 0.08, 30.0, 15);
                std::printf("== %s regular waves H = 0.08 m, T = %.2f s: RAO surge %.3f, heave %.3f, pitch/(ka) %.3f (mean surge offset %.4f m)\n", n, T, r.surge, r.heave, r.pitch, r.meanSurge);
                std::fflush(stdout);
            }
        }
        if (all || what == "decay") {
            std::printf("== %s surge decay: T_d = %.3f s\n", n, paredesDecayPeriod(cfg, 0, 0.1, 45.0, 2e-3));
        }
        std::fflush(stdout);
    }
    return 0;
}
