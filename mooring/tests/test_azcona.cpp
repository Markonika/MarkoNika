// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Validation against Azcona, Munduate, Gonzalez & Nygaard (2017), Ocean Engineering 129:415-427: submerged 21 m chain at ECN,
// fairlead driven horizontally, two anchor distances, three periods. Static tensions from Table 5; dynamic maxima read by eye from
// Figs 5 and 6 (mean markers, +-0.5 N). See docs/azcona_validation.md.
#include <cmath>
#include <cstdio>
#include <fstream>
#include "doctest.h"
#include "mooring/case_runner.hpp"

using namespace mooring;
using json = nlohmann::json;

namespace {
json load(const char* conf) {
    std::ifstream in(std::string(MOORING_SOURCE_DIR) + "/examples/azcona/" + conf + ".json");
    REQUIRE(in.good());
    json c = json::parse(in);
    c["output"]["write"] = false;
    return c;
}
// Mean of the per-step maxima of the last 4 complete cycles (the same definition as scripts/azcona_compare.py).
double lastCyclesMax(const CaseResult& r) {
    const int n = int(r.cycles.size()) - 1;                // the last entry is the single step at t = nCycles*T
    double s = 0; for (int k = n - 4; k < n; ++k) s += r.cycles[k].topMax;
    return s / 4.0;
}
double dyn(const char* conf, double T, int cycles) {
    json c = load(conf); c["motion"]["period_s"] = T; c["motion"]["cycles"] = cycles;
    return lastCyclesMax(runCase(c));
}
}  // namespace

TEST_CASE("Azcona 2017: static fairlead tension vs Table 5 (A priori criterion: within 5 %)") {
    struct S { const char* conf; double meas, paperModel; };
    for (S s : {S{"conf1", 8.13, 8.10}, S{"conf2", 14.48, 14.70}}) {
        json c = load(s.conf); c["motion"]["type"] = "none";
        const CaseResult r = runCase(c);
        std::printf("\n  %s: static top tension %.3f N; measured %.2f, the paper's own model %.2f (%+.2f %% vs measured)\n", s.conf,
                    r.staticTopTension, s.meas, s.paperModel, 100 * (r.staticTopTension - s.meas) / s.meas);
        CHECK(r.staticConverged);
        CHECK(std::fabs(r.staticTopTension - s.meas) / s.meas < 0.05);
    }
}

TEST_CASE("Azcona 2017: configuration 1, T = 4.74 s (harmonic regime) maximum tension vs Fig. 5") {
    const double m = dyn("conf1", 4.74, 8);
    std::printf("\n  conf1 T=4.74: model max %.2f N, measured ~9.2 N (%+.1f %%)\n", m, 100 * (m - 9.2) / 9.2);
    CHECK(std::fabs(m - 9.2) / 9.2 < 0.08);
}

TEST_CASE("Azcona 2017: all six dynamic cases, maximum fairlead tension vs Figs 5, 6 [validation run, ~1 min]" * doctest::skip()) {
    struct C { const char* conf; double T, meas; bool snap; };
    // A priori criterion: maximum within 8 % of the value read from the figure. Configuration 2 at 1.58 s (total loss of tension) is the
    // paper's own acknowledged sensitive case; its asserted envelope is the documented spread of the model (docs/azcona_validation.md).
    for (C c : {C{"conf1", 1.58, 14.5, false}, C{"conf1", 3.16, 9.8, false}, C{"conf1", 4.74, 9.2, false},
                C{"conf2", 1.58, 44.5, true}, C{"conf2", 3.16, 23.5, false}, C{"conf2", 4.74, 19.5, false}}) {
        const double m = dyn(c.conf, c.T, 12);
        std::printf("\n  %s T=%.2f: model max %.2f N, measured ~%.1f N (%+.1f %%)%s\n", c.conf, c.T, m, c.meas, 100 * (m - c.meas) / c.meas, c.snap ? "  [snap case]" : "");
        if (!c.snap) CHECK(std::fabs(m - c.meas) / c.meas < 0.08);
        else { CHECK(m > 0.9 * c.meas); CHECK(m < 1.40 * c.meas); }
    }
}
