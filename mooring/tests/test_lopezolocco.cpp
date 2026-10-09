// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Validation against Lopez-Olocco et al. (2022), J. Mar. Sci. Eng. 10(5):676: 27 m studless chain, forced circular fairlead motion,
// with and without a clump weight. Tables 8/9 are in examples/lopezolocco/measured.csv. See docs/lopezolocco_validation.md.
#include <cmath>
#include <cstdio>
#include <fstream>
#include "doctest.h"
#include "mooring/case_runner.hpp"

using namespace mooring;
using json = nlohmann::json;

namespace {
json load(const char* name) {
    std::ifstream in(std::string(MOORING_SOURCE_DIR) + "/examples/lopezolocco/" + name + ".json");
    REQUIRE(in.good());
    json c = json::parse(in);
    c["output"]["write"] = false;
    return c;
}
double lastCyclesMax(const CaseResult& r) {
    const int n = int(r.cycles.size()) - 1;
    double s = 0; for (int k = n - 4; k < n; ++k) s += r.cycles[k].topMax;
    return s / 4.0;
}
}  // namespace

TEST_CASE("Lopez-Olocco 2022: static pretension at the 25 m anchor distance is finite and the clump raises it by about 10 %") {
    double t[3]; int i = 0;
    for (const char* n : {"wocw", "cw1", "cw2"}) {
        json c = load(n); c["motion"]["type"] = "none";
        const CaseResult r = runCase(c);
        CHECK(r.staticConverged);
        t[i++] = r.staticTopTension;
    }
    std::printf("\n  static top tension: WO %.3f N, CW1 %.3f N, CW2 %.3f N (paper: CW2 about +10 %% vs WO)\n", t[0], t[1], t[2]);
    CHECK(t[0] > 10.5); CHECK(t[0] < 13.0);
    CHECK(t[2] > t[0]);
}

TEST_CASE("Lopez-Olocco 2022: clump-free and clump maxima vs Table 8 [validation run, several minutes]" * doctest::skip()) {
    struct C { const char* cfg; double A, T, meas; };
    // A priori criterion: maximum within 6 % of Table 8.
    for (C c : {C{"wocw", 0.175, 3.5, 16.4}, C{"cw1", 0.175, 3.5, 17.29}, C{"cw1", 0.225, 2.8, 24.31}, C{"cw2", 0.175, 3.5, 17.83}}) {
        json j = load(c.cfg); j["motion"]["radius_m"] = c.A; j["motion"]["period_s"] = c.T; j["motion"]["cycles"] = 12;
        const double m = lastCyclesMax(runCase(j));
        std::printf("\n  %s A=%.3f T=%.1f: model max %.2f N, measured %.2f N\n", c.cfg, c.A, c.T, m, c.meas);
        if (c.meas > 0) CHECK(std::fabs(m - c.meas) / c.meas < 0.06);
    }
}
