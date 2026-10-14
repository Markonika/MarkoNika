// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Adaptive (error-controlled) steps of the implicit generalised-alpha integrator (DynOptions::adaptive).
#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>
#include "doctest.h"
#include "mooring/case_runner.hpp"
#include "mooring/lumped_mass_cable.hpp"

using namespace mooring;
using json = nlohmann::json;

namespace {
const double kPi = 3.14159265358979323846;

struct ModeResult { double err; DynStats stats; };
// Undamped first longitudinal mode of a taut fixed-fixed chain: the discrete mode is exactly u(t) = u0 sin(pi s/L) cos(omega_d t) (omega_d of the chain).
// Returns the max node error at t = tEnd relative to u0 for the given options (max step = factor x the CFL step).
ModeResult mode(const DynOptions& o, double tEnd) {
    const double L = 10.0, EA = 1e5, ml = 0.1, e0 = 1e-3, u0 = 1e-5; const int N = 40;
    CableParams p; p.L = L; p.N = N; p.EA = EA; p.m_l = ml; p.w = 0.0; p.g = 0.0;
    const double l0 = p.l0();
    LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(L * (1 + e0), 0, 0));
    std::vector<Vec3> r(N + 1), v(N + 1);
    for (int i = 0; i <= N; ++i) r[i] = Vec3(i * l0 * (1 + e0) + u0 * std::sin(kPi * i * l0 / L), 0, 0);
    c.setInitialState(r, v);
    c.setDynOptions(o);
    c.advanceTo(tEnd);
    const double om = 2.0 * std::sqrt(EA / ml) / l0 * std::sin(kPi * l0 / (2.0 * L));
    double e = 0;
    for (int i = 0; i <= N; ++i) {
        const double ex = i * l0 * (1 + e0) + u0 * std::sin(kPi * i * l0 / L) * std::cos(om * tEnd);
        e = std::max(e, std::fabs(c.nodes()[i].x - ex));
    }
    return {e / u0, c.stats()};
}

json loadLopez() {
    std::ifstream in(std::string(MOORING_SOURCE_DIR) + "/examples/lopezolocco/wocw.json");
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

TEST_CASE("Adaptive: the error shrinks as the tolerance is tightened (longitudinal mode riding on a larger static strain, maximum step 40 x CFL)") {
    DynOptions o; o.scheme = Scheme::GenAlpha; o.rhoInf = 1.0; o.implicitDtFactor = 40.0; o.adaptive = true;
    double prev = 0;
    for (double rtol : {1e-4, 1e-5, 1e-6}) {   // the mode strain (3e-6) is small against the static strain (1e-3), the scale of the tolerance
        o.adaptRtol = rtol;
        const ModeResult m = mode(o, 0.2);
        std::printf("\n  rtol %.0e: max error %.3e u0, %ld steps (%ld rejected), dt %.2e..%.2e", rtol, m.err, m.stats.steps, m.stats.rejectedSteps, m.stats.dtMinUsed, m.stats.dtMaxUsed);
        if (prev > 0) CHECK(m.err < 0.5 * prev);
        prev = m.err;
    }
    // Adaptive stepping respects the maximum step.
    CHECK(mode(o, 0.2).stats.dtMaxUsed <= 40.0 * 0.5 * 0.25 / 1000.0 * 1.0001);
}

TEST_CASE("Adaptive: Lopez-Olocco clump-free chain converges to the explicit result and refines steps around slack events") {
    json c = loadLopez();
    c["motion"]["radius_m"] = 0.225; c["motion"]["period_s"] = 2.8; c["motion"]["cycles"] = 12;
    c["numerics"]["implicit_damping"] = true;
    const double ref = lastCyclesMax(runCase(c));
    json g = c; g["numerics"]["scheme"] = "genalpha"; g["numerics"]["implicit_dt_factor"] = 20.0; g["numerics"]["adaptive"] = true;
    double prevd = 1e9;
    for (double rtol : {1e-3, 1e-4}) {
        g["numerics"]["adapt_rtol"] = rtol;
        const CaseResult r = runCase(g);
        const double m = lastCyclesMax(r), d = std::fabs(m - ref) / ref;
        std::printf("\n  rtol %.0e: max %.4f N vs RK4 %.4f N (%+.4f %%), %ld steps, %ld rejected, dt %.1e..%.1e s\n", rtol, m, ref, 100 * (m - ref) / ref, r.steps, r.rejectedSteps, r.dtMinUsed, r.dtMaxUsed);
        CHECK(r.finite);
        CHECK(d < 0.001);
        CHECK(r.dtMinUsed < 0.05 * r.dtMaxUsed);          // the step is refined around slack / snap events
        CHECK(r.newtonFailures == 0);
        prevd = d;
    }
    CHECK(prevd < 5e-4);
}
