// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Non-uniform lines (LineSection): static equilibrium of a hanging two-section line, transverse eigenfrequency of a two-medium string,
// equivalence of a uniform line split into identical sections, and the JSON runner.
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
}

TEST_CASE("Sections: vertical two-section line (different EA, mass, weight, segment length) hangs exactly as the recurrence T_i = T_0 + sum W") {
    CableParams p; p.g = 9.81;
    LineSection a, b;
    a.length = 4.0; a.segments = 8;  a.EA = 2e4; a.m_l = 0.5; a.w = 4.0;
    b.length = 6.0; b.segments = 20; b.EA = 6e4; b.m_l = 0.2; b.w = 1.5;
    p.sections = {a, b};
    // Choose the anchor force T0 and derive the fairlead height from the nodal equilibrium recurrence (linear elastic segments).
    const double T0 = 60.0;
    const int N = 28;
    std::vector<double> l0(N), EA(N), wseg(N);
    for (int j = 0; j < N; ++j) { const LineSection& q = j < 8 ? a : b; l0[j] = q.length / q.segments; EA[j] = q.EA; wseg[j] = q.w; }
    double T = T0, h = 0.0;
    for (int j = 0; j < N; ++j) {
        h += l0[j] * (1.0 + T / EA[j]);
        if (j + 1 < N) T += 0.5 * (wseg[j] * l0[j] + wseg[j + 1] * l0[j + 1]);   // node j+1 weight
    }
    const double Ttop = T;
    LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(0, 0, h));
    CHECK(c.params().N == 28);
    CHECK(c.params().L == doctest::Approx(10.0));
    std::vector<Vec3> r(N + 1), v(N + 1);
    double z = 0.0;
    for (int j = 0; j <= N; ++j) { r[j] = Vec3(0, 0, z); if (j < N) z += l0[j]; }
    for (Vec3& q : r) q.z *= h / z;                                                  // spread the nodes between the fixed ends; relaxation finds the strains
    c.setInitialState(r, v);
    RelaxOptions ro; ro.forceTol = 1e-8;
    const RelaxResult rr = c.relaxStatic(ro);
    CHECK(rr.converged);
    CHECK(c.endSegmentTension(false) == doctest::Approx(T0).epsilon(1e-6));
    CHECK(c.endSegmentTension(true) == doctest::Approx(Ttop).epsilon(1e-6));
    // End-node forces (tension + own half weight): vertical sum equals minus the total weight, sum_j w_j l0_j (as in test_hydro).
    double Wtot = 0.0; for (int j = 0; j < N; ++j) Wtot += wseg[j] * l0[j];
    CHECK(c.endForce(false).z + c.endForce(true).z == doctest::Approx(-Wtot).epsilon(1e-6));
    CHECK(c.nodeAtArclength(4.0) == 8);
    CHECK(c.nodeAtArclength(4.1) == 8);
    CHECK(c.nodeAtArclength(4.4) == 9);
    CHECK(c.nodeAtArclength(10.0) == 28);
}

TEST_CASE("Sections: lowest transverse frequency of a taut two-medium string matches the exact eigenvalue equation") {
    // Fixed-fixed, equal tension T, stretched lengths S_k = L_k (1 + T/EA_k), line density rho_k = m_k / (1 + T/EA_k):
    // k1 cot(k1 S1) + k2 cot(k2 S2) = 0 with k_k = w sqrt(rho_k / T).
    CableParams p; p.g = 0.0;
    LineSection a, b;
    a.length = 4.0; a.segments = 80;  a.EA = 1e5; a.m_l = 0.1; a.w = 0.0;
    b.length = 6.0; b.segments = 120; b.EA = 3e5; b.m_l = 0.3; b.w = 0.0;
    p.sections = {a, b};
    const double T = 200.0;
    const double e1 = T / a.EA, e2 = T / b.EA, S1 = a.length * (1 + e1), S2 = b.length * (1 + e2);
    const double rho1 = a.m_l / (1 + e1), rho2 = b.m_l / (1 + e2);
    auto F = [&](double w) {
        const double k1 = w * std::sqrt(rho1 / T), k2 = w * std::sqrt(rho2 / T);
        return k1 / std::tan(k1 * S1) + k2 / std::tan(k2 * S2);
    };
    // Bracket the lowest root: scan upward from small w for the first sign change of F that is not a pole.
    double w0 = 0.5, step = 0.01, lo = 0, hi = 0;
    for (double w = w0; w < 200.0; w += step) {
        if (F(w) * F(w + step) < 0.0 && std::fabs(F(w)) < 50.0 && std::fabs(F(w + step)) < 50.0) { lo = w; hi = w + step; break; }
    }
    REQUIRE(hi > 0.0);
    for (int i = 0; i < 80; ++i) { const double m = 0.5 * (lo + hi); (F(lo) * F(m) <= 0.0 ? hi : lo) = m; }
    const double wExact = 0.5 * (lo + hi);
    // Simulation: pretension by stretching between fixed ends, small transverse perturbation, time from first to last upward zero crossing.
    const double h = S1 + S2;
    LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(h, 0, 0));
    const int N = 200;
    std::vector<Vec3> r(N + 1), v(N + 1);
    const double k1 = wExact * std::sqrt(rho1 / T), k2 = wExact * std::sqrt(rho2 / T), Bamp = std::sin(k1 * S1) / std::sin(k2 * S2);
    double x = 0.0;
    for (int j = 0; j <= N; ++j) {
        const double y = x <= S1 ? std::sin(k1 * x) : Bamp * std::sin(k2 * (h - x));   // exact first mode (amplitude 1e-4)
        r[j] = Vec3(x, 1e-4 * y, 0.0);
        if (j < N) x += (j < 80 ? a.length / 80 * (1 + e1) : b.length / 120 * (1 + e2));
    }
    c.setInitialState(r, v);
    const double dt = c.stableDt();
    std::vector<double> crossings; double yPrev = c.nodes()[60].y, tPrev = 0.0;
    while (c.time() < 3.0 * 2.0 * kPi / wExact && crossings.size() < 4) {
        c.advanceTo(c.time() + dt);
        const double y = c.nodes()[60].y;
        if (yPrev < 0.0 && y >= 0.0) crossings.push_back(tPrev + (c.time() - tPrev) * (-yPrev) / (y - yPrev));
        yPrev = y; tPrev = c.time();
    }
    REQUIRE(crossings.size() >= 3);
    const double period = (crossings.back() - crossings.front()) / double(crossings.size() - 1);
    const double wSim = 2.0 * kPi / period;
    std::printf("\n  two-medium string: omega exact %.4f, simulated %.4f (%+.3f %%)\n", wExact, wSim, 100 * (wSim - wExact) / wExact);
    CHECK(wSim == doctest::Approx(wExact).epsilon(3e-3));
}

TEST_CASE("Sections: a uniform line split into identical sections reproduces the uniform line (JSON runner, Azcona set-up)") {
    std::ifstream in(std::string(MOORING_SOURCE_DIR) + "/examples/azcona/conf1.json");
    REQUIRE(in.good());
    json c = json::parse(in);
    c["output"]["write"] = false;
    c["motion"]["period_s"] = 3.16; c["motion"]["cycles"] = 4;
    json cs = c;
    const json line = c.at("line");
    json s1 = json::object(), s2 = json::object();
    s1["length_m"] = 0.5 * line.at("length_m").get<double>(); s1["segments"] = line.at("segments").get<int>() / 2;
    s2 = s1;
    cs["line"].erase("length_m"); cs["line"].erase("segments");
    cs["line"]["sections"] = json::array({s1, s2});
    const CaseResult r0 = runCase(c), r1 = runCase(cs);
    std::printf("\n  static top tension uniform %.9f, sectioned %.9f\n", r0.staticTopTension, r1.staticTopTension);
    CHECK(r1.staticTopTension == doctest::Approx(r0.staticTopTension).epsilon(1e-8));
    CHECK(r1.cycles.back().topMax == doctest::Approx(r0.cycles.back().topMax).epsilon(1e-6));
    CHECK(r1.cycles[2].topMax == doctest::Approx(r0.cycles[2].topMax).epsilon(1e-6));
}
