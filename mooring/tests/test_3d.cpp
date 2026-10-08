// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Milestone 5: planar (y = 0) vs 3D regression, rotation invariance, mirror symmetry, out-of-plane behaviour.
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

json chalmers(double r, double T, int cycles = 5) {
    std::ifstream in(std::string(MOORING_SOURCE_DIR) + "/examples/chalmers/chalmers_config.json");
    REQUIRE(in.good());
    json c = json::parse(in);
    c["output"]["write"] = false;
    c["motion"]["radius_m"] = r; c["motion"]["period_s"] = T; c["motion"]["cycles"] = cycles;
    c["statistics"]["first_cycle"] = 2;
    return c;
}
struct Case { double r, T, rotTol; };
}  // namespace

TEST_CASE("planar mode reproduces the 3D solution exactly for in-plane problems (Chalmers)") {
    std::printf("\n  planar vs 3D (mean max top tension, N)\n  %5s %5s %14s %14s %10s\n", "r", "T", "3D", "planar", "diff");
    for (Case c : {Case{0.2, 3.5, 0}, Case{0.1, 2.0, 0}, Case{0.2, 1.25, 0}}) {
        json a = chalmers(c.r, c.T); json b = a; b["numerics"]["planar"] = true;
        const CaseResult ra = runCase(a), rb = runCase(b);
        std::printf("  %5.2f %5.2f %14.9f %14.9f %10.1e\n", c.r, c.T, ra.meanMax, rb.meanMax, ra.meanMax - rb.meanMax);
        CHECK(std::fabs(ra.meanMax - rb.meanMax) < 1e-9);
        CHECK(std::fabs(ra.staticTopTension - rb.staticTopTension) < 1e-9);
        CHECK(rb.maxOutOfPlane == 0.0);
        CHECK(ra.maxOutOfPlane < 1e-12);               // 3D run with in-plane data stays in plane
    }
}

TEST_CASE("3D rotation invariance: whole set-up rotated 37 deg about the vertical axis") {
    const double al = 37.0 * kPi / 180.0;
    std::printf("\n  rotated vs unrotated (mean max top tension, N)\n  %5s %5s %14s %14s %10s\n", "r", "T", "x-z plane", "rotated", "rel diff");
    for (Case c : {Case{0.2, 3.5, 1e-4}, Case{0.1, 2.0, 1e-4}, Case{0.2, 1.25, 1e-2}}) {
        // Snap case: round-off differences between the two runs are amplified by slack/snap events; the tolerance is
        // the documented snap-case scatter of milestone 4 (about +/-1 %), observed 0.1-0.5 % depending on the window.
        json a = chalmers(c.r, c.T); json b = a;
        b["fairlead_rest_m"] = {32.554 * std::cos(al), 32.554 * std::sin(al), 3.3};
        b["motion"]["plane_angle_deg"] = 37.0;
        const CaseResult ra = runCase(a), rb = runCase(b);
        const double rel = std::fabs(ra.meanMax - rb.meanMax) / ra.meanMax;
        std::printf("  %5.2f %5.2f %14.6f %14.6f %10.1e\n", c.r, c.T, ra.meanMax, rb.meanMax, rel);
        CHECK(rel < c.rotTol);
        CHECK(std::fabs(ra.staticTopTension - rb.staticTopTension) < 1e-5);           // limited by the relaxation tolerance
        CHECK(rb.maxOutOfPlane < 1e-9);                                                 // stays in the rotated plane
    }
}

TEST_CASE("mirror symmetry: lateral current +U and -U give mirrored deflection and equal tension") {
    json a = chalmers(0.05, 3.5, 5);
    a["environment"]["current_m_s"] = {0.0, 0.3, 0.0};
    json b = a; b["environment"]["current_m_s"] = {0.0, -0.3, 0.0};
    json n = a; n["environment"]["current_m_s"] = {0.0, 0.0, 0.0};                      // no current
    const CaseResult ra = runCase(a), rb = runCase(b), rn = runCase(n);
    std::printf("\n  lateral current 0.3 m/s: max tension +U %.9f, -U %.9f, none %.6f N; max |y| %.4f / %.4f m\n",
                ra.meanMax, rb.meanMax, rn.meanMax, ra.maxOutOfPlane, rb.maxOutOfPlane);
    CHECK(std::fabs(ra.meanMax - rb.meanMax) < 1e-9);
    CHECK(std::fabs(ra.maxOutOfPlane - rb.maxOutOfPlane) < 1e-12);
    CHECK(ra.maxOutOfPlane > 0.05);                                                      // the current really pushes it out of plane
    CHECK(rn.maxOutOfPlane < 1e-12);
    CHECK(ra.meanMax > rn.meanMax);                                                      // sideways drag adds tension
}

TEST_CASE("out-of-plane perturbation of the Chalmers solution does not grow") {
    std::printf("\n  1 mm half-sine out-of-plane perturbation, 6 cycles\n");
    for (Case c : {Case{0.2, 3.5, 0}, Case{0.1, 2.0, 0}, Case{0.2, 1.25, 0}}) {
        json a = chalmers(c.r, c.T, 6); json p = a; p["initial"]["perturbation_y_m"] = 1e-3;
        const CaseResult ra = runCase(a), rp = runCase(p);
        std::printf("  r=%.2f T=%.2f: max |y| %.3e m, final |y| %.3e m (initial 1.0e-3), mean max tension %.4f vs %.4f N\n",
                    c.r, c.T, rp.maxOutOfPlane, rp.finalOutOfPlane, rp.meanMax, ra.meanMax);
        CHECK(rp.finalOutOfPlane <= 1.0e-3);
        CHECK(rp.maxOutOfPlane <= 1.0e-3 * 1.001);                                       // never exceeds the initial amplitude
        CHECK(std::fabs(rp.meanMax - ra.meanMax) / ra.meanMax < 0.015);                  // within the snap-case scatter
    }
}

TEST_CASE("circularly polarised string: constant radial amplitude, uniform phase rotation (3D analytic)") {
    // Same string as the milestone-2 test (L0 = 0.5, eps = 1, EA = 1 N, m_l = 1 kg/m, no gravity), but displaced as
    // y = A cos(wt) s(x), z = A sin(wt) s(x). In the linear limit each node moves on a circle of radius A s(x).
    const int N = 40; const double L0 = 0.5, eps = 1.0, A = 1e-4, tEnd = 1.3;
    CableParams p; p.L = L0; p.N = N; p.EA = 1.0; p.m_l = 1.0; p.w = 0.0;
    LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(L0 * (1 + eps), 0, 0));
    std::vector<Vec3> r(N + 1), v(N + 1);
    const double cT = std::sqrt(1.0 / (p.m_l * (1 + eps)));
    const double om = 2.0 / p.l0() * cT * std::sin(kPi / (2.0 * N));                  // discrete frequency of the LM string
    for (int i = 0; i <= N; ++i) {
        const double sh = std::sin(kPi * i / N);
        r[i] = Vec3(i * p.l0() * (1 + eps), A * sh, 0.0);
        v[i] = Vec3(0.0, 0.0, A * om * sh);                                            // released with the matching velocity
    }
    c.setInitialState(r, v);
    DynOptions o; o.dt = 1e-4; c.setDynOptions(o);
    c.advanceTo(tEnd);
    double maxRadErr = 0, maxPhErr = 0;
    for (int i = 1; i < N; ++i) {
        const double sh = std::sin(kPi * i / N);
        const Vec3& q = c.nodes()[i];
        maxRadErr = std::max(maxRadErr, std::fabs(std::hypot(q.y, q.z) - A * sh) / (A * sh));
        double ph = std::atan2(q.z, q.y) - om * tEnd;
        ph = std::remainder(ph, 2.0 * kPi);
        maxPhErr = std::max(maxPhErr, std::fabs(ph));
    }
    std::printf("\n  polarised string, t = %.1f s (%.2f periods): max radial error %.2e, max phase error %.2e rad\n",
                tEnd, om * tEnd / (2 * kPi), maxRadErr, maxPhErr);
    CHECK(maxRadErr < 1e-4);
    CHECK(maxPhErr < 1e-3);
}
