// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Milestone 6: point elements (clump weights, floaters) on nodes.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>
#include "doctest.h"
#include "mooring/case_runner.hpp"
#include "mooring/catenary.hpp"
#include "mooring/lumped_mass_cable.hpp"

using namespace mooring;
using json = nlohmann::json;

namespace {
const double kPi = 3.14159265358979323846;

// Exact elastic catenary of one half of a symmetric line carrying a point load W (down positive) at mid-span.
// Vertical force balance at the kink fixes V0 = -(w L + W)/2 for the half starting at the support (s = 0);
// H follows from x(L/2) = span/2.
ElasticCatenary halfCatenary(double L, double EA, double w, double span, double W) {
    ElasticCatenary c; c.L = 0.5 * L; c.EA = EA; c.w = w; c.V0 = -0.5 * (w * L + W);
    double lo = 1e-6, hi = 1e9;
    for (int i = 0; i < 300; ++i) {
        c.H = 0.5 * (lo + hi);
        double x, z; c.position(c.L, x, z);
        (x > 0.5 * span ? hi : lo) = c.H;
    }
    c.H = 0.5 * (lo + hi);
    return c;
}

struct PointRun { double l2, Herr, reactionErr, symErr; bool conv; long steps; };

PointRun symmetricRun(int N, double W, bool floater, double* Hth = nullptr) {
    const double L = 100.5, span = 100.0, EA = 200e3, ml = 1.738, g = 9.81;
    CableParams p; p.L = L; p.N = N; p.EA = EA; p.m_l = ml; p.g = g;
    Environment env;
    if (floater) {                                   // submerged: net weight from Eq. 3.26, buoyancy from the point volume
        env.hydro = true; env.surfaceZ = 1e3;
        p.w = CableParams::submergedWeight(ml, 7800.0, 1000.0, g); p.D0 = 0.01;
    } else {
        p.w = ml * g;                                // in air
    }
    LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(span, 0, 0));
    c.setEnvironment(env);
    // W > 0: clump of net (downward) weight W; W < 0: floater with net upward force |W| when added to the cable weight.
    if (floater) c.addPointElement(PointElement::floater(N / 2, 0.0, -W, 0.1, 1.0, 1000.0, g));
    else         c.addPointElement({N / 2, W / g, 0.0, 0.0, 0.0, 0.0});
    RelaxOptions ro; ro.forceTol = 1e-8; ro.maxSteps = 8000000;
    const RelaxResult rr = c.relaxStatic(ro);
    const ElasticCatenary cat = halfCatenary(L, EA, p.w, span, W);
    double num = 0, den = 0;
    for (int i = 0; i <= N; ++i) {
        const bool left = i <= N / 2;
        const double s = (left ? i : N - i) * p.l0();
        double x, z; cat.position(s, x, z);
        if (!left) x = span - x;
        const Vec3& r = c.nodes()[i];
        num += (r.x - x) * (r.x - x) + (r.z - z) * (r.z - z); den += x * x + z * z;
    }
    const Vec3 Ft = c.endForce(true), Fa = c.endForce(false);
    PointRun o;
    o.l2 = std::sqrt(num / den); o.conv = rr.converged; o.steps = rr.steps;
    o.Herr = std::fabs(-Ft.x - cat.H) / cat.H;
    o.reactionErr = std::fabs((Ft.z + Fa.z) + (p.w * L + W)) / std::fabs(p.w * L + W);   // total support load = total weight
    o.symErr = std::fabs(Ft.z - Fa.z) / std::fabs(Ft.z);
    if (Hth) *Hth = cat.H;
    return o;
}

void study(const char* label, double W, bool floater) {
    std::printf("\n  %s (net point load %+.1f N): L2 error of the shape vs the exact half-catenaries\n  %6s %12s %12s %12s %10s\n",
                label, W, "N", "L2 error", "H error", "reaction err", "steps");
    std::vector<double> l2;
    for (int N : {20, 40, 80, 160}) {
        const PointRun r = symmetricRun(N, W, floater);
        std::printf("  %6d %12.4e %12.4e %12.2e %10ld\n", N, r.l2, r.Herr, r.reactionErr, r.steps);
        CHECK(r.conv);
        CHECK(r.reactionErr < 1e-6);                 // sum of support loads = cable weight + point load (to relaxation accuracy)
        CHECK(r.symErr < 1e-6);                      // symmetric: both supports carry half
        l2.push_back(r.l2);
    }
    for (size_t i = 1; i < l2.size(); ++i) CHECK(l2[i] < l2[i - 1]);
    const double order = std::log(l2[0] / l2[3]) / std::log(8.0);
    std::printf("  observed order (N=20..160): %.2f\n", order);
    CHECK(order > 1.5);
    CHECK(l2.back() < 1e-4);
}
}  // namespace

TEST_CASE("clump weight at mid-span of a catenary vs exact elastic half-catenaries") { study("clump 300 N, in air", 300.0, false); }

TEST_CASE("floater at mid-span of a catenary (arch) vs exact elastic half-catenaries") { study("floater 2500 N buoyancy, submerged", -2500.0 + 0.0, true); }

TEST_CASE("point-element drag: horizontal support load equals Morison drag of the element") {
    const int N = 20; const double D = 0.1, Cd = 1.2, U = 1.5, rho = 1000.0, g = 9.81;
    CableParams p; p.L = 10.2; p.N = N; p.EA = 1e6; p.m_l = 0.5; p.g = g;
    p.w = CableParams::submergedWeight(p.m_l, 7800.0, rho, g); p.D0 = 0.01; p.Cdn = 0.0; p.Cdt = 0.0;   // cable drag off
    Environment env; env.hydro = true; env.surfaceZ = 1e3;
    env.water = [U](const Vec3&, double, Vec3& vw, Vec3& aw) { vw = Vec3(U, 0, 0); aw = Vec3(); };
    LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(10, 0, 0)); c.setEnvironment(env);
    c.addPointElement(PointElement::clump(N / 2, 1.0, 5.0, D, Cd, rho, g));
    RelaxOptions ro; ro.forceTol = 1e-8; REQUIRE(c.relaxStatic(ro).converged);
    const double A = kPi * D * D / 4.0, Fd = 0.5 * rho * Cd * A * U * U;
    const double Fx = c.endForce(true).x + c.endForce(false).x;
    std::printf("\n  point drag: expected %.6f N, support sum %.6f N, weight check z: %.6f vs %.6f N\n", Fd, Fx,
                c.endForce(true).z + c.endForce(false).z, -(p.w * p.L + 5.0));
    CHECK(Fx == doctest::Approx(Fd).epsilon(1e-6));
    CHECK(c.endForce(true).z + c.endForce(false).z == doctest::Approx(-(p.w * p.L + 5.0)).epsilon(1e-6));
}

TEST_CASE("point-element added mass: two-segment oscillator frequency") {
    // Gravity-free taut string, N = 2, point element on the middle node: omega^2 = 2 T / (l M), M = m_node + m_pt + Cm rho V.
    auto period = [](double Cm) {
        CableParams p; p.L = 2.0; p.N = 2; p.EA = 100.0; p.m_l = 0.05; p.g = 0.0; p.w = 0.0;
        Environment env; env.hydro = true; env.surfaceZ = 1e3;
        LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(2.02, 0, 0)); c.setEnvironment(env);
        c.addPointElement({1, 0.4, 1e-4, 0.0, 0.0, Cm});
        std::vector<Vec3> r = {Vec3(0, 0, 0), Vec3(1.01, 1e-4, 0), Vec3(2.02, 0, 0)}, v(3);
        c.setInitialState(r, v);
        DynOptions o; o.dt = 1e-4; c.setDynOptions(o);
        double tPrev = 0, yPrev = 1e-4;
        for (double t = 1e-3; t < 5.0; t += 1e-3) {
            c.advanceTo(t);
            const double y = c.nodes()[1].y;
            if (yPrev > 0 && y <= 0) return 4.0 * (tPrev + (t - tPrev) * yPrev / (yPrev - y));   // quarter period -> period
            tPrev = t; yPrev = y;
        }
        return 0.0;
    };
    // Tension T = EA*eps = 1 N (eps = 0.01), chord length per segment l = 1.01 m.
    const double l = 1.01, rhoV = 1000.0 * 1e-4, mNode = 0.05;
    const double w0 = std::sqrt(2.0 * 1.0 / (l * (mNode + 0.4))), w1 = std::sqrt(2.0 * 1.0 / (l * (mNode + 0.4 + rhoV)));
    const double T0 = period(0.0), T1 = period(1.0);
    std::printf("\n  added mass oscillator: period Cm=0 %.5f s (exact %.5f), Cm=1 %.5f s (exact %.5f), ratio %.5f (exact %.5f)\n",
                T0, 2 * kPi / w0, T1, 2 * kPi / w1, T1 / T0, w0 / w1);
    CHECK(T0 == doctest::Approx(2 * kPi / w0).epsilon(2e-3));
    CHECK(T1 == doctest::Approx(2 * kPi / w1).epsilon(2e-3));
    CHECK(T1 / T0 == doctest::Approx(w0 / w1).epsilon(1e-3));
}

TEST_CASE("point elements through the JSON runner (arc-length placement)") {
    std::ifstream in(std::string(MOORING_SOURCE_DIR) + "/examples/chalmers/chalmers_config.json");
    REQUIRE(in.good());
    json base = json::parse(in);
    base["output"]["write"] = false; base["motion"]["type"] = "none";
    const CaseResult r0 = runCase(base);
    json withClump = base;
    withClump["point_elements"] = json::array({ {{"type", "clump"}, {"arclength_m", 30.0}, {"mass_kg", 0.05},
                                                  {"submerged_weight_N", 0.3}, {"diameter_m", 0.02}, {"Cd", 1.0}} });
    const CaseResult r1 = runCase(withClump);
    json withFloat = base;
    withFloat["point_elements"] = json::array({ {{"type", "floater"}, {"arclength_m", 30.0}, {"mass_kg", 0.01},
                                                  {"buoyancy_N", 0.3}, {"diameter_m", 0.02}, {"Cd", 1.0}} });
    const CaseResult r2 = runCase(withFloat);
    std::printf("\n  static top tension: none %.4f, clump(+0.3 N) %.4f, floater(-0.3 N) %.4f N\n", r0.staticTopTension,
                r1.staticTopTension, r2.staticTopTension);
    CHECK(r1.staticConverged); CHECK(r2.staticConverged);
    CHECK(r1.staticTopTension > r0.staticTopTension + 0.05);                            // clump near the top adds tension
    CHECK(r2.staticTopTension < r0.staticTopTension - 0.05);                            // floater reduces it
}
