// Copyright (c) [OWNER NAME]. All rights reserved. Proprietary and confidential.
// Test 1 (static catenary): 100.5 m cable hung between points 100 m apart at equal height,
// m_l = 1.738 kg/m (Paredes 2016, sec. 3.5.1). Settled LM shape vs. Eqs. 3.54-3.56.
#include <cmath>
#include <cstdio>
#include <vector>
#include "doctest.h"
#include "mooring/catenary.hpp"
#include "mooring/lumped_mass_cable.hpp"

using namespace mooring;

namespace {
constexpr double kL = 100.5, kSpan = 100.0, kML = 1.738, kG = 9.81;

// Relative L2 position error over the nodes (Lagrangian coordinate s_i = i*l0) and the
// error of the end tension vs. the analytic catenary.
struct Err { double l2, tensionErr; RelaxResult rr; };

Err runCase(int N, double EA) {
    const double w = kML * kG;
    CableParams p; p.L = kL; p.N = N; p.EA = EA; p.m_l = kML; p.w = w;
    LumpedMassCable cable(p, Vec3(0, 0, 0), Vec3(kSpan, 0, 0));
    RelaxOptions opt; opt.forceTol = 1e-7; opt.maxSteps = 3000000;
    Err e; e.rr = cable.relaxStatic(opt);
    const ElasticCatenary cat = ElasticCatenary::solve(kL, EA, w, kSpan, 0.0);
    REQUIRE(cat.converged);
    double num = 0, den = 0;
    for (int i = 0; i <= N; ++i) {
        double x, z; cat.position(i * p.l0(), x, z);
        const Vec3& r = cable.nodes()[i];
        num += (r.x - x) * (r.x - x) + (r.z - z) * (r.z - z);
        den += x * x + z * z;
    }
    e.l2 = std::sqrt(num / den);
    e.tensionErr = std::fabs(norm(cable.endTension(false)) - cat.tension(0)) / cat.tension(0);
    return e;
}
}  // namespace

TEST_CASE("elastic catenary solution is self-consistent") {
    const double w = kML * kG;
    for (double EA : {200e3, 2.0e7}) {
        const auto c = ElasticCatenary::solve(kL, EA, w, kSpan, 0.0);
        CHECK(c.converged);
        double x, z; c.position(kL, x, z);
        CHECK(x == doctest::Approx(kSpan).epsilon(1e-12));
        CHECK(std::fabs(z) < 1e-9);
        CHECK(c.V0 == doctest::Approx(-0.5 * w * kL).epsilon(1e-9));  // symmetric: V0 = -wL/2
        double xm, zm; c.position(0.5 * kL, xm, zm);                  // mid-point at x = span/2
        CHECK(xm == doctest::Approx(0.5 * kSpan).epsilon(1e-9));
    }
}

TEST_CASE("catenary with unequal end heights satisfies both end conditions") {
    const auto c = ElasticCatenary::solve(60.0, 5e6, 20.0, 40.0, 15.0);
    REQUIRE(c.converged);
    double x, z; c.position(60.0, x, z);
    CHECK(x == doctest::Approx(40.0).epsilon(1e-12));
    CHECK(z == doctest::Approx(15.0).epsilon(1e-10));
}

TEST_CASE("static LM cable converges to the elastic catenary; L2 error vs N") {
    // K = 200 kN/m in the brief is interpreted as EA = 200 kN (see docs/assumptions.md); the stiff
    // case EA = K*L is also run to separate discretisation error from elastic effects.
    for (double EA : {200e3, 200e3 * kL}) {
        std::printf("\n  EA = %.4g N\n  %6s %12s %12s %10s %8s\n", EA, "N", "L2 error", "T_end err", "steps", "conv");
        std::vector<double> l2;
        const std::vector<int> Ns = {10, 20, 40, 80, 160};
        for (int N : Ns) {
            const Err e = runCase(N, EA);
            std::printf("  %6d %12.4e %12.4e %10ld %8s\n", N, e.l2, e.tensionErr, e.rr.steps,
                        e.rr.converged ? "yes" : "NO");
            CHECK(e.rr.converged);
            l2.push_back(e.l2);
        }
        for (size_t i = 1; i < l2.size(); ++i) CHECK(l2[i] < l2[i - 1]);       // monotone convergence
        const double order = std::log(l2[1] / l2[4]) / std::log(8.0);           // N = 20 -> 160
        std::printf("  observed order (N=20..160): %.2f\n", order);
        CHECK(order > 1.5);
        CHECK(l2.back() < 1e-4);
    }
}
