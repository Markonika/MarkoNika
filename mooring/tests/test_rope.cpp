// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Nonlinear / viscoelastic rope law (include/mooring/rope.hpp): unit tests against closed forms and a cable-level dynamic-stiffness test.
#include <cmath>
#include <cstdio>
#include "doctest.h"
#include "mooring/lumped_mass_cable.hpp"

using namespace mooring;

namespace {
const double kPi = 3.14159265358979323846;
}

TEST_CASE("Rope law: piecewise-linear static curve, slope, energy, extrapolation") {
    RopeLaw r; r.curve = {{0.02, 200.0}, {0.05, 800.0}};   // slopes 1e4 then 2e4 N
    const double EA = 1e4;
    CHECK(r.staticTension(EA, -0.1) == doctest::Approx(0.0).scale(1));
    CHECK(r.staticTension(EA, 0.01) == doctest::Approx(100.0));
    CHECK(r.staticTension(EA, 0.02) == doctest::Approx(200.0));
    CHECK(r.staticTension(EA, 0.035) == doctest::Approx(500.0));
    CHECK(r.staticTension(EA, 0.08) == doctest::Approx(1400.0));        // last slope extrapolated
    CHECK(r.staticSlope(EA, 0.01) == doctest::Approx(1e4));
    CHECK(r.staticSlope(EA, 0.04) == doctest::Approx(2e4));
    CHECK(r.maxStaticStiffness(EA) == doctest::Approx(2e4));
    // Energy = integral of the tension (trapezoid rule is exact for piecewise-linear T)
    CHECK(r.staticEnergy(EA, 0.02) == doctest::Approx(0.5 * 200.0 * 0.02));
    CHECK(r.staticEnergy(EA, 0.05) == doctest::Approx(0.5 * 200.0 * 0.02 + 0.5 * (200.0 + 800.0) * 0.03));
    CHECK(r.staticEnergy(EA, 0.08) == doctest::Approx(0.5 * 200.0 * 0.02 + 0.5 * (200.0 + 800.0) * 0.03 + 0.5 * (800.0 + 1400.0) * 0.03));
    RopeLaw bad; bad.curve = {{0.05, 800.0}, {0.02, 900.0}};
    CHECK_THROWS(bad.validate());
}

TEST_CASE("Rope law: Maxwell branch matches stress relaxation and the complex modulus") {
    // Ramp to strain e1 at rate e1/t1 then hold: alpha(t) closed form for t > t1 is e1 - (e1 tau/t1)(1-exp(-t1/tau)) exp(-(t-t1)/tau).
    const double K = 3e4, tau = 0.4, e1 = 0.01, t1 = 0.2, dt = 1e-3;
    double a = 0.0, eps = 0.0;
    for (double t = 0.0; t < 1.0 - 1e-12; t += dt) {
        const double e0 = eps; eps = t + dt <= t1 ? e1 * (t + dt) / t1 : e1;
        a = RopeLaw::advanceAlpha(a, e0, eps, dt, tau);
    }
    const double ex = e1 - (e1 * tau / t1) * (1.0 - std::exp(-t1 / tau)) * std::exp(-(1.0 - t1) / tau);
    CHECK(a == doctest::Approx(ex).epsilon(1e-9));
    // Sinusoidal strain: storage K w^2 tau^2/(1+w^2 tau^2), loss K w tau/(1+w^2 tau^2) from a Fourier fit over whole cycles after a transient.
    for (double w : {0.5, 2.5, 12.0}) {
        const double ea = 1e-3, T = 2.0 * kPi / w;
        double al = 0.0, e = 0.0, s = 0.0, c = 0.0; const int nc = 40; const double h = T / 2000.0;
        for (int k = 0; k < nc * 2000; ++k) {
            const double t = k * h, e0 = e; e = ea * std::sin(w * (t + h));
            al = RopeLaw::advanceAlpha(al, e0, e, h, tau);
            if (k >= 20 * 2000) { const double F = K * (e - al); s += F * std::sin(w * (t + h)) * h; c += F * std::cos(w * (t + h)) * h; }
        }
        const double span = 20.0 * T, Es = 2.0 * s / span / ea, Ec = 2.0 * c / span / ea;
        const double wt = w * tau;
        std::printf("\n  w=%.1f: E' %.1f (exact %.1f), E'' %.1f (exact %.1f)\n", w, Es, K * wt * wt / (1 + wt * wt), Ec, K * wt / (1 + wt * wt));
        CHECK(Es == doctest::Approx(K * wt * wt / (1 + wt * wt)).epsilon(2e-3));
        CHECK(Ec == doctest::Approx(K * wt / (1 + wt * wt)).epsilon(2e-3));
    }
}

TEST_CASE("Rope law in the cable: static tension follows the nonlinear curve; relaxed state carries no branch force") {
    CableParams p; p.L = 10.0; p.N = 10; p.EA = 0.0; p.m_l = 0.1; p.w = 0.0; p.g = 0.0;
    p.rope.curve = {{0.02, 200.0}, {0.05, 800.0}};
    p.rope.branches = {{5e3, 0.5}};
    for (double eps : {0.01, 0.03, 0.07}) {
        LumpedMassCable cable(p, Vec3(0, 0, 0), Vec3(10.0 * (1 + eps), 0, 0));
        std::vector<Vec3> r(11), v(11);
        for (int i = 0; i <= 10; ++i) r[i] = Vec3(i * (1 + eps), 0, 0);
        cable.setInitialState(r, v);
        RelaxOptions ro; ro.forceTol = 1e-6;   // weightless taut line: the default tolerance is below round-off
        const RelaxResult rr = cable.relaxStatic(ro);
        CHECK(rr.converged);
        const double expect = p.rope.staticTension(1e4, eps);
        CHECK(norm(cable.endTension(true)) == doctest::Approx(expect).epsilon(1e-6));
        CHECK(cable.endSegmentTension(true) == doctest::Approx(expect).epsilon(1e-6));
    }
}

TEST_CASE("Rope law in the cable: dynamic stiffness of a Maxwell branch under slow sinusoidal fairlead motion") {
    CableParams p; p.L = 10.0; p.N = 10; p.EA = 1e5; p.m_l = 0.1; p.w = 0.0; p.g = 0.0;
    p.rope.branches = {{5e4, 0.5}};
    const double e0 = 0.02, L = p.L;
    LumpedMassCable cable(p, Vec3(0, 0, 0), Vec3(L * (1 + e0), 0, 0));
    std::vector<Vec3> r(11), v(11);
    for (int i = 0; i <= 10; ++i) r[i] = Vec3(i * (1 + e0), 0, 0);
    cable.setInitialState(r, v);
    cable.relaxStatic();
    const double w = 2.0 * kPi / 2.0, a = 0.02, ea = a / L, ramp = 3.0 * 2.0;   // 2 s period, 3-cycle cosine ramp
    cable.setTopMotion([&](double t, Vec3& pos, Vec3& vel) {
        const double s = t < ramp ? 0.5 * (1.0 - std::cos(kPi * t / ramp)) : 1.0;
        const double ds = t < ramp ? 0.5 * kPi / ramp * std::sin(kPi * t / ramp) : 0.0;
        pos = Vec3(L * (1 + e0) + a * s * std::sin(w * t), 0, 0);
        vel = Vec3(a * (ds * std::sin(w * t) + s * w * std::cos(w * t)), 0, 0);
    });
    const double dt = cable.stableDt(), tFit0 = ramp + 2.0 * 2.0, tEnd = tFit0 + 4.0 * 2.0;
    double sS = 0, sC = 0, n = 0;
    while (cable.time() < tEnd) {
        cable.advanceTo(cable.time() + dt);
        if (cable.time() >= tFit0) { const double T = cable.endSegmentTension(true); sS += T * std::sin(w * cable.time()); sC += T * std::cos(w * cable.time()); ++n; }
    }
    const double Ts = 2.0 * sS / n, Tc = 2.0 * sC / n, wt = w * 0.5;
    const double Es = 1e5 + 5e4 * wt * wt / (1 + wt * wt), Ec = 5e4 * wt / (1 + wt * wt);
    std::printf("\n  dynamic stiffness: E' %.0f (exact %.0f), E'' %.0f (exact %.0f), dt %.2e s\n", Ts / ea, Es, Tc / ea, Ec, dt);
    CHECK(Ts / ea == doctest::Approx(Es).epsilon(0.02));
    CHECK(Tc / ea == doctest::Approx(Ec).epsilon(0.03));
}
