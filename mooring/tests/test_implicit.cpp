// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Implicit treatment of the internal axial damping (DynOptions::implicitDamping): longitudinal mode of a taut chain with stiffness-proportional damping.
#include <cmath>
#include <cstdio>
#include <vector>
#include "doctest.h"
#include "mooring/lumped_mass_cable.hpp"

using namespace mooring;

namespace {
const double kPi = 3.14159265358979323846;

// Fixed-fixed taut chain, first longitudinal mode as initial displacement. Returns the mode energy E - E0 at the requested times.
struct ModeRun { std::vector<double> E; double dt; double omega; bool finite; };
ModeRun run(double beta, bool implicitDamping, const std::vector<double>& times, int N = 40) {
    const double L = 10.0, EA = 1e5, ml = 0.1, e0 = 1e-3, u0 = 1e-5;
    CableParams p; p.L = L; p.N = N; p.EA = EA; p.m_l = ml; p.w = 0.0; p.g = 0.0; p.c_int = beta * EA;
    const double l0 = p.l0();
    LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(L * (1 + e0), 0, 0));
    std::vector<Vec3> r(N + 1), v(N + 1);
    for (int i = 0; i <= N; ++i) r[i] = Vec3(i * l0 * (1 + e0) + u0 * std::sin(kPi * i * l0 / L), 0, 0);
    c.setInitialState(r, v);
    DynOptions o; o.implicitDamping = implicitDamping; c.setDynOptions(o);
    const double E0 = N * 0.5 * EA * l0 * e0 * e0;
    ModeRun out; out.dt = c.stableDt(); out.finite = true;
    out.omega = 2.0 * std::sqrt(EA / ml) / l0 * std::sin(kPi * l0 / (2.0 * L));
    for (double t : times) {
        c.advanceTo(t);
        const double E = c.energy() - E0;
        if (!std::isfinite(E)) out.finite = false;
        out.E.push_back(E);
    }
    return out;
}
}  // namespace

TEST_CASE("Implicit damping: mode energy decays at the exact stiffness-proportional rate beta*omega^2, like the explicit scheme") {
    const double beta = 1e-4;
    const std::vector<double> ts = {0.04, 0.08, 0.12};
    const ModeRun ex = run(beta, false, ts), im = run(beta, true, ts), ref0 = run(0.0, false, {1e-9});
    const double E0 = ref0.E[0];
    for (size_t k = 0; k < ts.size(); ++k) {
        const double exact = E0 * std::exp(-beta * im.omega * im.omega * ts[k]);
        std::printf("\n  t=%.2f: E/E0 explicit %.5f, implicit %.5f, exact %.5f\n", ts[k], ex.E[k] / E0, im.E[k] / E0, exact / E0);
        CHECK(ex.E[k] == doctest::Approx(exact).epsilon(0.03));
        CHECK(im.E[k] == doctest::Approx(exact).epsilon(0.03));
        CHECK(im.E[k] == doctest::Approx(ex.E[k]).epsilon(0.01));
    }
    CHECK(im.dt >= ex.dt);
}

TEST_CASE("Implicit damping: stable far beyond the explicit damping limit and monotonically dissipative") {
    const double beta = 5e-3;                                     // explicit step limit ~40 times below the wave CFL step
    const std::vector<double> ts = {0.01, 0.02, 0.05, 0.1, 0.2, 0.4};
    const ModeRun ex = run(beta, false, {1e-9}), im = run(beta, true, ts);
    std::printf("\n  explicit dt %.2e s, implicit dt %.2e s\n", ex.dt, im.dt);
    CHECK(im.dt > 20.0 * ex.dt);
    CHECK(im.finite);
    // The static energy is ~0.5 J, so differences below ~1e-12 J are round-off.
    double prev = 1e300;
    for (double E : im.E) { CHECK(E <= prev + 1e-12); CHECK(E >= -1e-12); prev = E; }
    CHECK(im.E.back() < 1e-3 * im.E.front() + 1e-12);              // mode 1 is strongly damped (zeta 0.78) after 0.4 s
}
