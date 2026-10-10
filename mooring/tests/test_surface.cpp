// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Submergence at the instantaneous free surface (Environment::elevation).
#include <cmath>
#include <vector>
#include "doctest.h"
#include "mooring/lumped_mass_cable.hpp"

using namespace mooring;

namespace {
LumpedMassCable makeLine(const Environment& env) {
    CableParams p; p.L = 10.0; p.N = 20; p.EA = 1e5; p.m_l = 0.5; p.w = 3.0; p.D0 = 0.02; p.Cdn = 1.0; p.Cdt = 0.5;
    LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(8.0, 0, 5.0));
    c.setEnvironment(env);
    c.addPointElement(PointElement::floater(10, 0.2, 5.0, 0.1, 0.5));
    return c;
}
}  // namespace

TEST_CASE("Instantaneous surface: elevation eta shifts the blend exactly like surfaceZ + eta (weights, buoyancy, drag)") {
    const double A = 0.4, om = 2.0;
    Environment still; still.hydro = true; still.surfaceZ = 2.5;
    Environment wavy = still;
    wavy.elevation = [&](const Vec3& x, double t) { return A * std::sin(om * t - 0.3 * x.x); };
    for (double t : {0.0, 0.7, 1.9, 3.3}) {
        const LumpedMassCable cw = makeLine(wavy);
        for (int i = 0; i <= 20; i += 3) {
            Environment shifted = still; shifted.surfaceZ = 2.5 + wavy.elevation(cw.nodes()[i], t);
            const LumpedMassCable cs = makeLine(shifted);
            CHECK(cw.submergedFraction(cw.nodes(), i, t) == doctest::Approx(cs.submergedFraction(cs.nodes(), i, t)).scale(1));
            CHECK(cw.nodeWeight(cw.nodes(), i, t) == doctest::Approx(cs.nodeWeight(cs.nodes(), i, t)));
        }
        const LumpedMassCable c10 = makeLine(wavy);
        Environment sh = still; sh.surfaceZ = 2.5 + wavy.elevation(c10.nodes()[10], t);
        CHECK(c10.pointNetWeight(c10.nodes(), 10, t) == doctest::Approx(makeLine(sh).pointNetWeight(c10.nodes(), 10, t)));
    }
}

TEST_CASE("Instantaneous surface: the floater's buoyancy tracks the crest and trough; static solve ignores the waves") {
    Environment env; env.hydro = true; env.surfaceZ = 2.5;
    env.elevation = [](const Vec3&, double t) { return 0.4 * std::sin(2.0 * t); };
    LumpedMassCable c = makeLine(env);
    const double zf = c.nodes()[10].z;
    // The floater node is at the interpolated height; at crest/trough the submerged fraction differs and the buoyant share follows it.
    const double tc = std::acos(0.0) / 2.0 * 1.0, tt = 3.0 * std::acos(0.0) / 2.0;   // sin(2 t) = +1 / -1
    const double phiC = c.submergedFraction(c.nodes(), 10, tc), phiT = c.submergedFraction(c.nodes(), 10, tt);
    const double bc = 0.2 * 9.81 - c.pointNetWeight(c.nodes(), 10, tc), bt = 0.2 * 9.81 - c.pointNetWeight(c.nodes(), 10, tt);
    CHECK(bc == doctest::Approx(phiC * 5.0).epsilon(1e-9));
    CHECK(bt == doctest::Approx(phiT * 5.0).epsilon(1e-9));
    CHECK(phiC >= phiT);
    (void)zf;
    // Static relaxation uses the still-water level: same shape with and without the elevation function.
    Environment e0 = env; e0.elevation = nullptr;
    LumpedMassCable a = makeLine(env), b = makeLine(e0);
    a.relaxStatic(); b.relaxStatic();
    for (int i = 0; i <= 20; ++i) CHECK(norm(a.nodes()[i] - b.nodes()[i]) == doctest::Approx(0.0).scale(1).epsilon(1e-12));
}
