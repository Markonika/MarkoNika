// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Bending stiffness: discrete curvature energy U = sum EI |t_i - t_{i-1}|^2 / (2 lbar).
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>
#include "doctest.h"
#include "mooring/lumped_mass_cable.hpp"

using namespace mooring;

namespace {
const double kPi = 3.14159265358979323846;

// Weightless-axially line between fixed ends with bending stiffness EI, tiny axial stiffness (negligible tension), gravity as the lateral load.
LumpedMassCable beamLine(double L, int N, double EI, double w, bool clamped, double EA = 1.0) {
    CableParams p; p.L = L; p.N = N; p.EA = EA; p.m_l = 0.1; p.g = 9.81; p.w = w; p.EI = EI; p.planar = true;
    if (clamped) { p.clampAnchor = p.clampFairlead = true; p.clampAnchorDir = p.clampFairleadDir = Vec3(1, 0, 0); }
    LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(L, 0, 0));
    std::vector<Vec3> r(N + 1), v(N + 1);
    for (int i = 0; i <= N; ++i) r[i] = Vec3(L * i / N, 0, 0);
    c.setInitialState(r, v);
    return c;
}
}  // namespace

TEST_CASE("Bending: the force is minus the gradient of the bending energy (random 3-D shape, clamped ends)") {
    const int N = 12;
    CableParams p; p.L = 6.0; p.N = N; p.EA = 1e5; p.m_l = 0.1; p.w = 0.0; p.g = 0.0; p.EI = 3.5;
    p.clampAnchor = p.clampFairlead = true; p.clampAnchorDir = Vec3(0.3, 0.2, -1.0); p.clampFairleadDir = Vec3(1, 0.4, 0.2);
    LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(5.5, 0, 0));
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> u(-0.12, 0.12);
    std::vector<Vec3> r(N + 1), v(N + 1);
    for (int i = 0; i <= N; ++i) r[i] = Vec3(5.5 * i / N + u(rng), u(rng), u(rng));
    r[0] = Vec3(0, 0, 0); r[N] = Vec3(5.5, 0, 0);
    c.setInitialState(r, v);
    CableParams p0 = p; p0.EI = 0.0; p0.clampAnchor = p0.clampFairlead = false;
    LumpedMassCable c0(p0, Vec3(0, 0, 0), Vec3(5.5, 0, 0));
    c0.setInitialState(r, v);
    std::vector<Vec3> f, f0;
    c.computeForces(r, v, 0.0, f);
    c0.computeForces(r, v, 0.0, f0);
    // bending force = difference; compare with -dU/dr by central differences of energy() differences
    double worst = 0.0, fmax = 0.0;
    for (int i = 0; i <= N; ++i) {
        const Vec3 fb = f[i] - f0[i];
        for (int d = 0; d < 3; ++d) {
            const double h = 1e-6;
            auto Ub = [&](double shift) {
                std::vector<Vec3> rp = r;
                (d == 0 ? rp[i].x : d == 1 ? rp[i].y : rp[i].z) += shift;
                c.setInitialState(rp, v); c0.setInitialState(rp, v);
                return c.energy() - c0.energy();
            };
            const double g = -(Ub(h) - Ub(-h)) / (2 * h);
            const double fv = d == 0 ? fb.x : d == 1 ? fb.y : fb.z;
            worst = std::max(worst, std::fabs(g - fv)); fmax = std::max(fmax, std::fabs(fv));
        }
    }
    std::printf("\n  max |F - (-dU/dr)| = %.2e (largest bending force component %.3f N)\n", worst, fmax);
    CHECK(fmax > 0.1);
    CHECK(worst < 1e-5 * fmax + 1e-8);
    // the bending forces on the whole line sum to zero (translation invariance) and have no net moment about the origin
    Vec3 sum, mom;
    for (int i = 0; i <= N; ++i) { const Vec3 fb = f[i] - f0[i]; sum += fb; mom += Vec3(r[i].y * fb.z - r[i].z * fb.y, r[i].z * fb.x - r[i].x * fb.z, r[i].x * fb.y - r[i].y * fb.x); }
    // (the clamp directions break rotation invariance, so only the net force is tested)
    CHECK(norm(sum) < 1e-10 * fmax + 1e-10);
    (void)mom;
}

TEST_CASE("Bending: beam statics - pinned (5qL^4/384EI) and clamped (qL^4/384EI) mid-span deflection under uniform load") {
    const double L = 10.0, EI = 100.0, q = 0.003;       // small load: the membrane tension (EA times the arc-length gain) and the axial softness stay negligible
    const int N = 40;
    for (bool clamped : {false, true}) {
        // w acts in -z; EA = 1e4 N keeps the axial extension negligible, the small load keeps the membrane tension (EA x arc-length gain) negligible. Planar x-z.
        LumpedMassCable c = beamLine(L, N, EI, q, clamped, 1e4);
        RelaxOptions ro; ro.forceTol = 1e-6; ro.maxSteps = 5000000;
        const RelaxResult rr = c.relaxStatic(ro);
        REQUIRE(rr.converged);
        const double mid = -c.nodes()[N / 2].z;
        const double exact = clamped ? q * std::pow(L, 4) / (384.0 * EI) : 5.0 * q * std::pow(L, 4) / (384.0 * EI);
        std::printf("\n  %s beam: mid-span deflection %.5f m (beam theory %.5f m, %+.2f %%)", clamped ? "clamped" : "pinned", mid, exact, 100 * (mid - exact) / exact);
        CHECK(mid == doctest::Approx(exact).scale(0).epsilon(0.01));      // (the nodal load omits the half-segment loads at the supports: -2.5 %/-5 % of the load)
    }
    std::printf("\n");
}

TEST_CASE("Bending: natural frequencies of a taut pinned beam, omega_k^2 = (T (k pi/L)^2 + EI (k pi/L)^4) / rho (finite-difference Jacobian, five-node stencil)") {
    const double L = 10.0, EA = 1e5, ml = 0.1, EI = 50.0, e0 = 1e-3; const int N = 100;
    CableParams p; p.L = L; p.N = N; p.EA = EA; p.m_l = ml; p.w = 0.0; p.g = 0.0; p.EI = EI; p.planar = true;
    const double T = EA * e0, rho = ml / (1.0 + e0), Ls = L * (1.0 + e0);
    LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(Ls, 0, 0));
    std::vector<Vec3> r(N + 1), v(N + 1);
    for (int i = 0; i <= N; ++i) r[i] = Vec3(Ls * i / N, 0, 0);
    c.setInitialState(r, v);
    const LineModes m = c.naturalModes(4);
    for (int k = 1; k <= 3; ++k) {
        const double kk = k * kPi / Ls, ex = std::sqrt((T * kk * kk + EI * kk * kk * kk * kk) / rho);
        std::printf("\n  mode %d: omega %.4f (beam-string theory %.4f, %+.3f %%)", k, m.omega[k - 1], ex, 100 * (m.omega[k - 1] - ex) / ex);
        CHECK(m.omega[k - 1] == doctest::Approx(ex).scale(0).epsilon(2e-3));
    }
    // without bending the same line gives the string frequencies (the bending terms are what changed)
    CableParams p0 = p; p0.EI = 0.0;
    LumpedMassCable c0(p0, Vec3(0, 0, 0), Vec3(Ls, 0, 0));
    c0.setInitialState(r, v);
    const LineModes m0 = c0.naturalModes(3);
    CHECK(m.omega[2] > 1.15 * m0.omega[2]);
    std::printf("\n  mode 3: with EI %.4f rad/s, without %.4f rad/s\n", m.omega[2], m0.omega[2]);
}

TEST_CASE("Bending: the implicit scheme (pentadiagonal Jacobian) and the explicit scheme agree on a vibrating beam; explicit step respects the bending limit") {
    const double L = 10.0, EA = 1e5, ml = 0.1, EI = 50.0, e0 = 1e-3; const int N = 60;
    CableParams p; p.L = L; p.N = N; p.EA = EA; p.m_l = ml; p.w = 0.0; p.g = 0.0; p.EI = EI; p.planar = true;
    const double Ls = L * (1.0 + e0);
    auto run = [&](Scheme sch, double factor, double* dtOut) {
        LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(Ls, 0, 0));
        std::vector<Vec3> r(N + 1), v(N + 1);
        for (int i = 0; i <= N; ++i) r[i] = Vec3(Ls * i / N, 0, 1e-3 * std::sin(kPi * i / N) + 4e-4 * std::sin(3 * kPi * i / N));
        c.setInitialState(r, v);
        DynOptions o; o.scheme = sch; o.implicitDtFactor = factor; o.rhoInf = 1.0; c.setDynOptions(o);
        *dtOut = c.stableDt();
        c.advanceTo(0.4);
        return c.nodes();
    };
    double dtE, dtI;
    const auto re = run(Scheme::RK4, 1.0, &dtE);
    const auto ri = run(Scheme::GenAlpha, 5.0, &dtI);
    double d = 0, nrm = 0;
    for (int i = 0; i <= N; ++i) { d = std::max(d, std::fabs(re[i].z - ri[i].z)); nrm = std::max(nrm, std::fabs(re[i].z)); }
    std::printf("\n  explicit dt %.2e s (bending limit), implicit dt %.2e s: max difference %.2e of amplitude %.2e\n", dtE, dtI, d, nrm);
    CHECK(dtE < 0.5 * 0.5 * (L / N) * (L / N) * std::sqrt(ml / EI) * 1.0001);           // the bending limit is the active one
    CHECK(d < 0.02 * nrm);
}
