// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Linearisation: Jacobi eigen solver, natural modes of a line, linearised mooring stiffness and natural frequencies of the body.
#include <cmath>
#include <cstdio>
#include <vector>
#include "doctest.h"
#include "mooring/linearize.hpp"
#include "mooring/lumped_mass_cable.hpp"
#include "mooring/rigid_body.hpp"

using namespace mooring;

namespace {
const double kPi = 3.14159265358979323846;
}

TEST_CASE("Linearize: Jacobi eigen-decomposition (1-D Laplacian eigenvalues, orthonormal eigenvectors)") {
    const int n = 30;
    std::vector<double> a(n * n, 0.0);
    for (int i = 0; i < n; ++i) { a[i * n + i] = 2.0; if (i + 1 < n) a[i * n + i + 1] = a[(i + 1) * n + i] = -1.0; }
    const SymEigen e = symmetricEigen(a, n);
    for (int k = 0; k < n; ++k) CHECK(e.w[k] == doctest::Approx(2.0 - 2.0 * std::cos((k + 1) * kPi / (n + 1))).epsilon(1e-12));
    double maxdev = 0;
    for (int p = 0; p < n; ++p) for (int q = 0; q < n; ++q) {
        double s = 0; for (int i = 0; i < n; ++i) s += e.v[p * n + i] * e.v[q * n + i];
        maxdev = std::max(maxdev, std::fabs(s - (p == q ? 1.0 : 0.0)));
    }
    CHECK(maxdev < 1e-12);
    // A v = w v
    double res = 0;
    for (int k = 0; k < n; ++k) for (int i = 0; i < n; ++i) { double s = 0; for (int j = 0; j < n; ++j) s += a[i * n + j] * e.v[k * n + j]; res = std::max(res, std::fabs(s - e.w[k] * e.v[k * n + i])); }
    CHECK(res < 1e-12);
}

TEST_CASE("Linearize: natural modes of a taut chain match the discrete transverse and longitudinal formulas; modal damping = beta omega / 2") {
    const double L = 10.0, EA = 1e5, ml = 0.1, e0 = 1e-3, beta = 1e-4; const int N = 40;
    CableParams p; p.L = L; p.N = N; p.EA = EA; p.m_l = ml; p.w = 0.0; p.g = 0.0; p.planar = true; p.c_int = beta * EA;
    const double l0 = p.l0(), T = EA * e0;
    LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(L * (1 + e0), 0, 0));
    std::vector<Vec3> r(N + 1), v(N + 1);
    for (int i = 0; i <= N; ++i) r[i] = Vec3(i * l0 * (1 + e0), 0, 0);
    c.setInitialState(r, v);
    const LineModes m = c.naturalModes(2 * (N - 1));
    REQUIRE(int(m.omega.size()) == 2 * (N - 1));
    CHECK(m.dof == 2 * (N - 1));
    // transverse (z) modes first: omega_k = (2/l0) sqrt(T / (ml (1 + e0))) sin(k pi / (2N)); longitudinal: (2/l0) sqrt(EA/ml) sin(k pi/(2N))
    for (int k = 1; k <= 3; ++k) {
        const double wt = 2.0 / l0 * std::sqrt(T / (ml * (1 + e0))) * std::sin(k * kPi / (2 * N));
        std::printf("\n  transverse mode %d: omega %.5f (analytic %.5f)", k, m.omega[k - 1], wt);
        CHECK(m.omega[k - 1] == doctest::Approx(wt).epsilon(1e-4));
    }
    const double wl1 = 2.0 / l0 * std::sqrt(EA / ml) * std::sin(kPi / (2 * N));
    int kl = -1;
    for (int k = 0; k < int(m.omega.size()); ++k) if (std::fabs(m.omega[k] - wl1) < 0.05 * wl1) { kl = k; break; }
    REQUIRE(kl >= 0);
    std::printf("\n  longitudinal mode 1: omega %.4f (analytic %.4f), zeta %.5f (beta omega / 2 = %.5f)\n", m.omega[kl], wl1, m.zeta[kl], 0.5 * beta * wl1);
    CHECK(m.omega[kl] == doctest::Approx(wl1).epsilon(2e-4));
    CHECK(m.zeta[kl] == doctest::Approx(0.5 * beta * wl1).epsilon(2e-2));
    CHECK(m.asymmetry < 1e-3);
}

TEST_CASE("Linearize: two-section string, lowest natural frequency vs the exact two-medium eigenvalue equation") {
    CableParams p; p.g = 0.0; p.planar = true;
    LineSection a, b;
    a.length = 4.0; a.segments = 80;  a.EA = 1e5; a.m_l = 0.1; a.w = 0.0;
    b.length = 6.0; b.segments = 120; b.EA = 3e5; b.m_l = 0.3; b.w = 0.0;
    p.sections = {a, b};
    const double T = 200.0;
    const double e1 = T / a.EA, e2 = T / b.EA, S1 = a.length * (1 + e1), S2 = b.length * (1 + e2);
    const double rho1 = a.m_l / (1 + e1), rho2 = b.m_l / (1 + e2);
    auto F = [&](double w) { const double k1 = w * std::sqrt(rho1 / T), k2 = w * std::sqrt(rho2 / T); return k1 / std::tan(k1 * S1) + k2 / std::tan(k2 * S2); };
    double lo = 0, hi = 0;
    for (double w = 0.5, st = 0.01; w < 200.0; w += st)
        if (F(w) * F(w + st) < 0.0 && std::fabs(F(w)) < 50.0 && std::fabs(F(w + st)) < 50.0) { lo = w; hi = w + st; break; }
    REQUIRE(hi > 0.0);
    for (int i = 0; i < 80; ++i) { const double mid = 0.5 * (lo + hi); (F(lo) * F(mid) <= 0.0 ? hi : lo) = mid; }
    const double wExact = 0.5 * (lo + hi);
    LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(S1 + S2, 0, 0));
    const int N = 200;
    std::vector<Vec3> r(N + 1), v(N + 1);
    double x = 0.0;
    for (int j = 0; j <= N; ++j) { r[j] = Vec3(x, 0, 0); if (j < N) x += (j < 80 ? a.length / 80 * (1 + e1) : b.length / 120 * (1 + e2)); }
    c.setInitialState(r, v);
    const LineModes m = c.naturalModes(3);
    std::printf("\n  two-medium string: omega exact %.4f, linearised %.4f (%+.4f %%)\n", wExact, m.omega[0], 100 * (m.omega[0] - wExact) / wExact);
    CHECK(m.omega[0] == doctest::Approx(wExact).epsilon(3e-3));
}

namespace {
struct Hang { std::unique_ptr<CoupledSystem> sys; LumpedMassCable* cable{nullptr}; double M, L0, EA; };
Hang makeHang() {
    Hang h; h.M = 2.0; h.L0 = 10.0; h.EA = 2000.0;
    const double ml = 0.001; const int N = 20;
    BodyParams bp; bp.mass = h.M; bp.cgRef = Vec3(0, 0, -h.L0 * 1.005);
    bp.F0[2] = -h.M * 9.81;
    bp.inertia = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    h.sys = std::make_unique<CoupledSystem>(bp);
    CableParams p; p.L = h.L0; p.N = N; p.EA = h.EA; p.m_l = ml; p.w = ml * 9.81;
    auto cab = std::make_unique<LumpedMassCable>(p, Vec3(0, 0, 0), Vec3(0, 0, -h.L0 * 1.005));
    h.cable = cab.get();
    h.sys->addLine(std::move(cab), Vec3(0, 0, 0), "hang");
    REQUIRE(h.sys->solveEquilibrium(1e-8).converged);
    return h;
}
}  // namespace

TEST_CASE("Linearize: mooring stiffness of a body hanging on a line (axial EA/L, geometric T/L) and the resulting natural frequencies") {
    Hang h = makeHang();
    const Mat6 K = h.sys->stiffness();
    const double Ltot = norm(h.sys->fairleadPosition(0));
    const double T = h.cable->endSegmentTension(true);
    const double Kz = h.EA / h.L0, Kx = T / Ltot;
    std::printf("\n  K(heave) %.4f (EA/L0 = %.4f), K(surge) %.5f (T/L = %.5f), T = %.4f N, asymmetry:", K(2, 2), Kz, K(0, 0), Kx, T);
    CHECK(K(2, 2) == doctest::Approx(Kz).scale(0).epsilon(1e-4));
    CHECK(K(0, 0) == doctest::Approx(Kx).scale(0).epsilon(5e-3));     // T/L ignores the weight of the line itself (0.24 %)
    CHECK(K(1, 1) == doctest::Approx(Kx).scale(0).epsilon(5e-3));
    const BodyModes bm = bodyModes(h.sys->body().totalMass(), K);
    double asym = 0; for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) asym = std::max(asym, std::fabs(K(i, j) - K(j, i)));
    std::printf(" max |K - K^T| %.2e\n", asym);
    // heave: omega = sqrt(K/M); pendulum: omega = sqrt(g/L) with the lateral (geometric) stiffness T/L
    bool foundHeave = false, foundPend = false;
    for (int k = 0; k < 6; ++k) {
        if (std::fabs(bm.omega[k] - std::sqrt(Kz / h.M)) < 1e-3 * std::sqrt(Kz / h.M)) foundHeave = true;
        if (std::fabs(bm.omega[k] - std::sqrt(Kx / h.M)) < 5e-3 * std::sqrt(Kx / h.M)) foundPend = true;
    }
    CHECK(foundHeave);
    CHECK(foundPend);
    CHECK(h.sys->lineStiffness()(2, 2) == doctest::Approx(K(2, 2)).scale(0).epsilon(1e-9));       // no hydrostatic C in this set-up
}
