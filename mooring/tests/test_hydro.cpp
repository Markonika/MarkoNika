// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Milestone 3: Morison hydrodynamics (Eqs. 3.26-3.31) and seabed (Eqs. 3.32-3.36).
#include <cmath>
#include <cstdio>
#include <vector>
#include "doctest.h"
#include "mooring/hydro.hpp"
#include "mooring/lumped_mass_cable.hpp"

using namespace mooring;

namespace { const double kPi = 3.14159265358979323846; }

TEST_CASE("added-mass node solve: M a = f for oblique tangent") {
    const Vec3 t = Vec3(1, 2, 2) / 3.0;                 // unit vector
    const Vec3 f(0.3, -1.1, 0.7);
    const double m = 0.8, c = 0.35;
    const Vec3 a = addedMassSolve(m, c, t, f);
    const Vec3 Ma = a * m + (a - t * dot(a, t)) * c;     // M a = m a + c (I - t t^T) a
    CHECK(norm(Ma - f) < 1e-14);
    // Pure tangential force sees only m; pure normal force sees m + c.
    CHECK(addedMassSolve(m, c, Vec3(0, 0, 1), Vec3(0, 0, 2)).z == doctest::Approx(2.0 / m));
    CHECK(addedMassSolve(m, c, Vec3(0, 0, 1), Vec3(2, 0, 0)).x == doctest::Approx(2.0 / (m + c)));
}

TEST_CASE("Morison drag: tangential and normal closed forms") {
    const double rho = 1000, D = 0.0022, len = 0.5;
    const Vec3 t(1, 0, 0);
    Vec3 F = morisonDrag(Vec3(-2, 0, 0), t, 0.5, 2.5, rho, D, len);   // cable moves +x at 2 m/s relative to water
    CHECK(F.x == doctest::Approx(-0.5 * 0.5 * rho * D * 4.0 * len));
    CHECK(F.y == 0.0);
    F = morisonDrag(Vec3(0, 3, 4), t, 0.5, 2.5, rho, D, len);         // |vn| = 5
    CHECK(F.y == doctest::Approx(0.5 * 2.5 * rho * D * 5.0 * 3.0 * len));
    CHECK(F.z == doctest::Approx(0.5 * 2.5 * rho * D * 5.0 * 4.0 * len));
    CHECK(F.x == 0.0);
}

TEST_CASE("seabed force law: spring, one-sided damping, ramped Coulomb friction") {
    SoilParams s; s.Ks = 3e9; s.zeta = 1.0; s.mu = 0.3; s.vlim = 0.01;
    const double D1 = 0.0022, ml = 0.0818, wSub = 0.699, len = 0.33;
    CHECK(norm(seabedForce(s, D1, ml, wSub, -1e-3, Vec3(1, 0, -1), len)) == 0.0);        // no contact
    Vec3 F = seabedForce(s, D1, ml, wSub, 1e-4, Vec3(), len);                             // static spring
    CHECK(F.z == doctest::Approx(3e9 * D1 * 1e-4 * len));
    CHECK(F.x == 0.0);
    const double cd = 2 * 1.0 * std::sqrt(3e9 * D1 * ml);
    F = seabedForce(s, D1, ml, wSub, 1e-4, Vec3(0, 0, -0.1), len);                        // penetrating: damped
    CHECK(F.z == doctest::Approx((3e9 * D1 * 1e-4 + cd * 0.1) * len));
    F = seabedForce(s, D1, ml, wSub, 1e-4, Vec3(0, 0, +0.1), len);                        // lifting: no damping
    CHECK(F.z == doctest::Approx(3e9 * D1 * 1e-4 * len));
    F = seabedForce(s, D1, ml, wSub, 1e-4, Vec3(0.005, 0, 0), len);                       // half of v_lim
    CHECK(F.x == doctest::Approx(-wSub * 0.3 * 0.5 * len));
    F = seabedForce(s, D1, ml, wSub, 1e-4, Vec3(0.5, 0.5, 0), len);                       // above v_lim, oblique
    CHECK(F.x == doctest::Approx(-wSub * 0.3 * len / std::sqrt(2.0)));
    CHECK(F.y == doctest::Approx(F.x));
}

namespace {
// Taut string, L0 = 1 m, eps = 0.01, EA = 100 N (T = 1 N), m_l = 0.05 kg/m, submerged, no gravity.
struct StringRun { double omegaNum, ampEnd; };

// Modal projection of the y-displacement on sin(pi s/L0) -> (q, qdot).
void modal(const LumpedMassCable& c, double& q, double& qd) {
    const int N = c.params().N; q = qd = 0;
    for (int i = 1; i < N; ++i) {
        const double sh = std::sin(kPi * i / N);
        q += c.nodes()[i].y * sh; qd += c.velocities()[i].y * sh;
    }
    q *= 2.0 / N; qd *= 2.0 / N;
}

CableParams stringParams(int N) {
    CableParams p; p.L = 1.0; p.N = N; p.EA = 100.0; p.m_l = 0.05; p.w = 0.0; p.g = 9.81;
    return p;
}
LumpedMassCable makeString(const CableParams& p, double A, const Environment& env) {
    const double eps = 0.01;
    LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(p.L * (1 + eps), 0, 0));
    c.setEnvironment(env);
    std::vector<Vec3> r(p.N + 1), v(p.N + 1);
    for (int i = 0; i <= p.N; ++i)
        r[i] = Vec3(i * p.l0() * (1 + eps), A * std::sin(kPi * i / p.N), 0);
    c.setInitialState(r, v);
    return c;
}
}  // namespace

TEST_CASE("added mass: string frequency in water follows m_l -> m_l + Cm rho A1 (1+eps)") {
    const int N = 40; const double eps = 0.01, A = 1e-4;
    CableParams p = stringParams(N);
    p.D0 = p.D1 = 0.01; p.Cm = 1.0; p.Cdn = 0; p.Cdt = 0;
    p.A1 = p.m_l / (1000.0 * p.Cm * (1 + eps));          // makes Cm rho A1 (1+eps) = m_l  -> added mass ratio 1
    Environment env; env.hydro = true; env.surfaceZ = 10.0;   // fully submerged
    LumpedMassCable c = makeString(p, A, env);
    DynOptions o; o.dt = 2e-4; c.setDynOptions(o);
    // Discrete frequency of the LM string with the added mass: omega_N = (2/l0) c sin(pi/2N), c^2 = T / ((m_l+ca)(1+eps))
    const double T = p.EA * eps, mEff = p.m_l + p.Cm * 1000.0 * p.A1 * (1 + eps);
    const double cN = std::sqrt(T / (mEff * (1 + eps)));
    const double omN = 2.0 / p.l0() * cN * std::sin(kPi / (2.0 * N));
    const double omDry = 2.0 / p.l0() * std::sqrt(T / (p.m_l * (1 + eps))) * std::sin(kPi / (2.0 * N));
    // Time of first zero crossing of q(t) = quarter period.
    double q, qd, tPrev = 0, qPrev; modal(c, qPrev, qd);
    double tq = 0;
    for (double t = 2e-3; t < 2.0; t += 2e-3) {
        c.advanceTo(t); modal(c, q, qd);
        if (qPrev > 0 && q <= 0) { tq = tPrev + (t - tPrev) * qPrev / (qPrev - q); break; }
        tPrev = t; qPrev = q;
    }
    REQUIRE(tq > 0);
    const double omMeas = 0.5 * kPi / tq;
    std::printf("\n  added-mass string: omega meas %.4f, predicted %.4f, dry %.4f (ratio %.4f, expected %.4f)\n",
                omMeas, omN, omDry, omMeas / omDry, 1.0 / std::sqrt(2.0));
    CHECK(omMeas == doctest::Approx(omN).epsilon(2e-3));
    CHECK(omMeas / omDry == doctest::Approx(1.0 / std::sqrt(2.0)).epsilon(1e-3));
}

TEST_CASE("quadratic normal drag: modal amplitude decay 1/a = 1/a0 + kappa t") {
    // Slowly-varying-amplitude solution of m_l q_tt = -1/2 Cdn rho D (1+eps) |q_t| q_t + elastic restoring:
    // da/dt = -kappa a^2 with kappa = 16 Cdn rho D (1+eps) omega / (9 pi^2 m_l).
    const int N = 40; const double eps = 0.01, a0 = 0.02;
    CableParams p = stringParams(N);
    p.D0 = p.D1 = 0.0002; p.Cm = 0; p.Cdn = 2.5; p.Cdt = 0.5;
    Environment env; env.hydro = true; env.surfaceZ = 10.0;
    LumpedMassCable c = makeString(p, a0, env);
    DynOptions o; o.dt = 2e-4; c.setDynOptions(o);
    const double T = p.EA * eps;
    const double omN = 2.0 / p.l0() * std::sqrt(T / (p.m_l * (1 + eps))) * std::sin(kPi / (2.0 * N));
    const double kappa = 16.0 * p.Cdn * 1000.0 * p.D0 * (1 + eps) * omN / (9.0 * kPi * kPi * p.m_l);
    std::printf("\n  drag decay: omega %.3f rad/s, kappa %.3f 1/(m s)\n  %8s %12s %12s %10s\n", omN, kappa, "t [s]", "a numeric", "a analytic", "rel err");
    // The law is a first-order averaging result: its error is O(kappa a0 / omega) at the start
    // (here ~3.6 %) and decays with the amplitude, so t = 0.5 s gets that bound and later times 1 %.
    double worst = 0, early = 0;
    for (double t : {0.5, 1.0, 2.0, 3.0, 4.0}) {
        c.advanceTo(t);
        double q, qd; modal(c, q, qd);
        const double a = std::sqrt(q * q + (qd / omN) * (qd / omN));   // envelope from (q, qdot/omega)
        const double aTh = 1.0 / (1.0 / a0 + kappa * t);
        (t < 0.75 ? early : worst) = std::max(t < 0.75 ? early : worst, std::fabs(a - aTh) / aTh);
        std::printf("  %8.2f %12.5e %12.5e %10.2e\n", t, a, aTh, (a - aTh) / aTh);
    }
    CHECK(early < 1.5 * kappa * a0 / omN);
    CHECK(worst < 0.01);
}

TEST_CASE("partly submerged catenary: dry/wet weights and support reactions balance") {
    const int N = 40;
    CableParams p; p.L = 20.5; p.N = N; p.EA = 5e5; p.m_l = 0.5; p.g = 9.81;
    p.w = CableParams::submergedWeight(p.m_l, 7800.0, 1000.0);        // Eq. 3.26
    p.D0 = p.D1 = 0.01;
    Environment env; env.hydro = true; env.surfaceZ = -1.5;            // lowest ~part of the sag is wet
    LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(20, 0, 0)); c.setEnvironment(env);
    CHECK(p.w == doctest::Approx(0.5 * 9.81 * (7800.0 - 1000.0) / 7800.0));
    RelaxOptions ro; ro.forceTol = 1e-8; const RelaxResult rr = c.relaxStatic(ro);
    CHECK(rr.converged);
    double wsum = 0; int wet = 0, dry = 0;
    for (int i = 0; i <= N; ++i) {
        wsum += c.nodeWeight(c.nodes(), i);
        const double ph = c.submergedFraction(c.nodes(), i);
        if (ph >= 1.0) ++wet; else if (ph <= 0.0) ++dry;
    }
    const Vec3 Ft = c.endForce(true), Fa = c.endForce(false);
    std::printf("\n  partly submerged: wet nodes %d, dry nodes %d, sum weight %.4f N, support z sum %.4f N\n",
                wet, dry, wsum, Ft.z + Fa.z);
    CHECK(wet > 0); CHECK(dry > 0);
    CHECK(wsum < p.dryWeight() * p.L);              // buoyancy lightens the wet part
    CHECK(wsum > p.w * p.L);
    CHECK(Ft.z + Fa.z == doctest::Approx(-wsum).epsilon(1e-6));
    CHECK(Ft.x + Fa.x == doctest::Approx(0.0).epsilon(1e-6).scale(std::fabs(Ft.x)));
}

TEST_CASE("chain with seabed touchdown: horizontal tension and suspended weight vs closed form") {
    // Inextensible-chain geometry from a chosen catenary parameter a = H/w: lifted to height h, suspended
    // length s = sqrt(h^2 + 2 a h), touchdown distance a acosh(1 + h/a), plus a lying length ell.
    const double a = 5.0, h = 1.5, ell = 6.0;
    const double s = std::sqrt(h * h + 2 * a * h), xs = a * std::acosh(1.0 + h / a);
    const double L = ell + s, D = ell + xs;
    std::printf("\n  touchdown geometry: L = %.4f m, span = %.4f m, h = %.2f m, a = H/w = %.1f m\n", L, D, h, a);
    for (int N : {30, 60}) {
        CableParams p; p.L = L; p.N = N; p.EA = 1e7; p.m_l = 0.5; p.g = 9.81;
        p.w = CableParams::submergedWeight(p.m_l, 7800.0, 1000.0);
        p.D0 = p.D1 = 0.02;
        p.soil.Ks = 3e9; p.soil.zeta = 1.0; p.soil.mu = 0.3; p.soil.vlim = 0.01;
        Environment env; env.hydro = true; env.seabed = true; env.surfaceZ = 100.0; env.seabedZ = 0.0;
        LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(D, 0, h)); c.setEnvironment(env);
        RelaxOptions ro; ro.forceTol = 1e-6; ro.maxSteps = 6000000;
        const RelaxResult rr = c.relaxStatic(ro);
        double zmin = 0; int onBed = 0;
        for (const Vec3& r : c.nodes()) { zmin = std::min(zmin, r.z); if (r.z < 1e-3) ++onBed; }
        const Vec3 Ft = c.endForce(true);
        std::printf("  N=%3d conv=%d steps=%ld  H num %.4f th %.4f  Vtop num %.4f th %.4f  zmin %.2e  nodes on bed %d\n",
                    N, rr.converged, rr.steps, -Ft.x, p.w * a, -Ft.z, p.w * s, zmin, onBed);
        CHECK(rr.converged);
        CHECK(zmin > -1e-3);                                   // penetration stays tiny
        CHECK(-Ft.x == doctest::Approx(p.w * a).epsilon(0.02));
        CHECK(-Ft.z == doctest::Approx(p.w * s).epsilon(0.02));
    }
}

TEST_CASE("dynamic seabed contact at true Ks: chain dropped onto the bed settles to the static solution") {
    const double a = 5.0, h = 1.5, ell = 6.0;
    const double s = std::sqrt(h * h + 2 * a * h), xs = a * std::acosh(1.0 + h / a);
    const double L = ell + s, D = ell + xs;
    const int N = 40;
    CableParams p; p.L = L; p.N = N; p.EA = 1e7; p.m_l = 0.5; p.g = 9.81;
    p.w = CableParams::submergedWeight(p.m_l, 7800.0, 1000.0);
    p.D0 = p.D1 = 0.02; p.Cdn = 2.5; p.Cdt = 0.5; p.Cm = 0.0;
    p.soil.Ks = 3e9; p.soil.zeta = 1.0; p.soil.mu = 0.3; p.soil.vlim = 0.01;
    p.c_int = 50.0;                                              // damps axial ringing from the stiff soil
    Environment env; env.hydro = true; env.seabed = true; env.surfaceZ = 100.0; env.seabedZ = 0.0;
    LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(D, 0, h)); c.setEnvironment(env);
    RelaxOptions ro; ro.forceTol = 1e-6; ro.maxSteps = 6000000; REQUIRE(c.relaxStatic(ro).converged);
    std::vector<Vec3> r = c.nodes(), v(N + 1);
    for (int i = 1; i < N; ++i) if (r[i].z < 1e-3) r[i].z += 0.3;      // lift the lying part, release from rest
    c.setInitialState(r, v);
    DynOptions o; c.setDynOptions(o);
    c.advanceTo(20.0);
    double zmin = 1e9;
    for (const Vec3& q : c.nodes()) zmin = std::min(zmin, q.z);
    const Vec3 Ft = c.endForce(true);
    std::printf("\n  dynamic drop: dt = %.3e s, soil contact evals %ld, zmin %.2e, H %.3f (static %.3f)\n",
                c.stats().dtUsed, c.stats().soilContactEvals, zmin, -Ft.x, p.w * a);
    CHECK(std::isfinite(zmin));
    CHECK(zmin > -1e-3);
    CHECK(c.stats().soilContactEvals > 0);
    CHECK(-Ft.x == doctest::Approx(p.w * a).epsilon(0.05));
}
