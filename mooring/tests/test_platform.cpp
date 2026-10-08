// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Milestone 7: 6-DOF rigid-body platform, multi-line fairlead coupling, sub-stepping.
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>
#include "doctest.h"
#include "mooring/coupled_runner.hpp"
#include "mooring/lumped_mass_cable.hpp"
#include "mooring/rigid_body.hpp"

using namespace mooring;

namespace { const double kPi = 3.14159265358979323846; }

TEST_CASE("rotation kinematics: orthonormal exponential map and rigid-body fairlead velocity") {
    const Vec3 th(0.3, -0.2, 0.5);
    const Mat3 R = rotationFromVector(th);
    // R R^T = I, det R = 1
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
        double s = 0; for (int k = 0; k < 3; ++k) s += R[3 * i + k] * R[3 * j + k];
        CHECK(s == doctest::Approx(i == j ? 1.0 : 0.0).epsilon(1e-14).scale(1.0));
    }
    const double det = R[0] * (R[4] * R[8] - R[5] * R[7]) - R[1] * (R[3] * R[8] - R[5] * R[6]) + R[2] * (R[3] * R[7] - R[4] * R[6]);
    CHECK(det == doctest::Approx(1.0).epsilon(1e-14));
    // 90 deg about z maps x -> y
    const Vec3 y = rotate(rotationFromVector(Vec3(0, 0, 0.5 * kPi)), Vec3(1, 0, 0));
    CHECK(norm(y - Vec3(0, 1, 0)) < 1e-14);
    // Small-angle series branch is continuous with the closed form.
    const Mat3 Ra = rotationFromVector(Vec3(0.99e-4, 0, 0)), Rb = rotationFromVector(Vec3(1.01e-4, 0, 0));
    CHECK(std::fabs(Ra[5] - Rb[5]) < 3e-6);
    // Velocity of a body point: v = u' + w x (R a), to first order in the rotation rate: compare with a central difference.
    BodyParams bp; bp.cgRef = Vec3(1, 2, 3);
    RigidBody6DOF b(bp);
    b.xi = {0.1, -0.2, 0.05, 0.02, -0.03, 0.04};
    b.xiDot = {0.3, 0.1, -0.2, 0.05, 0.02, -0.04};
    const Vec3 a(0.4, -0.3, 0.2);
    const double h = 1e-6;
    RigidBody6DOF bp1 = b, bm = b;
    for (int k = 0; k < 6; ++k) { bp1.xi[k] += h * b.xiDot[k]; bm.xi[k] -= h * b.xiDot[k]; }
    const Vec3 vfd = (bp1.pointPosition(a) - bm.pointPosition(a)) / (2 * h);
    // rotation vector rate vs angular velocity differ at O(theta * theta') (~1e-3 here): documented small-angle model
    CHECK(norm(vfd - b.pointVelocity(a)) < 5e-3 * norm(b.pointVelocity(a)) + 1e-9);
}

TEST_CASE("free heave decay with constant A, B, C matches the damped-oscillator solution") {
    BodyParams p; p.mass = 10.0; p.cgRef = Vec3(0, 0, 0);
    p.A(2, 2) = 5.0; p.B(2, 2) = 3.0; p.C(2, 2) = 600.0;
    CoupledSystem sys(p);
    sys.body().xi[2] = 0.02;
    sys.begin(0.0);
    const double M = 15.0, w0 = std::sqrt(600.0 / M), z = 3.0 / (2.0 * std::sqrt(600.0 * M)), wd = w0 * std::sqrt(1 - z * z);
    std::vector<double> errs;
    for (double dt : {2e-3, 1e-3, 5e-4}) {
        CoupledSystem s(p); s.body().xi[2] = 0.02; s.begin(0.0);
        s.advanceTo(5.0, dt);
        const double t = s.time();
        const double ex = 0.02 * std::exp(-z * w0 * t) * (std::cos(wd * t) + z * w0 / wd * std::sin(wd * t));
        const double err = std::fabs(s.body().xi[2] - ex) / 0.02;
        std::printf("  free heave decay, dt = %.0e s: error %.2e (relative to amplitude)\n", dt, err);
        errs.push_back(err);
    }
    // Second order in dt: the error drops by ~4 when dt is halved.
    CHECK(errs[0] / errs[1] == doctest::Approx(4.0).epsilon(0.1));
    CHECK(errs[1] / errs[2] == doctest::Approx(4.0).epsilon(0.1));
    CHECK(errs[2] < 3e-6);
}

namespace {
// Body of mass M hanging from a fixed point on one nearly vertical line (anchor above, fairlead = body CG).
struct Hang {
    CoupledSystem sys;
    LumpedMassCable* cable{nullptr};
    explicit Hang(const BodyParams& bp) : sys(bp) {}
};

std::unique_ptr<Hang> makeHang(double M, double L0, double EA, double ml, int N, double dtLine, double* Ltot) {
    BodyParams bp; bp.mass = M; bp.cgRef = Vec3(0, 0, -L0 * 1.005);               // reference pose slightly taut (Newton needs a taut line)
    bp.F0[2] = -M * 9.81;                                                           // gravity as a constant load (no hydrostatics)
    bp.inertia = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    auto h = std::make_unique<Hang>(bp);
    CableParams p; p.L = L0; p.N = N; p.EA = EA; p.m_l = ml; p.w = ml * 9.81;
    auto cab = std::make_unique<LumpedMassCable>(p, Vec3(0, 0, 0), Vec3(0, 0, -L0 * 1.005));
    if (dtLine > 0.0) { DynOptions o; o.dt = dtLine; cab->setDynOptions(o); }   // 0 = automatic CFL step
    h->cable = cab.get();
    h->sys.addLine(std::move(cab), Vec3(0, 0, 0), "hang");
    const EquilibriumResult er = h->sys.solveEquilibrium(1e-8);
    REQUIRE(er.converged);
    if (Ltot) *Ltot = norm(h->sys.fairleadPosition(0));
    return h;
}
}  // namespace

TEST_CASE("body on a vertical line: heave frequency vs exact end-mass spring; sub-step convergence") {
    const double M = 2.0, L0 = 10.0, EA = 2000.0, ml = 0.001; const int N = 20;
    const double dtLine = 1e-4;
    // Exact: beta tan(beta) = m_s / M, omega = beta c / L, c = sqrt(EA/m_l).
    const double ms = ml * L0, c = std::sqrt(EA / ml);
    double beta = std::sqrt(ms / M);
    for (int i = 0; i < 50; ++i) { const double f = beta * std::tan(beta) - ms / M, d = std::tan(beta) + beta / (std::cos(beta) * std::cos(beta)); beta -= f / d; }
    const double omExact = beta * c / L0;
    std::printf("\n  heave: omega exact (spring with distributed mass) %.5f rad/s, rigid-spring estimate %.5f\n", omExact, std::sqrt(EA / L0 / M));

    auto run = [&](int ratio, double* omOut, double* zEnd) {
        double Ltot; auto h = makeHang(M, L0, EA, ml, N, dtLine, &Ltot);
        const double z0 = h->sys.body().xi[2];
        h->sys.body().xi[2] += 0.01;
        h->sys.staticResidual();                                   // lines follow the displaced fairlead (static shape)
        h->sys.begin(0.0);
        const double dt = ratio * dtLine;
        // measure the half period between the first two zero crossings of (z - z0) by interpolation
        double tPrev = 0, zPrev = 0.01, tc[2] = {0, 0}; int nc = 0;
        while (h->sys.time() < 3.0 && nc < 2) {
            h->sys.step(dt);
            const double zz = h->sys.body().xi[2] - z0 - (h->sys.body().xi[2] > 1e9 ? 0 : 0);
            const double zr = zz;
            if ((zPrev > 0) != (zr > 0)) { tc[nc++] = tPrev + (h->sys.time() - tPrev) * zPrev / (zPrev - zr); }
            tPrev = h->sys.time(); zPrev = zr;
        }
        *omOut = kPi / (tc[1] - tc[0]);
        *zEnd = h->sys.body().xi[2];
        return h->sys.report().subStepRatio;
    };
    double om1, ze;
    run(1, &om1, &ze);
    std::printf("  ratio  measured omega    rel error vs exact   rel diff vs ratio 1\n");
    double prev = 0;
    for (int ratio : {1, 2, 5, 10, 20, 50}) {
        double om, z;
        const int rep = run(ratio, &om, &z);
        std::printf("  %5d  %14.5f  %18.2e  %18.2e   (reported ratio %d)\n", ratio, om, (om - omExact) / omExact, (om - om1) / om1, rep);
        CHECK(rep >= ratio - 1);
        if (ratio == 1) CHECK(om == doctest::Approx(omExact).epsilon(2e-3));
        if (ratio == 50) CHECK(std::fabs(om - om1) / om1 < 5e-3);
        prev = om;
    }
    (void)prev;
}

TEST_CASE("body on a vertical line: pendulum swing, omega = sqrt(g/L); sub-step convergence") {
    const double M = 2.0, L0 = 10.0, EA = 2000.0, ml = 2e-4; const int N = 20;
    double Ltot;
    auto run = [&](int ratio, double th0) {
        auto h = makeHang(M, L0, EA, ml, N, 0.0, &Ltot);
        const double dtLine = h->cable->stableDt();
        const Vec3 anchor(0, 0, 0);
        const Vec3 eq = h->sys.body().cg();
        // displace along the circle of radius |cg - anchor|
        const double Lp = norm(eq - anchor);
        h->sys.body().xi[0] = Lp * std::sin(th0);
        h->sys.body().xi[2] += Lp * (1 - std::cos(th0));
        h->sys.staticResidual();
        h->sys.begin(0.0);
        const double dt = ratio * dtLine;
        double tPrev = 0, xPrev = Lp * std::sin(th0), tc[3] = {0, 0, 0}; int nc = 0;
        while (h->sys.time() < 12.0 && nc < 3) {
            h->sys.step(dt);
            const double x = h->sys.body().xi[0];
            if ((xPrev > 0) != (x > 0)) tc[nc++] = tPrev + (h->sys.time() - tPrev) * xPrev / (xPrev - x);
            tPrev = h->sys.time(); xPrev = x;
        }
        return std::make_pair(2.0 * (tc[2] - tc[0]) / 2.0 * 1.0, Lp);       // full period = tc[2]-tc[0]
    };
    const double th0 = 0.02;
    const auto base = run(1, th0);
    const double Texact = 2.0 * kPi * std::sqrt(base.second / 9.81) * (1.0 + th0 * th0 / 16.0);
    std::printf("\n  pendulum: L = %.5f m, exact period %.5f s (incl. finite-amplitude correction)\n  ratio   period   rel error vs exact\n", base.second, Texact);
    CHECK(base.first == doctest::Approx(Texact).epsilon(1e-3));
    for (int ratio : {1, 5, 25, 100}) {
        const auto r = run(ratio, th0);
        std::printf("  %5d  %8.5f  %12.2e\n", ratio, r.first, (r.first - Texact) / Texact);
        CHECK(r.first == doctest::Approx(Texact).epsilon(1.5e-3));
    }
}

TEST_CASE("three lines at 120 deg: static equilibrium with independent force and moment balance") {
    BodyParams bp; bp.mass = 35.5; bp.cgRef = Vec3(0, 0, 0);
    bp.C(2, 2) = 2040.0; bp.C(3, 3) = 150.0; bp.C(4, 4) = 150.0;                  // heave, roll, pitch hydrostatics
    CoupledSystem sys(bp);
    const double R = 4.0, depth = 2.5;
    const double pretL[3] = {0.0, 0.03, -0.02};                                      // unequal rest lengths -> asymmetric equilibrium
    for (int k = 0; k < 3; ++k) {
        const double ang = 2.0 * kPi * k / 3.0;
        const Vec3 anchor(R * std::cos(ang), R * std::sin(ang), -depth);
        const Vec3 fair(0.3 * std::cos(ang), 0.3 * std::sin(ang), -0.15);            // body coordinates, relative to the CG
        const double chord = norm(Vec3(0, 0, 0) + fair - anchor);
        CableParams p; p.L = chord - 0.02 + pretL[k]; p.N = 16; p.EA = 3000.0; p.m_l = 0.01; p.w = 0.0;
        sys.addLine(std::make_unique<LumpedMassCable>(p, anchor, fair), fair, "leg" + std::to_string(k + 1));
    }
    const EquilibriumResult er = sys.solveEquilibrium(1e-7);
    std::printf("\n  equilibrium: converged %d in %d Newton iterations, |F| %.2e N, |M| %.2e N m\n", er.converged, er.iterations, er.forceResidual, er.momentResidual);
    CHECK(er.converged);
    const Vec6& xi = sys.body().xi;
    std::printf("  xi = [% .5f % .5f % .5f | % .5f % .5f % .5f]\n", xi[0], xi[1], xi[2], xi[3], xi[4], xi[5]);
    sys.begin(0.0);                                                                  // evaluates lines at the equilibrium pose
    // Independent balance: sum of line forces + restoring force = 0; sum of moments about the CG (+ restoring moments) = 0.
    Vec3 F; Vec3 M;
    for (int k = 0; k < 3; ++k) {
        F += sys.lineForce(k);
        M += cross(sys.fairleadPosition(k) - sys.body().cg(), sys.lineForce(k));
    }
    const Vec6 Cx = bp.C * xi;
    std::printf("  line force sum (%.6f, %.6f, %.6f) N vs restoring (%.6f, %.6f, %.6f)\n", F.x, F.y, F.z, Cx[0], Cx[1], Cx[2]);
    CHECK(F.x - Cx[0] == doctest::Approx(0.0).scale(1.0).epsilon(1e-6));
    CHECK(F.y - Cx[1] == doctest::Approx(0.0).scale(1.0).epsilon(1e-6));
    CHECK(F.z - Cx[2] == doctest::Approx(0.0).scale(1.0).epsilon(1e-6));
    CHECK(M.x - Cx[3] == doctest::Approx(0.0).scale(1.0).epsilon(1e-6));
    CHECK(M.y - Cx[4] == doctest::Approx(0.0).scale(1.0).epsilon(1e-6));
    CHECK(M.z - Cx[5] == doctest::Approx(0.0).scale(1.0).epsilon(1e-6));
    CHECK(xi[2] < 0.0);                                                              // taut legs pull the body down
    // Symmetric legs: no horizontal drift, no tilt, no yaw.
    CoupledSystem sym(bp);
    for (int k = 0; k < 3; ++k) {
        const double ang = 2.0 * kPi * k / 3.0;
        const Vec3 anchor(R * std::cos(ang), R * std::sin(ang), -depth), fair(0.3 * std::cos(ang), 0.3 * std::sin(ang), -0.15);
        CableParams p; p.L = norm(fair - anchor) - 0.02; p.N = 16; p.EA = 3000.0; p.m_l = 0.01; p.w = 0.0;
        sym.addLine(std::make_unique<LumpedMassCable>(p, anchor, fair), fair);
    }
    REQUIRE(sym.solveEquilibrium(1e-7).converged);
    const Vec6& xs = sym.body().xi;
    std::printf("  symmetric: xi = [% .2e % .2e % .6f | % .2e % .2e % .2e]\n", xs[0], xs[1], xs[2], xs[3], xs[4], xs[5]);
    CHECK(std::fabs(xs[0]) < 1e-7); CHECK(std::fabs(xs[1]) < 1e-7);
    CHECK(std::fabs(xs[3]) < 1e-7); CHECK(std::fabs(xs[4]) < 1e-7); CHECK(std::fabs(xs[5]) < 1e-7);
    CHECK(xs[2] < 0.0);
}

TEST_CASE("coupling stability: velocity-Verlet limit dt_body ~ 2/omega_s with the static line stiffness") {
    // omega_s^2 = EA / (L (M + m_line/3)); measured limits for several (EA, N, M) are in docs/assumptions.md.
    auto blowsUp = [](double M, double dtBody) {
        const double L0 = 10.0, EA = 2000.0, ml = 0.001; const int N = 20;
        BodyParams bp; bp.mass = M; bp.cgRef = Vec3(0, 0, -L0 * 1.005); bp.F0[2] = -M * 9.81;
        CoupledSystem sys(bp);
        CableParams p; p.L = L0; p.N = N; p.EA = EA; p.m_l = ml; p.w = ml * 9.81;
        sys.addLine(std::make_unique<LumpedMassCable>(p, Vec3(0, 0, 0), Vec3(0, 0, -L0 * 1.005)), Vec3(), "l");
        REQUIRE(sys.solveEquilibrium(1e-8).converged);
        sys.body().xi[2] += 1e-3; sys.staticResidual(); sys.begin(0.0);
        for (double t = 0; t < 2.0; t += dtBody) {
            sys.step(dtBody);
            if (!(std::fabs(sys.body().xi[2]) < 0.1)) return true;
        }
        return false;
    };
    std::printf("\n  stability of the explicit partitioned coupling (line: EA = 2000 N, L = 10 m, N = 20)\n  %6s %12s %12s %10s\n", "M [kg]", "2/omega_s", "dt_body", "unstable");
    for (double M : {2.0, 0.2}) {
        const double lim = 2.0 / std::sqrt(2000.0 / (10.0 * (M + 0.01 / 3.0)));
        const bool lo = blowsUp(M, 0.4 * lim), hi = blowsUp(M, 2.0 * lim);
        std::printf("  %6.1f %12.4f %12.4f %10d\n  %6s %12s %12.4f %10d\n", M, lim, 0.4 * lim, int(lo), "", "", 2.0 * lim, int(hi));
        CHECK(!lo);
        CHECK(hi);
    }
}

TEST_CASE("JSON platform runner: body on a line reproduces the pendulum period") {
    const nlohmann::json cfg = nlohmann::json::parse(R"({
      "body": { "mass_kg": 2.0, "inertia_diag_kg_m2": [1, 1, 1], "cg_ref_m": [0, 0, -10.05], "F0": [0, 0, -19.62, 0, 0, 0] },
      "lines": [ { "name": "string", "anchor_m": [0, 0, 0], "fairlead_body_m": [0, 0, 0],
                   "line": { "length_m": 10.0, "segments": 20, "EA_N": 2000.0, "mass_per_length_kg_m": 2e-4, "weight_per_length_N_m": 0.00196 } } ],
      "initial": { "equilibrium": true, "xi_offset": [0.202, 0, 0.002, 0, 0, 0] },
      "numerics": { "dt_body_s": 0.004, "t_end_s": 10.0 },
      "output": { "dt_out_s": 0.004, "write": false } })");
    const CoupledResult r = runCoupledCase(cfg);
    REQUIRE(r.equilibrium.converged);
    REQUIRE(r.finite);
    // pendulum length: anchor to CG at equilibrium (CG = fairlead)
    const double Lp = norm(Vec3(0, 0, -10.05) + Vec3(r.xiStart[0] - 0.202, r.xiStart[1], r.xiStart[2] - 0.002));
    std::vector<double> tc;
    for (size_t i = 1; i < r.samples.size() && tc.size() < 3; ++i) {
        const double a = r.samples[i - 1].xi[0], b = r.samples[i].xi[0];
        if ((a > 0) != (b > 0)) tc.push_back(r.samples[i - 1].t + (r.samples[i].t - r.samples[i - 1].t) * a / (a - b));
    }
    REQUIRE(tc.size() == 3);
    const double th0 = 0.202 / Lp;
    const double Texact = 2.0 * kPi * std::sqrt(Lp / 9.81) * (1.0 + th0 * th0 / 16.0);
    std::printf("\n  JSON platform pendulum: L = %.5f m, period %.5f s, exact %.5f s, error %.2e; sub-step ratio %d (dt_body %.1e, shortest line step %.1e)\n",
                Lp, tc[2] - tc[0], Texact, (tc[2] - tc[0] - Texact) / Texact, r.report.subStepRatio, r.report.dtBody, r.report.dtLineMin);
    CHECK((tc[2] - tc[0]) == doctest::Approx(Texact).epsilon(1.5e-3));
    CHECK(r.report.subStepRatio > 1);
}
