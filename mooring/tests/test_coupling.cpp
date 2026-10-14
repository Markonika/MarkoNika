// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Implicit body integration with iterated coupling (ImplicitCoupling): a stiff line on a light body, where the explicit partitioned scheme is
// limited to dt_body < 2/omega_s.
#include <cmath>
#include <cstdio>
#include <memory>
#include "doctest.h"
#include "mooring/lumped_mass_cable.hpp"
#include "mooring/rigid_body.hpp"

using namespace mooring;

namespace {
const double kPi = 3.14159265358979323846;

struct StiffHang {
    std::unique_ptr<CoupledSystem> sys;
    LumpedMassCable* cable{nullptr};
};
// Body of mass M hanging on a stiff line (EA = 2e6 N: omega_s = sqrt(EA / (L0 M)) = 316 rad/s, explicit limit dt_body < 6.3 ms).
StiffHang makeStiffHang(bool implicitCoupling, double M = 2.0, double L0 = 10.0, double EA = 2e6, int N = 10) {
    StiffHang h;
    BodyParams bp; bp.mass = M; bp.cgRef = Vec3(0, 0, -L0 * 1.005);
    bp.F0[2] = -M * 9.81;
    bp.inertia = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    h.sys = std::make_unique<CoupledSystem>(bp);
    CableParams p; p.L = L0; p.N = N; p.EA = EA; p.m_l = 0.001; p.w = 0.001 * 9.81;
    auto cab = std::make_unique<LumpedMassCable>(p, Vec3(0, 0, 0), Vec3(0, 0, -L0 * 1.005));
    h.cable = cab.get();
    h.sys->addLine(std::move(cab), Vec3(0, 0, 0), "hang");
    REQUIRE(h.sys->solveEquilibrium(1e-8).converged);
    if (implicitCoupling) { ImplicitCoupling ic; ic.enabled = true; h.sys->setImplicitCoupling(ic); }
    return h;
}
}  // namespace

TEST_CASE("Implicit coupling: stiff line on a light body, slow heave forcing at 8 x the explicit stability limit") {
    const double Fa = 5.0, T = 2.0, w = 2.0 * kPi / T, tEnd = 1.5;
    auto run = [&](bool impl, double dt, double* zEnd, bool* finite, long* iters) {
        StiffHang h = makeStiffHang(impl);
        const double z0 = h.sys->body().xi[2];
        h.sys->setExcitation([&](double t) { Vec6 F{}; F[2] = Fa * std::sin(w * t); return F; });
        h.sys->begin(0.0);
        *finite = true;
        while (h.sys->time() < tEnd - 1e-9 && *finite) {
            h.sys->step(std::min(dt, tEnd - h.sys->time()));
            if (!std::isfinite(h.sys->body().xi[2]) || std::fabs(h.sys->body().xi[2] - z0) > 1.0) *finite = false;
        }
        *zEnd = h.sys->body().xi[2] - z0;
        *iters = h.sys->report().couplingIterations;
    };
    double zRef, zImp, zExp; bool fRef, fImp, fExp; long iRef, iImp, iExp;
    run(false, 2e-3, &zRef, &fRef, &iRef);                    // explicit, below the limit (reference)
    run(true, 5e-2, &zImp, &fImp, &iImp);                     // implicit, 25 x larger
    run(false, 5e-2, &zExp, &fExp, &iExp);                    // explicit at the same step: unstable
    std::printf("\n  heave at t = %.1f s: explicit dt 2 ms %.6e m; implicit dt 50 ms %.6e m (%ld line advances); explicit dt 50 ms %s\n", tEnd, zRef, zImp, iImp, fExp ? "finite" : "diverged");
    REQUIRE(fRef);
    REQUIRE(fImp);
    CHECK(zImp == doctest::Approx(zRef).scale(0).epsilon(0.02));
    CHECK_FALSE(fExp);
    CHECK(iImp < 4 * (long)(tEnd / 5e-2) + 20);               // about two line advances per step or fewer
}

TEST_CASE("Implicit coupling: pendulum swing with large body steps; accuracy limited by the line's own vibration period, converging for smaller steps") {
    const double L0 = 10.0;
    auto run = [&](bool impl, double dt, double tEnd) {
        StiffHang h = makeStiffHang(impl, 2.0, L0, 2e4, 10);   // moderately stiff: explicit limit 63 ms; the lowest transverse line mode has a period of about 0.14 s
        h.sys->body().xi[0] = 0.2;
        h.sys->staticResidual();
        h.sys->begin(0.0);
        h.sys->advanceTo(tEnd, dt);
        return h.sys->body().xi[0];
    };
    const double ref = run(false, 5e-3, 3.0);                  // explicit reference
    std::printf("\n  pendulum x(3 s): explicit dt 5 ms %.6f", ref);
    double err[3]; int k = 0;
    for (double dt : {0.2, 0.1, 0.025}) {
        const double x = run(true, dt, 3.0);
        err[k++] = std::fabs(x - ref);
        std::printf("\n  implicit dt %.3f: %.6f (diff %.2e = %.2f %% of the swing amplitude)", dt, x, std::fabs(x - ref), 100 * std::fabs(x - ref) / 0.2);
    }
    std::printf("\n");
    CHECK(err[0] < 0.015 * 0.2);              // dt = 0.2 s is 1.4 x the line mode period: swing still right to 1.5 %
    CHECK(err[1] < 0.015 * 0.2);
    CHECK(err[2] < 0.002 * 0.2);              // once the line vibration is resolved the second-order Newmark error is small
}
