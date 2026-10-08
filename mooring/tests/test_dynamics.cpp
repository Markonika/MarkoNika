// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Test 2 (vibrating string, Paredes 2016 sec. 3.5.1, Eqs. 3.57-3.58) and test 3 (energy).
#include <cmath>
#include <cstdio>
#include <vector>
#include "doctest.h"
#include "mooring/catenary.hpp"
#include "mooring/lumped_mass_cable.hpp"

using namespace mooring;

namespace {
const double kPi = 3.14159265358979323846;

// 0.5 m cable stretched to 1 m (eps = 1), m_l = 1 kg/m, pretension 1 N => EA = 1 N, no gravity.
// Linear Lagrangian wave equation: c^2 = tau / (m_l (1 + eps)), mode sin(pi s / L0),
// omega = pi c / L0.  (Eq. 3.58 as printed uses L(1+eps); the physically consistent form with the
// unstretched length L0 is used here, see docs/assumptions.md.)
double stringError(int N, Scheme scheme, double tEnd, double dt) {
    const double L0 = 0.5, eps = 1.0, tau = 1.0, A = 1e-4;
    CableParams p; p.L = L0; p.N = N; p.EA = tau / eps; p.m_l = 1.0; p.w = 0.0;
    LumpedMassCable cable(p, Vec3(0, 0, 0), Vec3(L0 * (1 + eps), 0, 0));
    std::vector<Vec3> r(N + 1), v(N + 1);
    for (int i = 0; i <= N; ++i) {
        const double s = i * p.l0();
        r[i] = Vec3(s * (1 + eps), A * std::sin(kPi * s / L0), 0);
    }
    cable.setInitialState(r, v);
    DynOptions o; o.scheme = scheme; o.dt = dt; cable.setDynOptions(o);
    cable.advanceTo(tEnd);
    const double c = std::sqrt(tau / (p.m_l * (1 + eps))), om = kPi * c / L0;
    double num = 0, den = 0;
    for (int i = 0; i <= N; ++i) {
        const double ex = A * std::cos(om * tEnd) * std::sin(kPi * i * p.l0() / L0);
        const double d = cable.nodes()[i].y - ex;
        num += d * d; den += A * A * std::sin(kPi * i * p.l0() / L0) * std::sin(kPi * i * p.l0() / L0);
    }
    return std::sqrt(num / den);
}
}  // namespace

TEST_CASE("vibrating string: L2 error vs N (RK4 and velocity-Verlet), 2 s = 1.41 periods") {
    std::printf("\n  string: period = %.4f s, run to t = 2 s, dt = 1e-4\n  %6s %14s %14s\n",
                2 * kPi / (kPi * std::sqrt(0.5) / 0.5), "N", "RK4", "Verlet");
    const std::vector<int> Ns = {5, 10, 20, 40, 80};
    std::vector<double> eR, eV;
    for (int N : Ns) {
        eR.push_back(stringError(N, Scheme::RK4, 2.0, 1e-4));
        eV.push_back(stringError(N, Scheme::Verlet, 2.0, 1e-4));
        std::printf("  %6d %14.4e %14.4e\n", N, eR.back(), eV.back());
    }
    for (size_t i = 1; i < Ns.size(); ++i) { CHECK(eR[i] < eR[i - 1]); CHECK(eV[i] < eV[i - 1]); }
    const double ordR = std::log(eR[1] / eR[4]) / std::log(8.0);
    const double ordV = std::log(eV[1] / eV[4]) / std::log(8.0);
    std::printf("  observed order (N=10..80): RK4 %.2f, Verlet %.2f\n", ordR, ordV);
    CHECK(ordR > 1.8);
    CHECK(ordV > 1.8);
    CHECK(eR.back() < 1e-3);
}

namespace {
// Relaxed catenary (EA = 200 kN, Paredes static case) kicked by an out-of-equilibrium
// perturbation, then integrated with no damping and no dissipation.
struct EnergyRun { double rel_drift_max, E_range; double dt; };

EnergyRun energyRun(Scheme sch, double tEnd, double cfl) {
    const int N = 40;
    CableParams p; p.L = 100.5; p.N = N; p.EA = 200e3; p.m_l = 1.738; p.w = 1.738 * 9.81; p.c_int = 0.0;
    LumpedMassCable cable(p, Vec3(0, 0, 0), Vec3(100, 0, 0));
    RelaxOptions ro; ro.forceTol = 1e-9; cable.relaxStatic(ro);
    const double Eeq = cable.energy();
    std::vector<Vec3> r = cable.nodes(), v(N + 1);
    for (int i = 1; i < N; ++i) {
        const double s = std::sin(kPi * i / N);
        r[i].z += 1.0 * s; r[i].y += 0.5 * s;          // 3D perturbation, released from rest
    }
    cable.setInitialState(r, v);
    DynOptions o; o.scheme = sch; o.cfl = cfl; cable.setDynOptions(o);
    const double E0 = cable.energy();
    EnergyRun out{0, E0 - Eeq, cable.stableDt()};
    for (double t = 5.0; t <= tEnd + 1e-9; t += 5.0) {
        cable.advanceTo(t);
        out.rel_drift_max = std::max(out.rel_drift_max, std::fabs(cable.energy() - E0) / (E0 - Eeq));
    }
    return out;
}
}  // namespace

TEST_CASE("energy conservation: frictionless undamped catenary, 200 s (about 40 periods)") {
    const EnergyRun rk5 = energyRun(Scheme::RK4, 200.0, 0.5), vl5 = energyRun(Scheme::Verlet, 200.0, 0.5);
    const EnergyRun rk = energyRun(Scheme::RK4, 200.0, 0.25), vl = energyRun(Scheme::Verlet, 200.0, 0.25);
    std::printf("\n  excitation energy E0-Eeq = %.4f J\n", rk.E_range);
    std::printf("  max |E-E0|/(E0-Eeq), CFL 0.50 (default): RK4 %.3e  Verlet %.3e\n", rk5.rel_drift_max, vl5.rel_drift_max);
    std::printf("  max |E-E0|/(E0-Eeq), CFL 0.25          : RK4 %.3e  Verlet %.3e  (dt = %.3e s)\n",
                rk.rel_drift_max, vl.rel_drift_max, rk.dt);
    // Tolerances chosen after a first run (RK4 error ~ dt^4, Verlet ~ dt^2) with margin; stated, not tuned away.
    CHECK(rk.rel_drift_max < 1e-4);
    CHECK(vl.rel_drift_max < 2e-4);
}

TEST_CASE("forceOnBody returns the static end tension for a fixed fairlead") {
    const int N = 40;
    CableParams p; p.L = 100.5; p.N = N; p.EA = 200e3; p.m_l = 1.738; p.w = 1.738 * 9.81; p.c_int = 5e3;
    LumpedMassCable cable(p, Vec3(0, 0, 0), Vec3(100, 0, 0));
    RelaxOptions ro; ro.forceTol = 1e-9; cable.relaxStatic(ro);
    const auto cat = ElasticCatenary::solve(p.L, p.EA, p.w, 100.0, 0.0);
    Vec3 F;
    for (double t = 0.0; t <= 3.0 + 1e-9; t += 0.25) F = cable.forceOnBody(Vec3(100, 0, 0), Vec3(), t);
    // The line pulls the fairlead toward itself and down: (-H, -wL/2) in static equilibrium.
    CHECK(F.x == doctest::Approx(-cat.H).epsilon(2e-3));
    CHECK(F.z == doctest::Approx(-0.5 * p.w * p.L).epsilon(1e-3));
    const Vec3 Fa = cable.endForce(false);
    CHECK(F.z + Fa.z == doctest::Approx(-p.w * p.L).epsilon(1e-6));   // vertical balance: both supports carry the weight
    CHECK(F.x + Fa.x == doctest::Approx(0.0).epsilon(1e-6).scale(cat.H));
    CHECK(cable.stats().steps > 0);
}
