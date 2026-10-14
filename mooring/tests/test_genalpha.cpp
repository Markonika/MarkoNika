// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Fully implicit generalised-alpha / Newmark integration of the line (Scheme::GenAlpha): temporal order, damped-mode decay at steps far above the
// explicit limits, independence from the Jacobian refresh interval, and agreement with the explicit scheme on a Lopez-Olocco case.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>
#include "doctest.h"
#include "mooring/case_runner.hpp"
#include "mooring/lumped_mass_cable.hpp"

using namespace mooring;
using json = nlohmann::json;

namespace {
const double kPi = 3.14159265358979323846;

// The vibrating string of test_dynamics.cpp (0.5 m cable stretched to 1 m, no gravity): final positions after tEnd for the given scheme and step.
std::vector<Vec3> stringRun(Scheme sch, double dt, double rho, double tEnd, int N = 40) {
    const double L0 = 0.5, eps = 1.0, tau = 1.0, A = 1e-4;
    CableParams p; p.L = L0; p.N = N; p.EA = tau / eps; p.m_l = 1.0; p.w = 0.0;
    LumpedMassCable cable(p, Vec3(0, 0, 0), Vec3(L0 * (1 + eps), 0, 0));
    std::vector<Vec3> r(N + 1), v(N + 1);
    for (int i = 0; i <= N; ++i) { const double s = i * p.l0(); r[i] = Vec3(s * (1 + eps), A * std::sin(kPi * s / L0), 0); }
    cable.setInitialState(r, v);
    DynOptions o; o.scheme = sch; o.dt = dt; o.rhoInf = rho; cable.setDynOptions(o);
    cable.advanceTo(tEnd);
    return cable.nodes();
}
double diff(const std::vector<Vec3>& a, const std::vector<Vec3>& b) {
    double s = 0, n = 0;
    for (size_t i = 0; i < a.size(); ++i) { s += dot(a[i] - b[i], a[i] - b[i]); n += 1e-8; }
    return std::sqrt(s / n);
}

// Longitudinal mode of a taut chain with stiffness-proportional damping beta (as in test_implicit.cpp); mode energy at the given times.
std::vector<double> modeEnergy(double beta, const DynOptions& o, const std::vector<double>& times, double* dtOut = nullptr) {
    const double L = 10.0, EA = 1e5, ml = 0.1, e0 = 1e-3, u0 = 1e-5; const int N = 40;
    CableParams p; p.L = L; p.N = N; p.EA = EA; p.m_l = ml; p.w = 0.0; p.g = 0.0; p.c_int = beta * EA;
    const double l0 = p.l0();
    LumpedMassCable c(p, Vec3(0, 0, 0), Vec3(L * (1 + e0), 0, 0));
    std::vector<Vec3> r(N + 1), v(N + 1);
    for (int i = 0; i <= N; ++i) r[i] = Vec3(i * l0 * (1 + e0) + u0 * std::sin(kPi * i * l0 / L), 0, 0);
    c.setInitialState(r, v);
    c.setDynOptions(o);
    if (dtOut) *dtOut = c.stableDt();
    const double E0 = N * 0.5 * EA * l0 * e0 * e0;
    std::vector<double> E;
    for (double t : times) { c.advanceTo(t); E.push_back(c.energy() - E0); }
    return E;
}

json loadLopez(const char* name) {
    std::ifstream in(std::string(MOORING_SOURCE_DIR) + "/examples/lopezolocco/" + name + ".json");
    REQUIRE(in.good());
    json c = json::parse(in);
    c["output"]["write"] = false;
    return c;
}
double lastCyclesMax(const CaseResult& r) {
    const int n = int(r.cycles.size()) - 1;
    double s = 0; for (int k = n - 4; k < n; ++k) s += r.cycles[k].topMax;
    return s / 4.0;
}
}  // namespace

TEST_CASE("GenAlpha: Newmark (rho_inf = 1) is second order in time on the vibrating string") {
    const double tEnd = 2.0;
    const std::vector<Vec3> ref = stringRun(Scheme::RK4, 2e-4, 1.0, tEnd);          // temporal error of the reference is negligible
    double prev = 0;
    for (double dt : {0.04, 0.02, 0.01}) {
        const double e = diff(stringRun(Scheme::GenAlpha, dt, 1.0, tEnd), ref);
        std::printf("\n  Newmark dt=%.3f: rel. difference to the RK4 reference %.3e%s", dt, e, prev > 0 ? (", ratio " + std::to_string(prev / e)).c_str() : "");
        if (prev > 0) { CHECK(prev / e > 3.3); CHECK(prev / e < 4.7); }
        prev = e;
    }
}

TEST_CASE("GenAlpha: rho_inf < 1 stays second order and damps only high frequencies") {
    const double tEnd = 2.0;
    const std::vector<Vec3> ref = stringRun(Scheme::RK4, 2e-4, 1.0, tEnd);
    const double e1 = diff(stringRun(Scheme::GenAlpha, 0.02, 0.8, tEnd), ref), e2 = diff(stringRun(Scheme::GenAlpha, 0.01, 0.8, tEnd), ref);
    std::printf("\n  generalised-alpha rho=0.8: error %.3e (dt 0.02), %.3e (dt 0.01)\n", e1, e2);
    CHECK(e1 / e2 > 3.0);
    CHECK(e1 / e2 < 4.7);
}

TEST_CASE("GenAlpha: damped longitudinal mode decays at the exact rate with steps far above the explicit damping limit") {
    const double beta = 2e-4;                       // zeta = 0.031 for mode 1
    const std::vector<double> ts = {0.04, 0.08, 0.12};
    DynOptions rk; rk.implicitDamping = false;
    DynOptions ga; ga.scheme = Scheme::GenAlpha; ga.rhoInf = 1.0;
    double dtRk = 0, dtGa = 0;
    const std::vector<double> E0 = modeEnergy(0.0, rk, {1e-9}), Ex = modeEnergy(beta, rk, ts, &dtRk), Eg = modeEnergy(beta, ga, ts, &dtGa);
    const double omega = 2.0 * std::sqrt(1e5 / 0.1) / 0.25 * std::sin(kPi * 0.25 / 20.0);
    std::printf("\n  explicit dt %.2e s, GenAlpha dt %.2e s\n", dtRk, dtGa);
    CHECK(dtGa > 8.0 * dtRk);
    for (size_t k = 0; k < ts.size(); ++k) {
        const double exact = E0[0] * std::exp(-beta * omega * omega * ts[k]);
        std::printf("  t=%.2f: E/E0 explicit %.5f, GenAlpha %.5f, exact %.5f\n", ts[k], Ex[k] / E0[0], Eg[k] / E0[0], exact / E0[0]);
        CHECK(Eg[k] == doctest::Approx(exact).epsilon(0.04));
    }
}

TEST_CASE("GenAlpha: unconditionally stable and dissipative for a strongly damped mode at 25 times the explicit step") {
    DynOptions ga; ga.scheme = Scheme::GenAlpha; ga.rhoInf = 0.8;
    double dtGa = 0;
    const std::vector<double> E = modeEnergy(5e-3, ga, {0.01, 0.02, 0.05, 0.1, 0.2, 0.4}, &dtGa);
    double prev = 1e300;
    for (double e : E) { CHECK(std::isfinite(e)); CHECK(e <= prev + 1e-12); CHECK(e >= -1e-12); prev = e; }
    CHECK(E.back() < 1e-3 * E.front() + 1e-12);
}

TEST_CASE("GenAlpha: Lopez-Olocco clump-free chain (seabed, Morison, internal damping) agrees with the explicit scheme; result independent of the Jacobian interval") {
    json c = loadLopez("wocw");
    c["motion"]["radius_m"] = 0.225; c["motion"]["period_s"] = 2.8; c["motion"]["cycles"] = 12;
    c["numerics"]["implicit_damping"] = true;
    const double ref = lastCyclesMax(runCase(c));                // explicit RK4 + implicit damping (same as the explicit damping to 1e-4)
    json g = c; g["numerics"]["scheme"] = "genalpha"; g["numerics"]["implicit_dt_factor"] = 10.0;
    const CaseResult r10 = runCase(g);
    const double m10 = lastCyclesMax(r10);
    g["numerics"]["jacobian_interval"] = 1;
    const double m1 = lastCyclesMax(runCase(g));
    std::printf("\n  max tension: RK4 %.4f N, GenAlpha (dt %.1e s) %.4f N (%+.3f %%), Jacobian every step %.4f N\n", ref, r10.dt, m10, 100 * (m10 - ref) / ref, m1);
    CHECK(r10.finite);
    CHECK(r10.dt > 5.0 * 1.4e-5);
    CHECK(std::fabs(m10 - ref) / ref < 0.005);
    CHECK(std::fabs(m1 - m10) / m10 < 5e-4);
}
