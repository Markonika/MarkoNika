// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Irregular waves (JONSWAP / Pierson-Moskowitz spectra, frequency-dependent body load) and depth-varying currents.
#include <cmath>
#include <complex>
#include <cstdio>
#include "doctest.h"
#include "mooring/case_runner.hpp"
#include "mooring/coupled_runner.hpp"
#include "mooring/current.hpp"
#include "mooring/lumped_mass_cable.hpp"
#include "mooring/waves.hpp"

using namespace mooring;
using json = nlohmann::json;

namespace {
const double kPi = 3.14159265358979323846;
}

TEST_CASE("Spectrum: Pierson-Moskowitz closed form, exact m0, JONSWAP peak enhancement") {
    const double hs = 2.5, tp = 9.0, wp = 2 * kPi / tp;
    const JonswapSpectrum pm(hs, tp, 1.0);
    for (double r : {0.6, 0.8, 1.0, 1.3, 2.0}) {
        const double w = r * wp;
        const double exact = 5.0 / 16.0 * hs * hs * std::pow(wp, 4) * std::pow(w, -5) * std::exp(-1.25 * std::pow(wp / w, 4));
        CHECK(pm.density(w) == doctest::Approx(exact).scale(0).epsilon(1e-4));
    }
    // m0 = integral S = Hs^2 / 16 for both
    for (double gamma : {1.0, 3.3, 7.0}) {
        const JonswapSpectrum S(hs, tp, gamma);
        double m0 = 0; const double a = 0.1 * wp, b = 40 * wp; const int n = 200000; const double h = (b - a) / n;
        for (int i = 0; i < n; ++i) m0 += S.density(a + (i + 0.5) * h) * h;
        CHECK(m0 == doctest::Approx(hs * hs / 16.0).scale(0).epsilon(1e-5));
    }
    // the peak is enhanced and sharper for gamma = 3.3 (at the same Hs the peak density is higher than PM)
    CHECK(JonswapSpectrum(hs, tp, 3.3).density(wp) > 1.5 * pm.density(wp));
    CHECK_THROWS(JonswapSpectrum(-1.0, tp, 3.3));
}

TEST_CASE("Spectrum discretisation: exact Hs, reproducible seed, one component per bin, time-averaged variance, exact superposition of the kinematics") {
    SpectrumSpec sp; sp.hs = 3.0; sp.tp = 8.0; sp.gamma = 3.3; sp.n = 80; sp.seed = 42;
    WaveField a(50.0, 0.0, Stretching::Wheeler), b(50.0, 0.0, Stretching::Wheeler), c(50.0, 0.0, Stretching::Wheeler);
    a.addSpectrum(sp); b.addSpectrum(sp); sp.seed = 43; c.addSpectrum(sp);
    CHECK(a.significantHeight() == doctest::Approx(3.0).scale(0).epsilon(1e-12));
    REQUIRE(a.components().size() == 80);
    bool same = true, differs = false;
    for (size_t i = 0; i < 80; ++i) { same = same && a.components()[i].omega == b.components()[i].omega && a.components()[i].phase == b.components()[i].phase; differs = differs || a.components()[i].phase != c.components()[i].phase; }
    CHECK(same);
    CHECK(differs);
    const double wp = 2 * kPi / 8.0, dw = 3.5 * wp / 80;
    for (size_t i = 0; i < 80; ++i) {
        const double w = a.components()[i].omega;
        CHECK(w >= 0.5 * wp + i * dw);
        CHECK(w <= 0.5 * wp + (i + 1) * dw);
    }
    // surface elevation: sample variance over a long record ~ m0 = Hs^2/16 (random-phase process; the record has no repetition)
    double s2 = 0, s1 = 0; const int ns = 40000; const double dt = 0.5;
    for (int i = 0; i < ns; ++i) { const double e = a.elevation(0, 0, i * dt); s1 += e; s2 += e * e; }
    const double var = s2 / ns - (s1 / ns) * (s1 / ns);
    std::printf("\n  elevation variance %.4f m^2 (target m0 = %.4f), Hs from the record %.3f m\n", var, 9.0 / 16.0, 4 * std::sqrt(var));
    CHECK(var == doctest::Approx(9.0 / 16.0).scale(0).epsilon(0.06));
    // superposition: velocity and acceleration of the field = sum of single-component fields
    Vec3 u, acc, us, as;
    const Vec3 pos(3.0, -2.0, -4.0);
    a.kinematics(pos, 17.3, u, acc);
    for (const WaveComponent& w : a.components()) {
        WaveField one(50.0, 0.0, Stretching::Wheeler);
        one.addComponent(w.amplitude, 2 * kPi / w.omega, w.phase, 0.0);
        Vec3 uu, aa;
        one.kinematics(pos, 17.3, uu, aa);
        us += uu; as += aa;
    }
    // (Wheeler stretching depends on the total elevation, so the sum of separately stretched components differs slightly: the unstretched field is exact)
    WaveField ns1(50.0, 0.0, Stretching::None), ns2(50.0, 0.0, Stretching::None);
    ns1.addSpectrum(sp);
    Vec3 un, an, usum, asum;
    ns1.kinematics(pos, 17.3, un, an);
    for (const WaveComponent& w : ns1.components()) { WaveField one(50.0, 0.0, Stretching::None); one.addComponent(w.amplitude, 2 * kPi / w.omega, w.phase, 0.0); Vec3 uu, aa; one.kinematics(pos, 17.3, uu, aa); usum += uu; asum += aa; }
    CHECK(norm(un - usum) < 1e-12);
    CHECK(norm(an - asum) < 1e-12);
    CHECK(norm(u - us) < 0.1 * norm(u));            // stretched: close, not identical
    (void)ns2;
}

TEST_CASE("Current profile: table interpolation, clamping, power law, composition with wave kinematics") {
    CurrentProfile c(2.0);                                 // still-water level at z = 2
    c.addPoint(0.0, Vec3(1.0, 0, 0)); c.addPoint(10.0, Vec3(0.4, 0.2, 0)); c.addPoint(30.0, Vec3(0, 0, 0));
    CHECK(norm(c.velocity(5.0) - Vec3(1, 0, 0)) < 1e-14);                 // above the surface: surface value
    CHECK(norm(c.velocity(2.0) - Vec3(1, 0, 0)) < 1e-14);
    CHECK(norm(c.velocity(-3.0) - Vec3(0.7, 0.1, 0)) < 1e-14);            // depth 5: halfway between the first two points
    CHECK(norm(c.velocity(-8.0) - Vec3(0.4, 0.2, 0)) < 1e-14);            // depth 10
    CHECK(norm(c.velocity(-18.0) - Vec3(0.2, 0.1, 0)) < 1e-14);           // depth 20
    CHECK(norm(c.velocity(-100.0)) < 1e-14);                              // below the last point: last value (0)
    CHECK_THROWS(c.addPoint(30.0, Vec3()));
    CurrentProfile pw(0.0);
    pw.setPowerLaw(Vec3(0.8, 0, 0), 1.0 / 7.0, 100.0);
    CHECK(pw.velocity(0.0).x == doctest::Approx(0.8).scale(0).epsilon(1e-14));
    CHECK(pw.velocity(-50.0).x == doctest::Approx(0.8 * std::pow(0.5, 1.0 / 7.0)).scale(0).epsilon(1e-14));
    CHECK(norm(pw.velocity(-100.0)) < 1e-14);
    // JSON helper and composition with waves
    const json jenv = json::parse(R"({"water_surface_z_m": 2.0, "current_profile": {"type": "table", "points": [{"depth_m": 0, "velocity_m_s": [1, 0, 0]}, {"depth_m": 10, "velocity_m_s": [0.4, 0.2, 0]}]}})");
    const CurrentProfile cj = makeCurrentProfile(jenv, 2.0);
    CHECK(norm(cj.velocity(-3.0) - Vec3(0.7, 0.1, 0)) < 1e-14);
    WaveField w(50.0, 2.0, Stretching::Wheeler);
    w.addRegular(1.0, 8.0);
    const WaterField both = withCurrent(w.asWaterField(), cj);
    Vec3 u0, a0, u1, a1;
    w.kinematics(Vec3(1, 0, -3), 3.1, u0, a0);
    both(Vec3(1, 0, -3), 3.1, u1, a1);
    CHECK(norm(u1 - u0 - Vec3(0.7, 0.1, 0)) < 1e-14);
    CHECK(norm(a1 - a0) < 1e-14);
    const WaterField cur = withCurrent(WaterField(), cj);
    cur(Vec3(0, 0, -3), 0.0, u1, a1);
    CHECK(norm(u1 - Vec3(0.7, 0.1, 0)) < 1e-14);
    CHECK(norm(a1) < 1e-14);
}

TEST_CASE("Current profile on a line: Morison drag follows the local current at every node (vertical straight line in a sheared current)") {
    const int N = 20; const double L = 10.0, D = 0.1, Cdn = 1.2, rho = 1025.0;
    CableParams p; p.L = L; p.N = N; p.EA = 1e8; p.m_l = 1.0; p.w = 0.0; p.g = 0.0; p.D0 = D; p.D1 = D; p.Cdn = Cdn; p.Cdt = 0.0;
    LumpedMassCable c(p, Vec3(0, 0, -10.0), Vec3(0, 0, 0.0));
    Environment env; env.hydro = true; env.rho_w = rho; env.surfaceZ = 5.0;           // fully submerged
    CurrentProfile cur(5.0);
    cur.addPoint(0.0, Vec3(1.0, 0, 0)); cur.addPoint(8.0, Vec3(0.6, 0, 0)); cur.addPoint(15.0, Vec3(0.1, 0, 0));
    env.water = withCurrent(WaterField(), cur);
    c.setEnvironment(env);
    std::vector<Vec3> r(N + 1), v(N + 1);
    for (int i = 0; i <= N; ++i) r[i] = Vec3(0, 0, -10.0 + L * i / N);
    c.setInitialState(r, v);
    std::vector<Vec3> f;
    c.computeForces(r, v, 0.0, f);
    double worst = 0, ftot = 0, fexp = 0;
    for (int i = 0; i <= N; ++i) {
        const double Ls = (i == 0 || i == N ? 0.5 : 1.0) * L / N;
        const Vec3 U = cur.velocity(r[i].z);
        const Vec3 fe = U * (0.5 * Cdn * rho * D * norm(U) * Ls);
        worst = std::max(worst, norm(f[i] - fe));
        ftot += f[i].x; fexp += fe.x;
    }
    std::printf("\n  total drag %.4f N (analytic %.4f N), worst nodal difference %.2e N\n", ftot, fexp, worst);
    CHECK(worst < 1e-12);
    CHECK(ftot == doctest::Approx(fexp).scale(0).epsilon(1e-12));
    // against the uniform surface current the sheared profile gives a clearly smaller load
    CHECK(ftot < 0.8 * 0.5 * Cdn * rho * D * L * 1.0);
}

TEST_CASE("Irregular waves on a free body: heave response to a spectrum equals the sum of the steady-state responses of the components") {
    const double M = 10.0, A = 5.0, B = 12.0, C = 900.0;
    json cfg = json::parse(R"({
      "body": {"mass_kg": 10.0, "inertia_diag_kg_m2": [1,1,1], "cg_ref_m": [0,0,0], "A_diag": [0,0,5.0,0,0,0], "B_diag": [0,0,12.0,0,0,0], "C_diag": [0,0,900.0,0,0,0],
               "Dq": [0,0,0,0,0,0], "F0": [0,0,0,0,0,0],
               "wave_force_table": [ {"period_s": 0.4, "w": [0,0,60,0,0,0], "delta": [0,0,0.1,0,0,0]}, {"period_s": 1.0, "w": [0,0,140,0,0,0], "delta": [0,0,0.5,0,0,0]},
                                     {"period_s": 2.0, "w": [0,0,220,0,0,0], "delta": [0,0,-0.2,0,0,0]}, {"period_s": 3.6, "w": [0,0,300,0,0,0], "delta": [0,0,-0.6,0,0,0]} ]},
      "lines": [],
      "waves": {"depth_m": 30.0, "surface_z_m": 0.0, "spectrum": {"type": "jonswap", "hs_m": 0.1, "tp_s": 1.7, "gamma": 3.3, "n_components": 30, "seed": 5}},
      "numerics": {"dt_body_s": 0.002, "t_end_s": 60.0},
      "output": {"write": false, "dt_out_s": 0.05}
    })");
    const CoupledResult res = runCoupledCase(cfg);
    REQUIRE(res.finite);
    // analytic steady state
    WaveField wf(30.0, 0.0);
    SpectrumSpec sp; sp.hs = 0.1; sp.tp = 1.7; sp.gamma = 3.3; sp.n = 30; sp.seed = 5;
    wf.addSpectrum(sp);
    struct Row { double T, w, d; };
    const Row rows[4] = {{0.4, 60, 0.1}, {1.0, 140, 0.5}, {2.0, 220, -0.2}, {3.6, 300, -0.6}};
    auto interp = [&](double om, double* w, double* d) {
        const double oms[4] = {2 * kPi / 3.6, 2 * kPi / 2.0, 2 * kPi / 1.0, 2 * kPi / 0.4};
        const Row* rr[4] = {&rows[3], &rows[2], &rows[1], &rows[0]};
        int j = 0; while (j + 2 < 4 && om > oms[j + 1]) ++j;
        const double f = std::min(1.0, std::max(0.0, (om - oms[j]) / (oms[j + 1] - oms[j])));
        *w = rr[j]->w * (1 - f) + rr[j + 1]->w * f; *d = rr[j]->d * (1 - f) + rr[j + 1]->d * f;
    };
    double worst = 0, amp = 0;
    for (const CoupledSample& s : res.samples) {
        if (s.t < 40.0) continue;
        double z = 0;
        for (const WaveComponent& c : wf.components()) {
            double w, d; interp(c.omega, &w, &d);
            const std::complex<double> den(C - c.omega * c.omega * (M + A), c.omega * B);
            z += std::imag(w * c.amplitude * std::exp(std::complex<double>(0, c.omega * s.t + c.phase + d)) / den);
        }
        worst = std::max(worst, std::fabs(s.xi[2] - z)); amp = std::max(amp, std::fabs(z));
    }
    std::printf("\n  heave: max |simulated - analytic| = %.3e m for a response amplitude of %.3e m (%zu samples)\n", worst, amp, res.samples.size());
    CHECK(amp > 1e-3);
    CHECK(worst < 2e-3 * amp);
    // a component outside the table range, or a spectrum with the single-wave input, is rejected
    json bad = cfg; bad["waves"]["spectrum"]["tp_s"] = 6.0;                          // components down to ~2.4 s * ... and above the 3.6 s limit
    CHECK_THROWS(runCoupledCase(bad));
    json bad2 = cfg; bad2["body"].erase("wave_force_table"); bad2["body"]["wave_force"] = {{"w", {0,0,1,0,0,0}}, {"delta", {0,0,0,0,0,0}}};
    CHECK_THROWS(runCoupledCase(bad2));
}
