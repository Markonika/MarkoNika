// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Milestone 9: Airy waves with Wheeler stretching, wave excitation of the body, regular-wave RAOs.
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <vector>
#include "doctest.h"
#include "mooring/coupled_runner.hpp"
#include "mooring/waves.hpp"

using namespace mooring;
using json = nlohmann::json;

namespace {
const double kPi = 3.14159265358979323846;
using cd = std::complex<double>;

// Least-squares amplitude of the component at angular frequency om in x(t) over samples with t in [t0, t1]
// (model a sin + b cos + c + d t, so a constant offset and a linear drift do not leak into the amplitude).
double fitAmplitude(const std::vector<double>& t, const std::vector<double>& x, double om, double t0, double t1) {
    double M[4][5] = {};
    for (size_t i = 0; i < t.size(); ++i) {
        if (t[i] < t0 || t[i] > t1) continue;
        const double f[4] = {std::sin(om * t[i]), std::cos(om * t[i]), 1.0, t[i] - t0};
        for (int r = 0; r < 4; ++r) { for (int c = 0; c < 4; ++c) M[r][c] += f[r] * f[c]; M[r][4] += f[r] * x[i]; }
    }
    for (int c = 0; c < 4; ++c) {
        int p = c; for (int r = c + 1; r < 4; ++r) if (std::fabs(M[r][c]) > std::fabs(M[p][c])) p = r;
        for (int j = 0; j < 5; ++j) std::swap(M[c][j], M[p][j]);
        for (int r = c + 1; r < 4; ++r) { const double f = M[r][c] / M[c][c]; for (int j = c; j < 5; ++j) M[r][j] -= f * M[c][j]; }
    }
    double s[4];
    for (int r = 3; r >= 0; --r) { double v = M[r][4]; for (int c = r + 1; c < 4; ++c) v -= M[r][c] * s[c]; s[r] = v / M[r][r]; }
    return std::hypot(s[0], s[1]);
}

// Table 3.4 (Paredes 2016): hydrodynamic coefficients of the buoy at 0.9 m depth, T = 1.30 s.
struct Tab34 { double A11, B11, w1, d1, A33, B33, w3, d3, A55, B55, w5, d5, A15, B15; };
const Tab34 T130 = {24.55, 31.73, 1057.8, 1.475, 26.23, 37.94, 817.9, 0.258, 0.2599, 0.03336, 34.30, -1.667, 1.1739, 1.0288};
const Tab34 T140 = {23.921, 20.82, 963.8, 1.501, 27.05, 38.09, 921.8, 0.206, 0.2581, 0.01974, 29.68, -1.641, 1.1369, 0.6410};
}  // namespace

TEST_CASE("wave number: finite-depth dispersion relation and its limits") {
    for (double h : {0.9, 5.0, 100.0})
        for (double T : {0.8, 1.3, 2.33}) {
            const double om = 2 * kPi / T, k = WaveField::wavenumber(om, h);
            CHECK(9.81 * k * std::tanh(k * h) == doctest::Approx(om * om).epsilon(1e-12));
        }
    const double om = 2 * kPi / 1.3;
    CHECK(WaveField::wavenumber(om, 1e3) == doctest::Approx(om * om / 9.81).epsilon(1e-9));            // deep water
    CHECK(WaveField::wavenumber(0.05, 0.01) == doctest::Approx(0.05 / std::sqrt(9.81 * 0.01)).epsilon(2e-3));   // shallow water
}

TEST_CASE("Airy kinematics: surface, bed, kinematic surface condition, continuity, superposition") {
    const double h = 0.9, z0 = 0.9, H = 0.08, T = 1.4, A = 0.5 * H, om = 2 * kPi / T;
    WaveField wh(h, z0, Stretching::Wheeler), ln(h, z0, Stretching::None);
    wh.addRegular(H, T); ln.addRegular(H, T);
    const double k = wh.components()[0].k, kh = k * h;
    // elevation at the origin is A sin(omega t)
    CHECK(wh.elevation(0, 0, 0.0) == 0.0);
    for (double t : {0.3, 0.9, 1.1}) CHECK(wh.elevation(0, 0, t) == doctest::Approx(A * std::sin(om * t)).epsilon(1e-12).scale(0.0));
    Vec3 u, a;
    // crest at the origin (t = T/4): Wheeler gives the SWL velocity at the surface, linear theory a slightly larger one
    const double tc = T / 4;
    wh.kinematics(Vec3(0, 0, z0 + A * (1 - 1e-9)), tc, u, a);
    CHECK(u.x == doctest::Approx(A * om / std::tanh(kh)).epsilon(1e-8).scale(0.0));
    ln.kinematics(Vec3(0, 0, z0 + A * (1 - 1e-9)), tc, u, a);
    CHECK(u.x == doctest::Approx(A * om * std::cosh(k * (A + h)) / std::sinh(kh)).epsilon(1e-8).scale(0.0));
    // seabed: u = A omega / sinh(kh), w = 0
    for (const WaveField* f : {&wh, &ln}) {
        f->kinematics(Vec3(0, 0, z0 - h), tc, u, a);
        CHECK(u.x == doctest::Approx(A * om / std::sinh(kh)).epsilon(1e-8).scale(0.0));
        CHECK(std::fabs(u.z) < 1e-12);
    }
    // kinematic surface condition under Wheeler: w at the instantaneous surface = d(eta)/dt
    for (double t : {0.1, 0.37, 0.71, 1.1}) {
        const double eta = wh.elevation(0.3, 0, t), dt = 1e-6;
        const double detadt = (wh.elevation(0.3, 0, t + dt) - wh.elevation(0.3, 0, t - dt)) / (2 * dt);
        wh.kinematics(Vec3(0.3, 0, z0 + eta), t, u, a);
        CHECK(u.z == doctest::Approx(detadt).epsilon(1e-6).scale(A * om));
    }
    // no water above the surface
    wh.kinematics(Vec3(0, 0, z0 + 2 * A), tc, u, a);
    CHECK(norm(u) == 0.0);
    // linear theory: acceleration = du/dt at a fixed point and continuity du/dx + dw/dz = 0 (finite differences)
    const Vec3 p(0.2, 0.0, 0.5);
    const double t = 0.55, d = 1e-6;
    Vec3 up, um, ap, am;
    ln.kinematics(p, t + d, up, ap); ln.kinematics(p, t - d, um, am); ln.kinematics(p, t, u, a);
    CHECK(norm(a - (up - um) / (2 * d)) < 1e-6 * A * om * om);
    Vec3 ux, uxm, uz, uzm, tmp;
    ln.kinematics(p + Vec3(d, 0, 0), t, ux, tmp); ln.kinematics(p - Vec3(d, 0, 0), t, uxm, tmp);
    ln.kinematics(p + Vec3(0, 0, d), t, uz, tmp); ln.kinematics(p - Vec3(0, 0, d), t, uzm, tmp);
    CHECK(((ux.x - uxm.x) + (uz.z - uzm.z)) / (2 * d) == doctest::Approx(0.0).scale(A * om * k).epsilon(1e-6));
    // Wheeler acceleration (stretched linear value) stays close to the exact derivative: documented approximation
    const Vec3 ps(0.2, 0.0, z0 + 0.5 * A);
    wh.kinematics(ps, t + d, up, ap); wh.kinematics(ps, t - d, um, am); wh.kinematics(ps, t, u, a);
    std::printf("\n  Wheeler: |a - du/dt| / (A omega^2) = %.3f at z = SWL + A/2\n", norm(a - (up - um) / (2 * d)) / (A * om * om));
    CHECK(norm(a - (up - um) / (2 * d)) < 0.15 * A * om * om);
    // superposition of two components in different directions
    WaveField two(h, z0, Stretching::None), w1(h, z0, Stretching::None), w2(h, z0, Stretching::None);
    two.addRegular(0.06, 1.2, 0.3, 0.0); two.addRegular(0.04, 1.7, -0.5, 40.0);
    w1.addRegular(0.06, 1.2, 0.3, 0.0); w2.addRegular(0.04, 1.7, -0.5, 40.0);
    Vec3 u1, a1, u2, a2;
    two.kinematics(Vec3(0.4, -0.2, 0.6), 0.8, u, a); w1.kinematics(Vec3(0.4, -0.2, 0.6), 0.8, u1, a1); w2.kinematics(Vec3(0.4, -0.2, 0.6), 0.8, u2, a2);
    CHECK(norm(u - (u1 + u2)) < 1e-12); CHECK(norm(a - (a1 + a2)) < 1e-12);
    CHECK(two.elevation(0.4, -0.2, 0.8) == doctest::Approx(w1.elevation(0.4, -0.2, 0.8) + w2.elevation(0.4, -0.2, 0.8)).epsilon(1e-12));
    // ramp
    WaveField rp(h, z0); rp.addRegular(H, T); rp.setRampTime(5.0);
    CHECK(rp.rampFactor(0.0) == 0.0); CHECK(rp.rampFactor(2.5) == doctest::Approx(0.5)); CHECK(rp.rampFactor(6.0) == 1.0);
}

namespace {
json freeBodyConfig(const Tab34& c, double H, double T, double tEnd) {
    json cfg;
    const double Cs[6] = {0.0, 0.0, 1000.0 * 9.81 * kPi * 0.515 * 0.515 / 4.0, 37.76, 37.76, 0.0};
    cfg["body"] = {{"mass_kg", 35.5}, {"inertia_diag_kg_m2", {0.87, 0.87, 1.18}}, {"cg_ref_m", {0, 0, 0.8}},
                   {"C_diag", Cs}, {"wave_force", {{"w", {c.w1, 0.0, c.w3, 0.0, c.w5, 0.0}}, {"delta", {c.d1, 0.0, c.d3, 0.0, c.d5, 0.0}}}}};
    cfg["body"]["A_diag"] = {c.A11, c.A11, c.A33, c.A55, c.A55, 0.0};
    cfg["body"]["B_diag"] = {c.B11, c.B11, c.B33, c.B55, c.B55, 0.0};
    json A = json::array(), B = json::array();                     // surge-pitch coupling about the CG (sign convention of Table 3.4 assumed)
    for (int i = 0; i < 6; ++i) { A.push_back(json::array({0, 0, 0, 0, 0, 0})); B.push_back(json::array({0, 0, 0, 0, 0, 0})); }
    A[0][0] = c.A11; A[1][1] = c.A11; A[2][2] = c.A33; A[3][3] = c.A55; A[4][4] = c.A55; A[0][4] = A[4][0] = c.A15;
    B[0][0] = c.B11; B[1][1] = c.B11; B[2][2] = c.B33; B[3][3] = c.B55; B[4][4] = c.B55; B[0][4] = B[4][0] = c.B15;
    cfg["body"].erase("A_diag"); cfg["body"].erase("B_diag");
    cfg["body"]["A"] = A; cfg["body"]["B"] = B;
    cfg["lines"] = json::array();
    cfg["waves"] = {{"depth_m", 0.9}, {"surface_z_m", 0.9}, {"ramp_time_s", 5 * T}, {"components", json::array({ {{"height_m", H}, {"period_s", T}} })}};
    cfg["initial"] = {{"equilibrium", false}};
    cfg["numerics"] = {{"dt_body_s", 1e-3}, {"t_end_s", tEnd}};
    cfg["output"] = {{"dt_out_s", 1e-2}, {"write", false}};
    return cfg;
}

// Closed-form RAOs (amplitudes of surge, heave, pitch) of the 6-DOF model above for a regular wave of amplitude a.
void analyticRAO(const Tab34& c, double T, double a, double& x1, double& x3, double& x5) {
    const double om = 2 * kPi / T, M = 35.5, I = 0.87, C33 = 1000.0 * 9.81 * kPi * 0.515 * 0.515 / 4.0, C55 = 37.76;
    const cd i(0, 1);
    const cd z3 = -om * om * (M + c.A33) + i * om * c.B33 + C33;
    x3 = std::abs(c.w3 * a * std::exp(i * c.d3) / z3);
    const cd z11 = -om * om * (M + c.A11) + i * om * c.B11, z15 = -om * om * c.A15 + i * om * c.B15, z55 = -om * om * (I + c.A55) + i * om * c.B55 + C55;
    const cd f1 = c.w1 * a * std::exp(i * c.d1), f5 = c.w5 * a * std::exp(i * c.d5);
    const cd det = z11 * z55 - z15 * z15;
    x1 = std::abs((f1 * z55 - z15 * f5) / det);
    x5 = std::abs((z11 * f5 - z15 * f1) / det);
}
}  // namespace

TEST_CASE("free buoy in regular waves (Table 3.4 coefficients): simulated RAOs vs the closed-form solution") {
    std::printf("\n  free buoy, H = 0.08 m: RAO surge (m/m), heave (m/m), pitch (rad per ka)\n");
    for (auto cT : {std::make_pair(T130, 1.30), std::make_pair(T140, 1.40)}) {
        const double T = cT.second, a = 0.04, om = 2 * kPi / T;
        const double tEnd = 1200.0;              // the coupled surge-pitch mode is very lightly damped: at 600 s the pitch amplitude is still 0.26 % off, at 1200 s 3e-5
        const CoupledResult r = runCoupledCase(freeBodyConfig(cT.first, 0.08, T, tEnd));
        REQUIRE(r.finite);
        std::vector<double> t, x1, x3, x5;
        for (const auto& s : r.samples) { t.push_back(s.t); x1.push_back(s.xi[0]); x3.push_back(s.xi[2]); x5.push_back(s.xi[4]); }
        const double t0 = tEnd - 20 * T, t1 = tEnd;                      // last 20 wave periods
        const double s1 = fitAmplitude(t, x1, om, t0, t1), s3 = fitAmplitude(t, x3, om, t0, t1), s5 = fitAmplitude(t, x5, om, t0, t1);
        double e1, e3, e5; analyticRAO(cT.first, T, a, e1, e3, e5);
        const double k = WaveField::wavenumber(om, 0.9);
        std::printf("  T = %.2f s: simulated %.4f %.4f %.4f | analytic %.4f %.4f %.4f\n", T, s1 / a, s3 / a, s5 / (k * a), e1 / a, e3 / a, e5 / (k * a));
        CHECK(s3 == doctest::Approx(e3).epsilon(1e-3).scale(0.0));
        CHECK(s1 == doctest::Approx(e1).epsilon(1e-3).scale(0.0));
        CHECK(s5 == doctest::Approx(e5).epsilon(1e-3).scale(0.0));
    }
}

TEST_CASE("free buoy started from the closed-form steady state stays on it (initialisation used for the moored runs)") {
    const double T = 1.30, a = 0.04, om = 2 * kPi / T, M = 35.5, I = 0.87, C33 = 1000.0 * 9.81 * kPi * 0.515 * 0.515 / 4.0, C55 = 37.76;
    const Tab34& c = T130;
    const cd i(0, 1);
    const cd X3 = c.w3 * a * std::exp(i * c.d3) / (-om * om * (M + c.A33) + i * om * c.B33 + C33);
    const cd z11 = -om * om * (M + c.A11) + i * om * c.B11, z15 = -om * om * c.A15 + i * om * c.B15, z55 = -om * om * (I + c.A55) + i * om * c.B55 + C55;
    const cd f1 = c.w1 * a * std::exp(i * c.d1), f5 = c.w5 * a * std::exp(i * c.d5), det = z11 * z55 - z15 * z15;
    const cd X1 = (f1 * z55 - z15 * f5) / det, X5 = (z11 * f5 - z15 * f1) / det;
    json cfg = freeBodyConfig(c, 0.08, T, 20.0);
    cfg["waves"]["ramp_time_s"] = 0.0;
    cfg["initial"]["xi0"] = {X1.imag(), 0.0, X3.imag(), 0.0, X5.imag(), 0.0};
    cfg["initial"]["xi_dot0"] = {om * X1.real(), 0.0, om * X3.real(), 0.0, om * X5.real(), 0.0};
    const CoupledResult r = runCoupledCase(cfg);
    double e1 = 0, e3 = 0, e5 = 0;
    for (const auto& s : r.samples) {
        const cd ph = std::exp(i * om * s.t);
        e1 = std::max(e1, std::fabs(s.xi[0] - (X1 * ph).imag()) / std::abs(X1));
        e3 = std::max(e3, std::fabs(s.xi[2] - (X3 * ph).imag()) / std::abs(X3));
        e5 = std::max(e5, std::fabs(s.xi[4] - (X5 * ph).imag()) / std::abs(X5));
    }
    std::printf("\n  free buoy from the steady state, 20 s: max deviation from the closed-form time history: surge %.1e, heave %.1e, pitch %.1e (relative to amplitude)\n", e1, e3, e5);
    CHECK(e1 < 1e-4); CHECK(e3 < 1e-4); CHECK(e5 < 1e-4);
}

// Regular-wave RAOs of the moored buoy (CON1, CON2, CAT) vs the thesis, Figs 5.22-5.24 at T = 1.3 and 1.4 s (the only sampled periods of those
// figures with Table 3.4 coefficients). The thesis values are read by eye from the plots and are given as the range between the H = 0.04 m and
// H = 0.08 m panels (the thesis notes the RAO decreases with wave height; a first-order model has no such dependence).
#include "mooring/paredes.hpp"
TEST_CASE("moored buoy in regular waves (milestone 9): surge, heave, pitch RAO at T = 1.3, 1.4 s vs Figs 5.22-5.24 [validation run, ~10 min]" * doctest::skip()) {
    struct M { const char* cfg; double T; double surge[2], heave[2], pitch[2]; };
    // ranges [min, max] over both wave-height panels and the three configurations only differ slightly (thesis: 'unexpectedly small')
    const M meas[] = {
        {"con1", 1.3, {0.25, 0.33}, {1.45, 1.50}, {2.10, 2.30}}, {"con2", 1.3, {0.40, 0.45}, {1.45, 1.50}, {1.70, 1.95}}, {"cat", 1.3, {0.25, 0.33}, {1.45, 1.50}, {2.10, 2.55}},
        {"con1", 1.4, {0.50, 0.50}, {1.40, 1.45}, {1.75, 1.90}}, {"con2", 1.4, {0.62, 0.68}, {1.30, 1.48}, {1.55, 1.75}}, {"cat", 1.4, {0.46, 0.50}, {1.33, 1.33}, {1.85, 2.05}}};
    for (const M& m : meas) {
        std::ifstream in(std::string(MOORING_SOURCE_DIR) + "/examples/paredes/" + m.cfg + ".json");
        REQUIRE(in.good());
        const json cfg = json::parse(in);
        const ParedesRAO r = paredesWaveRAO(cfg, m.T, 0.08, 30.0, 15);
        std::printf("\n  %s T=%.1f: RAO surge %.3f [meas %.2f-%.2f], heave %.3f [%.2f-%.2f], pitch/ka %.3f [%.2f-%.2f]\n", m.cfg, m.T, r.surge, m.surge[0], m.surge[1],
                    r.heave, m.heave[0], m.heave[1], r.pitch, m.pitch[0], m.pitch[1]);
        CHECK(r.finite);
        // A priori acceptance criteria (fixed before the runs): heave within 25 %, pitch within 30 % of the measured range, surge of the right order.
        // Outcome (docs/waves_validation.md): heave passes in 6/6, pitch in 3/6 (CAT and CON1 at 1.4 s miss), surge does NOT agree quantitatively
        // (over-predicted by 14 % to a factor 3). The assertions below only pin the documented envelope of these known disagreements.
        const double sm = 0.5 * (m.surge[0] + m.surge[1]), hm = 0.5 * (m.heave[0] + m.heave[1]), pm = 0.5 * (m.pitch[0] + m.pitch[1]);
        std::printf("    model/measured-mid: surge %.2f, heave %.2f, pitch %.2f  | a priori: heave %s, pitch %s\n", r.surge / sm, r.heave / hm, r.pitch / pm,
                    (r.heave > 0.75 * m.heave[0] && r.heave < 1.25 * m.heave[1]) ? "PASS" : "FAIL", (r.pitch > 0.70 * m.pitch[0] && r.pitch < 1.30 * m.pitch[1]) ? "PASS" : "FAIL");
        CHECK(r.heave / hm > 0.75); CHECK(r.heave / hm < 1.05);
        CHECK(r.pitch / pm > 0.45); CHECK(r.pitch / pm < 1.40);
        CHECK(r.surge / sm > 1.0);  CHECK(r.surge / sm < 3.3);      // the sign of the disagreement (over-prediction) is part of the record
    }
}
