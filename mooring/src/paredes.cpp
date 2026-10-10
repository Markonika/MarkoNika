// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#include "mooring/paredes.hpp"
#include <cmath>
#include <complex>
#include <stdexcept>
#include "mooring/lumped_mass_cable.hpp"
#include "mooring/waves.hpp"

namespace mooring {
using json = nlohmann::json;
namespace {

Vec6 lineLoad(CoupledSystem& sys) {
    Vec6 F{};
    for (int k = 0; k < sys.numLines(); ++k) {
        const Vec3 Ra = sys.fairleadPosition(k) - sys.body().cg();
        RigidBody6DOF::addPointLoad(Ra, sys.lineForce(k), F);
    }
    return F;
}

// Solve with DOF 'fixedDof' held at 'value' and read the line loads.
Vec6 heldEquilibrium(CoupledSystem& sys, int fixedDof, double value, bool* ok) {
    std::array<bool, 6> fixed{};
    fixed[fixedDof] = true;
    sys.body().xi[fixedDof] = value;
    const EquilibriumResult r = sys.solveEquilibrium(1e-6, 40, fixed);
    if (ok) *ok = r.converged;
    sys.begin(0.0);                                   // lines at their static state, forces evaluated
    return lineLoad(sys);
}

}  // namespace

ParedesStatics paredesStatics(const json& cfg) {
    auto sysPtr = buildCoupledSystem(cfg);
    CoupledSystem& sys = *sysPtr;
    ParedesStatics s;
    s.eq = sys.solveEquilibrium(1e-6);
    sys.begin(0.0);
    s.xi = sys.body().xi;
    s.draftChange = -s.xi[2];
    for (int k = 0; k < sys.numLines(); ++k) {
        const auto* lm = dynamic_cast<const LumpedMassCable*>(&sys.line(k));
        s.topTension.push_back(lm ? lm->endSegmentTension(true) : 0.0);
        s.fairleadForce.push_back(norm(sys.lineForce(k)));
        s.netLineForce += sys.lineForce(k);
    }
    return s;
}

ParedesStiffness paredesSurgeStiffness(const json& cfg, const std::vector<double>& xs) {
    auto sysPtr = buildCoupledSystem(cfg);
    CoupledSystem& sys = *sysPtr;
    ParedesStiffness st;
    sys.solveEquilibrium(1e-6);
    sys.begin(0.0);
    st.F0 = lineLoad(sys)[0];
    const Vec6 xi0 = sys.body().xi;
    for (double x : xs) {
        sys.body().xi = xi0;
        bool ok = false;
        const Vec6 F = heldEquilibrium(sys, 0, x, &ok);
        st.x.push_back(x); st.Fx.push_back(F[0]); st.K.push_back(-(F[0] - st.F0) / x); st.converged.push_back(ok);
    }
    return st;
}

double paredesHeaveStiffness(const json& cfg, double delta) {
    auto sysPtr = buildCoupledSystem(cfg);
    CoupledSystem& sys = *sysPtr;
    sys.solveEquilibrium(1e-6);
    const Vec6 xi0 = sys.body().xi;
    sys.body().xi = xi0;
    const double Fp = heldEquilibrium(sys, 2, xi0[2] + delta, nullptr)[2];
    sys.body().xi = xi0;
    const double Fm = heldEquilibrium(sys, 2, xi0[2] - delta, nullptr)[2];
    return -(Fp - Fm) / (2.0 * delta);
}

double paredesPitchStiffness(const json& cfg, double d) {
    auto sysPtr = buildCoupledSystem(cfg);
    CoupledSystem& sys = *sysPtr;
    sys.solveEquilibrium(1e-6);
    const Vec6 xi0 = sys.body().xi;
    sys.body().xi = xi0;
    const double Mp = heldEquilibrium(sys, 4, xi0[4] + d, nullptr)[4];
    sys.body().xi = xi0;
    const double Mm = heldEquilibrium(sys, 4, xi0[4] - d, nullptr)[4];
    return -(Mp - Mm) / (2.0 * d);
}

double paredesDecayPeriod(json cfg, int dof, double offset, double tEnd, double dtBody) {
    cfg["numerics"]["t_end_s"] = tEnd;
    cfg["numerics"]["dt_body_s"] = dtBody;
    cfg["output"]["write"] = false;
    cfg["output"]["dt_out_s"] = std::max(dtBody, 0.01);
    if (cfg["initial"].value("equilibrium", false)) {
        Vec6 off{}; off[dof] = offset;
        cfg["initial"]["xi_offset"] = std::vector<double>(off.begin(), off.end());
    } else {
        std::vector<double> x0(6, 0.0); x0[dof] = offset;
        cfg["initial"]["xi0"] = x0;
    }
    const CoupledResult r = runCoupledCase(cfg);
    if (!r.finite) return 0.0;
    // Reference: mean of the last 20 % of the record (equilibrium after the decay), crossings counted about it.
    double ref = 0; int n = 0;
    for (size_t i = r.samples.size() * 4 / 5; i < r.samples.size(); ++i) { ref += r.samples[i].xi[dof]; ++n; }
    ref = n ? ref / n : 0.0;
    std::vector<double> up, down;
    for (size_t i = 1; i < r.samples.size(); ++i) {
        const double a = r.samples[i - 1].xi[dof] - ref, b = r.samples[i].xi[dof] - ref;
        if (a < 0 && b >= 0) up.push_back(r.samples[i - 1].t + (r.samples[i].t - r.samples[i - 1].t) * (-a) / (b - a));
        if (a > 0 && b <= 0) down.push_back(r.samples[i - 1].t + (r.samples[i].t - r.samples[i - 1].t) * a / (a - b));
    }
    auto per = [](const std::vector<double>& c) { return c.size() >= 3 ? (c.back() - c.front()) / (c.size() - 1) : 0.0; };
    const double pu = per(up), pd = per(down);
    if (pu > 0 && pd > 0) return 0.5 * (pu + pd);
    return pu > 0 ? pu : pd;
}

namespace {
const double kPi = 3.14159265358979323846;
// Amplitude of the component at omega in x(t) over [t0, t1]: least squares for a sin + b cos + c + d t.
double fitAmp(const std::vector<double>& t, const std::vector<double>& x, double om, double t0, double t1, double* mean) {
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
    if (mean) *mean = s[2];
    return std::hypot(s[0], s[1]);
}
struct Tab34 { double T, A11, B11, w1, d1, A33, B33, w3, d3, A55, B55, w5, d5, A15, B15; };
// Paredes (2016), Table 3.4 (hull alone, 0.9 m water depth)
const Tab34 kTab[] = {{1.30, 24.55, 31.73, 1057.8, 1.475, 26.23, 37.94, 817.9, 0.258, 0.2599, 0.03336, 34.30, -1.667, 1.1739, 1.0288},
                      {1.40, 23.921, 20.82, 963.8, 1.501, 27.05, 38.09, 921.8, 0.206, 0.2581, 0.01974, 29.68, -1.641, 1.1369, 0.6410},
                      {1.50, 23.17, 14.05, 880.2, 1.519, 27.83, 37.47, 1016.4, 0.168, 0.2564, 0.01213, 25.86, -1.623, 1.1004, 0.4129}};
}  // namespace

ParedesRAO paredesWaveRAO(json cfg, double period, double height, double tEnd, int fitCycles, double dtBody, bool steadyStart, bool instantaneousSurface) {
    const Tab34* c = nullptr;
    for (const Tab34& e : kTab) if (std::fabs(e.T - period) < 1e-9) c = &e;
    if (!c) throw std::invalid_argument("Table 3.4 has coefficients only for T = 1.30, 1.40, 1.50 s");
    json A = json::array(), B = json::array();
    for (int i = 0; i < 6; ++i) { A.push_back(json::array({0, 0, 0, 0, 0, 0})); B.push_back(json::array({0, 0, 0, 0, 0, 0})); }
    A[0][0] = A[1][1] = c->A11; A[2][2] = c->A33; A[3][3] = A[4][4] = c->A55; A[0][4] = A[4][0] = c->A15;
    B[0][0] = B[1][1] = c->B11; B[2][2] = c->B33; B[3][3] = B[4][4] = c->B55; B[0][4] = B[4][0] = c->B15;
    cfg["body"].erase("A_diag"); cfg["body"].erase("B_diag");
    cfg["body"]["A"] = A; cfg["body"]["B"] = B;
    cfg["body"]["wave_force"] = {{"w", {c->w1, 0.0, c->w3, 0.0, c->w5, 0.0}}, {"delta", {c->d1, 0.0, c->d3, 0.0, c->d5, 0.0}}};
    const double depth = 0.9;
    const double ramp = steadyStart ? 0.0 : 5.0 * period;
    cfg["waves"] = {{"depth_m", depth}, {"surface_z_m", cfg["environment"].value("water_surface_z_m", depth)}, {"ramp_time_s", ramp}, {"instantaneous_surface", instantaneousSurface},
                    {"components", json::array({ {{"height_m", height}, {"period_s", period}} })}};
    if (steadyStart) {
        // Start from the steady-state response of the FREE body (closed form of the 6-DOF model without lines), so that only the small
        // mooring-induced transient remains: the lightly damped surge-pitch mode would otherwise need ~1000 s to settle (see test_waves.cpp).
        const double om = 2.0 * kPi / period, a = 0.5 * height, M = cfg["body"].at("mass_kg");
        const double I55 = cfg["body"].at("inertia_diag_kg_m2").at(1);
        const double C33 = cfg["body"].at("C_diag").at(2), C55 = cfg["body"].at("C_diag").at(4);
        const std::complex<double> im(0, 1);
        const auto z3 = -om * om * (M + c->A33) + im * om * c->B33 + C33;
        const auto X3 = c->w3 * a * std::exp(im * c->d3) / z3;
        const auto z11 = -om * om * (M + c->A11) + im * om * c->B11, z15 = -om * om * c->A15 + im * om * c->B15;
        const auto z55 = -om * om * (I55 + c->A55) + im * om * c->B55 + C55;
        const auto f1 = c->w1 * a * std::exp(im * c->d1), f5 = c->w5 * a * std::exp(im * c->d5);
        const auto det = z11 * z55 - z15 * z15;
        const auto X1 = (f1 * z55 - z15 * f5) / det, X5 = (z11 * f5 - z15 * f1) / det;
        cfg["initial"]["xi_offset"] = {X1.imag(), 0.0, X3.imag(), 0.0, X5.imag(), 0.0};
        cfg["initial"]["xi_dot0"] = {om * X1.real(), 0.0, om * X3.real(), 0.0, om * X5.real(), 0.0};
    }
    cfg["numerics"]["t_end_s"] = tEnd; cfg["numerics"]["dt_body_s"] = dtBody;
    cfg["output"]["write"] = false; cfg["output"]["dt_out_s"] = 0.01;
    const CoupledResult r = runCoupledCase(cfg);
    ParedesRAO o; o.finite = r.finite;
    const double om = 2.0 * kPi / period;
    o.a = 0.5 * height; o.k = WaveField::wavenumber(om, depth);
    std::vector<double> t, x1, x3, x5;
    for (const auto& s : r.samples) { t.push_back(s.t); x1.push_back(s.xi[0]); x3.push_back(s.xi[2]); x5.push_back(s.xi[4]); }
    const double t0 = tEnd - fitCycles * period;
    o.surge = fitAmp(t, x1, om, t0, tEnd, &o.meanSurge) / o.a;
    o.heave = fitAmp(t, x3, om, t0, tEnd, nullptr) / o.a;
    o.pitch = fitAmp(t, x5, om, t0, tEnd, nullptr) / (o.k * o.a);
    return o;
}

}  // namespace mooring
