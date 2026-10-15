// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Frequency-dependent added mass and radiation damping (Cummins equation with a retardation function).
#include <cmath>
#include <complex>
#include <cstdio>
#include "doctest.h"
#include "mooring/coupled_runner.hpp"
#include "mooring/radiation.hpp"
#include "mooring/waves.hpp"

using namespace mooring;
using json = nlohmann::json;

namespace {
const double kPi = 3.14159265358979323846;

// Heave data with an exactly known retardation function K(t) = k0 exp(-a t):  B(w) = k0 a / (a^2 + w^2),  A(w) = Ainf - k0 / (a^2 + w^2).
const double kK0 = 200.0, kA = 5.0, kAinf = 8.0;
double Bex(double w) { return kK0 * kA / (kA * kA + w * w); }
double Aex(double w) { return kAinf - kK0 / (kA * kA + w * w); }

HydroTable lorentzTable(double wmax, double dw) {
    HydroTable t;
    for (double w = 0.001; w <= wmax; w += (w < 0.5 ? 0.05 : dw)) {
        t.omega.push_back(w);
        Mat6 A, B; A(2, 2) = Aex(w); B(2, 2) = Bex(w);
        t.A.push_back(A); t.B.push_back(B);
    }
    return t;
}
}  // namespace

TEST_CASE("Retardation function from a tabulated B(omega): K(t) = k0 exp(-a t); Ogilvie relations; A_inf from the table") {
    const HydroTable t = lorentzTable(300.0, 0.1);
    const RadiationModel rm(t, 4.0, 0.01, true);
    for (double tt : {0.0, 0.1, 0.5, 1.0, 2.0}) {
        const double ex = kK0 * std::exp(-kA * tt);
        std::printf("\n  K(%.1f) = %.5f (exact %.5f)", tt, rm.K(2, 2, tt), ex);
        CHECK(std::fabs(rm.K(2, 2, tt) - ex) < 1e-4 * kK0);        // absolute: the table (piecewise-linear B, spacing 0.1, up to 300 rad/s) limits K to ~1e-4 k0
    }
    std::printf("\n  A_inf from the Ogilvie relation at the highest frequency: %.5f (exact %.5f)\n", rm.Ainf()(2, 2), kAinf);
    CHECK(rm.Ainf()(2, 2) == doctest::Approx(kAinf).scale(0).epsilon(2e-3));
    for (double w : {1.0, 3.0, 6.0, 12.0}) {
        CHECK(std::fabs(rm.AFromKernel(2, 2, w) - Aex(w)) < 4e-3 * kAinf);
        CHECK(std::fabs(rm.BFromKernel(2, 2, w) - Bex(w)) < 4e-3 * Bex(0.0));
    }
    CHECK(rm.pairs().size() == 1);
    CHECK(rm.K(0, 0, 0.3) == 0.0);                     // channels without radiation damping have no kernel
    CHECK(rm.K(2, 2, 4.5) == 0.0);                     // beyond the kernel length
    CHECK_THROWS(RadiationModel(HydroTable(), 4.0, 0.01));
}

namespace {
// Heave of a body (M = 10, C = 900) with the Lorentzian radiation data forced by F sin(w t): maximum difference from the frequency-domain steady state.
double forcedHeaveError(double dt, double hK, double w, double F, double* ampOut) {
    HydroTable t = lorentzTable(300.0, 0.1);
    auto rm = std::make_shared<RadiationModel>(t, 3.0, hK, true);
    BodyParams bp; bp.mass = 10.0; bp.C(2, 2) = 900.0; bp.A(2, 2) = rm->Ainf()(2, 2); bp.cgRef = Vec3(0, 0, 0);
    CoupledSystem sys(bp);
    sys.setRadiation(rm);
    sys.setExcitation([&](double tt) { Vec6 f{}; f[2] = F * std::sin(w * tt); return f; });
    sys.begin(0.0);
    const std::complex<double> Z(-w * w * (10.0 + Aex(w)) + 900.0, w * Bex(w));
    double worst = 0, amp = 0;
    while (sys.time() < 40.0 - 1e-9) {
        sys.step(dt);
        if (sys.time() > 30.0) {
            const double ex = std::imag(F * std::exp(std::complex<double>(0, w * sys.time())) / Z);
            worst = std::max(worst, std::fabs(sys.body().xi[2] - ex)); amp = std::max(amp, std::fabs(ex));
        }
    }
    if (ampOut) *ampOut = amp;
    return worst;
}
}  // namespace

TEST_CASE("Cummins time-domain heave response equals the frequency-domain steady state with A(omega), B(omega)") {
    double amp;
    double prev = 1e300;
    for (double dt : {0.04, 0.02, 0.01}) {
        const double e = forcedHeaveError(dt, 0.0025, 3.0, 50.0, &amp);
        std::printf("\n  dt %.3f: max |x - x_fd| = %.3e m for an amplitude of %.3e m (%.1e relative)", dt, e, amp, e / amp);
        CHECK(e < 5e-4 * amp);
        CHECK(e <= prev * 1.01);               // the error does not grow when the body step is reduced (it ends at the floor set by the tabulated kernel, ~4e-5)
        prev = e;
    }
    std::printf("\n");
    // near resonance (the added mass drops: omega_0 ~ 7.6 rad/s instead of sqrt(900/(10 + A_inf)) = 7.07) the frequency dependence matters: compare
    const double w = 7.0;
    const double xCons = 50.0 / std::fabs(-w * w * (10.0 + kAinf) + 900.0 + 1e-12);                  // constant A_inf, no damping: far off
    const double xFd = 50.0 / std::abs(std::complex<double>(-w * w * (10.0 + Aex(w)) + 900.0, w * Bex(w)));
    std::printf("  at omega = 7: amplitude with A(omega), B(omega) %.4e m; with constant A_inf and no radiation damping %.4e m\n", xFd, xCons);
    CHECK(xCons > 2.0 * xFd);
    double amp7;
    const double e7 = forcedHeaveError(0.01, 0.0025, w, 50.0, &amp7);
    CHECK(e7 < 5e-3 * amp7);
}

TEST_CASE("Frequency-dependent hydrodynamics through the JSON runner in an irregular sea, and natural frequency with A(omega)") {
    json rows = json::array();
    const HydroTable t = lorentzTable(300.0, 0.5);
    for (size_t i = 0; i < t.omega.size(); ++i) rows.push_back({{"omega_rad_s", t.omega[i]}, {"A_diag", {0, 0, t.A[i](2, 2), 0, 0, 0}}, {"B_diag", {0, 0, t.B[i](2, 2), 0, 0, 0}}});
    json cfg = json::parse(R"({
      "body": {"mass_kg": 10.0, "inertia_diag_kg_m2": [1,1,1], "cg_ref_m": [0,0,0], "C_diag": [0,0,900.0,0,0,0], "Dq": [0,0,0,0,0,0], "F0": [0,0,0,0,0,0],
               "wave_force_table": [ {"period_s": 0.4, "w": [0,0,60,0,0,0], "delta": [0,0,0.1,0,0,0]}, {"period_s": 3.6, "w": [0,0,300,0,0,0], "delta": [0,0,-0.6,0,0,0]} ]},
      "lines": [],
      "waves": {"depth_m": 30.0, "surface_z_m": 0.0, "spectrum": {"type": "jonswap", "hs_m": 0.1, "tp_s": 1.7, "gamma": 3.3, "n_components": 30, "seed": 5}},
      "numerics": {"dt_body_s": 0.005, "t_end_s": 60.0},
      "output": {"write": false, "dt_out_s": 0.05}
    })");
    cfg["body"]["hydro_table"] = {{"rows", rows}, {"A_inf_diag", {0, 0, kAinf, 0, 0, 0}}, {"retardation_length_s", 3.0}, {"retardation_step_s", 0.0025}};
    const CoupledResult res = runCoupledCase(cfg);
    REQUIRE(res.finite);
    const RadiationModel rm(t, 3.0, 0.0025, true);
    WaveField wf(30.0, 0.0);
    SpectrumSpec sp; sp.hs = 0.1; sp.tp = 1.7; sp.gamma = 3.3; sp.n = 30; sp.seed = 5;
    wf.addSpectrum(sp);
    auto force = [](double om, double* w, double* d) {          // linear interpolation between the two table rows, as in the runner
        const double o1 = 2 * kPi / 3.6, o2 = 2 * kPi / 0.4, f = (om - o1) / (o2 - o1);
        *w = 300 * (1 - f) + 60 * f; *d = -0.6 * (1 - f) + 0.1 * f;
    };
    double worst = 0, amp = 0;
    for (const CoupledSample& s : res.samples) {
        if (s.t < 45.0) continue;
        double z = 0;
        for (const WaveComponent& c : wf.components()) {
            double w, d; force(c.omega, &w, &d);
            const std::complex<double> den(-c.omega * c.omega * (10.0 + rm.AFromKernel(2, 2, c.omega)) + 900.0, c.omega * rm.BFromKernel(2, 2, c.omega));
            z += std::imag(w * c.amplitude * std::exp(std::complex<double>(0, c.omega * s.t + c.phase + d)) / den);
        }
        worst = std::max(worst, std::fabs(s.xi[2] - z)); amp = std::max(amp, std::fabs(z));
    }
    std::printf("\n  irregular sea, heave with A(omega), B(omega): max difference %.3e m for an amplitude of %.3e m\n", worst, amp);
    CHECK(amp > 1e-4);
    CHECK(worst < 5e-3 * amp);
    // natural frequency with A(omega): omega^2 (M + A(omega)) = C solved by the fixed-point iteration, against a direct root-find on the same table
    Mat6 K; K(2, 2) = 900.0;
    Mat6 Mrb; for (int i = 0; i < 6; ++i) Mrb(i, i) = 10.0;
    const BodyModes bm = bodyModesRadiation(Mrb, rm, K);
    double lo = 3.0, hi = 20.0;
    for (int i = 0; i < 100; ++i) { const double m = 0.5 * (lo + hi); (m * m * (10.0 + rm.A(m)(2, 2)) > 900.0 ? hi : lo) = m; }
    double found = -1; for (int k = 0; k < 6; ++k) if (std::fabs(bm.omega[k] - 0.5 * (lo + hi)) < 1e-6) found = bm.omega[k];
    std::printf("  heave natural frequency: constant A_inf %.5f rad/s, A(omega) %.5f rad/s (root %.5f)\n", std::sqrt(900.0 / (10.0 + kAinf)), found, 0.5 * (lo + hi));
    CHECK(found > 0);
    CHECK(found > std::sqrt(900.0 / (10.0 + kAinf)));            // A(omega) < A_inf: higher frequency
    // the implicit coupling cannot be combined with the memory force
    json bad = cfg; bad["numerics"]["coupling"] = "implicit";
    CHECK_THROWS(runCoupledCase(bad));
}
