// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#include "mooring/coupled_runner.hpp"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include "mooring/case_runner.hpp"
#include "mooring/radiation.hpp"
#include "mooring/waves.hpp"
#include "mooring/lumped_mass_cable.hpp"

namespace mooring {
using json = nlohmann::json;
namespace {

Vec3 vec(const json& j) { return Vec3(j.at(0).get<double>(), j.at(1).get<double>(), j.at(2).get<double>()); }

// A 6x6 matrix from "<key>" (36 numbers or 6 rows) or "<key>_diag" (6 numbers); zero if absent.
Mat6 matrix6(const json& j, const std::string& key) {
    Mat6 m;
    if (j.contains(key + "_diag")) { for (int i = 0; i < 6; ++i) m(i, i) = j[key + "_diag"].at(i).get<double>(); }
    else if (j.contains(key)) {
        const json& a = j[key];
        if (a.size() == 6 && a[0].is_array()) { for (int i = 0; i < 6; ++i) for (int k = 0; k < 6; ++k) m(i, k) = a[i][k].get<double>(); }
        else if (a.size() == 36) { for (int i = 0; i < 36; ++i) m.a[i] = a[i].get<double>(); }
        else throw std::invalid_argument("matrix '" + key + "' must be 6x6 or 36 numbers");
    }
    return m;
}
Vec6 vector6(const json& j, const std::string& key) {
    Vec6 v{};
    if (j.contains(key)) { if (j[key].size() != 6) throw std::invalid_argument(key + " needs 6 numbers"); for (int i = 0; i < 6; ++i) v[i] = j[key][i].get<double>(); }
    return v;
}

}  // namespace

std::unique_ptr<CoupledSystem> buildCoupledSystem(const json& cfg, std::vector<std::string>* lineNames) {
    const json& jb = cfg.at("body");
    BodyParams bp;
    bp.mass = jb.at("mass_kg");
    if (jb.contains("inertia_diag_kg_m2")) { bp.inertia = Mat3{}; for (int i = 0; i < 3; ++i) bp.inertia[4 * i] = jb["inertia_diag_kg_m2"].at(i).get<double>(); }
    else if (jb.contains("inertia_kg_m2")) { for (int i = 0; i < 9; ++i) bp.inertia[i] = jb["inertia_kg_m2"].at(i).get<double>(); }
    bp.A = matrix6(jb, "A"); bp.B = matrix6(jb, "B"); bp.C = matrix6(jb, "C");
    bp.Dq = vector6(jb, "Dq"); bp.F0 = vector6(jb, "F0");
    bp.cgRef = vec(jb.at("cg_ref_m"));
    // Frequency-dependent added mass / radiation damping from a table (BEM output): the body mass matrix gets A_inf (plus the constant body.A as an
    // extra), the memory force acts through the retardation function K(t) built from B(omega) (docs/radiation.md).
    std::shared_ptr<RadiationModel> radiation;
    if (jb.contains("hydro_table")) {
        const json& jh = jb["hydro_table"];
        HydroTable ht;
        struct Row { double om; Mat6 A, B; };
        std::vector<Row> rows;
        for (const json& r : jh.at("rows")) {
            const double om = r.contains("omega_rad_s") ? r["omega_rad_s"].get<double>() : 2.0 * 3.14159265358979323846 / r.at("period_s").get<double>();
            rows.push_back({om, matrix6(r, "A"), matrix6(r, "B")});
        }
        std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return a.om < b.om; });
        for (const Row& r : rows) { ht.omega.push_back(r.om); ht.A.push_back(r.A); ht.B.push_back(r.B); }
        if (jh.contains("A_inf") || jh.contains("A_inf_diag")) { ht.hasAinf = true; ht.Ainf = matrix6(jh, "A_inf"); }
        radiation = std::make_shared<RadiationModel>(ht, jh.value("retardation_length_s", 30.0), jh.value("retardation_step_s", 0.01), jh.value("high_frequency_tail", true));
        bp.A = bp.A + radiation->Ainf();
    }
    auto sysPtr = std::make_unique<CoupledSystem>(bp);
    CoupledSystem& sys = *sysPtr;
    if (radiation) sys.setRadiation(radiation);

    // Waves: Airy components (+ Wheeler stretching), used for the line/point-element kinematics and, with
    // body.wave_force {w [N/m or N per m of amplitude], delta [rad]}, for the body excitation f_i = w_i A sin(omega t + delta_i) (Eq. 3.60).
    std::shared_ptr<WaveField> waves;
    if (cfg.contains("waves")) {
        const json& jw = cfg["waves"];
        const double depth = jw.at("depth_m");
        const double surfaceZ = jw.value("surface_z_m", cfg.value("environment", json::object()).value("water_surface_z_m", 0.0));
        waves = std::make_shared<WaveField>(depth, surfaceZ, jw.value("stretching", std::string("wheeler")) == "none" ? Stretching::None : Stretching::Wheeler);
        if (jw.contains("components"))
            for (const json& c : jw.at("components"))
                waves->addComponent(c.contains("amplitude_m") ? c["amplitude_m"].get<double>() : 0.5 * c.at("height_m").get<double>(),
                                    c.at("period_s"), c.value("phase_rad", 0.0), c.value("direction_deg", 0.0));
        if (jw.contains("spectrum")) {                                  // JONSWAP / Pierson-Moskowitz, random-frequency discretisation
            const json& js = jw["spectrum"];
            SpectrumSpec sp;
            sp.hs = js.at("hs_m"); sp.tp = js.at("tp_s");
            const std::string st = js.value("type", std::string("jonswap"));
            sp.gamma = st == "pm" ? 1.0 : js.value("gamma", 3.3);
            sp.n = js.value("n_components", 50); sp.seed = js.value("seed", 1u);
            sp.wMinFactor = js.value("omega_min_factor", 0.5); sp.wMaxFactor = js.value("omega_max_factor", 4.0);
            sp.directionDeg = js.value("direction_deg", 0.0);
            waves->addSpectrum(sp);
        }
        if (waves->components().empty()) throw std::invalid_argument("waves: give components and/or a spectrum");
        waves->setRampTime(jw.value("ramp_time_s", 0.0));
        if (jb.contains("wave_force_table")) {
            // Frequency-dependent first-order body load (e.g. from a BEM solution): per component f_i = w_i(omega) a sin(omega t + phi + delta_i(omega)),
            // w and delta interpolated linearly in omega between the table rows; a component outside the table range is an error.
            struct Row { double om; Vec6 w, dl; };
            std::vector<Row> rows;
            for (const json& r : jb["wave_force_table"]) rows.push_back({2.0 * 3.14159265358979323846 / r.at("period_s").get<double>(), vector6(r, "w"), vector6(r, "delta")});
            std::sort(rows.begin(), rows.end(), [](const Row& x, const Row& y) { return x.om < y.om; });
            if (rows.empty()) throw std::invalid_argument("body.wave_force_table is empty");
            struct Comp { double om, a, ph; Vec6 w, dl; };
            std::vector<Comp> cs;
            for (const WaveComponent& c : waves->components()) {
                if (std::fabs(c.direction) > 1e-12) throw std::invalid_argument("body.wave_force_table: the load coefficients are for waves along +x (direction_deg = 0)");
                if (c.omega < rows.front().om * (1.0 - 1e-9) || c.omega > rows.back().om * (1.0 + 1e-9))
                    throw std::invalid_argument("a wave component lies outside the period range of body.wave_force_table");
                Comp k; k.om = c.omega; k.a = c.amplitude; k.ph = c.phase;
                size_t j = 0;
                while (j + 2 < rows.size() && c.omega > rows[j + 1].om) ++j;
                if (rows.size() == 1) { k.w = rows[0].w; k.dl = rows[0].dl; }
                else {
                    const double f = std::min(1.0, std::max(0.0, (c.omega - rows[j].om) / (rows[j + 1].om - rows[j].om)));
                    for (int i = 0; i < 6; ++i) { k.w[i] = rows[j].w[i] * (1 - f) + rows[j + 1].w[i] * f; k.dl[i] = rows[j].dl[i] * (1 - f) + rows[j + 1].dl[i] * f; }
                }
                cs.push_back(k);
            }
            auto wp = waves;
            sys.setExcitation([cs, wp](double t) {
                Vec6 F{};
                const double r = wp->rampFactor(t);
                for (const Comp& k : cs) for (int i = 0; i < 6; ++i) F[i] += k.w[i] * k.a * std::sin(k.om * t + k.ph + k.dl[i]);
                for (int i = 0; i < 6; ++i) F[i] *= r;
                return F;
            });
        } else if (jb.contains("wave_force")) {
            if (waves->components().size() != 1) throw std::invalid_argument("body.wave_force is for one regular wave; use body.wave_force_table for several components or a spectrum");
            const Vec6 w = vector6(jb["wave_force"], "w"), dl = vector6(jb["wave_force"], "delta");
            const WaveComponent c0 = waves->components().at(0);          // the force coefficients belong to the first (regular) component
            auto wp = waves;
            sys.setExcitation([w, dl, c0, wp](double t) {
                Vec6 F{};
                const double r = wp->rampFactor(t);
                for (int k = 0; k < 6; ++k) F[k] = w[k] * c0.amplitude * r * std::sin(c0.omega * t + c0.phase + dl[k]);
                return F;
            });
        }
    }

    const json jn = cfg.value("numerics", json::object());
    const json jinit = cfg.value("initial", json::object());
    sys.body().xi = vector6(jinit, "xi0");
    sys.body().xiDot = vector6(jinit, "xi_dot0");

    for (const json& jl : cfg.at("lines")) {
        json lc = jl;                                              // shared defaults, overridden by the line's own entries
        for (const char* k : {"environment", "soil", "numerics"})
            if (cfg.contains(k) && !lc.contains(k)) lc[k] = cfg[k];
        if (cfg.contains("line_defaults")) {
            json merged = cfg["line_defaults"]; merged.update(lc.value("line", json::object())); lc["line"] = merged;
        }
        const Vec3 a = vec(jl.at("fairlead_body_m")), anchor = vec(jl.at("anchor_m"));
        std::unique_ptr<LumpedMassCable> cab = buildLine(lc, anchor, sys.body().pointPosition(a));
        if (waves) { Environment e = cab->environment(); e.water = withCurrent(waves->asWaterField(), makeCurrentProfile(lc.value("environment", json::object()), e.surfaceZ));
            if (cfg.at("waves").value("instantaneous_surface", false)) e.elevation = waves->asElevation();
            cab->setEnvironment(e); }
        if (jn.contains("cfl") || jn.contains("line_dt_s") || jn.contains("scheme")) {
            DynOptions o; parseDynOptions(jn, o); o.dt = jn.value("line_dt_s", 0.0);
            cab->setDynOptions(o);
        }
        const std::string nm = jl.value("name", "line" + std::to_string(sys.numLines() + 1));
        if (lineNames) lineNames->push_back(nm);
        sys.addLine(std::move(cab), a, nm);
    }
    sys.setThreads(jn.value("threads", 1));
    if (jn.value("coupling", std::string("explicit")) == "implicit") {
        ImplicitCoupling ic; ic.enabled = true;
        ic.newmarkGamma = jn.value("newmark_gamma", 0.5); ic.tol = jn.value("coupling_tol", 1e-9); ic.maxIter = jn.value("coupling_max_iter", 20);
        sys.setImplicitCoupling(ic);
    }
    return sysPtr;
}

CoupledResult runCoupledCase(const json& cfg) {
    const auto wall0 = std::chrono::steady_clock::now();
    CoupledResult res;
    std::unique_ptr<CoupledSystem> sysPtr = buildCoupledSystem(cfg, &res.lineNames);
    CoupledSystem& sys = *sysPtr;
    const json jn = cfg.value("numerics", json::object());
    const json jinit = cfg.value("initial", json::object());

    res.equilibriumRequested = jinit.value("equilibrium", false);
    if (res.equilibriumRequested) {
        res.equilibrium = sys.solveEquilibrium(jinit.value("equilibrium_tol", 1e-6));
        if (jinit.contains("xi_offset")) {                         // displace from equilibrium, e.g. for a decay test
            const Vec6 off = vector6(jinit, "xi_offset");
            for (int k = 0; k < 6; ++k) sys.body().xi[k] += off[k];
            sys.staticResidual();
        }
        if (jinit.contains("xi_dot0")) sys.body().xiDot = vector6(jinit, "xi_dot0");     // the equilibrium solve leaves the body at rest
    } else {
        sys.staticResidual();                                      // lines at their static shape for the given pose
    }
    res.xiStart = sys.body().xi;
    sys.begin(0.0);

    const double dt = jn.at("dt_body_s"), tEnd = jn.at("t_end_s");
    const json jo = cfg.value("output", json::object());
    const double dtOut = jo.value("dt_out_s", 0.01);
    const bool write = jo.value("write", false);
    const std::string dir = jo.value("directory", "out"), tag = jo.value("tag", "coupled");
    std::ofstream ts;
    if (write) {
        std::filesystem::create_directories(dir);
        ts.open(dir + "/" + tag + "_timeseries.csv");
        ts << "time_s,surge,sway,heave,roll,pitch,yaw";
        for (const auto& n : res.lineNames) ts << "," << n << "_force_N," << n << "_top_tension_N";
        ts << "\n";
    }
    auto sample = [&](double t) {
        CoupledSample s; s.t = t; s.xi = sys.body().xi;
        for (int k = 0; k < sys.numLines(); ++k) {
            s.lineForce.push_back(norm(sys.lineForce(k)));
            const auto* lm = dynamic_cast<const LumpedMassCable*>(&sys.line(k));
            s.topTension.push_back(lm ? lm->endSegmentTension(true) : 0.0);
        }
        if (write) {
            char buf[128]; std::snprintf(buf, sizeof buf, "%.6f,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g", t, s.xi[0], s.xi[1], s.xi[2], s.xi[3], s.xi[4], s.xi[5]);
            ts << buf;
            for (size_t k = 0; k < s.lineForce.size(); ++k) { std::snprintf(buf, sizeof buf, ",%.9g,%.9g", s.lineForce[k], s.topTension[k]); ts << buf; }
            ts << "\n";
        }
        for (double v : s.xi) if (!std::isfinite(v)) res.finite = false;
        res.samples.push_back(std::move(s));
    };
    sample(0.0);
    double nextOut = dtOut;
    while (sys.time() < tEnd - 1e-9 * dt && res.finite) {
        sys.step(std::min(dt, tEnd - sys.time()));
        if (sys.time() >= nextOut - 1e-9 * dt) { sample(sys.time()); nextOut += dtOut; }
    }
    res.report = sys.report();
    res.wallSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
    if (write) {
        json log = cfg;
        log["derived"] = {{"dt_body_s", res.report.dtBody}, {"dt_line_min_s", res.report.dtLineMin}, {"sub_step_ratio", res.report.subStepRatio},
                          {"coupling", "explicit partitioned (conventional serial staggered); body velocity-Verlet, implicit linear damping"},
                          {"body_steps", res.report.bodySteps}, {"equilibrium_converged", res.equilibrium.converged},
                          {"equilibrium_iterations", res.equilibrium.iterations}, {"xi_start", res.xiStart}, {"finite", res.finite}, {"wall_s", res.wallSeconds}};
        std::ofstream(dir + "/" + tag + "_params.json") << log.dump(2) << "\n";
    }
    return res;
}

}  // namespace mooring
