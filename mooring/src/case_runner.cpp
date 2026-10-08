// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#include "mooring/case_runner.hpp"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include "mooring/lumped_mass_cable.hpp"

namespace mooring {
using json = nlohmann::json;
namespace {
const double kPi = 3.14159265358979323846;

Vec3 vec(const json& j) { return Vec3(j.at(0).get<double>(), j.at(1).get<double>(), j.at(2).get<double>()); }

template <class T> T get(const json& j, const char* key, T def) { return j.contains(key) ? j[key].get<T>() : def; }
}  // namespace

void applyOverride(json& cfg, const std::string& a) {
    const auto eq = a.find('=');
    if (eq == std::string::npos) throw std::invalid_argument("override needs key=value: " + a);
    std::string path = a.substr(0, eq), val = a.substr(eq + 1);
    json v = json::parse(val, nullptr, false);
    if (v.is_discarded()) v = val;
    json* node = &cfg;
    size_t pos = 0;
    while (true) {
        const size_t dot = path.find('.', pos);
        const std::string key = path.substr(pos, dot == std::string::npos ? std::string::npos : dot - pos);
        if (dot == std::string::npos) { (*node)[key] = v; break; }
        node = &(*node)[key];
        pos = dot + 1;
    }
}

CaseResult runCase(const json& cfg) {
    const auto wall0 = std::chrono::steady_clock::now();
    const json& jl = cfg.at("line");
    const json jenv = cfg.value("environment", json::object());
    const json jsoil = cfg.value("soil", json::object());
    const json jm = cfg.value("motion", json::object());
    const json jn = cfg.value("numerics", json::object());
    const json jo = cfg.value("output", json::object());

    CableParams p;
    p.L = jl.at("length_m"); p.N = jl.at("segments"); p.EA = jl.at("EA_N");
    p.m_l = jl.at("mass_per_length_kg_m"); p.g = get(jl, "gravity_m_s2", 9.81);
    p.c_int = get(jl, "internal_damping_Ns", 0.0);
    p.D0 = get(jl, "hydro_diameter_m", 0.0); p.D1 = get(jl, "soil_diameter_m", p.D0);
    p.A1 = get(jl, "nominal_area_m2", 0.0);
    p.planar = get(jn, "planar", false);
    p.Cm = get(jl, "Cm", 0.0); p.Cdt = get(jl, "Cdt", 0.0); p.Cdn = get(jl, "Cdn", 0.0);
    Environment env;
    env.hydro = get(jenv, "hydro", false); env.seabed = get(jenv, "seabed", false);
    env.rho_w = get(jenv, "rho_w_kg_m3", 1000.0);
    env.surfaceZ = get(jenv, "water_surface_z_m", 0.0); env.seabedZ = get(jenv, "seabed_z_m", 0.0);
    const double rho_c = get(jl, "density_kg_m3", 0.0);
    p.w = rho_c > 0.0 ? CableParams::submergedWeight(p.m_l, rho_c, env.rho_w, p.g) : get(jl, "weight_per_length_N_m", p.m_l * p.g);
    p.soil.Ks = get(jsoil, "stiffness_Pa_per_m", 0.0); p.soil.zeta = get(jsoil, "damping_factor", 1.0);
    p.soil.mu = get(jsoil, "friction", 0.0); p.soil.vlim = get(jsoil, "v_lim_m_s", 0.01);

    const Vec3 anchor = vec(cfg.at("anchor_m")), rest = vec(cfg.at("fairlead_rest_m"));
    if (jenv.contains("current_m_s")) {              // uniform steady current (no water acceleration)
        const Vec3 U = vec(jenv["current_m_s"]);
        env.water = [U](const Vec3&, double, Vec3& vw, Vec3& aw) { vw = U; aw = Vec3(); };
    }
    LumpedMassCable cable(p, anchor, rest);
    cable.setEnvironment(env);

    CaseResult res;
    res.l0 = p.l0(); res.cWave = p.waveSpeed();
    const bool touchdown = get<std::string>(cfg.value("initial", json::object()), "shape", "touchdown") == "touchdown";
    if (touchdown && env.seabed) cable.initTouchdownCatenary();
    RelaxOptions ro;
    ro.forceTol = get(jn, "relax_force_tol", 1e-6); ro.maxSteps = get<long>(jn, "relax_max_steps", 20000000L);
    const RelaxResult rr = cable.relaxStatic(ro);
    res.staticConverged = rr.converged; res.staticSteps = rr.steps;
    res.staticTopTension = norm(cable.endTension(true));
    res.staticAnchorTension = norm(cable.endTension(false));

    const double pert = get(cfg.value("initial", json::object()), "perturbation_y_m", 0.0);
    std::vector<Vec3> rp = cable.nodes();
    if (pert != 0.0 && !p.planar) {                    // small out-of-plane half-sine, released from rest
        for (int i = 1; i < p.N; ++i) rp[i].y += pert * std::sin(kPi * i / p.N);
        cable.setInitialState(rp, std::vector<Vec3>(rp.size()), 0.0);
    }
    double maxY = 0.0;
    DynOptions dyn;
    dyn.scheme = get<std::string>(jn, "scheme", "rk4") == "verlet" ? Scheme::Verlet : Scheme::RK4;
    dyn.cfl = get(jn, "cfl", 0.5); dyn.dt = get(jn, "dt_s", 0.0);
    cable.setDynOptions(dyn);
    res.dt = cable.stableDt();

    // ---- prescribed top-end motion: circle in the x-z plane about the rest position, cosine ramp ----
    const std::string mtype = get<std::string>(jm, "type", "none");
    const double planeAng = get(jm, "plane_angle_deg", 0.0) * kPi / 180.0;   // rotation of the motion plane about z
    const Vec3 ex(std::cos(planeAng), std::sin(planeAng), 0.0);
    const double radius = get(jm, "radius_m", 0.0), period = get(jm, "period_s", 1.0);
    const double dir = get(jm, "direction", 1.0), phase = get(jm, "phase_deg", 0.0) * kPi / 180.0;
    const double rampCycles = get(jm, "ramp_cycles", 2.0), nCycles = get(jm, "cycles", 10.0);
    const double om = 2.0 * kPi / period, Tramp = rampCycles * period;
    const Vec3 centre = jm.contains("centre_m") ? vec(jm["centre_m"]) : rest;
    if (mtype == "circle_xz") {
        cable.setTopMotion([=](double t, Vec3& pos, Vec3& vel) {
            double rho = 1.0, drho = 0.0;
            if (t < Tramp) { rho = 0.5 * (1.0 - std::cos(kPi * t / Tramp)); drho = 0.5 * kPi / Tramp * std::sin(kPi * t / Tramp); }
            const double th = phase + dir * om * t;
            const Vec3 e = ex * std::cos(th) + Vec3(0, 0, std::sin(th)), de = ex * (-std::sin(th) * dir * om) + Vec3(0, 0, std::cos(th) * dir * om);
            pos = centre + e * (radius * rho);
            vel = e * (radius * drho) + de * (radius * rho);
        });
    } else if (mtype != "none") {
        throw std::invalid_argument("unknown motion type: " + mtype);
    }

    // ---- output ----
    const bool write = get(jo, "write", true);
    const std::string dir_out = get<std::string>(jo, "directory", "out"), tag = get<std::string>(jo, "tag", "case");
    const double dtOut = get(jo, "dt_out_s", 1e-2);
    std::ofstream ts;
    if (write) {
        std::filesystem::create_directories(dir_out);
        ts.open(dir_out + "/" + tag + "_timeseries.csv");
        ts << "time_s,top_tension_N,anchor_tension_N,top_x_m,top_y_m,top_z_m\n";
    }
    std::map<int, CycleMax> cyc;
    bool finite = true;
    const Vec3 nrm(-std::sin(planeAng), std::cos(planeAng), 0.0);
    cable.setStepObserver([&](const LumpedMassCable& c) {
        if ((c.stats().steps & 63) == 0)
            for (const Vec3& q : c.nodes()) maxY = std::max(maxY, std::fabs(dot(q - anchor, nrm)));
        const double T = c.endSegmentTension(true), A = c.endSegmentTension(false);
        if (!std::isfinite(T) || !std::isfinite(A)) { finite = false; return; }
        const int k = static_cast<int>(c.time() / period);
        auto it = cyc.find(k);
        if (it == cyc.end()) it = cyc.emplace(k, CycleMax{k, 0.0, 0.0}).first;
        it->second.topMax = std::max(it->second.topMax, T);
        it->second.anchorMax = std::max(it->second.anchorMax, A);
    });

    if (mtype != "none") {
        const double tEnd = nCycles * period;
        double nextOut = 0.0;
        auto emit = [&](double t) {
            if (!write) return;
            const Vec3& tp = cable.nodes().back();
            char buf[256];
            std::snprintf(buf, sizeof buf, "%.6f,%.9g,%.9g,%.9g,%.9g,%.9g\n", t, cable.endSegmentTension(true),
                          cable.endSegmentTension(false), tp.x, tp.y, tp.z);
            ts << buf;
        };
        emit(0.0);
        while (cable.time() < tEnd - 1e-12 && finite) {
            nextOut = std::min(nextOut + dtOut, tEnd);
            cable.advanceTo(nextOut);
            emit(cable.time());
        }
    }
    res.maxOutOfPlane = maxY;
    for (const Vec3& q : cable.nodes()) res.finalOutOfPlane = std::max(res.finalOutOfPlane, std::fabs(dot(q - anchor, nrm)));
    res.finite = finite;
    res.steps = cable.stats().steps;
    res.slackSegmentEvals = cable.stats().slackSegmentEvals;
    res.clippedTensionEvals = cable.stats().clippedTensionEvals;
    res.soilContactEvals = cable.stats().soilContactEvals;

    const int first = get<int>(cfg.value("statistics", json::object()), "first_cycle", 5);
    double sum = 0.0; int n = 0;
    const int lastComplete = static_cast<int>(nCycles) - 1;
    for (auto& kv : cyc) {
        res.cycles.push_back(kv.second);
        res.maxTopOverall = std::max(res.maxTopOverall, kv.second.topMax);
        if (kv.first >= first && kv.first <= lastComplete) { sum += kv.second.topMax; ++n; }
    }
    res.cyclesAveraged = n;
    res.meanMax = n ? sum / n : 0.0;
    res.wallSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();

    if (write) {
        std::ofstream cy(dir_out + "/" + tag + "_cycles.csv");
        cy << "cycle,top_max_N,anchor_max_N\n";
        for (const auto& c : res.cycles) cy << c.cycle << "," << c.topMax << "," << c.anchorMax << "\n";
        json log = cfg;
        log["derived"] = {{"submerged_weight_N_m", p.w}, {"dry_weight_N_m", p.dryWeight()}, {"l0_m", res.l0},
                          {"axial_wave_speed_m_s", res.cWave}, {"dt_used_s", res.dt}, {"steps", res.steps},
                          {"static_converged", res.staticConverged}, {"static_relax_steps", res.staticSteps},
                          {"static_top_tension_N", res.staticTopTension}, {"static_anchor_tension_N", res.staticAnchorTension},
                          {"slack_segment_evals", res.slackSegmentEvals}, {"clipped_tension_evals", res.clippedTensionEvals},
                          {"soil_contact_evals", res.soilContactEvals}, {"mean_cycle_max_N", res.meanMax},
                          {"cycles_averaged", res.cyclesAveraged}, {"finite", res.finite}, {"wall_s", res.wallSeconds}};
        std::ofstream(dir_out + "/" + tag + "_params.json") << log.dump(2) << "\n";
    }
    return res;
}

}  // namespace mooring
