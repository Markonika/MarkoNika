// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#include "mooring/paredes.hpp"
#include <cmath>
#include "mooring/lumped_mass_cable.hpp"

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

}  // namespace mooring
