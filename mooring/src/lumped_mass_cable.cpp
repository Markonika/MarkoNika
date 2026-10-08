// Copyright (c) [OWNER NAME]. All rights reserved. Proprietary and confidential.
#include "mooring/lumped_mass_cable.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace mooring {

double CableParams::waveSpeed() const { return std::sqrt(EA / m_l); }

LumpedMassCable::LumpedMassCable(const CableParams& p, const Vec3& anchor, const Vec3& fairlead) : p_(p) {
    if (p.N < 2 || p.L <= 0 || p.EA <= 0 || p.m_l <= 0) throw std::invalid_argument("invalid CableParams");
    r_.resize(p.N + 1);
    // Initial shape: straight chord with a parabolic sag matching the slack length.
    const Vec3 d = fairlead - anchor;
    const double chord = norm(d);
    const double slack = std::max(p.L - chord, 0.0);
    const double sag = std::sqrt(3.0 * chord * slack / 8.0);
    for (int i = 0; i <= p.N; ++i) {
        const double u = double(i) / p.N;
        r_[i] = anchor + d * u + Vec3(0, 0, -4.0 * sag * u * (1.0 - u));
    }
}

double LumpedMassCable::nodeMass(int i) const {
    const double m = p_.m_l * p_.l0();
    return (i == 0 || i == p_.N) ? 0.5 * m : m;
}

double LumpedMassCable::segmentTension(const std::vector<Vec3>& r, int seg) const {
    const double eps = norm(r[seg + 1] - r[seg]) / p_.l0() - 1.0;
    return eps > 0.0 ? p_.EA * eps : 0.0;
}

void LumpedMassCable::computeForces(const std::vector<Vec3>& r, const std::vector<Vec3>& /*v*/,
                                    std::vector<Vec3>& f) const {
    f.assign(r.size(), Vec3());
    for (int i = 0; i < p_.N; ++i) {
        const Vec3 d = r[i + 1] - r[i];
        const double len = norm(d);
        const double eps = len / p_.l0() - 1.0;
        if (eps <= 0.0 || len <= 0.0) continue;
        const Vec3 F = d * (p_.EA * eps / len);
        f[i] += F;
        f[i + 1] -= F;
    }
    for (int i = 0; i <= p_.N; ++i) f[i].z -= p_.w * p_.l0() * ((i == 0 || i == p_.N) ? 0.5 : 1.0);
}

Vec3 LumpedMassCable::endTension(bool top) const {
    const int seg = top ? p_.N - 1 : 0;
    const Vec3 d = top ? r_[p_.N - 1] - r_[p_.N] : r_[1] - r_[0];
    const double len = norm(d);
    return len > 0 ? d * (segmentTension(r_, seg) / len) : Vec3();
}

RelaxResult LumpedMassCable::relaxStatic(const RelaxOptions& opt) {
    // Dynamic relaxation with kinetic damping. Only the final equilibrium is physical, so the
    // inertia is fictitious: m_i = massFactor * dt^2 * EA / l0 keeps the explicit scheme stable for
    // any stiffness (massFactor = 2 is the linear limit and was found unstable; the default 32
    // converges robustly), so the pseudo-time step is independent of the axial wave speed.
    RelaxResult res;
    const int n = p_.N + 1;
    res.dt = 1.0;
    const double mFict = opt.massFactor * res.dt * res.dt * p_.EA / p_.l0();
    std::vector<Vec3> v(n), f(n);
    const double tol = opt.forceTol * std::max(p_.w, 1e-12) * p_.l0();
    double keOld = 0.0;
    for (long k = 0; k < opt.maxSteps; ++k) {
        computeForces(r_, v, f);
        double ke = 0.0, rmax = 0.0;
        for (int i = 1; i < p_.N; ++i) {          // end nodes fixed
            v[i] += f[i] * (res.dt / mFict);
            r_[i] += v[i] * res.dt;
            ke += 0.5 * mFict * dot(v[i], v[i]);
            rmax = std::max(rmax, norm(f[i]));
        }
        res.steps = k + 1;
        res.maxResidual = rmax;
        if (rmax < tol) { res.converged = true; break; }
        if (ke < keOld) {                          // kinetic-energy peak passed: restart from rest
            std::fill(v.begin(), v.end(), Vec3());
            ke = 0.0;
            ++res.kineticResets;
        }
        keOld = ke;
    }
    return res;
}

}  // namespace mooring
