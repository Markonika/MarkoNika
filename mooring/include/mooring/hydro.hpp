// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Per-node force laws: Morison hydrodynamics (Paredes 2016 Eqs. 3.27-3.31), seabed (Eqs. 3.32-3.36).
#pragma once
#include "mooring/vec3.hpp"

namespace mooring {

// Solve M a = f for M = m I + c (I - t t^T), t a unit vector (added mass acts normal to the cable only).
// M^-1 = (1/m) t t^T + 1/(m+c) (I - t t^T): exact, no explicit relative acceleration needed.
inline Vec3 addedMassSolve(double m, double c, const Vec3& t, const Vec3& f) {
    const Vec3 ft = t * dot(f, t);
    return ft / m + (f - ft) / (m + c);
}

// Morison drag force on a node [N] (Eqs. 3.28-3.29 with the (1+eps) factor folded into the
// stretched tributary length 'len'). vrel = water velocity - cable velocity. The tangential
// term is the sign-preserving quadratic 0.5 Cdt rho D (v.t)|v.t| t.
inline Vec3 morisonDrag(const Vec3& vrel, const Vec3& t, double Cdt, double Cdn, double rho,
                        double D, double len) {
    const double vt = dot(vrel, t);
    const Vec3 vn = vrel - t * vt;
    return (t * (0.5 * Cdt * rho * D * vt * std::fabs(vt)) +
            vn * (0.5 * Cdn * rho * D * norm(vn))) * len;
}

// Same with different tangential and normal masses: M^-1 = t t^T / mt + (I - t t^T) / mn.
inline Vec3 addedMassSolve2(double mt, double mn, const Vec3& t, const Vec3& f) {
    const Vec3 ft = t * dot(f, t);
    return ft / mt + (f - ft) / mn;
}

struct SoilParams {
    double Ks{0};      // soil stiffness per unit length and diameter [Pa/m = N/m^3]
    double zeta{1.0};  // damping factor
    double mu{0};      // Coulomb friction coefficient
    double vlim{0.01}; // speed over which the friction coefficient is ramped 0 -> mu [m/s]
};

// Seabed force on a node [N] (flat seabed with unit normal ns = +z). 'pen' is the penetration depth
// (seabedZ - z), 'len' the unstretched tributary length, 'wSub' the submerged weight per unit length.
// Normal: spring Ks D1 pen plus damping 2 zeta sqrt(Ks D1 m_l) only while penetrating (Eq. 3.32).
// Tangential: -wSub mu min(|vst|/vlim, 1) vst/|vst| (Eqs. 3.34-3.36), applied while in contact.
inline Vec3 seabedForce(const SoilParams& s, double D1, double m_l, double wSub, double pen,
                        const Vec3& v, double len) {
    if (pen <= 0.0) return Vec3();
    const double vsn = v.z;                                     // v . ns
    const double fn = s.Ks * D1 * pen - 2.0 * s.zeta * std::sqrt(s.Ks * D1 * m_l) * std::min(0.0, vsn);
    const Vec3 vst(v.x, v.y, 0.0);
    const double sp = norm(vst);
    Vec3 ft;
    if (sp > 0.0) ft = vst * (-wSub * s.mu * std::min(sp / s.vlim, 1.0) / sp);
    return (Vec3(0, 0, fn) + ft) * len;
}

}  // namespace mooring
