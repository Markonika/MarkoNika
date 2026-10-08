// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#include "mooring/catenary.hpp"
#include <algorithm>
#include <cmath>

namespace mooring {

void ElasticCatenary::position(double s, double& x, double& z) const {
    const double Vs = V0 + w * s;
    x = H * s / EA + (H / w) * (std::asinh(Vs / H) - std::asinh(V0 / H));
    z = (V0 * s + 0.5 * w * s * s) / EA + (std::sqrt(H * H + Vs * Vs) - std::sqrt(H * H + V0 * V0)) / w;
}

double ElasticCatenary::tension(double s) const {
    const double Vs = V0 + w * s;
    return std::sqrt(H * H + Vs * Vs);
}

ElasticCatenary ElasticCatenary::solve(double L, double EA, double w, double dx, double dz) {
    ElasticCatenary c;
    c.L = L; c.EA = EA; c.w = w;
    // Initial guess: parabolic sag estimate for the chord, then correct the vertical split.
    const double chord = std::sqrt(dx * dx + dz * dz);
    const double sag = std::sqrt(std::max(3.0 * chord * std::max(L - chord, 1e-3 * L) / 8.0, 1e-6));
    double H = std::max(w * dx * dx / (8.0 * sag), 1e-3 * w * L);
    double V0 = -0.5 * w * L + H * dz / dx;

    auto resid = [&](double h, double v, double& rx, double& rz) {
        ElasticCatenary t = c; t.H = h; t.V0 = v;
        double x, z; t.position(L, x, z);
        rx = x - dx; rz = z - dz;
    };
    double rx, rz; resid(H, V0, rx, rz);
    for (int it = 0; it < 200; ++it) {
        const double r0 = std::hypot(rx, rz);
        c.residual = r0;
        if (r0 < 1e-12 * std::max(1.0, chord)) break;
        const double eh = 1e-7 * H, ev = 1e-7 * std::max(std::fabs(V0), w * L);
        double ax, az, bx, bz;
        resid(H + eh, V0, ax, az); resid(H, V0 + ev, bx, bz);
        const double a11 = (ax - rx) / eh, a21 = (az - rz) / eh;
        const double a12 = (bx - rx) / ev, a22 = (bz - rz) / ev;
        const double det = a11 * a22 - a12 * a21;
        if (std::fabs(det) < 1e-300) break;
        double dH = -( a22 * rx - a12 * rz) / det;
        double dV = -(-a21 * rx + a11 * rz) / det;
        double lam = 1.0;                       // backtracking: keep H > 0 and reduce residual
        for (int k = 0; k < 40; ++k, lam *= 0.5) {
            const double h = H + lam * dH, v = V0 + lam * dV;
            if (h <= 0) continue;
            double tx, tz; resid(h, v, tx, tz);
            if (std::hypot(tx, tz) < r0 || k == 39) { H = h; V0 = v; rx = tx; rz = tz; break; }
        }
    }
    c.H = H; c.V0 = V0;
    c.residual = std::hypot(rx, rz);
    c.converged = c.residual < 1e-9 * std::max(1.0, chord);
    return c;
}

}  // namespace mooring
