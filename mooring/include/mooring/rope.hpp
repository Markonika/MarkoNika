// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace mooring {

// One Maxwell branch of the generalised Maxwell (Zener) rope law: a spring K [N] in series with a dashpot, parallel to the static curve.
// Branch force K (eps - alpha), internal strain alpha obeys d alpha/dt = (eps - alpha)/tau  (tau = c/K, relaxation time [s]).
struct MaxwellBranch { double K{0}; double tau{1}; };

// Phenomenological axial law of a synthetic rope, per segment of strain eps = l/l0 - 1 (tension only, eps <= 0 is slack):
//   T = T_static(eps) + sum_k K_k (eps - alpha_k)   [+ c_int d eps/dt, applied by the cable]
// T_static is linear (EA eps) or the piecewise-linear 'curve' through (0,0) and the given (strain, tension) points, extrapolated with
// the last slope. The relaxed (static) stiffness is the curve; the instantaneous (high-rate) stiffness adds sum K_k, so the dynamic
// stiffness exceeds the static one as measured for polyester and nylon ropes. This is a model form, not a fit: the parameters must come
// from tests of the rope concerned. It has no hysteresis in the static curve, no load-history (mean-load) dependence of the dynamic
// stiffness and no permanent set.
struct RopeLaw {
    std::vector<std::array<double, 2>> curve;   // (strain, tension [N]), strictly increasing, slopes > 0; empty = linear
    std::vector<MaxwellBranch> branches;

    bool nonlinear() const { return !curve.empty(); }
    bool viscoelastic() const { return !branches.empty(); }

    void validate() const {
        double e0 = 0.0, t0 = 0.0;
        for (const auto& p : curve) {
            if (!(p[0] > e0) || !(p[1] > t0)) throw std::invalid_argument("rope tension curve: strain and tension must be strictly increasing and positive");
            e0 = p[0]; t0 = p[1];
        }
        for (const MaxwellBranch& b : branches)
            if (!(b.K > 0.0) || !(b.tau > 0.0)) throw std::invalid_argument("rope Maxwell branch: K and tau must be > 0");
    }
    // Tangent stiffness of the static curve at the given strain [N].
    double staticSlope(double EA, double eps) const {
        if (curve.empty()) return EA;
        double e0 = 0.0, t0 = 0.0;
        for (size_t i = 0; i < curve.size(); ++i) {
            if (eps <= curve[i][0] || i + 1 == curve.size()) return (curve[i][1] - t0) / (curve[i][0] - e0);
            e0 = curve[i][0]; t0 = curve[i][1];
        }
        return EA;
    }
    double staticTension(double EA, double eps) const {
        if (eps <= 0.0) return 0.0;
        if (curve.empty()) return EA * eps;
        double e0 = 0.0, t0 = 0.0;
        for (size_t i = 0; i < curve.size(); ++i) {
            const double s = (curve[i][1] - t0) / (curve[i][0] - e0);
            if (eps <= curve[i][0] || i + 1 == curve.size()) return t0 + s * (eps - e0);
            e0 = curve[i][0]; t0 = curve[i][1];
        }
        return EA * eps;
    }
    // Strain energy of the static curve per unit unstretched length scaled by l0: integral of T_static d eps.
    double staticEnergy(double EA, double eps) const {
        if (eps <= 0.0) return 0.0;
        if (curve.empty()) return 0.5 * EA * eps * eps;
        double e0 = 0.0, t0 = 0.0, E = 0.0;
        for (size_t i = 0; i < curve.size(); ++i) {
            const double s = (curve[i][1] - t0) / (curve[i][0] - e0);
            const double e1 = (eps <= curve[i][0] || i + 1 == curve.size()) ? eps : curve[i][0];
            E += t0 * (e1 - e0) + 0.5 * s * (e1 - e0) * (e1 - e0);
            if (e1 == eps) return E;
            e0 = curve[i][0]; t0 = curve[i][1];
        }
        return E;
    }
    // Largest static tangent stiffness (governs the quasi-static relaxation) and the largest instantaneous one (governs the explicit step).
    double maxStaticStiffness(double EA) const {
        if (curve.empty()) return EA;
        double m = 0.0, e0 = 0.0, t0 = 0.0;
        for (const auto& p : curve) { m = std::max(m, (p[1] - t0) / (p[0] - e0)); e0 = p[0]; t0 = p[1]; }
        return m;
    }
    double maxStiffness(double EA) const {
        double m = maxStaticStiffness(EA);
        for (const MaxwellBranch& b : branches) m += b.K;
        return m;
    }
    // Exact update of one internal strain for a strain that varies linearly from eps0 to eps1 over dt.
    static double advanceAlpha(double alpha, double eps0, double eps1, double dt, double tau) {
        const double rate = (eps1 - eps0) / dt, ex = std::exp(-dt / tau);
        return eps1 - rate * tau + (alpha - eps0 + rate * tau) * ex;
    }
};

}  // namespace mooring
