// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
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
    if (p.planar) for (Vec3& q : r_) q.y = anchor.y;
}

double LumpedMassCable::nodeMass(int i) const {
    const double m = p_.m_l * p_.l0();
    return (i == 0 || i == p_.N) ? 0.5 * m : m;
}

double LumpedMassCable::segmentTension(const std::vector<Vec3>& r, int seg) const {
    const double eps = norm(r[seg + 1] - r[seg]) / p_.l0() - 1.0;
    return eps > 0.0 ? p_.EA * eps : 0.0;
}

double LumpedMassCable::submergedFraction(const std::vector<Vec3>& r, int i) const {
    if (!env_.hydro) return 1.0;
    // Linear blend over one segment length around the still-water level (avoids a force step).
    return std::min(1.0, std::max(0.0, 0.5 + (env_.surfaceZ - r[i].z) / p_.l0()));
}

double LumpedMassCable::nodeWeight(const std::vector<Vec3>& r, int i) const {
    const double half = (i == 0 || i == p_.N) ? 0.5 : 1.0;
    const double phi = submergedFraction(r, i);
    const double wl = env_.hydro ? phi * p_.w + (1.0 - phi) * p_.dryWeight() : p_.w;
    return wl * p_.l0() * half;
}

void LumpedMassCable::computeForces(const std::vector<Vec3>& r, const std::vector<Vec3>& v, double t,
                                    std::vector<Vec3>& f, std::vector<double>* ca,
                                    std::vector<Vec3>* aw) const {
    const int N = p_.N;
    f.assign(r.size(), Vec3());
    std::vector<double> len(N);
    for (int i = 0; i < N; ++i) {
        const Vec3 d = r[i + 1] - r[i];
        len[i] = norm(d);
        const double eps = len[i] / p_.l0() - 1.0;
        if (eps <= 0.0 || len[i] <= 0.0) { ++stats_.slackSegmentEvals; continue; }   // Eq. 3.25: T = 0
        double T = p_.EA * eps;
        if (p_.c_int != 0.0) {
            const double epsDot = dot(d, v[i + 1] - v[i]) / (len[i] * p_.l0());
            T += p_.c_int * epsDot;
            if (T < 0.0) { T = 0.0; ++stats_.clippedTensionEvals; }
        }
        const Vec3 F = d * (T / len[i]);
        f[i] += F;
        f[i + 1] -= F;
    }
    if (ca) ca->assign(r.size(), 0.0);
    if (aw) aw->assign(r.size(), Vec3());
    const double wSub = p_.w;
    for (int i = 0; i <= N; ++i) {
        const bool end = (i == 0 || i == N);
        f[i].z -= nodeWeight(r, i);
        if (env_.hydro) {
            const double phi = submergedFraction(r, i);
            if (phi > 0.0) {
                // Tangent from neighbours; stretched tributary length carries the (1+eps) of Eqs. 3.27-3.29.
                const Vec3 dt = i == 0 ? r[1] - r[0] : (i == N ? r[N] - r[N - 1] : r[i + 1] - r[i - 1]);
                const double ndt = norm(dt);
                const Vec3 tg = ndt > 0.0 ? dt / ndt : Vec3(1, 0, 0);
                const double Ls = 0.5 * ((i > 0 ? len[i - 1] : 0.0) + (i < N ? len[i] : 0.0));
                Vec3 vw, a_w;
                if (env_.water) env_.water(r[i], t, vw, a_w);
                f[i] += morisonDrag(vw - v[i], tg, p_.Cdt, p_.Cdn, env_.rho_w, p_.D0, Ls) * phi;
                if (ca) (*ca)[i] = phi * p_.Cm * env_.rho_w * p_.nominalArea() * Ls;
                if (aw) (*aw)[i] = a_w;
            }
        }
        if (env_.seabed) {
            const double pen = env_.seabedZ - r[i].z;
            if (pen > 0.0) {
                ++stats_.soilContactEvals;
                SoilParams sp = p_.soil;
                sp.Ks = std::min(sp.Ks, ksCap_);
                f[i] += seabedForce(sp, p_.D1, p_.m_l, wSub, pen, v[i], p_.l0() * (end ? 0.5 : 1.0));
            }
        }
    }
    if (p_.planar) for (Vec3& q : f) q.y = 0.0;      // 2D mode: no out-of-plane force
}

Vec3 LumpedMassCable::endTension(bool top) const {
    const int seg = top ? p_.N - 1 : 0;
    const Vec3 d = top ? r_[p_.N - 1] - r_[p_.N] : r_[1] - r_[0];
    const double len = norm(d);
    return len > 0 ? d * (segmentTension(r_, seg) / len) : Vec3();
}

double LumpedMassCable::endSegmentTension(bool top) const {
    const int seg = top ? p_.N - 1 : 0;
    const Vec3 d = r_[seg + 1] - r_[seg];
    const double len = norm(d), eps = len / p_.l0() - 1.0;
    if (eps <= 0.0 || len <= 0.0) return 0.0;
    double T = p_.EA * eps;
    if (p_.c_int != 0.0 && v_.size() == r_.size())
        T = std::max(0.0, T + p_.c_int * dot(d, v_[seg + 1] - v_[seg]) / (len * p_.l0()));
    return T;
}

bool LumpedMassCable::initTouchdownCatenary() {
    const Vec3 a0 = r_[0], a1 = r_[p_.N];
    const Vec3 dh(a1.x - a0.x, a1.y - a0.y, 0.0);
    const double D = norm(dh), h = a1.z - a0.z, L = p_.L;
    if (D <= 0.0 || h <= 0.0 || L <= std::sqrt(D * D + h * h)) return false;
    // g(a) = horizontal reach of (lying + catenary) - D, a = H/w; monotone increasing in a.
    auto g = [&](double a) { return (L - std::sqrt(h * h + 2.0 * a * h)) + a * std::acosh(1.0 + h / a) - D; };
    double lo = 1e-3, hi = 1e9;
    if (g(lo) > 0.0 || g(hi) < 0.0) return false;
    for (int i = 0; i < 300; ++i) { const double mid = 0.5 * (lo + hi); (g(mid) > 0.0 ? hi : lo) = mid; }
    const double a = 0.5 * (lo + hi), susp = std::sqrt(h * h + 2.0 * a * h), lying = L - susp;
    const Vec3 e = dh / D;
    for (int i = 0; i <= p_.N; ++i) {
        const double s = i * p_.l0();
        if (s <= lying) { r_[i] = a0 + e * s; }
        else {
            const double sg = s - lying;
            r_[i] = a0 + e * (lying + a * std::asinh(sg / a)) + Vec3(0, 0, std::sqrt(a * a + sg * sg) - a);
        }
    }
    r_[p_.N] = a1;
    return true;
}

Vec3 LumpedMassCable::endForce(bool top) const {
    std::vector<Vec3> f;
    const std::vector<Vec3> vz(r_.size());
    computeForces(r_, v_.size() == r_.size() ? v_ : vz, t_, f);
    return top ? f[p_.N] : f[0];
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
    const std::vector<Vec3> vz(n);                // fictitious velocities must not feed c_int damping
    const double tol = opt.forceTol * std::max(p_.w, 1e-12) * p_.l0();
    double keOld = 0.0;
    // Soil stiffness is capped so that the node contact stiffness equals the axial one (EA/l0); the
    // equilibrium penetration becomes w l0^2/EA instead of w/(Ks D1) - both negligible - and the
    // relaxation stays stable. Dynamic runs use the true Ks.
    ksCap_ = p_.D1 > 0.0 ? p_.EA / (p_.l0() * p_.l0() * p_.D1) : 1e300;
    for (long k = 0; k < opt.maxSteps; ++k) {
        computeForces(r_, vz, t_, f);
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
    ksCap_ = 1e300;
    return res;
}

// ---------------------------------------------------------------------------------------------
// Dynamics
// ---------------------------------------------------------------------------------------------

double LumpedMassCable::stableDt() const {
    if (dyn_.dt > 0.0) return dyn_.dt;
    double dt = dyn_.cfl * p_.l0() / p_.waveSpeed();                      // axial wave CFL
    // Explicit internal damping (c_int/l0 per segment, up to 4 c_int/(m_l l0^2) per node) adds its own limit.
    if (p_.c_int > 0.0) dt = std::min(dt, dyn_.cfl * p_.m_l * p_.l0() * p_.l0() / (2.0 * p_.c_int));
    if (env_.seabed && p_.soil.Ks > 0.0 && p_.D1 > 0.0) {
        // Contact oscillator: omega^2 = Ks D1 / m_l (independent of l0); largest eigenvalue magnitude
        // of the damped spring is omega (zeta + sqrt(zeta^2 - 1)) for zeta >= 1, omega otherwise.
        const double om = std::sqrt(p_.soil.Ks * p_.D1 / p_.m_l), z = p_.soil.zeta;
        const double lam = om * (z > 1.0 ? z + std::sqrt(z * z - 1.0) : 1.0);
        dt = std::min(dt, dyn_.cfl * 2.0 / lam);
    }
    return dt;
}

void LumpedMassCable::setInitialState(const std::vector<Vec3>& r, const std::vector<Vec3>& v, double t0) {
    if (r.size() != r_.size() || v.size() != r_.size()) throw std::invalid_argument("state size mismatch");
    r_ = r; v_ = v; t_ = t0;
    if (p_.planar) for (size_t i = 0; i < r_.size(); ++i) { r_[i].y = r_[0].y; v_[i].y = 0.0; }
}

void LumpedMassCable::acceleration(std::vector<Vec3>& r, std::vector<Vec3>& v, double t,
                                   std::vector<Vec3>& a) const {
    const int N = p_.N;
    r[0] = r_[0]; v[0] = Vec3();                   // anchor fixed
    if (top_) top_(t, r[N], v[N]);                 // otherwise the top end is held where it is
    std::vector<double> ca;
    std::vector<Vec3> aw;
    std::vector<Vec3>& f = a;
    computeForces(r, v, t, f, &ca, &aw);
    for (int i = 1; i < N; ++i) {
        const double m = nodeMass(i);
        if (ca[i] > 0.0) {
            // Added mass: M = m I + ca (I - t t^T), RHS carries ca (I - t t^T) a_w (Eq. 3.27 with a_rel = a_w - a).
            const Vec3 dt = r[i + 1] - r[i - 1];
            const double ndt = norm(dt);
            const Vec3 tg = ndt > 0.0 ? dt / ndt : Vec3(1, 0, 0);
            const Vec3 awn = aw[i] - tg * dot(aw[i], tg);
            f[i] = addedMassSolve(m, ca[i], tg, f[i] + awn * ca[i]);
        } else {
            f[i] = f[i] / m;
        }
    }
    f[0] = Vec3(); f[N] = Vec3();
}

void LumpedMassCable::step(double dt) {
    const int n = p_.N + 1, N = p_.N;
    if (v_.size() != r_.size()) v_.assign(n, Vec3());
    if (top_) top_(t_, r_[N], v_[N]);
    std::vector<Vec3> a(n);
    if (dyn_.scheme == Scheme::Verlet) {
        acceleration(r_, v_, t_, a);
        std::vector<Vec3> vh = v_;
        for (int i = 1; i < N; ++i) { vh[i] += a[i] * (0.5 * dt); r_[i] += vh[i] * dt; }
        acceleration(r_, vh, t_ + dt, a);
        for (int i = 1; i < N; ++i) v_[i] = vh[i] + a[i] * (0.5 * dt);
    } else {
        std::vector<Vec3> r1 = r_, v1 = v_, k1v(n), k2v(n), k3v(n), k4v(n), k2r(n), k3r(n), k4r(n);
        std::vector<Vec3> rs(n), vs(n);
        // stage 1
        acceleration(r1, v1, t_, k1v);
        const std::vector<Vec3>& k1r = v1;
        // stage 2
        rs = r_; vs = v_;
        for (int i = 1; i < N; ++i) { rs[i] = r_[i] + k1r[i] * (0.5 * dt); vs[i] = v_[i] + k1v[i] * (0.5 * dt); }
        acceleration(rs, vs, t_ + 0.5 * dt, k2v); k2r = vs;
        // stage 3
        for (int i = 1; i < N; ++i) { rs[i] = r_[i] + k2r[i] * (0.5 * dt); vs[i] = v_[i] + k2v[i] * (0.5 * dt); }
        acceleration(rs, vs, t_ + 0.5 * dt, k3v); k3r = vs;
        // stage 4
        for (int i = 1; i < N; ++i) { rs[i] = r_[i] + k3r[i] * dt; vs[i] = v_[i] + k3v[i] * dt; }
        acceleration(rs, vs, t_ + dt, k4v); k4r = vs;
        for (int i = 1; i < N; ++i) {
            r_[i] += (k1r[i] + k2r[i] * 2.0 + k3r[i] * 2.0 + k4r[i]) * (dt / 6.0);
            v_[i] += (k1v[i] + k2v[i] * 2.0 + k3v[i] * 2.0 + k4v[i]) * (dt / 6.0);
        }
    }
    t_ += dt;
    if (top_) top_(t_, r_[N], v_[N]);
    ++stats_.steps;
    if (observer_) observer_(*this);
}

void LumpedMassCable::advanceTo(double tEnd) {
    const double dt0 = stableDt();
    stats_.dtUsed = dt0;
    if (v_.size() != r_.size()) v_.assign(r_.size(), Vec3());
    while (t_ < tEnd - 1e-12 * dt0) step(std::min(dt0, tEnd - t_));
}

double LumpedMassCable::energy() const {
    double E = 0.0;
    for (int i = 0; i <= p_.N; ++i) {
        if (!v_.empty()) E += 0.5 * nodeMass(i) * dot(v_[i], v_[i]);
        E += p_.w * p_.l0() * ((i == 0 || i == p_.N) ? 0.5 : 1.0) * r_[i].z;
    }
    for (int i = 0; i < p_.N; ++i) {
        const double eps = norm(r_[i + 1] - r_[i]) / p_.l0() - 1.0;
        if (eps > 0.0) E += 0.5 * p_.EA * eps * eps * p_.l0();
    }
    return E;
}

Vec3 LumpedMassCable::forceOnBody(const Vec3& pos, const Vec3& vel, double t) {
    if (!fairInit_) {
        fairPrevPos_ = pos; fairPrevVel_ = vel; fairPrevT_ = t; fairInit_ = true;
        r_[p_.N] = pos; if (v_.size() == r_.size()) v_[p_.N] = vel;
        t_ = t;
    } else if (t > fairPrevT_) {
        const Vec3 p0 = fairPrevPos_, p1 = pos;
        const double t0 = fairPrevT_, t1 = t;
        top_ = [p0, p1, t0, t1](double tt, Vec3& r, Vec3& v) {   // linear interpolation of fairlead motion
            const double u = (tt - t0) / (t1 - t0);
            r = p0 + (p1 - p0) * u;
            v = (p1 - p0) / (t1 - t0);
        };
        advanceTo(t);
        fairPrevPos_ = pos; fairPrevVel_ = vel; fairPrevT_ = t;
    }
    return endForce(true);
}

}  // namespace mooring
