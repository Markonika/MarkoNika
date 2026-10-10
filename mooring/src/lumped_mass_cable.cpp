// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#include "mooring/lumped_mass_cable.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace mooring {

double CableParams::waveSpeed() const {
    if (sections.empty()) return std::sqrt(rope.maxStiffness(EA) / m_l);
    double c = 0.0;
    for (const LineSection& q : sections) c = std::max(c, std::sqrt(q.rope.maxStiffness(q.EA) / q.m_l));
    return c;
}

LumpedMassCable::LumpedMassCable(const CableParams& p, const Vec3& anchor, const Vec3& fairlead) : p_(p) {
    std::vector<LineSection> secs = p.sections;
    if (secs.empty()) {                                    // uniform line: one section from the top-level fields
        if (p.rope.nonlinear() && p_.EA <= 0) p_.EA = p.rope.curve[0][1] / p.rope.curve[0][0];   // initial slope as reference EA
        if (p_.N < 2 || p_.L <= 0 || p_.EA <= 0 || p_.m_l <= 0) throw std::invalid_argument("invalid CableParams");
        LineSection q; q.length = p_.L; q.segments = p_.N; q.EA = p_.EA; q.m_l = p_.m_l; q.w = p_.w; q.c_int = p_.c_int;
        q.D0 = p_.D0; q.D1 = p_.D1; q.A1 = p_.A1; q.Cm = p_.Cm; q.Cdt = p_.Cdt; q.Cdn = p_.Cdn; q.rope = p_.rope;
        secs.push_back(q);
    } else {
        double L = 0, mass = 0, wsum = 0, comp = 0; int N = 0;
        for (LineSection& q : secs) {
            if (q.rope.nonlinear() && q.EA <= 0) q.EA = q.rope.curve[0][1] / q.rope.curve[0][0];
            if (q.segments < 1 || q.length <= 0 || q.EA <= 0 || q.m_l <= 0) throw std::invalid_argument("invalid LineSection");
            L += q.length; N += q.segments; mass += q.m_l * q.length; wsum += q.w * q.length; comp += q.length / q.EA;
        }
        if (N < 2) throw std::invalid_argument("invalid CableParams: fewer than 2 segments");
        p_.L = L; p_.N = N; p_.m_l = mass / L; p_.w = wsum / L; p_.EA = L / comp;   // reporting values only
    }
    const int N = p_.N;
    for (size_t k = 0; k < secs.size(); ++k) {
        secs[k].rope.validate();
        ropes_.push_back(secs[k].rope);
        const LineSection& q = secs[k];
        for (int j = 0; j < q.segments; ++j)
            seg_.push_back({q.length / q.segments, q.EA, q.m_l, q.w, q.m_l * p_.g, q.D0, q.D1, q.A1 > 0.0 ? q.A1 : 0.7853981633974483 * q.D0 * q.D0,
                            q.Cm, q.Cdt, q.Cdn, q.c_int, int(k), 0});
    }
    size_t off = 0;
    for (SegProps& g : seg_) { g.aoff = off; off += ropes_[g.sec].branches.size(); }
    alpha_.assign(off, 0.0);
    s_.assign(N + 1, 0.0);
    {   // arc length = start of section + k l0 (identical to i l0 for a uniform line)
        int node = 0; double start = 0.0;
        for (const LineSection& q : secs) {
            const double l0 = q.length / q.segments;
            for (int k = 1; k <= q.segments; ++k) s_[node + k] = start + k * l0;
            node += q.segments; start += q.length;
        }
    }
    same_.assign(N + 1, 0);
    for (int i = 1; i < N; ++i) {
        const SegProps &a = seg_[i - 1], &b = seg_[i];
        same_[i] = a.l0 == b.l0 && a.EA == b.EA && a.ml == b.ml && a.w == b.w && a.D0 == b.D0 && a.D1 == b.D1 && a.A1 == b.A1 &&
                   a.Cm == b.Cm && a.Cdt == b.Cdt && a.Cdn == b.Cdn && a.cint == b.cint && a.sec == b.sec;
    }
    r_.resize(N + 1);
    // Initial shape: straight chord with a parabolic sag matching the slack length.
    const Vec3 d = fairlead - anchor;
    const double chord = norm(d);
    const double slack = std::max(p_.L - chord, 0.0);
    const double sag = std::sqrt(3.0 * chord * slack / 8.0);
    for (int i = 0; i <= N; ++i) {
        const double u = double(i) / N;
        r_[i] = anchor + d * u + Vec3(0, 0, -4.0 * sag * u * (1.0 - u));
    }
    if (p.planar) for (Vec3& q : r_) q.y = anchor.y;
}

int LumpedMassCable::nodeAtArclength(double s) const {
    int best = 0;
    for (int i = 1; i <= p_.N; ++i) if (std::fabs(s_[i] - s) < std::fabs(s_[best] - s)) best = i;
    return best;
}

double LumpedMassCable::nodeL0(int i) const {
    if (i == 0) return seg_[0].l0;
    if (i == p_.N) return seg_[p_.N - 1].l0;
    return 0.5 * (seg_[i - 1].l0 + seg_[i].l0);
}

double LumpedMassCable::nodeMass(int i) const {
    const double a = i > 0 ? seg_[i - 1].ml * seg_[i - 1].l0 : 0.0, b = i < p_.N ? seg_[i].ml * seg_[i].l0 : 0.0;
    return 0.5 * (a + b);
}

double LumpedMassCable::segmentTension(const std::vector<Vec3>& r, int seg) const {
    const double eps = norm(r[seg + 1] - r[seg]) / seg_[seg].l0 - 1.0;
    return eps > 0.0 ? std::max(0.0, segTension(seg, eps, 0.0)) : 0.0;
}

double LumpedMassCable::segTension(int seg, double eps, double epsDot) const {
    const SegProps& g = seg_[seg];
    const RopeLaw& rl = ropes_[g.sec];
    double T = rl.staticTension(g.EA, eps);
    if (!staticMode_) {
        const size_t nb = rl.branches.size();
        for (size_t k = 0; k < nb; ++k) T += rl.branches[k].K * (eps - alpha_[g.aoff + k]);
    }
    return T + g.cint * epsDot;
}

void LumpedMassCable::resetInternalState() {
    for (int s = 0; s < p_.N; ++s) {
        const size_t nb = ropes_[seg_[s].sec].branches.size();
        if (nb == 0) continue;
        const double eps = norm(r_[s + 1] - r_[s]) / seg_[s].l0 - 1.0;
        for (size_t k = 0; k < nb; ++k) alpha_[seg_[s].aoff + k] = eps;
    }
}

double LumpedMassCable::submergedFraction(const std::vector<Vec3>& r, int i, double t) const {
    if (!env_.hydro) return 1.0;
    double zs = env_.surfaceZ;
    if (env_.elevation && !staticMode_) zs += env_.elevation(r[i], t);
    // Linear blend over one segment length around the still-water level (avoids a force step).
    return std::min(1.0, std::max(0.0, 0.5 + (zs - r[i].z) / nodeL0(i)));
}

void LumpedMassCable::addPointElement(const PointElement& pe) {
    if (pe.node < 0 || pe.node > p_.N) throw std::invalid_argument("point element node out of range");
    if (pe.mass < 0.0 || pe.volume < 0.0) throw std::invalid_argument("point element mass/volume must be >= 0");
    points_.push_back(pe);
}

double LumpedMassCable::pointNetWeight(const std::vector<Vec3>& r, int i, double t) const {
    double w = 0.0;
    for (const PointElement& pe : points_) {
        if (pe.node != i) continue;
        w += pe.mass * p_.g;
        if (env_.hydro) w -= submergedFraction(r, i, t) * env_.rho_w * p_.g * pe.volume;
    }
    return w;
}

double LumpedMassCable::nodeWeight(const std::vector<Vec3>& r, int i, double t) const {
    const double phi = submergedFraction(r, i, t);
    auto part = [&](int j) {
        const SegProps& g = seg_[j];
        const double wl = env_.hydro ? phi * g.w + (1.0 - phi) * g.wd : g.w;
        return wl * g.l0;
    };
    const double a = i > 0 ? part(i - 1) : 0.0, b = i < p_.N ? part(i) : 0.0;
    return 0.5 * (a + b);
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
        const double eps = len[i] / seg_[i].l0 - 1.0;
        if (eps <= 0.0 || len[i] <= 0.0) { ++stats_.slackSegmentEvals; continue; }   // Eq. 3.25: T = 0
        const double epsDot = (seg_[i].cint != 0.0 && !skipDamping_) ? dot(d, v[i + 1] - v[i]) / (len[i] * seg_[i].l0) : 0.0;
        double T = segTension(i, eps, epsDot);
        if (T < 0.0) { T = 0.0; ++stats_.clippedTensionEvals; }
        const Vec3 F = d * (T / len[i]);
        f[i] += F;
        f[i + 1] -= F;
    }
    if (ca) ca->assign(r.size(), 0.0);
    if (aw) aw->assign(r.size(), Vec3());
    for (int i = 0; i <= N; ++i) {
        f[i].z -= nodeWeight(r, i, t);
        if (env_.hydro) {
            const double phi = submergedFraction(r, i, t);
            if (phi > 0.0) {
                // Tangent from neighbours; stretched tributary length carries the (1+eps) of Eqs. 3.27-3.29.
                const Vec3 dt = i == 0 ? r[1] - r[0] : (i == N ? r[N] - r[N - 1] : r[i + 1] - r[i - 1]);
                const double ndt = norm(dt);
                const Vec3 tg = ndt > 0.0 ? dt / ndt : Vec3(1, 0, 0);
                Vec3 vw, a_w;
                if (env_.water) env_.water(r[i], t, vw, a_w);
                double caSum = 0.0;
                // Tributary length of each adjacent segment (half of it, or the sum when both neighbours share their properties).
                auto addSeg = [&](int j, double Ls) {
                    const SegProps& g = seg_[j];
                    f[i] += morisonDrag(vw - v[i], tg, g.Cdt, g.Cdn, env_.rho_w, g.D0, Ls) * phi;
                    caSum += phi * g.Cm * env_.rho_w * g.A1 * Ls;
                };
                if (i > 0 && i < N && same_[i]) addSeg(i, 0.5 * (len[i - 1] + len[i]));
                else { if (i > 0) addSeg(i - 1, 0.5 * len[i - 1]); if (i < N) addSeg(i, 0.5 * len[i]); }
                if (ca) (*ca)[i] = caSum;
                if (aw) (*aw)[i] = a_w;
            }
        }
        if (env_.seabed) {
            const double pen = env_.seabedZ - r[i].z;
            if (pen > 0.0) {
                ++stats_.soilContactEvals;
                auto addSoil = [&](int j, double Ls) {
                    const SegProps& g = seg_[j];
                    SoilParams sp = p_.soil;
                    if (capSoil_ && g.D1 > 0.0) sp.Ks = std::min(sp.Ks, ropes_[g.sec].maxStaticStiffness(g.EA) / (g.l0 * g.l0 * g.D1));
                    f[i] += seabedForce(sp, g.D1, g.ml, g.w, pen, v[i], Ls);
                };
                if (i > 0 && i < N && same_[i]) addSoil(i, seg_[i].l0);
                else { if (i > 0) addSoil(i - 1, 0.5 * seg_[i - 1].l0); if (i < N) addSoil(i, 0.5 * seg_[i].l0); }
            }
        }
    }
    for (const PointElement& pe : points_) {
        const int i = pe.node;
        f[i].z -= pe.mass * p_.g;
        if (!env_.hydro) continue;
        const double phi = submergedFraction(r, i, t);
        if (phi <= 0.0) continue;
        f[i].z += phi * env_.rho_w * p_.g * pe.volume;
        const bool aniso = pe.CdT >= 0.0 || pe.CdN >= 0.0;
        if (aniso) {
            Vec3 vw, a_w;
            if (env_.water) env_.water(r[i], t, vw, a_w);
            const Vec3 vr = vw - v[i];
            const Vec3 dtv = i == 0 ? r[1] - r[0] : (i == N ? r[N] - r[N - 1] : r[i + 1] - r[i - 1]);
            const double ndt = norm(dtv);
            const Vec3 tg = ndt > 0.0 ? dtv / ndt : Vec3(1, 0, 0);
            const Vec3 vt = tg * dot(vr, tg), vn = vr - vt;
            f[i] += (vt * (0.5 * env_.rho_w * std::max(pe.CdT, 0.0) * pe.areaT * norm(vt)) +
                     vn * (0.5 * env_.rho_w * std::max(pe.CdN, 0.0) * pe.areaN * norm(vn))) * phi;
        } else if (pe.Cd > 0.0 && pe.area > 0.0) {
            Vec3 vw, a_w;
            if (env_.water) env_.water(r[i], t, vw, a_w);
            const Vec3 vr = vw - v[i];
            f[i] += vr * (0.5 * env_.rho_w * pe.Cd * pe.area * norm(vr) * phi);
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
    const double len = norm(d), eps = len / seg_[seg].l0 - 1.0;
    if (eps <= 0.0 || len <= 0.0) return 0.0;
    const double epsDot = (seg_[seg].cint != 0.0 && v_.size() == r_.size()) ? dot(d, v_[seg + 1] - v_[seg]) / (len * seg_[seg].l0) : 0.0;
    return std::max(0.0, segTension(seg, eps, epsDot));
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
        const double s = s_[i];
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
    // The relaxed (static) stiffness per segment k_j = EA_j / l0_j governs the quasi-static relaxation; each node gets the larger of its two.
    std::vector<double> kseg(p_.N), mFict(n, 0.0);
    double EAmax = 0.0, wl0 = 0.0;
    for (int j = 0; j < p_.N; ++j) {
        const SegProps& g = seg_[j];
        const double EAs = ropes_[g.sec].maxStaticStiffness(g.EA);
        kseg[j] = EAs / g.l0; EAmax = std::max(EAmax, EAs);
        wl0 = std::max({wl0, g.w * g.l0, g.ml * p_.g * g.l0});
    }
    for (int i = 1; i < p_.N; ++i) {
        const int j = kseg[i - 1] >= kseg[i] ? i - 1 : i;
        const SegProps& g = seg_[j];
        mFict[i] = opt.massFactor * res.dt * res.dt * ropes_[g.sec].maxStaticStiffness(g.EA) / g.l0;   // same operation order as the uniform-line code
    }
    staticMode_ = true;
    std::vector<Vec3> v(n), f(n);
    const std::vector<Vec3> vz(n);                // fictitious velocities must not feed c_int damping
    // Force scale for the convergence test: nodal weight, or (taut, weightless lines) a small fraction of the axial stiffness.
    const double fScale = std::max(wl0, 1e-7 * EAmax);
    const double tol = opt.forceTol * fScale;
    double keOld = 0.0;
    // Soil stiffness is capped so that the node contact stiffness equals the axial one (EA/l0); the
    // equilibrium penetration becomes w l0^2/EA instead of w/(Ks D1) - both negligible - and the
    // relaxation stays stable. Dynamic runs use the true Ks.
    capSoil_ = true;
    for (long k = 0; k < opt.maxSteps; ++k) {
        computeForces(r_, vz, t_, f);
        double ke = 0.0, rmax = 0.0;
        for (int i = 1; i < p_.N; ++i) {          // end nodes fixed
            v[i] += f[i] * (res.dt / mFict[i]);
            r_[i] += v[i] * res.dt;
            ke += 0.5 * mFict[i] * dot(v[i], v[i]);
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
    capSoil_ = false;
    staticMode_ = false;
    resetInternalState();
    v_.assign(r_.size(), Vec3());                  // the line is left at rest
    return res;
}

// ---------------------------------------------------------------------------------------------
// Dynamics
// ---------------------------------------------------------------------------------------------

double LumpedMassCable::stableDt() const {
    if (dyn_.dt > 0.0) return dyn_.dt;
    double dt = 1e300;
    for (const SegProps& g : seg_) {
        const double c = std::sqrt(ropes_[g.sec].maxStiffness(g.EA) / g.ml);
        dt = std::min(dt, dyn_.cfl * g.l0 / c);                                  // axial wave CFL
        // Explicit internal damping (c_int/l0 per segment, up to 4 c_int/(m_l l0^2) per node) adds its own limit.
        if (g.cint > 0.0 && !dyn_.implicitDamping) dt = std::min(dt, dyn_.cfl * g.ml * g.l0 * g.l0 / (2.0 * g.cint));
        if (env_.seabed && p_.soil.Ks > 0.0 && g.D1 > 0.0) {
            // Contact oscillator: omega^2 = Ks D1 / m_l (independent of l0); largest eigenvalue magnitude
            // of the damped spring is omega (zeta + sqrt(zeta^2 - 1)) for zeta >= 1, omega otherwise.
            const double om = std::sqrt(p_.soil.Ks * g.D1 / g.ml), z = p_.soil.zeta;
            const double lam = om * (z > 1.0 ? z + std::sqrt(z * z - 1.0) : 1.0);
            dt = std::min(dt, dyn_.cfl * 2.0 / lam);
        }
    }
    return dt;
}

void LumpedMassCable::setInitialState(const std::vector<Vec3>& r, const std::vector<Vec3>& v, double t0) {
    if (r.size() != r_.size() || v.size() != r_.size()) throw std::invalid_argument("state size mismatch");
    r_ = r; v_ = v; t_ = t0;
    if (p_.planar) for (size_t i = 0; i < r_.size(); ++i) { r_[i].y = r_[0].y; v_[i].y = 0.0; }
    resetInternalState();
}

void LumpedMassCable::acceleration(std::vector<Vec3>& r, std::vector<Vec3>& v, double t,
                                   std::vector<Vec3>& a) const {
    const int N = p_.N;
    r[0] = r_[0]; v[0] = Vec3();                   // anchor fixed
    if (top_) top_(t, r[N], v[N]);                 // otherwise the top end is held where it is
    std::vector<double> ca;
    std::vector<Vec3> aw;
    std::vector<Vec3>& f = a;
    skipDamping_ = dyn_.implicitDamping;
    computeForces(r, v, t, f, &ca, &aw);
    skipDamping_ = false;
    std::vector<double> mPt(N + 1, 0.0), caPt(N + 1, 0.0), maT(N + 1, 0.0), maN(N + 1, 0.0);   // point-element inertia and added masses
    for (const PointElement& pe : points_) {
        mPt[pe.node] += pe.mass;
        if (!env_.hydro) continue;
        const double phi = submergedFraction(r, pe.node, t);
        if (pe.Cm > 0.0) caPt[pe.node] += phi * pe.Cm * env_.rho_w * pe.volume;
        maT[pe.node] += phi * pe.maTan; maN[pe.node] += phi * pe.maNorm;
    }
    for (int i = 1; i < N; ++i) {
        const double m = nodeMass(i) + mPt[i] + caPt[i];
        if (caPt[i] > 0.0) f[i] += aw[i] * caPt[i];             // RHS of Cm rho V (a_w - a), isotropic
        if (ca[i] > 0.0 || maT[i] > 0.0 || maN[i] > 0.0) {
            // M = m I + (ca + maN) (I - t t^T) + maT t t^T; the right-hand side carries the water-acceleration parts.
            const Vec3 dt = r[i + 1] - r[i - 1];
            const double ndt = norm(dt);
            const Vec3 tg = ndt > 0.0 ? dt / ndt : Vec3(1, 0, 0);
            const Vec3 awt = tg * dot(aw[i], tg), awn = aw[i] - awt;
            f[i] = addedMassSolve2(m + maT[i], m + maN[i] + ca[i], tg, f[i] + awn * (ca[i] + maN[i]) + awt * maT[i]);
        } else {
            f[i] = f[i] / m;
        }
    }
    f[0] = Vec3(); f[N] = Vec3();
}

// Added mass [kg] of the line at node i (sum over the two adjacent segments of phi Cm rho A (stretched half length)).
double LumpedMassCable::nodeAddedMass(const std::vector<Vec3>& r, int i, double phi) const {
    const int N = p_.N;
    auto term = [&](int j, double Ls) { const SegProps& g = seg_[j]; return phi * g.Cm * env_.rho_w * g.A1 * Ls; };
    const double a = i > 0 ? norm(r[i] - r[i - 1]) : 0.0, b = i < N ? norm(r[i + 1] - r[i]) : 0.0;
    if (i > 0 && i < N && same_[i]) return term(i, 0.5 * (a + b));          // same expression as the force evaluation
    return (i > 0 ? term(i - 1, 0.5 * a) : 0.0) + (i < N ? term(i, 0.5 * b) : 0.0);
}

// Inertia of an interior node in the form M = mt t t^T + mn (I - t t^T), identical to the one solved in acceleration().
void LumpedMassCable::nodeInertia(const std::vector<Vec3>& r, double t, int i, double& mt, double& mn, Vec3& tg) const {
    double m = nodeMass(i), caPt = 0.0, maT = 0.0, maN = 0.0, ca = 0.0;
    for (const PointElement& pe : points_) {
        if (pe.node != i) continue;
        m += pe.mass;
        if (!env_.hydro) continue;
        const double phi = submergedFraction(r, i, t);
        if (pe.Cm > 0.0) caPt += phi * pe.Cm * env_.rho_w * pe.volume;
        maT += phi * pe.maTan; maN += phi * pe.maNorm;
    }
    m += caPt;
    if (env_.hydro) {
        const double phi = submergedFraction(r, i, t);
        if (phi > 0.0) {
            ca = nodeAddedMass(r, i, phi);
        }
    }
    const Vec3 dtv = r[i + 1] - r[i - 1];
    const double ndt = norm(dtv);
    tg = ndt > 0.0 ? dtv / ndt : Vec3(1, 0, 0);
    mt = m + maT; mn = m + maN + ca;
}

namespace {
struct Mat3 {
    double a[9]{};
    static Mat3 outer(const Vec3& u, double s) {
        Mat3 m; const double x[3] = {u.x, u.y, u.z};
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) m.a[3 * i + j] = s * x[i] * x[j];
        return m;
    }
    Mat3& operator+=(const Mat3& o) { for (int i = 0; i < 9; ++i) a[i] += o.a[i]; return *this; }
    Mat3& operator-=(const Mat3& o) { for (int i = 0; i < 9; ++i) a[i] -= o.a[i]; return *this; }
    Vec3 operator*(const Vec3& v) const {
        return Vec3(a[0] * v.x + a[1] * v.y + a[2] * v.z, a[3] * v.x + a[4] * v.y + a[5] * v.z, a[6] * v.x + a[7] * v.y + a[8] * v.z);
    }
    Mat3 operator*(const Mat3& o) const {
        Mat3 m;
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) m.a[3 * i + j] += a[3 * i + k] * o.a[3 * k + j];
        return m;
    }
    Mat3 inverse() const {
        const double c00 = a[4] * a[8] - a[5] * a[7], c01 = a[5] * a[6] - a[3] * a[8], c02 = a[3] * a[7] - a[4] * a[6];
        const double det = a[0] * c00 + a[1] * c01 + a[2] * c02;
        Mat3 m;
        m.a[0] = c00 / det; m.a[1] = (a[2] * a[7] - a[1] * a[8]) / det; m.a[2] = (a[1] * a[5] - a[2] * a[4]) / det;
        m.a[3] = c01 / det; m.a[4] = (a[0] * a[8] - a[2] * a[6]) / det; m.a[5] = (a[2] * a[3] - a[0] * a[5]) / det;
        m.a[6] = c02 / det; m.a[7] = (a[1] * a[6] - a[0] * a[7]) / det; m.a[8] = (a[0] * a[4] - a[1] * a[3]) / det;
        return m;
    }
};
}  // namespace

// Backward Euler for M dv/dt = -D(r) v, where the segment dashpots D = sum k P_j couple neighbouring nodes (block-tridiagonal 3x3 system,
// solved by block Thomas). k = c_int / l0 [N s/m] acts along the segment direction (T_damp = c_int d eps/dt); slack segments carry no damping.
// The anchor is fixed and the top end velocity is prescribed (taken from v_[N]).
void LumpedMassCable::implicitDamp(double h) {
    const int N = p_.N;
    std::vector<Mat3> P(N);
    for (int j = 0; j < N; ++j) {
        const Vec3 d = r_[j + 1] - r_[j];
        const double len = norm(d), eps = len / seg_[j].l0 - 1.0;
        if (eps > 0.0 && len > 0.0 && seg_[j].cint > 0.0) P[j] = Mat3::outer(d / len, seg_[j].cint / seg_[j].l0 * h);
    }
    const int n = N - 1;                               // interior nodes 1..N-1
    std::vector<Mat3> Cp(n);                           // A_i^{-1} U_i
    std::vector<Vec3> gp(n);                           // A_i^{-1} rhs_i (modified)
    for (int i = 1; i <= n; ++i) {
        double mt, mn; Vec3 tg;
        nodeInertia(r_, t_, i, mt, mn, tg);
        Mat3 M = Mat3::outer(tg, mt - mn);
        M.a[0] += mn; M.a[4] += mn; M.a[8] += mn;
        Mat3 A = M; A += P[i - 1]; A += P[i];
        Vec3 b = M * v_[i];
        if (i == n) b += P[i] * v_[N];                 // prescribed top velocity
        if (i > 1) {                                   // eliminate the sub-diagonal block L_i = -P[i-1]
            const Mat3 LC = P[i - 1] * Cp[i - 2];      // L_i * C_{i-1} = -P * C  =>  A -= L C  is  A += P C
            A += LC;
            b += P[i - 1] * gp[i - 2];
        }
        const Mat3 Ai = A.inverse();
        gp[i - 1] = Ai * b;
        if (i < n) { Mat3 U = P[i]; for (double& x : U.a) x = -x; Cp[i - 1] = Ai * U; }
    }
    Vec3 vnext;
    for (int i = n; i >= 1; --i) {
        vnext = i == n ? gp[i - 1] : gp[i - 1] - Cp[i - 1] * vnext;
        v_[i] = vnext;
    }
}

void LumpedMassCable::step(double dt) {
    const int n = p_.N + 1, N = p_.N;
    if (v_.size() != r_.size()) v_.assign(n, Vec3());
    if (top_) top_(t_, r_[N], v_[N]);
    std::vector<Vec3> a(n);
    // Maxwell internal strains are held at their start-of-step values during the (RK4 / Verlet) stages and advanced exactly afterwards
    // for a strain varying linearly over the step: first order in dt / tau (documented, docs/assumptions.md).
    const bool anyBranch = !alpha_.empty();
    std::vector<double> eps0;
    if (anyBranch) {
        eps0.resize(N);
        for (int s = 0; s < N; ++s) eps0[s] = norm(r_[s + 1] - r_[s]) / seg_[s].l0 - 1.0;
    }
    bool anyDamp = false;
    for (const SegProps& g : seg_) anyDamp = anyDamp || g.cint > 0.0;
    const bool impl = dyn_.implicitDamping && anyDamp;
    if (impl) implicitDamp(0.5 * dt);
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
    if (impl) implicitDamp(0.5 * dt);
    for (int s = 0; s < N && anyBranch; ++s) {
        const RopeLaw& rl = ropes_[seg_[s].sec];
        if (rl.branches.empty()) continue;
        const double eps1 = norm(r_[s + 1] - r_[s]) / seg_[s].l0 - 1.0;
        for (size_t k = 0; k < rl.branches.size(); ++k)
            alpha_[seg_[s].aoff + k] = RopeLaw::advanceAlpha(alpha_[seg_[s].aoff + k], eps0[s], eps1, dt, rl.branches[k].tau);
    }
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
        const double wa = i > 0 ? seg_[i - 1].w * seg_[i - 1].l0 : 0.0, wb = i < p_.N ? seg_[i].w * seg_[i].l0 : 0.0;
        E += 0.5 * (wa + wb) * r_[i].z;
    }
    for (int i = 0; i < p_.N; ++i) {
        const double eps = norm(r_[i + 1] - r_[i]) / seg_[i].l0 - 1.0;
        if (eps > 0.0) E += ropes_[seg_[i].sec].staticEnergy(seg_[i].EA, eps) * seg_[i].l0;
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

void LumpedMassCable::setFairlead(const Vec3& pos) {
    r_[p_.N] = pos;
    v_.assign(r_.size(), Vec3());
    resetInternalState();
    fairInit_ = false;                             // next forceOnBody() restarts the fairlead interpolation
}

Vec3 LumpedMassCable::staticForceOnBody(const Vec3& fairleadPos) {
    setFairlead(fairleadPos);
    RelaxOptions ro; ro.forceTol = 1e-7; ro.maxSteps = 5000000;   // 1e-7 of the force scale: ~1e-9 N here, above the round-off floor
    relaxStatic(ro);
    return endForce(true);
}

}  // namespace mooring
