// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#include "mooring/rigid_body.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace mooring {

Vec6 operator*(const Mat6& m, const Vec6& v) {
    Vec6 r{};
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) r[i] += m(i, j) * v[j];
    return r;
}
Mat6 operator+(const Mat6& a, const Mat6& b) { Mat6 r; for (int i = 0; i < 36; ++i) r.a[i] = a.a[i] + b.a[i]; return r; }
Mat6 operator*(const Mat6& m, double s) { Mat6 r; for (int i = 0; i < 36; ++i) r.a[i] = m.a[i] * s; return r; }

Vec6 solve6(Mat6 m, Vec6 b) {
    for (int c = 0; c < 6; ++c) {
        int piv = c;
        for (int r = c + 1; r < 6; ++r) if (std::fabs(m(r, c)) > std::fabs(m(piv, c))) piv = r;
        if (std::fabs(m(piv, c)) < 1e-300) throw std::runtime_error("solve6: singular matrix");
        if (piv != c) { for (int j = 0; j < 6; ++j) std::swap(m(c, j), m(piv, j)); std::swap(b[c], b[piv]); }
        for (int r = c + 1; r < 6; ++r) {
            const double f = m(r, c) / m(c, c);
            for (int j = c; j < 6; ++j) m(r, j) -= f * m(c, j);
            b[r] -= f * b[c];
        }
    }
    Vec6 x{};
    for (int i = 5; i >= 0; --i) {
        double s = b[i];
        for (int j = i + 1; j < 6; ++j) s -= m(i, j) * x[j];
        x[i] = s / m(i, i);
    }
    return x;
}

Mat3 rotationFromVector(const Vec3& th) {
    const double a = norm(th);
    double s1, s2;                                   // sin(a)/a and (1-cos(a))/a^2 with series for small a
    if (a < 1e-4) { s1 = 1.0 - a * a / 6.0; s2 = 0.5 - a * a / 24.0; }
    else { s1 = std::sin(a) / a; s2 = (1.0 - std::cos(a)) / (a * a); }
    const double x = th.x, y = th.y, z = th.z;
    // R = I + s1 K + s2 K^2, K = [th]x
    return {1.0 + s2 * (-y * y - z * z), -s1 * z + s2 * x * y,          s1 * y + s2 * x * z,
            s1 * z + s2 * x * y,         1.0 + s2 * (-x * x - z * z),  -s1 * x + s2 * y * z,
            -s1 * y + s2 * x * z,        s1 * x + s2 * y * z,           1.0 + s2 * (-x * x - y * y)};
}

Vec3 rotate(const Mat3& R, const Vec3& v) {
    return {R[0] * v.x + R[1] * v.y + R[2] * v.z, R[3] * v.x + R[4] * v.y + R[5] * v.z, R[6] * v.x + R[7] * v.y + R[8] * v.z};
}

RigidBody6DOF::RigidBody6DOF(const BodyParams& p) : xi{}, xiDot{}, p_(p) {
    Mat6 Mrb;
    for (int i = 0; i < 3; ++i) Mrb(i, i) = p.mass;
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) Mrb(3 + i, 3 + j) = p.inertia[3 * i + j];
    M_ = Mrb + p.A;
}

Vec3 RigidBody6DOF::pointVelocity(const Vec3& a) const {
    const Vec3 Ra = rotate(rotation(), a);
    return Vec3(xiDot[0], xiDot[1], xiDot[2]) + cross(Vec3(xiDot[3], xiDot[4], xiDot[5]), Ra);
}

void RigidBody6DOF::addPointLoad(const Vec3& Ra, const Vec3& f, Vec6& F) {
    const Vec3 m = cross(Ra, f);
    F[0] += f.x; F[1] += f.y; F[2] += f.z; F[3] += m.x; F[4] += m.y; F[5] += m.z;
}

int CoupledSystem::addLine(std::unique_ptr<CableModel> line, const Vec3& a, const std::string& name) {
    lines_.push_back({std::move(line), a, name, Vec3()});
    return int(lines_.size()) - 1;
}

Vec6 CoupledSystem::externalLoad(const Vec6& xi, const std::vector<Vec3>& lf) const {
    const BodyParams& p = body_.params();
    Vec6 F = p.F0;
    const Vec6 Cx = p.C * xi;
    for (int k = 0; k < 6; ++k) F[k] -= Cx[k];
    const Mat3 R = rotationFromVector(Vec3(xi[3], xi[4], xi[5]));
    for (size_t i = 0; i < lines_.size(); ++i) RigidBody6DOF::addPointLoad(rotate(R, lines_[i].a), lf[i], F);
    return F;
}

Vec6 CoupledSystem::staticResidual() {
    std::vector<Vec3> lf(lines_.size());
    for (size_t i = 0; i < lines_.size(); ++i) lf[i] = lines_[i].line->staticForceOnBody(body_.pointPosition(lines_[i].a));
    return externalLoad(body_.xi, lf);
}

EquilibriumResult CoupledSystem::solveEquilibrium(double tol, int maxIter, const std::array<bool, 6>& fixed) {
    EquilibriumResult res;
    body_.xiDot = Vec6{};
    auto norms = [&fixed](Vec6 G, double& f, double& m) {
        for (int k = 0; k < 6; ++k) if (fixed[k]) G[k] = 0.0;        // a held DOF carries a reaction, not a residual
        f = std::sqrt(G[0] * G[0] + G[1] * G[1] + G[2] * G[2]); m = std::sqrt(G[3] * G[3] + G[4] * G[4] + G[5] * G[5]);
    };
    Vec6 G = staticResidual();
    for (int it = 0; it < maxIter; ++it) {
        norms(G, res.forceResidual, res.momentResidual);
        res.iterations = it;
        if (res.forceResidual < tol && res.momentResidual < tol) { res.converged = true; break; }
        Mat6 J;                                           // J = dG/dxi by forward differences
        const Vec6 x0 = body_.xi;
        for (int j = 0; j < 6; ++j) {
            if (fixed[j]) { J(j, j) = 1.0; continue; }
            const double h = 1e-4;
            body_.xi = x0; body_.xi[j] += h;
            const Vec6 Gj = staticResidual();
            for (int i = 0; i < 6; ++i) J(i, j) = fixed[i] ? 0.0 : (Gj[i] - G[i]) / h;
        }
        body_.xi = x0;
        Vec6 rhs; for (int i = 0; i < 6; ++i) rhs[i] = fixed[i] ? 0.0 : -G[i];
        double jmax = 0; for (double v : J.a) jmax = std::max(jmax, std::fabs(v));
        for (int i = 0; i < 6; ++i) J(i, i) += 1e-10 * jmax + 1e-12;     // Tikhonov guard for a DOF without any restoring force (a slack line gives a zero Jacobian: start taut)
        const Vec6 d = solve6(J, rhs);
        double lam = 1.0, g0 = res.forceResidual + res.momentResidual;
        for (int k = 0; k < 30; ++k, lam *= 0.5) {
            for (int i = 0; i < 6; ++i) body_.xi[i] = x0[i] + lam * d[i];
            const Vec6 Gn = staticResidual();
            double f, m; norms(Gn, f, m);
            if (f + m < g0 || k == 29) { G = Gn; break; }
        }
    }
    norms(G, res.forceResidual, res.momentResidual);
    res.converged = res.forceResidual < tol * 10 && res.momentResidual < tol * 10;
    return res;
}

void CoupledSystem::begin(double t0) {
    t_ = t0;
    std::vector<Vec3> lf(lines_.size());
    for (size_t i = 0; i < lines_.size(); ++i) {
        Slot& s = lines_[i];
        s.force = s.line->forceOnBody(body_.pointPosition(s.a), body_.pointVelocity(s.a), t0);
        lf[i] = s.force;
    }
    Fnow_ = externalLoad(body_.xi, lf);
    if (excitation_) { const Vec6 e = excitation_(t0); for (int k = 0; k < 6; ++k) Fnow_[k] += e[k]; }
    Vec6 F = Fnow_;
    const Mat6 Bm = body_.params().B;
    const Vec6 Bv = Bm * body_.xiDot;
    for (int k = 0; k < 6; ++k) F[k] -= Bv[k] + body_.params().Dq[k] * std::fabs(body_.xiDot[k]) * body_.xiDot[k];
    acc_ = solve6(body_.totalMass(), F);
    report_ = CouplingReport{};
}

void CoupledSystem::step(double dt) {
    const BodyParams& p = body_.params();
    Vec6 vh;
    for (int k = 0; k < 6; ++k) vh[k] = body_.xiDot[k] + 0.5 * dt * acc_[k];
    for (int k = 0; k < 6; ++k) body_.xi[k] += dt * vh[k];
    const double tNew = t_ + dt;

    // Lines advanced to tNew with the fairlead moving linearly from its old to its new position.
    RigidBody6DOF probe = body_;                           // velocity estimate at the new pose uses the half-step velocity
    probe.xiDot = vh;
    std::vector<Vec3> lf(lines_.size());
    double dtLineMin = 1e300;
    for (size_t i = 0; i < lines_.size(); ++i) {
        Slot& s = lines_[i];
        s.force = s.line->forceOnBody(body_.pointPosition(s.a), probe.pointVelocity(s.a), tNew);
        lf[i] = s.force;
        dtLineMin = std::min(dtLineMin, s.line->internalTimeStep());
    }
    Fnow_ = externalLoad(body_.xi, lf);
    if (excitation_) { const Vec6 e = excitation_(tNew); for (int k = 0; k < 6; ++k) Fnow_[k] += e[k]; }

    // v_{n+1}: (M + dt/2 B) v = M vh + dt/2 (F - Dq|vh|vh)
    Vec6 rhs = body_.totalMass() * vh;
    for (int k = 0; k < 6; ++k) rhs[k] += 0.5 * dt * (Fnow_[k] - p.Dq[k] * std::fabs(vh[k]) * vh[k]);
    Mat6 lhs = body_.totalMass() + p.B * (0.5 * dt);
    const Vec6 vNew = solve6(lhs, rhs);
    body_.xiDot = vNew;
    Vec6 F = Fnow_;
    const Vec6 Bv = p.B * vNew;
    for (int k = 0; k < 6; ++k) F[k] -= Bv[k] + p.Dq[k] * std::fabs(vNew[k]) * vNew[k];
    acc_ = solve6(body_.totalMass(), F);
    t_ = tNew;
    ++report_.bodySteps;
    report_.dtBody = std::max(report_.dtBody, dt);                 // nominal (largest) body step; the last step may be shortened
    report_.dtLineMin = dtLineMin;
    report_.subStepRatio = static_cast<int>(std::ceil(report_.dtBody / dtLineMin - 1e-12));
}

void CoupledSystem::advanceTo(double tEnd, double dt) {
    while (t_ < tEnd - 1e-9 * dt) step(std::min(dt, tEnd - t_));
}

}  // namespace mooring
