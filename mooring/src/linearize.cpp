// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#include "mooring/linearize.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include "mooring/lumped_mass_cable.hpp"
#include "mooring/rigid_body.hpp"

namespace mooring {

SymEigen symmetricEigen(std::vector<double> a, int n) {
    std::vector<double> v(size_t(n) * n, 0.0);
    for (int i = 0; i < n; ++i) v[size_t(i) * n + i] = 1.0;       // v holds the eigenvectors as columns during the sweeps
    auto A = [&](int i, int j) -> double& { return a[size_t(i) * n + j]; };
    auto V = [&](int i, int j) -> double& { return v[size_t(i) * n + j]; };
    for (int sweep = 0; sweep < 60; ++sweep) {
        double off = 0.0, diag = 0.0;
        for (int i = 0; i < n; ++i) { diag += A(i, i) * A(i, i); for (int j = i + 1; j < n; ++j) off += A(i, j) * A(i, j); }
        if (off <= 1e-26 * (diag + off) || off == 0.0) break;
        for (int p = 0; p < n - 1; ++p) {
            for (int q = p + 1; q < n; ++q) {
                const double apq = A(p, q);
                if (std::fabs(apq) < 1e-300) continue;
                const double theta = (A(q, q) - A(p, p)) / (2.0 * apq);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < n; ++k) {                       // rotate columns p, q
                    const double akp = A(k, p), akq = A(k, q);
                    A(k, p) = c * akp - s * akq; A(k, q) = s * akp + c * akq;
                }
                for (int k = 0; k < n; ++k) {                       // rotate rows p, q
                    const double apk = A(p, k), aqk = A(q, k);
                    A(p, k) = c * apk - s * aqk; A(q, k) = s * apk + c * aqk;
                }
                for (int k = 0; k < n; ++k) {
                    const double vkp = V(k, p), vkq = V(k, q);
                    V(k, p) = c * vkp - s * vkq; V(k, q) = s * vkp + c * vkq;
                }
            }
        }
    }
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int x, int y) { return A(x, x) < A(y, y); });
    SymEigen r; r.n = n; r.w.resize(n); r.v.resize(size_t(n) * n);
    for (int k = 0; k < n; ++k) {
        r.w[k] = A(idx[k], idx[k]);
        for (int i = 0; i < n; ++i) r.v[size_t(k) * n + i] = V(i, idx[k]);
    }
    return r;
}

LineModes LumpedMassCable::naturalModes(int nModes) {
    const int N = p_.N;
    const std::vector<int> comps = p_.planar ? std::vector<int>{0, 2} : std::vector<int>{0, 1, 2};
    const int nc = int(comps.size()), n = nc * (N - 1);
    std::vector<Vec3> rf = r_, vf(r_.size());
    buildJacobian(rf, vf, t_);
    --stats_.jacobianBuilds;                                       // an analysis call, not a solver rebuild
    jacValid_ = false;
    std::vector<double> K(size_t(n) * n, 0.0), C(size_t(n) * n, 0.0);
    std::vector<std::array<double, 9>> Mh(N + 1);                  // M^{-1/2} per node (3x3)
    std::vector<std::array<double, 9>> Mm(N + 1);                  // M per node
    for (int i = 1; i < N; ++i) {
        double mt, mn; Vec3 tg;
        nodeInertia(r_, t_, i, mt, mn, tg);
        const double t3[3] = {tg.x, tg.y, tg.z};
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
            const double tt = t3[a] * t3[b], id = a == b ? 1.0 : 0.0;
            Mm[i][3 * a + b] = mt * tt + mn * (id - tt);
            Mh[i][3 * a + b] = tt / std::sqrt(mt) + (id - tt) / std::sqrt(mn);
        }
        for (int j = std::max(1, i - bw_); j <= std::min(N - 1, i + bw_); ++j) {
            const double* br = dAdr_[i].a[j - i + bw_];
            const double* bv = dAdv_[i].a[j - i + bw_];
            for (int ra = 0; ra < nc; ++ra) for (int cb = 0; cb < nc; ++cb) {
                double kr = 0.0, kv = 0.0;
                for (int c3 = 0; c3 < 3; ++c3) { kr += Mm[i][3 * comps[ra] + c3] * br[3 * c3 + comps[cb]]; kv += Mm[i][3 * comps[ra] + c3] * bv[3 * c3 + comps[cb]]; }
                K[size_t(nc * (i - 1) + ra) * n + nc * (j - 1) + cb] = -kr;
                C[size_t(nc * (i - 1) + ra) * n + nc * (j - 1) + cb] = -kv;
            }
        }
    }
    LineModes out; out.dof = n;
    double asym = 0.0, nrm = 0.0;
    for (int a = 0; a < n; ++a) for (int b = a + 1; b < n; ++b) {
        const double d = K[size_t(a) * n + b] - K[size_t(b) * n + a];
        asym += d * d;
        const double m = 0.5 * (K[size_t(a) * n + b] + K[size_t(b) * n + a]);
        K[size_t(a) * n + b] = K[size_t(b) * n + a] = m;
        const double mc = 0.5 * (C[size_t(a) * n + b] + C[size_t(b) * n + a]);
        C[size_t(a) * n + b] = C[size_t(b) * n + a] = mc;
    }
    for (double x : K) nrm += x * x;
    out.asymmetry = std::sqrt(2.0 * asym / std::max(nrm, 1e-300));
    // S = M^{-1/2} K M^{-1/2}: block rows / columns scaled by the 3x3 blocks restricted to the active components
    auto scale = [&](std::vector<double>& A) {
        std::vector<double> T(size_t(n) * n, 0.0);
        for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) {            // T = Mh A
            const int ia = a / nc + 1, ra = a % nc;
            double s = 0.0;
            for (int rb = 0; rb < nc; ++rb) s += Mh[ia][3 * comps[ra] + comps[rb]] * A[size_t(nc * (ia - 1) + rb) * n + b];
            T[size_t(a) * n + b] = s;
        }
        std::vector<double> R(size_t(n) * n, 0.0);
        for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) {            // R = T Mh
            const int ib = b / nc + 1, cb = b % nc;
            double s = 0.0;
            for (int rb = 0; rb < nc; ++rb) s += T[size_t(a) * n + nc * (ib - 1) + rb] * Mh[ib][3 * comps[rb] + comps[cb]];
            R[size_t(a) * n + b] = s;
        }
        A = R;
    };
    std::vector<double> S = K, Sc = C;
    scale(S); scale(Sc);
    for (int a = 0; a < n; ++a) for (int b = a + 1; b < n; ++b) { S[size_t(a) * n + b] = S[size_t(b) * n + a] = 0.5 * (S[size_t(a) * n + b] + S[size_t(b) * n + a]); }
    const SymEigen eg = symmetricEigen(S, n);
    const int m = std::min(nModes, n);
    for (int k = 0; k < m; ++k) {
        out.omega2.push_back(eg.w[k]);
        const double om = eg.w[k] > 0.0 ? std::sqrt(eg.w[k]) : 0.0;
        out.omega.push_back(om);
        double yCy = 0.0;                                           // y^T (M^{-1/2} C M^{-1/2}) y
        for (int a = 0; a < n; ++a) { double s = 0.0; for (int b = 0; b < n; ++b) s += Sc[size_t(a) * n + b] * eg.v[size_t(k) * n + b]; yCy += eg.v[size_t(k) * n + a] * s; }
        out.zeta.push_back(om > 0.0 ? yCy / (2.0 * om) : 0.0);
        std::vector<Vec3> shape(N + 1);
        for (int i = 1; i < N; ++i) {
            double y[3] = {0, 0, 0};
            for (int r = 0; r < nc; ++r) y[comps[r]] = eg.v[size_t(k) * n + nc * (i - 1) + r];
            Vec3 phi;
            double* pv[3] = {&phi.x, &phi.y, &phi.z};
            for (int a = 0; a < 3; ++a) { double s = 0.0; for (int b = 0; b < 3; ++b) s += Mh[i][3 * a + b] * y[b]; *pv[a] = s; }
            shape[i] = phi;
        }
        out.shape.push_back(shape);
    }
    return out;
}

namespace {
using cd = std::complex<double>;
inline double& comp(Vec3& v, int d) { return d == 0 ? v.x : (d == 1 ? v.y : v.z); }
// Solve A X = B for X (n x m, row-major) by Gaussian elimination with partial pivoting; A is overwritten.
void complexSolve(std::vector<cd>& A, std::vector<cd>& B, int n, int m) {
    for (int c = 0; c < n; ++c) {
        int piv = c;
        for (int r = c + 1; r < n; ++r) if (std::abs(A[size_t(r) * n + c]) > std::abs(A[size_t(piv) * n + c])) piv = r;
        if (std::abs(A[size_t(piv) * n + c]) < 1e-300) throw std::runtime_error("fairleadDynamicStiffness: singular interior matrix");
        if (piv != c) { for (int j = 0; j < n; ++j) std::swap(A[size_t(c) * n + j], A[size_t(piv) * n + j]); for (int j = 0; j < m; ++j) std::swap(B[size_t(c) * m + j], B[size_t(piv) * m + j]); }
        for (int r = c + 1; r < n; ++r) {
            const cd f = A[size_t(r) * n + c] / A[size_t(c) * n + c];
            if (f == cd(0)) continue;
            for (int j = c; j < n; ++j) A[size_t(r) * n + j] -= f * A[size_t(c) * n + j];
            for (int j = 0; j < m; ++j) B[size_t(r) * m + j] -= f * B[size_t(c) * m + j];
        }
    }
    for (int c = n - 1; c >= 0; --c)
        for (int j = 0; j < m; ++j) {
            cd s = B[size_t(c) * m + j];
            for (int k = c + 1; k < n; ++k) s -= A[size_t(c) * n + k] * B[size_t(k) * m + j];
            B[size_t(c) * m + j] = s / A[size_t(c) * n + c];
        }
}
}  // namespace

void LumpedMassCable::fairleadDynamicStiffness(double omega, std::complex<double> D[9]) {
    using cd = std::complex<double>;
    const int N = p_.N;
    const std::vector<int> comps = p_.planar ? std::vector<int>{0, 2} : std::vector<int>{0, 1, 2};
    const int nc = int(comps.size()), nI = nc * (N - 1), nU = nI + nc;         // interior dofs, then the top node
    std::vector<Vec3> r = r_, v = v_.size() == r_.size() ? v_ : std::vector<Vec3>(r_.size());
    double l0min = 1e300;
    for (const SegProps& g : seg_) l0min = std::min(l0min, g.l0);
    const double hr = 1e-6 * l0min, hv = 1e-6;
    auto unode = [&](int u) { return u < nI ? u / nc + 1 : N; };
    auto ucomp = [&](int u) { return comps[u % nc]; };
    std::vector<double> Kr(size_t(nU) * nU, 0.0), Kv(size_t(nU) * nU, 0.0);   // dF_row/d(r or v)_col
    std::vector<Vec3> fp, fm;
    auto fvec = [&](const Vec3& x, int comp_) { return comp_ == 0 ? x.x : (comp_ == 1 ? x.y : x.z); };
    for (int mode = 0; mode < 2; ++mode) {
        const double h = mode == 0 ? hr : hv;
        for (int c = 0; c < nU; ++c) {
            std::vector<Vec3> rp = r, vp = v, rm = r, vm = v;
            Vec3& xp = mode == 0 ? rp[unode(c)] : vp[unode(c)];
            Vec3& xm = mode == 0 ? rm[unode(c)] : vm[unode(c)];
            comp(xp, ucomp(c)) += h; comp(xm, ucomp(c)) -= h;
            computeForces(rp, vp, t_, fp);
            computeForces(rm, vm, t_, fm);
            for (int rr = 0; rr < nU; ++rr)
                (mode == 0 ? Kr : Kv)[size_t(rr) * nU + c] = (fvec(fp[unode(rr)], ucomp(rr)) - fvec(fm[unode(rr)], ucomp(rr))) / (2.0 * h);
        }
    }
    // Interior inertia (3x3 blocks restricted to the active components)
    std::vector<cd> Gii(size_t(nI) * nI), Git(size_t(nI) * nc), Gti(size_t(nc) * nI), Gtt(size_t(nc) * nc);
    const cd iw(0.0, omega);
    auto G = [&](int a, int b) { return cd(Kr[size_t(a) * nU + b], 0.0) + iw * Kv[size_t(a) * nU + b]; };
    for (int a = 0; a < nI; ++a) for (int b = 0; b < nI; ++b) Gii[size_t(a) * nI + b] = G(a, b);
    for (int a = 0; a < nI; ++a) for (int b = 0; b < nc; ++b) Git[size_t(a) * nc + b] = G(a, nI + b);
    for (int a = 0; a < nc; ++a) for (int b = 0; b < nI; ++b) Gti[size_t(a) * nI + b] = G(nI + a, b);
    for (int a = 0; a < nc; ++a) for (int b = 0; b < nc; ++b) Gtt[size_t(a) * nc + b] = G(nI + a, nI + b);
    for (int i = 1; i < N; ++i) {
        double mt, mn; Vec3 tg;
        nodeInertia(r_, t_, i, mt, mn, tg);
        const double t3[3] = {tg.x, tg.y, tg.z};
        for (int ra = 0; ra < nc; ++ra) for (int cb = 0; cb < nc; ++cb) {
            const int a = comps[ra], b = comps[cb];
            const double tt = t3[a] * t3[b], id = a == b ? 1.0 : 0.0;
            Gii[size_t(nc * (i - 1) + ra) * nI + nc * (i - 1) + cb] += omega * omega * (mt * tt + mn * (id - tt));
        }
    }
    // Schur complement: Gtt - Gti (Gii + w^2 M)^-1 Git
    std::vector<cd> X = Git;
    complexSolve(Gii, X, nI, nc);
    for (int a = 0; a < 9; ++a) D[a] = cd(0.0, 0.0);
    for (int ra = 0; ra < nc; ++ra) for (int cb = 0; cb < nc; ++cb) {
        cd s = Gtt[size_t(ra) * nc + cb];
        for (int k = 0; k < nI; ++k) s -= Gti[size_t(ra) * nI + k] * X[size_t(k) * nc + cb];
        D[3 * comps[ra] + comps[cb]] = -s;
    }
}

std::array<std::complex<double>, 36> CoupledSystem::dynamicLineStiffness(double omega) {
    using cd = std::complex<double>;
    std::array<cd, 36> K{};
    const Mat3 R = body_.rotation();
    for (int k = 0; k < numLines(); ++k) {
        auto* lm = dynamic_cast<LumpedMassCable*>(&line(k));
        if (!lm) throw std::invalid_argument("dynamicLineStiffness: lumped-mass lines only");
        cd D[9];
        lm->fairleadDynamicStiffness(omega, D);
        const Vec3 Ra = rotate(R, lines_[k].a);
        const double ra[3] = {Ra.x, Ra.y, Ra.z};
        // J (3 x 6): u = dx - [Ra]x dtheta;  [Ra]x = [[0,-z,y],[z,0,-x],[-y,x,0]]
        double J[3][6] = {};
        for (int i = 0; i < 3; ++i) J[i][i] = 1.0;
        const double cr[3][3] = {{0, -ra[2], ra[1]}, {ra[2], 0, -ra[0]}, {-ra[1], ra[0], 0}};
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) J[i][3 + j] = -cr[i][j];
        for (int a = 0; a < 6; ++a) for (int b = 0; b < 6; ++b) {
            cd s(0.0, 0.0);
            for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) s += J[i][a] * D[3 * i + j] * J[j][b];
            K[6 * a + b] += s;
        }
        // geometric stiffness: dQ_rot = (Ra f0^T - (f0.Ra) I) dtheta  ->  K_geo = -(Ra f0^T - (f0.Ra) I)
        const Vec3 f0 = lm->endForce(true);
        const double f3[3] = {f0.x, f0.y, f0.z};
        const double fr = f0.x * Ra.x + f0.y * Ra.y + f0.z * Ra.z;
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) K[6 * (3 + i) + (3 + j)] += cd(-(ra[i] * f3[j] - (i == j ? fr : 0.0)), 0.0);
    }
    return K;
}

Mat6 CoupledSystem::stiffness(double h) {
    const Vec6 x0 = body_.xi;
    Mat6 K;
    for (int j = 0; j < 6; ++j) {
        body_.xi = x0; body_.xi[j] += h;
        const Vec6 Gp = staticResidual();
        body_.xi = x0; body_.xi[j] -= h;
        const Vec6 Gm = staticResidual();
        for (int i = 0; i < 6; ++i) K(i, j) = -(Gp[i] - Gm[i]) / (2.0 * h);
    }
    body_.xi = x0;
    staticResidual();                                             // lines back at the original pose
    return K;
}

Mat6 CoupledSystem::lineStiffness(double h) {
    Mat6 K = stiffness(h);
    for (int i = 0; i < 36; ++i) K.a[i] -= body_.params().C.a[i];
    return K;
}

BodyModes bodyModes(const Mat6& M, const Mat6& Kin) {
    // Cholesky M = L L^T, S = L^-1 K L^-T (symmetric), Jacobi, x = L^-T y.
    double L[6][6] = {};
    for (int i = 0; i < 6; ++i) for (int j = 0; j <= i; ++j) {
        double s = 0.5 * (M(i, j) + M(j, i));
        for (int k = 0; k < j; ++k) s -= L[i][k] * L[j][k];
        if (i == j) { if (s <= 0.0) throw std::runtime_error("bodyModes: mass matrix not positive definite"); L[i][i] = std::sqrt(s); }
        else L[i][j] = s / L[j][j];
    }
    BodyModes out;
    double asym = 0.0, nrm = 0.0;
    double K[6][6];
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) { K[i][j] = 0.5 * (Kin(i, j) + Kin(j, i)); asym += std::pow(Kin(i, j) - Kin(j, i), 2); nrm += Kin(i, j) * Kin(i, j); }
    out.asymmetry = std::sqrt(asym / std::max(nrm, 1e-300));
    double Y[6][6];                                              // Y = L^-1 K
    for (int c = 0; c < 6; ++c) for (int i = 0; i < 6; ++i) { double s = K[i][c]; for (int k = 0; k < i; ++k) s -= L[i][k] * Y[k][c]; Y[i][c] = s / L[i][i]; }
    double S[6][6];                                              // S = Y L^-T  (solve S L^T = Y, row by row)
    for (int r = 0; r < 6; ++r) for (int j = 0; j < 6; ++j) { double s = Y[r][j]; for (int k = 0; k < j; ++k) s -= S[r][k] * L[j][k]; S[r][j] = s / L[j][j]; }
    std::vector<double> a(36);
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) a[6 * i + j] = 0.5 * (S[i][j] + S[j][i]);
    const SymEigen eg = symmetricEigen(a, 6);
    for (int k = 0; k < 6; ++k) {
        out.omega2[k] = eg.w[k];
        out.omega[k] = eg.w[k] > 0.0 ? std::sqrt(eg.w[k]) : 0.0;
        out.period[k] = out.omega[k] > 0.0 ? 6.283185307179586 / out.omega[k] : 0.0;
        double x[6];                                             // x = L^-T y (back substitution)
        for (int i = 5; i >= 0; --i) { double s = eg.v[size_t(k) * 6 + i]; for (int j = i + 1; j < 6; ++j) s -= L[j][i] * x[j]; x[i] = s / L[i][i]; }
        for (int i = 0; i < 6; ++i) out.shape[k][i] = x[i];
    }
    return out;
}

}  // namespace mooring
