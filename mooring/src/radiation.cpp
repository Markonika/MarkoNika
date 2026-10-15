// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#include "mooring/radiation.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace mooring {
namespace {
const double kPi = 3.14159265358979323846;

// Sine integral Si(x) = integral_0^x sin(u)/u du, x >= 0: power series up to 12, asymptotic expansion beyond.
double sineIntegral(double x) {
    if (x < 12.0) {
        double sum = 0.0, term = x;
        for (int k = 0; k < 80; ++k) {
            sum += term / (2 * k + 1);
            term *= -x * x / ((2.0 * k + 2.0) * (2.0 * k + 3.0));
            if (std::fabs(term) < 1e-18 * std::max(1.0, std::fabs(sum))) break;
        }
        return sum;
    }
    double f = 0.0, g = 0.0, tf = 1.0 / x, tg = 1.0 / (x * x);
    double ff = tf, gg = tg;                         // f = 1/x (1 - 2!/x^2 + 4!/x^4 ...), g = 1/x^2 (1 - 3!/x^2 + 5!/x^4 ...)
    f = ff; g = gg;
    for (int k = 1; k < 12; ++k) {
        ff *= -(2.0 * k) * (2.0 * k - 1.0) / (x * x);
        gg *= -(2.0 * k + 1.0) * (2.0 * k) / (x * x);
        if (std::fabs(ff) > std::fabs(f) && k > 3) break;
        f += ff; g += gg;
    }
    return kPi / 2.0 - f * std::cos(x) - g * std::sin(x);
}
}  // namespace

RadiationModel::RadiationModel(const HydroTable& t, double TK, double hK, bool tail) : t_(t), TK_(TK), hK_(hK) {
    const size_t n = t.omega.size();
    if (n < 2 || t.A.size() != n || t.B.size() != n) throw std::invalid_argument("hydro table needs at least two frequencies with A and B at each");
    for (size_t i = 0; i < n; ++i) if (t.omega[i] <= 0.0 || (i && t.omega[i] <= t.omega[i - 1])) throw std::invalid_argument("hydro table frequencies must be positive and ascending");
    if (TK <= 0.0 || hK <= 0.0 || hK > TK) throw std::invalid_argument("invalid retardation length / step");
    for (int i = 0; i < 36; ++i) idx_[i] = -1;
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
        double mx = 0.0;
        for (size_t k = 0; k < n; ++k) mx = std::max(mx, std::fabs(t.B[k](i, j)));
        if (mx > 0.0) { idx_[6 * i + j] = int(pairs_.size()); pairs_.push_back({i, j}); }
    }
    const int nK = int(std::ceil(TK / hK));
    kern_.assign(pairs_.size(), std::vector<double>(nK + 1, 0.0));
    for (size_t p = 0; p < pairs_.size(); ++p) {
        const int i = pairs_[p].first, j = pairs_[p].second;
        for (int m = 0; m <= nK; ++m) {
            const double tt = m * hK;
            double sum = 0.0;
            double wa = 0.0, Ba = 0.0;                          // segment [0, omega_1]: B(0) = 0
            for (size_t k = 0; k <= n; ++k) {
                double wb, Bb;
                if (k < n) { wb = t.omega[k]; Bb = t.B[k](i, j); } else break;
                const double d = wb - wa, slope = (Bb - Ba) / d;
                if (tt * d < 1e-6) sum += 0.5 * (Ba + Bb) * d;
                else sum += (Bb * std::sin(wb * tt) - Ba * std::sin(wa * tt)) / tt + slope / (tt * tt) * (std::cos(wb * tt) - std::cos(wa * tt));
                wa = wb; Ba = Bb;
            }
            if (tail) {                                          // B = B_N (w_N / w)^2 above the table
                const double wN = t.omega[n - 1], BN = t.B[n - 1](i, j);
                if (tt * wN < 1e-9) sum += BN * wN;
                else sum += BN * wN * wN * (std::cos(wN * tt) / wN - tt * (kPi / 2.0 - sineIntegral(wN * tt)));
            }
            kern_[p][m] = 2.0 / kPi * sum;
        }
    }
    // A_inf
    if (t.hasAinf) Ainf_ = t.Ainf;
    else {
        const double w = t.omega[n - 1];
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
            double integral = 0.0;
            const int p = idx_[6 * i + j];
            if (p >= 0) {
                for (int m = 0; m < nK; ++m) {
                    const double t0 = m * hK, t1 = (m + 1) * hK, K0v = kern_[p][m], K1v = kern_[p][m + 1], s = (K1v - K0v) / hK;
                    integral += (K0v * std::cos(w * t0) - K1v * std::cos(w * t1)) / w + s / (w * w) * (std::sin(w * t1) - std::sin(w * t0));
                }
            }
            Ainf_(i, j) = t.A[n - 1](i, j) + integral / w;
        }
    }
}

double RadiationModel::K(int i, int j, double s) const {
    const int p = idx_[6 * i + j];
    if (p < 0 || s >= TK_) return 0.0;
    const double x = std::max(0.0, s) / hK_;
    const size_t m = std::min<size_t>(size_t(x), kern_[p].size() - 2);
    const double f = x - double(m);
    return kern_[p][m] * (1.0 - f) + kern_[p][m + 1] * f;
}

static Mat6 interp(const HydroTable& t, const std::vector<Mat6>& V, double w) {
    const size_t n = t.omega.size();
    if (w <= t.omega.front()) return V.front();
    if (w >= t.omega.back()) return V.back();
    const size_t k = std::upper_bound(t.omega.begin(), t.omega.end(), w) - t.omega.begin();
    const double f = (w - t.omega[k - 1]) / (t.omega[k] - t.omega[k - 1]);
    return V[k - 1] * (1.0 - f) + V[k] * f;
    (void)n;
}
Mat6 RadiationModel::A(double w) const { return interp(t_, t_.A, w); }
Mat6 RadiationModel::B(double w) const { return interp(t_, t_.B, w); }

double RadiationModel::AFromKernel(int i, int j, double w) const {
    const int p = idx_[6 * i + j];
    double integral = 0.0;
    if (p >= 0) {
        const int nK = int(kern_[p].size()) - 1;
        for (int m = 0; m < nK; ++m) {
            const double t0 = m * hK_, t1 = (m + 1) * hK_, K0v = kern_[p][m], K1v = kern_[p][m + 1], s = (K1v - K0v) / hK_;
            integral += (K0v * std::cos(w * t0) - K1v * std::cos(w * t1)) / w + s / (w * w) * (std::sin(w * t1) - std::sin(w * t0));
        }
    }
    return Ainf_(i, j) - integral / w;
}

double RadiationModel::BFromKernel(int i, int j, double w) const {
    const int p = idx_[6 * i + j];
    double integral = 0.0;
    if (p >= 0) {
        const int nK = int(kern_[p].size()) - 1;
        for (int m = 0; m < nK; ++m) {
            const double t0 = m * hK_, t1 = (m + 1) * hK_, K0v = kern_[p][m], K1v = kern_[p][m + 1], s = (K1v - K0v) / hK_;
            integral += (K1v * std::sin(w * t1) - K0v * std::sin(w * t0)) / w + s / (w * w) * (std::cos(w * t1) - std::cos(w * t0));
        }
    }
    return integral;
}

void RadiationMemory::reset(double t0, const Vec6& v0) { t_.clear(); v_.clear(); t_.push_back(t0); v_.push_back(v0); }

void RadiationMemory::push(double t, const Vec6& v) {
    t_.push_back(t); v_.push_back(v);
    while (t_.size() > 2 && t - t_[1] >= m_->length()) { t_.pop_front(); v_.pop_front(); }     // intervals entirely beyond the kernel length
}

Vec6 RadiationMemory::known(double tNew) const {
    Vec6 I{};
    const size_t n = t_.size();
    if (n == 0) return I;
    const auto& pr = m_->pairs();
    for (const auto& pq : pr) {
        const int i = pq.first, j = pq.second;
        double sum = 0.0;
        for (size_t m = 0; m + 1 < n; ++m) {
            const double dtm = t_[m + 1] - t_[m];
            sum += 0.5 * dtm * (m_->K(i, j, tNew - t_[m]) * v_[m][j] + m_->K(i, j, tNew - t_[m + 1]) * v_[m + 1][j]);
        }
        const double dtl = tNew - t_[n - 1];
        sum += 0.5 * dtl * m_->K(i, j, dtl) * v_[n - 1][j];
        I[i] += sum;
    }
    return I;
}

Mat6 RadiationMemory::zeroLag(double dt) const {
    Mat6 C;
    for (const auto& pq : m_->pairs()) C(pq.first, pq.second) = 0.5 * dt * m_->K0(pq.first, pq.second);
    return C;
}

BodyModes bodyModesRadiation(const Mat6& Mrb, const RadiationModel& rad, const Mat6& K) {
    BodyModes out = bodyModes(Mrb + rad.Ainf(), K);
    for (int k = 0; k < 6; ++k) {
        double w = out.omega[k];
        if (w <= 0.0) continue;
        BodyModes cur;
        for (int it = 0; it < 100; ++it) {
            cur = bodyModes(Mrb + rad.A(w), K);
            const double wn = cur.omega[k];
            if (std::fabs(wn - w) < 1e-12 * std::max(1.0, w)) { w = wn; break; }
            w = wn;
        }
        out.omega[k] = w; out.omega2[k] = w * w; out.period[k] = 2.0 * kPi / w; out.shape[k] = cur.shape[k];
    }
    return out;
}

}  // namespace mooring
