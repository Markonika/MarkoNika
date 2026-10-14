// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#include "mooring/waves.hpp"
#include <cmath>
#include <random>
#include <stdexcept>

namespace mooring {
namespace { const double kPi = 3.14159265358979323846; }

double WaveField::wavenumber(double omega, double depth, double g) {
    // Newton on f(k) = g k tanh(kh) - omega^2, starting from the deep-water value.
    double k = omega * omega / g;
    for (int i = 0; i < 100; ++i) {
        const double th = std::tanh(k * depth);
        const double f = g * k * th - omega * omega;
        const double df = g * th + g * k * depth * (1.0 - th * th);
        const double dk = f / df;
        k -= dk;
        if (std::fabs(dk) < 1e-14 * k) break;
    }
    return k;
}

void WaveField::finalize(WaveComponent& w) const {
    w.cdir = std::cos(w.direction); w.sdir = std::sin(w.direction); w.q = std::exp(-2.0 * w.k * h_);
}

void WaveField::addComponent(double amplitude, double period, double phase, double directionDeg) {
    if (period <= 0.0 || amplitude < 0.0) throw std::invalid_argument("invalid wave component");
    WaveComponent w;
    w.amplitude = amplitude; w.omega = 2.0 * kPi / period; w.k = wavenumber(w.omega, h_, g_);
    w.phase = phase; w.direction = directionDeg * kPi / 180.0;
    finalize(w);
    c_.push_back(w);
}

void WaveField::addRegular(double height, double period, double phase, double directionDeg) {
    addComponent(0.5 * height, period, phase, directionDeg);
}

double WaveField::rampFactor(double t) const {
    if (ramp_ <= 0.0 || t >= ramp_) return 1.0;
    if (t <= 0.0) return 0.0;
    return 0.5 * (1.0 - std::cos(kPi * t / ramp_));
}

double WaveField::elevation(double x, double y, double t) const {
    double eta = 0.0;
    for (const WaveComponent& w : c_) {
        const double s = x * w.cdir + y * w.sdir;
        eta += w.amplitude * std::sin(w.omega * t - w.k * s + w.phase);
    }
    return eta * rampFactor(t);
}

void WaveField::kinematics(const Vec3& pos, double t, Vec3& u, Vec3& a) const {
    u = Vec3(); a = Vec3();
    const size_t n = c_.size();
    // one sincos per component, shared by the elevation and the velocity / acceleration (the cost scales with the number of components)
    thread_local std::vector<double> sn, cs;
    sn.resize(n); cs.resize(n);
    double eta = 0.0;
    for (size_t j = 0; j < n; ++j) {
        const WaveComponent& w = c_[j];
        const double ph = w.omega * t - w.k * (pos.x * w.cdir + pos.y * w.sdir) + w.phase;
        sincos(ph, &sn[j], &cs[j]);
        eta += w.amplitude * sn[j];
    }
    const double r = rampFactor(t);
    eta *= r;
    const double zeta = pos.z - z0_;                       // from the still-water level
    if (zeta > eta) return;                                // above the instantaneous surface: no water
    double zp = zeta;
    if (s_ == Stretching::Wheeler) zp = h_ * (zeta - eta) / (h_ + eta);
    for (size_t j = 0; j < n; ++j) {
        const WaveComponent& w = c_[j];
        // cosh(k(z'+h))/sinh(kh) and sinh(k(z'+h))/sinh(kh) from one exponential: with E = e^{k z'}, q = e^{-2kh}
        const double q = w.q, E = std::exp(w.k * zp), Ei = q / E, den = 1.0 - q;
        const double ch = (E + Ei) / den, sh = (E - Ei) / den;
        const double cdir = w.cdir, sdir = w.sdir;
        const double A = w.amplitude * r, om = w.omega;
        const double uh = A * om * ch * sn[j], ah = A * om * om * ch * cs[j];
        const double uz = A * om * sh * cs[j], az = -A * om * om * sh * sn[j];
        u += Vec3(uh * cdir, uh * sdir, uz);
        a += Vec3(ah * cdir, ah * sdir, az);
    }
}

double JonswapSpectrum::shape(double w) const {
    if (w <= 0.0) return 0.0;
    const double sigma = w <= wp_ ? 0.07 : 0.09;
    const double r = (w - wp_) / (sigma * wp_);
    return std::pow(w, -5.0) * std::exp(-1.25 * std::pow(wp_ / w, 4.0)) * std::pow(gamma_, std::exp(-0.5 * r * r));
}

JonswapSpectrum::JonswapSpectrum(double hs, double tp, double gamma) : hs_(hs), wp_(2.0 * kPi / tp), gamma_(gamma) {
    if (hs <= 0.0 || tp <= 0.0 || gamma < 1.0) throw std::invalid_argument("invalid JONSWAP parameters (need Hs > 0, Tp > 0, gamma >= 1)");
    // m0 = integral of the shape over [0.2 wp, 30 wp] (Simpson, 40000 intervals); the neglected tails are below 1e-9 of m0
    const double a = 0.2 * wp_, b = 30.0 * wp_;
    const int n = 40000;
    const double h = (b - a) / n;
    double sum = shape(a) + shape(b);
    for (int i = 1; i < n; ++i) sum += shape(a + i * h) * (i % 2 ? 4.0 : 2.0);
    norm_ = (hs * hs / 16.0) / (sum * h / 3.0);
}

double JonswapSpectrum::density(double w) const { return norm_ * shape(w); }

void WaveField::addSpectrum(const SpectrumSpec& sp) {
    if (sp.n < 1 || sp.wMinFactor <= 0.0 || sp.wMaxFactor <= sp.wMinFactor) throw std::invalid_argument("invalid spectrum discretisation");
    const JonswapSpectrum S(sp.hs, sp.tp, sp.gamma);
    const double wp = S.peakOmega(), w0 = sp.wMinFactor * wp, dw = (sp.wMaxFactor - sp.wMinFactor) * wp / sp.n;
    std::mt19937_64 rng(sp.seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<WaveComponent> comps;
    double m0 = 0.0;
    for (int i = 0; i < sp.n; ++i) {
        WaveComponent w;
        w.omega = w0 + (i + u(rng)) * dw;
        w.amplitude = std::sqrt(2.0 * S.density(w.omega) * dw);
        w.phase = 2.0 * kPi * u(rng);
        w.k = wavenumber(w.omega, h_, g_);
        w.direction = sp.directionDeg * kPi / 180.0;
        finalize(w);
        m0 += 0.5 * w.amplitude * w.amplitude;
        comps.push_back(w);
    }
    const double scale = std::sqrt(sp.hs * sp.hs / 16.0 / m0);        // exact variance of the discretised spectrum
    for (WaveComponent& w : comps) { w.amplitude *= scale; c_.push_back(w); }
}

double WaveField::significantHeight() const {
    double m0 = 0.0;
    for (const WaveComponent& w : c_) m0 += 0.5 * w.amplitude * w.amplitude;
    return 4.0 * std::sqrt(m0);
}

WaterField WaveField::asWaterField() const {
    WaveField copy = *this;
    return [copy](const Vec3& pos, double t, Vec3& vw, Vec3& aw) { copy.kinematics(pos, t, vw, aw); };
}

std::function<double(const Vec3&, double)> WaveField::asElevation() const {
    WaveField copy = *this;
    return [copy](const Vec3& pos, double t) { return copy.elevation(pos.x, pos.y, t); };
}

}  // namespace mooring
