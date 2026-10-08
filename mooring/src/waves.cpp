// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#include "mooring/waves.hpp"
#include <cmath>
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

void WaveField::addComponent(double amplitude, double period, double phase, double directionDeg) {
    if (period <= 0.0 || amplitude < 0.0) throw std::invalid_argument("invalid wave component");
    WaveComponent w;
    w.amplitude = amplitude; w.omega = 2.0 * kPi / period; w.k = wavenumber(w.omega, h_, g_);
    w.phase = phase; w.direction = directionDeg * kPi / 180.0;
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
        const double s = x * std::cos(w.direction) + y * std::sin(w.direction);
        eta += w.amplitude * std::sin(w.omega * t - w.k * s + w.phase);
    }
    return eta * rampFactor(t);
}

void WaveField::kinematics(const Vec3& pos, double t, Vec3& u, Vec3& a) const {
    u = Vec3(); a = Vec3();
    const double eta = elevation(pos.x, pos.y, t);
    const double zeta = pos.z - z0_;                       // from the still-water level
    if (zeta > eta) return;                                // above the instantaneous surface: no water
    double zp = zeta;
    if (s_ == Stretching::Wheeler) zp = h_ * (zeta - eta) / (h_ + eta);
    const double r = rampFactor(t);
    for (const WaveComponent& w : c_) {
        const double cs = std::cos(w.direction), sn = std::sin(w.direction);
        const double ph = w.omega * t - w.k * (pos.x * cs + pos.y * sn) + w.phase;
        const double kh = w.k * h_, ch = std::cosh(w.k * (zp + h_)) / std::sinh(kh), sh = std::sinh(w.k * (zp + h_)) / std::sinh(kh);
        const double A = w.amplitude * r, om = w.omega;
        const double uh = A * om * ch * std::sin(ph), ah = A * om * om * ch * std::cos(ph);
        const double uz = A * om * sh * std::cos(ph), az = -A * om * om * sh * std::sin(ph);
        u += Vec3(uh * cs, uh * sn, uz);
        a += Vec3(ah * cs, ah * sn, az);
    }
}

WaterField WaveField::asWaterField() const {
    WaveField copy = *this;
    return [copy](const Vec3& pos, double t, Vec3& vw, Vec3& aw) { copy.kinematics(pos, t, vw, aw); };
}

}  // namespace mooring
