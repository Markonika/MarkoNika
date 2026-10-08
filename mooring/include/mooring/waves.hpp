// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Linear (Airy) wave kinematics with optional Wheeler stretching (milestone 9).
#pragma once
#include <vector>
#include "mooring/lumped_mass_cable.hpp"
#include "mooring/vec3.hpp"

namespace mooring {

struct WaveComponent {
    double amplitude{0};     // [m]
    double omega{0};         // [rad/s]
    double k{0};             // wave number [rad/m], from the finite-depth dispersion relation
    double phase{0};         // [rad]
    double direction{0};     // propagation direction, angle from +x about +z [rad]
};

enum class Stretching { None, Wheeler };

// Sum of Airy components over constant depth h. Elevation of component j at (x, y, t) is
//   eta_j = A_j sin(omega_j t - k_j (x cos th + y sin th) + phi_j),
// so with phi = 0 the elevation at the origin is A sin(omega t) (the convention of Eq. 3.60: f = w A sin(omega t + delta)).
// Velocity potential: u_h = A omega cosh(k(z'+h))/sinh(kh) sin(.), w = A omega sinh(k(z'+h))/sinh(kh) cos(.),
// z' = z - surfaceZ measured from the still-water level, with Wheeler stretching z' -> h (z' - eta)/(h + eta).
// Optional cosine ramp of all amplitudes over the first 'rampTime' seconds (the ramp derivative is neglected
// in the accelerations, which is accurate for rampTime >> period).
class WaveField {
public:
    WaveField(double depth, double surfaceZ, Stretching s = Stretching::Wheeler, double g = 9.81)
        : h_(depth), z0_(surfaceZ), s_(s), g_(g) {}
    void addRegular(double height, double period, double phase = 0.0, double directionDeg = 0.0);
    void addComponent(double amplitude, double period, double phase, double directionDeg);
    void setRampTime(double t) { ramp_ = t; }

    static double wavenumber(double omega, double depth, double g = 9.81);   // omega^2 = g k tanh(k h)
    double elevation(double x, double y, double t) const;
    // Water velocity and acceleration at 'pos' (absolute z), zero above the instantaneous surface.
    void kinematics(const Vec3& pos, double t, Vec3& u, Vec3& a) const;
    WaterField asWaterField() const;       // for Environment::water (the field is copied)
    const std::vector<WaveComponent>& components() const { return c_; }
    double depth() const { return h_; }
    double surfaceZ() const { return z0_; }
    double rampFactor(double t) const;

private:
    double h_, z0_;
    Stretching s_;
    double g_;
    double ramp_{0.0};
    std::vector<WaveComponent> c_;
};

}  // namespace mooring
