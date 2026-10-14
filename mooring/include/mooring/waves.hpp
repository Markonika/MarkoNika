// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Linear (Airy) wave kinematics with optional Wheeler stretching (milestone 9).
#pragma once
#include <cstdint>
#include <vector>
#include "mooring/current.hpp"
#include "mooring/lumped_mass_cable.hpp"
#include "mooring/vec3.hpp"

namespace mooring {

struct WaveComponent {
    double amplitude{0};     // [m]
    double omega{0};         // [rad/s]
    double k{0};             // wave number [rad/m], from the finite-depth dispersion relation
    double phase{0};         // [rad]
    double direction{0};     // propagation direction, angle from +x about +z [rad]
    double cdir{1}, sdir{0}, q{0};   // cos / sin of the direction and exp(-2 k h): set by WaveField when the component is added (speed)
};

enum class Stretching { None, Wheeler };

// JONSWAP spectrum (gamma = 1 gives the Pierson-Moskowitz spectrum with the same Hs, Tp):
//   S(w) = norm w^-5 exp(-1.25 (wp/w)^4) gamma^exp(-(w - wp)^2 / (2 sigma^2 wp^2)),  sigma = 0.07 (w <= wp), 0.09 (w > wp),
// with 'norm' fixed numerically so that the zeroth moment is exactly m0 = Hs^2/16. Units m^2 s/rad.
class JonswapSpectrum {
public:
    JonswapSpectrum(double hs, double tp, double gamma = 3.3);
    double density(double omega) const;
    double hs() const { return hs_; }
    double peakOmega() const { return wp_; }
private:
    double shape(double omega) const;
    double hs_, wp_, gamma_, norm_{1.0};
};

// Discretisation of a spectrum into Airy components: n frequency bins of equal width over [wMinFactor, wMaxFactor] x wp, one component per bin
// at a uniformly random frequency inside the bin (deterministic for a given seed; avoids the artificial repetition of equally spaced lines), amplitude
// sqrt(2 S(w_i) dw), random phase. The amplitudes are rescaled so that sum(a^2/2) equals the target m0 exactly.
struct SpectrumSpec {
    double hs{1.0}, tp{8.0}, gamma{3.3};
    int n{50};
    double wMinFactor{0.5}, wMaxFactor{4.0};
    std::uint64_t seed{1};
    double directionDeg{0.0};
};

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
    void addSpectrum(const SpectrumSpec& spec);        // see SpectrumSpec
    double significantHeight() const;                   // 4 sqrt(sum a^2 / 2) of the components

    static double wavenumber(double omega, double depth, double g = 9.81);   // omega^2 = g k tanh(k h)
    double elevation(double x, double y, double t) const;
    // Water velocity and acceleration at 'pos' (absolute z), zero above the instantaneous surface.
    void kinematics(const Vec3& pos, double t, Vec3& u, Vec3& a) const;
    WaterField asWaterField() const;
    std::function<double(const Vec3&, double)> asElevation() const;   // surface elevation above surfaceZ for Environment::elevation (copied)       // for Environment::water (the field is copied)
    const std::vector<WaveComponent>& components() const { return c_; }
    double depth() const { return h_; }
    double surfaceZ() const { return z0_; }
    double rampFactor(double t) const;

private:
    double h_, z0_;
    Stretching s_;
    double g_;
    double ramp_{0.0};
    void finalize(WaveComponent& w) const;
    std::vector<WaveComponent> c_;
};

}  // namespace mooring
