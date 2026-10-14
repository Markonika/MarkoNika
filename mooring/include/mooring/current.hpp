// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include "mooring/vec3.hpp"

namespace mooring {

// Steady current velocity as a function of depth. Depth d = surfaceZ - z is measured down from the still-water level; above the surface and
// between the given points the profile is the surface value / linear interpolation, below the deepest point the last value is kept
// (use an explicit zero entry at the seabed for a bottom boundary layer). The current is horizontal or 3-D as given; no wave-current interaction.
class CurrentProfile {
public:
    explicit CurrentProfile(double surfaceZ = 0.0) : z0_(surfaceZ) {}
    // Table point: depth below the still-water level (>= 0, strictly increasing in the order added) and the velocity there.
    void addPoint(double depth, const Vec3& U) {
        if (!d_.empty() && depth <= d_.back()) throw std::invalid_argument("current profile: depths must be strictly increasing");
        d_.push_back(depth); u_.push_back(U);
    }
    // Power law U(d) = Us ((h - d)/h)^alpha for 0 <= d <= h (1/7 power law: alpha = 0.143); zero below the seabed depth h.
    void setPowerLaw(const Vec3& Us, double alpha, double waterDepth) { Us_ = Us; alpha_ = alpha; h_ = waterDepth; power_ = true; }
    bool empty() const { return d_.empty() && !power_; }
    double surfaceZ() const { return z0_; }

    Vec3 velocity(double z) const {
        const double d = z0_ - z;
        if (power_) {
            if (d <= 0.0) return Us_;
            if (d >= h_) return Vec3();
            return Us_ * std::pow((h_ - d) / h_, alpha_);
        }
        if (d_.empty()) return Vec3();
        if (d <= d_.front()) return u_.front();
        if (d >= d_.back()) return u_.back();
        const size_t k = std::upper_bound(d_.begin(), d_.end(), d) - d_.begin();
        const double f = (d - d_[k - 1]) / (d_[k] - d_[k - 1]);
        return u_[k - 1] * (1.0 - f) + u_[k] * f;
    }

private:
    double z0_;
    std::vector<double> d_;
    std::vector<Vec3> u_;
    bool power_{false};
    Vec3 Us_;
    double alpha_{0.143}, h_{1.0};
};

}  // namespace mooring
