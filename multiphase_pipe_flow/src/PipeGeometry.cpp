#include "mfs/PipeGeometry.hpp"
#include "mfs/Constants.hpp"

#include <algorithm>
#include <cmath>

namespace mfs {

PipeGeometry::PipeGeometry(double diameter)
    : diameter_(diameter), area_(constants::pi * diameter * diameter / 4.0) {}

double PipeGeometry::chordWidth(double h1) const {
    const double h = std::clamp(h1, 0.0, diameter_);
    const double w = 2.0 * std::sqrt(std::max(0.0, h * (diameter_ - h)));
    return w;
}

double PipeGeometry::heightFromArea(double A1) const {
    const double A1c = std::clamp(A1, 1.0e-9 * area_, area_ * (1.0 - 1.0e-9));

    // theta(h) = 2*acos(1 - 2h/D);  A1(h) = D^2/8 * (theta - sin(theta))
    // Monotonic increasing in h on [0, D]. Safeguarded Newton with a
    // bisection bracket fallback for robustness at the extremes, where
    // dA1/dh -> 0.
    double lo = 0.0, hi = diameter_;
    double h = diameter_ * (A1c / area_); // linear initial guess

    auto areaAt = [&](double hh) {
        const double arg = std::clamp(1.0 - 2.0 * hh / diameter_, -1.0, 1.0);
        const double theta = 2.0 * std::acos(arg);
        return diameter_ * diameter_ / 8.0 * (theta - std::sin(theta));
    };

    for (int it = 0; it < 60; ++it) {
        const double f = areaAt(h) - A1c;
        if (f > 0.0) hi = h; else lo = h;

        const double dfdh = chordWidth(h); // dA1/dh
        double hNewton = h;
        if (dfdh > 1.0e-9 * diameter_) {
            hNewton = h - f / dfdh;
        }
        if (hNewton > lo && hNewton < hi) {
            h = hNewton;
        } else {
            h = 0.5 * (lo + hi);
        }
        if (std::fabs(hi - lo) < 1.0e-10 * diameter_) break;
    }
    return std::clamp(h, 0.0, diameter_);
}

PipeGeometry::StratifiedGeometry PipeGeometry::fromAreaFraction(double e1) const {
    const double e1c = std::clamp(e1, constants::small_e, 1.0 - constants::small_e);
    const double A1 = e1c * area_;
    const double h1 = heightFromArea(A1);

    const double arg = std::clamp(1.0 - 2.0 * h1 / diameter_, -1.0, 1.0);
    const double theta = 2.0 * std::acos(arg);

    StratifiedGeometry g{};
    g.h1 = h1;
    g.A1 = diameter_ * diameter_ / 8.0 * (theta - std::sin(theta));
    g.A2 = std::max(area_ - g.A1, 1.0e-12 * area_);
    g.Swp1 = theta * diameter_ / 2.0;
    g.Swp2 = std::max(constants::pi * diameter_ - g.Swp1, 1.0e-9 * diameter_);
    g.Si = chordWidth(h1);
    g.Si = std::max(g.Si, 1.0e-6 * diameter_);

    g.D1 = 4.0 * g.A1 / std::max(g.Swp1, 1.0e-9 * diameter_);
    g.D2 = 4.0 * g.A2 / (g.Swp2 + g.Si);
    return g;
}

} // namespace mfs
