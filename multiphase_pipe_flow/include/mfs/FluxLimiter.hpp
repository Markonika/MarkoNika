#pragma once

#include <algorithm>
#include <cmath>

namespace mfs {

// TVD flux limiters, Sweby (1984) psi(r) form. r is the ratio of the
// upwind-side gradient to the local (across-face) gradient; each limiter
// satisfies psi(1)=1 (recovers second-order central/upwind-biased
// reconstruction for a locally linear profile) and psi(r)=0 for r<=0 (falls
// back to first-order upwind at a local extremum, which is exactly what
// keeps the scheme Total-Variation-Diminishing -- no new oscillation can be
// created at a front or a sign change in the gradient).
enum class FluxLimiterType { None, Minmod, VanLeer, Superbee, MC };

inline double fluxLimiterPsi(FluxLimiterType type, double r) {
    switch (type) {
        case FluxLimiterType::None:
            return 0.0; // always first-order upwind
        case FluxLimiterType::Minmod:
            return std::max(0.0, std::min(1.0, r));
        case FluxLimiterType::VanLeer:
            return (r > 0.0) ? (r + std::fabs(r)) / (1.0 + std::fabs(r)) : 0.0;
        case FluxLimiterType::Superbee:
            return std::max({0.0, std::min(2.0 * r, 1.0), std::min(r, 2.0)});
        case FluxLimiterType::MC:
            return std::max(0.0, std::min({2.0 * r, 0.5 * (1.0 + r), 2.0}));
    }
    return 0.0;
}

} // namespace mfs
