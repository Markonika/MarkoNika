#pragma once

namespace mfs {

// Physical / numerical constants shared across the solver.
namespace constants {
    inline constexpr double gravity   = 9.81;          // [m/s^2]
    inline constexpr double pi        = 3.14159265358979323846;
    inline constexpr double tiny      = 1.0e-12;        // guard against division by zero
    // Minimum clamp for volume fractions (and, via PipeGeometry, cross-
    // sectional area fractions). Kept well above machine epsilon: as a
    // layer's cross-section is pinched toward zero (e.g. the thin gas
    // layer just before a liquid slug fully bridges the pipe), several
    // closures (interfacial friction, drag, the pressure-equation
    // compressibility coefficient) become numerically stiff long before
    // the fraction is actually zero. 1e-3 keeps a "thin but resolvable"
    // layer -- physically consistent with the paper's own flow-regime
    // identification thresholds (c1, c2 ~ 0.02-0.15, Table 4) for what
    // counts as a vanishing gas layer -- while keeping the segregated
    // pressure-correction step well conditioned.
    inline constexpr double small_e   = 1.0e-3;
}

} // namespace mfs
