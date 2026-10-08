// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#pragma once
#include "mooring/vec3.hpp"

namespace mooring {

// Elastic catenary of a uniform cable hanging freely between two points in the vertical x-z plane
// (Paredes 2016, Eqs. 3.54-3.56). s is the unstretched (Lagrangian) arc length from end A.
// The tangent tension components along the cable are (H, V0 + w s).
struct ElasticCatenary {
    double L{0};    // unstretched length [m]
    double EA{0};   // axial stiffness [N]
    double w{0};    // weight per unit unstretched length [N/m] (acts in -z)
    double H{0};    // horizontal tension [N]
    double V0{0};   // vertical tension component at s = 0 [N]
    bool converged{false};
    double residual{0};   // final position residual [m]

    // Solve for (H, V0) so that the cable runs from A to B = A + (dx, dz), dx > 0.
    static ElasticCatenary solve(double L, double EA, double w, double dx, double dz);

    // Position relative to A at Lagrangian coordinate s in [0, L] -> (x, z).
    void position(double s, double& x, double& z) const;
    // Tension magnitude at s.
    double tension(double s) const;
};

}  // namespace mooring
