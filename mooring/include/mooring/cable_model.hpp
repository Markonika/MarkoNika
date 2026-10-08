// Copyright (c) [OWNER NAME]. All rights reserved. Proprietary and confidential.
#pragma once
#include "mooring/vec3.hpp"

namespace mooring {

// Abstract cable model. The lumped-mass solver implements it today; a spectral/hp
// discontinuous Galerkin model (Paredes 2016, Ch. 3) can be added behind the same interface.
class CableModel {
public:
    virtual ~CableModel() = default;
    // Advance the line to time t with the top end at (fairleadPos, fairleadVel) and return the
    // force the line exerts on the fairlead [N]. (Implemented from milestone 2.)
    virtual Vec3 forceOnBody(const Vec3& fairleadPos, const Vec3& fairleadVel, double t) = 0;
};

}  // namespace mooring
