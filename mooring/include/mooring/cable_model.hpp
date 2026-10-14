// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#pragma once
#include <stdexcept>
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
    // Static equilibrium of the line with the fairlead held at fairleadPos; returns the force on the fairlead.
    // Leaves the line at rest in that state (used for body equilibrium and as the start of a run).
    virtual Vec3 staticForceOnBody(const Vec3& fairleadPos) = 0;
    // Internal (stable) time step of the line [s]; the coupling reports dt_body / this value as the sub-step ratio.
    virtual double internalTimeStep() const = 0;
    // Snapshot / rollback of the physical state (needed by implicit coupling, which re-advances the line from the start of the step on every
    // iteration). The default implementation does not support it.
    virtual void saveState() { throw std::logic_error("this cable model does not support saveState()"); }
    virtual void restoreState() { throw std::logic_error("this cable model does not support restoreState()"); }
};

}  // namespace mooring
