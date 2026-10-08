// Copyright (c) [OWNER NAME]. All rights reserved. Proprietary and confidential.
#pragma once
#include <vector>
#include "mooring/vec3.hpp"

namespace mooring {

struct CableParams {
    double L{0};         // unstretched length [m]
    int    N{0};         // number of segments (N+1 nodes)
    double EA{0};        // axial stiffness [N]
    double m_l{0};       // mass per unit length [kg/m]
    double w{0};         // weight per unit length acting in -z [N/m] (submerged or dry, set by caller)
    double c_int{0};     // internal damping coefficient [N s], T = EA*eps + c_int*d(eps)/dt
    double g{9.81};
    double l0() const { return L / N; }
    double waveSpeed() const;   // c = sqrt(EA / m_l)
};

struct RelaxOptions {
    double forceTol{1e-8};     // stop when max nodal residual < forceTol * w * l0
    double massFactor{32.0};   // fictitious nodal mass = massFactor * EA / l0 (dt = 1)
    long   maxSteps{2000000};  // pseudo-time step limit
};

struct RelaxResult {
    bool converged{false};
    long steps{0};
    double dt{0};              // pseudo-time step (fictitious-mass relaxation)
    double maxResidual{0};     // [N]
    long kineticResets{0};
};

// Lumped-mass perfectly flexible cable (Paredes 2016, Ch. 3). Milestone 1: elastic + weight forces
// and static equilibrium by dynamic relaxation with kinetic damping.
class LumpedMassCable {
public:
    LumpedMassCable(const CableParams& p, const Vec3& anchor, const Vec3& fairlead);

    const CableParams& params() const { return p_; }
    const std::vector<Vec3>& nodes() const { return r_; }
    std::vector<Vec3>& nodes() { return r_; }

    // Net force on each node from segment tension and weight (end nodes included; ends are held
    // fixed by the caller). Tension is zero in compression (bilinear, Eq. 3.25).
    void computeForces(const std::vector<Vec3>& r, const std::vector<Vec3>& v, std::vector<Vec3>& f) const;
    double segmentTension(const std::vector<Vec3>& r, int seg) const;
    double nodeMass(int i) const;
    // Force exerted by the line on the end node's support: -(net force on the end node incl. weight)
    // is not what is wanted; this returns the segment tension vector pulling the end node into the line.
    Vec3 endTension(bool top) const;

    RelaxResult relaxStatic(const RelaxOptions& opt = {});

private:
    CableParams p_;
    std::vector<Vec3> r_;
};

}  // namespace mooring
