// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#pragma once
#include <functional>
#include <vector>
#include "mooring/cable_model.hpp"
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

enum class Scheme { RK4, Verlet };

// Prescribed top-end motion: fills position and velocity at time t.
using TopMotion = std::function<void(double t, Vec3& pos, Vec3& vel)>;

struct DynOptions {
    Scheme scheme{Scheme::RK4};
    double cfl{0.5};     // dt = cfl * l0 / c, c = sqrt(EA/m_l)
    double dt{0.0};      // if > 0, overrides the CFL value
};

// Counters for the slack-state regularisation (never hidden, see docs/assumptions.md).
struct DynStats {
    long steps{0};
    long slackSegmentEvals{0};   // segment force evaluations with eps <= 0 (tension set to zero)
    long clippedTensionEvals{0}; // evaluations where damping made T < 0 and T was clipped to 0
    double dtUsed{0};
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
class LumpedMassCable : public CableModel {
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
    // Tension force vector the end segment exerts on the end node (points into the line); this is the
    // force the line applies to the anchor / fairlead.
    Vec3 endTension(bool top) const;
    // Total force the line transmits to its support at that end: net force on the end node
    // (segment tension + the end node's half-weight, + hydrodynamics from milestone 3). Summed over
    // both ends in static equilibrium it equals the line's weight.
    Vec3 endForce(bool top) const;

    RelaxResult relaxStatic(const RelaxOptions& opt = {});

    // ---- dynamics (milestone 2) ----
    void setDynOptions(const DynOptions& o) { dyn_ = o; }
    double stableDt() const;                       // dt that will be used (CFL or override)
    void setTopMotion(TopMotion m) { top_ = std::move(m); }
    void setAnchor(const Vec3& a) { r_[0] = a; }
    void setInitialState(const std::vector<Vec3>& r, const std::vector<Vec3>& v, double t0 = 0.0);
    // Advance to time tEnd with fixed sub-steps (last step shortened to land exactly on tEnd).
    void advanceTo(double tEnd);
    double time() const { return t_; }
    const std::vector<Vec3>& velocities() const { return v_; }
    // Kinetic + elastic (eps>0) + gravitational energy [J]; excludes work done by the boundaries.
    double energy() const;
    const DynStats& stats() const { return stats_; }

    // CableModel: sub-steps to t with the fairlead moving linearly from its previous state, returns
    // the force the line exerts on the fairlead. Explicit/partitioned coupling.
    Vec3 forceOnBody(const Vec3& fairleadPos, const Vec3& fairleadVel, double t) override;

private:
    void acceleration(std::vector<Vec3>& r, std::vector<Vec3>& v, double t, std::vector<Vec3>& a) const;
    void step(double dt);

    CableParams p_;
    std::vector<Vec3> r_, v_;
    double t_{0.0};
    DynOptions dyn_;
    TopMotion top_;
    mutable DynStats stats_;
    Vec3 fairPrevPos_, fairPrevVel_;
    double fairPrevT_{0.0};
    bool fairInit_{false};
};

}  // namespace mooring
