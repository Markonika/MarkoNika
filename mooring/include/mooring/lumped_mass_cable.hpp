// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#pragma once
#include <functional>
#include <vector>
#include "mooring/cable_model.hpp"
#include "mooring/hydro.hpp"
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
    // Hydrodynamics (Eqs. 3.27-3.29). A1 = 0 means pi/4 D0^2.
    double D0{0}, D1{0}, A1{0};
    double Cm{0}, Cdt{0}, Cdn{0};
    SoilParams soil;
    bool planar{false};   // 2D mode: motion confined to the plane y = y_anchor (forces and velocities in y are dropped)
    double l0() const { return L / N; }
    double dryWeight() const { return m_l * g; }          // weight per length above the surface
    double nominalArea() const { return A1 > 0.0 ? A1 : 0.7853981633974483 * D0 * D0; }
    // Submerged weight per unit length, Eq. 3.26.
    static double submergedWeight(double m_l, double rho_c, double rho_w, double g = 9.81) {
        return (rho_c - rho_w) / rho_c * m_l * g;
    }
    double waveSpeed() const;   // c = sqrt(EA / m_l)
};

// Water kinematics at a point: velocity and acceleration of the water. Default is still water.
using WaterField = std::function<void(const Vec3& pos, double t, Vec3& vw, Vec3& aw)>;

struct Environment {
    bool hydro{false};          // Morison forces + submerged weight below the surface (else CableParams::w everywhere)
    bool seabed{false};         // flat seabed contact at z = seabedZ
    double rho_w{1000.0};
    double surfaceZ{0.0};       // still-water level: nodes above carry dry weight and no hydrodynamics
    double seabedZ{0.0};
    WaterField water;           // empty = still water
};

struct RelaxOptions {
    double forceTol{1e-8};     // stop when max nodal residual < forceTol * w * l0
    double massFactor{32.0};   // fictitious nodal mass = massFactor * EA / l0 (dt = 1)
    long   maxSteps{2000000};  // pseudo-time step limit
};

// Point element lumped on a node: clump weight, floater, instrument... Weight m g acts always; buoyancy
// rho_w g V, drag and added mass act on the submerged fraction of the node (requires Environment::hydro).
// Drag: 0.5 rho_w Cd A |v_rel| v_rel (isotropic). Added mass: Cm rho_w V (a_w - a), isotropic; the
// Froude-Krylov term is omitted, consistently with the cable's Eq. 3.27. No seabed contact of its own
// (the node's cable-soil contact still acts).
struct PointElement {
    int node{0};
    double mass{0};     // dry mass [kg]
    double volume{0};   // displaced volume [m^3]
    double Cd{0};
    double area{0};     // drag reference area [m^2]
    double Cm{0};

    // Floater from its buoyancy force [N] (net upward in water when only the volume is counted) and diameter (sphere).
    static PointElement floater(int node, double mass, double buoyancyN, double D, double Cd, double rho_w = 1000.0,
                                double g = 9.81, double Cm = 0.0) {
        return {node, mass, buoyancyN / (rho_w * g), Cd, 0.7853981633974483 * D * D, Cm};
    }
    // Clump weight from its mass and submerged weight [N]: V = (m g - W_sub) / (rho_w g).
    static PointElement clump(int node, double mass, double submergedWeightN, double D, double Cd, double rho_w = 1000.0,
                              double g = 9.81, double Cm = 0.0) {
        return {node, mass, (mass * g - submergedWeightN) / (rho_w * g), Cd, 0.7853981633974483 * D * D, Cm};
    }
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
    long soilContactEvals{0};    // node evaluations in seabed contact
    double dtUsed{0};
};

class LumpedMassCable;
using StepObserver = std::function<void(const LumpedMassCable&)>;

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
    // Elastic + weight + (optional) hydrodynamic + seabed forces. If 'ca' / 'aw' are given they receive
    // the per-node added-mass coefficient [kg] and water acceleration used by the added-mass solve.
    void computeForces(const std::vector<Vec3>& r, const std::vector<Vec3>& v, double t,
                       std::vector<Vec3>& f, std::vector<double>* ca = nullptr,
                       std::vector<Vec3>* aw = nullptr) const;
    void setEnvironment(const Environment& e) { env_ = e; }
    // Attach a point element to an interior or end node (0..N).
    void addPointElement(const PointElement& pe);
    const std::vector<PointElement>& pointElements() const { return points_; }
    // Net vertical force of the point elements on node i in the current configuration (weight - buoyancy), N, down positive.
    double pointNetWeight(const std::vector<Vec3>& r, int i) const;
    const Environment& environment() const { return env_; }
    // Weight [N] on node i for the given configuration (submergence-blended dry / submerged weight).
    double nodeWeight(const std::vector<Vec3>& r, int i) const;
    // Fraction of node i below the still-water level (0 = air, 1 = fully submerged).
    double submergedFraction(const std::vector<Vec3>& r, int i) const;
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
    // Start shape for a line that lies on a flat seabed and rises to the fairlead: inextensible
    // catenary with horizontal touchdown tangent (suspended part) plus a straight lying part. Needs
    // the anchor on the seabed and Environment::seabed. Returns false (keeps the current shape) when the
    // geometry has no such solution. Intended as the first guess for relaxStatic().
    bool initTouchdownCatenary();
    // Raw tension magnitude [N] in the end segment, including the c_int damping term (as in the force
    // evaluation, clipped at 0). Uses the current velocities.
    double endSegmentTension(bool top) const;
    void setStepObserver(StepObserver o) { observer_ = std::move(o); }

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
    Vec3 staticForceOnBody(const Vec3& fairleadPos) override;
    double internalTimeStep() const override { return stableDt(); }
    // Move the top node (e.g. to a new body position) without touching the other nodes; velocities are reset.
    void setFairlead(const Vec3& pos);

private:
    void acceleration(std::vector<Vec3>& r, std::vector<Vec3>& v, double t, std::vector<Vec3>& a) const;
    void step(double dt);

    CableParams p_;
    std::vector<Vec3> r_, v_;
    double t_{0.0};
    DynOptions dyn_;
    TopMotion top_;
    mutable DynStats stats_;
    Environment env_;
    std::vector<PointElement> points_;
    StepObserver observer_;
    double ksCap_{1e300};      // soil stiffness cap used only during static relaxation
    Vec3 fairPrevPos_, fairPrevVel_;
    double fairPrevT_{0.0};
    bool fairInit_{false};
};

}  // namespace mooring
