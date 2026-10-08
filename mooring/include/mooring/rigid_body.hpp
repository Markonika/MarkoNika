// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Rigid-body platform with 6 DOF and multi-line fairlead coupling (milestone 7).
#pragma once
#include <array>
#include <memory>
#include <string>
#include <vector>
#include "mooring/cable_model.hpp"
#include "mooring/vec3.hpp"

namespace mooring {

using Vec6 = std::array<double, 6>;
using Mat3 = std::array<double, 9>;      // row-major
struct Mat6 {
    std::array<double, 36> a{};
    double& operator()(int i, int j) { return a[6 * i + j]; }
    double operator()(int i, int j) const { return a[6 * i + j]; }
    static Mat6 diag(const Vec6& d) { Mat6 m; for (int i = 0; i < 6; ++i) m(i, i) = d[i]; return m; }
};
Vec6 operator*(const Mat6& m, const Vec6& v);
Mat6 operator+(const Mat6& a, const Mat6& b);
Mat6 operator*(const Mat6& m, double s);
// Solve m x = b by Gaussian elimination with partial pivoting; throws std::runtime_error if singular.
Vec6 solve6(Mat6 m, Vec6 b);

Mat3 rotationFromVector(const Vec3& theta);      // Rodrigues / exponential map, exactly orthonormal
Vec3 rotate(const Mat3& R, const Vec3& v);

// Linear (frequency-independent) potential-flow model, Paredes (2016) Eq. 3.59:
//   (M + A) xi'' + B xi' + C xi + Dq |xi'| xi' = F0 + F_lines(xi, xi', t)
// xi = (surge, sway, heave, rx, ry, rz) = CG translation and rotation vector measured from the reference
// (free-floating equilibrium) pose; the rotational part uses omega ~ d(theta)/dt (small angles). Fairlead
// kinematics use the exact exponential map R(theta), so fairlead positions are right for moderate rotations,
// while the dynamics are the linearised ones; valid up to roughly 15-20 deg.
struct BodyParams {
    double mass{1.0};
    Mat3 inertia{1, 0, 0, 0, 1, 0, 0, 0, 1};   // about the CG, reference axes [kg m^2]
    Mat6 A, B, C;                               // added mass, damping, hydrostatic (+ other linear) stiffness
    Vec6 Dq{};                                  // quadratic drag per DOF: F_k = -Dq_k |v_k| v_k [N s^2/m^2 or N m s^2/rad^2]
    Vec6 F0{};                                  // constant generalised external load (e.g. mean wave drift) [N, N m]
    Vec3 cgRef;                                 // global position of the CG at xi = 0
};

class RigidBody6DOF {
public:
    explicit RigidBody6DOF(const BodyParams& p);
    const BodyParams& params() const { return p_; }
    Mat6 totalMass() const { return M_; }       // M_rb + A
    Vec6 xi, xiDot;

    Vec3 cg() const { return p_.cgRef + Vec3(xi[0], xi[1], xi[2]); }
    // Global position / velocity of the body point at 'a' (body coordinates relative to the CG, reference axes).
    Vec3 pointPosition(const Vec3& a) const { return cg() + rotate(rotationFromVector(Vec3(xi[3], xi[4], xi[5])), a); }
    Vec3 pointVelocity(const Vec3& a) const;
    // Generalised load of a force f applied at body point a: (f, (R a) x f).
    static void addPointLoad(const Vec3& Ra, const Vec3& f, Vec6& F);
    Mat3 rotation() const { return rotationFromVector(Vec3(xi[3], xi[4], xi[5])); }
private:
    BodyParams p_;
    Mat6 M_;
};

struct EquilibriumResult { bool converged{false}; int iterations{0}; double forceResidual{0}, momentResidual{0}; };

struct CouplingReport {
    double dtBody{0};
    double dtLineMin{0};       // smallest internal line time step [s]
    int subStepRatio{0};       // ceil(dtBody / dtLineMin): line steps per body step
    long bodySteps{0};
};

// Explicit partitioned coupling: one body, any number of lines (conventional serial staggered scheme).
class CoupledSystem {
public:
    explicit CoupledSystem(const BodyParams& p) : body_(p) {}
    RigidBody6DOF& body() { return body_; }
    const RigidBody6DOF& body() const { return body_; }
    // Attach a line; 'a' is the fairlead in body coordinates relative to the CG (reference axes).
    // The line's anchor and its initial shape are the caller's responsibility.
    int addLine(std::unique_ptr<CableModel> line, const Vec3& a, const std::string& name = "");
    CableModel& line(int i) { return *lines_[i].line; }
    int numLines() const { return int(lines_.size()); }
    const Vec3& lineForce(int i) const { return lines_[i].force; }   // last force of line i on the fairlead
    Vec3 fairleadPosition(int i) const { return body_.pointPosition(lines_[i].a); }
    Vec3 fairleadBodyCoord(int i) const { return lines_[i].a; }

    // Static equilibrium: Newton on F0 - C xi + sum(line loads(xi)) = 0 with finite-difference Jacobian,
    // statics of each line solved exactly at every evaluation. Leaves body and lines at rest at the solution.
    EquilibriumResult solveEquilibrium(double tol = 1e-6, int maxIter = 40);
    // Residual generalised load at the current pose (lines solved statically), N and N m.
    Vec6 staticResidual();

    // Start the time loop at time t0 from the current pose and velocity (lines assumed at rest/consistent).
    void begin(double t0 = 0.0);
    // One body step of size dt: velocity-Verlet for the body with implicit linear damping, lines sub-stepped
    // inside forceOnBody with the fairlead motion linearly interpolated over the step.
    void step(double dt);
    void advanceTo(double tEnd, double dt);
    double time() const { return t_; }
    CouplingReport report() const { return report_; }

private:
    struct Slot { std::unique_ptr<CableModel> line; Vec3 a; std::string name; Vec3 force; };
    Vec6 externalLoad(const Vec6& xi, const std::vector<Vec3>& lineForces) const;
    RigidBody6DOF body_;
    std::vector<Slot> lines_;
    double t_{0.0};
    Vec6 acc_{};
    Vec6 Fnow_{};          // line + constant + stiffness load at the current state (no damping terms)
    CouplingReport report_;
};

}  // namespace mooring
