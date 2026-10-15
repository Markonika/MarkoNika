// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Rigid-body platform with 6 DOF and multi-line fairlead coupling (milestone 7).
#pragma once
#include <array>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include "mooring/cable_model.hpp"
#include "mooring/linearize.hpp"
#include "mooring/parallel.hpp"
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
Mat6 operator-(const Mat6& a, const Mat6& b);
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

// Implicit (monolithic-in-time) body integration with iterated coupling: Newmark-beta for the body, Newton iterations on the pose with the
// approximate Jacobian M/(beta dt^2) + gamma/(beta dt) B + K (K = static linearised stiffness incl. hydrostatics, refreshed when the iteration
// is slow), every iteration re-advancing all lines from the saved state at the start of the step with the trial fairlead motion. Removes the
// explicit stability limit dt_body < 2/omega_s of stiff lines on light bodies.
struct ImplicitCoupling {
    bool enabled{false};
    double newmarkGamma{0.5};      // 0.5 = average acceleration (no numerical damping); > 0.5 damps high frequencies (beta = (gamma + 1/2)^2 / 4)
    double tol{1e-9};              // converged when the Newton correction is below tol (m, rad) x (1 + |xi|)
    int maxIter{20};
};

struct CouplingReport {
    double dtBody{0};
    double dtLineMin{0};       // smallest internal line time step [s]
    int subStepRatio{0};       // ceil(dtBody / dtLineMin): line steps per body step
    long bodySteps{0};
    long couplingIterations{0};  // implicit coupling: line re-advances in total
    long stiffnessBuilds{0};     // implicit coupling: static stiffness evaluations
    long couplingFailures{0};    // implicit coupling: steps accepted without convergence
};

// Explicit partitioned coupling: one body, any number of lines (conventional serial staggered scheme).
// Natural frequencies of a body with total mass matrix M and stiffness K (see linearize.hpp).
BodyModes bodyModes(const Mat6& M, const Mat6& K);

class RadiationModel;
class RadiationMemory;

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
    // 'fixed[k]' holds DOF k at its current value (e.g. surge imposed) while the others are solved.
    EquilibriumResult solveEquilibrium(double tol = 1e-6, int maxIter = 40, const std::array<bool, 6>& fixed = {});
    // Residual generalised load at the current pose (lines solved statically), N and N m.
    Vec6 staticResidual();
    // Total linearised stiffness K = -d(residual)/d(xi) at the current pose (hydrostatic C included), by central differences of the static line
    // solutions with step h (rad / m); leaves the pose and the lines unchanged. 'lineStiffness' excludes the body's own C. Rows/columns: surge, sway,
    // heave, rx, ry, rz.
    Mat6 stiffness(double h = 1e-4);
    Mat6 lineStiffness(double h = 1e-4);

    // Time-dependent generalised load on the body (e.g. wave excitation, Eq. 3.60); added to the loads at every step.
    void setExcitation(std::function<Vec6(double t)> f) { excitation_ = std::move(f); }

    // Start the time loop at time t0 from the current pose and velocity (lines assumed at rest/consistent).
    void begin(double t0 = 0.0);
    // One body step of size dt: velocity-Verlet for the body with implicit linear damping, lines sub-stepped
    // inside forceOnBody with the fairlead motion linearly interpolated over the step.
    void step(double dt);                                  // explicit or implicit coupling, see setImplicitCoupling
    void setImplicitCoupling(const ImplicitCoupling& c) { ic_ = c; Kvalid_ = false; }
    void advanceTo(double tEnd, double dt);
    double time() const { return t_; }
    // Advance / solve the lines concurrently on n threads (1 = serial, the default). The lines share no state, so the results are identical to the
    // serial ones bit for bit; it pays for 2+ lines when the line work per body step is large (small time steps, fine discretisation, static solves).
    // Frequency-dependent added mass / radiation damping (Cummins): the body mass matrix must already contain A_inf (BodyParams::A); the memory
    // force is integrated with the explicit coupling (not available together with the implicit coupling).
    void setRadiation(std::shared_ptr<const RadiationModel> m);
    const RadiationModel* radiation() const { return radModel_.get(); }
    void setThreads(int n);
    int threads() const { return pool_ ? pool_->concurrency() : 1; }
    CouplingReport report() const { return report_; }

private:
    struct Slot { std::unique_ptr<CableModel> line; Vec3 a; std::string name; Vec3 force; };
    Vec6 externalLoad(const Vec6& xi, const std::vector<Vec3>& lineForces) const;
    std::function<Vec6(double)> excitation_;
    RigidBody6DOF body_;
    std::vector<Slot> lines_;
    std::unique_ptr<WorkerPool> pool_;
    ImplicitCoupling ic_;
    std::shared_ptr<const RadiationModel> radModel_;
    std::shared_ptr<RadiationMemory> radMem_;
    Mat6 Kapprox_;
    bool Kvalid_{false};
    void stepImplicit(double dt);
    template <class F> void forLines(F&& f) {
        const int n = int(lines_.size());
        if (pool_ && n > 1) pool_->run(n, std::function<void(int)>(std::forward<F>(f)));
        else for (int i = 0; i < n; ++i) f(i);
    }
    double t_{0.0};
    Vec6 acc_{};
    Vec6 Fnow_{};          // line + constant + stiffness load at the current state (no damping terms)
    CouplingReport report_;
};

}  // namespace mooring
