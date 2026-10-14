// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Linearisation about the static state: natural frequencies / mode shapes of a line, linearised mooring stiffness and natural frequencies of the body.
#pragma once
#include <array>
#include <vector>
#include "mooring/vec3.hpp"

namespace mooring {

// Eigen-decomposition of a real symmetric matrix (row-major n x n) by the cyclic Jacobi method. Eigenvalues ascending; v[k*n + i] is component i of
// eigenvector k (orthonormal).
struct SymEigen { int n{0}; std::vector<double> w, v; };
SymEigen symmetricEigen(std::vector<double> a, int n);

// Natural modes of a line about its current (static) shape: K phi = omega^2 M phi with K = -M dA/dr (symmetrised), M the nodal inertia (mass + added
// mass in the normal direction), ends fixed. 'zeta' is the modal damping ratio phi^T C phi / (2 omega) from the linearised velocity-dependent forces
// (internal damping and the linearised Morison drag at rest, which is zero for quadratic drag). A mode with omega^2 <= 0 (slack / buckled part) is
// reported with omega = 0.
struct LineModes {
    std::vector<double> omega2, omega, zeta;
    std::vector<std::vector<Vec3>> shape;      // per mode, nodes 0..N, mass-normalised (phi^T M phi = 1)
    double asymmetry{0};                       // |K - K^T| / |K| before symmetrisation (finite-difference noise + non-conservative terms)
    int dof{0};
};

// Natural frequencies of the body: K x = omega^2 M x (6 x 6, symmetrised K).
struct BodyModes {
    std::array<double, 6> omega2{}, omega{}, period{};
    std::array<std::array<double, 6>, 6> shape{};   // shape[k] = mode k, mass-normalised
    double asymmetry{0};
};

}  // namespace mooring
