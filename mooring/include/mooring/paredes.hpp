// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Benchmark helpers for the Paredes (2016) free-buoy and 3-leg moored-buoy experiments (milestone 8).
#pragma once
#include <vector>
#include "mooring/coupled_runner.hpp"
#include "nlohmann/json.hpp"

namespace mooring {

struct ParedesStatics {
    EquilibriumResult eq;
    Vec6 xi{};
    double draftChange{0};                 // positive = deeper than the free-floating draft [m]
    std::vector<double> topTension;        // end-segment tension of each leg at the fairlead [N]
    std::vector<double> fairleadForce;     // |force on the buoy| of each leg [N]
    Vec3 netLineForce;                     // sum of the line forces on the buoy [N]
};
ParedesStatics paredesStatics(const nlohmann::json& cfg);

struct ParedesStiffness {
    double F0{0};                          // net surge force of the lines at the free equilibrium [N]
    std::vector<double> x, Fx, K;          // imposed surge [m], line force in x [N], secant stiffness -(Fx - F0)/x [N/m]
    std::vector<bool> converged;
};
// Surge held at x (other 5 DOF free); the stiffness is the restoring force change over the displacement.
ParedesStiffness paredesSurgeStiffness(const nlohmann::json& cfg, const std::vector<double>& xs);

// Mooring contribution to the heave stiffness (+-delta, N/m) and pitch stiffness (+-delta rad, N m/rad), other DOF free.
double paredesHeaveStiffness(const nlohmann::json& cfg, double delta = 0.01);
double paredesPitchStiffness(const nlohmann::json& cfg, double deltaRad = 0.017453292519943295);

// Damped period [s] of DOF 'dof' from a free decay released from 'offset' (crest/trough method of the thesis,
// sec. 5.5): time between the first and last crossing divided by the cycles in between. 0 if fewer than 3 crossings.
double paredesDecayPeriod(nlohmann::json cfg, int dof, double offset, double tEnd, double dtBody);

// Regular-wave response of a Paredes configuration: first-order amplitudes (least-squares fit at the wave frequency over the last
// 'fitCycles' wave periods) of surge, heave and pitch, normalised as in the thesis (Figs 5.22-5.24): surge/a, heave/a, pitch/(k a).
// Body coefficients A, B, wave force w, delta are those of Table 3.4 for T = 1.30 or 1.40 s (including the surge-pitch coupling).
struct ParedesRAO { double surge{0}, heave{0}, pitch{0}; double a{0}, k{0}; bool finite{true}; double meanSurge{0}; };
// steadyStart: initialise pose and velocity with the steady response of the free body (closed form) instead of starting from rest with a ramp.
ParedesRAO paredesWaveRAO(nlohmann::json cfg, double period, double height, double tEnd, int fitCycles, double dtBody = 2e-3, bool steadyStart = true);

}  // namespace mooring
