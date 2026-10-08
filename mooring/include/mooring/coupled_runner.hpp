// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// JSON-configured platform + multi-line case (milestone 7). See docs/config.md.
#pragma once
#include <memory>
#include <string>
#include <vector>
#include "mooring/rigid_body.hpp"
#include "nlohmann/json.hpp"

namespace mooring {

struct CoupledSample {
    double t;
    Vec6 xi;
    std::vector<double> lineForce;     // |force of line k on its fairlead| [N]
    std::vector<double> topTension;    // end-segment tension of line k at the fairlead [N] (lumped-mass lines)
};

struct CoupledResult {
    EquilibriumResult equilibrium;
    bool equilibriumRequested{false};
    Vec6 xiStart{};                    // pose at t = 0
    std::vector<CoupledSample> samples;
    CouplingReport report;
    std::vector<std::string> lineNames;
    bool finite{true};
    double wallSeconds{0};
};

// Build body + lines from the config (no equilibrium solve, no time stepping). Optionally returns the line names.
std::unique_ptr<CoupledSystem> buildCoupledSystem(const nlohmann::json& cfg, std::vector<std::string>* lineNames = nullptr);

// Run a case: build body and lines, optionally solve the static equilibrium, integrate to numerics.t_end_s.
CoupledResult runCoupledCase(const nlohmann::json& cfg);

}  // namespace mooring
