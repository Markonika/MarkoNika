// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// JSON-configured single-line case: static relaxation, prescribed top-end motion, raw output.
#pragma once
#include <string>
#include <vector>
#include "nlohmann/json.hpp"

namespace mooring {

struct CycleMax { int cycle; double topMax, anchorMax; };

struct CaseResult {
    double staticTopTension{0};     // [N] after static relaxation
    double staticAnchorTension{0};
    bool staticConverged{false};
    long staticSteps{0};
    double dt{0};                   // time step used [s]
    double l0{0};
    double cWave{0};
    long steps{0};
    long slackSegmentEvals{0}, clippedTensionEvals{0}, soilContactEvals{0};
    std::vector<CycleMax> cycles;   // per excitation cycle, from every time step
    double meanMax{0};              // mean of cycle maxima over cycles >= statistics.first_cycle
    int cyclesAveraged{0};
    double maxTopOverall{0};
    bool finite{true};
    double wallSeconds{0};
};

// Run a case. Output files (time series, cycle maxima, parameter log) are written to
// cfg["output"]["directory"] with the prefix cfg["output"]["tag"] when "write" is true (default).
CaseResult runCase(const nlohmann::json& cfg);

// Apply "a.b.c=value" overrides (value parsed as JSON, falling back to a string).
void applyOverride(nlohmann::json& cfg, const std::string& assignment);

}  // namespace mooring
