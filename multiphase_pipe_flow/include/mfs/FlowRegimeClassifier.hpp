#pragma once

#include "mfs/FlowState.hpp"

#include <string>
#include <vector>

namespace mfs {

enum class FlowRegime { Stratified, Annular, Slug, Bubbly };

std::string toString(FlowRegime r);

// Section 5.1 identification criteria, used purely as a diagnostic / for
// comparison against published flow-regime maps. These criteria play NO
// role in the calculation itself (no closure in Closures.hpp/.cpp may
// depend on this classification) -- that separation is the central point
// of the paper.
struct RegimeClassifierOptions {
    double c1 = 0.05; // threshold on |d(eg)/dz|-based fluctuation amplitude (slug bridging)
    double c2 = 0.04; // threshold on eg itself (near-zero continuous-gas layer)
    double bubblyRunLengthFraction = 0.95; // fraction of pipe length continuously bubbly => "bubbly"
};

// Classifies each cell over a time window of eg(z,t) history. `egHistory`
// is a vector of eg snapshots (each of size N) taken at successive times;
// the fluctuation used for the slug/annular criteria is estimated from the
// spread of eg across that history at each cell.
std::vector<FlowRegime> classifyFromHistory(const std::vector<std::vector<double>>& egHistory,
                                             const RegimeClassifierOptions& opts = {});

// Single-snapshot classification using a supplied local fluctuation
// estimate deg(z) (e.g. a running standard deviation maintained by the
// caller) alongside the instantaneous eg(z).
std::vector<FlowRegime> classifyInstant(const std::vector<double>& eg,
                                         const std::vector<double>& degFluctuation,
                                         const RegimeClassifierOptions& opts = {});

} // namespace mfs
