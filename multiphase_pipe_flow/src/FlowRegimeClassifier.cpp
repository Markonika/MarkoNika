#include "mfs/FlowRegimeClassifier.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace mfs {

std::string toString(FlowRegime r) {
    switch (r) {
        case FlowRegime::Stratified: return "stratified";
        case FlowRegime::Annular:    return "annular";
        case FlowRegime::Slug:       return "slug";
        case FlowRegime::Bubbly:     return "bubbly";
    }
    return "unknown";
}

std::vector<FlowRegime> classifyInstant(const std::vector<double>& eg,
                                         const std::vector<double>& degFluctuation,
                                         const RegimeClassifierOptions& opts) {
    const int N = static_cast<int>(eg.size());
    std::vector<FlowRegime> regime(N, FlowRegime::Stratified);

    for (int i = 0; i < N; ++i) {
        const bool lowFluctuation = degFluctuation[i] <= opts.c1;
        const bool nearZeroGasLayer = eg[i] <= opts.c2;

        if (lowFluctuation) {
            regime[i] = FlowRegime::Stratified;
        } else if (nearZeroGasLayer) {
            // Large void-fraction fluctuations AND a near-vanishing
            // continuous-gas layer: the interface is bridging the pipe --
            // slug flow (Section 5.1).
            regime[i] = FlowRegime::Slug;
        } else {
            // Large fluctuations but the continuous gas layer does not
            // (locally) vanish -- large-amplitude waves that do not bridge
            // the pipe, i.e. annular / large-wave stratified flow.
            regime[i] = FlowRegime::Annular;
        }
    }

    // Bubbly flow: pipe fully bridged with no regions of (non-thin)
    // stratified gas layer, over a long contiguous run. Approximate this
    // by checking that eg is uniformly very small (below c2) over a long
    // contiguous run of cells with LOW fluctuation (a persistent, not
    // transient, near-zero gas layer -- distinguishing it from the
    // transient near-zero gas layer seen momentarily under a slug body).
    int runStart = -1;
    for (int i = 0; i <= N; ++i) {
        const bool bubblyLike = (i < N) && (eg[i] <= opts.c2) && (degFluctuation[i] <= opts.c1);
        if (bubblyLike) {
            if (runStart < 0) runStart = i;
        } else {
            if (runStart >= 0) {
                const int runLen = i - runStart;
                if (static_cast<double>(runLen) / N >= opts.bubblyRunLengthFraction) {
                    for (int k = runStart; k < i; ++k) regime[k] = FlowRegime::Bubbly;
                }
            }
            runStart = -1;
        }
    }

    return regime;
}

std::vector<FlowRegime> classifyFromHistory(const std::vector<std::vector<double>>& egHistory,
                                             const RegimeClassifierOptions& opts) {
    if (egHistory.empty()) return {};
    const int N = static_cast<int>(egHistory.front().size());
    const int T = static_cast<int>(egHistory.size());

    std::vector<double> mean(N, 0.0), fluct(N, 0.0);
    for (const auto& snap : egHistory)
        for (int i = 0; i < N; ++i) mean[i] += snap[i];
    for (int i = 0; i < N; ++i) mean[i] /= T;

    for (const auto& snap : egHistory)
        for (int i = 0; i < N; ++i) fluct[i] += (snap[i] - mean[i]) * (snap[i] - mean[i]);
    for (int i = 0; i < N; ++i) fluct[i] = std::sqrt(fluct[i] / std::max(1, T - 1));

    return classifyInstant(mean, fluct, opts);
}

} // namespace mfs
