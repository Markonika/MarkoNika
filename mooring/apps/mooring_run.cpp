// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Usage: mooring_run config.json [key.path=value ...]
#include <cstdio>
#include <fstream>
#include "mooring/case_runner.hpp"

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: %s config.json [key.path=value ...]\n", argv[0]); return 2; }
    std::ifstream in(argv[1]);
    if (!in) { std::fprintf(stderr, "cannot open %s\n", argv[1]); return 2; }
    nlohmann::json cfg = nlohmann::json::parse(in);
    for (int i = 2; i < argc; ++i) mooring::applyOverride(cfg, argv[i]);
    const mooring::CaseResult r = mooring::runCase(cfg);
    std::printf("static top tension %.4f N (anchor %.4f N), relax %s in %ld steps\n", r.staticTopTension,
                r.staticAnchorTension, r.staticConverged ? "converged" : "NOT converged", r.staticSteps);
    std::printf("dt %.3e s, steps %ld, wall %.1f s, finite=%d\n", r.dt, r.steps, r.wallSeconds, int(r.finite));
    std::printf("slack evals %ld, clipped %ld, soil contact evals %ld\n", r.slackSegmentEvals, r.clippedTensionEvals, r.soilContactEvals);
    std::printf("mean of cycle maxima (%d cycles): %.4f N\n", r.cyclesAveraged, r.meanMax);
    return r.finite ? 0 : 1;
}
