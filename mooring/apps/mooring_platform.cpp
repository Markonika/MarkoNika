// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Usage: mooring_platform config.json [key.path=value ...]
#include <cstdio>
#include <fstream>
#include "mooring/case_runner.hpp"
#include "mooring/coupled_runner.hpp"

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: %s config.json [key.path=value ...]\n", argv[0]); return 2; }
    std::ifstream in(argv[1]);
    if (!in) { std::fprintf(stderr, "cannot open %s\n", argv[1]); return 2; }
    nlohmann::json cfg = nlohmann::json::parse(in);
    for (int i = 2; i < argc; ++i) mooring::applyOverride(cfg, argv[i]);
    const mooring::CoupledResult r = mooring::runCoupledCase(cfg);
    if (r.equilibriumRequested)
        std::printf("equilibrium: %s in %d iterations (|F| %.2e N, |M| %.2e N m)\n", r.equilibrium.converged ? "converged" : "NOT converged",
                    r.equilibrium.iterations, r.equilibrium.forceResidual, r.equilibrium.momentResidual);
    std::printf("xi(0) = [%.5f %.5f %.5f | %.5f %.5f %.5f]\n", r.xiStart[0], r.xiStart[1], r.xiStart[2], r.xiStart[3], r.xiStart[4], r.xiStart[5]);
    std::printf("dt_body %.3e s, shortest line step %.3e s, sub-step ratio %d, %ld body steps, wall %.1f s, finite=%d\n", r.report.dtBody,
                r.report.dtLineMin, r.report.subStepRatio, r.report.bodySteps, r.wallSeconds, int(r.finite));
    return r.finite ? 0 : 1;
}
