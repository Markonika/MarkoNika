// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Linearisation about the static state. Usage: mooring_modes config.json [key=value ...] [--lines N]
//   single-line config (has "anchor_m"): natural frequencies of the line
//   platform config (has "lines"): static equilibrium, linearised mooring stiffness, body natural frequencies and the lowest modes of every line
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include "mooring/case_runner.hpp"
#include "mooring/coupled_runner.hpp"
#include "mooring/radiation.hpp"

using namespace mooring;

static void printLine(const char* name, LumpedMassCable& cab, int n) {
    const LineModes m = cab.naturalModes(n);
    std::printf("line %s: %d degrees of freedom, K asymmetry %.1e; lowest %zu natural modes (dynamic stiffness):\n", name, m.dof, m.asymmetry, m.omega.size());
    std::printf("  mode   omega [rad/s]   f [Hz]   period [s]   zeta [-]\n");
    for (size_t k = 0; k < m.omega.size(); ++k)
        std::printf("  %4zu   %12.4f   %7.4f   %10.4f   %8.5f%s\n", k + 1, m.omega[k], m.omega[k] / 6.283185307179586, m.omega[k] > 0 ? 6.283185307179586 / m.omega[k] : 0.0, m.zeta[k], m.omega2[k] <= 0 ? "  (slack / unstable)" : "");
}

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: mooring_modes config.json [key=value ...] [--lines N]\n"); return 2; }
    std::ifstream in(argv[1]);
    if (!in) { std::fprintf(stderr, "cannot open %s\n", argv[1]); return 2; }
    nlohmann::json cfg = nlohmann::json::parse(in);
    int nLine = 8;
    for (int i = 2; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--lines") && i + 1 < argc) nLine = std::atoi(argv[++i]);
        else applyOverride(cfg, argv[i]);
    }
    if (cfg.contains("lines")) {
        std::vector<std::string> names;
        auto sysPtr = buildCoupledSystem(cfg, &names);
        CoupledSystem& sys = *sysPtr;
        const nlohmann::json jinit = cfg.value("initial", nlohmann::json::object());
        if (jinit.value("equilibrium", false)) {
            const EquilibriumResult er = sys.solveEquilibrium(jinit.value("equilibrium_tol", 1e-6));
            std::printf("equilibrium: %s in %d iterations (|F| %.2e N, |M| %.2e N m)\n", er.converged ? "converged" : "NOT converged", er.iterations, er.forceResidual, er.momentResidual);
        } else {
            sys.staticResidual();
        }
        const Mat6 K = sys.stiffness();
        std::printf("total linearised stiffness K = -dG/dxi [N/m, N/rad, N m/m, N m/rad] (rows/cols: surge sway heave rx ry rz):\n");
        for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j) std::printf(" %11.4f", K(i, j)); std::printf("\n"); }
        const RadiationModel* rad = sys.radiation();
        const BodyModes bm = rad ? bodyModesRadiation(sys.body().totalMass() - rad->Ainf(), *rad, K) : bodyModes(sys.body().totalMass(), K);
        std::printf("body natural frequencies (undamped, %s), K asymmetry %.1e:\n", rad ? "frequency-dependent added mass A(omega) at each mode" : "total mass incl. added mass", bm.asymmetry);
        const char* dn[6] = {"surge", "sway", "heave", "roll", "pitch", "yaw"};
        for (int k = 0; k < 6; ++k) {
            int dom = 0; for (int i = 1; i < 6; ++i) if (std::fabs(bm.shape[k][i]) > std::fabs(bm.shape[k][dom])) dom = i;
            std::printf("  mode %d: omega %9.4f rad/s, period %9.4f s  (dominant %s)\n", k + 1, bm.omega[k], bm.period[k], dn[dom]);
        }
        for (int k = 0; k < sys.numLines(); ++k) {
            auto* lm = dynamic_cast<LumpedMassCable*>(&sys.line(k));
            if (lm) printLine(names[k].c_str(), *lm, nLine);
        }
    } else {
        auto cab = buildStaticLine(cfg);
        std::printf("static top tension %.4f N\n", norm(cab->endTension(true)));
        printLine("(single)", *cab, nLine);
    }
    return 0;
}
