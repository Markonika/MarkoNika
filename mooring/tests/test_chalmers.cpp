// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Test 4 (Chalmers 33 m chain, Bergdahl, Eskilsson & Palm 2016): static tension and mean maximum top tension
// against chalmers_table7_max_tension.csv (+/-5 % reading error). The full 30-case grid is a separate ctest entry.
#include <atomic>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <thread>
#include <vector>
#include "doctest.h"
#include "mooring/case_runner.hpp"

using namespace mooring;
using json = nlohmann::json;

namespace {
const std::string kDir = std::string(MOORING_SOURCE_DIR) + "/examples/chalmers/";

json baseConfig() {
    std::ifstream in(kDir + "chalmers_config.json");
    REQUIRE(in.good());
    json c = json::parse(in);
    c["output"]["write"] = false;
    return c;
}

struct Meas { double T, r, tension; };
std::vector<Meas> readTable() {
    std::ifstream in(kDir + "chalmers_table7_max_tension.csv");
    REQUIRE(in.good());
    std::string line; std::vector<double> radii; std::vector<Meas> out;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line); std::string tok; std::vector<std::string> f;
        while (std::getline(ss, tok, ',')) f.push_back(tok);
        if (radii.empty()) { for (size_t i = 1; i < f.size(); ++i) radii.push_back(std::stod(f[i].substr(2, f[i].size() - 4))); continue; }
        for (size_t i = 1; i < f.size(); ++i) out.push_back({std::stod(f[0]), radii[i - 1], std::stod(f[i])});
    }
    return out;
}

double simulate(double r, double T, int N = 66) {
    json c = baseConfig();
    c["motion"]["radius_m"] = r; c["motion"]["period_s"] = T; c["line"]["segments"] = N;
    return runCase(c).meanMax;
}
}  // namespace

TEST_CASE("Chalmers: table has 30 values") { CHECK(readTable().size() == 30); }

TEST_CASE("Chalmers: static top tension vs published 22.68 N, converging in N") {
    std::printf("\n  static top tension (published 22.68 N)\n  %6s %12s %10s\n", "N", "T_top [N]", "rel diff");
    double prev = 0, last = 0;
    for (int N : {16, 33, 66, 132}) {
        json c = baseConfig(); c["line"]["segments"] = N; c["motion"]["type"] = "none";
        const CaseResult r = runCase(c);
        std::printf("  %6d %12.4f %10.2e\n", N, r.staticTopTension, (r.staticTopTension - 22.68) / 22.68);
        CHECK(r.staticConverged);
        if (N > 16) CHECK(std::fabs(r.staticTopTension - last) < std::fabs(last - prev) + 1e-9);   // differences shrink
        prev = last; last = r.staticTopTension;
    }
    CHECK(last == doctest::Approx(22.68).epsilon(0.01));
}

TEST_CASE("Chalmers: representative cases within the +/-5 % reading error (slow, mid, snap)") {
    struct C { double r, T, ref; };
    const std::vector<C> cs = {{0.2, 3.5, 50.1}, {0.1, 2.0, 39.5}, {0.2, 1.25, 70.3}};
    for (const C& c : cs) {
        const double s = simulate(c.r, c.T);
        std::printf("  r=%.3f T=%.2f  measured %.1f  simulated %.2f  (%+.1f %%)\n", c.r, c.T, c.ref, s, 100 * (s - c.ref) / c.ref);
        CHECK(std::fabs(s - c.ref) / c.ref < 0.05);
    }
}

TEST_CASE("Chalmers full grid: 30 cases vs Table 7 [validation, run via ctest -R chalmers_grid]" * doctest::skip()) {
    const auto table = readTable();
    std::vector<double> sim(table.size());
    std::atomic<size_t> next{0};
    std::vector<std::thread> pool;
    for (int k = 0; k < 4; ++k) pool.emplace_back([&] {
        for (size_t i; (i = next++) < table.size();) sim[i] = simulate(table[i].r, table[i].T);
    });
    for (auto& t : pool) t.join();
    double mx = 0, my = 0; const double n = double(table.size());
    for (size_t i = 0; i < table.size(); ++i) { mx += table[i].tension; my += sim[i]; }
    mx /= n; my /= n;
    double sxx = 0, syy = 0, sxy = 0, ssr = 0, maxRel = 0; int within5 = 0;
    for (size_t i = 0; i < table.size(); ++i) {
        const double x = table[i].tension, y = sim[i];
        sxx += (x - mx) * (x - mx); syy += (y - my) * (y - my); sxy += (x - mx) * (y - my); ssr += (y - x) * (y - x);
        maxRel = std::max(maxRel, std::fabs(y - x) / x); within5 += std::fabs(y - x) / x <= 0.05;
    }
    const double r2 = sxy * sxy / (sxx * syy), r2_11 = 1.0 - ssr / sxx;
    std::printf("\n  grid: regression r2 = %.4f (published MooDy DG: 0.98), slope %.3f, r2 vs 1:1 = %.4f, RMSE %.2f N, "
                "mean bias %+.2f N, within 5%%: %d/30, max |rel| %.3f\n", r2, sxy / sxx, r2_11, std::sqrt(ssr / n), my - mx, within5, maxRel);
    CHECK(r2 > 0.97);
    CHECK(r2_11 > 0.95);
    CHECK(within5 >= 22);
    CHECK(maxRel < 0.10);
}
