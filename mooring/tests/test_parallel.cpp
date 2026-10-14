// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// WorkerPool and concurrent line solves (CoupledSystem::setThreads): identical results to the serial run.
#include <atomic>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <vector>
#include "doctest.h"
#include "mooring/coupled_runner.hpp"
#include "mooring/parallel.hpp"

using namespace mooring;
using json = nlohmann::json;

TEST_CASE("WorkerPool: every index exactly once, repeated runs, exception propagation") {
    WorkerPool pool(4);
    CHECK(pool.concurrency() == 4);
    for (int rep = 0; rep < 200; ++rep) {
        std::vector<int> hit(7, 0);
        pool.run(7, [&](int i) { ++hit[i]; });
        for (int h : hit) CHECK(h == 1);
    }
    std::atomic<long> sum{0};
    pool.run(1000, [&](int i) { sum += i; });
    CHECK(sum == 999L * 1000L / 2);
    CHECK_THROWS_AS(pool.run(5, [](int i) { if (i == 3) throw std::runtime_error("boom"); }), std::runtime_error);
    std::vector<int> again(3, 0);                         // the pool is still usable after an exception
    pool.run(3, [&](int i) { again[i] = 1; });
    CHECK(std::accumulate(again.begin(), again.end(), 0) == 3);
    WorkerPool serial(1);
    int n = 0; serial.run(4, [&](int) { ++n; });
    CHECK(n == 4);
}

TEST_CASE("Concurrent line solves give results identical to the serial run (three-leg platform, equilibrium + time loop)") {
    std::ifstream in(std::string(MOORING_SOURCE_DIR) + "/examples/platform/three_leg_example.json");
    REQUIRE(in.good());
    json c = json::parse(in);
    c["output"]["write"] = false;
    c["numerics"]["t_end_s"] = 1.0;
    c["numerics"]["threads"] = 1;
    const CoupledResult a = runCoupledCase(c);
    c["numerics"]["threads"] = 3;
    const CoupledResult b = runCoupledCase(c);
    REQUIRE(a.finite);
    REQUIRE(b.finite);
    REQUIRE(a.samples.size() == b.samples.size());
    bool same = true;
    for (size_t k = 0; k < a.samples.size(); ++k) {
        for (int i = 0; i < 6; ++i) same = same && a.samples[k].xi[i] == b.samples[k].xi[i];
        for (size_t j = 0; j < a.samples[k].lineForce.size(); ++j) same = same && a.samples[k].lineForce[j] == b.samples[k].lineForce[j];
    }
    CHECK(same);
    for (int i = 0; i < 6; ++i) CHECK(a.xiStart[i] == b.xiStart[i]);
}
