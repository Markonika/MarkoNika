// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
#pragma once
#include <condition_variable>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace mooring {

// Minimal persistent thread pool for "run f(0..n-1) and wait" (used to advance independent mooring lines concurrently). The calling thread takes part,
// so a pool of concurrency c spawns c-1 workers. run() returns after every index has been processed; the first exception thrown by f is rethrown.
class WorkerPool {
public:
    explicit WorkerPool(int concurrency) {
        for (int i = 1; i < concurrency; ++i) workers_.emplace_back([this] { loop(); });
    }
    ~WorkerPool() {
        { std::lock_guard<std::mutex> lk(m_); stop_ = true; ++gen_; }
        cv_.notify_all();
        for (std::thread& t : workers_) t.join();
    }
    WorkerPool(const WorkerPool&) = delete;
    WorkerPool& operator=(const WorkerPool&) = delete;
    int concurrency() const { return int(workers_.size()) + 1; }

    void run(int n, const std::function<void(int)>& f) {
        if (n <= 0) return;
        if (workers_.empty() || n == 1) { for (int i = 0; i < n; ++i) f(i); return; }
        {
            std::lock_guard<std::mutex> lk(m_);
            f_ = &f; n_ = n; next_ = 0; active_ = int(workers_.size()); err_ = nullptr; ++gen_;
        }
        cv_.notify_all();
        work();                                              // the caller helps
        std::unique_lock<std::mutex> lk(m_);
        done_.wait(lk, [this] { return active_ == 0; });
        f_ = nullptr;
        if (err_) std::rethrow_exception(err_);
    }

private:
    void work() {
        for (;;) {
            int i;
            { std::lock_guard<std::mutex> lk(m_); if (next_ >= n_) return; i = next_++; }
            try { (*f_)(i); } catch (...) { std::lock_guard<std::mutex> lk(m_); if (!err_) err_ = std::current_exception(); }
        }
    }
    void loop() {
        unsigned seen = 0;
        for (;;) {
            { std::unique_lock<std::mutex> lk(m_); cv_.wait(lk, [&] { return gen_ != seen; }); seen = gen_; if (stop_) return; }
            work();
            { std::lock_guard<std::mutex> lk(m_); --active_; }
            done_.notify_one();
        }
    }
    std::vector<std::thread> workers_;
    std::mutex m_;
    std::condition_variable cv_, done_;
    const std::function<void(int)>* f_{nullptr};
    int n_{0}, next_{0}, active_{0};
    unsigned gen_{0};
    bool stop_{false};
    std::exception_ptr err_;
};

}  // namespace mooring
