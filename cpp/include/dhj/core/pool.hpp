// pool.hpp - persistent worker pool with a chunked parallel_for.
//
// The Python package uses multiprocessing.Pool; here one process with N
// persistent threads plays the same role. Work is split into one contiguous
// chunk per worker.
#pragma once

#include <condition_variable>
#include <cstddef>
#include <functional>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

namespace dhj {

class ThreadPool {
public:
    explicit ThreadPool(std::size_t n_workers) : n_(n_workers ? n_workers : 1) {
        ranges_.resize(n_);
        workers_.reserve(n_ - 1);
        for (std::size_t i = 1; i < n_; ++i) workers_.emplace_back([this, i] { worker_loop(i); });
    }

    ~ThreadPool() {
        {
            std::lock_guard<std::mutex> lk(m_);
            stop_ = true;
            ++generation_;
        }
        cv_start_.notify_all();
        for (auto& t : workers_) if (t.joinable()) t.join();
    }

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    std::size_t size() const { return n_; }

    // Calls fn(begin, end) on disjoint contiguous ranges covering [0, n) and
    // blocks until every range has been processed.
    void parallel_for(std::size_t n, const std::function<void(std::size_t, std::size_t)>& fn) {
        if (n == 0) return;
        const std::size_t nw = n_ < n ? n_ : n;
        if (nw <= 1) { fn(0, n); return; }

        const std::size_t base = n / nw, rem = n % nw;
        {
            std::lock_guard<std::mutex> lk(m_);
            std::size_t start = 0;
            for (std::size_t i = 0; i < nw; ++i) {
                const std::size_t len = base + (i < rem ? 1 : 0);
                ranges_[i] = {start, start + len};
                start += len;
            }
            fn_ = &fn;
            active_ = nw;
            remaining_ = nw - 1;  // the calling thread handles range 0 itself
            ++generation_;
        }
        cv_start_.notify_all();

        fn(ranges_[0].first, ranges_[0].second);

        std::unique_lock<std::mutex> lk(m_);
        cv_done_.wait(lk, [this] { return remaining_ == 0; });
        fn_ = nullptr;
    }

private:
    void worker_loop(std::size_t id) {
        std::size_t seen = 0;
        for (;;) {
            std::pair<std::size_t, std::size_t> range;
            const std::function<void(std::size_t, std::size_t)>* fn = nullptr;
            {
                std::unique_lock<std::mutex> lk(m_);
                cv_start_.wait(lk, [&] { return generation_ != seen; });
                seen = generation_;
                if (stop_) return;
                if (id >= active_) continue;
                range = ranges_[id];
                fn = fn_;
            }
            if (fn && range.first < range.second) (*fn)(range.first, range.second);
            {
                std::lock_guard<std::mutex> lk(m_);
                if (--remaining_ == 0) cv_done_.notify_one();
            }
        }
    }

    const std::size_t n_;
    std::vector<std::thread> workers_;
    std::vector<std::pair<std::size_t, std::size_t>> ranges_;
    std::mutex m_;
    std::condition_variable cv_start_, cv_done_;
    const std::function<void(std::size_t, std::size_t)>* fn_ = nullptr;
    std::size_t generation_ = 0;
    std::size_t remaining_ = 0;
    std::size_t active_ = 0;
    bool stop_ = false;
};

}  // namespace dhj
