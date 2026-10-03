// timing.hpp - wall clock, timestamps and the per-stage timers.
#pragma once

#include <chrono>
#include <cstdio>
#include <ctime>
#include <string>

namespace dhj {

inline double now_seconds() {
    using clock = std::chrono::steady_clock;
    static const clock::time_point t0 = clock::now();
    return std::chrono::duration<double>(clock::now() - t0).count();
}

// ISO-8601 local timestamp, seconds resolution (matches the Python logs).
inline std::string timestamp_now() {
    const std::time_t t = std::time(nullptr);
    char buf[32];
    std::strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%S", std::localtime(&t));
    return buf;
}

inline std::string timestamp_compact() {
    const std::time_t t = std::time(nullptr);
    char buf[32];
    std::strftime(buf, sizeof buf, "%Y%m%d_%H%M%S", std::localtime(&t));
    return buf;
}

// Cumulative wall time per stage. The refinement loop reads the ode/queries/
// value_iteration counters before and after each phase for phases.csv.
struct Timers {
    double grid_build = 0.0;
    double spatial_index = 0.0;
    double cell_init = 0.0;
    double ode = 0.0;
    double queries = 0.0;
    double value_iteration = 0.0;
    double refine_split = 0.0;
    double plotting = 0.0;
    double checkpoint = 0.0;
    double total = 0.0;
    std::size_t vi_iterations = 0;
    std::size_t ode_solves = 0;
    std::size_t final_cells = 0;

    void report(const char* label) const {
        std::printf("\n%s\n", std::string(70, '=').c_str());
        std::printf("TIMING BREAKDOWN (%s)\n", label);
        std::printf("%s\n", std::string(70, '=').c_str());
        std::printf("  grid construction      : %10.3f s\n", grid_build);
        std::printf("  spatial index builds   : %10.3f s\n", spatial_index);
        std::printf("  cell initialization    : %10.3f s\n", cell_init);
        std::printf("  ODE integration        : %10.3f s  (%zu solves)\n", ode, ode_solves);
        std::printf("  successor set queries  : %10.3f s\n", queries);
        std::printf("  value iteration        : %10.3f s  (%zu sweeps)\n", value_iteration, vi_iterations);
        std::printf("  cell splitting         : %10.3f s\n", refine_split);
        std::printf("  plotting / output      : %10.3f s\n", plotting);
        std::printf("  checkpoint writes      : %10.3f s\n", checkpoint);
        std::printf("  %s\n", std::string(48, '-').c_str());
        std::printf("  TOTAL                  : %10.3f s\n", total);
        std::printf("  final leaf cells       : %10zu\n", final_cells);
    }
};

}  // namespace dhj
