// run.hpp - Algorithm 1 / Algorithm 2 drivers for a concrete dynamics object.
#pragma once

#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>

#include "../dynamics/base.hpp"
#include "args.hpp"
#include "cell.hpp"
#include "checkpoint.hpp"
#include "params.hpp"
#include "pool.hpp"
#include "reach.hpp"
#include "refine.hpp"
#include "timing.hpp"
#include "types.hpp"
#include "vi.hpp"

namespace dhj {

inline std::string derived_json(double eta_min, const GronwallReachabilityAnalyzer& reach, int workers,
                                const std::string& out_dir, const std::string& indent) {
    const std::string in = indent + "  ";
    std::string j = "{\n";
    j += in + "\"eta_min\": " + py_repr(eta_min) + ",\n";
    j += in + "\"n_checkpoints\": " + std::to_string(reach.n_checkpoints()) + ",\n";
    j += in + "\"checkpoint_times\": [";
    for (std::size_t i = 0; i < reach.checkpoint_times().size(); ++i) j += std::string(i ? ", " : "") + py_repr(reach.checkpoint_times()[i]);
    j += "],\n" + in + "\"growth_factors\": [";
    for (std::size_t i = 0; i < reach.growth_factors().size(); ++i) j += std::string(i ? ", " : "") + py_repr(reach.growth_factors()[i]);
    j += "],\n" + in + "\"workers\": " + std::to_string(workers) + ",\n";
    j += in + "\"out_dir\": \"" + json_escape(out_dir) + "\"\n";
    j += indent + "}";
    return j;
}

template <std::size_t N>
void run_algorithm_1(const Args& args, const Dynamics<N>& dyn, const std::string& out_dir, ThreadPool& pool, Timers& timers) {
    std::printf("%s\n", std::string(70, '=').c_str());
    std::printf("ALGORITHM 1: Discretization Routine\n");
    std::printf("%s\n", std::string(70, '=').c_str());
    std::printf("Environment: %s\n", dyn.name());
    std::printf("Initializing grid with resolution %d^%zu...\n", args.resolution, N);

    const double grid0 = now_seconds();
    CellTree<N> tree(dyn.state_bounds(), args.resolution, dyn.periodic_mask());
    timers.grid_build += now_seconds() - grid0;
    timers.spatial_index += tree.index_build_seconds();

    GronwallReachabilityAnalyzer reach(dyn.L_f(), dyn.dt(), dyn.tau());
    write_params(out_dir, args, dyn.describe_json("  "), derived_json(0.0, reach, args.workers, out_dir, "  "));
    SafetyValueIterator<N> vi(args, dyn, tree, reach, out_dir, pool, timers);

    const double t0 = now_seconds();
    vi.value_iteration(args.iterations, args.tolerance, args.plot_freq, args.conservative, args.delta_max);
    vi.print_statistics();
    std::printf("\nALGORITHM 1 COMPLETE\n");
    std::printf("Total time: %.2f seconds\n", now_seconds() - t0);
    std::printf("Results saved to: %s/\n", out_dir.c_str());
    timers.final_cells = tree.num_leaves();
}

template <std::size_t N>
void run_algorithm_2(const Args& args, const Dynamics<N>& dyn, const std::string& out_dir, ThreadPool& pool, Timers& timers) {
    std::printf("%s\n", std::string(70, '=').c_str());
    std::printf("ALGORITHM 2/3: Adaptive Refinement\n");
    std::printf("%s\n", std::string(70, '=').c_str());
    std::printf("Environment: %s\n", dyn.name());

    const std::string ckpt_dir = out_dir + "/checkpoints";
    std::unique_ptr<CellTree<N>> tree;
    CheckpointMeta<N> resume_meta;
    bool resuming = false;

    if (!args.resume.empty()) {
        std::string path = args.resume;
        if (path == "auto") {
            path = latest_checkpoint(ckpt_dir);
            if (path.empty()) throw std::runtime_error("--resume auto: no checkpoint found in " + ckpt_dir);
            std::printf("--resume auto selected %s\n", path.c_str());
        }
        const double load0 = now_seconds();
        resume_meta = load_checkpoint<N>(path, args.mode, dyn.periodic_mask(), tree);
        timers.grid_build += now_seconds() - load0;
        timers.spatial_index += tree->index_build_seconds();
        resuming = true;
        std::printf("Loaded checkpoint %s: phase %u, %zu leaves\n", path.c_str(), resume_meta.phase, tree->num_leaves());
        if (resume_meta.initial_resolution != static_cast<std::uint32_t>(args.initial_resolution))
            std::fprintf(stderr, "  warning: checkpoint was made with --initial-resolution %u but this run passes %d; "
                                 "the checkpoint's grid wins\n", resume_meta.initial_resolution, args.initial_resolution);
    } else {
        std::printf("Initializing grid with resolution %d^%zu...\n", args.initial_resolution, N);
        const double grid0 = now_seconds();
        tree.reset(new CellTree<N>(dyn.state_bounds(), args.initial_resolution, dyn.periodic_mask()));
        timers.grid_build += now_seconds() - grid0;
        timers.spatial_index += tree->index_build_seconds();
    }

    GronwallReachabilityAnalyzer reach(dyn.L_f(), dyn.dt(), dyn.tau());
    write_params(out_dir, args, dyn.describe_json("  "),
                 derived_json(args.epsilon / (2.0 * dyn.L_l()), reach, args.workers, out_dir, "  "));
    AdaptiveRefinement<N> adaptive(args, dyn, *tree, reach, out_dir, pool, timers);

    const double t0 = now_seconds();
    adaptive.refine(args.epsilon, args.refinements, args.vi_iterations, resuming ? &resume_meta : nullptr);
    std::printf("ALGORITHM 2 COMPLETE\n");
    std::printf("Total time: %.2f seconds\n", now_seconds() - t0);
    timers.final_cells = tree->num_leaves();
}

// Entry point once the dynamics type is known (called from the registry).
template <std::size_t N>
int run_with_dynamics(const Args& args, const Dynamics<N>& dyn, const std::string& out_dir) {
    dyn.validate();
    std::printf("Dynamics: %s (%zu-D)\n", dyn.name(), N);
    std::printf("  Time discretization: %d steps of dt=%g s over tau=%g s\n", dyn.n_steps(), dyn.dt(), dyn.tau());

    ThreadPool pool(static_cast<std::size_t>(args.workers));
    Timers timers;
    const double t0 = now_seconds();
    try {
        if (args.algorithm == 1) run_algorithm_1<N>(args, dyn, out_dir, pool, timers);
        else run_algorithm_2<N>(args, dyn, out_dir, pool, timers);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    timers.total = now_seconds() - t0;
    timers.report(mode_script(args.mode));
    std::printf("TOTAL WALL TIME: %.2f seconds\n", timers.total);
    return 0;
}

}  // namespace dhj
