// refine.hpp - Algorithm 2/3: adaptive refinement of the boundary cells.
//
// Each phase: identify boundary cells larger than eta_min = eps / (2 L_l),
// bisect them along their widest dimension, rebuild the spatial index, run
// the local value iteration on the whole leaf set, save the figure, write a
// checkpoint and a row in phases.csv.
#pragma once

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "../dynamics/base.hpp"
#include "args.hpp"
#include "cell.hpp"
#include "checkpoint.hpp"
#include "params.hpp"
#include "reach.hpp"
#include "timing.hpp"
#include "types.hpp"
#include "vi.hpp"

namespace dhj {

template <std::size_t N>
class AdaptiveRefinement {
public:
    AdaptiveRefinement(const Args& args, const Dynamics<N>& dyn, CellTree<N>& tree,
                       const GronwallReachabilityAnalyzer& reach, const std::string& out_dir, ThreadPool& pool,
                       Timers& timers)
        : args_(args), dyn_(dyn), tree_(tree), timers_(timers), output_dir_(out_dir),
          ckpt_dir_(out_dir + "/checkpoints"), vi_(args, dyn, tree, reach, out_dir, pool, timers),
          phase_log_(out_dir), t_start_(now_seconds()) {}

    const std::string& output_dir() const { return output_dir_; }
    SafetyValueIterator<N>& vi() { return vi_; }

    // `resume` (may be null) is the loop state read from a checkpoint; the tree
    // then already holds converged values for `resume->phase` phases and phase 0
    // is skipped. `max_refinements` counts total phases, not additional ones.
    void refine(double epsilon, int max_refinements, int vi_iterations_per_refinement,
                const CheckpointMeta<N>* resume = nullptr) {
        const double eta_min = epsilon / (2.0 * dyn_.L_l());
        std::printf("\n%s\n", std::string(70, '=').c_str());
        std::printf("ADAPTIVE REFINEMENT CONFIGURATION\n");
        std::printf("%s\n", std::string(70, '=').c_str());
        std::printf("  Error tolerance eps: %g\n", epsilon);
        std::printf("  Minimum cell size eta_min: %.6f\n", eta_min);
        std::printf("  Maximum refinements: %d\n", max_refinements);
        std::printf("  VI iterations per refinement: %d\n", vi_iterations_per_refinement);
        std::printf("  Conservative mode: %s\n", py_bool(args_.conservative).c_str());
        if (args_.conservative) std::printf("  delta_max: %g\n", args_.delta_max);

        VIResult conv;
        int refinement_iter = 0;
        std::size_t total_refined = 0;

        if (resume) {
            refinement_iter = static_cast<int>(resume->phase);
            total_refined = static_cast<std::size_t>(resume->total_refined);
            vi_.set_refinement_phase(refinement_iter);
            std::printf("\n%s\n", std::string(70, '=').c_str());
            std::printf("RESUMING FROM CHECKPOINT\n");
            std::printf("%s\n", std::string(70, '=').c_str());
            std::printf("  Phases already completed: %d\n", refinement_iter);
            std::printf("  Cumulative refined:       %zu\n", total_refined);
            std::printf("  Leaf cells:               %zu\n", tree_.num_leaves());
            std::printf("  Successor cache will be rebuilt in full on the first local sweep (exact).\n");
        } else {
            std::printf("\n%s\n", std::string(70, '=').c_str());
            std::printf("PHASE 0: INITIAL VALUE ITERATION\n");
            std::printf("%s\n", std::string(70, '=').c_str());
            std::printf("Grid: %zu cells\n", tree_.num_leaves());
            vi_.set_refinement_phase(0);
            const double phase0_start = now_seconds();
            const double succ0 = timers_.ode + timers_.queries, vi0 = timers_.value_iteration;
            conv = vi_.value_iteration(vi_iterations_per_refinement, args_.phase0_tol, args_.plot_freq,
                                       args_.conservative, args_.delta_max);
            const double phase_time = now_seconds() - phase0_start;
            std::printf("total time for running value_iteration the first time (phase 0) on the initial grid including "
                        "grid initialization, setting up cpu parallelism, successor set computation and value iteration "
                        "is %g seconds\n", phase_time);
            if (!conv.conv_upper.empty()) {
                std::printf("Initial VI completed in %zu iterations\n", conv.conv_upper.size());
                std::printf("Final convergence: ||V_upper||_inf = %.8e, ||V_lower||_inf = %.8e\n", conv.conv_upper.back(),
                            conv.conv_lower.back());
            }
            std::printf("Saved initial visualization: %s/value_function_phase_0_complete.png\n", output_dir_.c_str());
            const auto boundary = vi_.identify_boundary_cells();
            const auto refinable = filter_refinable(boundary, eta_min);
            phase_log_.append(0, tree_.num_leaves(), boundary.size(), refinable.size(), 0, conv.iterations, conv.converged,
                              timers_.value_iteration - vi0, timers_.ode + timers_.queries - succ0, phase_time,
                              now_seconds() - t_start_);
            save_phase_checkpoint(0, 0, epsilon, false);
        }

        std::vector<std::uint32_t> boundary = vi_.identify_boundary_cells();
        std::vector<std::uint32_t> refinable = filter_refinable(boundary, eta_min);
        std::printf("\n%s\nINITIAL QUEUE STATE\n%s\n", std::string(70, '=').c_str(), std::string(70, '=').c_str());
        std::printf("Total boundary cells: %zu\n", boundary.size());
        std::printf("  Refinable (>eta_min=%.6f): %zu\n", eta_min, refinable.size());
        std::printf("  Below threshold: %zu\n", boundary.size() - refinable.size());

        while (refinement_iter < max_refinements && !refinable.empty()) {
            boundary = vi_.identify_boundary_cells();
            refinable = filter_refinable(boundary, eta_min);
            std::printf("\n%s\nREFINEMENT PHASE %d\n%s\n", std::string(70, '=').c_str(), refinement_iter + 1,
                        std::string(70, '=').c_str());
            const double phase_start = now_seconds();
            const double succ0 = timers_.ode + timers_.queries, vi0 = timers_.value_iteration;
            std::printf("Boundary cells: %zu\n", boundary.size());
            std::printf("  Refinable (>eta_min=%.6f): %zu\n", eta_min, refinable.size());
            std::printf("  Below threshold: %zu\n", boundary.size() - refinable.size());
            std::printf("  Cumulative refined: %zu\n", total_refined);
            if (refinable.empty()) {
                std::printf("  No refinable cells remaining - stopping refinement\n");
                break;
            }
            print_size_stats(refinable);

            std::printf("  Refining %zu cells...\n", refinable.size());
            const double split_start = now_seconds();
            const std::vector<std::uint32_t> new_cells = tree_.refine_cells(refinable);
            const double split_time = now_seconds() - split_start;
            timers_.refine_split += split_time;
            const std::size_t n_refined = refinable.size();
            total_refined += n_refined;
            std::printf("    Refined: %zu parent cells\n", n_refined);
            std::printf("    Created: %zu child cells\n", new_cells.size());
            std::printf("    Total cells: %zu\n", tree_.num_leaves());
            std::printf("    Refinement time: %.3fs\n", split_time);
            std::printf("    Cumulative refined: %zu\n", total_refined);

            std::printf("  Rebuilding spatial index...\n");
            const double index_start = now_seconds();
            tree_.rebuild_spatial_index();
            const double index_time = now_seconds() - index_start;
            timers_.spatial_index += index_time;
            std::printf("    Spatial index rebuilt in %.3fs\n", index_time);

            if (!new_cells.empty()) {
                std::printf("  Initializing %zu new cells...\n", new_cells.size());
                const double init_start = now_seconds();
                vi_.initialize_new_cells(new_cells);
                std::printf("    Initialization completed in %.3fs\n", now_seconds() - init_start);
            }

            vi_.set_refinement_phase(refinement_iter + 1);
            std::printf("  Starting local value iteration...\n");
            const double local_start = now_seconds();
            conv = vi_.local_value_iteration(new_cells, vi_iterations_per_refinement, args_.tolerance, args_.conservative,
                                             args_.delta_max);
            std::printf("  Local VI completed in %zu iterations, time: %.3fs\n", conv.conv_upper.size(),
                        now_seconds() - local_start);
            if (!conv.conv_upper.empty()) {
                std::printf("    Final ||V_upper^k - V_upper^k-1||_inf = %.20e\n", conv.conv_upper.back());
                std::printf("    Final ||V_lower^k - V_lower^k-1||_inf = %.20e\n", conv.conv_lower.back());
            }

            char buf[128];
            std::snprintf(buf, sizeof buf, "/value_function_phase_%d_complete.png", refinement_iter + 1);
            vi_.emit_output(output_dir_ + buf, refinement_iter + 1, refinement_iter + 1);
            std::printf("  Saved visualization: %s%s\n", output_dir_.c_str(), buf);

            const double phase_time = now_seconds() - phase_start;
            std::printf("  Phase %d completed in %.2fs with %zu leaf cells\n", refinement_iter + 1, phase_time,
                        tree_.num_leaves());
            ++refinement_iter;

            const auto boundary_after = vi_.identify_boundary_cells();
            const auto refinable_after = filter_refinable(boundary_after, eta_min);
            phase_log_.append(refinement_iter, tree_.num_leaves(), boundary_after.size(), refinable_after.size(),
                              n_refined, conv.iterations, conv.converged, timers_.value_iteration - vi0,
                              timers_.ode + timers_.queries - succ0, phase_time, now_seconds() - t_start_);
            save_phase_checkpoint(refinement_iter, total_refined, epsilon, false);
            // Like the paper scripts, `refinable` is re-evaluated at the top of
            // the next phase (which prints "No refinable cells remaining").
        }

        save_phase_checkpoint(refinement_iter, total_refined, epsilon, true);

        std::printf("\n%s\nADAPTIVE REFINEMENT COMPLETE\n%s\n", std::string(70, '=').c_str(), std::string(70, '=').c_str());
        std::printf("Refinement Summary:\n");
        std::printf("  Phases completed: %d\n", refinement_iter);
        std::printf("  Parent cells refined: %zu\n", total_refined);
        std::printf("  Child cells created: %zu\n", total_refined * 2);
        std::printf("  Final cells: %zu\n", tree_.num_leaves());
        const auto final_boundary = vi_.identify_boundary_cells();
        const auto final_refinable = filter_refinable(final_boundary, eta_min);
        std::printf("Final Boundary State:\n");
        std::printf("  Boundary cells: %zu\n", final_boundary.size());
        std::printf("  Still refinable: %zu\n", final_refinable.size());
        std::printf("  Below threshold: %zu\n", final_boundary.size() - final_refinable.size());
        vi_.print_statistics();
        std::printf("\n All results saved to: %s/\n", output_dir_.c_str());
    }

private:
    void save_phase_checkpoint(int phase, std::size_t total_refined, double epsilon, bool force) {
        if (args_.checkpoint_every <= 0) return;
        if (!force && phase % args_.checkpoint_every != 0) return;
        if (phase == last_saved_phase_) return;
        CheckpointMeta<N> meta;
        meta.phase = static_cast<std::uint32_t>(phase);
        meta.total_refined = static_cast<std::uint64_t>(total_refined);
        meta.initial_resolution = static_cast<std::uint32_t>(args_.initial_resolution);
        meta.epsilon = epsilon;
        meta.root_bounds = tree_.root_bounds();
        const std::string path = checkpoint_path(ckpt_dir_, phase);
        const double t0 = now_seconds();
        try {
            save_checkpoint<N>(path, args_.mode, tree_, meta);
        } catch (const std::exception& e) {
            std::fprintf(stderr, "  warning: checkpoint failed: %s\n", e.what());
            return;
        }
        const double el = now_seconds() - t0;
        timers_.checkpoint += el;
        last_saved_phase_ = phase;
        std::error_code ec;
        const auto bytes = std::filesystem::file_size(path, ec);
        std::printf("  Checkpoint written: %s (%.1f MB, %.2fs)\n", path.c_str(),
                    ec ? 0.0 : static_cast<double>(bytes) / (1024.0 * 1024.0), el);
        saved_phases_.push_back(phase);
        prune_checkpoints();
    }

    // Keeps the newest K checkpoint files written by THIS run; files left by an
    // earlier run (including the one resumed from) are not touched.
    void prune_checkpoints() {
        if (args_.keep_checkpoints <= 0) return;
        std::error_code ec;
        while (static_cast<int>(saved_phases_.size()) > args_.keep_checkpoints) {
            std::filesystem::remove(checkpoint_path(ckpt_dir_, saved_phases_.front()), ec);
            saved_phases_.erase(saved_phases_.begin());
        }
    }

    std::vector<std::uint32_t> filter_refinable(const std::vector<std::uint32_t>& boundary, double eta_min) const {
        std::vector<std::uint32_t> out;
        out.reserve(boundary.size());
        for (std::uint32_t id : boundary) if (tree_.cell(id).max_range() > eta_min) out.push_back(id);
        return out;
    }

    void print_size_stats(const std::vector<std::uint32_t>& ids) const {
        std::vector<double> sizes;
        sizes.reserve(ids.size());
        for (std::uint32_t id : ids) sizes.push_back(tree_.cell(id).max_range());
        std::sort(sizes.begin(), sizes.end());
        double mean = 0.0;
        for (double v : sizes) mean += v;
        mean /= static_cast<double>(sizes.size());
        const std::size_t m = sizes.size() / 2;
        const double median = sizes.size() % 2 ? sizes[m] : 0.5 * (sizes[m - 1] + sizes[m]);
        std::printf("  Refinable cell size statistics:\n");
        std::printf("    Max:    %.6f\n", sizes.back());
        std::printf("    Mean:   %.6f\n", mean);
        std::printf("    Min:    %.6f\n", sizes.front());
        std::printf("    Median: %.6f\n", median);
    }

    const Args& args_;
    const Dynamics<N>& dyn_;
    CellTree<N>& tree_;
    Timers& timers_;
    std::string output_dir_, ckpt_dir_;
    SafetyValueIterator<N> vi_;
    PhaseLog phase_log_;
    double t_start_;
    int last_saved_phase_ = -1;
    std::vector<int> saved_phases_;
};

}  // namespace dhj
