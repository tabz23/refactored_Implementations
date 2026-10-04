// vi.hpp - successor-set computation and the interval value iteration.
//
// Port of SafetyValueIterator from the paper scripts. The three modes differ
// only in the Bellman backup:
//   RANoDiscount     V <- min(l, max(r, best))
//   RADiscount       V <- min(l, max(r, gamma * best))
//   AvoidNoDiscount  V <- min(l, best)
// where best = max over actions of (min / max over successors of V) and an
// empty successor set contributes -inf.
#pragma once

#include <algorithm>
#include <atomic>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "../dynamics/base.hpp"
#include "args.hpp"
#include "cell.hpp"
#include "pool.hpp"
#include "reach.hpp"
#include "svg.hpp"
#include "timing.hpp"
#include "types.hpp"

namespace dhj {

// CSR successor graph over leaf positions for one VI run; index n is the
// out-of-bounds sink.
struct SuccGraph {
    std::vector<std::uint32_t> off;
    std::vector<std::uint32_t> idx;
};

struct VIResult {
    std::vector<double> conv_upper, conv_lower;
    int iterations = 0;
    bool converged = false;
};

template <std::size_t N>
class SafetyValueIterator {
public:
    using State = dhj::State<N>;
    using Tree = CellTree<N>;

    SafetyValueIterator(const Args& args, const Dynamics<N>& dyn, Tree& tree, const GronwallReachabilityAnalyzer& reach,
                        std::string output_dir, ThreadPool& pool, Timers& timers)
        : mode_(args.mode), args_(args), dyn_(dyn), tree_(tree), reach_(reach), output_dir_(std::move(output_dir)),
          pool_(pool), timers_(timers), gamma_(args.gamma), n_actions_(dyn.num_actions()),
          periodic_(dyn.periodic_mask()),
          plotter_(dyn, PlotStyle{true, false, mode_has_target(args.mode), args.plot_dpi, args.plot_slice}) {
        if (args.plot_slice >= 0 && static_cast<std::size_t>(args.plot_slice) >= dyn.slices().size()) {
            std::fprintf(stderr, "error: --plot-slice %d is out of range (%zu slices)\n", args.plot_slice, dyn.slices().size());
            std::exit(2);
        }
        std::error_code ec;
        std::filesystem::create_directories(output_dir_, ec);
        std::printf("Initialized with gamma=%g, L_f=%g, L_l=%g,  L_r=%g\n", gamma_, dyn.L_f(), dyn.L_l(), dyn.L_r());
        std::printf("Contraction factor: gamma*L_f = %.4f\n", gamma_ * dyn.L_f());
        std::printf("Using %zu parallel workers\n", pool_.size());
        if (args_.precompute) precompute_all_successors();
    }

    const std::string& output_dir() const { return output_dir_; }
    void set_refinement_phase(int p) { refinement_phase_ = p; }
    bool successor_cache_empty() const { return cache_entries_ == 0; }

    // ------------------------------------------------------------ cells --

    void initialize_cells() { initialize(tree_.leaves()); }

    void initialize_new_cells(const std::vector<std::uint32_t>& new_cells) {
        if (new_cells.empty()) return;
        std::printf("  Initializing %zu new cells in parallel...\n", new_cells.size());
        const double t0 = now_seconds();
        initialize(new_cells);
        const double el = now_seconds() - t0;
        std::printf("   Initialized in %.2fs (%.1f cells/s)\n", el, el > 0 ? new_cells.size() / el : 0.0);
    }

    // ------------------------------------------------------- successors --

    void precompute_all_successors() {
        std::printf("\nPrecomputing successor sets with HYBRID approach...\n");
        const auto& leaves = tree_.leaves();
        std::printf("  Step 1: Computing %zu trajectories in parallel...\n", leaves.size() * n_actions_);
        compute_successors(leaves);
    }

    // Drop every cached entry naming a now-split parent, then recompute every
    // leaf with a missing action (port of _update_successor_cache_for_new_cells).
    void update_successor_cache(const std::vector<std::uint32_t>& new_cells) {
        if (new_cells.empty()) return;
        std::printf("    Updating successor cache for %zu new cells (HYBRID)...\n", new_cells.size());
        const double t0 = now_seconds();
        ensure_cache_capacity();
        std::vector<bool> refined_parent(tree_.num_cells(), false);
        for (std::uint32_t id : new_cells) {
            const std::int64_t p = tree_.cell(id).parent;
            if (p >= 0) refined_parent[static_cast<std::size_t>(p)] = true;
        }
        for (std::size_t key = 0; key < cache_.size(); ++key) {
            if (!cache_present_[key]) continue;
            const std::size_t cid = key / n_actions_;
            bool drop = refined_parent[cid];
            if (!drop)
                for (std::int32_t sid : cache_[key])
                    if (sid != kOOB && refined_parent[static_cast<std::size_t>(sid)]) { drop = true; break; }
            if (drop) {
                cache_[key].clear();
                cache_[key].shrink_to_fit();
                cache_present_[key] = 0;
                --cache_entries_;
            }
        }
        std::vector<std::uint32_t> affected;
        for (std::uint32_t id : tree_.leaves()) {
            bool missing = false;
            for (std::size_t a = 0; a < n_actions_; ++a)
                if (!cache_present_[static_cast<std::size_t>(id) * n_actions_ + a]) { missing = true; break; }
            if (missing) affected.push_back(id);
        }
        std::printf("      Affected cells: %zu\n", affected.size());
        compute_successors(affected);
        const double el = now_seconds() - t0;
        std::printf("     Cache updated in %.2fs (%.1f tasks/s)\n", el, el > 0 ? affected.size() * n_actions_ / el : 0.0);
    }

    // ---------------------------------------------------------------- VI --

    // Algorithm 1 / phase 0. The refinement loop passes the same residual
    // tolerance as later phases (--delta-min), unless --phase0-tol is set.
    VIResult value_iteration(int max_iterations, double convergence_tol, int plot_freq, bool conservative_mode,
                             double delta_max) {
        initialize_cells();
        if (successor_cache_empty()) {
            std::printf("  Precomputing successor sets...\n");
            precompute_all_successors();
        }
        std::printf("\nStarting OPTIMIZED PARALLEL value iteration (max %d iterations)...\n", max_iterations);
        std::printf("Conservative mode: %s\n", py_bool(conservative_mode).c_str());
        if (conservative_mode) {
            std::printf("Conservative tolerance delta_max: %g\n", delta_max);
            std::printf("Conservative margin will be: eps_cons = (gamma * delta_actual) / (1 - gamma)\n");
        } else {
            std::printf("Convergence tolerance: %g\n", convergence_tol);
        }
        std::printf("Number of cells: %zu\n", tree_.num_leaves());

        VIResult res = run_sweeps(max_iterations, convergence_tol, conservative_mode, delta_max, plot_freq, false);
        emit_output(output_dir_ + "/value_function_phase_0_complete.png", res.iterations, 0);
        std::printf("\nValue iteration completed in %d iterations\n", res.iterations);
        return res;
    }

    // Local VI after a refinement: reinitialise every leaf, refresh the
    // successor cache, sweep.
    VIResult local_value_iteration(const std::vector<std::uint32_t>& new_cells, int max_iterations,
                                   double convergence_tol, bool conservative_mode, double delta_max) {
        const std::size_t n = tree_.num_leaves();
        std::printf("  Reinitializing ALL %zu cells (%s)\n", n,
                    mode_discounted(mode_) ? "V = l" : "V = min(l, r)");
        const double t0 = now_seconds();
        initialize_cells();
        const double reinit = now_seconds() - t0;
        std::printf("   Reinitialized %zu cells in %.2fs (%.1f cells/s)\n", n, reinit, reinit > 0 ? n / reinit : 0.0);
        std::printf("  Local VI: updating all %zu cells\n", n);
        std::printf("    Conservative mode: %s\n", py_bool(conservative_mode).c_str());
        if (conservative_mode) std::printf("    delta_max: %g\n", delta_max);
        else std::printf("    Convergence tolerance: %g\n", convergence_tol);
        std::printf("    Max iterations: %d\n", max_iterations);
        if (n == 0) { std::printf("    No cells to update!\n"); return VIResult{}; }
        if (!new_cells.empty()) {
            std::printf("    Updating successor cache for %zu updated cells...\n", new_cells.size());
            update_successor_cache(new_cells);
        } else if (successor_cache_empty()) {
            // resume path: nothing split since the checkpoint, but no cache either
            std::printf("    Successor cache empty - computing for all %zu cells...\n", n);
            compute_successors(tree_.leaves());
        }
        VIResult res = run_sweeps(max_iterations, convergence_tol, conservative_mode, delta_max, 0, true);
        std::printf("    Local VI completed in %d iterations\n", res.iterations);
        if (!res.conv_upper.empty()) {
            std::printf("    Final convergence values:\n");
            std::printf("      ||V_upper^final - V_upper^prev||_inf = %.8e\n", res.conv_upper.back());
            std::printf("      ||V_lower^final - V_lower^prev||_inf = %.8e\n", res.conv_lower.back());
        }
        return res;
    }

    std::vector<std::uint32_t> identify_boundary_cells() const {
        std::vector<std::uint32_t> out;
        for (std::uint32_t id : tree_.leaves()) {
            const auto& c = tree_.cell(id);
            if (c.V_upper > 0.0 && c.V_lower <= 0.0) out.push_back(id);
        }
        return out;
    }

    void emit_output(const std::string& png_filename, int iteration, int phase) {
        const double t0 = now_seconds();
        if (!args_.no_plot) plotter_.plot(tree_, png_filename, iteration);
        if (args_.dump_csv) {
            char buf[64];
            std::snprintf(buf, sizeof buf, "/value_function_phase_%d.csv", phase);
            SlicePlotter<N>::dump_csv(tree_, output_dir_ + buf);
        }
        timers_.plotting += now_seconds() - t0;
    }

    // Same rule as the classification panel: safe if V_lower > 0, unsafe if
    // V_upper <= 0, otherwise unclassified. Cell-count percentages treat every
    // leaf equally; volume percentages weight each leaf by the product of its
    // side lengths (the fraction of the state box). One pass over the leaves,
    // negligible next to a value-iteration sweep. Printed to stdout, which the
    // run logger copies into run.log.
    void print_statistics(const char* heading = "Final Cell Classification:") const {
        std::size_t safe = 0, unsafe = 0, unclassified = 0;
        double safe_v = 0.0, unsafe_v = 0.0, unclassified_v = 0.0;
        for (std::uint32_t id : tree_.leaves()) {
            const auto& c = tree_.cell(id);
            double vol = 1.0;
            for (std::size_t d = 0; d < N; ++d) vol *= c.range(d);
            if (c.V_lower > 0.0) { ++safe; safe_v += vol; }
            else if (c.V_upper <= 0.0) { ++unsafe; unsafe_v += vol; }
            else { ++unclassified; unclassified_v += vol; }
        }
        const double n = tree_.num_leaves() ? static_cast<double>(tree_.num_leaves()) : 1.0;
        const double v = (safe_v + unsafe_v + unclassified_v) > 0.0 ? (safe_v + unsafe_v + unclassified_v) : 1.0;
        std::printf("\n%s\n", heading);
        std::printf("  Safe:         %8zu cells (%5.1f%%),  %5.1f%% of state space\n", safe, 100.0 * safe / n, 100.0 * safe_v / v);
        std::printf("  Unsafe:       %8zu cells (%5.1f%%),  %5.1f%% of state space\n", unsafe, 100.0 * unsafe / n, 100.0 * unsafe_v / v);
        std::printf("  Unclassified: %8zu cells (%5.1f%%),  %5.1f%% of state space\n", unclassified, 100.0 * unclassified / n,
                    100.0 * unclassified_v / v);
    }

private:
    // l, r bounds from the centre value and the cell half-width r (inf norm).
    void initialize(const std::vector<std::uint32_t>& ids) {
        if (ids.empty()) return;
        const double t0 = now_seconds();
        const double L_l = dyn_.L_l(), L_r = dyn_.L_r();
        pool_.parallel_for(ids.size(), [&](std::size_t b, std::size_t e) {
            for (std::size_t i = b; i < e; ++i) {
                auto& c = tree_.cell(ids[i]);
                const double l_center = dyn_.failure_function(c.center);
                const double r_center = dyn_.reward_function(c.center);
                double mx = 0.0;
                for (std::size_t d = 0; d < N; ++d) mx = std::max(mx, c.range(d));
                const double r = 0.5 * mx;
                c.l_lower = l_center - L_l * r;
                c.l_upper = l_center + L_l * r;
                c.r_lower = r_center - L_r * r;
                c.r_upper = r_center + L_r * r;
                // Discounted runs decrease from l. Undiscounted runs increase from min(l, r).
                if (mode_discounted(mode_)) {
                    c.V_lower = c.l_lower;
                    c.V_upper = c.l_upper;
                } else {
                    c.V_lower = std::min(c.l_lower, c.r_lower);
                    c.V_upper = std::min(c.l_upper, c.r_upper);
                }
            }
        });
        timers_.cell_init += now_seconds() - t0;
    }

    void ensure_cache_capacity() {
        const std::size_t need = tree_.num_cells() * n_actions_;
        if (cache_.size() < need) { cache_.resize(need); cache_present_.resize(need, 0); }
    }

    void compute_successors(const std::vector<std::uint32_t>& cell_ids) {
        if (cell_ids.empty()) return;
        ensure_cache_capacity();
        const std::size_t n_tasks = cell_ids.size() * n_actions_;
        const double wall0 = now_seconds();
        timers_.ode_solves += n_tasks;

        if (args_.seq_queries) {
            // Python structure: parallel ODE integration, then sequential queries.
            std::vector<std::vector<State>> trajectories(n_tasks);
            pool_.parallel_for(n_tasks, [&](std::size_t b, std::size_t e) {
                for (std::size_t t = b; t < e; ++t) {
                    const std::uint32_t cid = cell_ids[t / n_actions_];
                    dyn_.dynamics_multi_step(tree_.cell(cid).center, t % n_actions_, reach_.tau(), reach_.dt(), trajectories[t]);
                }
            });
            const double ode_wall = now_seconds() - wall0;
            timers_.ode += ode_wall;
            std::printf("   Trajectories computed in %.2fs (%.1f tasks/s)\n", ode_wall, ode_wall > 0 ? n_tasks / ode_wall : 0.0);
            std::printf("  Step 2: Computing successors using spatial index...\n");
            const double q0 = now_seconds();
            std::vector<std::uint32_t> stamp(tree_.num_leaves(), 0);
            std::uint32_t epoch = 0;
            std::vector<std::int32_t> out;
            for (std::size_t t = 0; t < n_tasks; ++t)
                store_successors(cell_ids[t / n_actions_], t % n_actions_, trajectories[t], stamp, epoch, out);
            const double q_wall = now_seconds() - q0;
            timers_.queries += q_wall;
            std::printf("   Spatial queries completed in %.2fs\n", q_wall);
        } else {
            // Trajectory and query fused per task; the index is read-only here.
            const std::size_t nw = pool_.size();
            std::vector<double> ode_cpu(nw, 0.0), query_cpu(nw, 0.0);
            std::vector<std::vector<std::uint32_t>> stamps(nw);
            std::vector<std::uint32_t> epochs(nw, 0);
            for (auto& s : stamps) s.assign(tree_.num_leaves(), 0);
            std::atomic<std::size_t> next_worker{0};
            pool_.parallel_for(n_tasks, [&](std::size_t b, std::size_t e) {
                const std::size_t w = next_worker.fetch_add(1) % nw;
                std::vector<State> traj;
                std::vector<std::int32_t> out;
                for (std::size_t t = b; t < e; ++t) {
                    const std::uint32_t cid = cell_ids[t / n_actions_];
                    const std::size_t ai = t % n_actions_;
                    double t0 = now_seconds();
                    dyn_.dynamics_multi_step(tree_.cell(cid).center, ai, reach_.tau(), reach_.dt(), traj);
                    ode_cpu[w] += now_seconds() - t0;
                    t0 = now_seconds();
                    store_successors(cid, ai, traj, stamps[w], epochs[w], out);
                    query_cpu[w] += now_seconds() - t0;
                }
            });
            const double wall = now_seconds() - wall0;
            double ode_sum = 0.0, q_sum = 0.0;
            for (std::size_t w = 0; w < nw; ++w) { ode_sum += ode_cpu[w]; q_sum += query_cpu[w]; }
            const double denom = (ode_sum + q_sum) > 0 ? (ode_sum + q_sum) : 1.0;
            timers_.ode += wall * ode_sum / denom;
            timers_.queries += wall * q_sum / denom;
            std::printf("   Trajectories + spatial queries computed in %.2fs (%.1f tasks/s)\n", wall, wall > 0 ? n_tasks / wall : 0.0);
            std::printf("    Breakdown: %.1f%% ODE, %.1f%% queries\n", 100.0 * ode_sum / denom, 100.0 * q_sum / denom);
        }
        recount_cache();
    }

    // Leaves whose boxes meet the Gronwall-inflated reach boxes at every
    // checkpoint, plus the OOB sink if a reach box leaves the domain along a
    // non-periodic axis.
    void store_successors(std::uint32_t cell_id, std::size_t action_idx, const std::vector<State>& traj,
                          std::vector<std::uint32_t>& stamp, std::uint32_t& epoch, std::vector<std::int32_t>& out) {
        const auto& c = tree_.cell(cell_id);
        double r;
        if (reach_.use_infinity_norm()) {
            r = 0.5 * c.max_range();
        } else {
            double sq = 0.0;
            for (std::size_t d = 0; d < N; ++d) sq += c.range(d) * c.range(d);
            r = 0.5 * std::sqrt(sq);
        }
        const auto& gb = dyn_.state_bounds();
        out.clear();
        bool oob = false;
        ++epoch;
        for (std::size_t i = 0; i < traj.size(); ++i) {
            const double expansion = r * reach_.growth_factors()[i];
            Bounds<N> rb{};
            for (std::size_t d = 0; d < N; ++d) rb[d] = {{traj[i][d] - expansion, traj[i][d] + expansion}};
            tree_.for_each_intersecting(rb, stamp, epoch, [&](std::uint32_t id) { out.push_back(static_cast<std::int32_t>(id)); });
            for (std::size_t d = 0; d < N; ++d) {
                if (periodic_[d]) continue;
                if (rb[d][0] < gb[d][0] || rb[d][1] > gb[d][1]) { oob = true; break; }
            }
        }
        if (oob) out.push_back(kOOB);
        const std::size_t key = static_cast<std::size_t>(cell_id) * n_actions_ + action_idx;
        cache_[key] = out;
        cache_present_[key] = 1;
    }

    void recount_cache() {
        cache_entries_ = 0;
        for (std::size_t i = 0; i < cache_present_.size(); ++i) if (cache_present_[i]) ++cache_entries_;
    }

    SuccGraph build_graph() const {
        const auto& leaves = tree_.leaves();
        const std::size_t n = leaves.size();
        std::vector<std::int32_t> id_to_pos(tree_.num_cells(), -1);
        for (std::size_t i = 0; i < n; ++i) id_to_pos[leaves[i]] = static_cast<std::int32_t>(i);
        SuccGraph g;
        g.off.resize(n * n_actions_ + 1);
        std::size_t total = 0;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t a = 0; a < n_actions_; ++a) {
                const std::size_t key = static_cast<std::size_t>(leaves[i]) * n_actions_ + a;
                if (key < cache_present_.size() && cache_present_[key]) total += cache_[key].size();
            }
        g.idx.resize(total);
        std::size_t w = 0;
        const std::uint32_t oob_idx = static_cast<std::uint32_t>(n);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t a = 0; a < n_actions_; ++a) {
                g.off[i * n_actions_ + a] = static_cast<std::uint32_t>(w);
                const std::size_t key = static_cast<std::size_t>(leaves[i]) * n_actions_ + a;
                if (key >= cache_present_.size() || !cache_present_[key]) continue;
                for (std::int32_t sid : cache_[key]) {
                    if (sid == kOOB) g.idx[w++] = oob_idx;
                    else {
                        const std::int32_t p = id_to_pos[static_cast<std::size_t>(sid)];
                        if (p >= 0) g.idx[w++] = static_cast<std::uint32_t>(p);
                    }
                }
            }
        g.off[n * n_actions_] = static_cast<std::uint32_t>(w);
        g.idx.resize(w);
        return g;
    }

    VIResult run_sweeps(int max_iterations, double convergence_tol, bool conservative_mode, double delta_max,
                        int plot_freq, bool local) {
        const double t_start = now_seconds();
        const auto& leaves = tree_.leaves();
        const std::size_t n = leaves.size();
        const SuccGraph g = build_graph();

        std::vector<double> Vu(n + 1), Vl(n + 1), lu(n + 1), ll(n + 1), ru(n + 1), rl(n + 1);
        for (std::size_t i = 0; i < n; ++i) {
            const auto& c = tree_.cell(leaves[i]);
            Vu[i] = c.V_upper; Vl[i] = c.V_lower; lu[i] = c.l_upper; ll[i] = c.l_lower; ru[i] = c.r_upper; rl[i] = c.r_lower;
        }
        Vu[n] = Vl[n] = lu[n] = ll[n] = ru[n] = rl[n] = -1.0;   // OOB sink
        std::vector<double> newVu(n + 1), newVl(n + 1);
        newVu[n] = newVl[n] = -1.0;

        const std::size_t nw = pool_.size();
        std::vector<double> part_du(nw), part_dl(nw), part_delta(nw);
        VIResult res;
        const bool discount = mode_discounted(mode_);
        const bool avoid = (mode_ == Mode::AvoidNoDiscount);
        const double gamma = gamma_;

        int iteration = 0;
        bool converged = false;
        double diff_upper = 0.0, diff_lower = 0.0, delta_k = 0.0;
        for (iteration = 0; iteration < max_iterations; ++iteration) {
            std::atomic<std::size_t> chunk{0};
            std::fill(part_du.begin(), part_du.end(), 0.0);
            std::fill(part_dl.begin(), part_dl.end(), 0.0);
            std::fill(part_delta.begin(), part_delta.end(), kInf);

            pool_.parallel_for(n, [&](std::size_t b, std::size_t e) {
                const std::size_t w = chunk.fetch_add(1) % nw;
                double du = 0.0, dl = 0.0, dk = kInf;
                for (std::size_t i = b; i < e; ++i) {
                    double best_min = -kInf, best_max = -kInf;
                    for (std::size_t a = 0; a < n_actions_; ++a) {
                        const std::uint32_t s = g.off[i * n_actions_ + a], t = g.off[i * n_actions_ + a + 1];
                        if (s == t) continue;
                        double lo = kInf, hi = -kInf;
                        for (std::uint32_t k = s; k < t; ++k) {
                            const std::uint32_t j = g.idx[k];
                            lo = std::min(lo, Vl[j]);
                            hi = std::max(hi, Vu[j]);
                        }
                        if (discount) { lo *= gamma; hi *= gamma; }
                        best_min = std::max(best_min, lo);
                        best_max = std::max(best_max, hi);
                    }
                    const double nl = avoid ? std::min(ll[i], best_min) : std::min(ll[i], std::max(rl[i], best_min));
                    const double nu = avoid ? std::min(lu[i], best_max) : std::min(lu[i], std::max(ru[i], best_max));
                    newVl[i] = nl;
                    newVu[i] = nu;
                    du = std::max(du, std::fabs(nu - Vu[i]));
                    dl = std::max(dl, std::fabs(nl - Vl[i]));
                    dk = std::min(dk, nl - Vl[i]);
                }
                part_du[w] = std::max(part_du[w], du);
                part_dl[w] = std::max(part_dl[w], dl);
                part_delta[w] = std::min(part_delta[w], dk);
            });

            diff_upper = 0.0; diff_lower = 0.0; delta_k = kInf;
            for (std::size_t w = 0; w < nw; ++w) {
                diff_upper = std::max(diff_upper, part_du[w]);
                diff_lower = std::max(diff_lower, part_dl[w]);
                delta_k = std::min(delta_k, part_delta[w]);
            }
            std::swap(Vu, newVu);
            std::swap(Vl, newVl);
            res.conv_upper.push_back(diff_upper);
            res.conv_lower.push_back(diff_lower);
            ++timers_.vi_iterations;

            const char* pad = local ? "    Local Iteration" : "Iteration";
            if (conservative_mode) {
                std::printf("%s %3d: ||V_upper^k - V_upper^k-1||_inf = %12.20f, ||V_lower^k - V_lower^k-1||_inf = %12.20f, delta^k = %12.20f\n",
                            pad, iteration + 1, diff_upper, diff_lower, delta_k);
                if (delta_k >= -delta_max && diff_upper < convergence_tol) {
                    apply_conservative_correction(Vl, n, delta_k, convergence_tol, diff_upper);
                    converged = true;
                }
            } else {
                std::printf("%s %3d: ||V_upper^k - V_upper^k-1||_inf = %12.20f, ||V_lower^k - V_lower^k-1||_inf = %12.20f\n",
                            pad, iteration + 1, diff_upper, diff_lower);
                if (diff_upper < convergence_tol && diff_lower < convergence_tol) {
                    std::printf("\n%s\n", std::string(50, '=').c_str());
                    std::printf(" %s CONVERGENCE ACHIEVED\n", local ? "LOCAL" : "STANDARD");
                    std::printf("%s\n", std::string(50, '=').c_str());
                    std::printf("  ||V_upper^k - V_upper^k-1||_inf = %.20e < tolerance = %g\n", diff_upper, convergence_tol);
                    std::printf("  ||V_lower^k - V_lower^k-1||_inf = %.20e < tolerance = %g\n", diff_lower, convergence_tol);
                    std::printf("  Converged in %d iterations\n", iteration + 1);
                    converged = true;
                }
            }
            if (converged) { ++iteration; break; }

            if (plot_freq > 0 && (iteration + 1) % plot_freq == 0) {
                write_back(Vu, Vl, n);
                char buf[128];
                std::snprintf(buf, sizeof buf, "/iteration_%04d_refinement_%02d.png", iteration + 1, refinement_phase_);
                emit_output(output_dir_ + buf, iteration + 1, refinement_phase_);
                std::printf("  [Plot saved at iteration %d]\n", iteration + 1);
            }
        }

        if (!converged) {
            iteration = max_iterations;
            std::printf("\n%s\n", std::string(50, '!').c_str());
            std::printf("  MAXIMUM ITERATIONS REACHED: %d\n", max_iterations);
            std::printf("%s\n", std::string(50, '!').c_str());
            if (conservative_mode) {
                std::printf("  Final delta^k = %.20e\n", delta_k);
                const double eps_cons = (gamma * std::fabs(delta_k)) / (1.0 - gamma);
                std::printf("  Using delta_actual = %.10e for correction\n", std::fabs(delta_k));
                std::printf("  Applying conservative margin: eps_cons = %.10e\n", eps_cons);
                for (std::size_t i = 0; i < n; ++i) Vl[i] -= eps_cons;
                std::printf("   Conservative correction applied to %zu cells\n", n);
            }
            std::printf("  Final ||V_upper^k - V_upper^k-1||_inf = %.20e\n", diff_upper);
            std::printf("  Final ||V_lower^k - V_lower^k-1||_inf = %.20e\n", diff_lower);
        }

        write_back(Vu, Vl, n);
        res.iterations = iteration;
        res.converged = converged;
        timers_.value_iteration += now_seconds() - t_start;
        return res;
    }

    // Algorithm 3: V_lower <- V_lower - eps_cons - eps_machine.
    void apply_conservative_correction(std::vector<double>& Vl, std::size_t n, double delta_k, double convergence_tol,
                                       double diff_upper) {
        std::printf("\n%s\n", std::string(60, '=').c_str());
        std::printf(" CONSERVATIVE STOPPING CONDITION MET (Algorithm 3)\n");
        std::printf("%s\n", std::string(60, '=').c_str());
        std::printf("  delta^k = %.10e >= -delta_max\n", delta_k);
        std::printf("  ||V_upper^k - V_upper^k-1||_inf = %.10e < %g\n", diff_upper, convergence_tol);
        const double delta_to_use = std::fabs(delta_k);
        const double eps_cons = (gamma_ * delta_to_use) / (1.0 - gamma_);
        std::printf("  eps_cons = (%g * %.10e) / %g = %.10e\n", gamma_, delta_to_use, 1.0 - gamma_, eps_cons);
        std::printf("  Correcting V_lower for %zu cells: V_lower <- V_lower - eps_cons - eps_machine\n", n);
        const double eps_machine = DBL_EPSILON;
        std::size_t near_zero = 0, flipped = 0;
        for (std::size_t i = 0; i < n; ++i) {
            const double old = Vl[i];
            Vl[i] = old - eps_cons - eps_machine;
            if (std::fabs(old - eps_cons) < eps_machine) ++near_zero;
            if (old - eps_cons > 0.0 && Vl[i] <= 0.0) ++flipped;
        }
        std::printf("   Conservative correction applied successfully\n");
        std::printf("  Machine epsilon eps_machine = %.3e\n", eps_machine);
        if (near_zero) std::printf("    %zu cells experienced near-cancellation\n", near_zero);
        if (flipped) std::printf("    %zu cells changed from (+) to (<=0) due to eps_machine\n", flipped);
    }

    void write_back(const std::vector<double>& Vu, const std::vector<double>& Vl, std::size_t n) {
        const auto& leaves = tree_.leaves();
        for (std::size_t i = 0; i < n; ++i) {
            auto& c = tree_.cell(leaves[i]);
            c.V_upper = Vu[i];
            c.V_lower = Vl[i];
        }
    }

    Mode mode_;
    const Args& args_;
    const Dynamics<N>& dyn_;
    Tree& tree_;
    const GronwallReachabilityAnalyzer& reach_;
    std::string output_dir_;
    ThreadPool& pool_;
    Timers& timers_;
    double gamma_;
    std::size_t n_actions_;
    std::array<bool, N> periodic_;
    SlicePlotter<N> plotter_;
    int refinement_phase_ = 0;

    std::vector<std::vector<std::int32_t>> cache_;
    std::vector<std::uint8_t> cache_present_;  // not vector<bool>: written concurrently
    std::size_t cache_entries_ = 0;
};

}  // namespace dhj
