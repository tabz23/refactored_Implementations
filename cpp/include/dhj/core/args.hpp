// args.hpp - command line of the dhj solver.
#pragma once

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "types.hpp"

namespace dhj {

struct Args {
    // problem
    Mode mode = Mode::RANoDiscount;
    std::string mode_name = "ra_nodiscount";
    std::string dynamics = "dubins";
    double velocity = 1.0;   // dubins / evasion
    double a_max = 1.0;      // double integrator
    double v_max = 1.0;      // double integrator
    double dt = 0.1;
    double tau = 0.1;
    bool tau_set = false;    // false => tau = dt (one reach checkpoint per step)
    double gamma = 0.96;     // ra_discount only; the other modes use gamma = 1
    bool gamma_set = false;
    // value iteration and Algorithm 2 (Algorithm 1 of the paper is Algorithm 2 with --refinements 0)
    double delta_lower = 1e-13;   // delta_lower*: stop once the lower bound changes by less than this ...
    double delta_upper = 1e-13;   // ... and the upper bound by less than delta_upper* (sup norm, one sweep)
    int plot_freq = 1000;
    double eps_min = 0.01;   // threshold of Algorithm 2 on the cell radius
    int initial_resolution = 15;
    int refinements = 100;
    int vi_iterations = 100000;
    // Upper bound uses the control-cell tube; lower bound keeps the discrete action.
    bool two_abstractions = false;
    // Drop the part of each reach box outside the state box instead of adding the
    // out-of-bounds sink, as if f mapped X into X. For comparison only: not sound for
    // systems that can leave X.
    bool restrict_to_grid = false;
    // With restrict_to_grid: an action whose reach box lies entirely outside the state box
    // (no successor) contributes only the cell's own bounds, as if the trajectory ended
    // there: min{l, r} in reach-avoid (the same as skipping the action) and l in the avoid
    // case. Without this flag such an action is skipped (counts as -inf).
    bool if_leave_use_cell_bounds = false;
    // Experimental (binary built with `make capd`): replace the integrator's tolerance budget
    // in each reach box with a CAPD-validated bound on the integration error of the centre.
    bool validated = false;
    // execution
    int workers = 0;         // 0 => hardware threads - 1
    bool precompute = false;
    bool seq_queries = false;
    // output
    std::string results_root;  // empty => <refactored_Implementations>/results
    std::string out_dir;
    std::string tag;
    bool no_plot = false;
    int plot_dpi = 800;    // 800: vector SVG (Python: 800 dpi). Smaller: raster PNG
    int plot_slice = -1;   // -1: every slice the dynamics declares. Otherwise one index.
    bool dump_csv = false;
    // checkpoints
    int checkpoint_every = 1;
    int keep_checkpoints = 2;
    std::string resume;

    std::vector<std::string> raw;  // argv, for params.json / run.log
};

[[noreturn]] inline void usage_and_exit(const char* prog, int code) {
    std::fprintf(code == 0 ? stdout : stderr,
        "usage: %s --mode MODE [--dynamics NAME] [options]\n\n"
        "Discrete Hamilton-Jacobi reach-avoid / avoid solver (C++).\n\n"
        "problem:\n"
        "  --mode {ra_nodiscount,ra_discount,avoid_nodiscount}   (required)\n"
        "  --dynamics NAME               dubins | evasion | double_integrator | bicycle |\n"
        "                                van_der_pol (default: dubins)\n"
        "  --velocity FLOAT              constant speed for dubins/evasion (default: 1.0)\n"
        "  --a-max FLOAT                 acceleration bound for double_integrator and bicycle (default: 1.0)\n"
        "  --v-max FLOAT                 speed bound for double_integrator and bicycle (default: 1.0)\n"
        "  --dt FLOAT                    sampling time t_s (default: 0.1)\n"
        "  --tau FLOAT                   control duration, if it differs from --dt; --dt is then\n"
        "                                the spacing of the reach checkpoints (default: --dt)\n"
        "  --gamma FLOAT                 discount factor in (0, 1), for --mode ra_discount only\n"
        "                                (default: 0.96); the other modes use gamma = 1\n"
        "value iteration and Algorithm 2 (Algorithm 1 is Algorithm 2 with --refinements 0):\n"
        "  --initial-resolution INT      initial uniform grid, cells per dimension (default: 15)\n"
        "  --refinements INT             max refinement iterations of Algorithm 2 (default: 100)\n"
        "  --eps-min FLOAT               Algorithm 2 refines a cell only while its radius (half its\n"
        "                                largest side) exceeds eps_min (default: 0.01)\n"
        "  --delta-lower FLOAT           stopping tolerance of the lower bound: value iteration\n"
        "                                stops once a sweep changes V_lower by less than this in the\n"
        "                                sup norm (default: 1e-13) ...\n"
        "  --delta-upper FLOAT           ... and V_upper by less than this (default: 1e-13)\n"
        "  --vi-iterations INT           max VI sweeps per phase (default: 100000)\n"
        "  --plot-freq INT               also plot every N sweeps of a phase (default: 1000)\n"
        "  --two-abstractions            V_upper successors cover each control cell;\n"
        "                                V_lower successors keep the discrete action\n"
        "  --restrict-reach-set-to-grid  ignore the part of each reach box outside the state\n"
        "                                box (no sink state), as if f mapped X into X; for\n"
        "                                comparison only, not sound if the system can leave X\n"
        "  --if-leave-use-cellwise-bounds  with --restrict-reach-set-to-grid: an action with no\n"
        "                                successor uses only the cell's own bounds (as if the\n"
        "                                trajectory ended) instead of counting as -inf\n"
        "execution:\n"
        "  --workers INT                 worker threads (default: hardware threads - 1)\n"
        "  --precompute                  precompute successor sets up front\n"
        "  --seq-queries                 single-threaded spatial queries (Python structure)\n"
        "output:\n"
        "  --results-root PATH           (default: <refactored_Implementations>/results)\n"
        "  --out-dir PATH                override the derived output directory\n"
        "  --tag STR                     suffix appended to the derived directory name\n"
        "  --no-plot                     skip figures\n"
        "  --plot-dpi INT                800 = vector SVG / 800 dpi (default). A smaller\n"
        "                                value writes a raster PNG at that dpi\n"
        "  --plot-slice INT              draw only this slice (default: all; 0 is θ=0 for Dubins)\n"
        "  --dump-csv                    write value_function_phase_N.csv per phase\n"
        "checkpoints:\n"
        "  --checkpoint-every N          checkpoint every N phases (default: 1; 0 = never)\n"
        "  --keep-checkpoints K          keep only the newest K (default: 2; 0 = keep all)\n"
        "  --resume PATH|auto            resume from a checkpoint ('auto' = latest in <out-dir>/checkpoints)\n",
        prog);
    std::exit(code);
}

inline double parse_double(const char* s, const char* flag) {
    char* end = nullptr;
    const double v = std::strtod(s, &end);
    if (end == s || *end != '\0') throw std::invalid_argument(std::string("bad float for ") + flag);
    return v;
}

inline long parse_long(const char* s, const char* flag) {
    char* end = nullptr;
    const long v = std::strtol(s, &end, 10);
    if (end == s || *end != '\0') throw std::invalid_argument(std::string("bad int for ") + flag);
    return v;
}

// Prefer <refactored_Implementations>/results when the binary lives in cpp/build/.
inline std::string default_results_root(const char* argv0) {
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::path p = fs::absolute(argv0 ? argv0 : "", ec);
    if (!ec) {
        p = fs::weakly_canonical(p, ec);
        const fs::path parent = p.parent_path();
        if (parent.filename() == "build" && parent.parent_path().filename() == "cpp")
            return (parent.parent_path().parent_path() / "results").string();
    }
    return "results";
}

inline Args parse_args(int argc, char** argv) {
    Args a;
    for (int i = 0; i < argc; ++i) a.raw.emplace_back(argv[i]);
    bool have_mode = false;
    auto need = [&](int i) -> const char* {
        if (i + 1 >= argc) { std::fprintf(stderr, "error: %s needs a value\n", argv[i]); std::exit(2); }
        return argv[i + 1];
    };
    for (int i = 1; i < argc; ++i) {
        const std::string f = argv[i];
        if (f == "-h" || f == "--help") usage_and_exit(argv[0], 0);
        else if (f == "--mode") {
            a.mode_name = need(i); ++i;
            if (!parse_mode(a.mode_name, a.mode)) { std::fprintf(stderr, "error: unknown --mode %s\n", a.mode_name.c_str()); std::exit(2); }
            have_mode = true;
        }
        else if (f == "--dynamics") { a.dynamics = need(i); ++i; }
        else if (f == "--velocity") { a.velocity = parse_double(need(i), "--velocity"); ++i; }
        else if (f == "--a-max") { a.a_max = parse_double(need(i), "--a-max"); ++i; }
        else if (f == "--v-max") { a.v_max = parse_double(need(i), "--v-max"); ++i; }
        else if (f == "--dt") { a.dt = parse_double(need(i), "--dt"); ++i; }
        else if (f == "--tau") { a.tau = parse_double(need(i), "--tau"); a.tau_set = true; ++i; }
        else if (f == "--gamma") { a.gamma = parse_double(need(i), "--gamma"); a.gamma_set = true; ++i; }
        else if (f == "--delta-lower") { a.delta_lower = parse_double(need(i), "--delta-lower"); ++i; }
        else if (f == "--delta-upper") { a.delta_upper = parse_double(need(i), "--delta-upper"); ++i; }
        else if (f == "--plot-freq") { a.plot_freq = static_cast<int>(parse_long(need(i), "--plot-freq")); ++i; }
        else if (f == "--eps-min") { a.eps_min = parse_double(need(i), "--eps-min"); ++i; }
        else if (f == "--initial-resolution") { a.initial_resolution = static_cast<int>(parse_long(need(i), "--initial-resolution")); ++i; }
        else if (f == "--refinements") { a.refinements = static_cast<int>(parse_long(need(i), "--refinements")); ++i; }
        else if (f == "--vi-iterations") { a.vi_iterations = static_cast<int>(parse_long(need(i), "--vi-iterations")); ++i; }
        else if (f == "--two-abstractions") { a.two_abstractions = true; }
        else if (f == "--restrict-reach-set-to-grid") { a.restrict_to_grid = true; }
        else if (f == "--if-leave-use-cellwise-bounds") { a.if_leave_use_cell_bounds = true; }
        else if (f == "--validated") { a.validated = true; }
        else if (f == "--workers") { a.workers = static_cast<int>(parse_long(need(i), "--workers")); ++i; }
        else if (f == "--precompute") { a.precompute = true; }
        else if (f == "--seq-queries") { a.seq_queries = true; }
        else if (f == "--results-root") { a.results_root = need(i); ++i; }
        else if (f == "--out-dir") { a.out_dir = need(i); ++i; }
        else if (f == "--tag") { a.tag = need(i); ++i; }
        else if (f == "--no-plot") { a.no_plot = true; }
        else if (f == "--plot-dpi") { a.plot_dpi = static_cast<int>(parse_long(need(i), "--plot-dpi")); ++i; }
        else if (f == "--plot-slice") { a.plot_slice = static_cast<int>(parse_long(need(i), "--plot-slice")); ++i; }
        else if (f == "--dump-csv") { a.dump_csv = true; }
        else if (f == "--checkpoint-every") { a.checkpoint_every = static_cast<int>(parse_long(need(i), "--checkpoint-every")); ++i; }
        else if (f == "--keep-checkpoints") { a.keep_checkpoints = static_cast<int>(parse_long(need(i), "--keep-checkpoints")); ++i; }
        else if (f == "--resume") { a.resume = need(i); ++i; }
        else { std::fprintf(stderr, "error: unknown argument %s\n", argv[i]); usage_and_exit(argv[0], 2); }
    }
    if (!have_mode) { std::fprintf(stderr, "error: --mode is required\n"); usage_and_exit(argv[0], 2); }
    if (!(a.eps_min > 0.0)) { std::fprintf(stderr, "error: --eps-min must be positive\n"); std::exit(2); }
    if (mode_discounted(a.mode)) {
        if (!(a.gamma > 0.0 && a.gamma < 1.0)) { std::fprintf(stderr, "error: --gamma must be in (0, 1)\n"); std::exit(2); }
    } else {
        if (a.gamma_set) {
            std::fprintf(stderr, "error: --gamma applies only to --mode ra_discount; %s uses gamma = 1\n",
                         a.mode_name.c_str());
            std::exit(2);
        }
        a.gamma = 1.0;
    }
    if (a.if_leave_use_cell_bounds && !a.restrict_to_grid) {
        std::fprintf(stderr, "error: --if-leave-use-cellwise-bounds needs --restrict-reach-set-to-grid "
                             "(otherwise every action has a successor)\n");
        std::exit(2);
    }
    if (!a.tau_set) a.tau = a.dt;
    if (a.workers <= 0) {
        const unsigned hc = std::thread::hardware_concurrency();
        a.workers = static_cast<int>(hc > 1 ? hc - 1 : 1);
    }
    if (a.results_root.empty()) a.results_root = default_results_root(argv[0]);
    return a;
}

// <results_root>/<mode>/<dynamics>_res<R>[_epsmin<eps_min>]_dlower<delta_lower>_dupper<delta_upper>
// [_gamma<gamma>]_dt<dt>[_tau<tau>][_twoabs][_<tag>]: _epsmin only if the run refines, _gamma only for
// ra_discount, _tau only if --tau differs from --dt. Runs made before 2026-11-04 were named
// <dynamics>_res<R>_eps<...>_dmin<delta>_dmax<...>_gamma<gamma>_dt<dt>_tau<tau>_std.
inline std::string derive_out_dir(const Args& a) {
    if (!a.out_dir.empty()) return a.out_dir;
    std::string name = a.dynamics + "_res" + std::to_string(a.initial_resolution);
    if (a.refinements > 0) name += "_epsmin" + py_repr(a.eps_min);
    name += "_dlower" + py_exp(a.delta_lower, 0) + "_dupper" + py_exp(a.delta_upper, 0);
    if (mode_discounted(a.mode)) name += "_gamma" + py_repr(a.gamma);
    name += "_dt" + py_repr(a.dt);
    if (a.tau != a.dt) name += "_tau" + py_repr(a.tau);
    if (a.two_abstractions) name += "_twoabs";
    if (a.restrict_to_grid) name += "_restrictgrid";
    if (a.if_leave_use_cell_bounds) name += "_cellbounds";
    if (a.validated) name += "_validated";
    if (!a.tag.empty()) name += "_" + a.tag;
    return a.results_root + "/" + a.mode_name + "/" + name;
}

// JSON object with every argument (sorted by name, like the Python's vars(args)).
inline std::string args_json(const Args& a, const std::string& indent) {
    const std::string in = indent + "  ";
    auto s = [](const std::string& v) { return "\"" + v + "\""; };
    auto b = [](bool v) { return v ? "true" : "false"; };
    std::string j = "{\n";
    j += in + "\"a_max\": " + py_repr(a.a_max) + ",\n";
    j += in + "\"checkpoint_every\": " + std::to_string(a.checkpoint_every) + ",\n";
    j += in + "\"delta_lower\": " + py_repr(a.delta_lower) + ",\n";
    j += in + "\"delta_upper\": " + py_repr(a.delta_upper) + ",\n";
    j += in + "\"dt\": " + py_repr(a.dt) + ",\n";
    j += in + "\"dump_csv\": " + b(a.dump_csv) + ",\n";
    j += in + "\"dynamics\": " + s(a.dynamics) + ",\n";
    j += in + "\"eps_min\": " + py_repr(a.eps_min) + ",\n";
    j += in + "\"gamma\": " + py_repr(a.gamma) + ",\n";
    j += in + "\"if_leave_use_cellwise_bounds\": " + b(a.if_leave_use_cell_bounds) + ",\n";
    j += in + "\"initial_resolution\": " + std::to_string(a.initial_resolution) + ",\n";
    j += in + "\"keep_checkpoints\": " + std::to_string(a.keep_checkpoints) + ",\n";
    j += in + "\"mode\": " + s(a.mode_name) + ",\n";
    j += in + "\"no_plot\": " + b(a.no_plot) + ",\n";
    j += in + "\"out_dir\": " + (a.out_dir.empty() ? std::string("null") : s(a.out_dir)) + ",\n";
    j += in + "\"plot_dpi\": " + std::to_string(a.plot_dpi) + ",\n";
    j += in + "\"plot_slice\": " + std::to_string(a.plot_slice) + ",\n";
    j += in + "\"plot_freq\": " + std::to_string(a.plot_freq) + ",\n";
    j += in + "\"precompute\": " + b(a.precompute) + ",\n";
    j += in + "\"refinements\": " + std::to_string(a.refinements) + ",\n";
    j += in + "\"restrict_reach_set_to_grid\": " + b(a.restrict_to_grid) + ",\n";
    j += in + "\"results_root\": " + s(a.results_root) + ",\n";
    j += in + "\"resume\": " + (a.resume.empty() ? std::string("null") : s(a.resume)) + ",\n";
    j += in + "\"seq_queries\": " + b(a.seq_queries) + ",\n";
    j += in + "\"tag\": " + s(a.tag) + ",\n";
    j += in + "\"tau\": " + py_repr(a.tau) + ",\n";
    j += in + "\"two_abstractions\": " + b(a.two_abstractions) + ",\n";
    j += in + "\"validated\": " + b(a.validated) + ",\n";
    j += in + "\"v_max\": " + py_repr(a.v_max) + ",\n";
    j += in + "\"velocity\": " + py_repr(a.velocity) + ",\n";
    j += in + "\"vi_iterations\": " + std::to_string(a.vi_iterations) + ",\n";
    j += in + "\"workers\": " + std::to_string(a.workers) + "\n";
    j += indent + "}";
    return j;
}

}  // namespace dhj
