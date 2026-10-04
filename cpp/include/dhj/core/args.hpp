// args.hpp - command line, mirroring refactored_Implementations/python/run.py.
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
    int algorithm = 2;
    std::string dynamics = "dubins";
    double velocity = 1.0;   // dubins / evasion
    double a_max = 1.0;      // double integrator
    double v_max = 1.0;      // double integrator
    double dt = 0.1;
    double tau = 1.0;
    double gamma = 0.96;
    // algorithm 1
    int resolution = 10;
    int iterations = -1;     // -1 => per-mode default
    double tolerance = 1e-13;
    int plot_freq = 1000;
    // algorithm 2 / 3
    double epsilon = 0.1;
    int initial_resolution = 15;
    int refinements = 100;
    int vi_iterations = 20000;
    bool conservative = false;
    bool conservative_set = false;   // true if --conservative or --no-conservative was passed
    double delta_max = 1e-6;
    double phase0_tol = 1e-3;   // unused unless --phase0-tol is passed; otherwise phase 0 uses tolerance
    bool phase0_tol_set = false;
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
        "  --algorithm {1,2}             1 = fixed grid VI, 2 = adaptive refinement (default: 2)\n"
        "  --dynamics NAME               dubins | evasion | double_integrator | bicycle |\n"
        "                                van_der_pol | van_der_pol_avoid (default: dubins)\n"
        "  --velocity FLOAT              constant speed for dubins/evasion (default: 1.0)\n"
        "  --a-max FLOAT                 acceleration bound for double_integrator and bicycle (default: 1.0)\n"
        "  --v-max FLOAT                 speed bound for double_integrator and bicycle (default: 1.0)\n"
        "  --dt FLOAT                    checkpoint interval (default: 0.1)\n"
        "  --tau FLOAT                   control duration (default: 1.0)\n"
        "  --gamma FLOAT                 discount factor (default: 0.96)\n"
        "algorithm 1:\n"
        "  --resolution INT              grid resolution (default: 10)\n"
        "  --iterations INT              max VI sweeps (default: 200 / 2000 / 20000 per mode)\n"
        "  --tolerance, --delta-min F    convergence tolerance (default: 1e-13)\n"
        "  --plot-freq INT               plot every N sweeps (default: 1000)\n"
        "algorithm 2/3:\n"
        "  --epsilon FLOAT               refinement error tolerance (default: 0.1)\n"
        "  --initial-resolution INT      initial coarse grid (default: 15)\n"
        "  --refinements INT             max refinement phases (default: 100)\n"
        "  --vi-iterations INT           VI sweeps per phase (default: 20000)\n"
        "  --conservative                conservative stopping (Algorithm 3);\n"
        "                                default on for ra_discount, off otherwise\n"
        "  --no-conservative             force conservative stopping off\n"
        "  --delta-max FLOAT             delta_max for conservative stopping (default: 1e-6)\n"
        "  --phase0-tol FLOAT            phase-0 residual tol (default: same as --delta-min)\n"
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
        "  --dump-csv                    write value_function_phase_N.csv per phase\n"
        "checkpoints (algorithm 2):\n"
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
        else if (f == "--algorithm") { a.algorithm = static_cast<int>(parse_long(need(i), "--algorithm")); ++i; }
        else if (f == "--dynamics") { a.dynamics = need(i); ++i; }
        else if (f == "--velocity") { a.velocity = parse_double(need(i), "--velocity"); ++i; }
        else if (f == "--a-max") { a.a_max = parse_double(need(i), "--a-max"); ++i; }
        else if (f == "--v-max") { a.v_max = parse_double(need(i), "--v-max"); ++i; }
        else if (f == "--dt") { a.dt = parse_double(need(i), "--dt"); ++i; }
        else if (f == "--tau") { a.tau = parse_double(need(i), "--tau"); ++i; }
        else if (f == "--gamma") { a.gamma = parse_double(need(i), "--gamma"); ++i; }
        else if (f == "--resolution") { a.resolution = static_cast<int>(parse_long(need(i), "--resolution")); ++i; }
        else if (f == "--iterations") { a.iterations = static_cast<int>(parse_long(need(i), "--iterations")); ++i; }
        else if (f == "--tolerance" || f == "--delta-min") { a.tolerance = parse_double(need(i), "--tolerance"); ++i; }
        else if (f == "--plot-freq") { a.plot_freq = static_cast<int>(parse_long(need(i), "--plot-freq")); ++i; }
        else if (f == "--epsilon") { a.epsilon = parse_double(need(i), "--epsilon"); ++i; }
        else if (f == "--initial-resolution") { a.initial_resolution = static_cast<int>(parse_long(need(i), "--initial-resolution")); ++i; }
        else if (f == "--refinements") { a.refinements = static_cast<int>(parse_long(need(i), "--refinements")); ++i; }
        else if (f == "--vi-iterations") { a.vi_iterations = static_cast<int>(parse_long(need(i), "--vi-iterations")); ++i; }
        else if (f == "--conservative") {
            if (a.conservative_set && !a.conservative) {
                std::fprintf(stderr, "error: pass only one of --conservative and --no-conservative\n");
                std::exit(2);
            }
            a.conservative = true; a.conservative_set = true;
        }
        else if (f == "--no-conservative") {
            if (a.conservative_set && a.conservative) {
                std::fprintf(stderr, "error: pass only one of --conservative and --no-conservative\n");
                std::exit(2);
            }
            a.conservative = false; a.conservative_set = true;
        }
        else if (f == "--delta-max") { a.delta_max = parse_double(need(i), "--delta-max"); ++i; }
        else if (f == "--phase0-tol") { a.phase0_tol = parse_double(need(i), "--phase0-tol"); a.phase0_tol_set = true; ++i; }
        else if (f == "--workers") { a.workers = static_cast<int>(parse_long(need(i), "--workers")); ++i; }
        else if (f == "--precompute") { a.precompute = true; }
        else if (f == "--seq-queries") { a.seq_queries = true; }
        else if (f == "--results-root") { a.results_root = need(i); ++i; }
        else if (f == "--out-dir") { a.out_dir = need(i); ++i; }
        else if (f == "--tag") { a.tag = need(i); ++i; }
        else if (f == "--no-plot") { a.no_plot = true; }
        else if (f == "--plot-dpi") { a.plot_dpi = static_cast<int>(parse_long(need(i), "--plot-dpi")); ++i; }
        else if (f == "--dump-csv") { a.dump_csv = true; }
        else if (f == "--checkpoint-every") { a.checkpoint_every = static_cast<int>(parse_long(need(i), "--checkpoint-every")); ++i; }
        else if (f == "--keep-checkpoints") { a.keep_checkpoints = static_cast<int>(parse_long(need(i), "--keep-checkpoints")); ++i; }
        else if (f == "--resume") { a.resume = need(i); ++i; }
        else { std::fprintf(stderr, "error: unknown argument %s\n", argv[i]); usage_and_exit(argv[0], 2); }
    }
    if (!have_mode) { std::fprintf(stderr, "error: --mode is required\n"); usage_and_exit(argv[0], 2); }
    if (a.algorithm != 1 && a.algorithm != 2) { std::fprintf(stderr, "error: --algorithm must be 1 or 2\n"); std::exit(2); }
    if (a.iterations < 0) a.iterations = mode_default_iterations(a.mode);
    if (a.workers <= 0) {
        const unsigned hc = std::thread::hardware_concurrency();
        a.workers = static_cast<int>(hc > 1 ? hc - 1 : 1);
    }
    if (a.results_root.empty()) a.results_root = default_results_root(argv[0]);
    if (!a.phase0_tol_set) a.phase0_tol = a.tolerance;
    if (!a.conservative_set) a.conservative = mode_discounted(a.mode);
    return a;
}

// Same scheme as run.py::derive_out_dir, so Python and C++ runs of the same
// configuration land in sibling directories.
inline std::string derive_out_dir(const Args& a) {
    if (!a.out_dir.empty()) return a.out_dir;
    const std::string cons = a.conservative ? "cons" : "std";
    std::string name;
    if (a.algorithm == 1)
        name = "algorithm1_" + a.dynamics + "_res" + std::to_string(a.resolution) + "_dmin" + py_exp(a.tolerance, 0) +
               "_gamma" + py_repr(a.gamma) + "_dt" + py_repr(a.dt) + "_tau" + py_repr(a.tau) + "_" + cons;
    else
        name = a.dynamics + "_res" + std::to_string(a.initial_resolution) + "_eps" + py_repr(a.epsilon) +
               "_dmin" + py_exp(a.tolerance, 0) + "_dmax" + py_exp(a.delta_max, 0) + "_gamma" + py_repr(a.gamma) +
               "_dt" + py_repr(a.dt) + "_tau" + py_repr(a.tau) + "_" + cons;
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
    j += in + "\"algorithm\": " + std::to_string(a.algorithm) + ",\n";
    j += in + "\"checkpoint_every\": " + std::to_string(a.checkpoint_every) + ",\n";
    j += in + "\"conservative\": " + b(a.conservative) + ",\n";
    j += in + "\"delta_max\": " + py_repr(a.delta_max) + ",\n";
    j += in + "\"dt\": " + py_repr(a.dt) + ",\n";
    j += in + "\"dump_csv\": " + b(a.dump_csv) + ",\n";
    j += in + "\"dynamics\": " + s(a.dynamics) + ",\n";
    j += in + "\"epsilon\": " + py_repr(a.epsilon) + ",\n";
    j += in + "\"gamma\": " + py_repr(a.gamma) + ",\n";
    j += in + "\"initial_resolution\": " + std::to_string(a.initial_resolution) + ",\n";
    j += in + "\"iterations\": " + std::to_string(a.iterations) + ",\n";
    j += in + "\"keep_checkpoints\": " + std::to_string(a.keep_checkpoints) + ",\n";
    j += in + "\"mode\": " + s(a.mode_name) + ",\n";
    j += in + "\"no_plot\": " + b(a.no_plot) + ",\n";
    j += in + "\"out_dir\": " + (a.out_dir.empty() ? std::string("null") : s(a.out_dir)) + ",\n";
    j += in + "\"phase0_tol\": " + py_repr(a.phase0_tol) + ",\n";
    j += in + "\"plot_dpi\": " + std::to_string(a.plot_dpi) + ",\n";
    j += in + "\"plot_freq\": " + std::to_string(a.plot_freq) + ",\n";
    j += in + "\"precompute\": " + b(a.precompute) + ",\n";
    j += in + "\"refinements\": " + std::to_string(a.refinements) + ",\n";
    j += in + "\"resolution\": " + std::to_string(a.resolution) + ",\n";
    j += in + "\"results_root\": " + s(a.results_root) + ",\n";
    j += in + "\"resume\": " + (a.resume.empty() ? std::string("null") : s(a.resume)) + ",\n";
    j += in + "\"seq_queries\": " + b(a.seq_queries) + ",\n";
    j += in + "\"tag\": " + s(a.tag) + ",\n";
    j += in + "\"tau\": " + py_repr(a.tau) + ",\n";
    j += in + "\"tolerance\": " + py_repr(a.tolerance) + ",\n";
    j += in + "\"v_max\": " + py_repr(a.v_max) + ",\n";
    j += in + "\"velocity\": " + py_repr(a.velocity) + ",\n";
    j += in + "\"vi_iterations\": " + std::to_string(a.vi_iterations) + ",\n";
    j += in + "\"workers\": " + std::to_string(a.workers) + "\n";
    j += indent + "}";
    return j;
}

}  // namespace dhj
