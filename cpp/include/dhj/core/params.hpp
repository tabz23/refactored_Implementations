// params.hpp - params.json and phases.csv, the machine-readable run record.
//
// params.json  all CLI arguments, the dynamics description, derived quantities
//              (eta_min, growth factors, workers) and the environment.
// phases.csv   one row per completed refinement phase: |cells|, boundary
//              cells, VI sweeps, timings. Same columns as the Python package.
#pragma once

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <unistd.h>

#include "args.hpp"
#include "timing.hpp"
#include "types.hpp"

namespace dhj {

inline std::string json_escape(const std::string& s) {
    std::string o;
    for (char c : s) {
        if (c == '"' || c == '\\') { o += '\\'; o += c; }
        else if (c == '\n') o += "\\n";
        else o += c;
    }
    return o;
}

// Writes <out_dir>/params.json; a previous copy is renamed params_<stamp>.json
// so a resumed run never loses its history.
inline std::string write_params(const std::string& out_dir, const Args& args, const std::string& dynamics_json,
                                const std::string& derived_json) {
    namespace fs = std::filesystem;
    const std::string path = out_dir + "/params.json";
    std::error_code ec;
    if (fs::exists(path, ec)) {
        const auto t = fs::last_write_time(path, ec);
        (void)t;
        fs::rename(path, out_dir + "/params_" + timestamp_compact() + ".json", ec);
    }
    std::string cmd;
    for (std::size_t i = 0; i < args.raw.size(); ++i) { if (i) cmd += ' '; cmd += args.raw[i]; }
    char host[256] = "?";
    ::gethostname(host, sizeof host);

    std::ofstream f(path);
    f << "{\n";
    f << "  \"command\": \"" << json_escape(cmd) << "\",\n";
    f << "  \"started\": \"" << timestamp_now() << "\",\n";
    f << "  \"args\": " << args_json(args, "  ") << ",\n";
    f << "  \"dynamics\": " << dynamics_json << ",\n";
    f << "  \"derived\": " << derived_json << ",\n";
    f << "  \"environment\": {\n";
    f << "    \"host\": \"" << json_escape(host) << "\",\n";
    f << "    \"implementation\": \"cpp\",\n";
#if defined(__clang__)
    f << "    \"compiler\": \"clang " << __clang_major__ << "." << __clang_minor__ << "\",\n";
#elif defined(__GNUC__)
    f << "    \"compiler\": \"gcc " << __GNUC__ << "." << __GNUC_MINOR__ << "\",\n";
#else
    f << "    \"compiler\": \"unknown\",\n";
#endif
    f << "    \"build\": \"" << __DATE__ << " " << __TIME__ << "\"\n";
    f << "  }\n}\n";
    return path;
}

inline const char* kPhaseColumns =
    "phase,n_leaves,boundary_cells,refinable_cells,refined_parents,vi_iterations,converged,"
    "time_vi_s,time_successors_s,time_phase_s,wall_since_start_s,timestamp";

// Appends one CSV row per phase to <out_dir>/phases.csv.
class PhaseLog {
public:
    explicit PhaseLog(const std::string& out_dir) : path_(out_dir + "/phases.csv") {
        std::error_code ec;
        if (!std::filesystem::exists(path_, ec)) {
            std::ofstream f(path_);
            f << kPhaseColumns << "\n";
        }
    }

    void append(int phase, std::size_t n_leaves, std::size_t boundary, std::size_t refinable,
                std::size_t refined_parents, int vi_iterations, bool converged, double time_vi,
                double time_succ, double time_phase, double wall_since_start) {
        std::ofstream f(path_, std::ios::app);
        f << phase << ',' << n_leaves << ',' << boundary << ',' << refinable << ',' << refined_parents << ','
          << vi_iterations << ',' << (converged ? "True" : "False") << ',' << py_fixed(time_vi, 3) << ','
          << py_fixed(time_succ, 3) << ',' << py_fixed(time_phase, 3) << ',' << py_fixed(wall_since_start, 3)
          << ',' << timestamp_now() << "\n";
    }

private:
    std::string path_;
};

}  // namespace dhj
