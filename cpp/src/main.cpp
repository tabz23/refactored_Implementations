// main.cpp - single binary: `dhj --mode <mode> --dynamics <name> [options]`.
//
// Parses the command line, derives the output directory, starts the run.log
// tee, then hands over to the dynamics selected in dynamics/registry.hpp.
#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <string>

#include "dhj/core/args.hpp"
#include "dhj/core/log.hpp"
#include "dhj/core/run.hpp"
#include "dhj/dynamics/registry.hpp"

int main(int argc, char** argv) {
    using namespace dhj;
    Args args;
    try {
        args = parse_args(argc, argv);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }

    const std::string out_dir = derive_out_dir(args);
    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    if (ec) {
        std::fprintf(stderr, "error: cannot create %s: %s\n", out_dir.c_str(), ec.message().c_str());
        return 1;
    }

    int rc = 0;
    {
        RunLogger logger(out_dir + "/run.log", argc, argv);
        std::printf("DISCRETE HJ REACH-AVOID SOLVER (C++)\n");
        std::printf("Mode: %s (%s.py semantics)\n", args.mode_name.c_str(), mode_script(args.mode));
        std::printf("Algorithm: %d, Workers: %d\n", args.algorithm, args.workers);
        std::printf("Output directory: %s\n", out_dir.c_str());
        std::printf("Conservative mode: %s\n", py_bool(args.conservative).c_str());

        bool found = false;
        try {
            found = with_dynamics(args, [&](const auto& dyn) { rc = run_with_dynamics(args, dyn, out_dir); });
        } catch (const std::exception& e) {
            std::fprintf(stderr, "error: %s\n", e.what());
            rc = 1;
        }
        if (!found) {
            std::fprintf(stderr, "error: unknown --dynamics '%s'; known:", args.dynamics.c_str());
            for (const auto& n : registered_dynamics()) std::fprintf(stderr, " %s", n.c_str());
            std::fprintf(stderr, "\n");
            rc = 2;
        }
    }
    return rc;
}
