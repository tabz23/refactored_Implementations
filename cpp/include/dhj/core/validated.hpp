// validated.hpp - optional validated integration with CAPD (experimental).
//
// The default build (`make`) does not use CAPD: available() is false and the solver
// widens each reach box by the integrator's tolerance budget (ode.hpp). `make capd`
// links src/capd_bridge.cpp, which encloses the flow of the cell centre with CAPD's
// rigorous Taylor integrator, and `--validated` then replaces the budget with the
// distance from the numerical centre to the far edge of that enclosure.
//
// The bridge is compiled separately with CAPD's flags (-frounding-math); this header
// keeps CAPD out of the main translation unit, so the numerical integrator is compiled
// exactly as before.
#pragma once

#include <atomic>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace dhj::validated {

// Number of transitions for which CAPD failed in this run. Their reach boxes keep the
// integrator's budget, so they are not validated (run_validated.sh requires 0).
inline std::atomic<std::size_t> failures{0};
// Number of CAPD calls that needed the fixed-step retry (see capd_bridge.cpp).
inline std::atomic<std::size_t> retries{0};

// A vector field in CAPD's formula syntax, e.g.
//   "par:v,u;var:x,y,th;fun:v*cos(th),v*sin(th),u;"
// with the values of its parameters.
struct Spec {
    std::string formula;
    std::vector<std::pair<std::string, double>> params;
};

#ifdef DHJ_WITH_CAPD
// True when this binary was linked with the CAPD bridge.
bool available();
// Rigorous enclosures [lo, hi] of the solution of x' = f(x) from the point x0 (dimension n)
// at the times t[0] < ... < t[nt-1], each > 0. lo and hi receive n * nt values (checkpoint
// major). Returns false, with a message in *err, if CAPD fails. Leaves the calling thread
// in round-to-nearest mode.
bool enclose(const Spec& spec, const double* x0, std::size_t n, const double* t, std::size_t nt, double* lo,
             double* hi, std::string* err = nullptr);
#else
inline bool available() { return false; }
inline bool enclose(const Spec&, const double*, std::size_t, const double*, std::size_t, double*, double*,
                    std::string* err = nullptr) {
    if (err) *err = "built without CAPD (use `make capd`)";
    return false;
}
#endif

}  // namespace dhj::validated
