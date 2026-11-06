// capd_bridge.cpp - validated integration of a single point with CAPD (see validated.hpp).
//
// Compiled only by `make capd`, with the flags of capd-config (CAPD native intervals,
// -frounding-math). Each thread keeps one CAPD solver per vector-field formula, because
// CAPD objects are not thread-safe. The parameters are set through the solver's own
// vector field before every call.
#include "dhj/core/validated.hpp"

#include <capd/capdlib.h>

#include <cfenv>
#include <exception>
#include <memory>
#include <mutex>
#include <unordered_map>

namespace dhj::validated {

namespace {

constexpr int kTaylorOrder = 15;  // fastest order in a test on Dubins steps; enclosures ~2e-15 wide

struct Solver {
    capd::IMap field;
    capd::IOdeSolver solver;
    capd::ITimeMap time_map;
    explicit Solver(const std::string& formula)
        : field(formula), solver(field, kTaylorOrder), time_map(solver) {}
};

// Building a CAPD map from a formula is not thread-safe (concurrent construction corrupted the
// heap in DagIndexer::createIndexArray), so construction is serialized. Each thread builds one
// solver per formula, so the lock is taken only a few times per run.
std::mutex construction_mutex;

Solver& solver_for(const std::string& formula) {
    thread_local std::unordered_map<std::string, std::unique_ptr<Solver>> cache;
    auto it = cache.find(formula);
    if (it == cache.end()) {
        std::lock_guard<std::mutex> lock(construction_mutex);
        it = cache.emplace(formula, std::make_unique<Solver>(formula)).first;
    }
    return *it->second;
}

struct NearestOnExit {
    ~NearestOnExit() { std::fesetround(FE_TONEAREST); }
};

}  // namespace

bool available() { return true; }

bool enclose(const Spec& spec, const double* x0, std::size_t n, const double* t, std::size_t nt, double* lo,
             double* hi, std::string* err) {
    NearestOnExit restore;
    try {
        Solver& s = solver_for(spec.formula);
        auto& field = s.solver.getVectorField();
        for (const auto& p : spec.params) field.setParameter(p.first, capd::interval(p.second));
        capd::IVector start(static_cast<int>(n));
        for (std::size_t d = 0; d < n; ++d) start[static_cast<int>(d)] = capd::interval(x0[d]);
        capd::C0Rect2Set set(start);
        for (std::size_t i = 0; i < nt; ++i) {
            const capd::IVector y = s.time_map(capd::interval(t[i]), set);
            for (std::size_t d = 0; d < n; ++d) {
                lo[i * n + d] = y[static_cast<int>(d)].leftBound();
                hi[i * n + d] = y[static_cast<int>(d)].rightBound();
            }
        }
        return true;
    } catch (const std::exception& e) {
        if (err) *err = e.what();
        return false;
    }
}

}  // namespace dhj::validated
