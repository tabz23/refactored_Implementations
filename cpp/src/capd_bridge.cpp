// capd_bridge.cpp - validated integration of a single point with CAPD (see validated.hpp).
//
// Compiled only by `make capd`, with the flags of capd-config (CAPD native intervals,
// -frounding-math). Each thread keeps its own CAPD solvers per vector-field formula, because
// CAPD objects are not thread-safe. The parameters are set through the solver's own vector
// field before every call.
//
// The first attempt uses CAPD's automatic step control. In some degenerate cases it stops with
// "minimal time step reached" although the solution is simple (the bicycle with steering angle
// exactly 0, steering rate 0 and nonzero acceleration, whose solution is polynomial in time).
// Such a call is retried with fixed steps of t/1, t/2, t/4, t/8 and t/16, which retries()
// counts. Both attempts give rigorous enclosures. If every attempt fails, enclose() returns
// false.
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
    capd::IMap field, fixed_field;
    capd::IOdeSolver solver, fixed_solver;
    capd::ITimeMap time_map, fixed_time_map;
    explicit Solver(const std::string& formula)
        : field(formula), fixed_field(formula), solver(field, kTaylorOrder), fixed_solver(fixed_field, kTaylorOrder),
          time_map(solver), fixed_time_map(fixed_solver) {
        fixed_solver.turnOffStepControl();
    }
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

// Encloses the solution from x0 at the times t[0..nt-1] with the given solver; throws on failure.
void integrate(const Spec& spec, capd::IOdeSolver& solver, capd::ITimeMap& time_map, const double* x0,
               std::size_t n, const double* t, std::size_t nt, double* lo, double* hi) {
    auto& field = solver.getVectorField();
    for (const auto& p : spec.params) field.setParameter(p.first, capd::interval(p.second));
    capd::IVector start(static_cast<int>(n));
    for (std::size_t d = 0; d < n; ++d) start[static_cast<int>(d)] = capd::interval(x0[d]);
    capd::C0Rect2Set set(start);
    for (std::size_t i = 0; i < nt; ++i) {
        const capd::IVector y = time_map(capd::interval(t[i]), set);
        for (std::size_t d = 0; d < n; ++d) {
            lo[i * n + d] = y[static_cast<int>(d)].leftBound();
            hi[i * n + d] = y[static_cast<int>(d)].rightBound();
        }
    }
}

}  // namespace

bool available() { return true; }

bool enclose(const Spec& spec, const double* x0, std::size_t n, const double* t, std::size_t nt, double* lo,
             double* hi, std::string* err) {
    NearestOnExit restore;
    Solver* s = nullptr;
    try {
        s = &solver_for(spec.formula);
        integrate(spec, s->solver, s->time_map, x0, n, t, nt, lo, hi);
        return true;
    } catch (const std::exception& e) {
        if (err) *err = e.what();
        if (s == nullptr) return false;
    }
    for (int steps : {1, 2, 4, 8, 16}) {
        try {
            s->fixed_solver.setStep(capd::interval(t[nt - 1]) / steps);
            integrate(spec, s->fixed_solver, s->fixed_time_map, x0, n, t, nt, lo, hi);
            retries.fetch_add(1, std::memory_order_relaxed);
            return true;
        } catch (const std::exception& e) {
            if (err) *err = e.what();
        }
    }
    return false;
}

}  // namespace dhj::validated
