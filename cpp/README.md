# C++ implementation

Header-only C++17 port of the Python package, same algorithm, same CLI, same
output files; one binary for all modes and all dynamics. No dependencies
beyond a C++17 compiler and pthreads.

```
cpp/
├── Makefile / CMakeLists.txt
├── src/main.cpp                  parse args -> out dir -> run.log tee -> dispatch on --dynamics
├── tools/compare_checkpoints.cpp cell-by-cell comparison of two checkpoints
└── include/dhj/
    ├── core/
    │   ├── types.hpp      State<N>, Bounds<N>, Mode, Python-style number formatting
    │   ├── ode.hpp        Dormand-Prince 5(4), atol = rtol = 1e-12 (stands in for odeint)
    │   ├── rtree.hpp      BoxIndex<N>: bulk-built BVH (stands in for libspatialindex)
    │   ├── cell.hpp       Cell<N>, CellTree<N> (periodic wrap-around queries, bisection)
    │   ├── reach.hpp      Gronwall growth factors
    │   ├── vi.hpp         SafetyValueIterator<N>: successor cache, parallel Bellman sweeps
    │   ├── refine.hpp     AdaptiveRefinement<N>: phase loop, checkpoints, phases.csv
    │   ├── checkpoint.hpp binary checkpoints, --resume auto
    │   ├── svg.hpp        SlicePlotter<N> (SVG figure) + CSV dump
    │   ├── params.hpp     params.json, phases.csv
    │   ├── log.hpp        RunLogger: tee stdout/stderr -> run.log
    │   ├── args.hpp       CLI, output-directory naming (same scheme as run.py)
    │   ├── pool.hpp       thread pool
    │   └── run.hpp        algorithm 1 / 2 drivers
    └── dynamics/
        ├── base.hpp              Dynamics<N> interface + Overlay/Slice + dynamics_multi_step
        ├── dubins.hpp            Dubins car (3-D, theta periodic)
        ├── evasion.hpp           pursuit/evasion (3-D, theta periodic)
        ├── double_integrator.hpp planar double integrator (4-D)
        └── registry.hpp          --dynamics name -> class
```

## Build and run

```bash
make                      # -> build/dhj          (or: cmake -B build && cmake --build build)
make tools                # -> build/dhj_compare_checkpoints
./build/dhj --help
```

```bash
./build/dhj --mode ra_nodiscount --dynamics dubins \
    --dt 0.3 --tau 0.3 --gamma 0.9 --delta-min 1e-8 --delta-max 1e-8 \
    --epsilon 0.05 --initial-resolution 20 --vi-iterations 100033 --conservative
```

Flags are identical to `python/run.py` (`--mode`, `--algorithm`, `--dynamics`,
`--velocity`, `--a-max`, `--v-max`, `--dt`, `--tau`, `--gamma`, `--resolution`,
`--iterations`, `--tolerance|--delta-min`, `--plot-freq`, `--epsilon`,
`--initial-resolution`, `--refinements`, `--vi-iterations`, `--conservative`,
`--delta-max`, `--phase0-tol`, `--workers`, `--precompute`, `--results-root`,
`--out-dir`, `--tag`, `--no-plot`, `--dump-csv`, `--checkpoint-every`,
`--keep-checkpoints`, `--resume`). Extra: `--seq-queries` (single-threaded
spatial queries, the Python's structure; same results, slower).

Output directory and files are as described in the top-level README. Figures are
written as **SVG** (`value_function_phase_N_complete.svg`), one column per slice
declared by the dynamics, three rows (upper bound, lower bound, classification),
with colour bars and overlays. For publication-quality matplotlib figures
identical to the Python ones, add `--dump-csv` and run
`python ../python/plot_from_csv.py --dynamics <name> --mode <mode> value_function_phase_N.csv`.

### Checkpoints / resume

`checkpoints/checkpoint_phase_NNNN.bin` is written after phase 0 and after every
phase (`--checkpoint-every`, `--keep-checkpoints`). Resume with
`--resume auto` or `--resume PATH` plus the same mode/dynamics/parameters. The
format is the one used by the earlier `ComputingHJ/cpp` port, so those
checkpoints can also be resumed and compared.

```bash
./build/dhj_compare_checkpoints A/checkpoints/checkpoint_phase_0005.bin B/checkpoints/checkpoint_phase_0005.bin
```

prints the maximum differences in bounds / V_lower / V_upper / l / r and the
number of classification mismatches (exit 0 when within `--tol`, default 1e-12).

## Adding a new dynamics

Create `include/dhj/dynamics/my_system.hpp`:

```cpp
#pragma once
#include <cmath>
#include "base.hpp"

namespace dhj {

class MySystem : public Dynamics<3> {            // <N> = state dimension
public:
    MySystem(double dt, double tau, double speed = 1.0)
        : Dynamics<3>(dt, tau), speed_(speed) {
        bounds_ = {{{{-3.0, 3.0}}, {{-3.0, 3.0}}, {{-kPi, kPi}}}};   // periodic dims MUST be [-pi, pi]
        actions_ = {-1.0, 0.0, 1.0};
        L_f_ = speed;            // ||f(x,u) - f(y,u)||_inf <= L_f ||x - y||_inf
        L_l_ = std::sqrt(2.0);   // |l(x) - l(y)| <= L_l ||x - y||_inf
        L_r_ = std::sqrt(2.0);
    }
    const char* name() const override { return "my_system"; }
    const Bounds& state_bounds() const override { return bounds_; }
    std::size_t num_actions() const override { return actions_.size(); }
    std::vector<double> action(std::size_t i) const override { return {actions_[i]}; }  // for params.json
    std::vector<std::size_t> periodic_dims() const override { return {2}; }
    std::string state_name(std::size_t d) const override { return d == 0 ? "x" : d == 1 ? "y" : "θ"; }

    State ode(const State& s, std::size_t a) const override {        // f(x, u_a)
        return {speed_ * std::cos(s[2]), speed_ * std::sin(s[2]), actions_[a]};
    }
    double failure_function(const State& s) const override { return std::hypot(s[0], s[1]) - 1.3; }
    double reward_function(const State& s) const override { return -(std::hypot(s[0] - 2.5, s[1]) - 0.5); }

    // ---- plotting ----------------------------------------------------------
    std::vector<Overlay> overlays() const override {
        Overlay obs; obs.cx = 0; obs.cy = 0; obs.radius = 1.3; obs.color = "darkblue";
        Overlay tgt; tgt.cx = 2.5; tgt.cy = 0; tgt.radius = 0.5; tgt.color = "orange"; tgt.dashed = true;
        tgt.reach_avoid_only = true;                 // hidden in avoid mode
        return {obs, tgt};
    }
    std::vector<Slice> slices() const override {     // one figure column per slice
        std::vector<Slice> out;
        for (double th : {0.0, kPi, kPi / 2, -kPi / 2})
            out.push_back(Slice{{{2, th}}, "θ=" + py_fixed(th, 2) + " rad"});
        return out;
    }
    std::string describe_extra_json() const override { return "\"speed\": " + py_repr(speed_); }

private:
    double speed_;
    Bounds bounds_{};
    std::vector<double> actions_;
};

}  // namespace dhj
```

Register it in `include/dhj/dynamics/registry.hpp`:

```cpp
#include "my_system.hpp"
...
if (args.dynamics == "my_system") { MySystem dyn(args.dt, args.tau, args.velocity); fn(dyn); return true; }
```

and add the name to `registered_dynamics()`. Rebuild; `--dynamics my_system` works.

Notes:

* Actions are addressed by index: `ode(state, action_index)` looks the action up
  in the class's own table, so vector-valued actions (e.g. `(ax, ay)` for the
  double integrator) need no support in the core.
* `plot_dims()` (default `{0, 1}`) selects the two state dimensions of the
  panels; `slices()` pins all the others. For a 4-D system pin two coordinates
  per slice (see `double_integrator.hpp`).
* Periodic dimensions need nothing beyond `periodic_dims()`: wrapping after
  integration, ±2π queries, exclusion from the out-of-bounds test and modulo-2π
  slice matching are handled by `Dynamics<N>::normalize`, `CellTree<N>` and
  `SlicePlotter<N>`.
* `validate()` (called at start-up) checks the periodic boxes, the plot
  dimensions and that the three Lipschitz constants were set.

## Differences from the Python (all numerically negligible)

* ODE: Dormand–Prince 5(4) with atol = rtol = 1e-12 instead of LSODA (`odeint`);
  checkpoint states agree to ~1e-15, cell sets are identical in every validated run.
* Spatial index: BVH with inclusive intersection instead of libspatialindex;
  returns the same leaf sets.
* Figures: SVG instead of matplotlib PNG (use `plot_from_csv.py` for the latter).
