# C++ solver

Header-only C++17 implementation of the paper's algorithms: value iteration on a grid
abstraction for lower and upper bounds on the reach-avoid and avoid value functions
(Algorithm 1), and adaptive refinement of the cells left unclassified (Algorithm 2). One
binary covers all modes and all dynamics. It needs only a C++17 compiler, pthreads and zlib.
All results in `../paper_results/` were computed with it.

```
cpp/
├── Makefile / CMakeLists.txt
├── src/main.cpp                  parse args -> out dir -> run.log tee -> dispatch on --dynamics
└── include/dhj/
    ├── core/
    │   ├── types.hpp      State<N>, Bounds<N>, Mode, number formatting
    │   ├── ode.hpp        Dormand-Prince 5(4), atol = rtol = 1e-12
    │   ├── rtree.hpp      BoxIndex<N>: bulk-built bounding-volume hierarchy
    │   ├── cell.hpp       Cell<N>, CellTree<N> (periodic wrap-around queries, bisection)
    │   ├── reach.hpp      Gronwall growth factors (reach-box radius)
    │   ├── vi.hpp         SafetyValueIterator<N>: successor cache, parallel Bellman sweeps
    │   ├── refine.hpp     AdaptiveRefinement<N>: Algorithm 2's loop, checkpoints, phases.csv
    │   ├── checkpoint.hpp binary checkpoints, --resume auto
    │   ├── svg.hpp        SlicePlotter<N> (figures) + CSV dump
    │   ├── params.hpp     params.json, phases.csv
    │   ├── log.hpp        RunLogger: tee stdout/stderr -> run.log
    │   ├── args.hpp       command line, output-directory naming
    │   ├── pool.hpp       thread pool
    │   └── run.hpp        Algorithm 2 driver (Algorithm 1 = Algorithm 2 with --refinements 0)
    └── dynamics/
        ├── base.hpp              Dynamics<N> interface + Overlay/Slice + dynamics_multi_step
        ├── dubins.hpp            Dubins car (3-D, theta periodic)
        ├── evasion.hpp           3-D evasion (theta periodic)
        ├── double_integrator.hpp planar double integrator (4-D)
        ├── bicycle.hpp           kinematic bicycle (5-D)
        ├── van_der_pol.hpp       van der Pol oscillator (2-D)
        └── registry.hpp          --dynamics name -> class
```

## Build and run

```bash
make                      # -> build/dhj          (or: cmake -B build && cmake --build build)
./build/dhj --help
```

The run of Fig. 1 in the paper (Dubins car, γ = 0.96, 13 refinement iterations):

```bash
./build/dhj --mode ra_discount --gamma 0.96 --dynamics dubins --dt 0.3 \
    --initial-resolution 20 --refinements 13 --eps-min 1e-6 --delta-lower 1e-8 --delta-upper 1e-8
```

Main options (`--help` lists all of them):

* `--mode`: `ra_discount` (discounted reach-avoid with discount factor `--gamma` in (0, 1)),
  `ra_nodiscount` (reach-avoid with γ = 1, Proposition 3) or `avoid_nodiscount` (avoid with
  γ = 1, Section 8). Only `ra_discount` takes `--gamma`.
* `--dynamics`: `dubins`, `evasion`, `double_integrator`, `bicycle`, `van_der_pol`.
* `--dt`: the sampling time t_s. `--tau` (the time a control is held) defaults to it.
* `--initial-resolution`: cells per dimension of the initial uniform grid.
* `--refinements`: number of refinement iterations of Algorithm 2; 0 runs Algorithm 1 only.
* `--eps-min`: Algorithm 2 splits a cell only while its radius (half its largest side)
  exceeds ε_min.
* `--delta-lower`, `--delta-upper`: the stopping tolerances δ̲\* and δ̄\* of Algorithm 1.
  Value iteration stops once a sweep changes the lower bound by less than δ̲\* and the upper
  bound by less than δ̄\* (sup norm). `--vi-iterations` caps the sweeps per phase (default
  100000).
* `--two-abstractions`: upper-bound successors cover each control cell (Proposition 2);
  the paper's experiments do not use it.
* `--restrict-reach-set-to-grid`: ignores the part of each reach box outside the state box
  and adds no sink state, as if f mapped X into X. With `--if-leave-use-cellwise-bounds`, an
  action whose reach box lies entirely outside the box (no successor) uses only the cell's
  own bounds (min{l, r} in reach-avoid, l in the avoid case) instead of counting as -inf.
  Run folders get `_restrictgrid` and `_cellbounds`. For comparison only: the results are not
  sound for systems that can leave X, and the paper's experiments use neither flag (see
  `../old_results/results_stay_in_grid*`). Every value iteration logs how many reach boxes
  leave the state box.

Each run writes to `<--results-root>/<mode>/<dynamics>_res<R>[_epsmin<ε_min>]_dlower<δ̲*>_dupper<δ̄*>[_gamma<γ>]_dt<dt>`
(`_epsmin` only if the run refines, `_gamma` only for `ra_discount`, plus `_tau<τ>` if `--tau`
differs from `--dt`): `run.log`, `params.json`, `phases.csv`
and one figure per phase. Figures are SVG at the default `--plot-dpi 800` and PNG at a
smaller value (the paper's scripts use 300). `--dump-csv` also writes every cell's bounds and
values to `value_function_phase_N.csv`.

### Checkpoints / resume

`checkpoints/checkpoint_phase_NNNN.bin` is written after phase 0 and after every
phase (`--checkpoint-every`, `--keep-checkpoints`). Resume with
`--resume auto` or `--resume PATH` plus the same mode/dynamics/parameters.

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

A Python port of the original paper scripts, with which this solver was cross-checked
during development, is kept in `../old_results/python_port/`; no result in the paper uses it.
