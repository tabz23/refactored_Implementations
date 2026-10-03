# Discrete Hamilton–Jacobi reach-avoid solver — modular implementations

Modular Python and C++ implementations of the adaptive interval value iteration
used in the paper (`Discrete_HJ_submission/{RA_nodiscount,RA_discount,avoid_nodiscount}.py`).
Both reproduce the paper scripts **exactly** (same cells, same values; see
[Validation](#validation)) and add:

* one code base for the three problem types (`--mode ra_nodiscount | ra_discount | avoid_nodiscount`),
* dynamics as plug-in files — each file holds the ODE, the reach/avoid sets, the
  Lipschitz constants **and** the plotting (slices, overlays) for one system,
  including correct handling of periodic (angle) coordinates,
* a 4-D example (planar double integrator) next to the paper's two 3-D systems,
* per-phase checkpoints with `--resume`, `run.log`, `params.json`, `phases.csv` in
  every results directory,
* parallelism (multiprocessing pool in Python, thread pool in C++).

```
refactored_Implementations/
├── README.md               <- this file
├── python/                 <- Python package (reference semantics, matplotlib figures)
│   ├── run.py              <- CLI
│   ├── plot_from_csv.py    <- render a CSV dump (e.g. from the C++) with matplotlib
│   └── dhj/
│       ├── core/           <- cell tree, successors, value iteration, refinement, checkpoints, logging
│       └── dynamics/       <- base.py, dubins.py, evasion.py, double_integrator.py, plotting.py
└── cpp/                    <- C++17 port, header-only core, single binary `dhj`
    ├── src/main.cpp
    ├── include/dhj/core/   <- same decomposition as the Python core
    ├── include/dhj/dynamics/ <- base.hpp, dubins.hpp, evasion.hpp, double_integrator.hpp, registry.hpp
    └── tools/compare_checkpoints.cpp
```

Language-specific details: [python/README.md](python/README.md), [cpp/README.md](cpp/README.md).

## Quick start

```bash
# Python
cd python && pip install -r requirements.txt
python run.py --mode ra_nodiscount --dynamics dubins --dt 0.3 --tau 0.3 --gamma 0.9 \
    --delta-min 1e-8 --delta-max 1e-8 --epsilon 0.05 --initial-resolution 20 \
    --vi-iterations 100033 --conservative

# C++
cd cpp && make
./build/dhj --mode ra_nodiscount --dynamics dubins --dt 0.3 --tau 0.3 --gamma 0.9 \
    --delta-min 1e-8 --delta-max 1e-8 --epsilon 0.05 --initial-resolution 20 \
    --vi-iterations 100033 --conservative
```

Both write to `refactored_Implementations/results/<mode>/<dynamics>_res<R>_eps<ε>_dmin<δ>_dmax<δmax>_gamma<γ>_dt<dt>_tau<τ>_<cons|std>/`
(override with `--out-dir`, `--results-root`, `--tag`). The directory contains

| file | content |
|---|---|
| `run.log` | everything printed to the terminal (appended across resumes), with a banner holding the command, cwd and timestamps |
| `params.json` | all CLI arguments, the dynamics description (bounds, actions, Lipschitz constants, sets), derived quantities (η_min, Grönwall growth factors), environment |
| `phases.csv` | one row per phase: `|cells|`, boundary/refinable cells, refined parents, VI sweeps, converged flag, VI time, successor time, phase time, wall time |
| `value_function_phase_N_complete.png` / `.svg` | the figure after phase N (Python: matplotlib PNG; C++: SVG) |
| `value_function_phase_N.csv` | (with `--dump-csv`) every leaf: `cell_id, d<k>_lo, d<k>_hi, …, V_lower, V_upper, l_lower, l_upper, r_lower, r_upper` |
| `checkpoints/` | one checkpoint per phase (`--checkpoint-every`, `--keep-checkpoints`) |

The algorithm, parameters and output format are the same in both languages; the
C++ is roughly two to three orders of magnitude faster at large cell counts.

## Modes

| `--mode` | Bellman backup | original script |
|---|---|---|
| `ra_nodiscount` | V ← min(l, max(r, best)) | `RA_nodiscount.py` |
| `ra_discount` | V ← min(l, max(r, γ·best)) | `RA_discount.py` |
| `avoid_nodiscount` | V ← min(l, best) | `avoid_nodiscount.py` |

`best = max_u (min / max over successors of V)`; an empty successor set is −∞;
leaving a non-periodic state bound leads to the out-of-bounds sink with value −1.
`--conservative` enables the Algorithm-3 stopping rule (δ^k ≥ −δ_max ∧ ‖ΔV̄‖ < δ_min,
then V_lower −= ε_cons + ε_machine).

## Dynamics included

| `--dynamics` | state | actions | L_f | L_l = L_r | periodic |
|---|---|---|---|---|---|
| `dubins` | (x, y, θ) ∈ [−3,3]²×[−π,π] | ω ∈ {−1,0,1} | v | √2 | θ |
| `evasion` | (x, y, θ) | u ∈ linspace(−1,1,5) | 1+v | √2 | θ |
| `double_integrator` | (x, y, vx, vy) ∈ [−3,3]²×[−v_max,v_max]² | (ax, ay) ∈ {−a,0,a}² | 1 | √2 | – |
| `van_der_pol` | (x, y) ∈ [−2.5,2.5]×[−2.9,2.9], y = x′ | none (the singleton {0}) | 20.75 | √2 | – |

Dubins, evasion, and the double integrator use an obstacle disc at the origin
(radius 1.3; 1.0 for evasion) and a target disc of radius 0.5 at (2.5, 0).
The double integrator treats the velocity box as a state constraint (leaving it
is out-of-bounds, like leaving the position box). Its L_f = 1 because
‖f(x)−f(y)‖_∞ = max(|Δvx|, |Δvy|) ≤ ‖x−y‖_∞. For all four systems l and r are
Euclidean signed distances, so L_l = L_r = √2 in the infinity norm.

Van der Pol is autonomous, μ = 1. Its limit cycle (period ≈ 6.67) lies in
[−2.01, 2.01] × [−2.68, 2.68] and stays at least 1.53 from the origin. The
state box is only slightly larger than that cycle. The obstacle is a disc of
radius 0.7 at the origin (inside the cycle); the target is a disc of radius
0.4 at (2, 0), which the cycle passes through. L_f = 20.75 is the maximum
∞-norm of Df on this box. Use a short horizon, `--dt 0.02 --tau 0.02`,
so the Gronwall tube stays thin enough for the positive set to grow along the
orbit instead of filling the box in one step.

## Adding a new dynamics (summary)

1. Copy `python/dhj/dynamics/dubins.py` (or `cpp/include/dhj/dynamics/dubins.hpp`)
   to a new file and fill in: `name`, `dim`, `periodic_dims`, `state_names`,
   state bounds, action list, `ode`, `failure_function`, `reward_function`,
   the Lipschitz constants `L_f`, `L_l`, `L_r`, and the plotting hooks
   `overlays()` and `slices()`.
2. Register it: Python — add the class to `REGISTRY` in `dhj/dynamics/__init__.py`;
   C++ — add an `if` block to `with_dynamics` in `dynamics/registry.hpp`.
3. Run with `--dynamics <name>`.

Periodic (angle) coordinates: list them in `periodic_dims` and give them the box
[−π, π]. The core then (a) wraps them with `arctan2(sin, cos)` after integration,
(b) queries the spatial index three times (shift 0, ±2π) so reach boxes crossing
±π find the cells on the other side, (c) never treats them as out-of-bounds, and
(d) matches plot slices modulo 2π. Nothing else is needed. Full instructions with
code templates are in the language READMEs.

## Checkpoints and resume

A checkpoint is written after phase 0 and after every refinement phase
(`--checkpoint-every N`, `--keep-checkpoints K`, default: every phase, keep 2).
Resume with `--resume auto` (latest checkpoint in the run's `checkpoints/`
directory) or `--resume <path>`, using the **same** mode/dynamics/parameters.
Resuming reproduces the uninterrupted run exactly: the checkpoint holds the full
cell tree and all interval values, and the successor sets — a pure function of
the leaf set — are recomputed on the first local sweep. `run.log` and
`phases.csv` are appended, `params.json` is rotated to `params_<stamp>.json`.

## Faithfulness to the paper scripts — things to know

* **Phase 0 tolerance.** Phase 0 uses `--delta-min`, the same residual tolerance as later phases. `--phase0-tol` overrides that for phase 0 only.
* **Avoid-mode figures.** `avoid_nodiscount.py` left the upper/lower value
  panels blank (the drawing loop is commented out). The new code draws them;
  the classification panel and all numbers are unchanged. No target circle is
  drawn in avoid mode.
* **Numerics.** Python uses `scipy.integrate.odeint` (atol = rtol = 1e-12) and
  the `rtree` index, exactly as the scripts; the Python package is bit-identical
  to them. The C++ uses a Dormand–Prince 5(4) integrator at the same tolerances
  and a BVH; it agrees with the Python to ~1e-15 in the values, with identical
  cell sets.
* Everything else (grid construction order, cell ids, refinement order,
  Grönwall growth factors, Bellman backups, stopping rules, printed log lines)
  follows the scripts line by line.

## Validation

The three runs stored in `Discrete_HJ_submission/results_all` (dubins res 20,
evasion res 30, avoid res 80; ε = 0.05, δ_min = δ_max = 1e-8, dt = τ = 0.3,
conservative) were re-run with the new code and compared with the reference
logs (`|cells|`, boundary cells, VI sweeps per phase) and, cell by cell, with
the original Python / the earlier C++ port. See
[VALIDATION.md](VALIDATION.md) for the tables and the exact commands.
