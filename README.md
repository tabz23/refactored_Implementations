# Sound bounds on reach-avoid and avoid value functions — code for the paper

C++ solver for the paper *A Quantitative Approach to Reach–Avoid Reachability Analysis for
Discrete-Time Nonlinear Control Systems*: value iteration on a grid abstraction that computes
lower and upper bounds on the reach-avoid and avoid value functions (Algorithm 1), with
adaptive refinement of the cells it leaves unclassified (Algorithm 2).

**Paper results:** everything needed to reproduce the paper's results, including one script
that reruns and compares all of them, is in [`paper_results/`](paper_results/README.md).

```
refactored_Implementations/
├── README.md        <- this file
├── cpp/             <- the solver (header-only C++17, one binary `dhj`); see cpp/README.md
├── paper_results/   <- the paper's results, their scripts, README.md and reproduce_all.sh
└── old_results/     <- earlier experiments and code that the paper does not use
                        (including a Python port of the original paper scripts)
```

## Quick start

```bash
cd cpp && make
./build/dhj --mode ra_discount --dynamics dubins --dt 0.3 --gamma 0.96 \
    --delta-min 1e-8 --eps-min 1e-6 --initial-resolution 20 --refinements 13
```

This is the run of Fig. 2. It writes to
`refactored_Implementations/results/<mode>/<dynamics>_res<R>_epsmin<ε_min>_dmin<δ>_gamma<γ>_dt<dt>/`
(override with `--results-root`, `--out-dir` or `--tag`), which contains

| file | content |
|---|---|
| `run.log` | everything printed to the terminal (appended across resumes), with a banner holding the command, cwd and timestamps |
| `params.json` | all arguments, the dynamics description (bounds, actions, Lipschitz constants, sets), derived quantities (Grönwall growth factors), environment |
| `phases.csv` | one row per phase: cells, boundary/refinable cells, refined parents, value-iteration sweeps, converged flag, VI time, successor time, phase time, wall time |
| `value_function_phase_N_complete.svg` / `.png` | the figure after phase N (SVG at the default `--plot-dpi 800`, PNG at a smaller value) |
| `value_function_phase_N.csv` | (with `--dump-csv`) every cell: `cell_id, d<k>_lo, d<k>_hi, …, V_lower, V_upper, l_lower, l_upper, r_lower, r_upper` |
| `checkpoints/` | one checkpoint per phase (`--checkpoint-every`, `--keep-checkpoints`) |

## Modes

| `--mode` | value function | Bellman backup of the bounds |
|---|---|---|
| `ra_discount` | discounted reach-avoid, γ = `--gamma` (Algorithm 1) | V ← min(l, max(r, γ·best)) |
| `ra_nodiscount` | reach-avoid with γ = 1 (Proposition 3) | V ← min(l, max(r, best)) |
| `avoid_nodiscount` | avoid with γ = 1 (Section 8) | V ← min(l, best) |

`best` is the maximum over actions of the minimum (lower bound) or maximum (upper bound)
of V over the successor cells. Leaving a non-periodic state bound leads to an absorbing
failure state with value −1. Value iteration stops when the sup-norm changes of both
bounds are below `--delta-min`.

## Dynamics included

| `--dynamics` | state | actions | periodic |
|---|---|---|---|
| `dubins` | (x, y, θ) ∈ [−3,3]²×[−π,π] | ω ∈ {−1, 0, 1} | θ |
| `evasion` | (x, y, θ) ∈ [−3,3]²×[−π,π] | u ∈ 5 evenly spaced values in [−1, 1] | θ |
| `double_integrator` | (x, y, vx, vy) ∈ [−3,3]²×[−v_max,v_max]² | (ax, ay) ∈ {−a, 0, a}² | – |
| `bicycle` | (x, y, θ, v, δ) ∈ [−3,3]²×[−π,π]×[−v_max,v_max]×[−0.6,0.6] | (a, ω) ∈ {−a_max, 0, a_max}×{−1, 0, 1} | θ |
| `van_der_pol` | (x, y) ∈ [−2.5,2.5]×[−2.9,2.9], y = x′ | none | – |

Dubins, evasion, the double integrator and the bicycle use an obstacle disc at the origin
(radius 1.3; 1.0 for evasion) and a target disc of radius 0.5 at (2.5, 0). Van der Pol
(μ = 1) uses an obstacle disc of radius 0.7 at the origin and a target disc of radius 0.4
at (2, 0). The failure and reward functions are Euclidean signed distances, so
L_l = L_r = √2 in the infinity norm. The Lipschitz constants L_f of the dynamics are in
each class and in every run's `params.json`. To add a system, see `cpp/README.md`.

## Checkpoints and resume

A checkpoint is written after phase 0 and after every refinement phase
(`--checkpoint-every N`, `--keep-checkpoints K`, default: every phase, keep 2).
Resume with `--resume auto` (latest checkpoint in the run's `checkpoints/`
directory) or `--resume <path>`, using the **same** mode/dynamics/parameters.
Resuming reproduces the uninterrupted run exactly: the checkpoint holds the full
cell tree and all interval values, and the successor sets — a pure function of
the leaf set — are recomputed on the first local sweep. `run.log` and
`phases.csv` are appended, `params.json` is rotated to `params_<stamp>.json`.
