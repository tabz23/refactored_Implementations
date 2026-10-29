# Quantitative verification of Algorithm 1 (Dubins car)

This directory holds two runs of Algorithm 1 (no refinement) for the Dubins car on
a uniform 120 × 120 × 120 grid, their verification, the tables and figures, and the
scripts that produce and replicate all of it:

* reach-avoid with γ = 0.96: Proposition 4, labelled `prop:greedy` (under the
  greedy policy, safe arrival in T within ⌈log(V_lower(s0)/r*)/log γ⌉ steps) and the
  soundness of both value bounds;
* avoid with γ = 1: Proposition 6, labelled `prop:avoid` (under its policy π*,
  l(x_t) ≥ V_lower(s0) for all t) and the soundness of both value bounds.

## Reproduce

```bash
./reproduce.sh            # replicates everything in ../results_quantitative_verification_rerun
                          # and compares the numbers with the ones here
./reproduce.sh DIR        # the same, into DIR
./reproduce.sh "$PWD"     # regenerates the results in this directory
```

This takes about three minutes and needs a C++17 compiler, Python with numpy, scipy
and matplotlib, and optionally LaTeX for the figure text (`reproduce.sh` lists
the recorded versions). `provenance.txt` records the commit (d8463b7, no
uncommitted solver changes), the compiler and the solver binary hash; the binary
embeds its build time, so the hash changes with every build.

## Contents

| path | content |
|---|---|
| `run_alg1_res120.sh` | builds `cpp/build/dhj_head` and runs the solver for both cases |
| `verify_quantitative.py` | the verifier; writes the tables, figures and reports |
| `reproduce.sh` | runs both scripts and compares the numbers |
| `ra_discount/dubins_res120_…/`, `avoid_nodiscount/dubins_res120_…/` | solver output: `value_function_phase_0.csv` (every cell: bounds, V_lower, V_upper, l and r bounds), `params.json`, `phases.csv`, `run.log`, figure |
| `verification_ra/` | `table_ra.tex`, figure `ra_verification.{pdf,png}`, per-panel figures, `summary.json`, `dubins/report.json` |
| `verification_avoid/` | `table_avoid.tex`, figure `avoid_verification.{pdf,png}`, per-panel figures, `summary.json`, `dubins/report.json` |
| `solver.log`, `verification_ra.log`, `verification_avoid.log` | console output of the last run |

The CSV, NPY and NPZ files are excluded from git by `.gitignore` (hundreds of MB).
`*/dubins/trajectories.npz` holds every simulated initial state and its outcome.

## Method

**Solver (C++, `cpp/build/dhj_head`).** The state box [−3, 3]² × [−π, π] is split
into 120³ = 1,728,000 equal cells. For each cell, l_lower, l_upper, r_lower, r_upper
are the center values ∓ √2·r, where r is the largest half-width of the cell
(∞-norm). For each cell and action, the center is integrated over τ = 0.3 s; the
successors are the cells meeting the box of ∞-norm radius r·exp(L_f τ) plus the
integrator's error budget around the result, and a sink with value −1 is added if
the box leaves the state box along x or y. The actions are {−1, 0, 1} = U, so A = U.

* Reach-avoid: γ = 0.96, V_lower starts at min(l_lower, r_lower) and V_upper at
  l_upper. The iteration stops when both sup-norm changes are below 1e-8
  (71 iterations).
* Avoid: γ = 1 with the avoid operators, both bounds start at l_lower and
  l_upper. The iterates take values in a finite set, so the tolerance 1e-300 stops
  the iteration only when nothing changes, as in the paper (89 iterations; the last
  change is exactly 0).

**Verifier (Python, `verify_quantitative.py`).** It reads only the dumped grid,
values and l/r bounds and recomputes everything else independently of the solver.

1. *Grid, flow and cell bounds.* It checks that the grid is uniform and complete,
   and recomputes the cell bounds from the centers (largest difference 1.8e-15).
   The flow over τ is the closed-form solution (checked against `solve_ivp`,
   largest difference 8.9e-16).
2. *Abstraction.* For each cell center and action, the successors are the cells
   whose closed boxes meet the box of radius r·exp(L_f τ) around the exact flow
   (θ wraps around), plus the sink if that box leaves [−3, 3]². This is a sound
   abstraction of the sampled system.

Reach-avoid checks (`verification_ra/`):

3. *Policy.* The greedy policy of the proposition:
   Q(s, a) = min{ l_lower(s), max{ r_lower(s), γ · min_{s′ ∈ Δ(s,a)} V_lower(s′) } },
   π̂(s) = argmax_a Q(s, a) (ties: first action).
4. *Lower certificate.* max_a Q(s, a) ≥ V_lower(s) for every cell. For cells with
   V_lower > 0 this is the inequality the proof of the proposition starts from, so
   the guarantee holds for every state in those cells, not only sampled ones.
5. *Upper certificate.* min{ l_upper(s), max{ r_upper(s), γ · max_a max_{s′ ∈ Δ(s,a)}
   V_upper(s′) } } ≤ V_upper(s) for every cell. With A = U this gives V_γ ≤ V_upper,
   so the cells with V_upper ≤ 0 are correctly certified.
6. *Simulation.* One state drawn uniformly in every cell with V_lower > 0
   (`numpy.random.default_rng(0)`), closed loop x_{t+1} = Φ(x_t, π̂(s_t), τ) with the
   exact flow, followed until t*, the first time with V_lower(s_t) ≤ r_lower(s_t)
   (the time used in the proof), or until it enters F (l ≤ 0) or leaves [−3, 3]².
   For every trajectory: no entry into F and no exit up to t*; the target (r > 0)
   reached no later than t*; t* ≤ B(s0) = ⌈log(V_lower(s0)/r*)/log γ⌉ with
   r* = max_s r_lower(s) = 0.4276; γ^t l(x_t) ≥ V_lower(s0) for t ≤ t*; and the achieved
   J = max_{t ≤ t*} min{ γ^t r(x_t), min_{τ ≤ t} γ^τ l(x_τ) } satisfies
   V_lower(s0) ≤ J ≤ V_upper(s0).

Avoid checks (`verification_avoid/`):

3. *Policy.* The policy of the proposition, π*(s) = argmax_a min_{s′ ∈ Δ(s,a)}
   V_lower(s′) (ties: first action).
4. *Lower certificate.* V_lower(s) ≤ l_lower(s) and min_{s′ ∈ Δ(s,π*(s))} V_lower(s′) ≥
   V_lower(s) for every cell with V_lower > 0. This is the invariance the proof uses,
   so every state of such a cell keeps l(x_t) ≥ V_lower(s0) for all t ≥ 0 under π*.
5. *Upper bound.* The decreasing iteration W ← min(l_upper, max_a max_{s′} W(s′))
   from W = l_upper is run on the verifier's abstraction until nothing changes, and
   V_upper ≥ W must hold in every cell. Every iterate bounds V_1 from above, so the
   cells with V_upper < 0 lie in the BRS. (A super-solution check would not suffice
   for the avoid problem: W ≡ −1 passes it.)
6. *Simulation of the safe cells.* One state drawn uniformly in every cell with
   V_lower > 0, closed loop under π* for 1000 steps (300 s). For every trajectory: no
   entry into F (l ≤ 0) and no exit from [−3, 3]²; l(x_t) ≥ V_lower(s0) at every step; and the cells
   visited never have V_lower below V_lower(s0).
7. *Simulation of the BRS cells.* One state drawn uniformly in every cell with
   V_upper < 0, under five controls: each constant action, the greedy action on W,
   and uniformly random actions. Every trajectory must reach l < 0 or leave the grid
   within k̄(s0) steps, the first iteration with W^(k)(s0) < 0 (W^(k)(s0) < 0 means
   that every control sequence does so within k steps).

## Results

Reach-avoid, γ = 0.96 (`verification_ra/table_ra.tex`):

| | Dubins car |
|---|---|
| cells certified in R / L / unclassified (%) | 56.9 / 26.5 / 16.5 |
| lower / upper certificate violations | 0 / 0 |
| simulated initial states | 983,988 |
| reached T without entering F or leaving the grid | 100 % |
| reach-time bound violations | 0 |
| median observed reach time / median bound (steps) | 19 / 37 |
| largest observed reach time / bound | 0.71 |
| min γ^t l(x_t) / V_lower(s0) (≥ 1 required) | 1.021 |
| min J − V_lower(s0) (≥ 0 required) / max J − V_upper(s0) (≤ 0 required) | 0.0018 / −0.0020 |
| abstraction + value iteration time (s) | 2.1 |

Avoid, γ = 1 (`verification_avoid/table_avoid.tex`):

| | Dubins car |
|---|---|
| cells with V_lower > 0 / V_upper < 0 / unclassified (%) | 56.3 / 26.9 / 16.8 |
| lower / upper violations | 0 / 0 |
| cells where the recomputed upper bound W equals V_upper | 1,728,000 (all) |
| simulated initial states (V_lower > 0) | 973,616 |
| outside F and in the grid for 1000 steps | 100 % |
| min_t l(x_t) / V_lower(s0) (≥ 1 required) | 1.002 |
| min_t l(x_t) − V_lower(s0) (≥ 0 required) | 0.0021 |
| sampled states (V_upper < 0) × controls | 464,032 × 5 |
| entered F or left the grid within k̄(s0) steps | 100 % (k̄ ≤ 10) |
| of the sampled BRS states: inside F at t = 0 / left the grid / entered F later | 244,402 / 189,912 / 29,718 |
| abstraction + value iteration time (s) | 2.3 |

## Scope

The checks concern the sampled system at the sampling instants. Leaving the grid
counts as failure (the solver's sink), so cells with V_upper ≤ 0 (reach-avoid) or
V_upper < 0 (avoid) are certified with respect to the failure set extended by the
region outside the grid; 189,912 of the sampled avoid BRS states fail by leaving
the grid. The certificates cover every state of a cell and all times; the
simulations cover one state per cell and, in the avoid case, 1000 steps.
