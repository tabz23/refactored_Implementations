# Quantitative verification of Algorithm 1 (Dubins car)

This directory holds two runs of Algorithm 1 (no refinement) for the Dubins car on
a uniform 120 × 120 × 120 grid, their verification, the tables and figures, and the
scripts that produce and replicate all of it:

* reach-avoid with γ = 0.96: Proposition 4, labelled `prop:greedy`, and the soundness
  of both value bounds. Its first claim: under the greedy policy, safe arrival in T within
  ⌊log(V_lower(s0)/r*)/log γ⌋ steps. Its second claim (violation, checked by
  `verify_violation.py`): if V_upper(s0) < 0, every policy reaches F, before T, within
  ⌊log(V_upper(s0)/l_min)/log γ⌋ steps;
* avoid with γ = 1: Proposition 6, labelled `prop:avoid` (under its policy π*,
  l(x_t) ≥ V_lower(s0) for all t) and the soundness of both value bounds. Its second
  claim (failure-reaching time, checked by `verify_avoid_failure_time.py`): if A = U and
  V_upper^(k)(s0) ≤ 0, every trajectory from [[s0]] reaches F within k steps, under every
  control sequence.

## Reproduce

```bash
./reproduce.sh            # replicates everything in ../results_quantitative_verification_rerun
                          # and compares the numbers with the ones here
./reproduce.sh DIR        # the same, into DIR
./reproduce.sh "$PWD"     # regenerates the results in this directory
```

This takes about four minutes and needs a C++17 compiler, Python with numpy, scipy
and matplotlib, and optionally LaTeX for the figure text (`reproduce.sh` lists
the recorded versions). `provenance.txt` records the commit (4175b91, no
uncommitted solver changes; the original run in `../original_runs/` used d8463b7),
the compiler and the solver binary hash; the binary embeds its build time, so the
hash changes with every build. `reproduce.log` is the output of the last replication,
ending with the comparison of every report with the one here.

**Only the failure-reaching time of Proposition 6(2)** (paper, Section 9.2.4, avoid case:
"for each cell s with V̄₁^cons(s) ≤ 0, we simulated the system starting from one random
state in its interior and from its eight corners while following every possible action
sequence until it reaches F. In every resulting trajectory, F was reached at some time
τ ≤ k, where k is the first iteration index satisfying V̄₁^(k)(s) ≤ 0, as the second claim
of Proposition 6 states."),
in about a minute:

```bash
./run_alg1_res120.sh avoid DIR     # Algorithm 1, avoid case, 120^3 cells; also writes
                                   # upper_nonpos_iteration_phase_0.csv, the k(s) of every cell
python verify_avoid_failure_time.py \
    DIR/avoid_nodiscount/dubins_res120_dlower1e-300_dupper1e-300_dt0.3 \
    --out DIR/verification_avoid_failure_time --seed 0
```

The paper's numbers are, in `verification_avoid_failure_time/report.json`:
`cells_V_upper_nonpositive` (464,032) and `exhaustive.violations` (0, "always happens").
The report also gives `exhaustive.states` (4,176,288) and the share of start states with
k(s0) ≥ 1 for which some sequence first reaches F exactly at k(s0): 63 % =
`exhaustive.latest_hit_equals_k_(k>=1)` / (9 × `certificate.cells_k_positive`) =
1,280,735 / 2,024,928. The paper does not report these two.

## Contents

| path | content |
|---|---|
| `run_alg1_res120.sh` | builds `cpp/build/dhj_head` and runs the solver for both cases |
| `verify_quantitative.py` | the verifier; writes the tables, figures and reports |
| `verify_violation.py` | the check of the violation claim of Proposition 4; imports `verify_quantitative.py` |
| `verify_avoid_failure_time.py` | the check of the failure-reaching time of Proposition 6(2); imports `verify_quantitative.py` |
| `reproduce.sh` | runs the four scripts and compares the numbers |
| `ra_discount/dubins_res120_…/`, `avoid_nodiscount/dubins_res120_…/` | solver output: `value_function_phase_0.csv` (every cell: bounds, V_lower, V_upper, l and r bounds), `params.json`, `phases.csv`, `run.log`, figure; in the avoid run also `upper_nonpos_iteration_phase_0.csv` (every cell: k(s), the first iteration with V_upper^(k)(s) ≤ 0, or −1) |
| `verification_ra/` | `table_ra.tex`, figure `ra_verification.{pdf,png}`, per-panel figures, `summary.json`, `dubins/report.json` |
| `verification_avoid/` | `table_avoid.tex`, figure `avoid_verification.{pdf,png}`, per-panel figures, `summary.json`, `dubins/report.json` |
| `verification_ra_violation/` | `report.json` (every number), `sampled_states.npz` and `exhaustive.npz` (per start state) |
| `verification_avoid_failure_time/` | `report.json` (every number), `exhaustive.npz` (per start state: k(s0), latest first hit, number of sequences) |
| `solver.log`, `verification_ra.log`, `verification_avoid.log`, `verification_ra_violation.log`, `verification_avoid_failure_time.log` | console output of the last run |
| `reproduce.log` | console output of the last replication (`./reproduce.sh`), with the comparisons |

The CSV, NPY and NPZ files are excluded from git by `.gitignore` (hundreds of MB).
`*/dubins/trajectories.npz` holds every simulated initial state and its outcome.

## Method

**Solver (C++, `cpp/build/dhj_head`).** The state box [−3, 3]² × [−π, π] is split
into 120³ = 1,728,000 equal cells. For each cell, l_lower, l_upper, r_lower, r_upper
are the center values ∓ √2·r, where r is the largest half-width of the cell
(∞-norm). For each cell and action, the center is integrated over τ = 0.3 s; the
successors are the cells meeting the box of ∞-norm radius r·exp(L_f τ) plus the
integrator's error budget around the result, and a sink is added if the box leaves the
state box along x or y. The sink has upper value −1 and lower value −1 (reach-avoid) or
min(−1, −1.3) = −1.3 (avoid; paper, Section 8). The actions are {−1, 0, 1} = U, so A = U.

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
   π̂(s) = argmax_a min_{s′ ∈ Δ(s,a)} V_lower(s′) (ties: first action). With
   Q(s, a) = min{ l_lower(s), max{ r_lower(s), γ · min_{s′ ∈ Δ(s,a)} V_lower(s′) } },
   which is non-decreasing in that minimum, π̂(s) also maximizes Q(s, ·).
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
   reached no later than t*; t* ≤ B(s0) = ⌊log(V_lower(s0)/r*)/log γ⌋ with
   r* = max_s r_lower(s) = 0.4276; γ^t l(x_t) ≥ V_lower(s0) for t ≤ t*; and the achieved
   J = max_{t ≤ t*} min{ γ^t r(x_t), min_{τ ≤ t} γ^τ l(x_τ) } satisfies
   V_lower(s0) ≤ J ≤ V_upper(s0).
7. *Time bound.* B(s0) is computed by `time_bound`: the largest k with
   γ^k r* ≥ V_lower(s0), with the product formed by repeated multiplication, as in the
   solver. The bound is attained (t* = B(s0)) by 576 trajectories, and 77 of them have
   V_lower(s0) = γ^k r* exactly as the solver computed it, where `log` or `pow` can be
   off by one ulp. On this run, B(s0) equals `floor(log(V_lower(s0)/r*)/log γ)` for every
   trajectory.

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

Violation checks (`verification_ra_violation/`, script `verify_violation.py`, on the
reach-avoid run). The second claim of Proposition 4: if V_upper(s0) < 0, then for every
policy and every x0 in [[s0]], the trajectory reaches F at some time
t ≤ t̄_F(s0) = ⌊log(V_upper(s0)/l_min)/log γ⌋ with l(x_t) ≤ γ^(−t) V_upper(s0), and
r(x_τ) ≤ γ^(−τ) V_upper(s0) < 0 (so x_τ ∉ T) for τ < t. The solver gives the sink
outside the grid the fixed value −1, and its l and r bounds are −1 as well
(`cpp/include/dhj/core/vi.hpp`); only in the avoid run is its lower value −1.3, which
these reach-avoid checks do not use. The verifier models it as an absorbing state with
l = r = −1, whose value is −1. Hence l_min = −1.3, at the obstacle center. For a
trajectory, t_hit is the first time with γ^t l(x_t) ≤ V_upper(s0).

1. *Upper certificate.* As in check 5 of the reach-avoid checks.
2. *Start states.* One state drawn uniformly (`numpy.random.default_rng(0)`) and the
   8 corners (vertices) of every cell with V_upper < 0: 458,296 × 9 = 4,124,664 states.
3. *Policies.* Each constant action, uniformly random actions (3 seeds), π̂ of the
   first claim, and the one-step greedy policy on V_upper (argmax_a V_upper of the cell
   containing Φ(x, a, τ), with the sink counted as −1). For every trajectory: t_hit
   exists, t_hit ≤ t̄_F(s0), and γ^τ r(x_τ) ≤ V_upper(s0) for τ < t_hit.
4. *Exhaustive search.* From every start state, every action sequence is followed until
   t_hit. A branch is a violation if it reaches t̄_F(s0) without t_hit, enters the sink
   too late (γ^t < −V_upper(s0)), or has γ^τ r(x_τ) > V_upper(s0) before t_hit. Every
   policy, stationary or not, generates one of these sequences, so this covers every
   policy from these states.
5. *Rounding.* γ^t is formed by repeated multiplication, as in the solver. Some cells
   attain the bound exactly: V_upper(s0) = −γ^k, and the trajectory leaves the grid at
   step k. With `pow`, 5,924 of these ties show up as one-ulp misses. Hits within 1e-12
   would be reported as ties; there are none. t̄_F(s0) is computed by `time_bound`, as
   B(s0).

Failure-reaching time checks (`verification_avoid_failure_time/`, script
`verify_avoid_failure_time.py`, on the avoid run). The second claim of Proposition 6: if
A = U and V_upper^(k)(s0) ≤ 0 for some k ≥ 0, then under every control sequence every
trajectory starting in [[s0]] reaches F = {l ≤ 0} at some time τ ≤ k. The smallest such
k is k(s0), the first iteration of Algorithm 1 with V_upper^(k)(s0) ≤ 0. The solver
records it during value iteration, which is Jacobi (each sweep writes a new array), so
its iteration k is the paper's V_upper^(k), and writes it to
`upper_nonpos_iteration_phase_0.csv`. As in the solver, l = −1 outside [−3, 3]², so
leaving the grid reaches F.

1. *Recomputation.* The same decreasing iteration W^(k+1) = min(l_upper, max_a max_{s′ ∈
   Δ(s,a)} W^(k)(s′)) from W^(0) = l_upper on the verifier's own abstraction (sink value
   −1), with the first k at which W^(k)(s) ≤ 0. It must agree with the solver's k(s) and
   V_upper in every cell.
2. *Certificate.* The induction of the proof, checked cell by cell on the verifier's
   abstraction: k(s) = 0 only if l_upper(s) ≤ 0 (then the whole cell lies in F), and for
   k(s) ≥ 1, l_upper(s) > 0 and every successor s′ of every action has k(s′) ≤ k(s) − 1
   (the sink counts as k = 0). The abstraction is sound, so this proves the bound for
   every state of these cells, not only for the sampled ones.
3. *Start states.* One state drawn uniformly (`numpy.random.default_rng(0)`) and the
   8 corners of every cell with V_upper ≤ 0: 464,032 × 9 = 4,176,288 states.
4. *Exhaustive search.* From every start state, every action sequence is followed with
   the exact flow until it first reaches F (l(x_t) ≤ 0 or leaving the grid). A branch
   that has not reached F at time k(s0) is a violation. Every control sequence, from
   any policy, stationary or not, is one of these sequences. This replaces check 7 of the
   avoid checks, which follows five controls from one state per cell.

## Results

Reach-avoid, γ = 0.96 (`verification_ra/table_ra.tex`):

| | Dubins car |
|---|---|
| cells certified in R / L / unclassified (%) | 56.9 / 26.5 / 16.5 |
| lower / upper certificate violations | 0 / 0 |
| simulated initial states | 983,988 |
| reached T without entering F or leaving the grid | 100 % |
| reach-time bound violations | 0 |
| median observed reach time / median bound (steps) | 19 / 36 |
| largest observed reach time / bound | 0.73 (40/55) |
| trajectories with t* = B(s0) (bound attained) / with B(s0) = 0 | 576 / 480 |
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

Violation (second claim of Proposition 4), γ = 0.96 (`verification_ra_violation/report.json`):

| | Dubins car |
|---|---|
| cells with V_upper < 0 / with V_upper = 0 | 458,296 / 0 |
| min V_upper / l_min | −1.228 / −1.3 |
| t̄_F(s0): median / max (steps) | 15 / 152 |
| upper certificate violations | 0 |
| start states (uniform + vertices) | 4,124,664 |
| violations under the 8 policies | 0 |
| exhaustive search: action sequences followed / largest frontier | 40,024,654 / 555 |
| exhaustive search: violations | 0 |
| latest t_hit over all sequences: median / max (steps) | 1 / 8 |
| largest (latest t_hit) / t̄_F(s0) | 0.5 |
| start states with t_hit = 0 (deep enough in F at t = 0) | 1,580,225 |
| min V_upper(s0) − γ^t l(x_t) at t_hit / min V_upper(s0) − γ^τ r(x_τ) before t_hit | 0 (bound attained) / 0.0017 |
| min γ^τ r(x_τ) / V_upper(s0) before t_hit (≥ 1 required) | 1.006 |
| verifier time (s) | about 18 |

Failure-reaching time (second claim of Proposition 6), γ = 1
(`verification_avoid_failure_time/report.json`):

| | Dubins car |
|---|---|
| cells with V_upper ≤ 0 / with V_upper = 0 | 464,032 / 0 |
| cells with k = 0 (l_upper ≤ 0) / with k ≥ 1 | 239,040 / 224,992 |
| k(s0) over the start states: median / max (steps) | 0 / 10 |
| cells where the verifier's k and V_upper equal the solver's | 1,728,000 (all) |
| certificate violations | 0 |
| start states (uniform + vertices) | 4,176,288 |
| exhaustive search: action sequences followed / largest frontier | 37,922,424 / 555 |
| exhaustive search: violations (F not reached by k(s0)) | 0 |
| latest first time in F over all sequences: max (steps) | 8 |
| start states inside F at t = 0 | 2,198,073 |
| start states with k(s0) ≥ 1 where some sequence first reaches F exactly at k(s0) | 1,280,735 of 2,024,928 (63 %) |
| verifier time (s) | about 33 |

## Scope

The checks concern the sampled system at the sampling instants. Leaving the grid
counts as failure (the solver's sink), so cells with V_upper ≤ 0 are certified with
respect to the failure set extended by the region outside the grid; 189,912 of the
sampled avoid BRS states fail by leaving the grid. The certificates cover every state
of a cell and all times; the simulations cover one state per cell and, in the avoid
case, 1000 steps. The violation and failure-reaching time checks cover nine start
states per cell with V_upper < 0 (reach-avoid) or V_upper ≤ 0 (avoid) and every action
sequence from them; under the constant and random policies, 1,657,584 of the 4,124,664
reach-avoid trajectories reach the failure set by leaving the grid.
