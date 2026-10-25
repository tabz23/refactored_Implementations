# Verification of the quantitative guarantees

Checks two guarantees of the paper on Algorithm 1 runs (no refinement) for the
Dubins car over a uniform 120³ grid:

* reach-avoid with γ = 0.96: the proposition on the greedy policy (label
  `prop:greedy`), i.e. safe arrival within ⌈log(V_lower(s0)/r*)/log γ⌉ steps, and the
  soundness of both value bounds;
* avoid with γ = 1: the proposition on the avoid value function (label
  `prop:avoid`), i.e. dist(x_t, F) ≥ V_lower(s0) for all t under its policy, and the
  soundness of both value bounds.

The results, the tables and figures, and the copies of these scripts that produced
them are in `../results_quantitative_verification/`; its `reproduce.sh` regenerates
or replicates everything. The scripts here are identical to those copies.

## Scripts

* `run_alg1_res120.sh [ra|avoid|all] [ROOT]` builds `cpp/build/dhj_head` and runs
  the solver with `--refinements 0 --dump-csv`. The reach-avoid run stops when
  both sup-norm changes are below 1e-8. The avoid run stops only when the iterates
  no longer change (tolerance 1e-300, as in the paper; the iterates take values in
  a finite set).
* `verify_quantitative.py RUN_DIR ... --out DIR` checks runs of one mode
  (`ra_discount` or `avoid_nodiscount`) and writes the tables, figures and reports.
  Needs numpy, scipy and matplotlib; LaTeX is used for figure text when available.

## What is checked

The verifier reads only the dumped grid, values and l/r bounds. It recomputes
the abstraction independently of the solver: the exact flow of each cell center
over τ (closed form, checked against `solve_ivp`), a reach box of ∞-norm radius
r·exp(L_f τ) (Grönwall), the cells meeting it (θ periodic), and the out-of-bounds
sink (value −1).

Reach-avoid (`ra_discount`):

1. Cell bounds l_lower, l_upper, r_lower, r_upper recomputed from the centers.
2. Lower certificate: (T̂ V_lower)(s) ≥ V_lower(s) for every cell, with the greedy
   action of the proposition. For cells with V_lower > 0 this is the inequality the
   proof starts from, so the guarantee holds for every state in those cells.
3. Upper certificate: (T̄ V_upper)(s) ≤ V_upper(s) for every cell, which with A = U
   gives V_γ ≤ V_upper (cells with V_upper ≤ 0 are sound).
4. Closed-loop simulation with the exact flow from one uniformly random state in
   every cell with V_lower > 0, followed until t* (first time with
   V_lower(s_t) ≤ r_lower(s_t)): no entry into F and no exit from the grid, the
   target reached by t*, t* ≤ ⌈log(V_lower(s0)/r*)/log γ⌉, γ^t l(x_t) ≥ V_lower(s0)
   for t ≤ t*, and V_lower(s0) ≤ max_{t≤t*} g_γ ≤ V_upper(s0).

Avoid (`avoid_nodiscount`):

1. Cell bounds l_lower, l_upper recomputed from the centers.
2. Lower certificate: V_lower(s) ≤ l_lower(s) and min_{s′∈Δ(s,π(s))} V_lower(s′) ≥
   V_lower(s) for every cell with V_lower > 0, where π(s) maximizes
   min_{s′∈Δ(s,a)} V_lower(s′) (the policy of the proposition). This is the
   invariance the proof uses, so every state of such a cell keeps
   l(x_t) ≥ V_lower(s0) for all t under π.
3. Upper bound: the decreasing iteration W ← min(l_upper, max_a max_{s′} W(s′))
   from W = l_upper is run on the verifier's abstraction until nothing changes,
   and V_upper ≥ W must hold in every cell. Every iterate W^(k) bounds V_1 from
   above, so the cells with V_upper < 0 lie in the BRS. (A super-solution check
   would not suffice here: W ≡ −1 passes it.)
4. Closed-loop simulation with the exact flow from one uniformly random state in
   every cell with V_lower > 0 under π for 1000 steps: no entry into F, no exit
   from the grid, and l(x_t) ≥ V_lower(s0), i.e. dist(x_t, F) ≥ V_lower(s0), at
   every step.
5. Cells with V_upper < 0: from one uniformly random state in every such cell,
   each constant action, the greedy action on W and uniformly random actions must
   reach l < 0 or leave the grid within k̄(s0) steps, the first iteration with
   W^(k)(s0) < 0 (W^(k)(s0) < 0 means that every control sequence does so within
   k steps).

## Scope

The guarantees concern the sampled system x_{t+1} = Φ(x_t, u_t, τ) at the sampling
instants. Leaving the grid counts as failure (the solver's sink), so cells with
V_upper ≤ 0 (reach-avoid) or V_upper < 0 (avoid) are certified with respect to the
failure set extended by the region outside the grid.
