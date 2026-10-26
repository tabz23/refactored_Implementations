# Numerical check of the convergence theorem (Dubins car)

This directory checks the paper's convergence theorem (label `thm:convergence`):
if γ ∈ (0,1) with γ·L_f < 1, U is finite, A = U and all cells have radius η (written ε in
the paper), then for
every cell s and every x ∈ [[s]] the fixed points satisfy

  |V̄_γ(s) − V_γ(x)| ≤ 2Lη/(1−γ)  and  |V̲_γ(s) − V_γ(x)| ≤ 2Lη/(1−γ),  L = max{L_l, L_r},

and the outputs of Algorithm 1 with stopping tolerances δ̄\* and δ̲\* satisfy

  |V̄_γ^cons(s) − V_γ(x)| ≤ (2Lη + γδ̄\*)/(1−γ)  and  |V̲_γ^cons(s) − V_γ(x)| ≤ (2Lη + γδ̲\*)/(1−γ).

The check uses the second form, with δ̄\* = δ̲\* = δ\* = 10⁻¹².

## Reproduce

```bash
./reproduce.sh            # replicates everything in ../results_bound_verification_rerun
                          # and compares the numbers with summary.json here
./reproduce.sh DIR        # the same, into DIR
./reproduce.sh "$PWD"     # regenerates the results in this directory
```

Needs Python 3 with numpy and matplotlib (recorded: Python 3.12.8, numpy 2.3.3,
matplotlib 3.10.6); the sweep takes about 3.5 minutes.

## Setup

* **System.** The Dubins car of the paper (v = 1, actions {−1, 0, 1} = U, sampling time
  t_s = 0.3, obstacle disc of radius 1.3 at the origin, target disc of radius 0.5 at
  (2.5, 0)), on a periodic domain: x₁ and x₂ wrap around [−3, 3) and the heading wraps
  around [−π, π). The sampled dynamics then map the domain into itself, as the paper
  assumes (f : 𝒳 × 𝒰 → 𝒳). In the paper's own experiments, leaving the grid counts as
  failure, which makes V_γ discontinuous at the grid boundary, so the theorem does not
  apply to that setting.
* **Lipschitz constants.** l and r are the signed periodic distances to the obstacle and
  the target, so L_l = L_r = √2 in the ∞-norm. The vector field is 1-Lipschitz in the
  ∞-norm, so the sampled map is L_f = e^0.3 ≈ 1.350-Lipschitz (Grönwall), the constant
  the solver also uses.
* **Discount factors.** γ ∈ {0.3, 0.4, 0.5, 0.6, 0.7}, so γ·L_f ∈ {0.40, …, 0.94} < 1.
* **Grids.** Uniform grids of N³ cells, N ∈ {50, 60, …, 120}. The cell radius (largest
  half-width, ∞-norm) is η = π/N.
* **Algorithm 1.** Cell bounds l̲, l̄, r̲, r̄ = center value ∓ Lη; the successors of a cell
  under an action are the cells whose closed boxes meet the ball of radius L_f·η around
  the exact image of the cell center (closed-form flow); value iteration from
  V̲ = min(l̲, r̲) and V̄ = l̄ runs until both sup-norm changes are at most δ\* = 10⁻¹², so
  the computed bounds are within γδ\*/(1−γ) ≤ 2.4·10⁻¹² of the fixed points.

## What is checked

The theorem bounds |V̄_γ^cons(s) − V_γ(x)| and |V̲_γ^cons(s) − V_γ(x)|, which cannot be
computed because V_γ is unknown. The bounds are sound, V̲_γ^cons(s) ≤ V_γ(x) ≤ V̄_γ^cons(s)
for every x ∈ [[s]], so both errors are at most the gap V̄_γ^cons(s) − V̲_γ^cons(s). Hence

  max_s (V̄_γ^cons(s) − V̲_γ^cons(s)) ≤ (2Lη + γδ\*)/(1−γ)

implies both inequalities of the theorem for every state of every cell, without knowing
V_γ. This check is stronger than the theorem requires: the theorem only implies that the
gap is at most twice the bound. For every (γ, N) the script reports this largest gap, its ratio to the bound
(`ratio`, at most 1 means verified), the mean gap, the fractions of cells certified in 𝓡
(V̲ > 0), in 𝓛 (V̄ ≤ 0) and unclassified, and, per γ, the least-squares slope of
log(largest gap) against log η (1 for a linear rate).

## Files

| file | content |
|---|---|
| `bound_verification.py` | the check (self-contained) |
| `reproduce.sh` | runs it and compares the numbers |
| `summary.json` | every run and the per-γ summary |
| `table_bounds.tex` | booktabs table: one row per γ |
| `bound_verification_paper.{pdf,png}` | single-column figure for the paper (Appendix H, copied to the paper's `new_figures/`): max_s (V̄_γ^cons(s) − V̲_γ^cons(s)) of Algorithm 1 against η (solid) with the bound (2Lη + γδ\*)/(1−γ) (dashed), one color and marker per γ, LaTeX labels |
| `bound_verification.{pdf,png}` | largest gap against η with the bound, and their ratio |
| `value_functions/gamma<γ>_N<N>.png` | for every run: the lower bound, the upper bound and the classification computed by Algorithm 1, in six heading slices (0, π, ±π/4, ±π/2); no value tables are saved |
| `run.log` | console output of the last run |

## Results

All 40 pairs (γ, N) satisfy max_s (V̄_γ^cons(s) − V̲_γ^cons(s)) ≤ (2Lη + γδ\*)/(1−γ), so
the theorem holds for every state of every cell in every run. With δ\* = 10⁻¹² the term
γδ\*/(1−γ) is below 2.4·10⁻¹², so the bound is 2Lη/(1−γ) to all printed digits. Per discount factor (`table_bounds.tex`):

| γ | γ·L_f | largest gap / bound, max over N | rate of the largest gap | unclassified cells, N = 50 → 120 |
|---|---|---|---|---|
| 0.3 | 0.40 | 0.70 | 1.00 | 50.3% → 6.5% |
| 0.4 | 0.54 | 0.70 | 1.00 | 27.9% → 5.0% |
| 0.5 | 0.67 | 0.73 | 1.00 | 11.8% → 5.0% |
| 0.6 | 0.81 | 0.70 | 1.00 | 11.8% → 5.0% |
| 0.7 | 0.94 | 0.71 | 0.94 | 11.8% → 5.0% |

* The largest gap decreases linearly in η (fitted slopes 0.94 to 1.00), as the theorem
  predicts, and its ratio to the bound stays between 0.64 and 0.73 for every γ and N, so
  the factor 1/(1−γ) captures the dependence on γ.
* The curves for γ ≤ 0.6 are straight lines because G = max_s (V̄_γ^cons(s) − V̲_γ^cons(s)) is an
  exact multiple of η there. For γ = 0.3, G = 2Lη in every run, attained e.g. inside the obstacle,
  where V̄ = l̄ and V̲ = l̲. For γ ∈ {0.4, 0.5, 0.6}, G = γ(2L + 12√2/π)η in 18 of the 24 runs and
  less than 0.34% below it in the other 6, attained at cells whose bounds are γ·l̄(s') and γ·l̲(s'')
  for two successors whose centers are two cells apart in x₁ and x₂ on a ray from the obstacle
  center (l is affine along it, so l̄(s') − l̲(s'') = 2Lη + 12√2/N). For γ = 0.7, at least one bound
  at the maximizing cells is set by l or r two or three steps away, G/η varies by 10% across N
  (6.68 for N ∈ {60, 80, 100, 120}, 6.00 to 6.64 otherwise), and G is slightly larger for N = 120
  than for N = 110 (the grids are not nested).
* The unclassified cells are contained in {x : |V_γ(x)| ≤ 2Lη/(1−γ)}, as the theorem's last
  statement says, so for every κ > 0 they contain no state with |V_γ(x)| ≥ κ once 2Lη/(1−γ) < κ.
