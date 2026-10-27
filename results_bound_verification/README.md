# Numerical check of the convergence theorem (Dubins car)

This directory checks the paper's convergence theorem (label `thm:convergence`):
if γ ∈ (0,1) with γ·L_f < 1, U is finite, A = U and all cells have radius η (written ε in
the paper), then for
every cell s and every x ∈ [[s]] the fixed points satisfy

  |V̄_γ(s) − V_γ(x)| ≤ 2Lη/(1−γ)  and  |V̲_γ(s) − V_γ(x)| ≤ 2Lη/(1−γ),  L = max{L_l, L_r},

and the outputs of Algorithm 1 with stopping tolerances δ̄\* and δ̲\* satisfy

  |V̄_γ^cons(s) − V_γ(x)| ≤ (2Lη + γδ̄\*)/(1−γ)  and  |V̲_γ^cons(s) − V_γ(x)| ≤ (2Lη + γδ̲\*)/(1−γ).

The check uses the second form, with δ̄\* = δ̲\* = δ\* = 10⁻¹².

A second script, `unclassified_check.py`, checks the theorem's last claim on the same 40
runs: every cell s that the classification theorem (`thm:soundness`) leaves unclassified,
i.e., V̲_γ^cons(s) ≤ 0 < V̄_γ^cons(s), satisfies |V_γ(x)| ≤ β := (2Lη + γδ\*)/(1−γ) for all
x ∈ [[s]]. The paper reports this check in Appendix H, subsection *Bound on unclassified
cells* (`app:bv_unclassified`). Details are in the section *Unclassified cells* below.

## Reproduce

```bash
./reproduce.sh            # replicates everything in ../results_bound_verification_rerun
                          # and compares the numbers with summary.json here
./reproduce.sh DIR        # the same, into DIR
./reproduce.sh "$PWD"     # regenerates the results in this directory
```

Needs Python 3 with numpy and matplotlib (recorded: Python 3.12.8, numpy 2.3.3,
matplotlib 3.10.6); the sweep takes about 3.5 minutes.

The check of the last claim is a separate script. Run it with the same Python as
`reproduce.sh`. It reads `bound_verification.py` and `summary.json` from this folder, writes
nothing, prints one JSON line per run, and ends with `ALL CHECKS PASS` or `SOME CHECK FAILED`:

```bash
python unclassified_check.py | tee DIR/unclassified_check.log   # about 3.5 minutes
python unclassified_check.py 20 50 60   # quick test: 20 sampled cells per run, N = 50 and 60
```

The arguments are `[SAMPLES_PER_RUN] [N ...]`. The defaults (200 sampled cells per run, all
eight N, seed 0) reproduce `unclassified_check.log`. To compare a rerun with it, ignoring
timings:

```bash
python - unclassified_check.log DIR/unclassified_check.log <<'EOF'
import json, sys
def rows(path):
    return [{k: v for k, v in json.loads(line).items() if k != "seconds"}
            for line in open(path) if line.startswith("{")]
ref, new = (rows(p) for p in sys.argv[1:3])
print("identical to the reference results" if ref == new else "DIFFERENT from the reference results")
EOF
```

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

### Unclassified cells (`unclassified_check.py`)

The theorem's last claim bounds |V_γ(x)| by β = (2Lη + γδ\*)/(1−γ) on every unclassified
cell. The script reruns Algorithm 1 with the functions of `bound_verification.py`, checks
that the iteration counts and G match `summary.json` (`same_as_appendix_H`), and bounds
|V_γ(x)| in two ways.

1. **All unclassified cells.** Soundness gives V̲_γ^cons(s) ≤ V_γ(x) ≤ V̄_γ^cons(s), so on an
   unclassified cell |V_γ(x)| ≤ M(s) := max{V̄_γ^cons(s), −V̲_γ^cons(s)}. The script takes
   the largest M(s) over all unclassified cells and reports its ratio to β (`sound_ratio`).
2. **Sampled states, exhaustive search.** This bound on V_γ does not depend on Algorithm 1,
   whose outputs only identify the unclassified cells.
   * *Sampling.* In each run: 200 unclassified cells drawn uniformly at random without
     replacement, plus the two unclassified cells with the largest V̄_γ^cons and the largest
     −V̲_γ^cons, and one state drawn uniformly at random in each of these 202 cells. One
     generator, `numpy.random.default_rng(0)`, serves all runs, in the order N = 50, …, 120
     and, for each N, γ = 0.3, …, 0.7.
   * *Search.* V_H(x) is the largest discounted reach-avoid measure g_γ(ξ, t) over all 3^H
     action sequences in {−1, 0, 1}^H and all t ≤ H. The trajectories use the closed-form
     flow (`flow` in `bound_verification.py`). Then
     V_H(x) ≤ V_γ(x) ≤ max{V_H(x), γ^(H+1)·r_max} with r_max = max r = 0.5, because the
     Bellman equation gives V_γ(x) ≥ g_γ(ξ, t) for every action sequence, and every time
     t > H contributes at most γ^t·r(x_t) ≤ γ^(H+1)·r_max. H is the smallest horizon with
     γ^(H+1)·r_max ≤ 0.01: H = 3, 4, 5, 7, 10 for γ = 0.3, 0.4, 0.5, 0.6, 0.7.
   * *Checks.* V_H(x) ≥ −β and max{V_H(x), γ^(H+1)·r_max} ≤ β, which give |V_γ(x)| ≤ β
     (`direct_ok`, `direct_ratio`). In addition, the interval
     [V_H(x), max{V_H(x), γ^(H+1)·r_max}] must overlap [V̲_γ^cons(s), V̄_γ^cons(s)], because
     both contain V_γ(x) (`soundness_violations`, `overlap_min_slack`).

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
| `unclassified_check.py` | the check of the theorem's last claim on unclassified cells (imports `bound_verification.py`) |
| `unclassified_check.log` | its output for the paper: one JSON line per run, a summary line, and `ALL CHECKS PASS` |

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
* Unclassified cells (`unclassified_check.log`). The recomputed runs match `summary.json`.
  * Over all unclassified cells (14,756 to 112,996 per run), the largest M(s) is 0.58 to
    0.71 of β, so the last claim holds at every state of every unclassified cell.
  * At all 8,080 sampled states, the exhaustive-search interval lies in [−β, β]. Per run,
    the largest |V_H(x)| is 0.30 to 0.70 of β, and γ^(H+1)·r_max is at most 0.053·β.
  * The interval overlaps [V̲_γ^cons(s), V̄_γ^cons(s)] at every state, with the smallest
    overlap 0.0034.
