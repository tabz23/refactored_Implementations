# Reproducing the results in the paper

Paper: *A Quantitative Approach to Reach–Avoid Reachability Analysis for Discrete-Time
Nonlinear Control Systems*.

Three folders hold every result in the paper that this code produces. Each one has a script
that regenerates it from scratch. Result numbers below follow the paper as of 2026-10-26.
LaTeX labels are given too, because numbering can change.

| Folder | Paper results | Script | Time (Apple M3 Pro) |
|---|---|---|---|
| `results_hscc_2027_newbellmans_12workers_current_code/` | Figs. 2–5, the runtimes in the Fig. 3–5 captions and in §9.2.3 | `regenerate_hscc_newbellmans_tol1e-15_current_code.sh` | about 1 hour |
| `results_quantitative_verification/` | §9.2.5, the closed-loop validation of Proposition 3 (`prop:greedy`) and Proposition 4 (`prop:avoid`) | `reproduce.sh` | about 3 minutes |
| `results_bound_verification/` | Appendix H (`app:bound_verification`): the numerical check of Theorem 4 (`thm:convergence`), its figure, and the check of its last claim on unclassified cells (`app:bv_unclassified`) | `reproduce.sh` and `unclassified_check.py` | about 3.5 minutes each |

The other paper results do not come from this code. Fig. 1 and Fig. 4(e) come from the
MATLAB HJ reachability toolbox. Fig. 4(c,d) come from ROCS 2.0.

## Requirements

* **Solver:** a C++17 compiler and zlib. The recorded compiler is Apple clang 16.0.0. The
  scripts compile the solver with `-march=native`, so build it on the machine that runs it.
  The scripts do this for you.
* **Python:** Python 3 with numpy, scipy and matplotlib. The recorded versions are Python
  3.12.8, numpy 2.3.3, scipy 1.16.2 and matplotlib 3.10.6. The scripts use
  `~/miniconda3/bin/python` by default; set `PYTHON` to use another interpreter. LaTeX with
  dvipng is optional and only typesets figure text.
* **Hardware for the reported runtimes:** a MacBook with an Apple M3 Pro (12 cores) and 18 GB
  of RAM, with 12 worker threads.
  * Runtimes vary by up to about ±20% between runs on the same machine.
  * Classifications are deterministic.

## 1. Figures 2–5 and their runtimes

Folder: `results_hscc_2027_newbellmans_12workers_current_code/`

### Run

The script writes its results into the folder that contains it, and it builds the solver
from `../cpp`. Keep the script in a folder directly inside `refactored_Implementations/`. It
can be started from any directory.

**To reproduce the paper's results without touching the stored ones.** This is the
recommended way, including for artifact evaluation. Run a copy of the script in a new folder
next to the stored one:

```bash
cd refactored_Implementations
mkdir rerun_paper_results
cp results_hscc_2027_newbellmans_12workers_current_code/regenerate_hscc_newbellmans_tol1e-15_current_code.sh rerun_paper_results/
cd rerun_paper_results
nohup ./regenerate_hscc_newbellmans_tol1e-15_current_code.sh > regenerate.log 2>&1 &
tail -f regenerate.log          # each case prints START/END; the last line is "ALL DONE"
```

Then compare the new results with the stored ones, from `refactored_Implementations/`. Every
case should print "same", because classifications are deterministic. Runtimes depend on the
machine.

```bash
NEW=rerun_paper_results
REF=results_hscc_2027_newbellmans_12workers_current_code
# Stored folders carry "_eps<--epsilon>" (the option they were made with), new ones
# "_epsmin<--eps-min>"; runs are paired on the rest of the folder name.
key() { sed -E 's/_eps(min)?[^_]*//' <<< "$1"; }
for d in "$NEW"/{ra_nodiscount,ra_discount,avoid_nodiscount}/*/; do
  run="${d#"$NEW"/}"; run="${run%/}"; ref=""
  for r in "$REF"/{ra_nodiscount,ra_discount,avoid_nodiscount}/*/; do
    r="${r#"$REF"/}"; r="${r%/}"
    [ "$(key "$r")" = "$(key "$run")" ] && ref="$r"
  done
  if [ -n "$ref" ] && diff <(grep -A3 "Classification after phase" "$d/run.log") \
          <(grep -A3 "Classification after phase" "$REF/$ref/run.log") > /dev/null; then
    echo "same classification:      $run"
  else
    echo "DIFFERENT classification: $run"
  fi
done
grep -H "Total time" "$NEW"/*/*/run.log
```

The comparison is only meaningful with the same solver code. The new `provenance.txt` must
show no uncommitted solver changes. `git diff 8751e29 <new commit> -- cpp` may show only the
change to the refinement threshold described under "Changes after the paper run" at the end.
That change does not alter any result.

**To regenerate the stored results in place,** run the script inside
`results_hscc_2027_newbellmans_12workers_current_code/`. The outputs already there move to
`previous_<timestamp>/`.

The script does four things:

1. It compiles `cpp/src/main.cpp` into `cpp/build/dhj`, with the flags of `cpp/Makefile`, so
   every case runs the current solver code.
2. It writes `provenance.txt`, which records the commit, any uncommitted solver changes, the
   compiler and the binary hash. If the solver has uncommitted changes, it also writes them
   to `solver_changes.diff`.
3. It moves outputs of earlier runs in this folder to `previous_<timestamp>/`, so it never
   overwrites them.
4. It runs the eight cases below. Each uses 12 workers, writes PNG figures at 300 dpi, and
   writes no checkpoints.

### Cases

| Script name | Mode | Setup | Paper |
|---|---|---|---|
| `dubins120_std` | `ra_nodiscount` (γ = 1, Algorithm 1 with the modifications of Proposition 2, `prop:gamma_1`) | Dubins car, 120³ grid (cell radius π/120 ≈ 0.026), t_s = 0.3, no refinement | Fig. 4(a), §9.2.3 |
| `dubins120_discount` | `ra_discount` (γ = 0.96, δ\* = 1e-8) | the same grid | Fig. 4(b), §9.2.3 |
| `avoid80` | `avoid_nodiscount` (γ = 1, Section 8) | Dubins car, 80³ grid (cell radius π/80 ≈ 0.039), Algorithm 2, 6 refinements | Fig. 5 |
| `dubins20_discount` | `ra_discount` (γ = 0.96, δ\* = 1e-8) | Dubins car, 20³ grid (cell radius π/20 ≈ 0.157), Algorithm 2, 13 refinements | Fig. 2 |
| `vdp18` | `ra_nodiscount` | Van der Pol, 40² grid (cell radius 0.0725), t_s = 0.01, 18 refinements | Fig. 3(a) |
| `evasion12` | `ra_nodiscount` | 3D evasion, 30³ grid (cell radius π/30 ≈ 0.105), t_s = 0.3, 12 refinements | Fig. 3(b) |
| `di14` | `ra_nodiscount` | planar double integrator, 10⁴ grid (cell radius 0.3), t_s = 0.3, 14 refinements | Fig. 3(c) |
| `bicycle8` | `ra_nodiscount` | kinematic bicycle, 5⁵ grid (cell radius π/5 ≈ 0.628), t_s = 0.3, 8 refinements | Fig. 3(d) |

The cell radius is the paper's ε (Section 4): half the largest side of a cell over all
dimensions, so the heading θ ∈ [−π, π] counts. For the Dubins car and the evasion model, the
position is in [−3, 3]², so an N³ grid has half-sides 3/N, 3/N and π/N, and radius π/N.

Notes on the arguments:

* **`--gamma`:** it only affects `ra_discount`. In `ra_nodiscount` and `avoid_nodiscount`
  the solver uses γ = 1, and `--gamma` appears only in the folder name.
* **`--eps-min`:** the threshold ε_min of Algorithm 2. A boundary cell is refined only while
  its radius exceeds ε_min (`filter_refinable` in `refine.hpp`). The paper does not report it.
  * The script uses 1e-6 for every case. That is far below the smallest radius any case
    refines, 2.4e-4 in `vdp18`'s last iteration. So every case ends at its `--refinements`
    limit.
  * The stored outputs were made before this option existed, with `--epsilon` 0.05, 0.1 or
    0.001, which set ε_min = `--epsilon`/(4 L_l). In none of those runs was a cell skipped for
    being too small: `Below threshold` is 0 in every refinement phase. So they are exactly
    what `--eps-min 1e-6` produces.
  * The only visible difference is in the stored `vdp18` and `evasion12` runs. There, the cells
    left unclassified after the last iteration are below the old threshold, so `Final Boundary
    State` shows `Still refinable: 0`, and the last row of `phases.csv` has `refinable_cells`
    0. With `--eps-min 1e-6` both show the full count. This changes no result.
  * Folder names contain `eps<--epsilon>` in the stored runs and `epsmin<--eps-min>` in new
    ones.
* **Stopping tolerances:**
  * Discounted runs stop when both sup-norm changes are at most δ\* = 1e-8, as in the paper.
  * The γ = 1 runs use 1e-15. The paper says they stop when the change is 0. The changes that
    remain at that point are rounding effects of about 1e-16.
  * Both bounds stop on `--delta-min`, so it plays the role of both δ̲\* and δ̄\* of the paper.
  * `--delta-max` is used only by the "conservative" mode, which no case here uses
    (`Conservative mode: False` in every `run.log`). The scripts set it equal to
    `--delta-min`, and it appears only in the folder name.
  * Value iteration met the tolerance in every phase of every case; `--vi-iterations` was
    never reached.
* **`--dt` and `--tau`:** `--tau` is the sampling time t_s, and `--dt` the spacing of the
  intermediate reach checkpoints. Every case sets them equal, so each successor set is computed
  only at time t_s, as in Section 4 of the paper.
* **The bicycle:** it stops after 8 refinements, which is the output in Fig. 3(d). A run with
  13 refinements (`../results_bicycle_oom_test/`, 2026-10-26, same flags otherwise) split to
  1,435,125 cells in refinement 9 and exceeded the machine's 18 GB of RAM. Its log stops during
  refinement 9, and `time.txt` reports a peak memory footprint of 57.3 GB. §9.2.2 of the paper
  reports this.

### Where the paper's numbers come from

Each case writes `<mode>/<case folder>/` containing these files:

* **`run.log`:** the full log.
  * The runtime the paper reports is its last `Total time: … seconds` line.
  * The classified percentage is in the last `Classification after phase N` block: Safe plus
    Unsafe, as "% of state space".
* **`phases.csv`:** one row per phase, with cell counts, iterations and timings.
* **`value_function_phase_<N>_complete.png`:** the classification map after phase N. Green
  means $\underline V > 0$, red means $\overline V \le 0$, and gray means unclassified. The
  solver uses this test in every mode. For the avoid case, see known difference 4 at the
  end.
* **`params.json`:** the exact command and its arguments.

Quick extraction:

```bash
grep -H "Total time" */*/run.log
grep -A3 "Classification after phase" ra_nodiscount/van_der_pol*/run.log | tail -4
```

Values from the run used in the paper (2026-10-26, solver commit `8751e29`):

| Paper | Case | Total time | Classified (safe + unsafe) |
|---|---|---|---|
| Fig. 4(a), §9.2.3 | `dubins120_std` | 2.69 s | 83.4% |
| Fig. 4(b), §9.2.3 | `dubins120_discount` | 2.86 s | 83.4% |
| Fig. 5 | `avoid80` | 22.13 s (22.1 s in the caption) | 94.4% |
| Fig. 2 | `dubins20_discount` | 46.87 s | 94.5% |
| Fig. 3(a) | `vdp18` | 291.75 s | 97.7% |
| Fig. 3(b) | `evasion12` | 1301.23 s | 88.0% |
| Fig. 3(c) | `di14` | 1910.33 s | 94.3% |
| Fig. 3(d) | `bicycle8` | 63.74 s | 11.5% |

**Fig. 2b caveat:** it plots per-iteration runtimes, cell counts and the classified volume. It
was drawn with `python/plot_runtime.py` from an earlier run of `dubins20_discount`.
* The cell counts and classified volumes are identical to this run's.
* The per-iteration runtimes agree within about 10%. This run takes 0.52 s for iteration 1 and
  14.95 s for iteration 14 (`time_phase_s` in `phases.csv`).
* `python plot_runtime.py --run-dir <this run's folder>` redraws the figure from this run, into
  that folder. The paper's version labels the right axis "classified as in L/R" instead of
  "safe/unsafe".

**Comparison with the previous run:** `../results_hscc_2027_newbellmans_12workers/` was made
with older builds of the solver.
* Classifications are identical in all 79 phases of the eight cases.
* 78 of the 79 PNGs are pixel-identical. The remaining one, Fig. 4(b), differs only in its
  title line, because the old run used the "conservative" stopping mode.
* Only the runtimes differ.

## 2. Quantitative guarantees (§9.2.5)

Folder: `results_quantitative_verification/`. Its `README.md` has the details.

### Run

```bash
cd results_quantitative_verification
./reproduce.sh            # writes ../results_quantitative_verification_rerun and compares
                          # every number except timings with the stored summary.json
./reproduce.sh "$PWD"     # regenerates the results in this folder instead
```

### What it checks

1. **Solver runs.** `run_alg1_res120.sh` compiles the current solver into
   `cpp/build/dhj_head` and records `provenance.txt`. It then runs Algorithm 1 without
   refinement on the Dubins car with the 120³ grid:
   * reach–avoid with γ = 0.96 and δ\* = 1e-8;
   * avoid with γ = 1 and tolerance 1e-300, so it stops only when nothing changes, as in the
     paper.

   Both runs dump their value tables.
2. **Independent abstraction.** `verify_quantitative.py` reads only the dumped grid, the
   values, and the bounds on l and r. It rebuilds the abstraction itself:
   * the exact closed-form Dubins flow of each cell center;
   * the Grönwall reach box;
   * the successor cells and the sink outside the grid.
3. **Reach–avoid checks:**
   * the cell bounds;
   * the certificates $\hat T\underline V \ge \underline V$ and
     $\overline T\,\overline V \le \overline V$;
   * closed-loop simulation with the exact flow from one uniformly random initial state in
     each cell with $\underline V^{\mathrm{cons}} > 0$, under the policy $\hat\pi$ of
     Proposition 3. Each trajectory must reach the target within
     $\bar t = \lceil\log(\underline V(s_0)/r^*)/\log\gamma\rceil$ steps with
     $r^* = \max_s \underline r(s)$, never enter F or leave the grid, and keep
     $\gamma^\tau l(x_\tau) \ge \underline V(s_0)$ until it arrives.
4. **Avoid checks:**
   * the certificate $\underline V \le \min\{\underline l, \max_a \min_{s'} \underline V\}$;
   * closed-loop simulation for 1000 steps under the policy $\pi^*$ of Proposition 4. Each
     trajectory must keep $l(x_\tau) \ge \underline V_1(s_0)$.

### Results (stored `summary.json`, the numbers in §9.2.5)

* **Reach–avoid:** 71 iterations and 983,988 trajectories, with 0 violations.
  * The median arrival time is 19 steps against a median bound of 37.
  * $t \le 0.72\,\bar t$ for all trajectories (maximum 0.714). There are 76 trajectories with
    $\bar t = 0$, and all of them start inside T.
  * $\gamma^\tau l(x_\tau)/\underline V(s_0) \ge 1.02$.
* **Avoid:** 89 iterations and 973,616 trajectories, with 0 violations.
  * $l(x_\tau)/\underline V_1(s_0) \ge 1.002$ over all 1000 steps.

### Alignment with the paper

* **The policies are the propositions' definitions:**
  * $\hat\pi(s) \in \arg\max_a \min\{\underline l(s), \max\{\underline r(s), \gamma\min_{s'} \underline V(s')\}\}$ (Proposition 3);
  * $\pi^*(s) \in \arg\max_a \min_{s'} \underline V_1(s')$ (Proposition 4).
* **$t^*$ follows the proof:** it is the first τ with $\underline V(s_\tau) \le \underline r(s_\tau)$.
* **Naming in the outputs:** the summary key `distance_violations` counts trajectories with
  $\min_\tau l(x_\tau) < \underline V_1(s_0)$. It keeps its old name so that reruns still compare
  against the stored `summary.json`.

## 3. Convergence bound (Appendix H, Theorem 4)

Folder: `results_bound_verification/`. Its `README.md` has the details.

### Run

```bash
cd results_bound_verification
./reproduce.sh            # writes ../results_bound_verification_rerun and compares
                          # every number except timings with the stored summary.json
./reproduce.sh "$PWD"     # regenerates the results in this folder instead
python unclassified_check.py | tee DIR/unclassified_check.log
                          # the check of the last claim (unclassified cells); compare it
                          # with unclassified_check.log as described in the folder README
```

### What it checks

* **Implementation:** `bound_verification.py` is an independent numpy implementation of
  Algorithm 1 and does not use the C++ solver.
* **System:** a periodic Dubins car on the torus [−3, 3)² × [−π, π), so that f maps X into
  X, as Theorem 4 assumes.
* **Sweep:** γ ∈ {0.3, 0.4, 0.5, 0.6, 0.7} and N ∈ {50, 60, …, 120}.
  * Every γ satisfies γL_f < 1 with L_f = e^0.3.
  * Each grid has N³ cells with radius ε = π/N.
  * δ\* = 1e-12 and L = √2.
* **Quantity checked:** $G = \max_s(\overline V^{\mathrm{cons}}(s) - \underline V^{\mathrm{cons}}(s))$,
  against the bound $(2L\varepsilon + \gamma\delta^*)/(1-\gamma)$.
  * Theorem 2 gives $\underline V^{\mathrm{cons}} \le V_\gamma \le \overline V^{\mathrm{cons}}$,
    so both errors bounded in Theorem 4 are at most G.
  * G at or below the bound therefore confirms the theorem's inequalities for that run.
* **Unclassified cells:** `unclassified_check.py` checks Theorem 4's last claim,
  $|V_\gamma(x)| \le (2L\varepsilon + \gamma\delta^*)/(1-\gamma)$ on every cell with
  $\underline V^{\mathrm{cons}} \le 0 < \overline V^{\mathrm{cons}}$, in two ways.
  * Through Theorem 2, on all unclassified cells.
  * At 8,080 sampled states, by exhaustive search over action sequences with the
    closed-form flow, which bounds $V_\gamma$ without Algorithm 1.

### Results (stored `summary.json`)

* All 40 runs satisfy the bound. They take 24 to 41 iterations, and
  $(1-\gamma)G/(2L\varepsilon+\gamma\delta^*) \in [0.637, 0.727]$.
* The paper's figure is `bound_verification_paper.pdf`, copied to the paper repository as
  `new_figures/bound_verification_paper.pdf`.
* Unclassified cells (stored `unclassified_check.log`): all checks pass. The ratio to the
  bound is 0.58 to 0.71 over all unclassified cells, through Theorem 2. It is 0.30 to 0.70
  at the sampled states, by exhaustive search.

### Alignment with the paper

The implementation follows the paper exactly:
* Algorithm 1: start at $\min\{\underline l, \underline r\}$ and $\overline l$, take the min
  or max over successors, then the max over actions, and stop on the residual.
* The abstraction: a ball of radius $L_f\varepsilon$ around the exact image of each cell
  center, with closed cell boxes.
* The cell bounds: $\pm L\varepsilon$ around the center values.

This folder is not yet committed to git.

## Code version and alignment with the paper

### Which code produced which results

* Section 1 used solver commit `8751e29`, and Section 2 used `d8463b7`.
* `git diff d8463b7 8751e29 -- cpp` is empty, so both used the same solver code.
* Each folder's `provenance.txt` records the commit. To check that a rerun uses the same code,
  compare it with `git rev-parse HEAD` and `git status --short -- cpp`.

### Where the paper's algorithm lives in the solver (`cpp/include/dhj/core/`)

| Paper | Solver |
|---|---|
| Algorithm 1, discounted | `ra_discount`: $\underline V^{(0)} = \min\{\underline l, \underline r\}$, $\overline V^{(0)} = \overline l$; the operators of Theorem 1; stops when both sup-norm changes are at most δ\* (`vi.hpp`) |
| Proposition 2 (γ = 1, reach–avoid) | `ra_nodiscount`: both bounds start at $\min\{l, r\}$; γ is not applied |
| Section 8 (avoid) | `avoid_nodiscount`: both bounds start at $l$; the operators $\mathcal T^{\mathrm{av}}$ |
| Section 4 (abstraction) | reach box of radius $L_f$ × cell radius around the integrated image of the center, with $L_f = e^{\lambda t_s}$ (Grönwall); successors are the cells meeting the box; outside the grid is one absorbing sink with value −1 |
| Section 5 (cell bounds) | $l(x_c) \pm L_l\,$radius and $r(x_c) \pm L_r\,$radius, with $L_l = L_r = \sqrt2$ |
| Algorithm 2 | refines boundary cells ($\overline V > 0 \ge \underline V$) whose radius $\max_j \eta_{x,j}(s)$ exceeds ε_min (`--eps-min`), bisects them along the widest dimension, then re-initializes and recomputes every cell; stops when no cell qualifies or after `--refinements` iterations (`refine.hpp`) |
| §9.2 action sets | Dubins {−1, 0, 1}; evasion 5 evenly spaced values in [−1, 1]; double integrator {−1, 0, 1}²; bicycle {−0.2, 0, 0.2} × {−1, 0, 1} |

### Known differences from the paper's text

None of these is stated in the paper.

1. **Stopping at 1e-15:** the γ = 1 runs stop when the change is below 1e-15 instead of
   exactly 0. See Section 1.
2. **Integration error:** the solver enlarges each reach box by the Runge–Kutta integrator's
   tolerance (about 1e-12) to cover integration error. This makes the bounds slightly more
   conservative.
3. **Lipschitz constants:** the vector-field constants of Van der Pol (λ = 20.75) and the
   bicycle (λ = 2.152) hold on the state box. A trajectory from an edge cell can leave the
   box within a step, where the valid constants are about 21.38 and 2.28.
4. **Unsafe test in the avoid case:** the solver labels a cell unsafe when
   $\overline V \le 0$. The paper's avoid certificate (Proposition 4(2), Fig. 5 caption) uses
   $\overline V < 0$. The two agree unless a cell has $\overline V = 0$ exactly. Upper values
   come only from the cell bounds $\overline l(s)$ and the sink value −1, so this does not
   happen in practice; the §9.2.5 avoid run has no such cell (`cells_with_V_upper_zero` = 0).
5. **Algorithm 2 stops at an iteration limit:** the paper's Algorithm 2 runs until its queue is
   empty. With `--eps-min 1e-6`, every run stops after `--refinements` iterations while the
   queue is still nonempty, and the figures show that output.
   * Theorem 3 applies after every iteration, because each iteration ends with Algorithm 1 on
     the current grid.
   * The captions and §9.2.2 give the iteration counts for Fig. 5 (six) and the bicycle
     (eight).

### Changes after the paper run

* **2026-10-26, refinement threshold (`args.hpp`, `refine.hpp`, `cell.hpp`, `run.hpp`,
  `params.hpp`, `checkpoint.hpp`):**
  * **Before:** `--epsilon` set the threshold, and a cell was refined while its largest side
    exceeded `--epsilon`/(2 L_l). The paper does not use this formula.
  * **Now:** `--eps-min` gives ε_min directly. A cell is refined while its radius
    (`Cell::max_radius()`, half its largest side) exceeds ε_min, which is the condition
    $\max_j \eta_{x,j}(s) > \varepsilon_{\min}$ of Algorithm 2.
    * `--epsilon` no longer exists, and passing it is an error.
    * `params.json` records `eps_min`, and folder names contain `epsmin`.
    * The log prints `Refinable (radius > eps_min=…)`.
  * **Scripts:** they now pass `--eps-min 1e-6`, which never stops a refinement in these runs;
    see the `--eps-min` note in Section 1. Results are unchanged because no cell in the stored
    runs was ever skipped for size.
  * **Check:** rerunning `vdp18` with the new code and `--eps-min 1e-6` reproduced the stored
    run.
    * All 19 classification blocks are identical, and all 19 PNGs are byte-identical.
    * `phases.csv` is identical except `refinable_cells` in the last row, which is 9,539,610
      instead of 0, as explained in the `--eps-min` note.
    * The runtime was 317.82 s, against 291.75 s stored.
  * **Smoke test of the other cases** (same flags as the script):
    * `dubins120_std`, `dubins120_discount`, `avoid80`, `dubins20_discount` and `bicycle8`
      ran in full. The first 10 refinements of `di14` and the first 9 of `evasion12` ran too.
    * Every compared phase has identical `phases.csv` counts and classifications, and
      byte-identical PNGs.
    * `results_quantitative_verification/reproduce.sh` reported both summaries identical to
      the stored ones.
    * The logs differ only in the threshold lines, the radius statistics (half the old side
      statistics), output paths and timings.
