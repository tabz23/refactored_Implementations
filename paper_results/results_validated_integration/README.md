# Integration error checked with CAPD (paper, Appendix H.1)

Appendix H.1 of the paper (label `app:bv_setup`) states:

> In Section 9, f is computed numerically, and each ball of Section 4 is enlarged to cover the
> integration error, which we checked with CAPD [18].

This folder holds that check. **It does not change how the paper's results are produced:** the
reproduction scripts (`../reproduce_all.sh` and the scripts in the other folders) build
`cpp/build/dhj` without CAPD and never pass `--validated`, so they use the Dormand–Prince
integrator and its budget exactly as in the paper. The normal build of the current code
reproduces the stored per-cell values of the Section 9.2.4 run byte for byte (checked
2026-11-06).

## What is checked

The solver builds one transition per cell s and action a: it integrates the cell centre x_c(s)
over one sampling period t_s with the adaptive Dormand–Prince 5(4) method
(`cpp/include/dhj/core/ode.hpp`, tolerances 1e-12), which gives the numerical end point ŷ and
the budget e(s,a) (the tolerance of each accepted step, grown by exp(λh) over each step of length
h, where λ is the Lipschitz constant of the vector field). The reach box is ŷ ± (L_f r + e(s,a)),
where r is the cell radius, and the successors are the cells it meets (plus the out-of-domain
state if it leaves the state box).

`run_validated.sh` repeats each paper run with `build/dhj_capd --validated`. For every
transition of the run, in every refinement phase:

1. CAPD (Taylor order 15, automatic step control, CAPD's own interval arithmetic) encloses the
   exact solution from x_c(s) at t_s in a box [lo, hi] about 1e-15 wide. If the automatic step
   control fails, the call is retried with fixed steps t_s/1, ..., t_s/16 (also rigorous). This
   happens only for bicycle cells with steering angle exactly 0 under steering rate 0 and
   nonzero acceleration, whose solution is polynomial in time.
2. The validated error c(s,a) is the distance from ŷ to the far edge of [lo, hi] in the
   infinity norm, rounded up. So the exact end point lies within c(s,a) of ŷ.
3. The run uses the box ŷ ± (L_f r + c(s,a)), and also builds the box ŷ ± (L_f r + e(s,a)) of
   the normal build from the same ŷ and compares the successors of the two.

`compare_validated.py` then checks, for each run:

| check | meaning | needed to pass |
|---|---|---|
| CAPD failures | transitions for which CAPD gave no enclosure, even with fixed steps (their box keeps the budget) | 0 |
| validated error > budget | transitions with c(s,a) > e(s,a). If there are none, the box of the normal build (used in the paper) contains the validated box, which contains the exact image of the cell. This is the statement of Appendix H.1 | 0 |
| successor cells added / dropped | cells that only the validated box, or only the normal box, meets. Since the validated box lies inside the normal box, none can be added. A dropped cell is one that the normal box meets only within its last ~1e-11, so the paper's abstraction is slightly more conservative there | 0 added (dropped cells are reported) |
| same results as the paper run | `phases.csv` rows and the classification counts in `run.log`, and for Section 9.2.4 the per-cell values (`csv_check.txt`) | identical |

Not covered (these are not integration errors): the Lipschitz constant behind the Grönwall
factor L_f, and the rounding (about 1e-16) when the box edges are computed.

## Results

`summary.md` (written by `compare_validated.py`):

| case | paper | transitions | result |
|---|---|---|---|
| qv_ra | Sec. 9.2.4 reach-avoid, Fig. 3(b) | 5,184,000 | PASS |
| qv_avoid | Sec. 9.2.4 avoid | 5,184,000 | PASS |
| dubins120_std | Fig. 3(a) | 5,184,000 | PASS |
| avoid80 | Fig. 4 | 31,566,084 | PASS |
| dubins20 | Fig. 1 | 45,819,318 | PASS |
| vdp18 | Fig. 2(a) | 72,135,085 | PASS (4 transitions drop 12 cells, results identical) |
| bicycle8 | Fig. 2(d) | (running) | |
| evasion12 | Fig. 2(b) | (queued) | |
| di14 | Fig. 2(c) | (queued) | |

In every finished run, CAPD never failed and the validated error never exceeded the budget: the
largest validated error is about 1e-13 and the largest budget about 5e-11 (the validated error
is at most 0.85% of the budget). Every transition count equals the number of ODE solves in the
paper run's log. Fig. 3(b) is the Section 9.2.4 reach-avoid run without the value dump, so it has
the same transitions. Appendix H uses the closed-form Dubins flow and involves no numerical
integration.

## How to reproduce

1. Build CAPD with its native intervals (its default interval library, filib, needs x86):

   ```bash
   git clone https://github.com/CAPDGroup/CAPD && cd CAPD && mkdir build && cd build
   cmake .. -DCMAKE_BUILD_TYPE=Release -DCAPD_INTERVAL_TYPE=NATIVE && make -j8
   ```

   The runs here used CAPD 6.1.0 at commit 37057e5 on an Apple M3 Pro.

2. Run the cases and compare:

   ```bash
   CAPD_CONFIG=<CAPD>/build/bin/capd-config ./run_validated.sh            # all cases, or name some
   python3 compare_validated.py                                           # writes summary.md
   ```

   `run_validated.sh` builds `cpp/build/dhj_capd` (`make capd`), records `provenance.txt`, and
   writes `runs/` and `logs/` here. CAPD adds about 0.05–0.2 ms per transition (about 2–3 ms for
   the five-dimensional bicycle), so with 12 workers the Dubins 120³ runs take about 80 s each,
   Fig. 1 and Fig. 4 about 8–11 min, Fig. 2(a) about 22 min, and Fig. 2(b)–(d) about 1–2 hours
   each. Fig. 2(b) and 2(c) need 25–30 GB of memory, so macOS swaps to disk. The script stops a
   run if the free disk space drops below `MIN_FREE_GB` (default 8).

## Files

| file | content |
|---|---|
| `run_validated.sh` | runs the cases with `--validated` (same arguments as the paper runs, without plots) |
| `compare_validated.py` | checks each run, writes `summary.md` |
| `summary.md` | the results table |
| `csv_check.txt` | byte comparison of the Section 9.2.4 per-cell values with the stored paper runs |
| `logs/<case>.log` | solver output, including the "Validated integration (CAPD)" counts after every successor computation |
| `runs/<mode>/<run>_validated/` | `run.log`, `phases.csv`, `params.json` (the 400 MB value dumps are not kept, see `.gitignore`) |
| `provenance.txt` | code and CAPD commits, flags, compiler, binary hashes |

The code: `cpp/src/capd_bridge.cpp` (the only file that uses CAPD, compiled with CAPD's flags),
`cpp/include/dhj/core/validated.hpp`, `Dynamics::validated_errors` and `capd_spec` in
`cpp/include/dhj/dynamics/`, and `compare_with_budget_boxes` in `cpp/include/dhj/core/vi.hpp`.
See also "Validated integration with CAPD" in `cpp/README.md`.
