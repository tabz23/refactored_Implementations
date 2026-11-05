# Stay-in-grid runs: reach sets restricted to the grid

Exploratory reruns of the eight cases of Figs. 2–5, made on 2026-11-04 to see what the
paper's figures would be if f mapped X into X. **Not sound for our systems, which can leave
X; the paper does not use them.**

**What changes.** The solver flag `--restrict-reach-set-to-grid` ignores the part of each
reach box that lies outside the state box and adds no sink state. An action whose reach box
lies entirely outside the box then has no successor, and the solver counts it as −∞, so
cells whose every action leaves the grid still come out red.

**How they were generated.**

```bash
./run_stay_in_grid.sh              # all eight cases, each under 5 minutes (about 8 min total)
FULL=1 ./run_stay_in_grid.sh       # the paper's refinement counts, no time limit
python compare_with_paper.py > comparison.txt
python make_comparison_images.py   # comparison_images/, numpy and Pillow needed
```

The script compiles `cpp/build/dhj` from the current source and runs each case with the
paper's flags plus `--restrict-reach-set-to-grid` (12 workers, PNGs at 300 dpi). To keep
every case under 5 minutes, `vdp18`, `di14` and `evasion12` stop after 16, 10 and 9
refinements instead of 18, 14 and 12. `provenance.txt` records the commit. Run folders end
in `_restrictgrid`.

**Results** (`comparison.txt`, shares of the state space, paper run → this run):

* Only the certified-safe share grows, by up to 8 points: 56.9 → 60.2% for Fig. 4(a,b),
  63.9 → 66.9% for Fig. 5, 64.3 → 67.3% for Fig. 2, 71.3 → 79.0% for the double integrator
  (36% of its reach boxes cross the velocity limits). The bicycle does not change.
* On the starting grid of every case the unsafe share is identical; later phases differ by
  up to 1.5 points because the two runs refine different cells.
* The changed cells are thin strips along the inner edge of the red bands
  (`comparison_images/`: paper run on top, this run in the middle, changed cells below).
