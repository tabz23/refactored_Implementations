# Stay-in-grid runs with cell-wise bounds for actions that leave

A variant of `../results_stay_in_grid`, made on 2026-11-04. **Not sound for our systems;
the paper does not use it.**

**What changes.** As in `../results_stay_in_grid`, the part of each reach box outside the
state box is ignored and no sink state is added (`--restrict-reach-set-to-grid`). In
addition, `--if-leave-use-cellwise-bounds` lets an action whose reach box lies entirely
outside the box (no successor) use only the cell's own bounds, as if the trajectory ended
there: min{l, r} in reach–avoid, which is what `../results_stay_in_grid` already did, and l
in the avoid case, where `../results_stay_in_grid` used −∞.

**How they were generated.**

```bash
./run_stay_in_grid_cellwise.sh          # all eight cases, each under 5 minutes
FULL=1 ./run_stay_in_grid_cellwise.sh   # the paper's refinement counts, no time limit
python compare_with_paper.py > comparison.txt
python make_comparison_images.py
```

Same setup as `../results_stay_in_grid` (16, 10 and 9 refinements for `vdp18`, `di14` and
`evasion12`). Run folders end in `_restrictgrid_cellbounds`.

**Results** (`comparison.txt`):

* The seven reach–avoid cases are identical in every phase to `../results_stay_in_grid`.
* The avoid case (Fig. 5, phase 6) changes from 63.9% safe / 30.5% unsafe in the paper's run
  (66.9 / 29.0% in `../results_stay_in_grid`) to 81.1 / 17.1%. All red edge bands turn safe,
  because leaving the grid then counts as staying safe forever
  (`comparison_images/fig5_avoid80_phase6.png`).
