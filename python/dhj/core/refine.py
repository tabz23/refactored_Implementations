"""Algorithm 2/3: adaptive refinement of the undetermined (boundary) cells.

Each phase: identify boundary cells larger than eta_min = eps / (2 L_l),
bisect them along their widest dimension, rebuild the spatial index, run the
local value iteration on the whole leaf set, save the figure, write a
checkpoint and a row in phases.csv.
"""
import os
import time
from typing import Optional

import numpy as np

from .cell import CellTree
from .checkpoint import prune_checkpoints, save_checkpoint
from .params import PhaseLog
from .vi import SafetyValueIterator


class AdaptiveRefinement:
    def __init__(self, vi: SafetyValueIterator, cell_tree: CellTree, out_dir: str, args):
        self.vi = vi
        self.cell_tree = cell_tree
        self.out_dir = out_dir
        self.args = args
        self.L_l = vi.dyn.L_l
        self.phase_log = PhaseLog(out_dir)
        self.t_start = time.time()
        self._last_saved_phase = -1

    # ---- helpers ------------------------------------------------------------
    def _boundary_and_refinable(self, eta_min: float):
        boundary = self.vi.identify_boundary_cells()
        refinable = [c for c in boundary if c.get_max_range() > eta_min]
        return boundary, refinable

    def _checkpoint(self, phase: int, total_refined: int, epsilon: float, force: bool = False):
        every = self.args.checkpoint_every
        if every <= 0:
            return
        if not force and phase % every != 0:
            return
        if phase == self._last_saved_phase:
            return
        t0 = time.time()
        meta = dict(mode=str(self.vi.mode.value), dynamics=self.vi.dyn.name, epsilon=epsilon,
                    initial_resolution=self.args.initial_resolution)
        path = save_checkpoint(self.out_dir, self.cell_tree, phase, total_refined, meta)
        prune_checkpoints(self.out_dir, self.args.keep_checkpoints)
        self._last_saved_phase = phase
        mb = os.path.getsize(path) / (1024 * 1024)
        print(f"  Checkpoint written: {path} ({mb:.1f} MB, {time.time() - t0:.2f}s)")

    # ---- main loop ------------------------------------------------------------
    def refine(self, epsilon: float, max_refinements: int, vi_iterations_per_refinement: int,
               resume: Optional[dict] = None):
        eta_min = epsilon / (2 * self.L_l)
        a = self.args
        print("\n" + "=" * 70)
        print("ADAPTIVE REFINEMENT CONFIGURATION")
        print("=" * 70)
        print(f"  Error tolerance ε: {epsilon}")
        print(f"  Minimum cell size η_min: {eta_min:.6f}")
        print(f"  Maximum refinements: {max_refinements}")
        print(f"  VI iterations per refinement: {vi_iterations_per_refinement}")
        print(f"  Conservative mode: {a.conservative}")
        if a.conservative:
            print(f"  δ_max: {a.delta_max}")

        refinement_iter = 0
        total_refined = 0
        if resume is not None:
            refinement_iter = int(resume["phase"])
            total_refined = int(resume["total_refined"])
            self.vi.refinement_phase = refinement_iter
            print("\n" + "=" * 70)
            print("RESUMING FROM CHECKPOINT")
            print("=" * 70)
            print(f"  Phases already completed: {refinement_iter}")
            print(f"  Cumulative refined:       {total_refined}")
            print(f"  Leaf cells:               {self.cell_tree.get_num_leaves()}")
            print("  Successor cache will be rebuilt in full on the first local sweep (exact).")
        else:
            print("\n" + "=" * 70)
            print("PHASE 0: INITIAL VALUE ITERATION")
            print("=" * 70)
            print(f"Grid: {self.cell_tree.get_num_leaves()} cells")
            self.vi.refinement_phase = 0
            t0 = time.time()
            succ0 = self.vi.timers["ode"] + self.vi.timers["queries"]
            vi0 = self.vi.timers["vi"]
            conv_upper, conv_lower, iters, converged = self.vi.value_iteration(
                max_iterations=vi_iterations_per_refinement, convergence_tol=a.phase0_tol,
                plot_freq=a.plot_freq, conservative_mode=a.conservative, delta_max=a.delta_max)
            phase_time = time.time() - t0
            print(f"total time for running value_iteration the first time (phase 0) on the initial grid "
                  f"including grid initialization, setting up cpu parallelism, successor set computation "
                  f"and value iteration is {phase_time} seconds")
            if len(conv_upper) > 0:
                print(f"Initial VI completed in {len(conv_upper)} iterations")
                print(f"Final convergence: ||V̄||_∞ = {conv_upper[-1]:.8e}, ||V_||_∞ = {conv_lower[-1]:.8e}")
            print(f"Saved initial visualization: {os.path.join(self.out_dir, 'value_function_phase_0_complete.png')}")
            self.vi.print_statistics("Classification after phase 0:")
            boundary, refinable = self._boundary_and_refinable(eta_min)
            self.phase_log.append(phase=0, n_leaves=self.cell_tree.get_num_leaves(),
                                  boundary_cells=len(boundary), refinable_cells=len(refinable),
                                  refined_parents=0, vi_iterations=iters, converged=converged,
                                  time_vi_s=round(self.vi.timers["vi"] - vi0, 3),
                                  time_successors_s=round(self.vi.timers["ode"] + self.vi.timers["queries"] - succ0, 3),
                                  time_phase_s=round(phase_time, 3),
                                  wall_since_start_s=round(time.time() - self.t_start, 3))
            self._checkpoint(0, 0, epsilon)

        boundary, refinable = self._boundary_and_refinable(eta_min)
        print("\n" + "=" * 70)
        print("INITIAL QUEUE STATE")
        print("=" * 70)
        print(f"Total boundary cells: {len(boundary)}")
        print(f"  Refinable (>η_min={eta_min:.6f}): {len(refinable)}")
        print(f"  Below threshold: {len(boundary) - len(refinable)}")

        while refinement_iter < max_refinements and len(refinable) > 0:
            boundary, refinable = self._boundary_and_refinable(eta_min)
            print("\n" + "=" * 70)
            print(f"REFINEMENT PHASE {refinement_iter + 1}")
            print("=" * 70)
            phase_start = time.time()
            succ0 = self.vi.timers["ode"] + self.vi.timers["queries"]
            vi0 = self.vi.timers["vi"]
            print(f"Boundary cells: {len(boundary)}")
            print(f"  Refinable (>η_min={eta_min:.6f}): {len(refinable)}")
            print(f"  Below threshold: {len(boundary) - len(refinable)}")
            print(f"  Cumulative refined: {total_refined}")
            if len(refinable) == 0:
                print("  No refinable cells remaining - stopping refinement")
                break

            sizes = [c.get_max_range() for c in refinable]
            print("  Refinable cell size statistics:")
            print(f"    Max:    {max(sizes):.6f}")
            print(f"    Mean:   {np.mean(sizes):.6f}")
            print(f"    Min:    {min(sizes):.6f}")
            print(f"    Median: {np.median(sizes):.6f}")

            print(f"  Refining {len(refinable)} cells...")
            t0 = time.time()
            new_cells = self.cell_tree.refine_cells(refinable)
            refinement_time = time.time() - t0
            n_refined = len(refinable)
            total_refined += n_refined
            print(f"    Refined: {n_refined} parent cells")
            print(f"    Created: {len(new_cells)} child cells")
            print(f"    Total cells: {self.cell_tree.get_num_leaves()}")
            print(f"    Refinement time: {refinement_time:.3f}s")
            print(f"    Cumulative refined: {total_refined}")

            print("  Rebuilding spatial index...")
            t0 = time.time()
            self.cell_tree.rebuild_spatial_index()
            print(f"    Spatial index rebuilt in {time.time() - t0:.3f}s")

            if new_cells:
                print(f"  Initializing {len(new_cells)} new cells...")
                t0 = time.time()
                self.vi.initialize_new_cells(new_cells)
                print(f"    Initialization completed in {time.time() - t0:.3f}s")

            self.vi.refinement_phase = refinement_iter + 1
            print("  Starting local value iteration...")
            t0 = time.time()
            conv_upper, conv_lower, iters, converged = self.vi.local_value_iteration(
                new_cells, max_iterations=vi_iterations_per_refinement, convergence_tol=a.tolerance,
                conservative_mode=a.conservative, delta_max=a.delta_max)
            local_vi_time = time.time() - t0
            print(f"  Local VI completed in {len(conv_upper)} iterations, time: {local_vi_time:.3f}s")
            if len(conv_upper) > 0:
                print(f"    Final ||V̄^k - V̄^k-1||_∞ = {conv_upper[-1]:.20e}")
                print(f"    Final ||V_^k - V_^k-1||_∞ = {conv_lower[-1]:.20e}")

            fn = os.path.join(self.out_dir, f"value_function_phase_{refinement_iter + 1}_complete.png")
            self.vi.emit_output(fn, refinement_iter + 1, refinement_iter + 1)
            print(f"  Saved visualization: {fn}")

            phase_time = time.time() - phase_start
            n_leaves = self.cell_tree.get_num_leaves()
            print(f"  Phase {refinement_iter + 1} completed in {phase_time:.2f}s with {n_leaves} leaf cells")
            self.vi.print_statistics(f"Classification after phase {refinement_iter + 1}:")
            refinement_iter += 1

            boundary_after, refinable_after = self._boundary_and_refinable(eta_min)
            self.phase_log.append(phase=refinement_iter, n_leaves=n_leaves,
                                  boundary_cells=len(boundary_after), refinable_cells=len(refinable_after),
                                  refined_parents=n_refined, vi_iterations=iters, converged=converged,
                                  time_vi_s=round(self.vi.timers["vi"] - vi0, 3),
                                  time_successors_s=round(self.vi.timers["ode"] + self.vi.timers["queries"] - succ0, 3),
                                  time_phase_s=round(phase_time, 3),
                                  wall_since_start_s=round(time.time() - self.t_start, 3))
            self._checkpoint(refinement_iter, total_refined, epsilon)
            # NOTE: like the paper scripts, `refinable` is re-evaluated at the top
            # of the next phase (which prints "No refinable cells remaining").

        self._checkpoint(refinement_iter, total_refined, epsilon, force=True)

        print("\n" + "=" * 70)
        print("ADAPTIVE REFINEMENT COMPLETE")
        print("=" * 70)
        print("Refinement Summary:")
        print(f"  Phases completed: {refinement_iter}")
        print(f"  Parent cells refined: {total_refined}")
        print(f"  Child cells created: {total_refined * 2}")
        print(f"  Final cells: {self.cell_tree.get_num_leaves()}")
        final_boundary, final_refinable = self._boundary_and_refinable(eta_min)
        print("Final Boundary State:")
        print(f"  Boundary cells: {len(final_boundary)}")
        print(f"  Still refinable: {len(final_refinable)}")
        print(f"  Below threshold: {len(final_boundary) - len(final_refinable)}")
        self.vi.print_statistics()
        print(f"\n All results saved to: {self.out_dir}/")
