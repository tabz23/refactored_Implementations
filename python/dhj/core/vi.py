"""Interval value iteration (Algorithms 1 and 3 of the paper) over a CellTree.

Port of SafetyValueIterator from the paper scripts. Two implementation changes
that do not alter any number:

  * the Bellman sweep is vectorised with numpy (exact min/max reductions over a
    CSR successor graph) instead of one multiprocessing task per cell;
  * the dynamics object is installed once per worker instead of pickled into
    every ODE / initialisation task.

Everything else - successor sets, out-of-bounds handling, convergence tests,
the order of floating-point operations - is the same. Discounted runs start
V_lower at min(l, r) and stop on the residual test; they do not subtract the
Algorithm 3 lower-bound correction.
"""
import os
import time
from multiprocessing import Pool, cpu_count
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .cell import Cell, CellTree
from .modes import Mode
from .reach import GronwallReachabilityAnalyzer
from .workers import cell_init_worker, init_worker, trajectory_worker

OOB = -1  # successor id standing for "reach box leaves the (non-periodic) domain"


class SafetyValueIterator:
    def __init__(self, mode: Mode, dyn, gamma: float, cell_tree: CellTree,
                 reachability: GronwallReachabilityAnalyzer, output_dir: str,
                 n_workers: Optional[int] = None, plotter=None, dump_csv: bool = False):
        self.mode = mode
        self.dyn = dyn
        self.gamma = gamma
        self.cell_tree = cell_tree
        self.reachability = reachability
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        self.plotter = plotter          # callable(cells, filename, iteration) or None
        self.dump_csv = dump_csv
        self.L_f, self.L_l, self.L_r = dyn.get_lipschitz_constants()
        self.actions = list(dyn.get_action_space())
        self.n_actions = len(self.actions)
        self.refinement_phase = 0
        self.n_workers = n_workers or max(1, cpu_count() - 1)
        # (cell_id, action_idx) -> np.ndarray of successor cell ids (OOB = -1)
        self.successor_cache: Dict[Tuple[int, int], np.ndarray] = {}
        self.pool = None
        self.timers = {"cell_init": 0.0, "ode": 0.0, "queries": 0.0, "vi": 0.0, "plot": 0.0}
        self.counters = {"ode_solves": 0, "vi_sweeps": 0}
        print(f"Initialized with γ={gamma}, L_f={self.L_f}, L_l={self.L_l},  L_r={self.L_r}")
        print(f"Contraction factor: γL_f = {gamma * self.L_f:.4f}")
        print(f"Using {self.n_workers} parallel workers")

    # ---- pool ----------------------------------------------------------------
    def _get_pool(self):
        if self.pool is None:
            self.pool = Pool(self.n_workers, initializer=init_worker, initargs=(self.dyn,))
        return self.pool

    def close(self):
        if self.pool is not None:
            self.pool.close()
            self.pool.join()
            self.pool = None

    def __del__(self):
        try:
            self.close()
        except Exception:  # noqa: BLE001
            pass

    # ---- cell initialisation --------------------------------------------------
    def _initialize(self, cells: Sequence[Cell], label: str):
        if len(cells) == 0:
            return
        print(f"  Initializing {len(cells)} {label}cells in parallel...")
        t0 = time.time()
        tasks = [(c.cell_id, c.bounds, c.center, self.L_l, self.L_r) for c in cells]
        chunksize = max(1, len(tasks) // (self.n_workers * 4))
        results = self._get_pool().map(cell_init_worker, tasks, chunksize=chunksize)
        by_id = {c.cell_id: c for c in cells}
        for cell_id, l_lower, l_upper, r_lower, r_upper in results:
            c = by_id[cell_id]
            c.l_lower, c.l_upper, c.r_lower, c.r_upper = l_lower, l_upper, r_lower, r_upper
            # min(l, r) is the least value of the reach-avoid Bellman image, so
            # lower iterates that start there increase and stay below the fixed
            # point. Discounted upper values still start at l and decrease.
            if self.mode.discounted:
                c.V_lower = min(c.l_lower, c.r_lower)
                c.V_upper = c.l_upper
            else:
                c.V_lower = min(c.l_lower, c.r_lower)
                c.V_upper = min(c.l_upper, c.r_upper)
        el = time.time() - t0
        self.timers["cell_init"] += el
        print(f"   Initialized in {el:.2f}s ({len(cells) / el if el > 0 else float('inf'):.1f} cells/s)")

    def initialize_cells(self):
        self._initialize(self.cell_tree.get_leaves(), "")

    def initialize_new_cells(self, new_cells: Sequence[Cell]):
        self._initialize(new_cells, "new ")

    # ---- successor sets -------------------------------------------------------
    def _compute_successors(self, cells: Sequence[Cell], banner: bool):
        """Parallel ODE integration + sequential spatial-index queries (HYBRID)."""
        n_tasks = len(cells) * self.n_actions
        if banner:
            print(f"  Step 1: Computing {n_tasks} trajectories in parallel...")
        t0 = time.time()
        tasks = [(c.cell_id, c.center.copy(), ai, a, self.reachability.tau, self.reachability.dt)
                 for c in cells for ai, a in enumerate(self.actions)]
        chunksize = max(1, len(tasks) // (self.n_workers * 4))
        trajectories = self._get_pool().map(trajectory_worker, tasks, chunksize=chunksize)
        ode_time = time.time() - t0
        self.timers["ode"] += ode_time
        self.counters["ode_solves"] += n_tasks
        if banner:
            print(f"   Trajectories computed in {ode_time:.2f}s ({n_tasks / ode_time if ode_time > 0 else float('inf'):.1f} tasks/s)")
            print("  Step 2: Computing successors using spatial index...")

        q0 = time.time()
        by_id = {c.cell_id: c for c in cells}
        leaves = self.cell_tree.leaves
        leaf_ids = np.fromiter((c.cell_id for c in leaves), dtype=np.int64, count=len(leaves))
        grid = self.dyn.get_state_bounds()
        dim = self.dyn.dim
        periodic = set(self.dyn.periodic_dims)
        bounded_dims = [d for d in range(dim) if d not in periodic]
        growth = self.reachability.growth_factors
        tree = self.cell_tree

        for cell_id, action_idx, checkpoint_states in trajectories:
            cell = by_id[cell_id]
            if self.reachability.use_infinity_norm:
                r = 0.5 * cell.get_max_range()
            else:
                r = 0.5 * np.linalg.norm(np.array([cell.get_range(j) for j in range(dim)]))
            hits = set()
            oob = False
            reach_bounds = np.zeros((dim, 2))
            for i, s in enumerate(checkpoint_states):
                expansion = r * growth[i]
                for d in range(dim):
                    reach_bounds[d, 0] = s[d] - expansion
                    reach_bounds[d, 1] = s[d] + expansion
                hits.update(tree.get_intersecting_indices(reach_bounds))
                # a periodic dimension wraps rather than leaving the domain
                for d in bounded_dims:
                    if reach_bounds[d, 0] < grid[d, 0] or reach_bounds[d, 1] > grid[d, 1]:
                        oob = True
                        break
            succ = leaf_ids[np.fromiter(hits, dtype=np.int64, count=len(hits))] if hits \
                else np.empty(0, dtype=np.int64)
            if oob:
                succ = np.append(succ, OOB)
            self.successor_cache[(cell_id, action_idx)] = succ

        q_time = time.time() - q0
        self.timers["queries"] += q_time
        if banner:
            elapsed = time.time() - t0
            print(f"   Spatial queries completed in {q_time:.2f}s")
            print(f"   Total precomputation: {elapsed:.2f}s ({n_tasks / elapsed if elapsed > 0 else float('inf'):.1f} tasks/s)")
            print(f"    Breakdown: {ode_time / elapsed * 100 if elapsed > 0 else 0:.1f}% ODE, "
                  f"{q_time / elapsed * 100 if elapsed > 0 else 0:.1f}% queries")

    def precompute_all_successors(self):
        print("\nPrecomputing successor sets with HYBRID approach...")
        self._compute_successors(self.cell_tree.get_leaves(), banner=True)

    def update_successor_cache_for_new_cells(self, new_cells: Sequence[Cell]):
        """Drop every cached entry that names a now-split parent (as owner or as
        successor), then recompute every leaf with a missing entry."""
        if not new_cells:
            return
        print(f"    Updating successor cache for {len(new_cells)} new cells (HYBRID)...")
        t0 = time.time()
        refined_parent_ids = {c.parent.cell_id for c in new_cells if c.parent is not None}
        if refined_parent_ids:
            rp = np.fromiter(refined_parent_ids, dtype=np.int64, count=len(refined_parent_ids))
            keys_to_delete = [k for k, succ in self.successor_cache.items()
                              if k[0] in refined_parent_ids or (succ.size and np.isin(succ, rp).any())]
            for k in keys_to_delete:
                del self.successor_cache[k]
        affected = [c for c in self.cell_tree.get_leaves()
                    if any((c.cell_id, ai) not in self.successor_cache for ai in range(self.n_actions))]
        print(f"      Affected cells: {len(affected)}")
        self._compute_successors(affected, banner=False)
        el = time.time() - t0
        n_tasks = len(affected) * self.n_actions
        print(f"     Cache updated in {el:.2f}s ({n_tasks / el if el > 0 else float('inf'):.1f} tasks/s)")

    # ---- Bellman sweeps -------------------------------------------------------
    def _build_graph(self, leaves: List[Cell]):
        """CSR successor graph over leaf positions; position n is the OOB sink."""
        n = len(leaves)
        A = self.n_actions
        lookup = np.full(self.cell_tree.next_id + 1, -1, dtype=np.int64)  # index = cell_id + 1
        lookup[0] = n                                                    # OOB -> sink
        for i, c in enumerate(leaves):
            lookup[c.cell_id + 1] = i
        chunks, owners = [], []
        for i, c in enumerate(leaves):
            base = i * A
            for ai in range(A):
                succ = self.successor_cache.get((c.cell_id, ai))
                if succ is not None and succ.size:
                    chunks.append(succ)
                    owners.append(np.full(succ.size, base + ai, dtype=np.int64))
        if chunks:
            ids = np.concatenate(chunks)
            own = np.concatenate(owners)
            pos = lookup[ids + 1]
            keep = pos >= 0            # successors that are no longer leaves are skipped
            pos, own = pos[keep], own[keep]
        else:
            pos = np.empty(0, dtype=np.int64)
            own = np.empty(0, dtype=np.int64)
        lengths = np.bincount(own, minlength=n * A)
        starts = np.concatenate(([0], np.cumsum(lengths)[:-1]))
        return pos, lengths, starts

    def _run_sweeps(self, max_iterations: int, convergence_tol: float, conservative_mode: bool,
                    delta_max: float, plot_freq: int, local: bool):
        t_start = time.time()
        leaves = self.cell_tree.leaves
        n = len(leaves)
        A = self.n_actions
        idx, lengths, starts = self._build_graph(leaves)
        nz = lengths > 0
        starts_nz = starts[nz]

        Vu = np.fromiter((c.V_upper for c in leaves), dtype=float, count=n)
        Vl = np.fromiter((c.V_lower for c in leaves), dtype=float, count=n)
        lu = np.fromiter((c.l_upper for c in leaves), dtype=float, count=n)
        ll = np.fromiter((c.l_lower for c in leaves), dtype=float, count=n)
        ru = np.fromiter((c.r_upper for c in leaves), dtype=float, count=n)
        rl = np.fromiter((c.r_lower for c in leaves), dtype=float, count=n)
        Vu_ext = np.append(Vu, -1.0)   # index n: out-of-bounds sink, fixed at -1
        Vl_ext = np.append(Vl, -1.0)

        conv_upper, conv_lower = [], []
        converged = False
        iteration = 0
        diff_upper = diff_lower = delta_k = float("nan")
        min_lower_step = float("inf")
        conservative_mode = self._use_lower_correction(conservative_mode)
        pad = "    Local Iteration" if local else "Iteration"
        gamma = self.gamma

        for iteration in range(max_iterations):
            mins = np.full(n * A, -np.inf)
            maxs = np.full(n * A, -np.inf)
            if idx.size:
                mins[nz] = np.minimum.reduceat(Vl_ext[idx], starts_nz)
                maxs[nz] = np.maximum.reduceat(Vu_ext[idx], starts_nz)
            if self.mode.discounted:
                mins[nz] *= gamma
                maxs[nz] *= gamma
            best_min = mins.reshape(n, A).max(axis=1)
            best_max = maxs.reshape(n, A).max(axis=1)
            if self.mode is Mode.AVOID_NODISCOUNT:
                new_Vl = np.minimum(ll, best_min)
                new_Vu = np.minimum(lu, best_max)
            else:
                new_Vl = np.minimum(ll, np.maximum(rl, best_min))
                new_Vu = np.minimum(lu, np.maximum(ru, best_max))

            diff_upper = float(np.max(np.abs(new_Vu - Vu)))
            diff_lower = float(np.max(np.abs(new_Vl - Vl)))
            delta_k = float(np.min(new_Vl - Vl))
            if self.mode.discounted:
                min_lower_step = min(min_lower_step, delta_k)
            Vu, Vl = new_Vu, new_Vl
            Vu_ext[:n] = Vu
            Vl_ext[:n] = Vl
            conv_upper.append(diff_upper)
            conv_lower.append(diff_lower)
            self.counters["vi_sweeps"] += 1

            if conservative_mode:
                print(f"{pad} {iteration + 1:3d}: "
                      f"||V̄^k - V̄^k-1||_∞ = {diff_upper:12.20f}, "
                      f"||V_^k - V_^k-1||_∞ = {diff_lower:12.20f}, "
                      f"δ^k = {delta_k:12.20f}")
                if delta_k >= -delta_max and diff_upper < convergence_tol:
                    Vl = self._conservative_correction(Vl, delta_k, diff_upper, convergence_tol)
                    converged = True
            else:
                print(f"{pad} {iteration + 1:3d}: "
                      f"||V̄^k - V̄^k-1||_∞ = {diff_upper:12.20f}, "
                      f"||V_^k - V_^k-1||_∞ = {diff_lower:12.20f}")
                if diff_upper < convergence_tol and diff_lower < convergence_tol:
                    print("\n" + "=" * 50)
                    print(f" {'LOCAL' if local else 'STANDARD'} CONVERGENCE ACHIEVED")
                    print("=" * 50)
                    print(f"  ||V̄^k - V̄^k-1||_∞ = {diff_upper:.20e} < tolerance = {convergence_tol}")
                    print(f"  ||V_^k - V_^k-1||_∞ = {diff_lower:.20e} < tolerance = {convergence_tol}")
                    print(f"  Converged in {iteration + 1} iterations")
                    converged = True

            if converged:
                iteration += 1
                break

            if plot_freq > 0 and (iteration + 1) % plot_freq == 0:
                self._write_back(leaves, Vu, Vl)
                fn = os.path.join(self.output_dir,
                                  f"iteration_{iteration + 1:04d}_refinement_{self.refinement_phase:02d}.png")
                self.emit_output(fn, iteration + 1, self.refinement_phase, csv_name=None)
                print(f"  [Plot saved at iteration {iteration + 1}]")

        if not converged:
            iteration = max_iterations
            print("\n" + "!" * 50)
            print(f"  MAXIMUM ITERATIONS REACHED: {max_iterations}")
            print("!" * 50)
            if conservative_mode:
                print(f"  Final δ^k = {delta_k:.20e}")
            print(f"  Final ||V̄^k - V̄^k-1||_∞ = {diff_upper:.20e}")
            print(f"  Final ||V_^k - V_^k-1||_∞ = {diff_lower:.20e}")
            if conservative_mode:
                delta_to_use = abs(delta_k)
                epsilon_cons = (gamma * delta_to_use) / (1 - gamma)
                print(f"  Using δ_actual = {delta_to_use:.10e} for correction")
                print(f"  Applying conservative margin: ε_cons = {epsilon_cons:.10e}")
                Vl = Vl - epsilon_cons
                print(f"   Conservative correction applied to {n} cells")

        if self.mode.discounted and iteration > 0:
            print(f"  V_lower increased from min(l_lower, r_lower); smallest step {min_lower_step:.3e}")
        self._write_back(leaves, Vu, Vl)
        self.timers["vi"] += time.time() - t_start
        return np.array(conv_upper), np.array(conv_lower), iteration, converged

    def _use_lower_correction(self, conservative_mode: bool) -> bool:
        """Discounted lower iterates increase from min(l, r), already a lower
        bound of every Bellman image, so the Algorithm 3 margin is not applied."""
        if self.mode.discounted and conservative_mode:
            print("Discounted V_lower starts at min(l_lower, r_lower); residual stopping, no lower-bound correction.")
            return False
        return conservative_mode

    def _conservative_correction(self, Vl, delta_k, diff_upper, convergence_tol):
        print("\n" + "=" * 60)
        print(" CONSERVATIVE STOPPING CONDITION MET (Algorithm 3)")
        print("=" * 60)
        print(f"  δ^k = {delta_k:.10e} ≥ -δ_max")
        print(f"  ||V̄^k - V̄^k-1||_∞ = {diff_upper:.10e} < {convergence_tol}")
        delta_to_use = abs(delta_k)
        epsilon_cons = (self.gamma * delta_to_use) / (1 - self.gamma)
        print(f"  ε_cons = ({self.gamma} * {delta_to_use:.10e}) / {1 - self.gamma} = {epsilon_cons:.10e}")
        print(f"  Correcting V_lower for {len(Vl)} cells: V_lower ← V_lower - ε_cons - ε_machine")
        eps_machine = np.finfo(float).eps
        old = Vl
        new = old - epsilon_cons - eps_machine
        near_zero = int(np.sum(np.abs(old - epsilon_cons) < eps_machine))
        flipped = int(np.sum((old - epsilon_cons > 0) & (new <= 0)))
        print("   Conservative correction applied successfully")
        print(f"  Machine epsilon ε_machine = {eps_machine:.3e}")
        if near_zero > 0:
            print(f"    {near_zero} cells experienced near-cancellation")
        if flipped > 0:
            print(f"    {flipped} cells changed from (+) to (≤0) due to ε_machine")
        return new

    @staticmethod
    def _write_back(leaves: List[Cell], Vu: np.ndarray, Vl: np.ndarray):
        for i, c in enumerate(leaves):
            c.V_upper = float(Vu[i])
            c.V_lower = float(Vl[i])

    # ---- public drivers -------------------------------------------------------
    def value_iteration(self, max_iterations: int, convergence_tol: float, plot_freq: int,
                        conservative_mode: bool, delta_max: float):
        """Phase-0 / Algorithm-1 value iteration on the current leaves."""
        self.initialize_cells()
        if not self.successor_cache:
            print("  Precomputing successor sets...")
            self.precompute_all_successors()
        print(f"\nStarting OPTIMIZED PARALLEL value iteration (max {max_iterations} iterations)...")
        conservative_mode = self._use_lower_correction(conservative_mode)
        print(f"Conservative mode: {conservative_mode}")
        if conservative_mode:
            print(f"Conservative tolerance δ_max: {delta_max}")
            print("Conservative margin will be: ε_cons = (γ * δ_actual) / (1 - γ)")
        else:
            print(f"Convergence tolerance: {convergence_tol}")
        print(f"Number of cells: {self.cell_tree.get_num_leaves()}")

        cu, cl, iters, converged = self._run_sweeps(max_iterations, convergence_tol, conservative_mode,
                                                    delta_max, plot_freq, local=False)
        fn = os.path.join(self.output_dir, "value_function_phase_0_complete.png")
        self.emit_output(fn, iters, 0)
        print(f"\nValue iteration completed in {iters} iterations")
        return cu, cl, iters, converged

    def local_value_iteration(self, new_cells: Sequence[Cell], max_iterations: int, convergence_tol: float,
                              conservative_mode: bool, delta_max: float):
        """Re-initialise every leaf, refresh the successor cache, sweep to convergence."""
        leaves = self.cell_tree.get_leaves()
        kind = "V_lower = min(l, r), V_upper = l" if self.mode.discounted else "V = min(l, r)"
        print(f"  Reinitializing ALL {len(leaves)} cells ({kind})")
        t0 = time.time()
        self.initialize_cells()
        reinit = time.time() - t0
        print(f"   Reinitialized {len(leaves)} cells in {reinit:.2f}s ({len(leaves) / reinit if reinit > 0 else float('inf'):.1f} cells/s)")
        conservative_mode = self._use_lower_correction(conservative_mode)
        print(f"  Local VI: updating all {len(leaves)} cells")
        print(f"    Conservative mode: {conservative_mode}")
        if conservative_mode:
            print(f"    δ_max: {delta_max}")
        else:
            print(f"    Convergence tolerance: {convergence_tol}")
        print(f"    Max iterations: {max_iterations}")
        if len(leaves) == 0:
            print("    No cells to update!")
            return np.array([]), np.array([]), 0, True
        if new_cells:
            print(f"    Updating successor cache for {len(new_cells)} updated cells...")
            self.update_successor_cache_for_new_cells(new_cells)
        elif not self.successor_cache:
            # resumed from a checkpoint: the cache is rebuilt in full (exact)
            print("    Successor cache is empty (resume) - computing all successor sets...")
            self.precompute_all_successors()

        cu, cl, iters, converged = self._run_sweeps(max_iterations, convergence_tol, conservative_mode,
                                                    delta_max, 0, local=True)
        print(f"    Local VI completed in {iters} iterations")
        if len(cu) > 0:
            print("    Final convergence values:")
            print(f"      ||V̄^final - V̄^prev||_∞ = {cu[-1]:.8e}")
            print(f"      ||V_^final - V_^prev||_∞ = {cl[-1]:.8e}")
        return cu, cl, iters, converged

    # ---- analysis / output -------------------------------------------------------
    def identify_boundary_cells(self) -> List[Cell]:
        """Cells with V_upper > 0 and V_lower <= 0 (undetermined)."""
        return [c for c in self.cell_tree.get_leaves()
                if c.V_upper is not None and c.V_lower is not None and c.V_upper > 0 and c.V_lower <= 0]

    def classification_counts(self):
        safe = unsafe = boundary = 0
        for c in self.cell_tree.get_leaves():
            if c.V_lower is not None and c.V_upper is not None:
                if c.V_lower > 0:
                    safe += 1
                elif c.V_upper <= 0:
                    unsafe += 1
                else:
                    boundary += 1
        return safe, unsafe, boundary

    def print_statistics(self, heading="Final Cell Classification:"):
        """Safe / unsafe / unclassified, same rule as the classification panel.

        Reports both the share of leaves and the share of state-space volume
        (product of each leaf's side lengths). One pass over the leaves.
        stdout is copied into run.log.
        """
        safe = unsafe = unclassified = 0
        safe_v = unsafe_v = unclassified_v = 0.0
        dim = self.dyn.dim
        for c in self.cell_tree.get_leaves():
            vol = 1.0
            for d in range(dim):
                vol *= c.get_range(d)
            if c.V_lower is not None and c.V_upper is not None and c.V_lower > 0:
                safe += 1
                safe_v += vol
            elif c.V_lower is not None and c.V_upper is not None and c.V_upper <= 0:
                unsafe += 1
                unsafe_v += vol
            else:
                unclassified += 1
                unclassified_v += vol
        n = max(1, safe + unsafe + unclassified)
        v = safe_v + unsafe_v + unclassified_v
        if v <= 0.0:
            v = 1.0
        print(f"\n{heading}")
        print(f"  Safe:         {safe:8d} cells ({100 * safe / n:5.1f}%),  {100 * safe_v / v:5.1f}% of state space")
        print(f"  Unsafe:       {unsafe:8d} cells ({100 * unsafe / n:5.1f}%),  {100 * unsafe_v / v:5.1f}% of state space")
        print(f"  Unclassified: {unclassified:8d} cells ({100 * unclassified / n:5.1f}%),  {100 * unclassified_v / v:5.1f}% of state space")

    def emit_output(self, png_filename: str, iteration: int, phase: int, csv_name: Optional[str] = ""):
        """Figure (if a plotter is set) and optional CSV dump of the leaves."""
        t0 = time.time()
        if self.plotter is not None:
            self.plotter(self.cell_tree.get_leaves(), png_filename, iteration)
        if self.dump_csv and csv_name is not None:
            name = csv_name or f"value_function_phase_{phase}.csv"
            write_cells_csv(self.dyn, self.cell_tree.get_leaves(), os.path.join(self.output_dir, name))
        self.timers["plot"] += time.time() - t0


def write_cells_csv(dyn, cells: Sequence[Cell], path: str):
    """Generic per-leaf dump: cell_id, d<k>_lo, d<k>_hi (k < dim), values."""
    dim = dyn.dim
    header = ["cell_id"] + [f"d{k}_{s}" for k in range(dim) for s in ("lo", "hi")] + \
             ["V_lower", "V_upper", "l_lower", "l_upper", "r_lower", "r_upper"]
    with open(path, "w", encoding="utf-8") as f:
        f.write(",".join(header) + "\n")
        for c in cells:
            b = c.bounds
            row = [str(c.cell_id)] + [repr(float(b[k, s])) for k in range(dim) for s in (0, 1)] + \
                  [repr(float(v)) for v in (c.V_lower, c.V_upper, c.l_lower, c.l_upper, c.r_lower, c.r_upper)]
            f.write(",".join(row) + "\n")
