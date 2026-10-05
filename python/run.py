#!/usr/bin/env python
"""Discrete HJ reach-avoid / avoid solver - command line entry point.

One script replaces RA_nodiscount.py, RA_discount.py and avoid_nodiscount.py:
pick the Bellman backup with --mode and the system with --dynamics. All other
flags keep the names and defaults of the original scripts.

Examples (parameters of the paper runs, see README.md):

  python run.py --mode ra_nodiscount --dynamics dubins  --algorithm 2 --dt 0.3 --tau 0.3 \
      --gamma 0.9 --delta-min 1e-8 --delta-max 1e-8 --epsilon 0.05 --initial-resolution 20 --conservative

  python run.py --mode ra_nodiscount --dynamics evasion --algorithm 2 --dt 0.3 --tau 0.3 \
      --gamma 0.96 --delta-min 1e-8 --delta-max 1e-8 --epsilon 0.05 --initial-resolution 30 --conservative

  python run.py --mode avoid_nodiscount --dynamics dubins --algorithm 2 --dt 0.3 --tau 0.3 \
      --gamma 0.9 --delta-min 1e-8 --delta-max 1e-8 --epsilon 0.05 --initial-resolution 80 --conservative

  python run.py --resume auto --out-dir results/ra_nodiscount/<run>      # continue a run
"""
import argparse
import os
import sys
import time
from functools import partial
from multiprocessing import cpu_count

from dhj.core.cell import CellTree
from dhj.core.checkpoint import latest_checkpoint, load_checkpoint
from dhj.core.logging_utils import RunLogger
from dhj.core.modes import Mode
from dhj.core.params import write_params
from dhj.core.reach import GronwallReachabilityAnalyzer
from dhj.core.refine import AdaptiveRefinement
from dhj.core.vi import SafetyValueIterator
from dhj.dynamics import REGISTRY, make_dynamics
from dhj.dynamics.plotting import PlotStyle, plot_value_function


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Discrete HJ reachability (modular rewrite of the paper scripts)",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--mode", type=str, choices=[m.value for m in Mode], default=Mode.RA_NODISCOUNT.value,
                   help="Bellman backup: ra_nodiscount (RA_nodiscount.py), ra_discount (RA_discount.py), "
                        "avoid_nodiscount (avoid_nodiscount.py)")
    p.add_argument("--algorithm", type=int, choices=[1, 2], default=2,
                   help="1: value iteration on a fixed grid; 2: adaptive refinement (Alg. 2/3)")
    p.add_argument("--dynamics", type=str, default="dubins", choices=sorted(REGISTRY),
                   help="dynamics model (see dhj/dynamics/)")
    # dynamics parameters
    p.add_argument("--velocity", type=float, default=1.0, help="constant speed (dubins, evasion)")
    p.add_argument("--a-max", type=float, default=1.0, help="acceleration bound (double_integrator)")
    p.add_argument("--v-max", type=float, default=1.0, help="velocity bound (double_integrator)")
    p.add_argument("--dt", type=float, default=0.1, help="checkpoint interval for the reach sets")
    p.add_argument("--tau", type=float, default=1.0, help="control duration (action hold time)")
    p.add_argument("--gamma", type=float, default=0.96,
                   help="discount factor; also sets the conservative margin eps_cons = gamma*delta/(1-gamma)")
    # algorithm 1
    p.add_argument("--resolution", type=int, default=10, help="grid resolution per dimension (Alg. 1)")
    p.add_argument("--iterations", type=int, default=None,
                   help="max value iterations (Alg. 1); default per mode: 200 / 2000 / 20000")
    p.add_argument("--tolerance", "--delta-min", type=float, default=1e-13, dest="tolerance",
                   help="convergence tolerance delta_min")
    p.add_argument("--plot-freq", type=int, default=1000, help="save a figure every N sweeps (0 = never)")
    # algorithm 2
    p.add_argument("--epsilon", type=float, default=0.1, help="error tolerance eps; eta_min = eps / (2 L_l)")
    p.add_argument("--initial-resolution", type=int, default=15, help="initial coarse grid per dimension")
    p.add_argument("--refinements", type=int, default=100, help="maximum number of refinement phases")
    p.add_argument("--vi-iterations", type=int, default=20000, help="max VI sweeps per phase")
    p.add_argument("--conservative", action="store_true", default=False,
                   help="conservative stopping (Algorithm 3). Ignored for ra_discount, "
                        "which starts V_lower at min(l, r) and uses residual stopping")
    p.add_argument("--no-conservative", action="store_true", default=False,
                   help="force conservative stopping off")
    p.add_argument("--delta-max", type=float, default=1e-6, help="delta_max for the conservative stop")
    p.add_argument("--phase0-tol", type=float, default=None,
                   help="residual tolerance for phase 0 (default: same as --delta-min)")
    p.add_argument("--workers", type=int, default=None, help="worker processes (default: CPU count - 1)")
    p.add_argument("--precompute", action="store_true", help="precompute all successor sets up front")
    # output / checkpointing
    p.add_argument("--results-root", type=str, default=None,
                   help="root of the results tree (default: <refactored_Implementations>/results)")
    p.add_argument("--out-dir", type=str, default=None, help="explicit output directory (overrides the derived name)")
    p.add_argument("--tag", type=str, default="", help="suffix appended to the derived output directory name")
    p.add_argument("--no-plot", action="store_true", help="skip figures")
    p.add_argument("--plot-dpi", type=int, default=800,
                   help="savefig dpi. 800 matches the paper scripts. A smaller value writes a smaller PNG")
    p.add_argument("--dump-csv", action="store_true", help="write value_function_phase_N.csv per phase")
    p.add_argument("--checkpoint-every", type=int, default=1, help="checkpoint every N phases (0 = never)")
    p.add_argument("--keep-checkpoints", type=int, default=2, help="keep only the newest K checkpoints (0 = all)")
    p.add_argument("--resume", type=str, default=None,
                   help="path to a checkpoint .npz, or 'auto' for the newest one in the output directory")
    return p


# refactored_Implementations/  (parent of this file's directory)
_IMPLEMENTATIONS_ROOT = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_RESULTS_ROOT = os.path.join(os.path.dirname(_IMPLEMENTATIONS_ROOT), "results")


def derive_out_dir(args, mode: Mode) -> str:
    if args.out_dir:
        return args.out_dir
    root = args.results_root if args.results_root else _DEFAULT_RESULTS_ROOT
    cons = "cons" if args.conservative else "std"
    if args.algorithm == 1:
        name = (f"algorithm1_{args.dynamics}_res{args.resolution}_dmin{args.tolerance:.0e}"
                f"_gamma{args.gamma}_dt{args.dt}_tau{args.tau}_{cons}")
    else:
        name = (f"{args.dynamics}_res{args.initial_resolution}_eps{args.epsilon}_dmin{args.tolerance:.0e}"
                f"_dmax{args.delta_max:.0e}_gamma{args.gamma}_dt{args.dt}_tau{args.tau}_{cons}")
    if args.tag:
        name += f"_{args.tag}"
    return os.path.join(root, mode.value, name)


def make_plotter(args, mode: Mode, dyn):
    if args.no_plot:
        return None
    # A dynamics class may provide its own figure: def plot(self, cells, filename, iteration)
    custom = getattr(dyn, "plot", None)
    if callable(custom):
        return custom
    if args.plot_dpi == 800:
        style = PlotStyle(draw_value_cells=True, cell_labels=mode.discounted, draw_target=mode.has_target)
    else:
        style = PlotStyle(draw_value_cells=True, cell_labels=mode.discounted, draw_target=mode.has_target,
                          dpi=args.plot_dpi, save_dpi=args.plot_dpi)
    return partial(plot_value_function, dyn, style=style)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.conservative and args.no_conservative:
        parser.error("pass only one of --conservative and --no-conservative")
    if args.phase0_tol is None:
        args.phase0_tol = args.tolerance
    mode = Mode(args.mode)
    if mode.discounted and args.conservative:
        print("warning: --conservative is ignored for ra_discount; "
              "V_lower starts at min(l_lower, r_lower) and residual stopping is used",
              file=sys.stderr)
        args.conservative = False
    if args.workers is None:
        args.workers = max(1, cpu_count() - 1)
    if args.iterations is None:
        args.iterations = mode.default_vi_iterations_alg1
    out_dir = derive_out_dir(args, mode)
    os.makedirs(out_dir, exist_ok=True)

    with RunLogger(out_dir):
        print("DISCRETE HJ REACHABILITY - MODULAR IMPLEMENTATION")
        print(f"Mode: {mode.value} (reproduces {mode.script_name})")
        print(f"Algorithm: {args.algorithm}, Workers: {args.workers}")
        print(f"Dynamics: {args.dynamics}")
        print(f"Conservative mode: {args.conservative}")

        dyn = make_dynamics(args.dynamics, args)
        dyn.validate()
        t_total = time.time()

        if args.algorithm == 1:
            run_algorithm_1(args, mode, dyn, out_dir)
        else:
            run_algorithm_2(args, mode, dyn, out_dir)
        print(f"TOTAL WALL TIME: {time.time() - t_total:.2f} seconds")
    return 0


def run_algorithm_1(args, mode, dyn, out_dir):
    print("=" * 70)
    print("ALGORITHM 1: Discretization Routine")
    print("=" * 70)
    print(f"Environment: {type(dyn).__name__}")
    print(f"Initializing grid with resolution {args.resolution}^{dyn.dim}...")
    tree = CellTree(dyn.get_state_bounds(), initial_resolution=args.resolution, periodic_dims=dyn.periodic_dims)
    reach = GronwallReachabilityAnalyzer(dyn)
    write_params(out_dir, args, dyn, {"reachability": reach.describe(), "algorithm": 1})
    vi = SafetyValueIterator(mode, dyn, args.gamma, tree, reach, out_dir, n_workers=args.workers,
                             plotter=make_plotter(args, mode, dyn), dump_csv=args.dump_csv)
    try:
        t0 = time.time()
        vi.value_iteration(max_iterations=args.iterations, convergence_tol=args.tolerance,
                           plot_freq=args.plot_freq, conservative_mode=args.conservative,
                           delta_max=args.delta_max)
        print("\nALGORITHM 1 COMPLETE")
        print(f"Total time: {time.time() - t0:.2f} seconds")
        print(f"Results saved to: {out_dir}/")
        vi.print_statistics()
    finally:
        vi.close()


def run_algorithm_2(args, mode, dyn, out_dir):
    print("=" * 70)
    print("ALGORITHM 2/3: Adaptive Refinement")
    print("=" * 70)
    print(f"Environment: {type(dyn).__name__}")

    resume_meta = None
    if args.resume:
        path = latest_checkpoint(out_dir) if args.resume == "auto" else args.resume
        if not path or not os.path.exists(path):
            raise SystemExit(f"--resume: no checkpoint found ({args.resume}) in {out_dir}/checkpoints")
        print(f"Loading checkpoint {path}")
        tree, resume_meta = load_checkpoint(path)
        if resume_meta.get("mode") != mode.value or resume_meta.get("dynamics") != dyn.name:
            raise SystemExit(f"checkpoint was written for mode={resume_meta.get('mode')}, "
                             f"dynamics={resume_meta.get('dynamics')}; this run is {mode.value}/{dyn.name}")
        if resume_meta.get("initial_resolution") != args.initial_resolution:
            print(f"  warning: checkpoint was made with --initial-resolution {resume_meta.get('initial_resolution')}"
                  f" but this run passes {args.initial_resolution}; the checkpoint's grid wins")
    else:
        tree = CellTree(dyn.get_state_bounds(), initial_resolution=args.initial_resolution,
                        periodic_dims=dyn.periodic_dims)

    reach = GronwallReachabilityAnalyzer(dyn)
    eta_min = args.epsilon / (2 * dyn.L_l)
    write_params(out_dir, args, dyn, {"reachability": reach.describe(), "algorithm": 2, "eta_min": eta_min,
                                       "resumed_from": resume_meta and resume_meta.get("phase")})
    vi = SafetyValueIterator(mode, dyn, args.gamma, tree, reach, out_dir, n_workers=args.workers,
                             plotter=make_plotter(args, mode, dyn), dump_csv=args.dump_csv)
    if args.precompute and resume_meta is None:
        vi.precompute_all_successors()
    adaptive = AdaptiveRefinement(vi, tree, out_dir, args)
    try:
        t0 = time.time()
        adaptive.refine(epsilon=args.epsilon, max_refinements=args.refinements,
                        vi_iterations_per_refinement=args.vi_iterations, resume=resume_meta)
        print("ALGORITHM 2 COMPLETE")
        print(f"Total time: {time.time() - t0:.2f} seconds")
    finally:
        vi.close()


if __name__ == "__main__":
    sys.exit(main())
