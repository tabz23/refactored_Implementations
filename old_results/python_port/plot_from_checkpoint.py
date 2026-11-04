#!/usr/bin/env python
"""Render matplotlib figures from a C++ (or Python) checkpoint.

  python plot_from_checkpoint.py --run-dir <results/.../run>
  python plot_from_checkpoint.py --checkpoint path/to/checkpoint_phase_0003.bin \\
      --dynamics dubins --mode ra_nodiscount --out phase3.png
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import struct
import sys

import numpy as np

from dhj.core.modes import Mode
from dhj.dynamics import REGISTRY, make_dynamics
from dhj.dynamics.base import Dynamics
from dhj.dynamics.plotting import PlotStyle, _draw_overlays, _finish_axes, _safe_ticks, cell_in_slice


class _Cell:
    __slots__ = ("bounds", "V_lower", "V_upper")

    def __init__(self, bounds, vl, vu):
        self.bounds = bounds
        self.V_lower = vl
        self.V_upper = vu

    def get_range(self, d):
        return self.bounds[d, 1] - self.bounds[d, 0]


def load_cpp_checkpoint(path: str):
    """Return (phase, dim, mode_int, list[_Cell] of leaves)."""
    with open(path, "rb") as f:
        if f.read(8) != b"DHJCKPT1":
            raise RuntimeError(f"{path}: not a DHJ checkpoint")
        version, dim, mode, phase = struct.unpack("<IIII", f.read(16))
        if version != 1:
            raise RuntimeError(f"{path}: unsupported version {version}")
        f.read(8)  # total_refined
        f.read(8)  # init_res + pad
        f.read(8)  # epsilon
        f.read(16 * dim)  # root bounds
        num_cells, num_leaves = struct.unpack("<QQ", f.read(16))
        rec = 16 * dim + 56
        raw = f.read(num_cells * rec)
        if len(raw) != num_cells * rec:
            raise RuntimeError(f"{path}: truncated cell records")
        leaf_ids = np.frombuffer(f.read(num_leaves * 4), dtype="<u4")
        if leaf_ids.size != num_leaves:
            raise RuntimeError(f"{path}: truncated leaf list")

    dt = np.dtype([
        ("bounds", "<f8", (dim, 2)),
        ("parent", "<i4"),
        ("child0", "<i4"),
        ("Vu", "<f8"), ("Vl", "<f8"),
        ("lu", "<f8"), ("ll", "<f8"),
        ("ru", "<f8"), ("rl", "<f8"),
    ])
    recs = np.frombuffer(raw, dtype=dt)
    leaves = recs[leaf_ids]
    cells = [_Cell(leaves["bounds"][i].copy(), float(leaves["Vl"][i]), float(leaves["Vu"][i]))
             for i in range(len(leaves))]
    return phase, dim, mode, cells


def _quads(cells, px, py):
    verts = np.empty((len(cells), 4, 2))
    for i, c in enumerate(cells):
        ax, bx = c.bounds[px]
        ay, by = c.bounds[py]
        verts[i] = ((ax, ay), (bx, ay), (bx, by), (ax, by))
    return verts


def plot_cells(dyn: Dynamics, cells, filename: str, iteration: int, style: PlotStyle):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize
    from matplotlib.ticker import FuncFormatter

    slices = dyn.slices()
    n = len(slices)
    fig, axes = plt.subplots(3, n, figsize=(5 * n, 14), dpi=style.dpi)
    if n == 1:
        axes = axes.reshape(3, 1)
    px, py = dyn.plot_dims
    cmap = plt.cm.RdYlGn
    C_UNSAFE, C_SAFE, C_BOUND = "#d62728", "#2ca02c", "#7f7f7f"

    for idx, sl in enumerate(slices):
        relevant = [c for c in cells if cell_in_slice(c.bounds, sl, dyn.periodic_dims)]
        verts = _quads(relevant, px, py) if relevant else np.zeros((0, 4, 2))

        for row, upper, label, title in (
            (0, True, "V̄_γ", f"Upper Bound V̄_γ ({sl.label})"),
            (1, False, "V_γ", f"Lower Bound V_γ ({sl.label})"),
        ):
            ax = axes[row, idx]
            if not relevant:
                ax.text(0.5, 0.5, "No data for this slice", transform=ax.transAxes,
                        ha="center", va="center")
            elif style.draw_value_cells:
                values = np.array([c.V_upper if upper else c.V_lower for c in relevant])
                vmin, vmax = float(values.min()), float(values.max())
                if np.isclose(vmin, vmax):
                    vmin, vmax = vmin - 1e-12, vmax + 1e-12
                if vmax <= 0:
                    norm = Normalize(vmin=vmin, vmax=0)
                elif vmin >= 0:
                    norm = Normalize(vmin=0, vmax=vmax)
                else:
                    norm = Normalize(vmin=vmin, vmax=vmax)
                colors = cmap(norm(values))
                lw = 0.005 if style.save_dpi >= 800 else 0.0
                ax.add_collection(PolyCollection(verts, facecolors=colors, edgecolors="black",
                                                 linewidths=lw, antialiased=False))
                sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
                sm.set_array([])
                cbar = plt.colorbar(sm, ax=ax, label=label)
                cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.1e}"))
                cbar.ax.yaxis.get_offset_text().set_visible(False)
                ticks = _safe_ticks(vmin, vmax)
                cbar.set_ticks(ticks)
                cbar.ax.set_yticklabels([f"{v:.1e}" for v in ticks])
            _draw_overlays(ax, dyn, style)
            _finish_axes(ax, dyn)
            ax.set_title(title)

        ax = axes[2, idx]
        if relevant:
            colors = []
            for c in relevant:
                if c.V_lower > 0:
                    colors.append(C_SAFE)
                elif c.V_upper <= 0:
                    colors.append(C_UNSAFE)
                else:
                    colors.append(C_BOUND)
            lw = 0.005 if style.save_dpi >= 800 else 0.0
            ax.add_collection(PolyCollection(verts, facecolors=colors, edgecolors="black",
                                             linewidths=lw, antialiased=False))
        _draw_overlays(ax, dyn, style)
        _finish_axes(ax, dyn)
        ax.set_title(f"Cell Classification ({sl.label})")

    fig.suptitle(f"Safety Value Function - Iteration {iteration}", fontsize=16, y=0.995)
    plt.tight_layout()
    plt.savefig(filename, dpi=style.save_dpi, bbox_inches="tight")
    plt.close(fig)


def _style_for_n(n: int, mode: Mode, plot_dpi: int = 800) -> PlotStyle:
    # 800 keeps the paper scripts (figure dpi 600, savefig dpi 800). A smaller
    # value rasterizes at that dpi. At a low dpi a 0.005 pt edge snaps to one
    # pixel and covers the small cells, so the edge is omitted below 800.
    del n
    if plot_dpi == 800:
        return PlotStyle(draw_value_cells=True, cell_labels=mode.discounted,
                         draw_target=mode.has_target, dpi=600, save_dpi=800)
    return PlotStyle(draw_value_cells=True, cell_labels=mode.discounted,
                     draw_target=mode.has_target, dpi=plot_dpi, save_dpi=plot_dpi)


class _Args:
    def __init__(self, dt, tau, velocity=1.0, a_max=1.0, v_max=1.0):
        self.dt, self.tau = dt, tau
        self.velocity, self.a_max, self.v_max = velocity, a_max, v_max


def plot_one(ckpt: str, dyn: Dynamics, mode: Mode, out: str, iteration: int = None, plot_dpi: int = 800):
    phase, dim, _mode_i, cells = load_cpp_checkpoint(ckpt)
    if dim != dyn.dim:
        raise SystemExit(f"{ckpt}: dim {dim} != dynamics dim {dyn.dim}")
    iteration = phase if iteration is None else iteration
    style = _style_for_n(len(cells), mode, plot_dpi)
    plot_cells(dyn, cells, out, iteration, style)
    print(f"wrote {out}  (phase {phase}, {len(cells)} leaves)", flush=True)


def plot_run_dir(run_dir: str, plot_dpi: int = 800):
    params_path = os.path.join(run_dir, "params.json")
    with open(params_path, encoding="utf-8") as f:
        params = json.load(f)
    args = params["args"]
    mode = Mode(args["mode"])
    dyn = make_dynamics(args["dynamics"], _Args(args["dt"], args["tau"],
                                               args.get("velocity", 1.0),
                                               args.get("a_max", 1.0),
                                               args.get("v_max", 1.0)))
    ckpts = sorted(glob.glob(os.path.join(run_dir, "checkpoints", "checkpoint_phase_*.bin")))
    if not ckpts:
        raise SystemExit(f"no checkpoint_phase_*.bin in {run_dir}/checkpoints")
    for ckpt in ckpts:
        phase = int(os.path.basename(ckpt).split("_")[-1].split(".")[0])
        out = os.path.join(run_dir, f"value_function_phase_{phase}_complete.png")
        if os.path.exists(out) and os.path.getsize(out) > 0:
            print(f"skip {out}", flush=True)
            continue
        plot_one(ckpt, dyn, mode, out, iteration=phase, plot_dpi=plot_dpi)


def main(argv=None):
    p = argparse.ArgumentParser(description="Plot C++ DHJ checkpoints as matplotlib PNGs")
    p.add_argument("--run-dir", help="results directory containing params.json and checkpoints/")
    p.add_argument("--results-root", help="plot every run under this tree")
    p.add_argument("--checkpoint")
    p.add_argument("--dynamics", choices=sorted(REGISTRY))
    p.add_argument("--mode", default=Mode.RA_NODISCOUNT.value, choices=[m.value for m in Mode])
    p.add_argument("--out")
    p.add_argument("--dt", type=float, default=0.3)
    p.add_argument("--tau", type=float, default=0.3)
    p.add_argument("--plot-dpi", type=int, default=800,
                   help="800 matches the paper scripts. A smaller value writes a smaller PNG")
    args = p.parse_args(argv)


    if args.results_root:
        for dirpath, dirnames, filenames in os.walk(args.results_root):
            if "params.json" in filenames and os.path.isdir(os.path.join(dirpath, "checkpoints")):
                print(f"=== {dirpath}", flush=True)
                plot_run_dir(dirpath, args.plot_dpi)
        return 0
    if args.run_dir:
        plot_run_dir(args.run_dir, args.plot_dpi)
        return 0
    if not args.checkpoint or not args.dynamics:
        p.error("pass --results-root, --run-dir, or --checkpoint plus --dynamics")
    dyn = make_dynamics(args.dynamics, _Args(args.dt, args.tau))
    out = args.out or os.path.splitext(args.checkpoint)[0] + ".png"
    plot_one(args.checkpoint, dyn, Mode(args.mode), out, plot_dpi=args.plot_dpi)
    return 0


if __name__ == "__main__":
    sys.exit(main())
