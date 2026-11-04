"""Generic slice plotter shared by all dynamics.

Reproduces the figure of the paper scripts: for each slice a column with
three panels (upper bound, lower bound, cell classification) in which every
leaf cell intersecting the slice is drawn as a coloured rectangle on the
(plot_dims) plane. Periodic dimensions (angles) are handled by also matching
the slice value shifted by +-2*pi, so a request for theta = 3*pi/2 finds the
cells around -pi/2 and theta = pi finds both the cells ending at pi and those
starting at -pi.
"""
import math
from typing import Iterable, List, Tuple

import numpy as np

from .base import Dynamics, Slice


class PlotStyle:
    def __init__(self, draw_value_cells: bool = True, cell_labels: bool = False,
                 draw_target: bool = True, dpi: int = 600, save_dpi: int = 800):
        self.draw_value_cells = draw_value_cells
        self.cell_labels = cell_labels
        self.draw_target = draw_target
        self.dpi = dpi
        self.save_dpi = save_dpi


def _candidates(value: float, period: float) -> Tuple[float, ...]:
    return (value, value - period, value + period)


def cell_in_slice(bounds: np.ndarray, sl: Slice, periodic_dims: Iterable[int]) -> bool:
    """True if the cell box contains every pinned coordinate of the slice
    (inclusive, like the paper scripts), with wrap-around on periodic dims."""
    periodic = set(periodic_dims)
    for d, v in sl.fixed.items():
        lo, hi = bounds[d, 0], bounds[d, 1]
        if d in periodic:
            if not any(lo <= c <= hi for c in _candidates(v, 2.0 * math.pi)):
                return False
        else:
            if not (lo <= v <= hi):
                return False
    return True


def _safe_ticks(vmin, vmax, zero_threshold=0.1):
    ticks = [vmin, vmax]
    if abs(vmin) >= zero_threshold and abs(vmax) >= zero_threshold:
        ticks.insert(1, 0.0)
    return sorted(set(ticks))


def _draw_overlays(ax, dyn: Dynamics, style: PlotStyle):
    from matplotlib.patches import Circle  # noqa: PLC0415

    for ov in dyn.overlays():
        if ov.when == "reach_avoid" and not style.draw_target:
            continue
        if ov.kind == "circle":
            ax.add_patch(Circle(ov.center, ov.radius, facecolor="none", edgecolor=ov.color,
                                linewidth=ov.linewidth, zorder=10, linestyle=ov.linestyle))


def _finish_axes(ax, dyn: Dynamics):
    b = dyn.get_state_bounds()
    px, py = dyn.plot_dims
    ax.set_xlabel(dyn.state_name(px), fontsize=16)
    ax.set_ylabel(dyn.state_name(py), fontsize=16)
    ax.tick_params(axis="both", which="major", labelsize=14)
    ax.set_xlim(b[px, 0], b[px, 1])
    ax.set_ylim(b[py, 0], b[py, 1])
    ax.set_aspect("equal")
    ax.grid(False)


def _plot_value_panel(ax, dyn: Dynamics, cells: List, sl: Slice, upper: bool, label: str,
                      style: PlotStyle):
    import matplotlib.pyplot as plt  # noqa: PLC0415
    from matplotlib.colors import Normalize  # noqa: PLC0415
    from matplotlib.patches import Rectangle  # noqa: PLC0415
    from matplotlib.ticker import FuncFormatter  # noqa: PLC0415

    px, py = dyn.plot_dims
    b = dyn.get_state_bounds()
    relevant = [(c, c.V_upper if upper else c.V_lower) for c in cells
                if cell_in_slice(c.bounds, sl, dyn.periodic_dims)
                and (c.V_upper if upper else c.V_lower) is not None]
    if not relevant:
        ax.text(0.5, 0.5, "No data for this slice", transform=ax.transAxes, ha="center", va="center")
        ax.set_xlim(b[px, 0], b[px, 1])
        ax.set_ylim(b[py, 0], b[py, 1])
        return

    values = np.array([v for _, v in relevant])
    cmap = plt.cm.RdYlGn
    vmin, vmax = float(np.min(values)), float(np.max(values))
    if np.isclose(vmin, vmax):
        vmin, vmax = vmin - 1e-12, vmax + 1e-12
    if vmax <= 0:
        norm = Normalize(vmin=vmin, vmax=0)
    elif vmin >= 0:
        norm = Normalize(vmin=0, vmax=vmax)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)

    if style.draw_value_cells:
        for cell, value in relevant:
            a_x, b_x = cell.bounds[px]
            a_y, b_y = cell.bounds[py]
            edge = 0.005 if style.save_dpi >= 800 else 0.0
            ax.add_patch(Rectangle((a_x, a_y), b_x - a_x, b_y - a_y, facecolor=cmap(norm(value)),
                                   edgecolor="black", linewidth=edge, antialiased=False, alpha=1.0))
            if style.cell_labels:
                w, h = b_x - a_x, b_y - a_y
                if w > 0.05 and h > 0.05:
                    ax.text(0.5 * (a_x + b_x), 0.5 * (a_y + b_y), f"{value:.1e}", ha="center",
                            va="center", fontsize=min(w * 15, h * 15, 6),
                            color="black" if norm(value) > 0.5 else "white", clip_on=True)

    _draw_overlays(ax, dyn, style)
    _finish_axes(ax, dyn)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, label=label)
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.1e}"))
    cbar.ax.yaxis.get_offset_text().set_visible(False)
    ticks = _safe_ticks(vmin, vmax)
    cbar.set_ticks(ticks)
    cbar.ax.set_yticklabels([f"{v:.1e}" for v in ticks])


def _plot_classification_panel(ax, dyn: Dynamics, cells: List, sl: Slice, style: PlotStyle):
    from matplotlib.patches import Rectangle  # noqa: PLC0415

    C_UNSAFE, C_SAFE, C_BOUND, C_UNKNOWN = "#d62728", "#2ca02c", "#7f7f7f", "#ffffff"
    px, py = dyn.plot_dims
    for cell in cells:
        if not cell_in_slice(cell.bounds, sl, dyn.periodic_dims):
            continue
        if cell.V_upper is None or cell.V_lower is None:
            color = C_UNKNOWN
        elif cell.V_lower > 0:
            color = C_SAFE
        elif cell.V_upper <= 0:
            color = C_UNSAFE
        else:
            color = C_BOUND
        edge = 0.005 if style.save_dpi >= 800 else 0.0
        ax.add_patch(Rectangle((cell.bounds[px, 0], cell.bounds[py, 0]), cell.get_range(px),
                               cell.get_range(py), facecolor=color, edgecolor="black",
                               linewidth=edge, antialiased=False, alpha=1.0))
    _draw_overlays(ax, dyn, style)
    _finish_axes(ax, dyn)


def plot_value_function(dyn: Dynamics, cells: List, filename: str, iteration: int,
                        style: PlotStyle = None, slices: List[Slice] = None):
    """Draw the 3 x n_slices figure for `cells` (any objects with .bounds,
    .V_upper, .V_lower, .get_range) and save it to `filename`."""
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415

    style = style or PlotStyle()
    slices = slices or dyn.slices()
    n = len(slices)
    fig, axes = plt.subplots(3, n, figsize=(5 * n, 14), dpi=style.dpi)
    if n == 1:
        axes = axes.reshape(3, 1)
    for idx, sl in enumerate(slices):
        _plot_value_panel(axes[0, idx], dyn, cells, sl, True, "V̄_γ", style)
        axes[0, idx].set_title(f"Upper Bound V̄_γ ({sl.label})")
        _plot_value_panel(axes[1, idx], dyn, cells, sl, False, "V_γ", style)
        axes[1, idx].set_title(f"Lower Bound V_γ ({sl.label})")
        _plot_classification_panel(axes[2, idx], dyn, cells, sl, style)
        axes[2, idx].set_title(f"Cell Classification ({sl.label})")
    fig.suptitle(f"Safety Value Function - Iteration {iteration}", fontsize=16, y=0.995)
    plt.tight_layout()
    plt.savefig(filename, dpi=style.save_dpi, bbox_inches="tight")
    plt.close(fig)
