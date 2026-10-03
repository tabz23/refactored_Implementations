#!/usr/bin/env python
"""Render a value_function_phase_N.csv (written by the Python run with
--dump-csv, or by the C++ binary) with the matplotlib plotter of a dynamics.

  python plot_from_csv.py --dynamics dubins --mode ra_nodiscount \
      path/to/value_function_phase_3.csv --out phase3.png --iteration 3
"""
import argparse
import os
import sys

import numpy as np

from dhj.core.modes import Mode
from dhj.dynamics import REGISTRY, make_dynamics
from dhj.dynamics.plotting import PlotStyle, plot_value_function


class _CsvCell:
    __slots__ = ("bounds", "V_lower", "V_upper")

    def __init__(self, bounds, vl, vu):
        self.bounds = bounds
        self.V_lower = vl
        self.V_upper = vu

    def get_range(self, d):
        return self.bounds[d, 1] - self.bounds[d, 0]


def load_cells(path: str, dim: int):
    with open(path, encoding="utf-8") as f:
        header = f.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    col = {name: i for i, name in enumerate(header)}
    lo = [col[f"d{k}_lo"] for k in range(dim)]
    hi = [col[f"d{k}_hi"] for k in range(dim)]
    cells = []
    for row in data:
        b = np.empty((dim, 2))
        for k in range(dim):
            b[k, 0], b[k, 1] = row[lo[k]], row[hi[k]]
        cells.append(_CsvCell(b, float(row[col["V_lower"]]), float(row[col["V_upper"]])))
    return cells


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("csv")
    p.add_argument("--dynamics", required=True, choices=sorted(REGISTRY))
    p.add_argument("--mode", default=Mode.RA_NODISCOUNT.value, choices=[m.value for m in Mode])
    p.add_argument("--out", default=None, help="output image (default: <csv>.png)")
    p.add_argument("--iteration", type=int, default=0, help="number shown in the figure title")
    p.add_argument("--dt", type=float, default=0.3)
    p.add_argument("--tau", type=float, default=0.3)
    p.add_argument("--velocity", type=float, default=1.0)
    p.add_argument("--a-max", type=float, default=1.0)
    p.add_argument("--v-max", type=float, default=1.0)
    args = p.parse_args(argv)

    mode = Mode(args.mode)
    dyn = make_dynamics(args.dynamics, args)
    cells = load_cells(args.csv, dyn.dim)
    out = args.out or os.path.splitext(args.csv)[0] + ".png"
    style = PlotStyle(draw_value_cells=True, cell_labels=mode.discounted, draw_target=mode.has_target)
    plot_value_function(dyn, cells, out, args.iteration, style)
    print(f"wrote {out} ({len(cells)} cells)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
