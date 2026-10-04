#!/usr/bin/env python
"""Runtime figures for one Algorithm 2 run.

Reads phases.csv and run.log in a results directory. Iteration 1 is phase 0.
The time is the wall time of that iteration, including its plot. The second
figure adds the share of state-space volume classified safe or unsafe.

  python plot_runtime.py --run-dir ../results_hscc_2027/ra_nodiscount/dubins_res20_..._std
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def last_run(rows):
    starts = [i for i, r in enumerate(rows) if r["phase"] == "0"]
    if not starts:
        raise SystemExit("phases.csv has no phase 0")
    return rows[starts[-1] :]


def g_label(n):
    v = n / 1000.0
    text = f"{v:.1f}".rstrip("0").rstrip(".")
    return text + "k"


def load(run: Path):
    rows = last_run(list(csv.DictReader((run / "phases.csv").open())))
    iters = [int(r["phase"]) + 1 for r in rows]
    times = [float(r["time_phase_s"]) for r in rows]
    cells = [int(r["n_leaves"]) for r in rows]

    text = (run / "run.log").read_text()
    pat = re.compile(
        r"Classification after phase (\d+):\n"
        r"  Safe:\s+\d+ cells \(\s*[0-9.]+%\),\s*([0-9.]+)% of state space\n"
        r"  Unsafe:\s+\d+ cells \(\s*[0-9.]+%\),\s*([0-9.]+)% of state space"
    )
    blocks = list(pat.finditer(text))
    if not blocks:
        raise SystemExit("run.log has no state-space classification lines")
    last0 = max(i for i, m in enumerate(blocks) if m.group(1) == "0")
    volume = {}
    for m in blocks[last0:]:
        volume[int(m.group(1))] = float(m.group(2)) + float(m.group(3))
    classified = [volume[i - 1] for i in iters]
    return iters, times, cells, classified


def style_axes(ax):
    ax.tick_params(axis="both", which="major", labelsize=16, width=1.2, length=6)
    ax.tick_params(axis="both", which="minor", width=0.8, length=3)
    ax.set_xlabel("Refinement Iteration", fontsize=20)
    ax.set_ylabel("Time (s)", fontsize=20)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(16)
    ax.xaxis.label.set_size(20)
    ax.yaxis.label.set_size(20)
    ax.grid(True, which="both", color="#b0b0b0", linewidth=0.7)
    ax.set_axisbelow(True)


def draw(run: Path, iters, times, cells, classified, with_class: bool):
    # The classified figure is 25% shorter so it fits the paper column.
    height = 4.2 if with_class else 5.6
    fig, ax = plt.subplots(figsize=(10.2, height), dpi=160)
    ax.plot(iters, times, color="#1f77b4", marker="o", markersize=7, linewidth=2.0)
    ax.set_yscale("log")
    style_axes(ax)
    ax.set_xticks(iters)
    ax.set_xlim(min(iters) - 0.5, max(iters) + 0.5)
    ax.tick_params(axis="y", labelsize=16)
    ymin, ymax = ax.get_ylim()
    number_labels = []
    for x, n in zip(iters, cells):
        number_labels.append(
            ax.text(x, ymax * 1.02, g_label(n), rotation=55, ha="left", va="bottom", fontsize=13, clip_on=False)
        )

    if with_class:
        ax2 = ax.twinx()
        class_label = "% of the full\nstate-space volume\nclassified as safe/unsafe"
        ax2.plot(iters, classified, color="#d62728", marker="s", markersize=5.5, linewidth=1.8)
        ax2.set_ylabel(class_label, fontsize=20, color="#d62728", labelpad=12)
        ax2.set_ylim(0, 100)
        ax2.tick_params(axis="y", labelsize=16, colors="#d62728", width=1.2, length=6)
        ax.set_ylabel("Time (s)", fontsize=20, color="#1f77b4")
        ax.tick_params(axis="y", labelsize=16, colors="#1f77b4", width=1.2, length=6)
        name = "runtime_algorithm2_classified.png"
    else:
        name = "runtime_algorithm2.png"

    fig.tight_layout()
    # Header only needs the rotated cell counts. |G| sits at their left.
    fig.subplots_adjust(top=0.80, right=0.80 if with_class else 0.97)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axbb = ax.get_window_extent(renderer)
    g_label_artist = fig.text(0, 0, r"$|\mathcal{G}|$", ha="right", va="bottom", fontsize=28)
    fig.canvas.draw()
    # Just left of the axes and just above the top spine, beside the upper tick.
    g_label_artist.set_position((
        (axbb.x0 - 2) / fig.bbox.width,
        (axbb.y1 + 8) / fig.bbox.height,
    ))
    path = run / name
    fig.savefig(path, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    print("wrote", path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", type=Path, required=True, help="directory with phases.csv and run.log")
    args = p.parse_args()
    run = args.run_dir.resolve()
    iters, times, cells, classified = load(run)
    draw(run, iters, times, cells, classified, False)
    draw(run, iters, times, cells, classified, True)
    for i, t, n, c in zip(iters, times, cells, classified):
        print(f"  iter {i:2d}  |G|={n:8d}  time={t:.3f}s  state space classified={c:.1f}%")


if __name__ == "__main__":
    main()
