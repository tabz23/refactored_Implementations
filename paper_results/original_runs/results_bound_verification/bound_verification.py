#!/usr/bin/env python3
"""Numerical check of the convergence theorem (label thm:convergence) on the Dubins car.

For every discount factor gamma with gamma * L_f < 1 and every resolution N, Algorithm 1
runs on a uniform grid of N^3 cells over a periodic Dubins domain: x1 and x2 wrap around
[-3, 3) and the heading wraps around [-pi, pi). The dynamics then map the domain into
itself, as the paper assumes (f : X x U -> X), and l and r use the periodic distance to
the obstacle and the target, so both are sqrt(2)-Lipschitz in the infinity norm.

The abstraction is the one of the paper: the cell bounds of l and r are the center values
-/+ L eta, where eta is the cell radius and L = max(L_l, L_r), and the successors of a
cell under an action are the cells whose closed boxes meet the ball of infinity-norm
radius L_f eta around the exact image of the cell center, with L_f = exp(tau) (the vector
field is 1-Lipschitz in the infinity norm). Value iteration is Algorithm 1 with the
stopping tolerances delta_lower^* = delta_upper^* = delta^* = 1e-12.

Check. The bounds are sound: V_lower(s) <= V_gamma(x) <= V_upper(s) for x in [[s]], so
|V_upper(s) - V_gamma(x)| and |V_lower(s) - V_gamma(x)| are at most V_upper(s) - V_lower(s).
Hence, if max_s (V_upper(s) - V_lower(s)) <= eps(eta) := (2 L eta + gamma delta^*)/(1 - gamma),
the bound of the theorem for the outputs of Algorithm 1, both inequalities of the theorem
hold for every state of every cell. The script reports this
ratio for every (gamma, N), the empirical rate of the largest gap, and the fractions of
cells certified in R, in L, and unclassified.

Usage: python bound_verification.py --out DIR [--gammas ...] [--resolutions ...]
"""
import argparse
import itertools
import json
import math
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

V_CONST, TAU = 1.0, 0.3
L_F = math.exp(1.0 * TAU)      # Lipschitz constant of the sampled map (Gronwall)
L = math.sqrt(2.0)             # L_l = L_r: Euclidean distances in the plane, infinity norm
ACTIONS = (-1.0, 0.0, 1.0)     # A = U
LO = np.array([-3.0, -3.0, -math.pi])
HI = np.array([3.0, 3.0, math.pi])
OBSTACLE, R_OBSTACLE = (0.0, 0.0), 1.3
TARGET, R_TARGET = (2.5, 0.0), 0.5
TOL = 1e-12                    # stopping tolerance of value iteration
IDX_EPS = 1e-9                 # index ranges are widened by this fraction of a cell (sound)


def flow(x, y, th, u, t):
    """Exact Dubins state after time t under the constant turn rate u."""
    if u == 0.0:
        return x + V_CONST * t * np.cos(th), y + V_CONST * t * np.sin(th), th
    th1 = th + u * t
    return (x + (V_CONST / u) * (np.sin(th1) - np.sin(th)),
            y - (V_CONST / u) * (np.cos(th1) - np.cos(th)), th1)


def periodic_distance(x, y, c):
    dx = np.mod(x - c[0] + 3.0, 6.0) - 3.0
    dy = np.mod(y - c[1] + 3.0, 6.0) - 3.0
    return np.hypot(dx, dy)


def successor_rows(N, C, h, eta, u):
    """Flat indices of the successors of every cell under u, one row per offset in the
    largest box; rows past a cell's own box repeat one of its successors."""
    x1, y1, th1 = flow(C[0], C[1], C[2], u, TAU)
    e = L_F * eta
    lo_i, w = [], []
    for d, c in enumerate((x1, y1, th1)):
        # cell k covers [LO + k h, LO + (k + 1) h] and meets [c - e, c + e] iff
        # (c - e - LO)/h - 1 <= k <= (c + e - LO)/h
        i0 = np.ceil((c - e - LO[d]) / h[d] - 1.0 - IDX_EPS).astype(np.int64)
        i1 = np.floor((c + e - LO[d]) / h[d] + IDX_EPS).astype(np.int64)
        lo_i.append(i0)
        w.append(i1 - i0 + 1)
    rows = []
    for off in itertools.product(*(range(int(wd.max())) for wd in w)):
        j = [np.mod(lo_i[d] + np.minimum(off[d], w[d] - 1), N) for d in range(3)]
        rows.append(((j[0] * N + j[1]) * N + j[2]).ravel().astype(np.int32))
    return np.stack(rows)


def reduce(rows, V, op):
    """op (np.minimum or np.maximum) of V over the successors of every cell."""
    flat = V.ravel()
    out = flat[rows[0]].copy()
    for r in rows[1:]:
        op(out, flat[r], out=out)
    return out.reshape(V.shape)


def grid(N):
    h = (HI - LO) / N
    eta = 0.5 * float(h.max())
    axes = [LO[d] + (np.arange(N) + 0.5) * h[d] for d in range(3)]
    C = np.meshgrid(*axes, indexing="ij")
    lc = periodic_distance(C[0], C[1], OBSTACLE) - R_OBSTACLE
    rc = -(periodic_distance(C[0], C[1], TARGET) - R_TARGET)
    bounds = (lc - L * eta, lc + L * eta, rc - L * eta, rc + L * eta)
    tables = [successor_rows(N, C, h, eta, u) for u in ACTIONS]
    return eta, bounds, tables


def value_iteration(gamma, bounds, tables):
    """Algorithm 1 run to the fixed points of both discrete Bellman operators."""
    ll, lu, rl, ru = bounds
    Vl, Vu = np.minimum(ll, rl), lu.copy()
    for k in range(1, 1000000):
        Ml = np.max([reduce(t, Vl, np.minimum) for t in tables], axis=0)
        Mu = np.max([reduce(t, Vu, np.maximum) for t in tables], axis=0)
        nl = np.minimum(ll, np.maximum(rl, gamma * Ml))
        nu = np.minimum(lu, np.maximum(ru, gamma * Mu))
        dl, du = float(np.abs(nl - Vl).max()), float(np.abs(nu - Vu).max())
        Vl, Vu = nl, nu
        if dl <= TOL and du <= TOL:
            return Vl, Vu, k
    raise RuntimeError("value iteration did not converge")


CLASS_CMAP = ListedColormap(["#c62828", "#9e9e9e", "#2e7d32"])   # red, gray, green
# Okabe-Ito colorblind-safe palette and one marker per discount factor.
PALETTE = ["#000000", "#E69F00", "#0072B2", "#009E73", "#D55E00", "#CC79A7", "#56B4E9"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X"]
HEADINGS = [(0.0, r"$\theta=0$"), (math.pi, r"$\theta=\pi$"), (math.pi / 4, r"$\theta=\pi/4$"),
            (-math.pi / 4, r"$\theta=-\pi/4$"), (math.pi / 2, r"$\theta=\pi/2$"), (-math.pi / 2, r"$\theta=-\pi/2$")]


def plot_values(gamma, N, eta, Vl, Vu, out):
    """Lower bound, upper bound and classification of one run in six heading slices."""
    ext = (LO[0], HI[0], LO[1], HI[1])
    norm = TwoSlopeNorm(vcenter=0.0, vmin=min(float(Vl.min()), -1e-6), vmax=max(float(Vu.max()), 1e-6))
    fig, axes = plt.subplots(3, len(HEADINGS), figsize=(16.0, 8.2))
    for j, (th, label) in enumerate(HEADINGS):
        k = int(np.mod(np.floor((th - LO[2]) / ((HI[2] - LO[2]) / N)), N))
        lo, up = Vl[:, :, k].T, Vu[:, :, k].T      # rows: x2, columns: x1
        cls = np.where(lo > 0, 2, np.where(up <= 0, 0, 1))
        im = axes[0, j].imshow(lo, origin="lower", extent=ext, cmap="RdYlGn", norm=norm, interpolation="nearest")
        axes[1, j].imshow(up, origin="lower", extent=ext, cmap="RdYlGn", norm=norm, interpolation="nearest")
        axes[2, j].imshow(cls, origin="lower", extent=ext, cmap=CLASS_CMAP, vmin=0, vmax=2, interpolation="nearest")
        axes[0, j].set_title(label, fontsize=10)
    for row, name in enumerate(["lower bound", "upper bound", "classification"]):
        axes[row, 0].set_ylabel(f"{name}\n$x_2$")
    for ax in axes.ravel():
        ax.add_patch(Circle(OBSTACLE, R_OBSTACLE, fill=False, color="#0d47a1", lw=1.0))
        ax.add_patch(Circle(TARGET, R_TARGET, fill=False, color="#f9a825", lw=1.0, ls="--"))
        ax.set_xlim(LO[0], HI[0])
        ax.set_ylim(LO[1], HI[1])
        ax.set_aspect("equal")
    for ax in axes[2]:
        ax.set_xlabel("$x_1$")
    fig.colorbar(im, ax=axes[:2, :].ravel().tolist(), fraction=0.015, pad=0.01, label="value")
    fig.suptitle(f"Algorithm 1 on the periodic domain, γ = {gamma}, N = {N} (η = {eta:.4f}); "
                 "classification: green R (lower bound > 0), red L (upper bound ≤ 0), gray unclassified",
                 fontsize=10)
    d = out / "value_functions"
    d.mkdir(exist_ok=True)
    fig.savefig(d / f"gamma{gamma}_N{N}.png", dpi=110)
    plt.close(fig)


def paper_figure(rows, gammas, out):
    """Single-column figure: largest gap of Algorithm 1 against eta, with the bound of the theorem."""
    import os
    import shutil
    from matplotlib.ticker import FixedLocator, NullLocator

    os.environ["PATH"] = "/Library/TeX/texbin:" + os.environ.get("PATH", "")
    usetex = bool(shutil.which("latex") and shutil.which("dvipng"))
    rc = {"text.usetex": usetex, "font.family": "serif", "font.size": 9, "axes.labelsize": 10,
          "legend.fontsize": 8, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5}
    if usetex:
        rc["text.latex.preamble"] = r"\usepackage{amsmath}"
        ylabel = (r"\shortstack{$\max_{s\in S}\bigl(\overline V^{\mathrm{cons}}_\gamma(s)"
                  r"-\underline V^{\mathrm{cons}}_\gamma(s)\bigr)$ (solid)\\[3pt]"
                  r"bound $(2L\varepsilon+\gamma\delta^*)/(1-\gamma)$ (dashed)}")
    else:
        ylabel = "max_s (V_upper^cons(s) - V_lower^cons(s)) (solid)\nbound (2Lε+γδ*)/(1-γ) (dashed)"
    with plt.rc_context(rc):
        fig, ax = plt.subplots(figsize=(3.4, 2.9))
        for g, col, mk in zip(gammas, PALETTE, MARKERS):
            rs = [r for r in rows if r["gamma"] == g]
            eta = [r["eta"] for r in rs]
            ax.loglog(eta, [r["gap_max"] for r in rs], "-", marker=mk, color=col, ms=3.2, lw=1.0,
                      label=rf"$\gamma={g}$")
            ax.loglog(eta, [r["eps"] for r in rs], "--", color=col, lw=0.9)
        ax.xaxis.set_major_locator(FixedLocator([0.03, 0.04, 0.05, 0.06]))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_xticklabels(["0.03", "0.04", "0.05", "0.06"])
        ax.yaxis.set_major_locator(FixedLocator([0.1, 0.2, 0.3, 0.4, 0.6]))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.set_yticklabels(["0.1", "0.2", "0.3", "0.4", "0.6"])
        ax.set_xlabel(r"cell radius $\varepsilon=\max_i\eta_{x,i}(s)$")
        ax.set_ylabel(ylabel)
        ax.legend(ncol=2, loc="upper left", frameon=False, handlelength=1.6, columnspacing=1.0)
        fig.tight_layout()
        for ext in ("pdf", "png"):
            fig.savefig(out / f"bound_verification_paper.{ext}", dpi=300)
        plt.close(fig)


def slope(x, y):
    """Least-squares slope of log y against log x."""
    return float(np.polyfit(np.log(x), np.log(y), 1)[0])


def write_outputs(rows, gammas, out):
    rows = sorted(rows, key=lambda r: (r["gamma"], r["N"]))
    summary = {"L": L, "L_f": L_F, "tolerance": TOL, "runs": rows, "per_gamma": {}}
    for g in gammas:
        rs = [r for r in rows if r["gamma"] == g]
        summary["per_gamma"][str(g)] = {
            "gamma_L_f": g * L_F,
            "largest_ratio": max(r["ratio"] for r in rs),
            "all_verified": all(r["verified"] for r in rs),
            "rate_gap_max": slope([r["eta"] for r in rs], [r["gap_max"] for r in rs]),
            "rate_gap_mean": slope([r["eta"] for r in rs], [r["gap_mean"] for r in rs]),
        }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))

    # Figure: largest gap against eta with the bound, and the ratio to the bound.
    paper_figure(rows, gammas, out)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.7))
    for g, col in zip(gammas, PALETTE):
        rs = [r for r in rows if r["gamma"] == g]
        eta = np.array([r["eta"] for r in rs])
        a1.loglog(eta, [r["gap_max"] for r in rs], "o-", color=col, ms=3.5, lw=1, label=rf"$\gamma={g}$")
        a1.loglog(eta, [r["eps"] for r in rs], "--", color=col, lw=0.9)
        a2.semilogx(eta, [r["ratio"] for r in rs], "o-", color=col, ms=3.5, lw=1, label=rf"$\gamma={g}$")
    a2.axhline(1.0, color="#c62828", lw=1)
    for ax in (a1, a2):
        ax.set_xlabel(r"cell radius $\eta$")
        ax.set_xticks([0.03, 0.04, 0.05, 0.06])
        ax.set_xticklabels(["0.03", "0.04", "0.05", "0.06"])
        ax.minorticks_off()
    a1.set_ylabel(r"largest gap (solid), $(2L\eta+\gamma\delta^*)/(1-\gamma)$ (dashed)", fontsize=8)
    a2.set_ylabel("largest gap / bound")
    a2.set_ylim(0.0, 1.1)
    a1.legend(fontsize=6, ncol=1)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out / f"bound_verification.{ext}", dpi=300)
    plt.close(fig)

    # LaTeX table (booktabs): one row per gamma.
    lines = [r"\begin{tabular}{rrrrr}", r"\toprule",
             r"$\gamma$ & $\gamma L_f$ & $\max_N$ gap/bound & rate & unclassified (\%) \\", r"\midrule"]
    for g in gammas:
        rs = [r for r in rows if r["gamma"] == g]
        p = summary["per_gamma"][str(g)]
        lines.append(f"{g} & {p['gamma_L_f']:.2f} & {p['largest_ratio']:.2f} & {p['rate_gap_max']:.2f} & "
                     f"{100 * rs[0]['frac_unclassified']:.1f} $\\to$ {100 * rs[-1]['frac_unclassified']:.1f} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    (out / "table_bounds.tex").write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.3, 0.4, 0.5, 0.6, 0.7])
    ap.add_argument("--resolutions", type=int, nargs="+", default=[50, 60, 70, 80, 90, 100, 110, 120])
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for g in a.gammas:
        if not (0.0 < g < 1.0 and g * L_F < 1.0):
            raise ValueError(f"gamma = {g} violates gamma * L_f < 1 (L_f = {L_F:.4f})")
    rows = []
    for N in a.resolutions:
        t0 = time.time()
        eta, bounds, tables = grid(N)
        t_grid = time.time() - t0
        for g in a.gammas:
            t1 = time.time()
            Vl, Vu, iters = value_iteration(g, bounds, tables)
            gap = Vu - Vl
            eps = (2.0 * L * eta + g * TOL) / (1.0 - g)     # theorem's bound for Algorithm 1, delta^* = TOL
            row = {"gamma": g, "N": N, "eta": eta, "cells": int(gap.size), "iterations": iters,
                   "eps": eps, "gap_max": float(gap.max()), "gap_mean": float(gap.mean()),
                   "gap_min": float(gap.min()), "ratio": float(gap.max() / eps),
                   "verified": bool(gap.max() <= eps and gap.min() >= -1e-12),
                   "frac_R": float((Vl > 0).mean()), "frac_L": float((Vu <= 0).mean()),
                   "frac_unclassified": float(((Vl <= 0) & (Vu > 0)).mean()),
                   "seconds": round(t_grid + time.time() - t1, 1)}
            rows.append(row)
            print(json.dumps(row), flush=True)
            plot_values(g, N, eta, Vl, Vu, out)
        del tables
    write_outputs(rows, a.gammas, out)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
