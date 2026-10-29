#!/usr/bin/env python3
"""Check the guarantees of the paper on Algorithm 1 runs over a uniform grid.

Result numbers follow the paper as of 2026-10-28: Proposition 4 is labelled
prop:greedy and Proposition 6 is labelled prop:avoid.

Usage:
    python verify_quantitative.py RUN_DIR [RUN_DIR ...] --out OUT_DIR

Each RUN_DIR is a results directory written by build/dhj_head with --dump-csv and
--refinements 0, and all of them use the same mode: ra_discount (Proposition 4
and the classification, checks below) or avoid_nodiscount (Proposition 6, see
verify_avoid). Only the dumped grid, values and l/r bounds are read. The
abstraction is recomputed independently of the solver:

  * the exact flow of each cell center over tau under each constant action
    (closed form for Dubins and evasion, checked against solve_ivp);
  * a reach box of infinity-norm radius r exp(L_f tau), where r is the largest
    half-width of the cell (Gronwall);
  * successors: every cell whose closed box meets the reach box (theta is
    periodic), plus the out-of-bounds sink with value -1 if the box leaves the
    state box along x or y.

Checks for ra_discount:
  1. the cell bounds l_lower, l_upper, r_lower, r_upper, recomputed from the centers;
  2. lower-bound certificate: (T_hat V_lower)(s) >= V_lower(s) for every cell, where
     T_hat uses the greedy action of Proposition 4. For cells with V_lower > 0 this
     is the inequality the proof of Proposition 4 starts from;
  3. upper-bound certificate: (T_bar V_upper)(s) <= V_upper(s) for every cell. With
     A = U this gives V_gamma <= V_upper, so the cells with V_upper <= 0 are sound;
  4. closed-loop simulation with the exact flow from random states in every cell
     with V_lower > 0, under the greedy policy: the trajectory must reach T within
     ceil(log(V_lower(s0)/r*)/log(gamma)) steps, never enter F or leave the grid,
     keep gamma^t l(x_t) >= V_lower(s0) up to the reach time, and achieve a
     measure max_t g_gamma(x, t) between V_lower(s0) and V_upper(s0).
"""
import argparse
import csv
import itertools
import json
import math
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap, LogNorm  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

TOL = 1e-12      # absolute tolerance for floating-point comparisons
IDX_EPS = 1e-9   # index ranges are widened by this fraction of a cell (keeps the abstraction sound)
SINK = -1.0      # value of the out-of-bounds sink, as in the solver
PAPER_SLICE = {"dubins": 0.0, "evasion": math.pi}
AVOID_PATH_STEPS = 60   # steps drawn of each example avoid trajectory (all steps are checked)


# ------------------------------------------------------------------- dynamics
def wrap(th):
    return np.arctan2(np.sin(th), np.cos(th))


def flow(name, v, x, y, th, u, t):
    """Exact state after time t under the constant control u (theta not wrapped)."""
    if name == "dubins":
        if u == 0.0:
            return x + v * t * np.cos(th), y + v * t * np.sin(th), th + 0.0 * (x + t)
        th1 = th + u * t
        return x + (v / u) * (np.sin(th1) - np.sin(th)), y - (v / u) * (np.cos(th1) - np.cos(th)), th1
    if name == "evasion":
        # z = x + iy satisfies z' = -iu z + v (exp(i theta) - 1), theta(t) = theta0 - u t.
        z = x + 1j * y
        if u == 0.0:
            z1 = z + v * t * (np.exp(1j * th) - 1.0)
            return z1.real, z1.imag, th + 0.0 * (x + t)
        e = np.exp(-1j * u * t)
        th1 = th - u * t
        z1 = e * z + v * t * np.exp(1j * th1) - (v / (1j * u)) * (1.0 - e)
        return z1.real, z1.imag, th1
    raise ValueError(f"no closed-form flow for {name}")


def vector_field(name, v, u, s):
    x, y, th = s
    if name == "dubins":
        return [v * np.cos(th), v * np.sin(th), u]
    return [-v + v * np.cos(th) + u * y, v * np.sin(th) - u * x, -u]


def check_flow(name, v, actions, tau, rng):
    """Largest difference between the closed form and a tight numerical solution."""
    from scipy.integrate import solve_ivp

    err = 0.0
    for u in actions:
        for _ in range(25):
            s0 = rng.uniform([-3.0, -3.0, -math.pi], [3.0, 3.0, math.pi])
            sol = solve_ivp(lambda t, s: vector_field(name, v, u, s), (0.0, tau), s0,
                            method="DOP853", rtol=1e-13, atol=1e-13)
            exact = np.array(flow(name, v, s0[0], s0[1], s0[2], u, tau), dtype=float)
            err = max(err, float(np.max(np.abs(sol.y[:, -1] - exact))))
    return err


def l_fun(dyn, x, y):
    ox, oy = dyn["obstacle_position"]
    return np.hypot(x - ox, y - oy) - dyn["obstacle_radius"]


def r_fun(dyn, x, y):
    tx, ty = dyn["target_position"]
    return -(np.hypot(x - tx, y - ty) - dyn["target_radius"])


# ----------------------------------------------------------------------- grid
class Grid:
    """Uniform grid rebuilt from the CSV dump; fields are arrays indexed by (i, j, k)."""

    def __init__(self, dyn, data):
        bounds = np.array(dyn["state_bounds"], dtype=float)
        D = bounds.shape[0]
        if D != 3:
            raise ValueError("only 3-D systems are supported")
        self.lo, self.hi = bounds[:, 0], bounds[:, 1]
        clo = data[:, 1:1 + 2 * D:2]
        chi = data[:, 2:2 + 2 * D:2]
        width = chi - clo
        self.n = np.rint((self.hi - self.lo) / width[0]).astype(int)
        self.h = (self.hi - self.lo) / self.n
        if not np.allclose(width, self.h, rtol=0.0, atol=1e-12):
            raise ValueError("the grid is not uniform (was the run refined?)")
        idx = np.rint((clo - self.lo) / self.h).astype(int)
        if (idx < 0).any() or (idx >= self.n).any():
            raise ValueError("cell index out of range")
        flat = np.ravel_multi_index(idx.T, tuple(self.n))
        if not (np.unique(flat).size == flat.size == int(np.prod(self.n))):
            raise ValueError("cells are missing or repeated")
        shape = tuple(self.n)

        def field(col):
            a = np.empty(shape)
            a.reshape(-1)[flat] = data[:, col]
            return a

        c0 = 1 + 2 * D
        self.Vl, self.Vu, self.ll, self.lu, self.rl, self.ru = (field(c0 + k) for k in range(6))
        axes = [self.lo[d] + (np.arange(self.n[d]) + 0.5) * self.h[d] for d in range(D)]
        self.axes = axes
        self.C = np.meshgrid(*axes, indexing="ij")
        self.r_half = 0.5 * float(self.h.max())
        self.periodic = set(dyn["periodic_dims"])

    def locate(self, pts):
        """Cell indices of points (theta wrapped); x, y must lie inside the box."""
        return [np.clip(np.floor((pts[d] - self.lo[d]) / self.h[d]).astype(int), 0, self.n[d] - 1) for d in range(3)]


def load_run(run_dir):
    run_dir = Path(run_dir)
    params = json.loads((run_dir / "params.json").read_text())
    csv_path = run_dir / "value_function_phase_0.csv"
    cache = run_dir / "value_table.npy"
    if cache.exists() and cache.stat().st_mtime >= csv_path.stat().st_mtime:
        data = np.load(cache)
    else:
        data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
        np.save(cache, data)
    phases = list(csv.DictReader((run_dir / "phases.csv").open()))
    return params, data, phases


# ------------------------------------------------------------------ abstraction
def successor_table(grid, dyn, u, tau):
    """Successors of every cell under the constant action u.

    "idx" holds flat cell indices, one row per offset in the largest reach box, so
    that a reduction over the rows is the min or max over the successors: in rows
    past a cell's own box the index is clamped to the box, which repeats one of its
    successors, and a cell whose box misses the grid along x or y points to the
    sink slot (index = number of cells). "oob" flags the boxes that leave the
    state box along x or y.
    """
    X, Y, TH = grid.C
    x1, y1, th1 = flow(dyn["name"], dyn["v_const"], X, Y, TH, u, tau)
    th1 = wrap(th1)
    e = grid.r_half * math.exp(dyn["L_f"] * tau)
    lo_i, w = [], []
    oob = np.zeros(X.shape, dtype=bool)
    for d, c in enumerate((x1, y1, th1)):
        a, b = c - e, c + e
        i0 = np.ceil((a - grid.lo[d]) / grid.h[d] - 1.0 - IDX_EPS).astype(int)
        i1 = np.floor((b - grid.lo[d]) / grid.h[d] + IDX_EPS).astype(int)
        if d not in grid.periodic:
            oob |= (a < grid.lo[d]) | (b > grid.hi[d])
            # Intersect with the grid; a box that misses the grid along d meets no cell.
            i0 = np.maximum(i0, 0)
            i1 = np.minimum(i1, grid.n[d] - 1)
        lo_i.append(i0)
        w.append(np.maximum(i1 - i0 + 1, 0))
    kmax = [int(wd.max()) for wd in w]
    empty = (w[0] == 0) | (w[1] == 0) | (w[2] == 0)
    rows = []
    for off in itertools.product(*(range(k) for k in kmax)):
        ii = []
        for d in range(3):
            j = lo_i[d] + np.clip(np.minimum(off[d], w[d] - 1), 0, None)
            ii.append(np.mod(j, grid.n[d]) if d in grid.periodic else np.clip(j, 0, grid.n[d] - 1))
        flat = np.ravel_multi_index(ii, tuple(grid.n))
        rows.append(np.where(empty, grid.Vl.size, flat).ravel().astype(np.int32))
    return {"idx": np.stack(rows), "oob": oob.ravel(), "e": e, "kmax": kmax}


def reduce_successors(table, V, op):
    """op (np.minimum or np.maximum) of V over the successors of every cell; the sink
    (value -1) counts as a successor when the reach box leaves the state box."""
    Vext = np.append(V.ravel(), SINK)
    out = Vext[table["idx"][0]]
    for row in table["idx"][1:]:
        op(out, Vext[row], out=out)
    return np.where(table["oob"], op(out, SINK), out).reshape(V.shape)


def certificates(grid, dyn, actions, gamma, tau):
    Q, Mmax, kmax_all = [], None, []
    for u in actions:
        tab = successor_table(grid, dyn, u, tau)
        mlo = reduce_successors(tab, grid.Vl, np.minimum)
        Mhi = reduce_successors(tab, grid.Vu, np.maximum)
        e = tab["e"]
        Q.append(np.minimum(grid.ll, np.maximum(grid.rl, gamma * mlo)))
        Mmax = Mhi if Mmax is None else np.maximum(Mmax, Mhi)
        kmax_all.append(tab["kmax"])
    Q = np.stack(Q)
    pi = np.argmax(Q, axis=0).astype(np.int8)
    T_hat = Q.max(axis=0)
    T_bar = np.minimum(grid.lu, np.maximum(grid.ru, gamma * Mmax))
    green = grid.Vl > 0.0
    deficit = grid.Vl - T_hat     # > 0: lower certificate violated
    excess = T_bar - grid.Vu      # > 0: upper certificate violated
    stats = {
        "reach_box_radius": e,
        "max_successor_cells_per_axis": [max(k[d] for k in kmax_all) for d in range(3)],
        "lower_cells_checked": int(green.sum()),
        "lower_violations": int((deficit[green] > TOL).sum()),
        "lower_max_deficit": float(max(0.0, deficit[green].max())),
        "lower_violations_all_cells": int((deficit > TOL).sum()),
        "upper_cells_checked": int(grid.Vu.size),
        "upper_violations": int((excess > TOL).sum()),
        "upper_max_excess": float(max(0.0, excess.max())),
    }
    return pi, stats


# ------------------------------------------------------------------- simulation
def simulate(grid, dyn, actions, gamma, tau, pi, x0, max_steps, record=None):
    """Closed loop under u_t = actions[pi(cell(x_t))] with the exact flow.

    A trajectory is followed until t*, the first time with
    V_lower(s_t) <= r_lower(s_t) (the time used in the proof of Proposition 4),
    or until it enters F (l <= 0) or leaves the grid along x or y.
    """
    m = x0.shape[0]
    X, Y, TH = x0[:, 0].copy(), x0[:, 1].copy(), x0[:, 2].copy()
    cells0 = grid.locate([X, Y, TH])
    v0 = grid.Vl[cells0[0], cells0[1], cells0[2]]
    active = np.ones(m, dtype=bool)
    reach_t = np.full(m, -1)        # first time with r(x_t) > 0
    star_t = np.full(m, -1)         # t*
    fail_t = np.full(m, -1)         # time of entering F or leaving the grid
    run_min = np.full(m, np.inf)    # min over tau <= t of gamma^tau l(x_tau)
    J = np.full(m, -np.inf)         # max over t <= t* of g_gamma(x, t)
    margin = np.full(m, np.inf)     # min over t <= t* of gamma^t l(x_t) / V_lower(s0)
    paths = None if record is None else {int(i): [(X[i], Y[i], TH[i])] for i in record}
    for t in range(max_steps + 1):
        ia = np.flatnonzero(active)
        if ia.size == 0:
            break
        x, y = X[ia], Y[ia]
        out = (x < grid.lo[0]) | (x > grid.hi[0]) | (y < grid.lo[1]) | (y > grid.hi[1])
        lx, rx = l_fun(dyn, x, y), r_fun(dyn, x, y)
        g = gamma ** t
        run_min[ia] = np.minimum(run_min[ia], g * lx)
        J[ia] = np.maximum(J[ia], np.minimum(g * rx, run_min[ia]))
        failed = out | (lx <= 0.0)
        margin[ia] = np.where(failed, margin[ia], np.minimum(margin[ia], g * lx / v0[ia]))
        first = ~failed & (rx > 0.0) & (reach_t[ia] < 0)
        reach_t[ia[first]] = t
        fail_t[ia[failed]] = t
        active[ia[failed]] = False
        keep = ia[~failed]
        ci = grid.locate([X[keep], Y[keep], wrap(TH[keep])])
        at_star = grid.Vl[ci[0], ci[1], ci[2]] <= grid.rl[ci[0], ci[1], ci[2]]
        star_t[keep[at_star]] = t
        active[keep[at_star]] = False
        go = keep[~at_star]
        if go.size == 0 or t == max_steps:
            break
        a = pi[ci[0][~at_star], ci[1][~at_star], ci[2][~at_star]]
        for k, u in enumerate(actions):
            sel = go[a == k]
            if sel.size == 0:
                continue
            if paths is not None:
                for i in np.intersect1d(sel, record):
                    ts = np.linspace(0.0, tau, 12)[1:]
                    xs, ys, ths = flow(dyn["name"], dyn["v_const"], X[i], Y[i], TH[i], u, ts)
                    paths[int(i)].extend(zip(xs, ys, ths))
            X[sel], Y[sel], th1 = flow(dyn["name"], dyn["v_const"], X[sel], Y[sel], TH[sel], u, tau)
            TH[sel] = wrap(th1)
    return {"v0": v0, "reach_t": reach_t, "star_t": star_t, "fail_t": fail_t, "unfinished": active, "J": J,
            "margin": margin, "cells0": np.stack(cells0, axis=1)}, paths


# ------------------------------------------------------------------- avoid case
def upper_avoid_iteration(grid, tables):
    """Decreasing iteration W <- min(l_upper, max_a max_{s' in Delta(s,a)} W(s')) from
    W = l_upper on the verifier's abstraction, run until nothing changes.

    Returns the fixed point, the number of iterations, the first iteration k at
    which W^(k)(s) < 0 for every cell (-1 if never) and the per-action maxima at
    the fixed point. Every W^(k) is an upper bound on V_1 (proof of Proposition 6).
    """
    W = grid.lu.copy()
    first_neg = np.where(W < 0.0, 0, -1)
    k = 0
    while True:
        k += 1
        M = np.stack([reduce_successors(t, W, np.maximum) for t in tables])
        Wn = np.minimum(grid.lu, M.max(axis=0))
        first_neg[(Wn < 0.0) & (first_neg < 0)] = k
        if np.array_equal(Wn, W):
            return W, k, first_neg, M
        W = Wn


def advance(dyn, actions, tau, X, Y, TH, idx, a):
    """One step of the exact flow for the states idx under the action indices a."""
    for k, u in enumerate(actions):
        sel = idx[a == k]
        if sel.size:
            X[sel], Y[sel], th1 = flow(dyn["name"], dyn["v_const"], X[sel], Y[sel], TH[sel], u, tau)
            TH[sel] = wrap(th1)


def outside_grid(grid, x, y):
    return (x < grid.lo[0]) | (x > grid.hi[0]) | (y < grid.lo[1]) | (y > grid.hi[1])


def simulate_avoid(grid, dyn, actions, tau, pi, x0, steps, record=None):
    """Closed loop under u_t = actions[pi(cell(x_t))] with the exact flow for `steps` steps.

    Records the first time a trajectory enters F (l <= 0) or leaves the grid, the
    smallest l(x_t), and the smallest V_lower of the cells it visits.
    """
    m = x0.shape[0]
    X, Y, TH = x0[:, 0].copy(), x0[:, 1].copy(), x0[:, 2].copy()
    cells0 = grid.locate([X, Y, TH])
    v0 = grid.Vl[cells0[0], cells0[1], cells0[2]]
    alive = np.ones(m, dtype=bool)
    fail_t = np.full(m, -1)          # time of entering F or leaving the grid
    l_min = np.full(m, np.inf)       # min over t of l(x_t)
    vl_min = np.full(m, np.inf)      # min over t of V_lower(s_t)
    paths = None if record is None else {int(i): [(X[i], Y[i], TH[i])] for i in record}
    for t in range(steps + 1):
        ia = np.flatnonzero(alive)
        if ia.size == 0:
            break
        x, y = X[ia], Y[ia]
        lx = l_fun(dyn, x, y)
        l_min[ia] = np.minimum(l_min[ia], lx)
        failed = outside_grid(grid, x, y) | (lx <= 0.0)
        fail_t[ia[failed]] = t
        alive[ia[failed]] = False
        keep = ia[~failed]
        ci = grid.locate([X[keep], Y[keep], TH[keep]])
        vl_min[keep] = np.minimum(vl_min[keep], grid.Vl[ci[0], ci[1], ci[2]])
        if keep.size == 0 or t == steps:
            break
        a = pi[ci[0], ci[1], ci[2]]
        if paths is not None and t < AVOID_PATH_STEPS:
            for i in np.intersect1d(keep, record):
                u = actions[int(a[np.searchsorted(keep, i)])]
                ts = np.linspace(0.0, tau, 12)[1:]
                xs, ys, ths = flow(dyn["name"], dyn["v_const"], X[i], Y[i], TH[i], u, ts)
                paths[int(i)].extend(zip(xs, ys, ths))
        advance(dyn, actions, tau, X, Y, TH, keep, a)
    return {"v0": v0, "fail_t": fail_t, "l_min": l_min, "vl_min": vl_min,
            "cells0": np.stack(cells0, axis=1)}, paths


def first_hits(grid, dyn, actions, tau, x0, steps, choose):
    """First time each trajectory has l(x_t) < 0 or leaves the grid (-1: neither within
    `steps` steps), with the action indices choose(cells, idx); also whether that
    first hit was leaving the grid."""
    m = x0.shape[0]
    X, Y, TH = x0[:, 0].copy(), x0[:, 1].copy(), x0[:, 2].copy()
    hit_t = np.full(m, -1)
    left = np.zeros(m, dtype=bool)
    alive = np.ones(m, dtype=bool)
    for t in range(steps + 1):
        ia = np.flatnonzero(alive)
        if ia.size == 0:
            break
        x, y = X[ia], Y[ia]
        out = outside_grid(grid, x, y)
        hit = out | (l_fun(dyn, x, y) < 0.0)
        hit_t[ia[hit]] = t
        left[ia[out]] = True
        alive[ia[hit]] = False
        keep = ia[~hit]
        if keep.size == 0 or t == steps:
            break
        ci = grid.locate([X[keep], Y[keep], TH[keep]])
        advance(dyn, actions, tau, X, Y, TH, keep, choose(ci, keep))
    return hit_t, left


def verify_avoid(run_dir, out_dir, rng, points_per_cell, n_paths, horizon):
    """Checks of Proposition 6 (avoid value function, gamma = 1).

    1. the cell bounds l_lower, l_upper, recomputed from the centers;
    2. lower-bound certificate: V_lower(s) <= l_lower(s) and
       min_{s' in Delta(s, pi(s))} V_lower(s') >= V_lower(s) for every cell with
       V_lower > 0, where pi(s) maximizes min_{s'} V_lower(s') (the policy of
       Proposition 6). This is the invariance used in the proof, so every state of
       such a cell keeps l(x_t) >= V_lower(s0) for all t under pi;
    3. upper bound: V_upper >= W in every cell, where W is the fixed point of
       upper_avoid_iteration, so the cells with V_upper < 0 lie in the BRS;
    4. closed-loop simulation with the exact flow from random states in every cell
       with V_lower > 0 under pi for `horizon` steps: no entry into F, no exit from
       the grid, and l(x_t) >= V_lower(s0) at every step;
    5. BRS cells: from a random state in every cell with V_upper < 0, every tested
       control (each constant action, the greedy action on W, uniformly random
       actions) reaches l < 0 or leaves the grid within kbar(s0) steps, the first
       iteration at which W^(k)(s0) < 0.
    Leaving the grid counts as failure, as in the solver (sink with value -1).
    """
    t_start = time.time()
    params, data, phases = load_run(run_dir)
    dyn, args = params["dynamics"], params["args"]
    name = dyn["name"]
    tau, dt = float(dyn["tau"]), float(dyn["dt"])
    if args["mode"] != "avoid_nodiscount":
        raise ValueError("Proposition 6 needs an avoid run")
    if abs(tau - dt) > 1e-12:
        raise ValueError("expects one checkpoint per step (dt == tau)")
    actions = [float(a) for a in dyn["actions"]]
    grid = Grid(dyn, data)
    report = {"system": name, "run_dir": str(run_dir), "cells": int(grid.Vl.size), "grid": grid.n.tolist(),
              "gamma": 1.0, "tau": tau, "actions": actions, "L_f": dyn["L_f"], "L_l": dyn["L_l"]}
    p0 = phases[0]
    report["solver"] = {"vi_iterations": int(p0["vi_iterations"]), "stop_tolerance": float(args["phase0_tol"]),
                        "time_vi_s": float(p0["time_vi_s"]), "time_successors_s": float(p0["time_successors_s"]),
                        "time_phase_s": float(p0["time_phase_s"])}
    green, red = grid.Vl > 0, grid.Vu < 0
    report["classification_percent"] = {"safe": 100.0 * green.mean(), "BRS": 100.0 * red.mean(),
                                        "unclassified": 100.0 * (~green & ~red).mean()}
    report["cells_with_V_upper_zero"] = int((grid.Vu == 0).sum())

    # 1. exact flow and cell bounds of l
    report["flow_check_max_error"] = check_flow(name, dyn["v_const"], actions, tau, rng)
    X, Y, _ = grid.C
    lc = l_fun(dyn, X, Y)
    report["cell_bound_max_error"] = float(max(np.abs(grid.ll - (lc - dyn["L_l"] * grid.r_half)).max(),
                                               np.abs(grid.lu - (lc + dyn["L_l"] * grid.r_half)).max()))

    # 2. lower certificate on the independently computed abstraction
    tables = [successor_table(grid, dyn, u, tau) for u in actions]
    m = np.stack([reduce_successors(t, grid.Vl, np.minimum) for t in tables])
    pi = np.argmax(m, axis=0).astype(np.int8)
    deficit = np.maximum(grid.Vl - grid.ll, grid.Vl - m.max(axis=0))    # > 0: certificate violated
    del m
    cert = {
        "reach_box_radius": tables[0]["e"],
        "max_successor_cells_per_axis": [max(t["kmax"][d] for t in tables) for d in range(3)],
        "lower_cells_checked": int(green.sum()),
        "lower_violations": int((deficit[green] > TOL).sum()),
        "lower_max_deficit": float(max(0.0, deficit[green].max())) if green.any() else 0.0,
        "lower_violations_all_cells": int((deficit > TOL).sum()),
    }

    # 3. upper bound: recomputed decreasing iteration
    W, iters, first_neg, M = upper_avoid_iteration(grid, tables)
    pibar = np.argmax(M, axis=0).astype(np.int8)
    del M, tables
    excess = W - grid.Vu    # > 0: the solver's V_upper is below an upper bound on V_1
    cert.update({
        "upper_iterations": iters,
        "upper_violations": int((excess > TOL).sum()),
        "upper_max_excess": float(max(0.0, excess.max())),
        "upper_cells_equal": int((W == grid.Vu).sum()),
        "BRS_cells_recomputed": int((W < 0).sum()),
    })
    report["certificates"] = cert

    # 4. closed-loop simulation from random states in every cell with V_lower > 0
    cells = np.repeat(np.argwhere(green), points_per_cell, axis=0)
    x0 = grid.lo + (cells + rng.uniform(size=cells.shape)) * grid.h
    sim = {"initial_states": int(cells.shape[0]), "points_per_cell": points_per_cell, "horizon_steps": horizon}
    res, paths = None, None
    if cells.size:
        k = int(np.clip(np.floor((PAPER_SLICE[name] - grid.lo[2]) / grid.h[2]), 0, grid.n[2] - 1))
        in_slice = np.flatnonzero(cells[:, 2] == k)
        record = rng.choice(in_slice, size=min(n_paths, in_slice.size), replace=False) if in_slice.size else np.array([], int)
        x0[record] = grid.lo + (cells[record] + 0.5) * grid.h
        res, paths = simulate_avoid(grid, dyn, actions, tau, pi, x0, horizon, record=record)
        v0, ok = res["v0"], res["fail_t"] < 0
        sim.update({
            "stayed_safe": int(ok.sum()),
            "entered_F_or_left_grid": int((~ok).sum()),
            "distance_violations": int((res["l_min"] < v0 - TOL).sum()),
            "margin_min": float((res["l_min"] / v0).min()),
            "l_min_minus_Vlower_min": float((res["l_min"] - v0).min()),
            "visited_cells_below_Vlower_s0": int((res["vl_min"] < v0 - TOL).sum()),
            "Vlower_s0_median": float(np.median(v0)),
            "Vlower_s0_max": float(v0.max()),
            "l_min_median": float(np.median(res["l_min"])),
        })
        res["x0"] = x0
    report["simulation"] = sim

    # 5. random states in the cells with V_upper < 0 under several controls
    rcells = np.argwhere(red)
    brs = {"samples": int(rcells.shape[0])}
    hits = None
    if rcells.size:
        y0 = grid.lo + (rcells + rng.uniform(size=rcells.shape)) * grid.h
        kbar = first_neg[rcells[:, 0], rcells[:, 1], rcells[:, 2]]
        strategies = [(f"constant u={u:g}", lambda ci, idx, k=k: np.full(idx.size, k, dtype=np.int8))
                      for k, u in enumerate(actions)]
        strategies.append(("greedy on W", lambda ci, idx: pibar[ci[0], ci[1], ci[2]]))
        strategies.append(("uniformly random", lambda ci, idx: rng.integers(len(actions), size=idx.size)))
        latest = np.zeros(len(y0), dtype=int)
        per = {}
        for label, choose in strategies:
            hit_t, left = first_hits(grid, dyn, actions, tau, y0, int(kbar.max()), choose)
            bad = (kbar < 0) | (hit_t < 0) | (hit_t > kbar)
            per[label] = {"violations": int(bad.sum()), "left_grid": int(left.sum()),
                          "entered_F": int(((hit_t >= 0) & ~left).sum()),
                          "hit_time_median": float(np.median(hit_t[hit_t >= 0])), "hit_time_max": int(hit_t.max())}
            latest = np.maximum(latest, np.where(hit_t < 0, np.iinfo(int).max, hit_t))
        brs.update({
            "start_outside_F": int((l_fun(dyn, y0[:, 0], y0[:, 1]) >= 0).sum()),
            "strategies": per,
            "violations": int(sum(p["violations"] for p in per.values())),
            "trajectories": int(len(y0) * len(strategies)),
            "bound_median": float(np.median(kbar)),
            "bound_max": int(kbar.max()),
            "ratio_latest_hit_to_bound_max": float((latest / np.maximum(kbar, 1)).max()),
        })
        hits = {"y0": y0, "kbar": kbar, "latest": latest}
    report["brs"] = brs
    report["verifier_seconds"] = time.time() - t_start

    sys_dir = out_dir / name
    sys_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(sys_dir / "trajectories.npz", **(res or {}), **{f"brs_{k}": v for k, v in (hits or {}).items()})
    (sys_dir / "report.json").write_text(json.dumps(report, indent=2))
    return report, grid, dyn, res, hits, paths


# ---------------------------------------------------------------------- figures
CLASS_CMAP = ListedColormap(["#c62828", "#9e9e9e", "#2e7d32"])   # red, gray, green


def classification_panel(ax, grid, dyn, theta, paths, title, avoid=False):
    """Green: V_lower > 0. Red: V_upper <= 0 (reach-avoid) or V_upper < 0 (avoid)."""
    k = int(np.clip(np.floor((theta - grid.lo[2]) / grid.h[2]), 0, grid.n[2] - 1))
    unsafe = grid.Vu[:, :, k] < 0 if avoid else grid.Vu[:, :, k] <= 0
    cls = np.where(grid.Vl[:, :, k] > 0, 2, np.where(unsafe, 0, 1))
    ax.imshow(cls.T, origin="lower", cmap=CLASS_CMAP, vmin=0, vmax=2, interpolation="nearest",
              extent=(grid.lo[0], grid.hi[0], grid.lo[1], grid.hi[1]))
    ax.add_patch(Circle(dyn["obstacle_position"], dyn["obstacle_radius"], fill=False, color="#0d47a1", lw=1.5))
    if not avoid:
        ax.add_patch(Circle(dyn["target_position"], dyn["target_radius"], fill=False, color="#f9a825", lw=1.5,
                            ls="--"))
    for pts in (paths or {}).values():
        p = np.array(pts)
        ax.plot(p[:, 0], p[:, 1], color="black", lw=0.8)
        ax.plot(p[0, 0], p[0, 1], "o", color="white", mec="black", ms=3.5, mew=0.8)
    ax.set_xlim(grid.lo[0], grid.hi[0])
    ax.set_ylim(grid.lo[1], grid.hi[1])
    ax.set_aspect("equal")
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.set_title(title, fontsize=9)


def bound_panel(ax, res, bound):
    ok = res["reach_t"] >= 0
    t, b = res["reach_t"][ok], bound[ok]
    hi = int(max(t.max(), b.max())) + 1
    H, xe, ye = np.histogram2d(t, b, bins=[np.arange(-0.5, hi + 0.5), np.arange(-0.5, hi + 0.5)])
    ax.pcolormesh(xe, ye, H.T, norm=LogNorm(vmin=1, vmax=max(1, H.max())), cmap="viridis", rasterized=True)
    ax.plot([0, hi], [0, hi], color="#c62828", lw=1)
    ax.set_xlim(-0.5, t.max() + 1)
    ax.set_ylim(-0.5, b.max() + 1)
    ax.set_xlabel("observed reach time (steps)")
    ax.set_ylabel("guaranteed bound (steps)")


def value_panel(ax, res):
    ok = res["reach_t"] >= 0
    v0, J = res["v0"][ok], res["J"][ok]
    hb = ax.hexbin(v0, J, gridsize=60, bins="log", cmap="viridis", mincnt=1, linewidths=0)
    hb.set_rasterized(True)
    top = max(v0.max(), J.max())
    ax.plot([0, top], [0, top], color="#c62828", lw=1)
    ax.set_xlabel(LABELS["v_lower"])
    ax.set_ylabel(LABELS["achieved"])


def margin_panel(ax, res, horizon):
    v0, lm = res["v0"], res["l_min"]
    hb = ax.hexbin(v0, lm, gridsize=60, bins="log", cmap="viridis", mincnt=1, linewidths=0)
    hb.set_rasterized(True)
    top = max(v0.max(), lm.max())
    ax.plot([0, top], [0, top], color="#c62828", lw=1)
    ax.set_xlabel(LABELS["v1_lower"])
    ax.set_ylabel(rf"$\min_{{t\leq {horizon}}} l(x_t)$")


def brs_panel(ax, hits):
    t, b = hits["latest"], hits["kbar"]
    hi = int(max(t.max(), b.max())) + 1
    H, xe, ye = np.histogram2d(t, b, bins=[np.arange(-0.5, hi + 0.5), np.arange(-0.5, hi + 0.5)])
    ax.pcolormesh(xe, ye, H.T, norm=LogNorm(vmin=1, vmax=max(1, H.max())), cmap="viridis", rasterized=True)
    ax.plot([0, hi], [0, hi], color="#c62828", lw=1)
    ax.set_xlim(-0.5, t.max() + 1)
    ax.set_ylim(-0.5, b.max() + 1)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("latest failure time (steps)")
    ax.set_ylabel("guaranteed bound (steps)")


def note_panel(ax, text):
    ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes, fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])


LABELS = {"v_lower": r"$V_{\gamma,\mathrm{lower}}^{\mathrm{cons}}(s_0)$",
          "v1_lower": r"$V_{1,\mathrm{lower}}^{\mathrm{cons}}(s_0)$",
          "v1_safe": r"$V_{1,\mathrm{lower}}^{\mathrm{cons}}>0$",
          "v1_brs": r"$V_{1,\mathrm{upper}}^{\mathrm{cons}}<0$",
          "achieved": r"achieved $\max_{t\leq t^*} g_\gamma(\xi,t)$"}


def use_latex():
    """Typeset figure text with LaTeX when it is installed (matches the paper)."""
    import os
    import shutil

    os.environ["PATH"] = "/Library/TeX/texbin:" + os.environ.get("PATH", "")
    if shutil.which("latex") and shutil.which("dvipng"):
        plt.rcParams.update({"text.usetex": True, "font.family": "serif",
                             "text.latex.preamble": r"\usepackage{amsmath}"})
        LABELS["v_lower"] = r"$\underline{V}_\gamma^{\mathrm{cons}}(s_0)$"
        LABELS["v1_lower"] = r"$\underline{V}_1^{\mathrm{cons}}(s_0)$"
        LABELS["v1_safe"] = r"$\underline{V}_1^{\mathrm{cons}}>0$"
        LABELS["v1_brs"] = r"$\overline{V}_1^{\mathrm{cons}}<0$"


# ------------------------------------------------------------------------- main
def verify(run_dir, out_dir, rng, points_per_cell, n_paths):
    t_start = time.time()
    params, data, phases = load_run(run_dir)
    dyn, args = params["dynamics"], params["args"]
    name = dyn["name"]
    gamma, tau, dt = float(args["gamma"]), float(dyn["tau"]), float(dyn["dt"])
    if args["mode"] != "ra_discount" or not (0.0 < gamma < 1.0):
        raise ValueError("Proposition 4 needs a discounted reach-avoid run")
    if abs(tau - dt) > 1e-12:
        raise ValueError("expects one checkpoint per step (dt == tau)")
    actions = [float(a) for a in dyn["actions"]]
    grid = Grid(dyn, data)
    report = {"system": name, "run_dir": str(run_dir), "cells": int(grid.Vl.size), "grid": grid.n.tolist(),
              "gamma": gamma, "tau": tau, "actions": actions, "L_f": dyn["L_f"], "L_l": dyn["L_l"], "L_r": dyn["L_r"]}
    p0 = phases[0]
    report["solver"] = {"vi_iterations": int(p0["vi_iterations"]), "time_vi_s": float(p0["time_vi_s"]),
                        "time_successors_s": float(p0["time_successors_s"]), "time_phase_s": float(p0["time_phase_s"])}
    green, red = grid.Vl > 0, grid.Vu <= 0
    report["classification_percent"] = {"R": 100.0 * green.mean(), "L": 100.0 * red.mean(),
                                        "unclassified": 100.0 * (~green & ~red).mean()}

    # 1. exact flow and cell bounds
    report["flow_check_max_error"] = check_flow(name, dyn["v_const"], actions, tau, rng)
    X, Y, _ = grid.C
    lc, rc = l_fun(dyn, X, Y), r_fun(dyn, X, Y)
    rh = grid.r_half
    report["cell_bound_max_error"] = float(max(
        np.abs(grid.ll - (lc - dyn["L_l"] * rh)).max(), np.abs(grid.lu - (lc + dyn["L_l"] * rh)).max(),
        np.abs(grid.rl - (rc - dyn["L_r"] * rh)).max(), np.abs(grid.ru - (rc + dyn["L_r"] * rh)).max()))

    # 2-3. certificates on the independently computed abstraction
    pi, cert = certificates(grid, dyn, actions, gamma, tau)
    report["certificates"] = cert

    # 4. closed-loop simulation from random states in every certified cell
    cells = np.repeat(np.argwhere(green), points_per_cell, axis=0)
    x0 = grid.lo + (cells + rng.uniform(size=cells.shape)) * grid.h
    v0 = grid.Vl[cells[:, 0], cells[:, 1], cells[:, 2]]
    r_star = float(grid.rl.max())
    bound = np.ceil(np.log(v0 / r_star) / math.log(gamma)).astype(int)
    # Example paths: start at cell centers in the theta slice shown in the paper.
    k = int(np.clip(np.floor((PAPER_SLICE[name] - grid.lo[2]) / grid.h[2]), 0, grid.n[2] - 1))
    in_slice = np.flatnonzero(cells[:, 2] == k)
    record = rng.choice(in_slice, size=min(n_paths, in_slice.size), replace=False) if in_slice.size else np.array([], int)
    x0[record] = grid.lo + (cells[record] + 0.5) * grid.h
    res, paths = simulate(grid, dyn, actions, gamma, tau, pi, x0, int(bound.max()) + 5, record=record)
    res["bound"] = bound
    # ok: no failure before t*, t* reached, and the target reached no later than t*.
    ok = (res["fail_t"] < 0) & (res["star_t"] >= 0) & (res["reach_t"] >= 0) & (res["reach_t"] <= res["star_t"])
    vu0 = grid.Vu[cells[:, 0], cells[:, 1], cells[:, 2]]
    sim = {
        "initial_states": int(cells.shape[0]),
        "points_per_cell": points_per_cell,
        "r_star": r_star,
        "reached_safely": int(ok.sum()),
        "entered_F_or_left_grid": int((res["fail_t"] >= 0).sum()),
        "unfinished": int(res["unfinished"].sum()),
        "time_bound_violations": int((~ok | (res["star_t"] > bound)).sum()),
        "t_star_median": float(np.median(res["star_t"][ok])),
        "t_star_max": int(res["star_t"][ok].max()),
        "reach_time_median": float(np.median(res["reach_t"][ok])),
        "reach_time_max": int(res["reach_t"][ok].max()),
        "bound_median": float(np.median(bound[ok])),
        "bound_max": int(bound.max()),
        "slack_median": float(np.median((bound - res["reach_t"])[ok])),
        "slack_min": int((bound - res["reach_t"])[ok].min()),
        "ratio_time_to_bound_max": float((res["reach_t"][ok] / np.maximum(bound[ok], 1)).max()),
        "margin_min": float(res["margin"][ok].min()),
        "margin_violations": int((res["margin"][ok] < 1.0 - 1e-9).sum()),
        "J_minus_Vlower_min": float((res["J"] - v0)[ok].min()),
        "J_below_Vlower": int((res["J"][ok] < v0[ok] - TOL).sum()),
        "J_minus_Vupper_max": float((res["J"] - vu0)[ok].max()),
        "J_above_Vupper": int((res["J"][ok] > vu0[ok] + TOL).sum()),
    }
    report["simulation"] = sim
    report["verifier_seconds"] = time.time() - t_start

    sys_dir = out_dir / name
    sys_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(sys_dir / "trajectories.npz", x0=x0, vu0=vu0, **res)
    (sys_dir / "report.json").write_text(json.dumps(report, indent=2))
    return report, grid, dyn, res, bound, paths


def write_outputs(results, out_dir):
    use_latex()
    names = {"dubins": "Dubins car", "evasion": "3D evasion"}
    slice_label = {"dubins": r"$\theta=0$", "evasion": r"$\theta=\pi$"}
    # Per-system panels and a combined figure (rows: systems; columns: panels).
    fig, axes = plt.subplots(len(results), 3, figsize=(7.0, 2.35 * len(results)), squeeze=False,
                             gridspec_kw={"width_ratios": [1.0, 1.05, 1.05]})
    for row, (rep, grid, dyn, res, bound, paths) in enumerate(results):
        n = rep["system"]
        classification_panel(axes[row, 0], grid, dyn, PAPER_SLICE[n], paths, f"{names[n]}, {slice_label[n]}")
        bound_panel(axes[row, 1], res, bound)
        value_panel(axes[row, 2], res)
        for j, (f, fn) in enumerate([(classification_panel, "classification"), (bound_panel, "bound"),
                                     (value_panel, "value")]):
            f1, a1 = plt.subplots(figsize=(3.3, 2.6))
            if fn == "classification":
                f(a1, grid, dyn, PAPER_SLICE[n], paths, f"{names[n]}, {slice_label[n]}")
            elif fn == "bound":
                f(a1, res, bound)
            else:
                f(a1, res)
            f1.tight_layout()
            for ext in ("pdf", "png"):
                f1.savefig(out_dir / n / f"{fn}.{ext}", dpi=300)
            plt.close(f1)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out_dir / f"ra_verification.{ext}", dpi=300)
    plt.close(fig)

    # LaTeX table (booktabs), one column per system.
    reps = [r[0] for r in results]
    cols = " & ".join(names[r["system"]] for r in reps)

    def row(label, fmt):
        return f"{label} & " + " & ".join(fmt(r) for r in reps) + r" \\"

    lines = [
        r"\begin{table}[t]",
        r"\caption{Verification of Proposition~\ref{prop:greedy} for Algorithm~\ref{alg:alg1} with "
        r"$\gamma=0.96$ on a uniform grid of $120^3$ cells.}",
        r"\label{tab:ra_verification}",
        r"\small",
        r"\begin{tabular}{l" + "r" * len(reps) + "}",
        r"\toprule",
        f" & {cols} \\\\",
        r"\midrule",
        row(r"Cells certified in $\mathcal R$ / $\mathcal L$ (\%)",
            lambda r: f"{r['classification_percent']['R']:.1f} / {r['classification_percent']['L']:.1f}"),
        row(r"Certificate violations (lower / upper)",
            lambda r: f"{r['certificates']['lower_violations']} / {r['certificates']['upper_violations']}"),
        row(r"Simulated initial states", lambda r: f"{r['simulation']['initial_states']:,}".replace(",", "{,}")),
        row(r"Reached $\mathcal T$ without entering $\mathcal F$",
            lambda r: f"{100.0 * r['simulation']['reached_safely'] / r['simulation']['initial_states']:.1f}\\%"),
        row(r"Reach-time bound violations", lambda r: f"{r['simulation']['time_bound_violations']}"),
        row(r"Median observed / bound (steps)",
            lambda r: f"{r['simulation']['reach_time_median']:.0f} / {r['simulation']['bound_median']:.0f}"),
        row(r"$\min_t \gamma^t l(x_t)/\underline V_\gamma^{\mathrm{cons}}(s_0)$",
            lambda r: f"{r['simulation']['margin_min']:.2f}"),
        row(r"$\min\,(J-\underline V_\gamma^{\mathrm{cons}}(s_0))$",
            lambda r: f"{r['simulation']['J_minus_Vlower_min']:.3f}"),
        row(r"Abstraction and value iteration time (s)",
            lambda r: f"{r['solver']['time_successors_s'] + r['solver']['time_vi_s']:.1f}"),
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ]
    (out_dir / "table_ra.tex").write_text("\n".join(lines) + "\n")
    (out_dir / "summary.json").write_text(json.dumps(reps, indent=2))


def write_outputs_avoid(results, out_dir, horizon):
    use_latex()
    names = {"dubins": "Dubins car", "evasion": "3D evasion"}
    slice_label = {"dubins": r"$\theta=0$", "evasion": r"$\theta=\pi$"}

    def draw(panel, ax, rep, grid, dyn, res, hits, paths):
        n = rep["system"]
        if panel == "classification":
            classification_panel(ax, grid, dyn, PAPER_SLICE[n], paths, f"{names[n]}, {slice_label[n]}", avoid=True)
        elif panel == "margin" and res is not None:
            margin_panel(ax, res, horizon)
        elif panel == "margin":
            note_panel(ax, f"no cell with {LABELS['v1_safe']}")
        elif hits is not None:
            brs_panel(ax, hits)
        else:
            note_panel(ax, f"no cell with {LABELS['v1_brs']}")

    panels = ("classification", "margin", "brs")
    fig, axes = plt.subplots(len(results), 3, figsize=(7.0, 2.35 * len(results)), squeeze=False,
                             gridspec_kw={"width_ratios": [1.0, 1.05, 1.05]})
    for row, r in enumerate(results):
        for j, panel in enumerate(panels):
            draw(panel, axes[row, j], *r)
            f1, a1 = plt.subplots(figsize=(3.3, 2.6))
            draw(panel, a1, *r)
            f1.tight_layout()
            for ext in ("pdf", "png"):
                f1.savefig(out_dir / r[0]["system"] / f"{panel}.{ext}", dpi=300)
            plt.close(f1)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out_dir / f"avoid_verification.{ext}", dpi=300)
    plt.close(fig)

    # LaTeX table (booktabs), one column per system.
    reps = [r[0] for r in results]
    cols = " & ".join(names[r["system"]] for r in reps)

    def row(label, fmt):
        return f"{label} & " + " & ".join(fmt(r) for r in reps) + r" \\"

    def sim(fmt):
        return lambda r: fmt(r["simulation"]) if r["simulation"]["initial_states"] else "--"

    def brs(fmt):
        return lambda r: fmt(r["brs"]) if r["brs"]["samples"] else "--"

    lines = [
        r"\begin{table}[t]",
        r"\caption{Verification of Proposition~\ref{prop:avoid} for Algorithm~\ref{alg:alg1} with the avoid "
        r"operators ($\gamma=1$) on a uniform grid of $120^3$ cells.}",
        r"\label{tab:avoid_verification}",
        r"\small",
        r"\begin{tabular}{l" + "r" * len(reps) + "}",
        r"\toprule",
        f" & {cols} \\\\",
        r"\midrule",
        row(r"Cells with $\underline V_1^{\mathrm{cons}}>0$ / $\overline V_1^{\mathrm{cons}}<0$ (\%)",
            lambda r: f"{r['classification_percent']['safe']:.1f} / {r['classification_percent']['BRS']:.1f}"),
        row(r"Certificate violations (lower / upper)",
            lambda r: f"{r['certificates']['lower_violations']} / {r['certificates']['upper_violations']}"),
        row(r"Simulated initial states ($\underline V_1^{\mathrm{cons}}>0$)",
            sim(lambda s: f"{s['initial_states']:,}".replace(",", "{,}"))),
        row(rf"Outside $\mathcal F$ and in the grid for {horizon} steps",
            sim(lambda s: f"{100.0 * s['stayed_safe'] / s['initial_states']:.1f}\\%")),
        row(r"$\min_t l(x_t)/\underline V_1^{\mathrm{cons}}(s_0)$", sim(lambda s: f"{s['margin_min']:.3f}")),
        row(r"Sampled states ($\overline V_1^{\mathrm{cons}}<0$) $\times$ controls",
            brs(lambda b: f"{b['trajectories']:,}".replace(",", "{,}"))),
        row(r"Entered $\mathcal F$ or left the grid within the bound",
            brs(lambda b: f"{100.0 * (b['trajectories'] - b['violations']) / b['trajectories']:.1f}\\%")),
        row(r"Abstraction and value iteration time (s)",
            lambda r: f"{r['solver']['time_successors_s'] + r['solver']['time_vi_s']:.1f}"),
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ]
    (out_dir / "table_avoid.tex").write_text("\n".join(lines) + "\n")
    (out_dir / "summary.json").write_text(json.dumps(reps, indent=2))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dirs", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--points-per-cell", type=int, default=1)
    ap.add_argument("--paths", type=int, default=14, help="example trajectories drawn per system")
    ap.add_argument("--horizon", type=int, default=1000, help="simulated steps in the avoid case")
    a = ap.parse_args()
    modes = {json.loads((Path(rd) / "params.json").read_text())["args"]["mode"] for rd in a.run_dirs}
    if len(modes) != 1:
        raise ValueError(f"all runs must use the same mode, got {sorted(modes)}")
    avoid = modes.pop() == "avoid_nodiscount"
    out_dir = Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    results = []
    for rd in a.run_dirs:
        print(f"verifying {rd}", flush=True)
        if avoid:
            results.append(verify_avoid(Path(rd), out_dir, rng, a.points_per_cell, a.paths, a.horizon))
        else:
            results.append(verify(Path(rd), out_dir, rng, a.points_per_cell, a.paths))
        print(json.dumps(results[-1][0], indent=2), flush=True)
    if avoid:
        write_outputs_avoid(results, out_dir, a.horizon)
    else:
        write_outputs(results, out_dir)
    print(f"wrote {out_dir}", flush=True)


if __name__ == "__main__":
    main()
