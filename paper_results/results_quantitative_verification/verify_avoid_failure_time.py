#!/usr/bin/env python3
"""Empirical check of the failure-reaching time of Proposition 6(2) (avoid case, gamma = 1).

Claim (paper, 2026-11-04): if A = U and V_upper^(k)(s0) <= 0 for some k >= 0, then, under
every control sequence, every trajectory starting in [[s0]] reaches F = {l <= 0} at some time
tau <= k. The smallest such k is k(s0), the first iteration of Algorithm 1 with
V_upper^(k)(s0) <= 0. The solver writes it, for every cell, to
upper_nonpos_iteration_phase_0.csv (avoid runs with --dump-csv; -1 means never).

Data: the avoid run of Section 9.2.4 (Dubins car, 120^3 cells, gamma = 1, U = A = {-1, 0, 1}).
As in the solver and the paper's experiments, l = -1 outside the state box, so leaving the
box along x or y reaches F.

Checks:
  1. k(s) from the solver against the verifier's own recomputation: the decreasing iteration
     W^(k+1) = min(l_upper, max_a max_{s' in Delta(s,a)} W^(k)) from W^(0) = l_upper, on the
     abstraction that verify_quantitative.successor_table builds independently of the solver
     (sink value -1), with the first k at which W^(k)(s) <= 0;
  2. a certificate of the bound on that abstraction, i.e., the induction in the proof:
     k(s) = 0 only if l_upper(s) <= 0 (then the whole cell lies in F), and for k(s) >= 1,
     l_upper(s) > 0 and every successor s' of every action has k(s') <= k(s) - 1 (the sink
     counts as k = 0). Since the abstraction is sound, this proves the bound for every state
     of the cell, not only for the sampled ones;
  3. exhaustive search from one uniformly random state and the 8 vertices of every cell with
     V_upper <= 0: every action sequence is followed with the exact flow until it first reaches
     F (l(x_t) <= 0 or leaving the box). A branch that has not reached F at time k(s0) is a
     violation. The search covers every control sequence, hence every policy, stationary or
     not.

Usage (the run reported in the paper, about a minute):
    python verify_avoid_failure_time.py RUN_DIR --out verification_avoid_failure_time --seed 0
"""
import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_quantitative import (SINK, Grid, check_flow, flow, l_fun, load_run,  # noqa: E402
                                 outside_grid, reduce_successors, successor_table, wrap)

NEVER = 1 << 20   # stands for k = -1 (V_upper never <= 0) in the certificate


def load_k(run_dir, grid, data):
    """k(s) from the solver's dump, as an array over the grid (-1: never)."""
    kk = np.loadtxt(Path(run_dir) / "upper_nonpos_iteration_phase_0.csv", delimiter=",", skiprows=1,
                    dtype=np.int64, ndmin=2)
    # Join on cell_id with the value dump, whose rows give the cell boxes.
    order_v, order_k = np.argsort(data[:, 0], kind="stable"), np.argsort(kk[:, 0], kind="stable")
    if kk.shape[0] != data.shape[0] or not np.array_equal(data[order_v, 0].astype(np.int64), kk[order_k, 0]):
        raise ValueError("the k dump and the value dump list different cells")
    D = 3
    clo = data[:, 1:1 + 2 * D:2]
    idx = np.rint((clo - grid.lo) / grid.h).astype(int)
    flat = np.ravel_multi_index(idx.T, tuple(grid.n))
    K = np.full(int(np.prod(grid.n)), -1, dtype=np.int64)
    K[flat[order_v]] = kk[order_k, 1]
    return K.reshape(tuple(grid.n))


def recompute_k(grid, tables):
    """First k with W^(k) <= 0 on the verifier's abstraction, and the fixed point W."""
    W = grid.lu.copy()
    K = np.where(W <= 0.0, 0, -1)
    k = 0
    while True:
        k += 1
        Wn = np.minimum(grid.lu, np.stack([reduce_successors(t, W, np.maximum, SINK) for t in tables]).max(axis=0))
        K[(Wn <= 0.0) & (K < 0)] = k
        if np.array_equal(Wn, W):
            return K, W, k
        W = Wn


def certificate(grid, tables, K):
    """Induction of the proof on the verifier's abstraction (see the module docstring)."""
    Kf = np.where(K < 0, NEVER, K).astype(float)
    # Largest k over the successors of every cell and action; the sink has k = 0.
    Kmax = np.stack([reduce_successors(t, Kf, np.maximum, 0.0) for t in tables]).max(axis=0)
    zero, pos = K == 0, K >= 1
    return {
        "cells_k0": int(zero.sum()),
        "k0_with_l_upper_positive": int((zero & (grid.lu > 0.0)).sum()),
        "cells_k_positive": int(pos.sum()),
        "k_positive_with_l_upper_nonpositive": int((pos & (grid.lu <= 0.0)).sum()),
        "k_positive_with_a_successor_not_below": int((pos & (Kmax > K - 1)).sum()),
        "violations": int((zero & (grid.lu > 0.0)).sum() + (pos & ((grid.lu <= 0.0) | (Kmax > K - 1))).sum()),
    }


def exhaustive(dyn, grid, actions, tau, x0, kbar, cap):
    """All action sequences from each x0, each followed until it first reaches F."""
    m = x0.shape[0]
    out = {"latest_hit": np.full(m, -1), "leaves": np.zeros(m, np.int64), "max_frontier": np.zeros(m, np.int64),
           "no_hit_by_k": np.zeros(m, bool), "exhausted": np.ones(m, bool), "left_grid_first": np.zeros(m, bool)}

    def run(rows):
        X, Y, TH = x0[rows, 0].copy(), x0[rows, 1].copy(), x0[rows, 2].copy()
        own = np.arange(rows.size)
        acc = {k: v[rows].copy() for k, v in out.items()}
        kb = kbar[rows]
        t = 0
        while own.size:
            acc["max_frontier"] = np.maximum(acc["max_frontier"], np.bincount(own, minlength=rows.size))
            left = outside_grid(grid, X, Y)
            hit = left | (l_fun(dyn, X, Y) <= 0.0)
            if hit.any():
                oh = own[hit]
                np.maximum.at(acc["latest_hit"], oh, t)
                np.add.at(acc["leaves"], oh, 1)
                acc["left_grid_first"][own[hit & left]] = True
            expired = ~hit & (t >= kb[own])
            acc["no_hit_by_k"][own[expired]] = True
            go = ~hit & ~expired
            n = int(go.sum())
            if n == 0:
                break
            if n * len(actions) > cap:
                return None
            Xg, Yg, THg, og = X[go], Y[go], TH[go], own[go]
            parts = [flow(dyn["name"], dyn["v_const"], Xg, Yg, THg, u, tau) for u in actions]
            X = np.concatenate([p[0] for p in parts])
            Y = np.concatenate([p[1] for p in parts])
            TH = wrap(np.concatenate([p[2] for p in parts]))
            own = np.tile(og, len(actions))
            t += 1
        return acc

    def solve(rows):
        acc = run(rows)
        if acc is None:
            if rows.size == 1:
                out["exhausted"][rows] = False
                return
            half = rows.size // 2
            solve(rows[:half])
            solve(rows[half:])
            return
        for k in out:
            out[k][rows] = acc[k]

    batch = 4096
    for b0 in range(0, m, batch):
        solve(np.arange(b0, min(m, b0 + batch)))
        if (b0 // batch) % 100 == 0:
            print(f"  exhaustive: {min(m, b0 + batch)}/{m} states", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cap", type=int, default=12_000_000, help="largest frontier per search batch")
    a = ap.parse_args()
    t0 = time.time()
    out_dir = Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    params, data, phases = load_run(Path(a.run_dir))
    dyn, args = params["dynamics"], params["args"]
    if args["mode"] != "avoid_nodiscount":
        raise ValueError("needs an avoid run")
    tau = float(dyn["tau"])
    actions = [float(u) for u in dyn["actions"]]
    grid = Grid(dyn, data)
    K = load_k(a.run_dir, grid, data)
    red = grid.Vu <= 0.0
    report = {"run_dir": str(a.run_dir), "tau": tau, "actions": actions, "seed": a.seed,
              "solver_vi_iterations": int(phases[0]["vi_iterations"]), "cells": int(grid.Vu.size),
              "cells_V_upper_nonpositive": int(red.sum()), "cells_V_upper_zero": int((grid.Vu == 0.0).sum()),
              "cells_with_k": int((K >= 0).sum()),
              "cells_with_k_but_V_upper_positive_or_vice_versa": int(((K >= 0) != red).sum())}
    vals, counts = np.unique(K[K >= 0], return_counts=True)
    report["k_histogram_over_cells"] = {int(v): int(c) for v, c in zip(vals, counts)}
    report["flow_check_max_error"] = check_flow(dyn["name"], dyn["v_const"], actions, tau, rng)

    # 1 and 2: recomputation and certificate on the verifier's abstraction
    print("abstraction ...", flush=True)
    tables = [successor_table(grid, dyn, u, tau) for u in actions]
    K_ver, W, iters = recompute_k(grid, tables)
    report["recomputation"] = {"iterations": iters, "fixed_point_equals_solver_V_upper": int((W == grid.Vu).sum()),
                               "k_equal_cells": int((K_ver == K).sum()), "k_different_cells": int((K_ver != K).sum())}
    report["certificate"] = certificate(grid, tables, K)
    del tables
    print(json.dumps({k: report[k] for k in ("cells_with_k", "k_histogram_over_cells", "recomputation",
                                             "certificate")}), flush=True)

    # 3: exhaustive search
    cells = np.argwhere(red)
    x0 = grid.lo + (cells + rng.uniform(size=cells.shape)) * grid.h
    kind = np.zeros(cells.shape[0], dtype=np.int8)          # 0: uniform sample, 1: cell vertex
    offs = np.array(list(itertools.product((0, 1), repeat=3)))
    verts = (grid.lo + (cells[:, None, :] + offs[None, :, :]) * grid.h).reshape(-1, 3)
    x0 = np.concatenate([x0, verts])
    cells = np.concatenate([cells, np.repeat(cells, len(offs), axis=0)])
    kind = np.concatenate([kind, np.ones(verts.shape[0], dtype=np.int8)])
    kbar = K[cells[:, 0], cells[:, 1], cells[:, 2]]
    ts = time.time()
    ex = exhaustive(dyn, grid, actions, tau, x0, kbar, a.cap)
    done = ex["exhausted"]
    lh = ex["latest_hit"]
    bad = done & (ex["no_hit_by_k"] | (lh < 0) | (lh > kbar))
    pos = done & (kbar >= 1)
    report["exhaustive"] = {
        "states": int(x0.shape[0]), "uniform": int((kind == 0).sum()), "vertices": int((kind == 1).sum()),
        "exhausted": int(done.sum()), "not_exhausted_(frontier_cap)": int((~done).sum()),
        "violations": int(bad.sum()), "no_hit_by_k": int((done & ex["no_hit_by_k"]).sum()),
        "action_prefixes_checked": int(ex["leaves"][done].sum()), "max_frontier": int(ex["max_frontier"].max()),
        "states_in_F_at_time_0": int((done & (lh == 0)).sum()),
        "k_max": int(kbar.max()), "k_median": float(np.median(kbar)),
        "latest_hit_max": int(lh[done].max()),
        "latest_hit_equals_k_(k>=1)": int((pos & (lh == kbar)).sum()),
        "max_latest_hit_over_k_(k>=1)": float((lh[pos] / kbar[pos]).max()) if pos.any() else None,
        "min_k_minus_latest_hit": int((kbar[done] - lh[done]).min()),
        "seconds": round(time.time() - ts, 1),
    }
    print("exhaustive", json.dumps(report["exhaustive"]), flush=True)
    np.savez_compressed(out_dir / "exhaustive.npz", x0=x0, cells=cells, kind=kind, kbar=kbar, **ex)
    report["seconds"] = round(time.time() - t0, 1)
    (out_dir / "report.json").write_text(json.dumps(report, indent=2))
    print(f"wrote {out_dir}", flush=True)


if __name__ == "__main__":
    main()
