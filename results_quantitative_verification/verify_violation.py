#!/usr/bin/env python3
"""Empirical check of the violation claim of Proposition 4 (prop:greedy, second claim).

Claim (paper, 2026-11-03): if U is finite, A = U and V_upper(s0) < 0, then under every
policy the trajectory from any x0 in [[s0]] reaches F at some time
    t <= tbar(s0) := floor(log(V_upper(s0) / l_min) / log(gamma)),   l_min = min_X l.
At this time l(x_t) <= gamma^(-t) V_upper(s0) < 0, and before it the trajectory does not
reach T: the proof gives r(x_tau) <= gamma^(-tau) V_upper(s0) < 0 for tau < t.

Data: the saved reach-avoid run of Section 9.2.5 (Dubins car, 120^3 cells, gamma = 0.96,
U = A = {-1, 0, 1}). As in the solver and the paper's experiments, leaving the grid along
x or y moves the state to an absorbing failure state. The solver gives this sink the fixed
value -1 (its l and r bounds are -1 as well, cpp/include/dhj/core/vi.hpp); here it is an
absorbing state with l = r = -1, whose value is -1. Hence l_min = min(min_X l, -1) = -1.3,
attained at the obstacle center.

For a trajectory, t_hit is the first time with gamma^t l(x_t) <= V_upper(s0), the time t
of the proof. A trajectory satisfies the claim iff t_hit exists, t_hit <= tbar(s0), and
gamma^tau r(x_tau) <= V_upper(s0) for every tau < t_hit (this implies r(x_tau) < 0,
i.e. x_tau is not in T).

Checks:
  1. upper certificate T_bar V_upper <= V_upper on the verifier's own abstraction (as in
     verify_quantitative.py); with A = U it gives V_gamma <= V_upper;
  2. simulation from the start states (one state drawn uniformly in every cell with
     V_upper < 0 and, with --corners, the 8 vertices of every such cell) under each
     constant action, uniformly random actions (several seeds), the policy of the
     proposition's first claim (argmax_a min_{s' in Delta(s,a)} V_lower(s')), and the
     one-step greedy policy on V_upper (argmax_a V_upper of the cell containing f(x, a),
     with the sink counted as -1);
  3. exhaustive search from the start states (all of them with --exhaustive -1, otherwise
     N drawn at random): every action sequence is followed until it meets t_hit. A branch
     that reaches tbar(s0) without meeting it, enters the sink too late to meet it, or has
     gamma^tau r(x_tau) > V_upper(s0) before meeting it is a violation. From these states
     the search covers every control sequence, hence every policy, stationary or not.

gamma^t is formed by repeated multiplication, as in the solver: some cells attain the
bound exactly (V_upper(s0) = -gamma^k for trajectories that leave the grid at step k),
and pow() rounds differently by one ulp. Hits within TOL = 1e-12 are counted as ties and
reported; with repeated multiplication there are none.

Usage (the run reported in the paper and in README.md, about 20 s):
    python verify_violation.py RUN_DIR --out verification_ra_violation \
        --seed 0 --corners --exhaustive -1 --random-seeds 3
"""
import argparse
import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_quantitative import (SINK, Grid, certificates, check_flow, flow, l_fun,  # noqa: E402
                                 load_run, outside_grid, r_fun, wrap)

TOL = 1e-12   # ties: the solver forms gamma^k by repeated multiplication, so bounds can be
              # attained exactly; a hit within TOL is counted (and reported) as a tie


def gamma_powers(gamma, n):
    """gamma^0..gamma^n by repeated multiplication, the rounding the solver uses."""
    g = np.empty(n + 1)
    g[0] = 1.0
    for k in range(1, n + 1):
        g[k] = g[k - 1] * gamma
    return g


def l_ext(dyn, x, y, out):
    return np.where(out, SINK, l_fun(dyn, x, y))


def r_ext(dyn, x, y, out):
    return np.where(out, SINK, r_fun(dyn, x, y))


def step_all(dyn, actions, tau, X, Y, TH, a):
    """One step of the exact flow; a holds action indices (one per state)."""
    X1, Y1, TH1 = X.copy(), Y.copy(), TH.copy()
    for k, u in enumerate(actions):
        sel = np.flatnonzero(a == k)
        if sel.size:
            X1[sel], Y1[sel], th = flow(dyn["name"], dyn["v_const"], X[sel], Y[sel], TH[sel], u, tau)
            TH1[sel] = wrap(th)
    return X1, Y1, TH1


# ------------------------------------------------------------------ simulation
def simulate(grid, dyn, actions, gamma, tau, x0, vbar, tbar, choose, extra=60):
    """Follow each trajectory until t_hit (or until it is stuck in the sink without a hit,
    or tbar + extra steps). choose(t, X, Y, TH, idx) returns action indices for idx."""
    m = x0.shape[0]
    X, Y, TH = x0[:, 0].copy(), x0[:, 1].copy(), x0[:, 2].copy()
    sink = np.zeros(m, dtype=bool)
    hit_t = np.full(m, -1)
    hit_by_sink = np.zeros(m, dtype=bool)
    hit_slack = np.full(m, np.nan)       # V_upper(s0) - gamma^t l(x_t) at t_hit (>= 0)
    r_slack = np.full(m, np.inf)         # min_{tau < t_hit} V_upper(s0) - gamma^tau r(x_tau) (>= 0 required)
    in_T = np.zeros(m, dtype=bool)       # r(x_tau) > 0 for some tau < t_hit
    entered_F = np.full(m, -1)           # first time with l(x_t) <= 0 (F or sink)
    alive = np.ones(m, dtype=bool)
    tmax = int(tbar.max()) + extra
    gp = gamma_powers(gamma, tmax)
    for t in range(tmax + 1):
        ia = np.flatnonzero(alive)
        if ia.size == 0:
            break
        out = sink[ia] | outside_grid(grid, X[ia], Y[ia])
        sink[ia] = out
        lx, rx = l_ext(dyn, X[ia], Y[ia], out), r_ext(dyn, X[ia], Y[ia], out)
        g = gp[t]
        v = vbar[ia]
        first_F = (lx <= 0.0) & (entered_F[ia] < 0)
        entered_F[ia[first_F]] = t
        hit = g * lx <= v + TOL
        hit_t[ia[hit]] = t
        hit_by_sink[ia[hit]] = out[hit]
        hit_slack[ia[hit]] = (v - g * lx)[hit]
        nh = ~hit
        r_slack[ia[nh]] = np.minimum(r_slack[ia[nh]], (v - g * rx)[nh])
        in_T[ia[nh & (rx > 0.0)]] |= True
        alive[ia[hit]] = False
        alive[ia[nh & out]] = False          # in the sink without a hit: gamma^t * (-1) only increases
        keep = ia[nh & ~out]
        if keep.size == 0 or t == tmax:
            break
        a = choose(t, X, Y, TH, keep)
        X[keep], Y[keep], TH[keep] = step_all(dyn, actions, tau, X[keep], Y[keep], TH[keep], a)
    return {"hit_t": hit_t, "hit_by_sink": hit_by_sink, "hit_slack": hit_slack, "r_slack": r_slack,
            "in_T": in_T, "entered_F": entered_F}


def summarize(res, tbar):
    hit_t = res["hit_t"]
    ok_hit = hit_t >= 0
    ok_time = ok_hit & (hit_t <= tbar)
    ok_r = res["r_slack"] >= -TOL
    ok = ok_time & ok_r & ~res["in_T"]
    late = ok_hit & (hit_t > 0)
    s = {
        "trajectories": int(hit_t.size),
        "satisfy_claim": int(ok.sum()),
        "violations": int((~ok).sum()),
        "no_hit": int((~ok_hit).sum()),
        "hit_after_bound": int((ok_hit & (hit_t > tbar)).sum()),
        "r_condition_violated": int((~ok_r).sum()),
        "visited_T_before_hit": int(res["in_T"].sum()),
        "hit_at_t0": int((hit_t == 0).sum()),
        "hits_that_are_ties_(|slack|<=1e-12)": int((res["hit_slack"] < 0).sum()),
        "hit_in_obstacle": int((ok_hit & ~res["hit_by_sink"]).sum()),
        "hit_by_leaving_grid": int((ok_hit & res["hit_by_sink"]).sum()),
    }
    if late.any():
        s.update({
            "t_hit_median_(t_hit>0)": float(np.median(hit_t[late])),
            "t_hit_max": int(hit_t[ok_hit].max()),
            "tbar_median_(t_hit>0)": float(np.median(tbar[late])),
            "max_t_hit_over_tbar": float((hit_t[late] / tbar[late]).max()),
            "min_tbar_minus_t_hit": int((tbar - hit_t)[ok_hit].min()),
        })
    s["min_hit_slack"] = float(np.nanmin(res["hit_slack"])) if ok_hit.any() else None
    fin = np.isfinite(res["r_slack"])
    s["min_r_slack"] = float(res["r_slack"][fin].min()) if fin.any() else None
    return s


# ---------------------------------------------------------- exhaustive search
def exhaustive(dyn, grid, actions, gamma, tau, x0, vbar, tbar, cap):
    """All action sequences from each x0 (one search per row), each followed until t_hit.

    Returns per-state arrays; "exhausted" is False where the frontier exceeded cap."""
    m = x0.shape[0]
    out_arrays = {
        "latest_hit": np.full(m, -1), "earliest_hit": np.full(m, 1 << 30), "leaves": np.zeros(m, np.int64),
        "max_frontier": np.zeros(m, np.int64), "no_hit_by_tbar": np.zeros(m, bool),
        "sink_without_hit": np.zeros(m, bool), "r_violation": np.zeros(m, bool), "in_T": np.zeros(m, bool),
        "min_hit_slack": np.full(m, np.inf), "min_r_slack": np.full(m, np.inf), "exhausted": np.ones(m, bool),
    }

    def run(rows):
        X, Y, TH = x0[rows, 0].copy(), x0[rows, 1].copy(), x0[rows, 2].copy()
        own = np.arange(rows.size)
        acc = {k: v[rows].copy() for k, v in out_arrays.items()}
        gp = gamma_powers(gamma, int(tbar[rows].max()) + 1)
        t = 0
        while own.size:
            counts = np.bincount(own, minlength=rows.size)
            acc["max_frontier"] = np.maximum(acc["max_frontier"], counts)
            out = outside_grid(grid, X, Y)
            lx, rx = l_ext(dyn, X, Y, out), r_ext(dyn, X, Y, out)
            g = gp[t]
            v = vbar[rows][own]
            hit = g * lx <= v + TOL
            if hit.any():
                oh = own[hit]
                np.maximum.at(acc["latest_hit"], oh, t)
                np.minimum.at(acc["earliest_hit"], oh, t)
                np.add.at(acc["leaves"], oh, 1)
                np.minimum.at(acc["min_hit_slack"], oh, (v - g * lx)[hit])
            nh = ~hit
            if nh.any():
                on = own[nh]
                np.minimum.at(acc["min_r_slack"], on, (v - g * rx)[nh])
                acc["r_violation"][own[nh & (g * rx > v + TOL)]] = True
                acc["in_T"][own[nh & (rx > 0.0)]] = True
            dead = nh & out
            acc["sink_without_hit"][own[dead]] = True
            expired = nh & ~out & (t >= tbar[rows][own])
            acc["no_hit_by_tbar"][own[expired]] = True
            go = nh & ~out & ~expired
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
                out_arrays["exhausted"][rows] = False
                return
            half = rows.size // 2
            solve(rows[:half])
            solve(rows[half:])
            return
        for k in out_arrays:
            out_arrays[k][rows] = acc[k]

    batch = 4096
    for b0 in range(0, m, batch):
        solve(np.arange(b0, min(m, b0 + batch)))
        if (b0 // batch) % 20 == 0:
            print(f"  exhaustive: {min(m, b0 + batch)}/{m} states", flush=True)
    return out_arrays


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--random-seeds", type=int, default=3, help="runs with uniformly random actions")
    ap.add_argument("--exhaustive", type=int, default=2000, help="states searched exhaustively (-1: all)")
    ap.add_argument("--cap", type=int, default=12_000_000, help="largest frontier per search batch")
    ap.add_argument("--skip-sim", action="store_true")
    ap.add_argument("--corners", action="store_true", help="also start from the 8 vertices of every cell")
    a = ap.parse_args()
    t0 = time.time()
    out_dir = Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    params, data, phases = load_run(Path(a.run_dir))
    dyn, args = params["dynamics"], params["args"]
    gamma, tau = float(args["gamma"]), float(dyn["tau"])
    if args["mode"] != "ra_discount" or not (0.0 < gamma < 1.0):
        raise ValueError("needs a discounted reach-avoid run")
    actions = [float(u) for u in dyn["actions"]]
    grid = Grid(dyn, data)
    # l_min over the state box (obstacle center inside the box) and the sink.
    ox, oy = dyn["obstacle_position"]
    assert grid.lo[0] <= ox <= grid.hi[0] and grid.lo[1] <= oy <= grid.hi[1]
    l_min = min(-float(dyn["obstacle_radius"]), SINK)
    report = {"run_dir": str(a.run_dir), "gamma": gamma, "tau": tau, "actions": actions, "l_min": l_min,
              "seed": a.seed}
    report["flow_check_max_error"] = check_flow(dyn["name"], dyn["v_const"], actions, tau, rng)

    red = grid.Vu < 0.0
    report["cells"] = int(grid.Vu.size)
    report["cells_V_upper_negative"] = int(red.sum())
    report["cells_V_upper_zero"] = int((grid.Vu == 0.0).sum())
    report["min_V_upper"] = float(grid.Vu.min())
    if grid.Vu.min() < l_min:
        raise ValueError("V_upper < l_min somewhere; the time bound is undefined there")

    print("certificates ...", flush=True)
    pi_hat, cert = certificates(grid, dyn, actions, gamma, tau)
    report["upper_certificate"] = {k: cert[k] for k in ("upper_cells_checked", "upper_violations", "upper_max_excess")}

    cells = np.argwhere(red)
    x0 = grid.lo + (cells + rng.uniform(size=cells.shape)) * grid.h
    kind = np.zeros(cells.shape[0], dtype=np.int8)          # 0: uniform sample, 1: cell vertex
    if a.corners:
        offs = np.array(list(itertools.product((0, 1), repeat=3)))
        verts = (grid.lo + (cells[:, None, :] + offs[None, :, :]) * grid.h).reshape(-1, 3)
        x0 = np.concatenate([x0, verts])
        cells = np.concatenate([cells, np.repeat(cells, len(offs), axis=0)])
        kind = np.concatenate([kind, np.ones(verts.shape[0], dtype=np.int8)])
    report["start_states"] = {"uniform": int((kind == 0).sum()), "vertices": int((kind == 1).sum())}
    vbar = grid.Vu[cells[:, 0], cells[:, 1], cells[:, 2]]
    tbar = np.floor(np.log(vbar / l_min) / math.log(gamma)).astype(int)
    report["tbar"] = {"median": float(np.median(tbar)), "max": int(tbar.max()), "min": int(tbar.min()),
                      "quantiles_10_50_90_99": [float(q) for q in np.quantile(tbar, [0.1, 0.5, 0.9, 0.99])]}
    print(json.dumps({k: report[k] for k in ("cells_V_upper_negative", "min_V_upper", "tbar", "upper_certificate")}),
          flush=True)

    def const(k):
        return lambda t, X, Y, TH, idx: np.full(idx.size, k)

    def random_policy(seed):
        r = np.random.default_rng(seed)
        return lambda t, X, Y, TH, idx: r.integers(0, len(actions), size=idx.size)

    def prop4_policy(t, X, Y, TH, idx):
        ci = grid.locate([X[idx], Y[idx], TH[idx]])
        return pi_hat[ci[0], ci[1], ci[2]].astype(int)

    def greedy_upper(t, X, Y, TH, idx):
        vals = []
        for u in actions:
            x1, y1, th1 = flow(dyn["name"], dyn["v_const"], X[idx], Y[idx], TH[idx], u, tau)
            out = outside_grid(grid, x1, y1)
            ci = grid.locate([np.clip(x1, grid.lo[0], grid.hi[0]), np.clip(y1, grid.lo[1], grid.hi[1]), wrap(th1)])
            vals.append(np.where(out, SINK, grid.Vu[ci[0], ci[1], ci[2]]))
        return np.argmax(np.stack(vals), axis=0)

    policies = {f"constant u={u:g}": const(k) for k, u in enumerate(actions)}
    for j in range(a.random_seeds):
        policies[f"random actions (seed {1000 + j})"] = random_policy(1000 + j)
    policies["Prop. 4 policy on V_lower"] = prop4_policy
    policies["greedy on V_upper"] = greedy_upper

    sims = {}
    worst_t = np.full(cells.shape[0], -1)
    any_viol = np.zeros(cells.shape[0], bool)
    if not a.skip_sim:
        for name, pol in policies.items():
            ts = time.time()
            res = simulate(grid, dyn, actions, gamma, tau, x0, vbar, tbar, pol)
            sims[name] = summarize(res, tbar)
            sims[name]["seconds"] = round(time.time() - ts, 1)
            worst_t = np.maximum(worst_t, np.where(res["hit_t"] >= 0, res["hit_t"], 1 << 30))
            any_viol |= ~((res["hit_t"] >= 0) & (res["hit_t"] <= tbar) & (res["r_slack"] >= -TOL) & ~res["in_T"])
            print(name, json.dumps(sims[name]), flush=True)
        late = worst_t > 0
        report["simulation"] = sims
        report["simulation_all_policies"] = {
            "states": int(cells.shape[0]),
            "states_violated_by_some_policy": int(any_viol.sum()),
            "max_over_policies_t_hit_over_tbar": float((worst_t[late] / tbar[late]).max()) if late.any() else None,
            "median_worst_t_hit_(t_hit>0)": float(np.median(worst_t[late])) if late.any() else None,
        }

    # exhaustive search on a random subsample (or all states)
    n_ex = cells.shape[0] if a.exhaustive < 0 else min(a.exhaustive, cells.shape[0])
    if n_ex > 0:
        pick = np.sort(rng.choice(cells.shape[0], size=n_ex, replace=False))
        ts = time.time()
        ex = exhaustive(dyn, grid, actions, gamma, tau, x0[pick], vbar[pick], tbar[pick], a.cap)
        done = ex["exhausted"]
        bad = done & (ex["no_hit_by_tbar"] | ex["sink_without_hit"] | ex["r_violation"] | ex["in_T"])
        lh, tb = ex["latest_hit"][done], tbar[pick][done]
        pos = lh > 0
        report["exhaustive"] = {
            "states": int(n_ex),
            "exhausted": int(done.sum()),
            "not_exhausted_(frontier_cap)": int((~done).sum()),
            "violations": int(bad.sum()),
            "no_hit_by_tbar": int((done & ex["no_hit_by_tbar"]).sum()),
            "sink_without_hit": int((done & ex["sink_without_hit"]).sum()),
            "r_violation": int((done & ex["r_violation"]).sum()),
            "visited_T_before_hit": int((done & ex["in_T"]).sum()),
            "action_prefixes_checked": int(ex["leaves"][done].sum()),
            "max_frontier": int(ex["max_frontier"].max()),
            "latest_hit_max": int(lh.max()) if lh.size else None,
            "latest_hit_median_(>0)": float(np.median(lh[pos])) if pos.any() else None,
            "tbar_median_(latest_hit>0)": float(np.median(tb[pos])) if pos.any() else None,
            "max_latest_hit_over_tbar": float((lh[pos] / tb[pos]).max()) if pos.any() else None,
            "min_tbar_minus_latest_hit": int((tb - lh).min()) if lh.size else None,
            "min_hit_slack": float(ex["min_hit_slack"][done].min()) if done.any() else None,
            "states_with_a_tie_hit": int((done & (ex["min_hit_slack"] < 0)).sum()),
            "min_r_slack": float(ex["min_r_slack"][done][np.isfinite(ex["min_r_slack"][done])].min())
            if np.isfinite(ex["min_r_slack"][done]).any() else None,
            "seconds": round(time.time() - ts, 1),
        }
        np.savez_compressed(out_dir / "exhaustive.npz", x0=x0[pick], vbar=vbar[pick], tbar=tbar[pick],
                            cells=cells[pick], **ex)
        print("exhaustive", json.dumps(report["exhaustive"]), flush=True)

    np.savez_compressed(out_dir / "sampled_states.npz", x0=x0, vbar=vbar, tbar=tbar, cells=cells,
                        worst_t_hit=worst_t, any_violation=any_viol)
    report["seconds"] = round(time.time() - t0, 1)
    (out_dir / "report.json").write_text(json.dumps(report, indent=2))
    print(f"wrote {out_dir}", flush=True)


if __name__ == "__main__":
    main()
