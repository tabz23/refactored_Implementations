"""Numerical check of the last claim of Theorem 4 (paper Appendix H, subsection "Bound on unclassified cells", label app:bv_unclassified).

Claim: on a uniform grid with cell radius eta, every cell s with V_lower(s) <= 0 < V_upper(s)
(outputs of Algorithm 1) satisfies |V_gamma(x)| <= (2 L eta + gamma delta*)/(1 - gamma) for all
x in [[s]], where delta* = max of the two stopping tolerances (both 1e-12 here).

For every (gamma, N) of bound_verification.py, Algorithm 1 is rerun with that code, and
  (1) sound check: by Theorem 2, V_lower(s) <= V_gamma(x) <= V_upper(s), so on an unclassified
      cell |V_gamma(x)| <= max(V_upper(s), -V_lower(s)). The largest such value over all
      unclassified cells is compared with the bound. This covers every state of every cell.
  (2) direct check, independent of the abstraction: at random states in random unclassified
      cells (and in the two cells that attain the maxima of (1)), V_gamma(x) is bounded from
      both sides by enumerating all 3^H action sequences with the exact flow:
          V_H(x) <= V_gamma(x) <= max(V_H(x), gamma^(H+1) r_max),
      where V_H(x) is the best discounted RA measure over times 0..H and r_max = 0.5 = max r.
      H is the smallest horizon with gamma^(H+1) r_max <= 0.01.
Usage: python unclassified_check.py [SAMPLES_PER_RUN] [N ...]   (paper: 200 samples, all N)
Output of the paper run: unclassified_check.log in this folder.
"""
import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
SRC = Path(__file__).resolve().parent      # folder of bound_verification.py and summary.json
K = int(sys.argv[1]) if len(sys.argv) > 1 else 200
NS = [int(a) for a in sys.argv[2:]] or [50, 60, 70, 80, 90, 100, 110, 120]
spec = importlib.util.spec_from_file_location("bv", SRC / "bound_verification.py")
bv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bv)
ref = {(r["gamma"], r["N"]): r for r in json.load(open(SRC / "summary.json"))["runs"]}
GAMMAS = [0.3, 0.4, 0.5, 0.6, 0.7]
R_MAX = bv.R_TARGET          # r = -(d - 0.5) <= 0.5
TRUNC = 0.01


def l_fn(x, y):
    return bv.periodic_distance(x, y, bv.OBSTACLE) - bv.R_OBSTACLE


def r_fn(x, y):
    return -(bv.periodic_distance(x, y, bv.TARGET) - bv.R_TARGET)


def horizon(gamma):
    return math.ceil(math.log(TRUNC / R_MAX) / math.log(gamma)) - 1


def v_H(states, gamma, H, batch):
    """max over all action sequences of length H of max_{t<=H} min(g^t r(x_t), min_{tau<=t} g^tau l(x_tau))."""
    out = np.empty(len(states))
    for b0 in range(0, len(states), batch):
        s = states[b0:b0 + batch]
        x, y, th = s[:, 0:1], s[:, 1:2], s[:, 2:3]
        m = l_fn(x, y)
        best = np.minimum(r_fn(x, y), m).max(axis=1)
        for t in range(1, H + 1):
            nxt = [bv.flow(x, y, th, u, bv.TAU) for u in bv.ACTIONS]
            x = np.concatenate([n[0] for n in nxt], axis=1)
            y = np.concatenate([n[1] for n in nxt], axis=1)
            th = np.concatenate([n[2] for n in nxt], axis=1)
            g = gamma ** t
            m = np.minimum(np.concatenate([m] * len(bv.ACTIONS), axis=1), g * l_fn(x, y))
            best = np.maximum(best, np.minimum(g * r_fn(x, y), m).max(axis=1))
        out[b0:b0 + batch] = best
    return out


rng = np.random.default_rng(0)
rows = []
for N in NS:
    t0 = time.time()
    eta, bounds, tables = bv.grid(N)
    h = (bv.HI - bv.LO) / N
    for g in GAMMAS:
        Vl, Vu, iters = bv.value_iteration(g, bounds, tables)
        r0 = ref[(g, N)]
        same = iters == r0["iterations"] and abs(float((Vu - Vl).max()) - r0["gap_max"]) <= 1e-12
        bound = (2.0 * bv.L * eta + g * bv.TOL) / (1.0 - g)
        unc = (Vl <= 0) & (Vu > 0)
        idx = np.flatnonzero(unc)
        up, lo = Vu.ravel()[idx], -Vl.ravel()[idx]
        sound = float(max(up.max(), lo.max()))
        # (2) sample states: the two extreme cells of (1) plus K random unclassified cells
        pick = np.concatenate([[idx[up.argmax()], idx[lo.argmax()]],
                               rng.choice(idx, size=min(K, idx.size), replace=False)])
        ijk = np.stack(np.unravel_index(pick, (N, N, N)), axis=1)
        states = bv.LO + (ijk + rng.random(ijk.shape)) * h
        H = horizon(g)
        vh = v_H(states, g, H, batch=max(1, 2_000_000 // 3 ** H))
        tail = g ** (H + 1) * R_MAX
        vmax = np.maximum(vh, tail)                     # V_gamma(x) <= vmax, V_gamma(x) >= vh
        sl, su = Vl.ravel()[pick], Vu.ravel()[pick]
        row = {"gamma": g, "N": N, "same_as_appendix_H": bool(same), "unclassified": int(idx.size),
               "bound": bound, "sound_max_over_cells": sound, "sound_ratio": sound / bound,
               "H": H, "tail": tail, "samples": int(len(pick)),
               "direct_max_upper": float(vmax.max()), "direct_max_neg_lower": float((-vh).max()),
               "direct_ratio": float(max(vmax.max(), (-vh).max()) / bound),
               "direct_ok": bool((vmax <= bound).all() and (vh >= -bound).all()),
               "soundness_violations": int(((vh > su + 1e-9) | (sl > vmax + 1e-9)).sum()),
               "overlap_min_slack": float(min((su - vh).min(), (vmax - sl).min())),
               "direct_min_lower": float(vh.min()), "direct_max_vh": float(vh.max()),
               "seconds": round(time.time() - t0, 1)}
        rows.append(row)
        print(json.dumps(row), flush=True)
        t0 = time.time()
    del tables
print("summary:", json.dumps({
    "runs": len(rows), "states": sum(r["samples"] for r in rows),
    "unclassified_range": [min(r["unclassified"] for r in rows), max(r["unclassified"] for r in rows)],
    "sound_ratio_range": [min(r["sound_ratio"] for r in rows), max(r["sound_ratio"] for r in rows)],
    "direct_ratio_range": [min(r["direct_ratio"] for r in rows), max(r["direct_ratio"] for r in rows)],
    "overlap_min_slack": min(r["overlap_min_slack"] for r in rows),
    "all_same_as_appendix_H": all(r["same_as_appendix_H"] for r in rows)}), flush=True)
ok = all(r["same_as_appendix_H"] and r["sound_ratio"] <= 1 and r["direct_ok"] and r["soundness_violations"] == 0
         for r in rows)
print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED", flush=True)
