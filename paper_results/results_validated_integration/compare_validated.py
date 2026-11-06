#!/usr/bin/env python3
"""Compare the --validated runs of run_validated.sh with the stored paper runs.

For each case, reads the solver's own comparison of the validated reach boxes with the boxes of
the normal build (the last "Validated integration (CAPD)" block of logs/<case>.log), and checks
that the run reproduces the stored paper run: the phases.csv rows (phase, cells, unclassified
cells, refinable cells, refined parents, value-iteration sweeps, converged), the classification
counts printed in run.log, and, where both runs dumped them, the per-cell value files.

Usage: python compare_validated.py [--root ROOT] [case ...]   (writes ROOT/summary.md)
Exit status 1 if any check fails.
"""
import argparse
import filecmp
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PR = os.path.dirname(HERE)  # paper_results
QV = os.path.join(PR, "results_quantitative_verification")
FIG = os.path.join(PR, "results_hscc_2027_newbellmans_12workers_current_code")
CASES = {  # case: (paper result, stored run, mode, run-folder prefix)
    "qv_ra": ("Sec. 9.2.4 reach-avoid, Fig. 3(b)",
              QV + "/ra_discount/dubins_res120_dlower1e-08_dupper1e-08_gamma0.96_dt0.3", "ra_discount", "dubins_res120_"),
    "qv_avoid": ("Sec. 9.2.4 avoid",
                 QV + "/avoid_nodiscount/dubins_res120_dlower1e-300_dupper1e-300_dt0.3", "avoid_nodiscount", "dubins_res120_"),
    "dubins120_std": ("Fig. 3(a)",
                      FIG + "/ra_nodiscount/dubins_res120_dlower1e-15_dupper1e-15_dt0.3", "ra_nodiscount", "dubins_res120_"),
    "avoid80": ("Fig. 4",
                FIG + "/avoid_nodiscount/dubins_res80_epsmin1e-06_dlower1e-15_dupper1e-15_dt0.3", "avoid_nodiscount", "dubins_res80_"),
    "dubins20": ("Fig. 1",
                 FIG + "/ra_discount/dubins_res20_epsmin1e-06_dlower1e-08_dupper1e-08_gamma0.96_dt0.3", "ra_discount", "dubins_res20_"),
    "bicycle8": ("Fig. 2(d)",
                 FIG + "/ra_nodiscount/bicycle_res5_epsmin1e-06_dlower1e-15_dupper1e-15_dt0.3", "ra_nodiscount", "bicycle_res5_"),
    "vdp18": ("Fig. 2(a)",
              FIG + "/ra_nodiscount/van_der_pol_res40_epsmin1e-06_dlower1e-15_dupper1e-15_dt0.01", "ra_nodiscount", "van_der_pol_res40_"),
    "evasion12": ("Fig. 2(b)",
                  FIG + "/ra_nodiscount/evasion_res30_epsmin1e-06_dlower1e-15_dupper1e-15_dt0.3", "ra_nodiscount", "evasion_res30_"),
    "di14": ("Fig. 2(c)",
             FIG + "/ra_nodiscount/double_integrator_res10_epsmin1e-06_dlower1e-15_dupper1e-15_dt0.3", "ra_nodiscount", "double_integrator_res10_"),
}


# A run folder that was reused appends every run to phases.csv and run.log (the stored Section
# 9.2.4 folders hold three identical runs). Only the last run is compared.
def phases(d):
    rows = open(os.path.join(d, "phases.csv")).read().strip().split("\n")[1:]
    starts = [i for i, r in enumerate(rows) if r.split(",")[0] == "0"]
    rows = rows[starts[-1]:] if starts else rows
    return [",".join(r.split(",")[:7]) for r in rows]


def classification(d):
    text = open(os.path.join(d, "run.log"), errors="replace").read().split("# run started")[-1]
    return re.findall(r"^\s*(Safe|Unsafe|Unclassified):\s+(\d+) cells", text, re.M)


def capd_stats(log):
    text = open(log, errors="replace").read()
    blocks = text.split("Validated integration (CAPD), all transitions so far:")
    if len(blocks) < 2:
        return None
    b = blocks[-1]
    g = lambda pat: re.search(pat, b).group(1)
    return dict(transitions=int(g(r"transitions: (\d+)")), failures=int(g(r"CAPD failures: (\d+)")),
                sets_differ=int(g(r"differ from the budget boxes: (\d+)")),
                oob_differ=int(g(r"out-of-domain flags that differ: (\d+)")),
                max_cert=float(g(r"certified error: max ([0-9.e+-]+)")),
                max_budget=float(g(r"budget: max ([0-9.e+-]+)")),
                above=int(g(r"certified > budget at (\d+)")), ratio=float(g(r"ratio: ([0-9.e+-]+)")),
                retries=int(m.group(1)) if (m := re.search(r"fixed-step retries: (\d+)", b)) else 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=HERE)
    ap.add_argument("cases", nargs="*")
    a = ap.parse_args()
    rows, ok_all = [], True
    for name in a.cases or CASES:
        paper, stored, mode, prefix = CASES[name]
        runs = glob.glob(os.path.join(a.root, "runs", mode, prefix + "*_validated"))
        log = os.path.join(a.root, "logs", name + ".log")
        s = capd_stats(log) if os.path.exists(log) else None
        if not runs or s is None:
            rows.append(f"| {name} | {paper} | not run or not finished | | | | | |")
            ok_all = False
            continue
        v = runs[0]
        pv, ps = phases(v), phases(stored)
        cv, cs = classification(v), classification(stored)
        same_phases = len(pv) == len(ps) and pv == ps
        same_class = len(cv) > 0 and cv == cs[:len(cv)] and len(cv) == len(cs)
        files = []
        for f in ("value_function_phase_0.csv", "upper_nonpos_iteration_phase_0.csv"):
            x, y = os.path.join(stored, f), os.path.join(v, f)
            if os.path.exists(x) and os.path.exists(y):
                files.append(f"{f}: {'identical' if filecmp.cmp(x, y, shallow=False) else 'DIFFERENT'}")
        ok = (s["failures"] == 0 and s["sets_differ"] == 0 and s["oob_differ"] == 0 and s["above"] == 0
              and same_phases and same_class and all(f.endswith("identical") for f in files))
        ok_all &= ok
        rows.append(f"| {name} | {paper} | {s['transitions']:,} | {s['failures']} ({s['retries']:,} retried) | {s['above']} | "
                    f"{s['max_cert']:.1e} / {s['max_budget']:.1e} | {s['sets_differ']} + {s['oob_differ']} | "
                    f"{'yes' if same_phases and same_class else 'NO'} ({len(pv)} phases)"
                    + (f", {'; '.join(files)}" if files else "") + f" | {'PASS' if ok else 'FAIL'} |")
    header = ("| case | paper | transitions | CAPD failures | validated error > budget | "
              "max validated error / max budget | successor sets + out-of-domain flags changed | "
              "same phases and classification as the paper run | result |\n|---|---|---|---|---|---|---|---|---|")
    table = header + "\n" + "\n".join(rows)
    print(table)
    with open(os.path.join(a.root, "summary.md"), "w") as f:
        f.write(table + "\n")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
