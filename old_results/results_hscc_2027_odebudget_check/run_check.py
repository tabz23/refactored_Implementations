"""Rerun each modified-init experiment for at most 5 minutes and compare."""
import os
import re
import signal
import subprocess
import time
from pathlib import Path

CPP = Path("/Users/i.k.tabbara/Documents/python directory/HSCC_2027/refactored_Implementations/cpp")
OLD = Path("/Users/i.k.tabbara/Documents/python directory/HSCC_2027/refactored_Implementations/results_hscc_2027_modified_init")
NEW = Path("/Users/i.k.tabbara/Documents/python directory/HSCC_2027/refactored_Implementations/results_hscc_2027_odebudget_check")
TIMEOUT = 300

# One command per results folder. Dubins discount was launched twice with the
# same solver flags; only the plot settings differed, so one --no-plot run covers it.
RUNS = [
    ["--mode", "ra_discount", "--dynamics", "van_der_pol", "--dt", "0.01", "--tau", "0.01",
     "--initial-resolution", "40", "--epsilon", "0.001", "--delta-min", "1e-8", "--phase0-tol", "1e-8",
     "--delta-max", "1e-8", "--gamma", "0.96", "--refinements", "18", "--vi-iterations", "20000"],
    ["--mode", "ra_discount", "--dynamics", "evasion", "--dt", "0.3", "--tau", "0.3",
     "--initial-resolution", "30", "--epsilon", "0.05", "--delta-min", "1e-8", "--phase0-tol", "1e-8",
     "--delta-max", "1e-8", "--gamma", "0.96", "--refinements", "12", "--vi-iterations", "1000033"],
    ["--mode", "ra_discount", "--dynamics", "double_integrator", "--dt", "0.3", "--tau", "0.3",
     "--initial-resolution", "10", "--epsilon", "0.05", "--delta-min", "1e-8", "--phase0-tol", "1e-8",
     "--delta-max", "1e-8", "--gamma", "0.96", "--refinements", "14", "--vi-iterations", "20000"],
    ["--mode", "ra_discount", "--dynamics", "bicycle", "--a-max", "0.2", "--v-max", "1",
     "--dt", "0.3", "--tau", "0.3", "--initial-resolution", "5", "--epsilon", "0.1",
     "--delta-min", "1e-8", "--phase0-tol", "1e-8", "--delta-max", "1e-8", "--gamma", "0.96",
     "--refinements", "8", "--vi-iterations", "20000"],
    ["--mode", "ra_discount", "--dynamics", "dubins", "--dt", "0.3", "--tau", "0.3",
     "--initial-resolution", "20", "--epsilon", "0.05", "--delta-min", "1e-8", "--phase0-tol", "1e-8",
     "--delta-max", "1e-8", "--gamma", "0.96", "--refinements", "13", "--vi-iterations", "100033"],
    ["--mode", "ra_nodiscount", "--dynamics", "van_der_pol", "--dt", "0.01", "--tau", "0.01",
     "--initial-resolution", "40", "--epsilon", "0.001", "--delta-min", "1e-8", "--phase0-tol", "1e-8",
     "--delta-max", "1e-8", "--gamma", "0.96", "--refinements", "18", "--vi-iterations", "20000"],
    ["--mode", "ra_nodiscount", "--dynamics", "evasion", "--dt", "0.3", "--tau", "0.3",
     "--initial-resolution", "30", "--epsilon", "0.05", "--delta-min", "1e-8", "--phase0-tol", "1e-8",
     "--delta-max", "1e-8", "--gamma", "0.96", "--refinements", "12", "--vi-iterations", "1000033"],
    ["--mode", "ra_nodiscount", "--dynamics", "double_integrator", "--dt", "0.3", "--tau", "0.3",
     "--initial-resolution", "10", "--epsilon", "0.05", "--delta-min", "1e-8", "--phase0-tol", "1e-8",
     "--delta-max", "1e-8", "--gamma", "0.9", "--refinements", "14", "--vi-iterations", "20000"],
    ["--mode", "ra_nodiscount", "--dynamics", "bicycle", "--a-max", "0.2", "--v-max", "1",
     "--dt", "0.3", "--tau", "0.3", "--initial-resolution", "5", "--epsilon", "0.1",
     "--delta-min", "1e-8", "--phase0-tol", "1e-8", "--delta-max", "1e-8", "--gamma", "0.9",
     "--refinements", "8", "--vi-iterations", "20000"],
    ["--mode", "ra_nodiscount", "--dynamics", "dubins", "--dt", "0.3", "--tau", "0.3",
     "--initial-resolution", "20", "--epsilon", "0.05", "--delta-min", "1e-8", "--phase0-tol", "1e-8",
     "--delta-max", "1e-8", "--gamma", "0.9", "--refinements", "13", "--vi-iterations", "100033"],
    ["--mode", "avoid_nodiscount", "--dynamics", "dubins", "--dt", "0.3", "--tau", "0.3",
     "--initial-resolution", "80", "--epsilon", "0.05", "--delta-min", "1e-8", "--phase0-tol", "1e-8",
     "--delta-max", "1e-8", "--gamma", "0.9", "--refinements", "6", "--vi-iterations", "20000"],
]

CLS = re.compile(
    r"Classification after phase (\d+):\n"
    r"\s*Safe:\s+(\d+) cells.*\n"
    r"\s*Unsafe:\s+(\d+) cells.*\n"
    r"\s*Unclassified:\s+(\d+) cells",
    re.M,
)


def last_log(text):
    parts = text.split("# run started ")
    return parts[-1] if parts else text


def classifications(text):
    return {int(m.group(1)): tuple(int(x) for x in m.group(2, 3, 4)) for m in CLS.finditer(last_log(text))}


def phase_rows(path):
    lines = path.read_text().splitlines()
    header, rows = lines[0].split(","), lines[1:]
    parsed = [dict(zip(header, r.split(","))) for r in rows if r.strip()]
    starts = [i for i, r in enumerate(parsed) if r["phase"] == "0"]
    parsed = parsed[starts[-1]:] if starts else parsed
    keys = ["n_leaves", "boundary_cells", "refinable_cells", "refined_parents", "vi_iterations", "converged"]
    return {int(r["phase"]): tuple(r[k] for k in keys) for r in parsed}


def find_old(mode, dynamics):
    hits = list((OLD / mode).glob(dynamics + "_*"))
    if len(hits) != 1:
        raise SystemExit(f"expected one old folder for {mode}/{dynamics}, found {hits}")
    return hits[0]


def main():
    summaries = []
    for args in RUNS:
        mode = args[args.index("--mode") + 1]
        dynamics = args[args.index("--dynamics") + 1]
        label = f"{mode}/{dynamics}"
        cmd = ["./build/dhj", *args, "--checkpoint-every", "0", "--keep-checkpoints", "0",
               "--no-plot", "--workers", "12", "--results-root", str(NEW)]
        print(f"\n=== START {label} ===", flush=True)
        t0 = time.time()
        proc = subprocess.Popen(cmd, cwd=CPP, start_new_session=True)
        timed_out = False
        try:
            code = proc.wait(timeout=TIMEOUT)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(proc.pid, signal.SIGKILL)
            code = proc.wait()
        elapsed = time.time() - t0
        print(f"=== END {label} elapsed={elapsed:.1f}s timeout={timed_out} code={code} ===", flush=True)
        summaries.append((label, timed_out, elapsed, code))

    print("\n=== COMPARE ===", flush=True)
    any_diff = False
    for args in RUNS:
        mode = args[args.index("--mode") + 1]
        dynamics = args[args.index("--dynamics") + 1]
        label = f"{mode}/{dynamics}"
        old_dir = find_old(mode, dynamics)
        new_dirs = list((NEW / mode).glob(dynamics + "_*"))
        if len(new_dirs) != 1:
            print(f"{label}: MISSING new output {new_dirs}")
            any_diff = True
            continue
        new_dir = new_dirs[0]
        old_ph, new_ph = phase_rows(old_dir / "phases.csv"), phase_rows(new_dir / "phases.csv")
        old_cls = classifications((old_dir / "run.log").read_text(errors="replace"))
        new_cls = classifications((new_dir / "run.log").read_text(errors="replace"))
        common = sorted(set(old_ph) & set(new_ph))
        grid_bad, cls_bad, iter_bad = [], [], []
        for p in common:
            if old_ph[p][:4] != new_ph[p][:4]:
                grid_bad.append(p)
            if old_ph[p][4] != new_ph[p][4] or old_ph[p][5] != new_ph[p][5]:
                iter_bad.append((p, old_ph[p][4], new_ph[p][4]))
            if p in old_cls and p in new_cls and old_cls[p] != new_cls[p]:
                cls_bad.append((p, old_cls[p], new_cls[p]))
            elif p in new_cls and p not in old_cls:
                cls_bad.append((p, None, new_cls[p]))
        status = "SAME" if not grid_bad and not cls_bad else "DIFFERENT"
        if grid_bad or cls_bad:
            any_diff = True
        print(f"{label}: {status}  compared_phases={common[0]}..{common[-1] if common else '-'} "
              f"({len(common)})  new_only_phases={sorted(set(new_ph)-set(old_ph))}  "
              f"grid_mismatch={grid_bad}  class_mismatch={cls_bad}  "
              f"iteration_mismatch={iter_bad}")
    print("COMPARE DONE" + (" WITH DIFFERENCES" if any_diff else " ALL MATCHED"), flush=True)


if __name__ == "__main__":
    main()
