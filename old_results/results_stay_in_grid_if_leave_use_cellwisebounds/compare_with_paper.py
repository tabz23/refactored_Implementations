"""Compare the runs in this folder with the paper's runs, phase by phase, and say whether each
one matches the corresponding run in ../results_stay_in_grid (same classifications in every
phase).

For every case, prints the share of the state space certified Safe (V_lower > 0: in R, or
outside the BRS in the avoid case), Unsafe (V_upper <= 0: in L, or in the BRS) and
Unclassified, in the paper's run (with the sink state) and in the run here (reach sets
restricted to the grid), and how many reach boxes leave the state box.

    python compare_with_paper.py [PAPER_FIGS_FOLDER]
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
def find_up(name):
    """The nearest folder called `name` in this folder or above it."""
    for d in [HERE, *HERE.parents]:
        if (d / name).is_dir():
            return d / name
    raise SystemExit(f"{name}/ not found above {HERE}")


PAPER = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    find_up("paper_results") / "results_hscc_2027_newbellmans_12workers_current_code"
CASES = [  # (paper figure, mode, folder prefix)
    ("Fig. 4(a) dubins120_std", "ra_nodiscount", "dubins_res120_"),
    ("Fig. 4(b) dubins120_discount", "ra_discount", "dubins_res120_"),
    ("Fig. 5 avoid80", "avoid_nodiscount", "dubins_res80_"),
    ("Fig. 2 dubins20_discount", "ra_discount", "dubins_res20_"),
    ("Fig. 3(d) bicycle8", "ra_nodiscount", "bicycle_res5_"),
    ("Fig. 3(a) vdp18", "ra_nodiscount", "van_der_pol_res40_"),
    ("Fig. 3(c) di14", "ra_nodiscount", "double_integrator_res10_"),
    ("Fig. 3(b) evasion12", "ra_nodiscount", "evasion_res30_"),
]


def find(root, mode, prefix, restricted):
    runs = [d for d in sorted((root / mode).glob(prefix + "*"))
            if d.is_dir() and ("_restrictgrid" in d.name) == restricted]
    return runs[0] if runs else None


def parse(run):
    """Per phase: (safe, unsafe, unclassified) in % of the state space, and box counts."""
    text = (run / "run.log").read_text()
    cls = {}
    for m in re.finditer(r"Classification after phase (\d+):\n(.*)\n(.*)\n(.*)\n", text):
        cls[int(m.group(1))] = tuple(float(re.search(r"([\d.]+)% of state space", m.group(k)).group(1))
                                     for k in (2, 3, 4))
    boxes = [tuple(int(x) for x in m.groups()) for m in re.finditer(
        r"Reach boxes leaving the state box: (\d+) of (\d+) cell-action pairs \([^)]*\); "
        r"pairs with no successor: (\d+)", text)]
    total = re.findall(r"Total time: ([\d.]+) seconds", text)
    return cls, boxes, (float(total[-1]) if total else None)


print(f"paper runs: {PAPER}\n")
print("Safe / Unsafe / Unclassified, % of the state space (paper run with sink -> run here)")
for title, mode, prefix in CASES:
    mine, ref = find(HERE, mode, prefix, True), find(PAPER, mode, prefix, False)
    print(f"\n== {title}")
    if mine is None or not (mine / "run.log").exists():
        print("   not run yet")
        continue
    cls, boxes, total = parse(mine)
    rcls, _, rtotal = parse(ref)
    if not cls:
        print("   no completed phase yet")
        continue
    for k in sorted(cls):
        a, b = rcls.get(k), cls[k]
        d = " ".join(f"{x:+.1f}" for x in (b[0] - a[0], b[1] - a[1], b[2] - a[2])) if a else ""
        old = f"{a[0]:5.1f} {a[1]:5.1f} {a[2]:5.1f}" if a else "  (no such phase)"
        print(f"   phase {k:2d}: {old} -> {b[0]:5.1f} {b[1]:5.1f} {b[2]:5.1f}   change {d}")
    if boxes:
        left, pairs, empty = boxes[-1]
        print(f"   last phase: {left:,} of {pairs:,} cell-action reach boxes leave the state box "
              f"({100 * left / pairs:.1f}%), {empty:,} of them entirely (no successor)")
    last = max(cls)
    print(f"   phases run here: 0-{last}; paper: 0-{max(rcls)}; time here: "
          f"{f'{total:.1f} s' if total else 'stopped early'} (paper run: {rtotal:.1f} s)")
    first = find(HERE.parent / "results_stay_in_grid", mode, prefix, True)
    if first is not None and first.resolve() != mine.resolve():
        fcls = parse(first)[0]
        diff = [k for k in sorted(set(cls) | set(fcls)) if cls.get(k) != fcls.get(k)]
        print("   vs ../results_stay_in_grid: " +
              ("identical in every phase" if not diff else f"differs from phase {diff[0]} on"))
