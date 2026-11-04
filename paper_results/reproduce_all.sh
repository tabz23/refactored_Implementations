#!/bin/bash
# Reproduces every result of the paper that this code produces and compares it with the
# results stored in this folder (times on an Apple M3 Pro with 12 workers):
#   1. Figs. 2-5 and their runtimes   results_hscc_2027_newbellmans_12workers_current_code/   about 1 hour
#   2. Section 9.2.5                  results_quantitative_verification/                       about 3.5 minutes
#   3. Appendix H                     results_bound_verification/ (reproduce.sh and
#                                     unclassified_check.py)                                   about 7 minutes
#
#   ./reproduce_all.sh [--quick] [OUT]     OUT defaults to ./rerun_all (ignored by git)
#
# --quick is a smoke test of about 20 minutes. In step 1, vdp18, di14 and evasion12 stop
# after 16, 10 and 9 refinements instead of 18, 14 and 12, and their phases are compared
# with the same phases of the stored runs. The other cases and steps 2 and 3 run in full.
#
# Every comparison prints "same" or "identical" when the results match the stored ones.
# Runtimes are printed but not compared, since they depend on the machine. The script
# exits with status 1 if anything differs. Needs a C++17 compiler with zlib, and Python 3
# with numpy, scipy, matplotlib and Pillow (set PYTHON to choose the interpreter; the
# default is ~/miniconda3/bin/python). LaTeX with dvipng is optional (figure text only).
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
QUICK=0
if [[ "${1:-}" == "--quick" ]]; then QUICK=1; shift; fi
OUT="${1:-$HERE/rerun_all}"
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
export PYTHON="${PYTHON:-$HOME/miniconda3/bin/python}"
export CPP_DIR="${CPP_DIR:-$(cd "$HERE/../cpp" && pwd)}"
FIGS=results_hscc_2027_newbellmans_12workers_current_code
REGEN=regenerate_hscc_newbellmans_tol1e-15_current_code.sh
status=0
start=$(date +%s)
echo "reproduce_all: quick=$QUICK, output in $OUT, solver sources in $CPP_DIR"

# ---------------------------------------------------------------- 1. Figs. 2-5
echo "=== 1. Figs. 2-5 ($(date)) ==="
mkdir -p "$OUT/figs"
cp "$HERE/$FIGS/$REGEN" "$OUT/figs/"
QUICK=$QUICK "$OUT/figs/$REGEN" > "$OUT/figs/regenerate.log" 2>&1 || {
  echo "the step-1 script failed; see $OUT/figs/regenerate.log"; exit 1; }
"$PYTHON" -I - "$HERE/$FIGS" "$OUT/figs" > "$OUT/figs/comparison.txt" <<'EOF' || status=1
"""Compare each new run with the stored run of the same case, phase by phase."""
import csv
import re
import sys
from pathlib import Path

ref_root, new_root = map(Path, sys.argv[1:3])
MODES = ("ra_nodiscount", "ra_discount", "avoid_nodiscount")


def key(name):
    # Stored folders carry "_eps<--epsilon>", new ones "_epsmin<--eps-min>".
    return re.sub(r"_eps(min)?[^_]*", "", name)


def blocks(log):
    lines = log.read_text().splitlines()
    return {int(m.group(1)): lines[i:i + 4] for i, l in enumerate(lines)
            if (m := re.match(r"Classification after phase (\d+):", l))}


def classified(block):
    # Volume certified in R or L: the first two lines, "... y% of state space".
    return sum(float(re.search(r"([\d.]+)% of state space", l).group(1)) for l in block[1:3])


def same_png(a, b):
    if a.read_bytes() == b.read_bytes():
        return True
    try:
        from PIL import Image, ImageChops
    except ImportError:
        return False
    with Image.open(a) as x, Image.open(b) as y:
        return x.size == y.size and ImageChops.difference(x.convert("RGB"), y.convert("RGB")).getbbox() is None


refs = {(m, key(d.name)): d for m in MODES for d in (ref_root / m).glob("*") if d.is_dir()}
bad = 0
for m in MODES:
    for new in sorted(d for d in (new_root / m).glob("*") if d.is_dir()):
        ref = refs.get((m, key(new.name)))
        if ref is None:
            print(f"DIFFERENT: no stored run for {m}/{new.name}")
            bad += 1
            continue
        nb, rb = blocks(new / "run.log"), blocks(ref / "run.log")
        phases = sorted(nb)
        ok_cls = bool(phases) and all(nb[k] == rb.get(k) for k in phases)
        # phases.csv without timings. refinable_cells of a run's last phase depends on how
        # the run ended, so it is compared only for the phases before it.
        nr, rr = (list(csv.DictReader(open(d / "phases.csv"))) for d in (new, ref))
        cols = ("n_leaves", "boundary_cells", "refined_parents", "vi_iterations", "converged")
        ok_csv = (len(nr) <= len(rr) and all(a[c] == b[c] for a, b in zip(nr, rr) for c in cols)
                  and all(a["refinable_cells"] == b["refinable_cells"] for a, b in zip(nr[:-1], rr)))
        pngs = sorted(new.glob("*.png"))
        bad_png = [p.name for p in pngs if not (ref / p.name).exists() or not same_png(p, ref / p.name)]
        t = re.search(r"Total time: ([\d.]+)", (new / "run.log").read_text())
        ok = ok_cls and ok_csv and not bad_png
        bad += not ok
        last = phases[-1] if phases else -1
        print(f"{'same' if ok else 'DIFFERENT'}: {m}/{key(new.name)}: phases 0-{last} of 0-{max(rb)}, "
              f"{classified(nb[last]) if phases else float('nan'):.1f}% classified after phase {last}, "
              f"{len(pngs)} PNGs {'pixel-identical' if not bad_png else 'differ: ' + ', '.join(bad_png)}, "
              f"{t.group(1) if t else '?'} s")
        if not ok_cls:
            print("    the classification blocks in run.log differ")
        if not ok_csv:
            print("    phases.csv differs")
sys.exit(1 if bad else 0)
EOF
cat "$OUT/figs/comparison.txt"

# ---------------------------------------------------------------- 2. Section 9.2.5
echo "=== 2. Section 9.2.5 ($(date)) ==="
"$HERE/results_quantitative_verification/reproduce.sh" "$OUT/quantitative_verification" \
  > "$OUT/quantitative_verification.log" 2>&1 || status=1
grep -E "^verification_[a-z_]+: " "$OUT/quantitative_verification.log" || status=1
[ "$(grep -c "identical to the reference results" "$OUT/quantitative_verification.log")" = 3 ] || status=1

# ---------------------------------------------------------------- 3. Appendix H
echo "=== 3. Appendix H ($(date)) ==="
"$HERE/results_bound_verification/reproduce.sh" "$OUT/bound_verification" \
  > "$OUT/bound_verification.log" 2>&1 || status=1
printf 'error bounds (Theorem 3): '
tail -1 "$OUT/bound_verification.log"
grep -q "^identical to the reference results" "$OUT/bound_verification.log" || status=1
"$PYTHON" -I "$HERE/results_bound_verification/unclassified_check.py" \
  > "$OUT/bound_verification/unclassified_check.log" 2> "$OUT/bound_verification/unclassified_check.err" || status=1
printf 'unclassified cells (Theorem 3): '
"$PYTHON" -I - "$HERE/results_bound_verification/unclassified_check.log" \
  "$OUT/bound_verification/unclassified_check.log" <<'EOF' || status=1
import json
import sys


def lines(path):
    out = []
    for line in open(path).read().splitlines():
        if line.startswith("{"):
            row = json.loads(line)
            row.pop("seconds", None)
            out.append(row)
        else:
            out.append(line)
    return out


ref, new = (lines(p) for p in sys.argv[1:3])
ok = ref == new
print("identical to the reference results" if ok else "DIFFERENT from the reference results")
sys.exit(0 if ok else 1)
EOF

echo "=== done in $(( ($(date +%s) - start) / 60 )) min: $([ $status = 0 ] && echo "everything matches the stored results" || echo "SOMETHING DIFFERS, see above") ==="
exit $status
