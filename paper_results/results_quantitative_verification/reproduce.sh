#!/bin/bash
# Produces or replicates every result in this directory: two runs of Algorithm 1
# for the Dubins car on a uniform 120^3 grid (reach-avoid with gamma = 0.96, and
# the avoid value function with gamma = 1) and their verification with
# verify_quantitative.py, including the tables and figures, the check of the
# violation claim of Proposition 4 with verify_violation.py, and the check of the
# failure-reaching time of Proposition 6(2) with verify_avoid_failure_time.py.
#
#   ./reproduce.sh          writes everything to ../results_quantitative_verification_rerun
#                           and compares the numbers with verification_*/summary.json,
#                           verification_ra_violation/report.json and
#                           verification_avoid_failure_time/report.json here
#   ./reproduce.sh DIR      the same, into DIR
#   ./reproduce.sh "$PWD"   regenerates the results in this directory (no comparison)
#
# Recorded environment: see provenance.txt for the commit and the compiler;
# Python 3.12.8 with numpy 2.3.3, scipy 1.16.2 and matplotlib 3.10.6
# (~/miniconda3/bin/python). Set PYTHON to use another interpreter. Figure text
# is typeset with LaTeX when latex and dvipng are installed.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
OUT="${1:-$REPO/results_quantitative_verification_rerun}"
PY="${PYTHON:-$HOME/miniconda3/bin/python}"
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"

# 1. Solver: Algorithm 1 only (--refinements 0), value tables dumped as CSV.
"$HERE/run_alg1_res120.sh" all "$OUT" 2>&1 | tee "$OUT/solver.log"

# The solver names each run folder after its options (dubins_res120_..._dt0.3); use the
# newest one of each mode, i.e. the run just made.
run_dir() {
  local d
  d="$(ls -td "$OUT/$1"/dubins_res120_*/ 2>/dev/null | head -1)"
  [ -n "$d" ] || { echo "no solver run found in $OUT/$1" >&2; exit 1; }
  echo "${d%/}"
}
RA_RUN="$(run_dir ra_discount)"
AVOID_RUN="$(run_dir avoid_nodiscount)"

# 2. Verification: seed 0, one random state per cell, 14 example paths.
"$PY" -I "$HERE/verify_quantitative.py" "$RA_RUN" \
  --out "$OUT/verification_ra" --seed 0 --points-per-cell 1 --paths 14 \
  2>&1 | tee "$OUT/verification_ra.log"
"$PY" -I "$HERE/verify_quantitative.py" "$AVOID_RUN" \
  --out "$OUT/verification_avoid" --seed 0 --points-per-cell 1 --paths 14 --horizon 1000 \
  2>&1 | tee "$OUT/verification_avoid.log"

# Violation claim of Proposition 4: from one random state and the 8 vertices of every
# cell with V_upper < 0, eight policies and an exhaustive search over all action sequences.
"$PY" -I "$HERE/verify_violation.py" "$RA_RUN" \
  --out "$OUT/verification_ra_violation" --seed 0 --corners --exhaustive -1 --random-seeds 3 \
  2>&1 | tee "$OUT/verification_ra_violation.log"

# Failure-reaching time of Proposition 6(2): the solver's k(s) (first iteration with
# V_upper^(k)(s) <= 0), its recomputation and a certificate on the verifier's abstraction,
# and an exhaustive search over all action sequences from one random state and the
# 8 vertices of every cell with V_upper <= 0.
"$PY" -I "$HERE/verify_avoid_failure_time.py" "$AVOID_RUN" \
  --out "$OUT/verification_avoid_failure_time" --seed 0 \
  2>&1 | tee "$OUT/verification_avoid_failure_time.log"

# 3. Compare everything except timings with the reference numbers in this directory.
if [[ "$OUT" != "$HERE" ]]; then
  for v in verification_ra verification_avoid; do
    printf '%s: ' "$v"
    "$PY" -I - "$HERE/$v/summary.json" "$OUT/$v/summary.json" <<'EOF'
import json
import sys

def numbers(rep):
    """Everything except paths and wall-clock times."""
    rep = {k: v for k, v in rep.items() if k not in ("run_dir", "verifier_seconds")}
    rep["solver"] = {k: v for k, v in rep["solver"].items() if not k.startswith("time_")}
    return rep


ref, new = ([numbers(r) for r in json.load(open(p))] for p in sys.argv[1:3])
diffs = []
for a, b in zip(ref, new):
    for key in a:
        if a[key] != b.get(key):
            diffs.append(f"{a['system']}.{key}: {a[key]} -> {b.get(key)}")
if len(ref) != len(new):
    diffs.append(f"{len(ref)} systems in the reference, {len(new)} in the new results")
print("identical to the reference results" if not diffs else "\n".join(diffs))
EOF
  done
  for v in verification_ra_violation verification_avoid_failure_time; do
    printf '%s: ' "$v"
    "$PY" -I - "$HERE/$v/report.json" "$OUT/$v/report.json" <<'EOF'
import json
import sys

def numbers(o):
    """Everything except paths and wall-clock times."""
    if isinstance(o, dict):
        return {k: numbers(v) for k, v in o.items() if k not in ("run_dir", "seconds")}
    return o


ref, new = (numbers(json.load(open(p))) for p in sys.argv[1:3])
diffs = [f"{k}: {ref[k]} -> {new.get(k)}" for k in ref if ref[k] != new.get(k)]
print("identical to the reference results" if not diffs else "\n".join(diffs))
EOF
  done
fi
