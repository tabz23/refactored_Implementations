#!/bin/bash
# Produces or replicates every result in this directory: two runs of Algorithm 1
# for the Dubins car on a uniform 120^3 grid (reach-avoid with gamma = 0.96, and
# the avoid value function with gamma = 1) and their verification with
# verify_quantitative.py, including the tables and figures.
#
#   ./reproduce.sh          writes everything to ../results_quantitative_verification_rerun
#                           and compares the numbers with verification_*/summary.json here
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

# 2. Verification: seed 0, one random state per cell, 14 example paths.
"$PY" -I "$HERE/verify_quantitative.py" \
  "$OUT/ra_discount/dubins_res120_epsmin1e-06_dmin1e-08_dmax1e-08_gamma0.96_dt0.3_tau0.3_std" \
  --out "$OUT/verification_ra" --seed 0 --points-per-cell 1 --paths 14 \
  2>&1 | tee "$OUT/verification_ra.log"
"$PY" -I "$HERE/verify_quantitative.py" \
  "$OUT/avoid_nodiscount/dubins_res120_epsmin1e-06_dmin1e-300_dmax1e-300_gamma1.0_dt0.3_tau0.3_std" \
  --out "$OUT/verification_avoid" --seed 0 --points-per-cell 1 --paths 14 --horizon 1000 \
  2>&1 | tee "$OUT/verification_avoid.log"

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
fi
