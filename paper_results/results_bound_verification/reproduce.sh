#!/bin/bash
# Reproduces the numerical check of the convergence theorem (label thm:convergence):
# Algorithm 1 on a periodic Dubins domain for gamma in {0.3, 0.4, 0.5, 0.6, 0.7}, all
# with gamma * L_f < 1, on uniform grids of N^3 cells for N = 50, 60, ..., 120.
#
#   ./reproduce.sh          writes everything to ../results_bound_verification_rerun
#                           and compares the numbers with summary.json here
#   ./reproduce.sh DIR      the same, into DIR
#   ./reproduce.sh "$PWD"   regenerates the results in this directory (no comparison)
#
# Needs Python 3 with numpy and matplotlib. Recorded environment: Python 3.12.8,
# numpy 2.3.3, matplotlib 3.10.6 (~/miniconda3/bin/python); set PYTHON to use another
# interpreter.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${1:-$HERE/../results_bound_verification_rerun}"
PY="${PYTHON:-$HOME/miniconda3/bin/python}"
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"

"$PY" -I "$HERE/bound_verification.py" --out "$OUT" 2>&1 | tee "$OUT/run.log"

# Compare everything except wall-clock times with the reference numbers in this directory.
if [[ "$OUT" != "$HERE" ]]; then
  "$PY" -I - "$HERE/summary.json" "$OUT/summary.json" <<'EOF'
import json
import sys

ref, new = (json.load(open(p)) for p in sys.argv[1:3])


def numbers(s):
    return [{k: v for k, v in r.items() if k != "seconds"} for r in s["runs"]], s["per_gamma"]


print("identical to the reference results" if numbers(ref) == numbers(new)
      else "DIFFERENT from the reference results")
EOF
fi
