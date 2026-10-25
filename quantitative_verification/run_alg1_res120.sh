#!/bin/bash
# Runs Algorithm 1 once (no refinement) for the Dubins car on a uniform 120^3
# grid and dumps the value tables, so that verify_quantitative.py can check the
# guarantees of the paper:
#   ra     reach-avoid with gamma = 0.96 (Proposition 4 and the classification).
#          The value iteration stops when both sup-norm changes are below 1e-8.
#   avoid  avoid value function with gamma = 1 (Proposition 5). The iterates take
#          values in a finite set, so a nonzero change is far above 1e-300 and the
#          iteration stops only when nothing changes, as in the paper.
#
# Usage: run_alg1_res120.sh [ra|avoid|all] [ROOT]
# Results go to ROOT/<mode>/<run>; ROOT defaults to ../results_quantitative_verification.
#
# The solver is compiled from the current sources into build/dhj_head, so the
# existing build/dhj is left untouched. ROOT/provenance.txt records the commit,
# the compiler and the binary hash (the binary embeds its build time, so the hash
# differs between builds); uncommitted solver changes are saved to
# ROOT/solver_changes.diff.
set -euo pipefail

PART="${1:-all}"
case "$PART" in
  ra|avoid|all) ;;
  *) echo "usage: $0 [ra|avoid|all] [ROOT]" >&2; exit 2 ;;
esac
REPO="$(cd "$(dirname "$0")/.." && pwd)"
ROOT="${2:-$REPO/results_quantitative_verification}"
mkdir -p "$ROOT"

cd "$REPO/cpp"
c++ -std=c++17 -O3 -DNDEBUG -march=native -fno-math-errno -Iinclude -Wall -Wextra \
    src/main.cpp -o build/dhj_head -pthread -lz

rm -f "$ROOT/solver_changes.diff"
git -C "$REPO" diff --quiet -- cpp || git -C "$REPO" diff -- cpp > "$ROOT/solver_changes.diff"
{
  echo "parts: ${PART}"
  echo "git commit: $(git -C "$REPO" rev-parse HEAD)"
  echo "uncommitted solver changes: $(git -C "$REPO" status --short -- cpp | tr '\n' ' ')"
  echo "compiler: $(c++ --version | head -1)"
  echo "binary sha256: $(shasum -a 256 build/dhj_head | cut -d' ' -f1)"
  echo "date: $(date)"
} > "$ROOT/provenance.txt"

run() {
  local mode="$1" gamma="$2" tol="$3"
  echo "=== START ${mode} dubins $(date) ==="
  ./build/dhj_head --mode "$mode" --dynamics dubins \
    --dt 0.3 --tau 0.3 --initial-resolution 120 --epsilon 0.05 \
    --delta-min "$tol" --phase0-tol "$tol" --delta-max "$tol" --gamma "$gamma" \
    --refinements 0 --vi-iterations 100033 \
    --checkpoint-every 0 --keep-checkpoints 0 --dump-csv \
    --plot-dpi 300 --workers 12 --results-root "$ROOT"
  echo "=== END ${mode} dubins $(date) ==="
}

if [[ "$PART" == ra || "$PART" == all ]]; then
  run ra_discount 0.96 1e-8
fi
if [[ "$PART" == avoid || "$PART" == all ]]; then
  run avoid_nodiscount 1 1e-300
fi
echo "=== ALL DONE $(date) ==="
