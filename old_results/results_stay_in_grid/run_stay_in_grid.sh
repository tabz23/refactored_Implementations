#!/bin/bash
# Reruns the eight cases of Figs. 2-5 with --restrict-reach-set-to-grid: the part of each
# reach box outside the state box X is ignored and no sink state is added, as if f mapped X
# into X. Our systems can leave X, so these results are NOT sound certificates. They show what
# the paper's figures would be under that assumption. The paper's results, with the sink,
# are in ../paper_results/results_hscc_2027_newbellmans_12workers_current_code/, and
# compare_with_paper.py compares the two phase by phase.
#
# By default every case finishes within a few minutes: vdp18, di14 and evasion12 stop after
# 16, 10 and 9 refinements (as in reproduce_all.sh --quick), and any case still running after
# TIME_LIMIT seconds (default 300) is stopped; its completed phases stay in its run.log.
# FULL=1 runs the paper's refinements with no time limit.
#
# Outputs of earlier runs in this folder move to previous_<timestamp>/.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
CPP="${CPP_DIR:-}"
if [ -z "$CPP" ]; then
  d="$ROOT"
  while [ ! -f "$d/cpp/src/main.cpp" ]; do
    [ "$d" = / ] && { echo "cpp/src/main.cpp not found above $ROOT; set CPP_DIR" >&2; exit 1; }
    d="$(dirname "$d")"
  done
  CPP="$d/cpp"
fi
FULL="${FULL:-0}"
TIME_LIMIT="${TIME_LIMIT:-300}"
WORKERS=12

cd "$CPP"
c++ -std=c++17 -O3 -DNDEBUG -march=native -fno-math-errno -Iinclude -Wall -Wextra \
    src/main.cpp -o build/dhj -pthread -lz
{
  echo "date: $(date)"
  echo "full: $FULL, time limit per case: $([ "$FULL" = 1 ] && echo none || echo "${TIME_LIMIT} s")"
  echo "git commit: $(git rev-parse HEAD)"
  echo "uncommitted solver changes: $(git status --short -- . | tr '\n' ' ')"
  echo "binary sha256: $(shasum -a 256 build/dhj | cut -d' ' -f1)"
} | tee "$ROOT/provenance.txt"

dest="$ROOT/previous_$(date +%Y%m%d_%H%M%S)"
for m in ra_nodiscount ra_discount avoid_nodiscount; do
  if [ -d "$ROOT/$m" ]; then mkdir -p "$dest"; mv "$ROOT/$m" "$dest/"; echo "Moved earlier $m outputs to $dest"; fi
done

quick_refinements() {   # the paper uses 18, 14 and 12
  case "$1" in vdp18) echo 16 ;; di14) echo 10 ;; evasion12) echo 9 ;; esac
}

run() {
  local name="$1" q
  shift
  q="$(quick_refinements "$name")"
  local limit=()
  if [ "$FULL" != 1 ]; then
    [ -z "$q" ] || set -- "$@" --refinements "$q"   # the last --refinements counts
    limit=(perl -e 'alarm shift; exec @ARGV' "$TIME_LIMIT")
  fi
  echo "=== START ${name} $(date) ==="
  # ${limit[@]+...}: bash 3.2 rejects an empty "${limit[@]}" under set -u
  caffeinate -ims ${limit[@]+"${limit[@]}"} ./build/dhj "$@" --restrict-reach-set-to-grid \
    --checkpoint-every 0 --plot-dpi 300 --workers "$WORKERS" --results-root "$ROOT" \
    || echo "=== ${name} stopped (exit $?), e.g. by the time limit ==="
  echo "=== END ${name} $(date) ==="
}

# The same cases and flags as the paper (paper_results/.../regenerate_hscc_newbellmans_tol1e-15_current_code.sh).
run dubins120_std      --mode ra_nodiscount --dynamics dubins --dt 0.3 --initial-resolution 120 --refinements 0 --delta-lower 1e-15 --delta-upper 1e-15
run dubins120_discount --mode ra_discount --gamma 0.96 --dynamics dubins --dt 0.3 --initial-resolution 120 --refinements 0 --delta-lower 1e-8 --delta-upper 1e-8
run avoid80            --mode avoid_nodiscount --dynamics dubins --dt 0.3 --initial-resolution 80 --refinements 6 --eps-min 1e-6 --delta-lower 1e-15 --delta-upper 1e-15
run dubins20_discount  --mode ra_discount --gamma 0.96 --dynamics dubins --dt 0.3 --initial-resolution 20 --refinements 13 --eps-min 1e-6 --delta-lower 1e-8 --delta-upper 1e-8
run bicycle8           --mode ra_nodiscount --dynamics bicycle --a-max 0.2 --v-max 1 --dt 0.3 --initial-resolution 5 --refinements 8 --eps-min 1e-6 --delta-lower 1e-15 --delta-upper 1e-15
run vdp18              --mode ra_nodiscount --dynamics van_der_pol --dt 0.01 --initial-resolution 40 --refinements 18 --eps-min 1e-6 --delta-lower 1e-15 --delta-upper 1e-15
run di14               --mode ra_nodiscount --dynamics double_integrator --dt 0.3 --initial-resolution 10 --refinements 14 --eps-min 1e-6 --delta-lower 1e-15 --delta-upper 1e-15
run evasion12          --mode ra_nodiscount --dynamics evasion --dt 0.3 --initial-resolution 30 --refinements 12 --eps-min 1e-6 --delta-lower 1e-15 --delta-upper 1e-15

echo "=== ALL DONE $(date) ==="
