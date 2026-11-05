#!/bin/bash
# Regenerates Figs. 2-5 of the paper and their runtimes: rebuilds build/dhj from the current
# solver source (the flags of cpp/Makefile), records the commit and binary hash in
# provenance.txt, and runs the eight cases below. Results go into the folder that holds this
# script; outputs of earlier runs there are moved to previous_<timestamp>/.
#
#   Fig. 2       dubins20_discount       Algorithm 2, discounted reach-avoid
#   Fig. 3       vdp18, evasion12, di14, bicycle8   Algorithm 2, reach-avoid with gamma = 1
#   Fig. 4(a,b)  dubins120_std, dubins120_discount  Algorithm 1 only (--refinements 0)
#   Fig. 5       avoid80                 Algorithm 2, avoid with gamma = 1
#
# Parameters of each case (see `cpp/build/dhj --help`):
#   --mode                ra_discount (discounted, --gamma), ra_nodiscount or avoid_nodiscount
#                         (both with gamma = 1, so they take no --gamma)
#   --dt                  the sampling time t_s
#   --initial-resolution  cells per dimension of the initial uniform grid
#   --refinements         refinement iterations of Algorithm 2 (0: Algorithm 1 only)
#   --eps-min             Algorithm 2's threshold epsilon_min on the cell radius. 1e-6 is far
#                         below every cell these runs refine, so each case ends after its
#                         --refinements iterations. Runs without refinement do not use it.
#   --delta-lower/upper   stopping tolerances of value iteration for the lower and upper bound
#                         (the paper's underline-delta* and overline-delta*): 1e-8 in the
#                         discounted runs, as in Section 9.2.5; 1e-15 in the gamma = 1 runs,
#                         which stop once nothing changes beyond rounding.
# All runs write PNG figures (--plot-dpi 300), no checkpoints, and use 12 worker threads.
#
# QUICK=1 runs a smoke test (reproduce_all.sh --quick): the three slowest cases stop after
# fewer refinement iterations than in the paper, given by quick_refinements below, so every
# case takes at most a few minutes. Their phases are the first phases of the full runs.
#
# CASES="di14 evasion12" runs only the listed cases and keeps the outputs already in this
# folder (used by reproduce_all.sh --updated --cases to finish an interrupted run). The
# provenance of such a rerun is appended to provenance.txt.
#
# Memory: the last iterations of di14 (14.8 million cells) and evasion12 (24.8 million cells)
# need far more than 18 GB, so macOS swaps to disk. Close memory-heavy apps and keep about
# 40 GB of disk free, or the system kills the solver when the disk fills.
set -euo pipefail

WORKERS=12
QUICK="${QUICK:-0}"
CASES="${CASES:-}"
quick_refinements() {   # smoke test only; the cases below give the paper's values
  case "$1" in
    vdp18) echo 16 ;;      # instead of 18
    di14) echo 10 ;;       # instead of 14
    evasion12) echo 9 ;;   # instead of 12
  esac
}
# The solver sources are the nearest cpp/ above this folder (refactored_Implementations/cpp);
# set CPP_DIR to use another copy.
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
cd "$CPP"
mkdir -p "$ROOT"

# Rebuild the solver from the current source, so every run uses the latest code.
c++ -std=c++17 -O3 -DNDEBUG -march=native -fno-math-errno -Iinclude -Wall -Wextra \
    src/main.cpp -o build/dhj -pthread -lz
suffix=""
[ -z "$CASES" ] || suffix="_cases"
rm -f "$ROOT/solver_changes${suffix}.diff"
git diff --quiet -- . || git diff -- . > "$ROOT/solver_changes${suffix}.diff"
provenance() {
  echo "date: $(date)"
  echo "quick: $QUICK"
  [ -z "$CASES" ] || echo "cases: $CASES"
  echo "git commit: $(git rev-parse HEAD)"
  echo "uncommitted solver changes: $(git status --short -- . | tr '\n' ' ')"
  echo "compiler: $(c++ --version | head -1)"
  echo "binary sha256: $(shasum -a 256 build/dhj | cut -d' ' -f1)"
}
if [ -z "$CASES" ]; then
  provenance > "$ROOT/provenance.txt"
  cat "$ROOT/provenance.txt"
else
  { echo "--- rerun of: $CASES"; provenance; } | tee -a "$ROOT/provenance.txt"
fi

# Keep outputs of earlier runs instead of overwriting them (not when rerunning some cases).
if [ -z "$CASES" ]; then
  dest="$ROOT/previous_$(date +%Y%m%d_%H%M%S)"
  for m in ra_nodiscount ra_discount avoid_nodiscount; do
    if [ -d "$ROOT/$m" ]; then
      mkdir -p "$dest"
      mv "$ROOT/$m" "$dest/"
      echo "Moved earlier $m outputs to $dest"
    fi
  done
fi

run() {
  local name="$1" quick
  shift
  if [ -n "$CASES" ] && [[ " $CASES " != *" $name "* ]]; then
    echo "=== SKIP ${name} (not in CASES) ==="
    return
  fi
  quick="$(quick_refinements "$name")"
  # In the smoke test, a second --refinements overrides the case's own (the last one counts).
  if [ "$QUICK" = 1 ] && [ -n "$quick" ]; then set -- "$@" --refinements "$quick"; fi
  echo "=== START ${name} $(date) ==="
  caffeinate -ims ./build/dhj "$@" \
    --checkpoint-every 0 --plot-dpi 300 --workers "$WORKERS" --results-root "$ROOT"
  echo "=== END ${name} $(date) ==="
}

# Fig. 4(a): Algorithm 1, reach-avoid with gamma = 1 (Proposition 3), 120^3 grid.
run dubins120_std \
  --mode ra_nodiscount --dynamics dubins --dt 0.3 \
  --initial-resolution 120 --refinements 0 \
  --delta-lower 1e-15 --delta-upper 1e-15

# Fig. 4(b): Algorithm 1, discounted reach-avoid, 120^3 grid.
run dubins120_discount \
  --mode ra_discount --gamma 0.96 --dynamics dubins --dt 0.3 \
  --initial-resolution 120 --refinements 0 \
  --delta-lower 1e-8 --delta-upper 1e-8

# Fig. 5: Algorithm 2, avoid with gamma = 1 (Section 8), 80^3 initial grid.
run avoid80 \
  --mode avoid_nodiscount --dynamics dubins --dt 0.3 \
  --initial-resolution 80 --refinements 6 --eps-min 1e-6 \
  --delta-lower 1e-15 --delta-upper 1e-15

# Fig. 2: Algorithm 2, discounted reach-avoid, 20^3 initial grid.
run dubins20_discount \
  --mode ra_discount --gamma 0.96 --dynamics dubins --dt 0.3 \
  --initial-resolution 20 --refinements 13 --eps-min 1e-6 \
  --delta-lower 1e-8 --delta-upper 1e-8

# Fig. 3(d): Algorithm 2, reach-avoid with gamma = 1, kinematic bicycle, 5^5 initial grid.
run bicycle8 \
  --mode ra_nodiscount --dynamics bicycle --a-max 0.2 --v-max 1 --dt 0.3 \
  --initial-resolution 5 --refinements 8 --eps-min 1e-6 \
  --delta-lower 1e-15 --delta-upper 1e-15

# Fig. 3(a): Algorithm 2, reach-avoid with gamma = 1, van der Pol, 40^2 initial grid.
run vdp18 \
  --mode ra_nodiscount --dynamics van_der_pol --dt 0.01 \
  --initial-resolution 40 --refinements 18 --eps-min 1e-6 \
  --delta-lower 1e-15 --delta-upper 1e-15

# Fig. 3(c): Algorithm 2, reach-avoid with gamma = 1, planar double integrator, 10^4 initial grid.
run di14 \
  --mode ra_nodiscount --dynamics double_integrator --dt 0.3 \
  --initial-resolution 10 --refinements 14 --eps-min 1e-6 \
  --delta-lower 1e-15 --delta-upper 1e-15

# Fig. 3(b): Algorithm 2, reach-avoid with gamma = 1, 3D evasion, 30^3 initial grid.
run evasion12 \
  --mode ra_nodiscount --dynamics evasion --dt 0.3 \
  --initial-resolution 30 --refinements 12 --eps-min 1e-6 \
  --delta-lower 1e-15 --delta-upper 1e-15

echo "=== ALL DONE $(date) ==="
