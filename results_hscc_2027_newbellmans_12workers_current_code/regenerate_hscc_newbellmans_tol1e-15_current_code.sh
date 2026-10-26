#!/bin/bash
# Copy of regenerate_hscc_newbellmans.sh that first rebuilds build/dhj from the
# current solver source (same compiler flags as the Makefile) and records the commit and binary hash in provenance.txt.
# Outputs of earlier runs in this folder are moved to previous_<timestamp>/.
# Other change: the undiscounted runs
# (ra_nodiscount, avoid_nodiscount) use tolerance 1e-15 instead of 1e-8.
# Discounted runs are unchanged. Output: PNG figures only (--plot-dpi 300), no checkpoints.
# Runs only the cases reported in the paper, with the current solver:
#   Fig. 2       dubins20_discount (Algorithm 2, gamma = 0.96, delta* = 1e-8)
#   Fig. 3       vdp18, evasion12, di14, bicycle8 (gamma = 1)
#   Fig. 4(a,b)  dubins120_std (gamma = 1), dubins120_discount (gamma = 0.96)
#   Fig. 5       avoid80 (avoid, gamma = 1)
# Current solver: discounted runs start V_lower at min(l, r) and V_upper at l and
# stop on the residual; undiscounted reach-avoid runs start at min(l, r); avoid
# runs start at l. Each case stops at the last phase the original run completed.
# --eps-min 1e-6 (Algorithm 2's threshold on the cell radius) is far below every
# cell radius these runs refine, so it never stops a refinement; each case ends at
# its --refinements limit. (The stored outputs were made with the earlier
# --epsilon option, which also never stopped a refinement; see README_paper_results.md.)
set -euo pipefail

WORKERS=12
# Results go into the folder that holds this script; the solver sources are in
# ../cpp, so keep the script in a folder directly inside refactored_Implementations/.
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT/../cpp"
mkdir -p "$ROOT"

# Rebuild the solver from the current source, so every run uses the latest code.
c++ -std=c++17 -O3 -DNDEBUG -march=native -fno-math-errno -Iinclude -Wall -Wextra \
    src/main.cpp -o build/dhj -pthread -lz
rm -f "$ROOT/solver_changes.diff"
git diff --quiet -- . || git diff -- . > "$ROOT/solver_changes.diff"
{
  echo "date: $(date)"
  echo "git commit: $(git rev-parse HEAD)"
  echo "uncommitted solver changes: $(git status --short -- . | tr '\n' ' ')"
  echo "compiler: $(c++ --version | head -1)"
  echo "binary sha256: $(shasum -a 256 build/dhj | cut -d' ' -f1)"
} > "$ROOT/provenance.txt"
cat "$ROOT/provenance.txt"

# Keep outputs of earlier runs instead of overwriting them.
dest="$ROOT/previous_$(date +%Y%m%d_%H%M%S)"
for m in ra_nodiscount ra_discount avoid_nodiscount; do
  if [ -d "$ROOT/$m" ]; then
    mkdir -p "$dest"
    mv "$ROOT/$m" "$dest/"
    echo "Moved earlier $m outputs to $dest"
  fi
done

run() {
  local name="$1"
  shift
  echo "=== START ${name} $(date) ==="
  caffeinate -ims ./build/dhj "$@" \
    --checkpoint-every 0 --keep-checkpoints 0 \
    --plot-dpi 300 --workers "$WORKERS"  \
    --results-root "$ROOT"
  echo "=== END ${name} $(date) ==="
}

run dubins120_std \
  --mode ra_nodiscount --dynamics dubins \
  --dt 0.3 --tau 0.3 --initial-resolution 120 --eps-min 1e-6 \
  --delta-min 1e-15 --delta-max 1e-15 --gamma 0.96 \
  --refinements 0 --vi-iterations 100033

run dubins120_discount \
  --mode ra_discount --dynamics dubins \
  --dt 0.3 --tau 0.3 --initial-resolution 120 --eps-min 1e-6 \
  --delta-min 1e-8 --delta-max 1e-8 --gamma 0.96 \
  --refinements 0 --vi-iterations 100033

run avoid80 \
  --mode avoid_nodiscount --dynamics dubins \
  --dt 0.3 --tau 0.3 --initial-resolution 80 --eps-min 1e-6 \
  --delta-min 1e-15 --delta-max 1e-15 --gamma 0.9 \
  --refinements 6 --vi-iterations 20000

run dubins20_discount \
  --mode ra_discount --dynamics dubins \
  --dt 0.3 --tau 0.3 --initial-resolution 20 --eps-min 1e-6 \
  --delta-min 1e-8 --delta-max 1e-8 --gamma 0.96 \
  --refinements 13 --vi-iterations 100033

run bicycle8 \
  --mode ra_nodiscount --dynamics bicycle --a-max 0.2 --v-max 1 \
  --dt 0.3 --tau 0.3 --initial-resolution 5 --eps-min 1e-6 \
  --delta-min 1e-15 --delta-max 1e-15 --gamma 0.9 \
  --refinements 8 --vi-iterations 20000

run vdp18 \
  --mode ra_nodiscount --dynamics van_der_pol \
  --dt 0.01 --tau 0.01 --initial-resolution 40 --eps-min 1e-6 \
  --delta-min 1e-15 --delta-max 1e-15 --gamma 0.96 \
  --refinements 18 --vi-iterations 20000

run di14 \
  --mode ra_nodiscount --dynamics double_integrator \
  --dt 0.3 --tau 0.3 --initial-resolution 10 --eps-min 1e-6 \
  --delta-min 1e-15 --delta-max 1e-15 --gamma 0.9 \
  --refinements 14 --vi-iterations 20000

run evasion12 \
  --mode ra_nodiscount --dynamics evasion \
  --dt 0.3 --tau 0.3 --initial-resolution 30 --eps-min 1e-6 \
  --delta-min 1e-15 --delta-max 1e-15 --gamma 0.96 \
  --refinements 12 --vi-iterations 1000033

echo "=== ALL DONE $(date) ==="
