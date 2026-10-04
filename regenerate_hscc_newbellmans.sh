#!/bin/bash
# Regenerates the results_hscc_2027 cases with the current solver.
# Discounted runs start at l and use conservative stopping.
# Undiscounted runs start at min(l, r) and do not.
# Each case stops at the last phase the original run completed.
# WORKERS is 10. The output directory name says 12workers.
set -euo pipefail

WORKERS=12
ROOT="/Users/i.k.tabbara/Documents/python directory/HSCC_2027/refactored_Implementations/results_hscc_2027_newbellmans_12workers"
cd "/Users/i.k.tabbara/Documents/python directory/HSCC_2027/refactored_Implementations/cpp"
mkdir -p "$ROOT"

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
  --dt 0.3 --tau 0.3 --initial-resolution 120 --epsilon 0.05 \
  --delta-min 1e-8 --delta-max 1e-8 --gamma 0.96 \
  --refinements 0 --vi-iterations 100033

run dubins120_discount \
  --mode ra_discount --dynamics dubins \
  --dt 0.3 --tau 0.3 --initial-resolution 120 --epsilon 0.05 \
  --delta-min 1e-8 --delta-max 1e-8 --gamma 0.96 \
  --refinements 0 --vi-iterations 100033

run dubins120_discount_1e-3 \
  --mode ra_discount --dynamics dubins \
  --dt 0.3 --tau 0.3 --initial-resolution 120 --epsilon 0.05 \
  --delta-min 1e-3 --delta-max 1e-3 --gamma 0.96 \
  --refinements 0 --vi-iterations 100033

run avoid80 \
  --mode avoid_nodiscount --dynamics dubins \
  --dt 0.3 --tau 0.3 --initial-resolution 80 --epsilon 0.05 \
  --delta-min 1e-8 --delta-max 1e-8 --gamma 0.9 \
  --refinements 6 --vi-iterations 20000

run dubins20_std \
  --mode ra_nodiscount --dynamics dubins \
  --dt 0.3 --tau 0.3 --initial-resolution 20 --epsilon 0.05 \
  --delta-min 1e-8 --delta-max 1e-8 --gamma 0.9 \
  --refinements 13 --vi-iterations 100033

run dubins20_discount \
  --mode ra_discount --dynamics dubins \
  --dt 0.3 --tau 0.3 --initial-resolution 20 --epsilon 0.05 \
  --delta-min 1e-8 --delta-max 1e-8 --gamma 0.96 \
  --refinements 13 --vi-iterations 100033

run bicycle8 \
  --mode ra_nodiscount --dynamics bicycle --a-max 0.2 --v-max 1 \
  --dt 0.3 --tau 0.3 --initial-resolution 5 --epsilon 0.1 \
  --delta-min 1e-8 --delta-max 1e-8 --gamma 0.9 \
  --refinements 8 --vi-iterations 20000

run vdp18 \
  --mode ra_nodiscount --dynamics van_der_pol \
  --dt 0.01 --tau 0.01 --initial-resolution 40 --epsilon 0.001 \
  --delta-min 1e-8 --delta-max 1e-8 --gamma 0.96 \
  --refinements 18 --vi-iterations 20000

run di14 \
  --mode ra_nodiscount --dynamics double_integrator \
  --dt 0.3 --tau 0.3 --initial-resolution 10 --epsilon 0.05 \
  --delta-min 1e-8 --delta-max 1e-8 --gamma 0.9 \
  --refinements 14 --vi-iterations 20000

run evasion12 \
  --mode ra_nodiscount --dynamics evasion \
  --dt 0.3 --tau 0.3 --initial-resolution 30 --epsilon 0.05 \
  --delta-min 1e-8 --delta-max 1e-8 --gamma 0.96 \
  --refinements 12 --vi-iterations 1000033

echo "=== ALL DONE $(date) ==="
