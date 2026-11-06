#!/bin/bash
# Checks with CAPD that the integration-error budget bounds the true integration error in every
# run of the paper's Section 9 (paper, Appendix I). Each case below repeats a paper run with
# --validated: every reach box uses a CAPD-validated bound on the integration error of the cell
# centre instead of the integrator's tolerance budget, and the solver also builds the box the
# normal build uses and counts every transition whose successors differ (see README.md).
#
# Usage:  CAPD_CONFIG=<CAPD>/build/bin/capd-config ./run_validated.sh [case ...]
#   cases: qv_ra qv_avoid dubins120_std avoid80 dubins20 bicycle8 vdp18 evasion12 di14
#   (default: all). Outputs go to ROOT (default: this folder): runs/<mode>/<run>_validated/
#   and logs/<case>.log. Compare them with the stored paper runs with compare_validated.py.
#
# evasion12 and di14 need 25-30 GB of memory, so macOS swaps to disk. A run is stopped if the
# free disk space drops below MIN_FREE_GB (default 8).
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="${ROOT:-$HERE}"
MIN_FREE_GB="${MIN_FREE_GB:-8}"
WORKERS="${WORKERS:-12}"
: "${CAPD_CONFIG:?set CAPD_CONFIG to <CAPD>/build/bin/capd-config (see README.md)}"

REPO="$HERE"
while [ ! -f "$REPO/cpp/src/main.cpp" ]; do
  [ "$REPO" = / ] && { echo "cpp/src/main.cpp not found above $HERE" >&2; exit 1; }
  REPO="$(dirname "$REPO")"
done
make -C "$REPO/cpp" capd CAPD_CONFIG="$CAPD_CONFIG" || exit 1
BIN="$REPO/cpp/build/dhj_capd"
mkdir -p "$ROOT/runs" "$ROOT/logs"
{
  echo "date: $(date)"
  echo "cases: ${*:-all}"
  echo "git commit: $(git -C "$REPO" rev-parse HEAD) ($(git -C "$REPO" rev-parse --abbrev-ref HEAD))"
  echo "uncommitted solver changes: $(git -C "$REPO" status --short -- cpp | tr '\n' ' ')"
  echo "CAPD: $(dirname "$(dirname "$(dirname "$CAPD_CONFIG")")") at commit $(git -C "$(dirname "$(dirname "$(dirname "$CAPD_CONFIG")")")" rev-parse --short HEAD 2>/dev/null)"
  echo "CAPD flags: $("$CAPD_CONFIG" --cflags | tr ' ' '\n' | grep -v '^-I' | tr '\n' ' ')"
  echo "compiler: $(c++ --version | head -1)"
  echo "binary sha256: $(shasum -a 256 "$BIN" | cut -d' ' -f1)"
} | tee -a "$ROOT/provenance.txt"

common=(--validated --checkpoint-every 0 --workers "$WORKERS" --no-plot --results-root "$ROOT/runs")
run() {
  local name="$1"; shift
  echo "=== START $name $(date) ==="
  "$BIN" "$@" "${common[@]}" > "$ROOT/logs/$name.log" 2>&1 &
  local pid=$!
  caffeinate -ims -w "$pid" &
  while kill -0 "$pid" 2>/dev/null; do
    local free_gb
    free_gb=$(df -g "$ROOT" | awk 'NR==2 {print $4}')
    if [ "$free_gb" -lt "$MIN_FREE_GB" ]; then
      echo "free disk ${free_gb} GB < ${MIN_FREE_GB} GB: stopping $name" | tee -a "$ROOT/logs/$name.log"
      kill "$pid"
    fi
    sleep 30
  done
  wait "$pid"
  echo "=== END $name exit $? $(date) ==="
}

cases=("$@")
[ ${#cases[@]} -eq 0 ] && cases=(qv_ra qv_avoid dubins120_std avoid80 dubins20 bicycle8 vdp18 evasion12 di14)
for c in "${cases[@]}"; do
  case "$c" in
    # Section 9.2.4, reach-avoid (the same run as Fig. 3(b), which differs only by not dumping values)
    qv_ra) run qv_ra --mode ra_discount --gamma 0.96 --dynamics dubins --dt 0.3 --initial-resolution 120 \
             --refinements 0 --delta-lower 1e-8 --delta-upper 1e-8 --dump-csv ;;
    # Section 9.2.4, avoid
    qv_avoid) run qv_avoid --mode avoid_nodiscount --dynamics dubins --dt 0.3 --initial-resolution 120 \
             --refinements 0 --delta-lower 1e-300 --delta-upper 1e-300 --dump-csv ;;
    # Fig. 3(a)
    dubins120_std) run dubins120_std --mode ra_nodiscount --dynamics dubins --dt 0.3 --initial-resolution 120 \
             --refinements 0 --delta-lower 1e-15 --delta-upper 1e-15 ;;
    # Fig. 4
    avoid80) run avoid80 --mode avoid_nodiscount --dynamics dubins --dt 0.3 --initial-resolution 80 \
             --refinements 6 --eps-min 1e-6 --delta-lower 1e-15 --delta-upper 1e-15 ;;
    # Fig. 1
    dubins20) run dubins20 --mode ra_discount --gamma 0.96 --dynamics dubins --dt 0.3 --initial-resolution 20 \
             --refinements 13 --eps-min 1e-6 --delta-lower 1e-8 --delta-upper 1e-8 ;;
    # Fig. 2(d)
    bicycle8) run bicycle8 --mode ra_nodiscount --dynamics bicycle --a-max 0.2 --v-max 1 --dt 0.3 \
             --initial-resolution 5 --refinements 8 --eps-min 1e-6 --delta-lower 1e-15 --delta-upper 1e-15 ;;
    # Fig. 2(a)
    vdp18) run vdp18 --mode ra_nodiscount --dynamics van_der_pol --dt 0.01 --initial-resolution 40 \
             --refinements 18 --eps-min 1e-6 --delta-lower 1e-15 --delta-upper 1e-15 ;;
    # Fig. 2(b)
    evasion12) run evasion12 --mode ra_nodiscount --dynamics evasion --dt 0.3 --initial-resolution 30 \
             --refinements 12 --eps-min 1e-6 --delta-lower 1e-15 --delta-upper 1e-15 ;;
    # Fig. 2(c)
    di14) run di14 --mode ra_nodiscount --dynamics double_integrator --dt 0.3 --initial-resolution 10 \
             --refinements 14 --eps-min 1e-6 --delta-lower 1e-15 --delta-upper 1e-15 ;;
    *) echo "unknown case $c" >&2 ;;
  esac
done
echo "=== ALL DONE $(date) ==="
