#!/usr/bin/env bash
# Benchmark Phase 3: site-pattern likelihood parallelisation
# Runs bench-site-likelihoods.Rev at OMP_NUM_THREADS=1,2,4,8 and reports speedup.
#
# Usage (from benchmarks/):
#   ./bench-site-likelihoods.sh [data.nex] [alphabet] [ngens]
#
# Defaults: primates cytb, DNA, 3000 generations

set -euo pipefail

DATA="${1:-../tests/data/primates_and_galeopterus_cytb.nex}"
ALPHA="${2:-DNA}"
NGENS="${3:-3000}"
RB="${RB:-../projects/cmake/rb}"

if [[ ! -x "$RB" ]]; then
    echo "rb not found at $RB — set RB=/path/to/rb or run from benchmarks/ dir"
    exit 1
fi

echo "=== Phase 3 site-likelihood benchmark ==="
echo "data=$DATA  alphabet=$ALPHA  ngens=$NGENS"
echo ""

run_bench() {
    local threads=$1
    local t0 t1 elapsed
    t0=$(date +%s%N)
    OMP_NUM_THREADS=$threads "$RB" bench-site-likelihoods.Rev "$DATA" "$ALPHA" "$NGENS" \
        > /dev/null 2>&1
    t1=$(date +%s%N)
    echo $(( (t1 - t0) / 1000000 ))   # milliseconds
}

echo -n "threads=1 (baseline) ... "
T1=$(run_bench 1)
echo "${T1} ms"

for T in 2 4 8; do
    echo -n "threads=$T ... "
    TN=$(run_bench $T)
    SPEEDUP=$(awk "BEGIN {printf \"%.2f\", $T1/$TN}")
    echo "${TN} ms  (speedup vs 1 thread: ${SPEEDUP}x)"
done

echo ""
echo "Done."
