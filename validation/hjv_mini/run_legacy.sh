#!/bin/bash
# Usage: run_legacy.sh <stage 1|2> <N> [np]. Builds prob_mini_s<stage> for N and runs it with mpirun -np <np>.
set -u
stage="$1"; N="$2"; np="${3:-8}"
cd "$(dirname "$0")"
out="results/legacy_s${stage}_${N}"; mkdir -p "$out"
OUT="prob_mini_s${stage}_${N}" legacy_build/build_prob.sh prob_hjv_mini.cc mpi -DSTAGE="$stage" -DNOBS="$N" || exit 1
start=$(date +%s.%N)
mpirun -np "$np" "./prob_mini_s${stage}_${N}_mpi" > "$out/run.log" 2>&1
rc=$?
end=$(date +%s.%N)
printf "legacy stage %s N=%s: exit=%d secs=%.1f np=%s\n" "$stage" "$N" "$rc" "$(echo "$end - $start" | bc)" "$np" | tee "$out/timing.txt"
grep -E "Final Likelihood|Number of Iterations|EXIT|FAILED|ABORT|Error" "$out/run.log" | head -8
