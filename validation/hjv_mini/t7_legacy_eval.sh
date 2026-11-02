#!/bin/bash
# T7, legacy side: one likelihood, gradient and Hessian evaluation at the truth via stage 0 (TestCode), serial.
# Usage: t7_legacy_eval.sh <N> <repeats>. Needs results/legacy_s0_<N>/init_par.txt (same truth for every N).
cd "$(dirname "$0")"; N="$1"; reps="${2:-5}"
d="results/legacy_s0_$N"; mkdir -p "$d"; [ -f "$d/init_par.txt" ] || cp results/legacy_s0_2000/init_par.txt "$d/"
[ -x "prob_mini_s0_$N" ] || OUT="prob_mini_s0_$N" legacy_build/build_prob.sh prob_hjv_mini.cc serial -DSTAGE=0 -DNOBS="$N" >/dev/null
: > "$d/t7_times.txt"
for r in $(seq 1 "$reps"); do
  ./prob_mini_s0_$N > "$d/t7_run_$r.log" 2>&1
  awk '/Calculating the Likelihood took/{l=$5} /Calculating the Gradient/{g=$5} /Calculating the Hessian/{h=$5} END{print l, g, h}' "$d/t7_run_$r.log" >> "$d/t7_times.txt"
done
echo "legacy N=$N per-evaluation seconds (loglik, loglik+grad, loglik+grad+hess), $reps repeats:"; cat "$d/t7_times.txt"
