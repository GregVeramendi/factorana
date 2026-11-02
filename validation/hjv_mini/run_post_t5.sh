#!/bin/bash
# After the first T5 pass: T4 (legacy scores serially), T1a (measurement-only derivatives), then the T5 "L" leg again.
cd "$(dirname "$0")"; N=2000
{
echo "== T4 legacy scores ML (serial) $(date)"; s=$(date +%s); ./prob_mini_s5_$N > results/legacy_s5_$N/run.log 2>&1; echo "exit=$? secs=$(( $(date +%s) - s ))"
echo "== T4 legacy scores with prior (serial) $(date)"; s=$(date +%s); ./prob_mini_s5p_$N > results/legacy_s5p_$N/run.log 2>&1; echo "exit=$? secs=$(( $(date +%s) - s ))"
wc -l results/legacy_s5_$N/factor_predictions.txt results/legacy_s5p_$N/factor_predictions.txt
echo "== T4 compare"; Rscript t4_compare.R $N 2>&1 | grep -vE "^Test passed|Warning|max;|^[0-9]+: In"
echo "== T1a legacy stage 0, measurement only $(date)"; ./prob_mini_s0m_$N > results/legacy_s0m_$N/run.log 2>&1; grep -E "Likelihood at initial|took|Wrote" results/legacy_s0m_$N/run.log
echo "== T1a R eval, measurement system"; Rscript t1_r_eval.R $N meas 2>&1 | grep -E "^R T1|free parameters"
echo "== T1a compare"; Rscript t1_compare.R $N m 2>&1 | grep -vE "^Test passed"
echo "== T5 leg L: legacy adaptive with legacy ML scores $(date)"
d="results/legacy_s2aL_$N"; mkdir -p "$d"; cp results/legacy_s1_$N/meas_par.txt "$d/"; cp results/legacy_s5_$N/factor_predictions.txt "$d/"
s=$(date +%s); mpirun -np 8 ./prob_mini_s2aL_${N}_mpi > "$d/run.log" 2>&1; echo "exit=$? secs=$(( $(date +%s) - s ))"; grep -E "Number of Iterations|EXIT" "$d/run.log" | head -4
echo "== T5 leg L: R adaptive with legacy ML scores"; Rscript t5_r_adapt.R $N L 8 2>&1 | grep -vE "^Test passed|^$"
echo "== compare"; Rscript compare_stage.R s3aL s2aL $N 1 2>&1 | grep -vE "^Test passed" | head -8
Rscript t5_vs_full.R $N L 2>&1 | grep -vE "^Test passed"; Rscript t5_points.R $N L 2>&1 | grep -vE "^Test passed"
echo "== done $(date)"
} > results/post_t5.log 2>&1
touch results/post_t5.done
