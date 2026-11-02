#!/bin/bash
# T7 speed. Per-evaluation timing (single core) at N = 2000 and 5000 for both codes; likelihood+gradient only at 20000.
# Time to convergence with 8 cores, stage 1 and the education stage, N = 5000. Scaling 1/4/8 cores, stage 1, N = 2000.
cd "$(dirname "$0")"
{
echo "== per-evaluation timing $(date)"
./t7_legacy_eval.sh 2000 3; ./t7_legacy_eval.sh 5000 3
Rscript t7_r_timing.R 2000 3 2>&1 | grep -E "^R N="; Rscript t7_r_timing.R 5000 1 2>&1 | grep -E "^R N="
echo "== convergence timing, N = 5000, 8 cores $(date)"
./run_legacy.sh 1 5000 8 | grep -E "^legacy|Iterations"
mkdir -p results/legacy_s2_5000 && cp results/legacy_s1_5000/meas_par.txt results/legacy_s2_5000/ && ./run_legacy.sh 2 5000 8 | grep -E "^legacy|Iterations"
Rscript run_r.R 1 5000 8 2>&1 | grep -E "^R stage"; Rscript run_r.R 3 5000 8 2>&1 | grep -E "^R stage"
Rscript compare_stage.R 1 1 5000 2>&1 | grep -E "^(matched|estimates|SEs|loglik|pass)"; Rscript compare_stage.R 3 2 5000 1 2>&1 | grep -E "^(matched|estimates|SEs|loglik|pass)"
echo "== scaling, stage 1, N = 2000 $(date)"
for np in 1 4; do d=results/legacy_s1_2000_np$np; mkdir -p $d; s=$(date +%s); mpirun -np $np ./prob_mini_s1_2000_mpi > $d/run.log 2>&1; echo "legacy stage 1 N=2000 np=$np secs=$(( $(date +%s) - s ))"; done
for nc in 1 4; do Rscript run_r.R 1 2000 $nc 2>&1 | grep -E "^R stage"; done
echo "== memory, stage 1, N = 5000, 8 cores (peak RSS of the launcher and ranks)"
/usr/bin/time -l mpirun -np 8 ./prob_mini_s1_5000_mpi > /dev/null 2> results/t7_mem_legacy.txt; grep -E "maximum resident|real" results/t7_mem_legacy.txt
/usr/bin/time -l Rscript run_r.R 1 5000 8 > /dev/null 2> results/t7_mem_r.txt; grep -E "maximum resident|real" results/t7_mem_r.txt
echo "== T7 done $(date)"
} > results/t7_run.log 2>&1
touch results/t7.done
