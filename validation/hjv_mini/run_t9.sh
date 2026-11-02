#!/bin/bash
# T9: latent types and the joint final stage (production layout: education + labor market in one estimation group with
# T types; the type probability logit has an intercept and free factor loadings; measurement models carry no types).
# N = 2000, seed 1, T = 3 for the parity checks (T1-style evaluation, stage-1 and final-stage estimation), T = 8 for
# per-evaluation timing only. Binaries prob_t9_* must exist (see NOTES.md, T9). 8 cores for estimation, serial for TestCode.
cd "$(dirname "$0")"; N=2000; T=3
{
echo "== T9 start $(date)"
Rscript generate_data.R $N 1 $T; Rscript generate_data.R $N 1 8
# 1. legacy stage 1 as an untyped model on the typed data (prob_t9_s1u, -DDSUF="_t3"): types do not enter the
#    measurement system, and the legacy code aborts (std::out_of_range) when stage 1 runs with ntypes > 1 because the
#    type model sits in a skipped outcome group. Parameter indices match the typed build (type norms are not parameters).
d1=results/legacy_s1_t${T}_$N; rm -rf $d1; mkdir -p $d1
s=$(date +%s); mpirun -np 8 ./prob_t9_s1u_mpi > $d1/run.log 2>&1; echo "legacy stage 1 (untyped model, typed data): exit $? secs $(( $(date +%s) - s ))"
grep -E "Number of Iterations|Total seconds in IPOPT" $d1/run.log | head -2
# 2. legacy final stage, run A: legacy default start
d2=results/legacy_s2t${T}_$N; rm -rf $d2 results/legacy_s2t${T}A_$N results/legacy_s2t${T}B_$N; mkdir -p $d2; cp $d1/meas_par.txt $d2/
s=$(date +%s); mpirun -np 8 ./prob_t9_s2t${T}_mpi > $d2/run.log 2>&1; echo "legacy final stage A (default start): exit $? secs $(( $(date +%s) - s ))"
grep -E "Number of Iterations|Total seconds in IPOPT|Objective\.\.\." $d2/run.log | head -3
cp -r $d2 results/legacy_s2t${T}A_$N
# 3. index map and true-value init file from run A's parameter block
d0=results/legacy_s0t${T}_$N; rm -rf $d0; mkdir -p $d0
Rscript make_init_par.R $N $d2/run.log $d0 0 _t$T 2>&1 | grep -vE "^Test passed"
# 4. legacy final stage, run B: start at the truth (init_par.txt), same measurement system
rm -rf $d2; mkdir -p $d2; cp $d1/meas_par.txt $d2/; cp $d0/init_par.txt $d2/
s=$(date +%s); mpirun -np 8 ./prob_t9_s2t${T}_mpi > $d2/run.log 2>&1; echo "legacy final stage B (truth start): exit $? secs $(( $(date +%s) - s ))"
grep -E "Number of Iterations|Total seconds in IPOPT|Objective\.\.\." $d2/run.log | head -3
mv $d2 results/legacy_s2t${T}B_$N
# 5. T1 with types: legacy TestCode (serial) and R evaluation at the truth, then compare
s=$(date +%s); ./prob_t9_s0t${T} > $d0/run.log 2>&1; echo "legacy TestCode typed: exit $? secs $(( $(date +%s) - s ))"
grep -E "Calculating the (Likelihood|Gradient|Hessian)" $d0/run.log | tail -3
Rscript t1_r_eval.R $N full $T 2>&1 | grep -E "^R T1|^free"
Rscript t1_compare.R $N t$T 2>&1 | grep -vE "^Test passed"
# 6. R stage 1 (typed factor model) and the joint final stage from the default start (vs A) and from the truth (vs B)
Rscript run_r.R 1 $N 8 1 $T 2>&1 | grep -E "^R stage"
Rscript run_r.R 5 $N 8 1 $T 2>&1 | grep -E "^R stage"
Rscript run_r.R 5 $N 8 1 $T truth 2>&1 | grep -E "^R stage"
echo "-- stage 1: R (typed factor model) vs legacy (untyped model, typed data)"; Rscript compare_stage.R s1t$T s1_t$T $N 2>&1 | grep -E "^(matched|estimates|SEs|loglik|pass|  )"
echo "-- final stage, default starts: R s5t$T vs legacy A"; Rscript compare_stage.R s5t$T s2t${T}A $N 1 2>&1 | grep -E "^(matched|estimates|SEs|loglik|pass|  )"
echo "-- final stage, truth starts: R s5t${T}T vs legacy B"; Rscript compare_stage.R s5t${T}T s2t${T}B $N 1 2>&1 | grep -E "^(matched|estimates|SEs|loglik|pass|  )"
echo "-- final stage, cross: R truth start vs legacy default start"; Rscript compare_stage.R s5t${T}T s2t${T}A $N 1 2>&1 | grep -E "^(estimates|loglik)"
# 7. per-evaluation timing with 3 and 8 types, single core (legacy TestCode prints its own; R via t7b)
Rscript t7b_r_timing.R $N all 1 "" $T 2>&1 | grep "^R \["
d8=results/legacy_s0t8_$N; mkdir -p $d8; ./prob_t9_s0t8 > $d8/run.log 2>&1; echo "legacy TestCode 8 types: exit $?"
grep -E "Calculating the (Likelihood|Gradient|Hessian)" $d8/run.log | tail -3
Rscript t7b_r_timing.R $N all 1 "" 8 2>&1 | grep "^R \["
echo "== T9 done $(date)"
} > results/t9_run.log 2>&1
touch results/t9.done
