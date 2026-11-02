#!/bin/bash
# T6: replications. N = 2001 (distinct file names), seeds 1..REPS. Per replication: legacy stage 1 and full-quadrature
# education stage (reference), legacy ML scores, R posterior-mode scores, then four adaptive education fits:
# legacy integrator x {ML, R} scores and R integrator x {ML, R} scores. 8 cores per step, sequential.
cd "$(dirname "$0")"; N=2001; REPS="${1:-20}"
for v in s1 s2 s5 s2aL s2aR; do :; done
OUT=prob_t6_s1 legacy_build/build_prob.sh prob_hjv_mini.cc mpi -DSTAGE=1 -DNOBS=$N >/dev/null
OUT=prob_t6_s2 legacy_build/build_prob.sh prob_hjv_mini.cc mpi -DSTAGE=2 -DNOBS=$N >/dev/null
OUT=prob_t6_s5 legacy_build/build_prob.sh prob_hjv_mini.cc serial -DSTAGE=5 -DNOBS=$N >/dev/null
OUT=prob_t6_s2aL legacy_build/build_prob.sh prob_hjv_mini.cc mpi -DSTAGE=2 -DNOBS=$N '-DADAPT="L"' >/dev/null
OUT=prob_t6_s2aR legacy_build/build_prob.sh prob_hjv_mini.cc mpi -DSTAGE=2 -DNOBS=$N '-DADAPT="R"' >/dev/null
mkdir -p results/t6
{
for rep in $(seq 1 "$REPS"); do
  echo "== rep $rep $(date)"; out=results/t6/rep$rep; mkdir -p $out
  Rscript generate_data.R $N $rep > /dev/null 2>&1; cp data/truth_$N.rds $out/
  for d in legacy_s1_$N legacy_s2_$N legacy_s5_$N legacy_s2aL_$N legacy_s2aR_$N; do rm -rf results/$d; mkdir -p results/$d; done
  mpirun -np 8 ./prob_t6_s1_mpi > results/legacy_s1_$N/run.log 2>&1 || echo "legacy s1 failed"
  cp results/legacy_s1_$N/meas_par.txt results/legacy_s2_$N/; cp results/legacy_s1_$N/meas_par.txt results/legacy_s5_$N/
  mpirun -np 8 ./prob_t6_s2_mpi > results/legacy_s2_$N/run.log 2>&1 || echo "legacy s2 failed"
  ./prob_t6_s5 > results/legacy_s5_$N/run.log 2>&1 || echo "legacy s5 failed"
  [ $rep -eq 1 ] && Rscript t6_rep.R $N $rep fit1 2>&1 | grep -vE "^Test passed"
  Rscript t6_rep.R $N $rep scores 2>&1 | grep -vE "^Test passed"
  cp results/legacy_s1_$N/meas_par.txt results/legacy_s2aL_$N/; cp results/legacy_s5_$N/factor_predictions.txt results/legacy_s2aL_$N/
  cp results/legacy_s1_$N/meas_par.txt results/legacy_s2aR_$N/; cp results/t6_tmp/factor_predictions_R.txt results/legacy_s2aR_$N/factor_predictions.txt
  mpirun -np 8 ./prob_t6_s2aL_mpi > results/legacy_s2aL_$N/run.log 2>&1 || echo "legacy adaptive L failed"
  mpirun -np 8 ./prob_t6_s2aR_mpi > results/legacy_s2aR_$N/run.log 2>&1 || echo "legacy adaptive R failed"
  Rscript t6_rep.R $N $rep adapt L 2>&1 | grep -vE "^Test passed"; Rscript t6_rep.R $N $rep adapt R 2>&1 | grep -vE "^Test passed"
  for d in legacy_s1_$N legacy_s2_$N legacy_s5_$N legacy_s2aL_$N legacy_s2aR_$N; do mkdir -p $out/$d; cp results/$d/*.txt results/$d/run.log $out/$d/ 2>/dev/null; done
  cp results/t6_tmp/r_adapt_L.csv results/t6_tmp/r_adapt_R.csv results/t6_tmp/r_scores.rds $out/
done
echo "== T6 done $(date)"
} > results/t6_run.log 2>&1
touch results/t6.done
