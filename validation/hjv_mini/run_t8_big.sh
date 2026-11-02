#!/bin/bash
# T8 at scale: N = 20000 (data/mini_20000.csv, seed 1). Stage 1 in both codes, posterior-mode scores (R) and ML scores
# (legacy), then the adaptive education stage with 8 workers / ranks in both codes and single-core evaluation times.
# Sequential so that no two 8-core jobs overlap. Usage: run_t8_big.sh
cd "$(dirname "$0")"; N=20000
{
echo "== T8 big start $(date)"
# ---- R: stage 1, scores, adaptive fits ----
[ -f results/r_s1_$N.rds ] || Rscript run_r.R 1 $N 8 2>&1 | grep -E "^R stage"
Rscript t8_big_r.R $N scores 2>&1 | grep -E "^R "
rm -f results/t8big/r_fits_$N.csv
for cfg in "post 0.5 8" "post 1.0 8" "post 1.0 1" "post 0.5 1"; do Rscript t8_big_r.R $N adapt $cfg 2>&1 | grep -E "^R adaptive|Error"; done
# ---- legacy: build, stage 1, ML scores, adaptive fits ----
OUT=prob_big_s1 legacy_build/build_prob.sh prob_hjv_mini.cc mpi -DSTAGE=1 -DNOBS=$N > /dev/null 2>&1
OUT=prob_big_s5 legacy_build/build_prob.sh prob_hjv_mini.cc serial -DSTAGE=5 -DNOBS=$N > /dev/null 2>&1
OUT=prob_big_s2aL legacy_build/build_prob.sh prob_hjv_mini.cc mpi -DSTAGE=2 -DNOBS=$N '-DADAPT="L"' > /dev/null 2>&1
OUT=prob_big_s2aP05 legacy_build/build_prob.sh prob_hjv_mini.cc mpi -DSTAGE=2 -DNOBS=$N '-DADAPT="P05"' -DATHRESH=0.5 > /dev/null 2>&1
OUT=prob_big_s2aP10 legacy_build/build_prob.sh prob_hjv_mini.cc mpi -DSTAGE=2 -DNOBS=$N '-DADAPT="P10"' -DATHRESH=1.0 > /dev/null 2>&1
OUT=prob_big_s0aL legacy_build/build_prob.sh prob_hjv_mini.cc serial -DSTAGE=0 -DNOBS=$N '-DADAPT="L"' > /dev/null 2>&1
OUT=prob_big_s0aP05 legacy_build/build_prob.sh prob_hjv_mini.cc serial -DSTAGE=0 -DNOBS=$N '-DADAPT="P05"' -DATHRESH=0.5 > /dev/null 2>&1
OUT=prob_big_s0aP10 legacy_build/build_prob.sh prob_hjv_mini.cc serial -DSTAGE=0 -DNOBS=$N '-DADAPT="P10"' -DATHRESH=1.0 > /dev/null 2>&1
OUT=prob_big_s0f legacy_build/build_prob.sh prob_hjv_mini.cc serial -DSTAGE=0 -DNOBS=$N > /dev/null 2>&1
d1=results/legacy_s1_$N; mkdir -p $d1
if [ ! -f $d1/meas_par.txt ]; then s=$(date +%s); mpirun -np 8 ./prob_big_s1_mpi > $d1/run.log 2>&1; echo "legacy stage 1 N=$N: exit $? secs $(( $(date +%s) - s )) $(grep -m1 'Number of Iterations' $d1/run.log)"; fi
d5=results/legacy_s5_$N; mkdir -p $d5; cp $d1/meas_par.txt $d5/; s=$(date +%s); ./prob_big_s5 > $d5/run.log 2>&1; echo "legacy ML scores: exit $? secs $(( $(date +%s) - s ))"
for tag in L P05 P10; do d=results/legacy_s2a${tag}_$N; rm -rf $d; mkdir -p $d; cp $d1/meas_par.txt $d/
  if [ $tag = L ]; then cp $d5/factor_predictions.txt $d/; else cp results/t8big/factor_predictions_P_$N.txt $d/factor_predictions.txt; fi
  s=$(date +%s); mpirun -np 8 ./prob_big_s2a${tag}_mpi > $d/run.log 2>&1; echo "legacy adaptive $tag N=$N: exit $? wall $(( $(date +%s) - s )) s; $(grep -m1 'Number of Iterations' $d/run.log | tr -s ' '); $(grep -m1 'Total seconds in IPOPT' $d/run.log | tr -s ' ')"
done
# ---- single-core evaluation times: legacy TestCode (init_par at the truth from the s2aL log) and R ----
d0=results/legacy_s0_$N; mkdir -p $d0; Rscript make_init_par.R $N results/legacy_s2aL_$N/run.log $d0 294 2>&1 | grep -E "^wrote"
for tag in aL aP05 aP10; do d=results/legacy_s0${tag}_$N; mkdir -p $d; cp $d0/init_par.txt $d/
  if [ $tag = aL ]; then cp $d5/factor_predictions.txt $d/; else cp results/t8big/factor_predictions_P_$N.txt $d/factor_predictions.txt; fi
  ./prob_big_s0${tag} > $d/run.log 2>&1; echo "legacy eval N=$N $tag: $(grep -E 'Calculating the' $d/run.log | tr -s ' ' | tr '\n' '|')"
done
./prob_big_s0f > $d0/run.log 2>&1; echo "legacy eval N=$N full: $(grep -E 'Calculating the' $d0/run.log | tr -s ' ' | tr '\n' '|')"
Rscript t8_big_r.R $N evaltime 1 2>&1 | grep -E "^R eval"
Rscript t8_big_r.R $N adapt ml 0.5 8 2>&1 | grep -E "^R adaptive|Error"
echo "== T8 big done $(date)"
} > results/t8_big_run.log 2>&1
touch results/t8_big.done
