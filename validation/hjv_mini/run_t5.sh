#!/bin/bash
# T5: 2 x 2 (score source L = legacy ML scores, R = R posterior-mode scores) x (integrator legacy, R). Sequential, 8 cores each.
cd "$(dirname "$0")"; N=2000
{
for src in L R; do
  d="results/legacy_s2a${src}_${N}"; mkdir -p "$d"; cp results/legacy_s1_${N}/meas_par.txt "$d/"
  if [ "$src" = L ]; then cp results/legacy_s5_${N}/factor_predictions.txt "$d/"; else Rscript make_factor_predictions.R $N "$d/factor_predictions.txt"; fi
  echo "== legacy adaptive, scores $src  $(date)"; s=$(date +%s); mpirun -np 8 ./prob_mini_s2a${src}_${N}_mpi > "$d/run.log" 2>&1; echo "exit=$? secs=$(( $(date +%s) - s ))"
  grep -E "Number of Iterations|EXIT" "$d/run.log" | head -4
  echo "== R adaptive, scores $src  $(date)"; Rscript t5_r_adapt.R $N $src 8 2>&1 | grep -vE "^Test passed|^$"
  echo "== compare R adaptive vs legacy adaptive (group 1 objective)"; Rscript compare_stage.R s3a${src} s2a${src} $N 1 2>&1 | grep -vE "^Test passed" | head -8
  echo "== adaptive vs full quadrature"; Rscript t5_vs_full.R $N $src 2>&1 | grep -vE "^Test passed"
  echo "== points per observation"; Rscript t5_points.R $N $src 2>&1 | grep -vE "^Test passed"
done
echo "== T5 done $(date)"
} > results/t5_run.log 2>&1
touch results/t5.done
