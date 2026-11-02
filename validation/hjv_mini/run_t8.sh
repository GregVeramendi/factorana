#!/bin/bash
# T8: cost against accuracy of adaptive rules with posterior-mode scores, on the T6 replications (N = 2001, seeds 1..REPS).
# Per replication: regenerate the data (same seed), posterior SEs, then one R adaptive education fit per variant,
# each timed. The R full-quadrature fit is timed on replication 1 only. 8 cores per step, sequential.
# Variants "<se>_t<thresh>_m<max>": se = obs (observation-only SE, R per-factor fallback), obsP (observation-only SE,
# person-level fallback to full quadrature), post (posterior SE, no fallback can fire).
cd "$(dirname "$0")"; N=2001; REPS="${1:-20}"
VARIANTS="obs_t0.5_m8 obsP_t0.5_m8 post_t0.5_m8 post_t0.3_m8 obs_t1.0_m8 post_t1.0_m8"
mkdir -p results/t8
{
for rep in $(seq 1 "$REPS"); do
  echo "== rep $rep $(date)"
  Rscript generate_data.R $N $rep > /dev/null 2>&1
  Rscript t8_rep.R $N $rep postse 2>&1 | grep -vE "^Test passed"
  [ $rep -eq 1 ] && Rscript t8_rep.R $N $rep full 2>&1 | grep -vE "^Test passed"
  for v in $VARIANTS; do Rscript t8_rep.R $N $rep adapt $v 2>&1 | grep -vE "^Test passed"; done
done
echo "== T8 done $(date)"
} > results/t8_run.log 2>&1
touch results/t8.done
