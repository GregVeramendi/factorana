#!/bin/bash
# T3: legacy stage 2 (education, then labor, measurement fixed), then R stage 3 and stage 4. Sequential: 8 cores each.
cd "$(dirname "$0")"
N=2000
{ echo "== legacy stage 2 $(date)"; ./run_legacy.sh 2 $N 8;
  echo "== R stage 3 $(date)"; Rscript run_r.R 3 $N 8 2>&1 | grep -vE "^Test passed|^$";
  echo "== R stage 4 $(date)"; Rscript run_r.R 4 $N 8 2>&1 | grep -vE "^Test passed|^$";
  echo "== done $(date)"; } > results/t3_run.log 2>&1
touch results/t3.done
