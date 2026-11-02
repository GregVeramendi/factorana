#!/bin/bash
# T7 redo, R side only, on an optimized build. The first T7 run timed the shared object that pkgload::load_all()
# compiles, which carries pkgbuild's debug flags (-g -O0); rebuild first with
#   Rscript -e 'pkgbuild::clean_dll("../../factorana_R"); pkgbuild::compile_dll("../../factorana_R", debug = FALSE)'
# Same steps and sizes as run_t7.sh. Run on an otherwise idle machine.
cd "$(dirname "$0")"
{
echo "== R per-evaluation timing, single core $(date)"
Rscript t7_r_timing.R 2000 3 2>&1 | grep -E "^R N="; Rscript t7_r_timing.R 5000 1 2>&1 | grep -E "^R N="
echo "== R convergence timing, N = 5000, 8 cores $(date)"
Rscript run_r.R 1 5000 8 2>&1 | grep -E "^R stage"; Rscript run_r.R 3 5000 8 2>&1 | grep -E "^R stage"
echo "== R scaling, stage 1, N = 2000 $(date)"
for nc in 1 4 8; do Rscript run_r.R 1 2000 $nc 2>&1 | grep -E "^R stage"; done
echo "== R adaptive education stage, N = 2000, 8 cores $(date)"
Rscript t5_r_adapt.R 2000 R 8 2>&1 | grep -E "^R stage"
echo "== T7r done $(date)"
} > results/t7r_run.log 2>&1
touch results/t7r.done
