# Patches applied to the legacy sources for the macOS build

Sources are copied from `libraries/` and `include/`; the originals are untouched.

1. `TMinLkhd.cc:4585`: linear solver `ma57` -> `mumps`. HSL is not installed; Homebrew Ipopt 3.14.20 bundles MUMPS 5.6.2.
   Changes the linear solver inside each Ipopt iteration, not the objective or the optimizer.
2. `minlkhd_nlp.cc`: `eval_g` and `eval_jac_g` return `true` instead of `false`.
   The problem has no constraints (m = 0). Ipopt 3.14 calls `eval_jac_g` during setup and treats `false` as an invalid NLP;
   3.11 did not. No effect on the optimization.
3. Build flags: `-std=c++20` (required by ROOT 6.40 headers), ROOT 6 library list from `root-config --libs` plus `-lMinuit -lMathCore`.
   `-lCint` and the KNITRO/HSL/ASL libraries from `bin/makeprob` are dropped.
4. `TMinLkhd.cc` `TestCode()` (stage 0): after the timing block, write `t1_dump.txt` (objective, gradient, packed
   upper-triangle Hessian over free parameters, 17 significant digits) and `t1_freeparlist.txt` to the working
   directory. Diagnostic output only; no computation is changed.
5. `TMinLkhd.cc` `TestCode()`: `return 0` after the dump, so stage 0 does not continue into the legacy finite-difference
   gradient checks (which ran for minutes and printed `nan`). Stage 0 is only used for the T1 dump and T7 timing.

## 6. TestCode: likelihood-only shortcut (November 2, 2026)
With the environment variable `LKHD_ONLY` set, `TMinLkhd::TestCode()` writes the objective to `t1_dump.txt` right after the
likelihood evaluation and returns, skipping the gradient and Hessian. Used by `t9_bisect.R` to compare log-likelihoods at
many parameter vectors cheaply. No effect when the variable is unset.
