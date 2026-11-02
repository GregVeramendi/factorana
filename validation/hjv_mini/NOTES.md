# Mini-HJV validation: working notes

## Legacy build on macOS (done)
- Toolchain: Homebrew ROOT 6.40.04, Ipopt 3.14.20 with bundled MUMPS 5.6.2, Open MPI 5.0.11, Apple clang 17.
- Patches: see `legacy_build/PATCHES.md`. Build: `legacy_build/build.sh`, then `legacy_build/build_prob.sh prob_x.cc [mpi]`.

## Tiny check (1 factor, 3 linear measures, N = 2000), `tiny/`
- Under the HJV identification (factor variance fixed at 1, loadings free), legacy and R agree on every
  estimate and SE to the six printed decimals (`tiny/results2/`).
- Legacy objective = R negative log-likelihood minus N * n_factors * 0.5 * log(pi):
  7905.682566 - 6760.952680 = 1144.7299 = 2000 * 0.5723649. Legacy does not divide the Gauss-Hermite
  weights by sqrt(pi). T1 adds this constant back before comparing log-likelihoods.
- Fixing both the factor variance and an anchor loading (`tiny/results/`) gave a singular Hessian in legacy
  (`TDecompLU: matrix is singular`) and SEs of 40 to 70. Legacy fixes the factor SD at 1 by default, so
  passing a `norm` of 1 for a loading over-identifies the model. Use free loadings with legacy.
- `define_model_component()` warns when the intercept column is not named `intercept`. Name it `intercept`.

## Legacy pitfall: `LastModel_*` calls after a skipped `AddModel`
- With `stage = 1`, `AddModel` returns without adding any outcome model, but `LastModel_FixParamValue()` still
  binds to `models.back()`, which is then the last measurement model. In the first Stage 1 run this fixed
  `stresstol_x2` at 0 (the call was written for `HS_track`), costing 1.25 log-likelihood units and moving
  f2-related estimates by up to 1.2 SE. For a non-logit model the `choice` argument is ignored, so there was
  no error. Guard such calls with `if (STAGE > 1)`. `prob_hjv.cc` avoids this by chance: its fixed variable
  is not a regressor of the last measurement model.

## T2, Stage 1 (measurement system), N = 2000, 8 cores both sides: PASS
- 96 of 96 free parameters paired by exact name. Max |estimate gap| = 1.2e-8 (`g_gpahs_intercept`),
  max |gap|/SE = 4e-7, max |SE ratio - 1| = 3.8e-7. Log-likelihood gap after the sqrt(pi) constant = 9e-12.
- Legacy: 35 Ipopt iterations, 56 s on 8 MPI ranks. R: 13 nlminb iterations, 178 s on 8 cores.
  Timing is indicative only: another R test session was running at the same time. T7 measures this properly.

## Pipeline difference: what is in the likelihood when a later outcome group is estimated
- Legacy `Est_outcomes` estimates group g with the measurement models (fixed) and group g only. Earlier outcome
  groups are ignored (`SetIgnoreMod`), so the labor-market group is estimated without the education models.
- `estimate_hjv.R` Stage 4 passes `previous_stage = stage3_result`, so the education components stay in the
  likelihood (fixed) when the labor-market equations are estimated. The factor posterior per person then also
  uses the education outcomes. Both are consistent under the model; the R version uses more information.
- For the comparison, R stage 4 fixes stage 1 only (`run_r.R 4 N 8 1`, file `r_s4_N`), which mirrors legacy.
  `run_r.R 4 N 8 3` (file `r_s4p3_N`) is the `estimate_hjv.R` version. The choice for the re-estimation is the
  user's; see the report.

## T3, education stage (R stage 3 vs legacy group 1), N = 2000: PASS
- 272 parameters paired (96 measurement fixed at identical values + 176 education). Max |estimate gap| = 1.2e-8,
  max |gap|/SE = 4e-7, max |SE ratio - 1| = 3.8e-7, log-likelihood gap after the constant = 1.8e-7.
- Legacy group 1: 11 Ipopt iterations; both groups 69 s on 8 ranks. R stage 3: 6 nlminb iterations, 285 s on 8 cores.

## R package bug found by T3: previous_stage cannot be chained twice
- `define_model_system(previous_stage = )` writes `variance_fixed` and `variance_value` into the factor model it
  stores in the returned system (define_model_system.R lines 272-292). A later call that passes that fit as
  `previous_stage` fails the `identical(prev_ms$factor, factor)` check with "Factor model in previous_stage must
  be identical to current factor model". `estimate_hjv.R` Stage 4 (`previous_stage = stage3_result`) hits this.
- `variance_value` was `c(factor_var_1 = 1)` for a 3-factor model, so the stored value covers one factor only.
- Not fixed here (no package changes in the validation). Workaround used: R stage 4 fixes stage 1 only.

## T3, labor-market stage (R stage 4 with stage 1 fixed vs legacy group 2), N = 2000: PASS
- 114 parameters paired (96 measurement + 18 labor). Max |estimate gap| = 1.2e-8, max |gap|/SE = 4e-7,
  max |SE ratio - 1| = 5.7e-7, log-likelihood gap after the constant = 3.7e-8.
- Legacy group 2: 14 Ipopt iterations. R stage 4: 6 nlminb iterations, 24 s on 8 cores.

## T1 (likelihood, gradient, Hessian at the true parameters, N = 2000), first run
- Log-likelihood: R -53079.9842913558, legacy -53079.9842913557 after the constant; gap 3.8e-11, relative 7e-16.
- Gradient and Hessian: all 194 education and labor parameters agree to about 1e-15 relative (median 2.5e-15),
  including all sigma, loading and beta entries. For the 96 measurement parameters legacy stage 0 returned
  exact zeros, so TestCode did not compute derivatives for group-0 models. Being resolved (see below).
- Timing at N = 2000, single core: legacy likelihood 1.59 s, gradient 2.42 s, Hessian 49.7 s;
  R likelihood + gradient + Hessian 523 s in one call (the per-call split is in T7).
- Legacy stage 5 (PredictFactors) aborts under MPI with MPI_ERR_TRUNCATE; it is run serially instead.

## T5, leg R (R posterior-mode scores fed to both integrators), N = 2000
- Same scores, two integrators: 272 parameters paired, max |estimate gap|/SE = 0.42 (`swesat_loading_3`), SE ratios
  within 2 percent. FAIL against the 0.01 SE criterion: the integrators are not the same given the same scores.
- Cause, from the scores alone: the point-count rules differ for 653 of 2000 observations (32.6 percent). Both use
  1 + 2 floor(SE/SD/0.5). R then sends any factor with SE > SD to 8 points centered at the prior mean (955 factor-
  observations here), while legacy keeps 5 importance-sampling points centered at the score, and its observation-
  level fallback (mean SE/SD > 1.5) fired for nobody. Mean points per observation: R 37.8, legacy 20.2 (full: 512).
- Against full quadrature (176 education parameters): R adaptive mean |gap|/SE 0.20, max 1.07; legacy adaptive
  mean 0.25, max 1.15. So adaptive integration at max 8 points moves estimates by about a fifth of an SE on average.
- Speed: legacy adaptive 3 s (9 iterations, 8 ranks); R adaptive 23 s (8 cores). Full quadrature: 69 s and 285 s.
- The reported log-likelihoods differ by 3487 because R's importance weight carries the SE/sigma jacobian.

## T4 note: legacy `SetIncludePriorFactorScore()` path is not usable
- The legacy ML score run finished 2000 observations in 2 s (serial). With `SetIncludePriorFactorScore()` the run
  produced 3 observations in 10 minutes and was stopped. HJV never used this option (`includePriorFactorScore = 0`
  everywhere), so the T4 "legacy with prior vs R" criterion cannot be evaluated. T4 compares legacy ML scores with
  R posterior-mode scores, both against the true factors.

## T6 scope (change from the plan)
- The plan said 20 replications at N = 5000. R's stage 1 and full-quadrature education stage take 3 to 4 minutes each
  at N = 2000 on 8 cores, so N = 5000 would need about 9 hours. T6 runs 20 replications at N = 2001 (about 3 minutes
  each) with the full-quadrature reference and stage 1 estimated by legacy only; T2 and T3 showed the two codes
  reach the same optimum to 1e-8, so R's factor scores are computed from the legacy stage-1 values placed in an R
  fit object. N = 5000 can follow if the N = 2001 results are close calls.

## T1a (measurement-only build), N = 2000: PASS. T1 complete.
- 96 free parameters. Log-likelihood gap 1.1e-10 (rel 3e-15); gradient max |g_R + g_L| = 7e-13 (rel 7e-14);
  Hessian max |H_R + H_L| = 3.5e-11 (rel 8e-12). With T1b (194 outcome parameters at 1e-15), the likelihood,
  gradient and Hessian code agree for every parameter type.

## T4, factor scores, N = 2000
- Legacy ML scores and R posterior-mode scores differ for 97 percent of people by more than 0.01; mean |diff|
  0.31, 0.24, 0.50 by factor; max 3.0, 2.2, 4.3. RMSE against the true factors: legacy ML 0.65, 0.58, 0.92;
  R posterior mode 0.39, 0.47, 0.57. The R scores are closer to the truth on every factor.
- Legacy SEs average 0.33, 0.53, 0.49 against R observation-only SEs of 0.55, 0.56, 0.86. The legacy SE is shrunk by
  factors of 0.8 until the likelihood ratio at +/- SE is at most 1, so it is not the inverse-Hessian SE.

## T5, leg L (legacy ML scores fed to both integrators), N = 2000
- Given the same scores, R adaptive and legacy adaptive are identical: all 176 education estimates agree to 1e-8.
  No legacy SE reaches 1, so neither fallback fires. Legacy SEs are below 0.5 for 99.9, 59 and 77 percent of people
  on factors 1 to 3, giving one point there and three points otherwise: 3.0 points per person on average against 512.
  For the 206 people (10 percent) with neither cognitive nor leadership measures the legacy SEs average 0.42, 0.62,
  0.51 while R's observation-only SEs average 1.14, 0.78, 1.56: the legacy 0.8-shrink guard makes weakly measured
  factors look precise, and the adaptive rule then integrates them with a single point.
  Combined with leg R, the only integrator difference is the fallback rule.
- Against full quadrature, the legacy production combination (ML scores, adaptive) moves the education estimates by
  0.49 SE on average and up to 4.45 SE (`swesat_loading_1`). With R posterior-mode scores the adaptive fits move by
  0.20 to 0.25 SE on average and at most 1.15 SE. Adaptive integration with 1 point is plug-in estimation.

## T7, legacy per-evaluation timing (single core, full 28-component system at the truth, 8 quadrature points)
- N = 2000: likelihood 1.59 s, likelihood + gradient 2.43 s, likelihood + gradient + Hessian 49.2 s (3 repeats, spread 0.03 s).
- N = 5000: 3.97 s, 6.07 s, 122.8 s. Linear in N. The Hessian dominates: about 20 times the gradient call.
- R per-evaluation timings, convergence timings, scaling and memory: pending (T7 still running).

## T7, R per-evaluation timing (single core, same system, same point), N = 5000, one repeat
- R: likelihood 12.9 s, likelihood + gradient 29.2 s, likelihood + gradient + Hessian 1313 s (290 free parameters).
- Legacy at N = 5000: 4.0 s, 6.1 s, 122.8 s. Ratios R / legacy: 3.2, 4.8 and 10.7.
- The first R N = 2000 line in t7_run.log (all zeros) is a timing-script bug (promise evaluated once); the N = 2000
  rerun with the corrected script runs after the other T7 steps.

## T7, time to convergence, N = 5000, 8 cores (machine otherwise idle)
- Legacy stage 1: 131 s, 38 Ipopt iterations. Legacy stage 2 (education then labor groups): 118 s, 8 + 13 iterations.
- R stage 1: 368 s, 12 nlminb iterations (2.8 times legacy). R education stage: pending.
- R education stage, N = 5000, 8 cores: 674 s, 6 nlminb iterations, against legacy's 118 s for both outcome groups
  (5.7 times). T2 and T3 repeated at N = 5000: pass (estimates to 7e-9, SEs to 5e-7, log-likelihood to 4e-7).
- Scaling, stage 1, N = 2000: legacy 394 s on 1 rank against 56 s on 8 ranks (7.0 times). R: pending.

## R package fix (November 1, 2026): previous_stage chaining
- `define_model_system()` now compares factor models through `.factor_spec()`, which drops the bookkeeping fields
  `variance_fixed`, `variance_value`, `variance_partially_fixed`, `fixed_variance_params`, `fixed_variance_values`;
  `variance_value` stores one value per factor. Regression test `test-previous-stage-chain.R`: 2 tests, 0 failures,
  run in final form. The existing previous_stage test files run after T6 on the idle machine.
- Legacy stage 1, N = 2000: 394 s on 1 rank, 96 s on 4, 56 s on 8 (speedups 4.1 and 7.0). R: pending.
- R stage 1, N = 2000: 1195 s on 1 core against 178 s on 8 (speedup 6.7; the 8-core run had another job alongside).
  Single core, R is 3.0 times legacy on this stage. R 4-core and memory: pending.
- R stage 1, N = 2000: 308 s on 4 cores (speedups 3.9 at 4, 6.7 at 8). Legacy at the same N: 394, 96, 56 s.
- Memory, stage 1, N = 5000, 8 cores, via `/usr/bin/time -l`: the figures are the launcher's peak RSS only
  (mpirun 37 MB, the Rscript master 201 MB), not the ranks or workers, so they do not compare the two codes.
  Wall times from the same runs: legacy 134 s (1048 s CPU), R 372 s (2916 s CPU).
- R per-evaluation timing, N = 2000, corrected script, median of 3: likelihood 5.15 s, + gradient 11.6 s,
  + Hessian 521.8 s. Legacy: 1.59, 2.43, 49.2 s. Ratios 3.2, 4.8, 10.6, the same as at N = 5000.

## T6 run log
- First launch failed at the R scoring step in every replication: `merge(tab, par)` produced `est.x`/`est.y`, so the
  legacy stage-1 values were never written into the template fit. Fixed (`t6_rep.R`), smoke-tested, relaunched from
  replication 1 at 18:5x. Replication 1 complete with all five estimators; R adaptive fits converged.

## T6, 20 replications at N = 2001, education stage (176 parameters), bias and RMSE in units of the full-quadrature SE
| estimator | mean bias | mean |bias| | RMSE | mean |dev from full| | max |dev from full| |
|---|---|---|---|---|---|
| full quadrature (legacy, reference) | -0.039 | 0.234 | 1.096 | 0 | 0 |
| legacy adaptive, legacy ML scores (legacy production) | -0.103 | 0.546 | 1.378 | 0.434 | 5.40 |
| legacy adaptive, R posterior-mode scores | -0.078 | 0.349 | 1.140 | 0.222 | 1.80 |
| R adaptive, R posterior-mode scores (R pipeline) | -0.071 | 0.312 | 1.121 | 0.177 | 1.49 |
| R adaptive, legacy ML scores | -0.103 | 0.546 | 1.378 | 0.434 | 5.40 |
- The two ML-score rows are identical, as T5 leg L predicted. The score type drives the result: ML scores raise mean
  |bias| from 0.23 to 0.55 SE and RMSE from 1.10 to 1.38. With posterior-mode scores the R fallback rule (per factor)
  does slightly better than the legacy rule (per person): 0.312 against 0.349 mean |bias|, max deviation 1.5 against 1.8.
- Largest biases of the legacy production estimator (mean z over 20 replications): swesat_loading_1 -6.74,
  advmath_loading_1 -4.93, adveng_loading_1 -4.20, HS_track_c3_loading_1 -2.68, swesat_sigma +2.29, apply_loading_1
  -2.16. Factor-1 loadings of outcome equations are attenuated. The R pipeline keeps all of these within about 3 SE
  (swesat_loading_1 -2.94, advmath_loading_1 -2.06) and most within 1 SE.
- The full-quadrature reference itself shows mean |bias| 0.23 SE and RMSE 1.10 at this sample size, with
  advmath_loading_1 at -2.15 and swesat_loading_1 at -2.53: finite-sample bias of the ML estimator, present in both codes.
- R adaptive non-convergence: 0 of 40 fits.

## Regression check for the previous_stage fix (idle machine, after T6)
- test-two-stage.R 1 test, test-stale-previous-stage.R 3 tests (1 warning), test-two-stage-se-types.R 7 tests,
  test-previous-stage-chain.R 2 tests: 13 tests, 0 failures, 0 errors. All seven validation steps are complete.

## T7 correction (November 2, 2026): the R timings used an unoptimized build
- Every R timing so far (T1 note, T3/T5 convergence times, T7 table) ran on `src/factorana.so` as compiled by
  `pkgload::load_all()`. pkgload calls `pkgbuild::compile_dll(debug = TRUE)`, whose flags are
  `-UNDEBUG -Wall -pedantic -g -O0`. The legacy build is -O2. Found after the user recalled that in February the R
  package estimated the measurement system slightly faster than legacy.
- Check: the February package (v1.0.2, commit 1cb3991) and v1.4.0 installed with `R CMD INSTALL` (-O2) into scratch
  libraries, timed with `t7b_r_timing.R` (single core, N = 2000, at the truth, medians of 3 for stage 1, one run for
  the full system; other single-core jobs were running, so absolute times carry about 30 percent noise):

  | Build | Stage 1 (96 free): loglik / +grad / +hess | Full system (290 free): loglik / +grad / +hess |
  |---|---|---|
  | legacy, -O2, idle machine (Nov 1) | 0.71 / 1.07 / 7.53 s | 1.59 / 2.42 / 49.7 s |
  | legacy, -O2, loaded machine (Nov 2) | 0.98 / 1.50 / 10.55 s | not rerun |
  | R v1.0.2 (Feb), -O2 | 0.58 / 1.16 / 9.26 s | 1.47 / 2.98 / 71.1 s |
  | R v1.4.0, -O2 | 0.58 / 1.17 / 9.44 s | 1.49 / 3.00 / 70.8 s |
  | R working copy 1.8.0, load_all build (-O0) | 2.55 / 5.81 / 104.8 s | 5.2 / 11.6 / 522 s (T7 table) |

  Both R -O2 builds return the same log-likelihood as the working copy (-34120.425305 stage 1, -53079.984291 full).
  The -O0 build is 4.4, 5.0 and 11.3 times slower than -O2 on stage 1, which reproduces the T7 ratios exactly.
- Consequence: the T7 table, the "3 to 6 times" budget in the report, and the R convergence times quoted in T3/T5
  are wrong by the -O0 factor. Legacy and R at -O2 are close to parity on the likelihood (R faster on stage 1),
  R about 1.2 times slower on the gradient and 1.4 times on the Hessian of the full system.
- Fix applied: `pkgbuild::clean_dll(); pkgbuild::compile_dll(debug = FALSE)` on the working copy (factorana.so shrank
  from 1.55 MB to 0.53 MB). `run_t7r.sh` redoes the R side of T7 on that build; T8 runs afterwards on it.
- No C++ code regression: v1.0.2, v1.4.0 and the current source give the same -O2 times, and the only C++ change
  between v1.4.0 and the current source is the optional per-observation score storage (off by default) plus the
  exploded-logit rank test.
- Caution for anyone repeating the check: `R CMD INSTALL .` on the working copy reused the -O0 objects already in
  `src/` (no compile lines in its log) and timed like the -O0 build (2.55 / 5.77 / 105 s). Use `--preclean`, or
  `pkgbuild::clean_dll()` first. The rebuilt working copy loaded through `load_all(compile = FALSE)` gives
  0.585 / 1.18 / 9.48 s on stage 1, equal to the v1.0.2 and v1.4.0 -O2 installs.

## T9 design (November 2, 2026): latent types and the joint final stage
- Production (`bin/prob_hjv.cc`) uses `TMinLkhd(..., 3 factors, 8 types, ...)`: measurement models fix their type
  intercepts at 0 (`norm[3..9] = 0`), and ONE final estimation group holds `AddTypesModel({"none","constant"})` (a
  logit over types with an intercept and free factor loadings) followed by every education and labor-market model with
  free type intercepts. The mini model had no types and split education and labor into two groups, so T1 to T8 never
  exercised the type code. The labor-stage question in the first report is moot under the production layout.
- Mini T9: `generate_data.R <N> <seed> <T>` draws type with P(type = t | f) = softmax(a_t + L_t f), adds deterministic
  type intercepts to every education and labor component (per alternative for logits), writes `*_t<T>` files and a
  `true_type` column. T = 3 for parity (344 parameters, 50 type-related, in the R full system), T = 8 for timing only.
- Legacy: `prob_hjv_mini.cc -DNTYPES=T` (binaries `prob_t9_*`). R: `build_mini_system(dat, 5, prev, ntypes = T)`
  (stage 5 = education + labor jointly, `use_types = TRUE` on every outcome component). Name map: legacy
  "Type Probability (Choice c)" -> `typeprob_c_intercept`, `type_c_loading_k`; "Typ j loading" in model m -> `m_type_(j+1)_intercept`.
- Steps (`run_t9.sh`): legacy stage 1 (typed build), legacy final stage from its default start (A) and from the truth
  (B), TestCode at the truth vs R (T1 criterion), R stage 1 and R final stage from default and truth starts, compare
  R vs A and R vs B (T3 criterion), per-evaluation timing at 3 and 8 types. Two starts because a 3-type mixture can
  have several local optima; the truth-start pair checks parity, the default-start pair checks whether the codes'
  own initializations land on the same optimum.

## T7 redo, R side on the -O2 build (`run_t7r.sh`, November 2)
- Per evaluation, single core, full system at the truth: N = 2000 (median of 3) 1.489 / 3.026 / 72.9 s; N = 5000 (one
  run) 2.961 / 8.627 / 160.0 s. Legacy (medians of 3): 1.59 / 2.43 / 49.2 s and 4.0 / 6.1 / 123 s.
- To convergence, 8 cores, N = 5000: stage 1 56.5 s (12 nlminb iterations) against legacy 131 s wall (Ipopt 123 s,
  38 iterations); education stage 76.0 s (6 iterations) against legacy Ipopt 79 s (118 s wall with the wage group).
- Scaling, stage 1, N = 2000: 1 / 4 / 8 cores 123.2 / 37.6 / 20.3 s against legacy 394 / 96 / 56 s wall
  (Ipopt 374 / 89 / 52). Speedup 1 to 8 cores: R 6.1, legacy 7.0.
- Adaptive education stage, N = 2000, 8 cores, R scores: 3.9 s wall against legacy 3 s wall (Ipopt 1.1 s).
- REPORT.md T7 section, results row and recommendation 1 rewritten accordingly.

## T8 (November 2): cost against accuracy of adaptive rules, posterior-mode scores, 20 replications at N = 2001
Same replications, stage-1 estimates and full-quadrature reference as T6. All 121 R fits converged (code 0). Points per
person computed from the factor-score SEs by each rule (full quadrature = 512). Seconds: legacy Ipopt time (education
group), R optimizer time (`optimization_time`, 8 workers), both on -O2 builds. Accuracy in units of the full-quadrature SE,
176 education parameters.

| Estimator | Code | Points/person | Seconds | Time / full | Mean abs bias | RMSE | Mean dev from full | Max dev |
|---|---|---|---|---|---|---|---|---|
| Full quadrature | legacy / R | 512 | 38.6 / 30.9 | 1 | 0.234 | 1.096 | 0 | 0 |
| ML scores, legacy rule (production) | legacy | 3.1 | 0.3 | 0.008 | 0.546 | 1.378 | 0.434 | 5.40 |
| Posterior mode, obs-only SE, legacy rule (no fallback fires) | legacy | 17.0 | 1.5 | 0.039 | 0.349 | 1.140 | 0.222 | 1.80 |
| Posterior mode, obs-only SE, R per-factor fallback, t 0.5 (T6 R pipeline) | R | 34.1 | 3.3 | 0.107 | 0.312 | 1.121 | 0.177 | 1.49 |
| Posterior mode, obs-only SE, R per-factor fallback, t 1.0 | R | 14.9 | 1.8 | 0.058 | 0.339 | 1.137 | 0.211 | 1.86 |
| Posterior mode, obs-only SE, per-person fallback (any SE > SD -> 512) | R | 156.7 | 11.8 | 0.382 | 0.283 | 1.107 | 0.140 | 1.06 |
| Posterior mode, posterior SE (include_prior = TRUE), t 0.3 | R | 35.3 | 3.2 | 0.104 | 0.201 | 1.076 | 0.109 | 1.22 |
| Posterior mode, posterior SE, t 0.5 | R | 8.1 | 1.4 | 0.045 | 0.281 | 1.104 | 0.141 | 1.22 |
| Posterior mode, posterior SE, t 1.0 (one point per person) | R | 1.0 | 0.9 | 0.029 | 0.395 | 1.170 | 0.274 | 2.22 |

- The posterior SE (`estimate_factorscores_rcpp(include_prior = TRUE)`) is the decisive change. It is bounded by the
  prior SD, so no fallback rule ever fires: a person with no measures has score 0 and SE 1 and gets 1 + 2 floor(1/t)
  points per factor centered at the prior with unit spread, which is standard quadrature with the IS weight equal to 1.
  Everyone else gets points in proportion to posterior uncertainty.
- Posterior SE, t 0.5: 8 points per person, 4.5 percent of the full-quadrature time, accuracy equal to the strict
  per-person fallback that costs 157 points (mean dev from full 0.141 against 0.140; max 1.22 against 1.06).
- Posterior SE, t 0.3: 35 points per person, 10 percent of the full time, and no worse than full quadrature itself on
  bias and RMSE (0.201 / 1.076 against 0.234 / 1.096; the difference is replication noise).
- The obs-only SE variants are dominated: the T6 R pipeline uses 34 points for a mean deviation of 0.177.
- One point per person (posterior SE, t 1.0) is a Laplace-type plug-in: 3 percent of the full time, mean deviation
  0.27 SE, maximum 2.2 SE. Cheaper than the production setting and far less biased than it, but t 0.5 at 8 points is
  the better trade.
- Production (ML scores, 3 points) stays the fastest and the most biased: mean deviation 0.43 SE, maximum 5.4 SE.
- `t8_summary.R` writes results/t8_summary_2001.csv and results/t8_long_2001.csv; per-replication fits in results/t8/rep*/.

## T9 results so far (November 2)
- Legacy limitation: with `ntypes > 1` and the type model in a skipped outcome group, `Est_measurementsys` throws
  `std::out_of_range` (type_model indexes past the model vector). Stage 1 therefore runs as an untyped model on the
  typed data (`prob_t9_s1u`, `-DDSUF="_t3"`), which is the same estimation problem; 54 s on 8 ranks.
- R package bug, fixed in 1.8.0: `initialize_factor_model_cpp()` allocated the type-probability parameters whenever
  the factor model had `n_types > 1`, while R adds them only when a component has `use_types = TRUE`. The typed
  stage 1 (factor model with 3 types, no typed component) failed with "Fixed values size mismatch". Fix in
  `src/rcpp_interface.cpp`: the typed-component pre-scan now runs before the constructor and `n_types` collapses to 1
  when no component uses types (also avoids looping over indistinguishable types). Test `test-types-unused.R`, 8 passes.
- Stage 1, typed factor model in R vs legacy untyped model on typed data: PASS. 96 parameters, max |diff| 5.1e-9,
  SE ratios within 4e-7, log-likelihood gap 1.8e-12 after the constant. R 20.4 s (8 cores, 13 iterations).
- Legacy final stage (3 types, education + labor jointly, 245 free outcome parameters): default start converges in
  28 Ipopt iterations, 561 s Ipopt time on 8 ranks (601 s wall), objective 49444.125. From the true parameters
  (init_par.txt) Ipopt stops at iteration 0 with "Invalid number in NLP function or derivative detected"; TestCode
  at the same parameters returns finite values (objective 49660.60, 58311 finite dump entries), so the failure is in
  the derivative Ipopt requests first, not in the likelihood. Not yet understood.
- T1 with 3 types at the truth: FAIL. R log-likelihood -53202.597 against legacy -53094.792 after the constant (gap
  107.8; both finite). Legacy stage-0 gradients for the measurement parameters are zero again (thisflag = 1 when
  outcome models are present), so only the log-likelihood is informative here. `t9_bisect.R` zeroes subsets of the
  type parameters on both sides to locate the difference.
- Per-evaluation timing, single core, N = 2000, 3 types (340 free): legacy 5.89 / 9.05 / 204 s, R 3.20 / 6.70 / 256 s.
- Bisection of the typed T1 gap (`t9_bisect.R`, log-likelihoods at the truth with subsets of the type parameters set
  to zero): no gap without type parameters, with type-probability parameters only, with linear or probit type
  intercepts only; the whole gap came from the multinomial logits (HS_track -2.1, application -13.9, switch_major
  -0.3; -107.8 with everything on). A single HS_track choice-1 intercept left R's log-likelihood unchanged, a single
  choice-3 intercept moved it as if applied to every alternative.
- Cause: a second R package bug. The parameter layout defines one type intercept per non-reference alternative
  (`<comp>_c<k>_type_<t>_intercept`), but `Model::Eval()` took one scalar per (model, type) and added it to every
  alternative's index, and the gradient and Hessian treated the block as holding n_types - 1 intercepts. Binary models
  were right; every logit with more than two alternatives was wrong in likelihood, derivatives and
  `fix_type_intercepts(choice = )`. This affects the production final stage (HS track, application, major switching).
- Fix (1.8.0): `Model::Eval()` takes a pointer to per-alternative intercepts; `FactorModel` fills them choice-major
  (`GetTypeInterceptIndex(ityp, imod, choice)`, `TypeInterceptBlocks()`), accumulates one gradient entry per
  alternative (the alternative's intercept-regressor gradient) and the Hessian blocks (diagonal, between alternatives
  of the same type, with factor variances, with base parameters); `rcpp_interface.cpp` and `optimize_model.R` fix
  type intercepts per alternative. After the fix all nine bisection variants agree with legacy to 1e-6 (printed
  precision); test files test-types, test-use-types, test-fix-type-intercepts, test-types-unused,
  test-two-stage-se-types: 0 failures. New `test-types-mnl.R` (likelihood sensitivity per alternative, FD gradient and
  Hessian, recovery, fix_type_intercepts(choice)).
- Legacy start values: `Est_outcomes` multiplies every value read from `init_par.txt` by `loadingMultiplier`
  (= 0.001), so a file holding the true values puts sigmas near 0 and Ipopt stops at iteration 0 with "Invalid number
  in NLP function or derivative". The same happens for the untyped model. Run B therefore uses init values scaled by
  1000. TestCode reads the same file without the multiplier.
- Package tests after the type fixes (NOT_CRAN=true, local): test-fix-type-intercepts 91 expectations, test-types-mnl 11,
  test-use-types 30, test-types 95, test-types-unused 8, test-two-stage-se-types 28: 0 failures, 0 errors, 0 warnings.
  A third repair came out of the new test: `initialize_parameters()` dropped a fixed per-alternative type intercept
  from the parameter vector while the C++ layout keeps the slot; it now keeps it with value 0, as for other models.
- Per-evaluation timing, single core, N = 2000, 8 types (legacy): 19.6 / 29.8 / 1146 s (likelihood / gradient / Hessian).
- T1 with 3 types after the fix (`t1_compare.R 2000 t3 outcome`): log-likelihood gap 1.4e-10 (relative 2.7e-15) over
  the whole 344-parameter system. Over the 244 free outcome parameters (legacy stage 0 returns zero derivatives for the
  measurement models when outcome models are present, and those are fixed in the final stage) the gradients agree to
  1.9e-13 relative and 59,524 of the 59,536 Hessian elements agree to 1e-6. The 12 that differ (6 symmetric pairs) are
  cross terms between a type intercept and a regressor in the three typed LINEAR components (swesat, wage1, wage2),
  e.g. wage1_type_3_intercept x wage1_x1: R -145.36, legacy -19.79. Arbiter: a central finite difference of the
  gradient, which both codes share to 1e-13, reproduces R's values to all printed digits (-145.3604, -3080.7511,
  -68.0392, -2086.5527, 158.0062). So the legacy Hessian is wrong on those cross terms (type "loadings" in the linear
  model); it does not move the optimum, but it perturbs legacy's Newton steps and the inverse-Hessian SEs of the type
  intercepts and regressors in typed linear outcomes.
- First typed R final-stage fits (before the next two repairs): default and truth starts both converged (code 0) to
  the same point, log-likelihood -52910.876, which is 32.6 below legacy's optimum (-52878.315 = -(49444.125 + const)),
  with swesat_sigma on its 0.01 bound. R's likelihood evaluated AT legacy's optimum is -52878.3149 with gradient
  3e-5, so legacy's optimum is a stationary point of R's likelihood too; the R optimizer worked on a mis-set vector.
- Cause (R package, fixed in 1.8.0): with `previous_stage`, `initialize_parameters()` reused the stage-1 vector as is,
  so the stage-5 vector lacked the 8 type-probability parameters (336 instead of 344). `estimate_model_rcpp()` then
  reconciled by name, and because `build_parameter_metadata()` named multinomial logit parameters `<comp>_beta_<cov>_alt<k>`
  while `initialize_parameters()` uses `<comp>_c<k>_<cov>`, the 93 logit parameters lost their values (and the
  reported names switched to the `_alt` form). Fixes: the type-probability block is now inserted after the factor-level
  parameters in the previous_stage path, and the metadata names follow the `_c<k>_` convention; both paths now give
  identical vectors on the typed and untyped systems (344 / 344, 276 / 276, 99 / 99). Tests: previous-stage-chain,
  two-stage, stale-previous-stage, two-stage-se-types, types-mnl, types-unused, fix-type-intercepts, exploded logit
  (two files), mlogit coding: 239 expectations, 0 failures.
- Per-evaluation timing, single core, N = 2000, 8 types (465 free parameters): R 10.9 / 23.4 / 1383 s against legacy
  19.6 / 29.8 / 1146 s (likelihood / gradient / Hessian). Both runs overlapped with other jobs on the machine, so
  the absolute times carry some contention; the ratios (R 0.56, 0.79 and 1.21 times legacy) match the untyped and
  3-type pattern: R faster on the likelihood and gradient, slower on the Hessian.
- Second typed R rerun (after the typeprob insertion and the metadata rename) still stopped 32.2 below legacy, and the
  truth start 18.1 below, with wage2_sigma = -0.754 (impossible under a 0.01 lower bound). Cause: R and C++ disagreed
  on the FREE parameter set, so bounds, gradient and Hessian slots were shifted. Two sources: (i) my metadata rename
  broke the `fix_coefficient()` lookup in setup_parameter_constraints(), which parsed the old `_beta_<cov>_alt<k>`
  names (R left HS_track_c3_x2 free, C++ fixed it: 177 vs 176 free on the untyped education stage); (ii) under
  `previous_stage` the C++ fixed EVERY factor-level parameter, including the eight type-probability parameters that the
  stage-1 anchor never had (R leaves them free): 245 vs 236 on the typed final stage. Both fixed; R and C++ now agree
  on all five systems checked (typed stage 5: 244, typed stage 1: 96, untyped stages 1/3/4: 96/176/18).
- Lesson recorded: whenever R-side naming or constraint code changes, compare `setup_parameter_constraints()$free_idx`
  with `extract_free_params_cpp(fm_ptr, seq_along(p))` before trusting a fit; a mismatch converges with code 0 to a
  wrong point.
- Affected test files after all constraint fixes (NOT_CRAN=true): previous-stage-chain 12, two-stage 11,
  stale-previous-stage 12 (1 expected warning), two-stage-se-types 55, types-mnl 11, types-unused 8,
  fix-type-intercepts 91, types 95, use-types 30, exploded-logit-missing-ranks 7, exploded-nested-logit 27,
  fixed-coefficients 33, mlogit-category-coding 5 expectations: 0 failures, 0 errors.
- Third typed R final-stage run (all fixes in): default start converges in 15 iterations, 429.5 s on 8 cores, to
  log-likelihood -52878.314907; truth start in 12 iterations, 330.6 s, to the same value. Legacy: A (default start)
  28 Ipopt iterations, 561 s; B (truth start) 13 iterations, 305 s; both objective 49444.125249 = the same log-likelihood
  (gap 1.1e-7 after the constant).
- Parameters: R (truth start) vs legacy B: 340 free parameters matched, max |diff| 6.7e-9, max |diff|/SE 0.0000,
  SE ratios within 4.0e-7: PASS at the T3 criterion. R (default start) vs legacy B: max |diff| 1.8e-6. Legacy A differs
  from B only by relabeling types 2 and 3: after swapping the labels every type-2/type-3 parameter agrees to 0.0e+00 and
  the non-type parameters are identical (typeprob intercepts A = (-0.0028, 1.0338), B = (1.0338, -0.0028)). So all
  four fits are one optimum up to the label permutation; R's default initialization lands on the truth's labeling here.
- Per-evaluation times with types (single core, N = 2000): 3 types (340 free) legacy 5.89 / 9.05 / 204 s, R
  3.20 / 6.70 / 256 s; 8 types (465 free) legacy 19.6 / 29.8 / 1146 s, R 10.9 / 23.4 / 1383 s. Joint typed final stage
  to convergence (3 types, 8 cores): legacy 561 s (28 it.), R 430 s (15 it.).
- Full package test suite on the final code (NOT_CRAN=true, pkgload build at -O2): 46 files, 265 tests, 1325 expectations,
  0 failures, 0 errors, 2 skipped, 3 warnings (the expected "Skipping stale previous_stage names" messages).

## T8 addendum (November 2): the threshold curve with the posterior SE (20 replications, N = 2001)
| threshold t | points / person | optimizer s | time / full | mean abs bias | RMSE | mean dev from full | max dev |
|---|---|---|---|---|---|---|---|
| 0.3 | 35.3 | 3.2 | 0.104 | 0.201 | 1.076 | 0.109 | 1.22 |
| 0.5 | 8.1 | 1.4 | 0.045 | 0.281 | 1.104 | 0.141 | 1.22 |
| 0.7 | 1.2 | 0.9 | 0.029 | 0.394 | 1.169 | 0.273 | 2.22 |
| 0.85 | 1.0 | 0.9 | 0.029 | 0.395 | 1.170 | 0.274 | 2.22 |
| 1.0 | 1.0 | 0.9 | 0.029 | 0.395 | 1.170 | 0.274 | 2.22 |
- With the posterior SE a factor gets 1 + 2 floor(SE / t) points (SE <= 1): at t = 0.5 a factor with SE in [0.5, 1) gets 3
  points and the clamped SE = 1 of an unmeasured person gets 5; at t = 0.7 almost every factor is at one point; at
  t >= 1 everyone is at one point (the floor), so thresholds above 1 change nothing.
- Going from 0.5 to the floor saves a third of an already small time (1.4 s to 0.9 s, 4.5 to 2.9 percent of full
  quadrature) and doubles the deviation from full quadrature (0.14 to 0.27 SE on average, 1.2 to 2.2 SE at most).
  The fixed per-person cost (every component evaluated at least once) dominates at that end.
