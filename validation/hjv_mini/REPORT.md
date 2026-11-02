# Legacy C++ versus R package on a mini-HJV model

November 1 to 2, 2026. Nine steps: T1 to T7 on November 1 (T7 redone on an optimized build on November 2), T8 (cost
against accuracy of the adaptive rules) and T9 (latent types and the joint final stage) on November 2. T9 exposed
four defects in the R package's typed two-stage path, fixed in factorana 1.8.0 (uncommitted at the time of writing).
Numbers and run logs: `NOTES.md` and `results/`.

## Question

Can the HJV re-estimation move from the legacy C++ code to the R package? We answer this stage by stage on a
simulated model with the same components as HJV, estimated by both codes on the same machine with 8 cores.

## The mini model

Three factors with variances fixed at 1 and free loadings, as in HJV. 28 components and 294 parameters
(HJV: 40 components, 1,253 parameters). Stage 1 has 14 measurements: eight grades on all three factors,
including two sports grades that load mostly on factor 2, three cognitive measures on factor 1, leadership and
stress tolerance on factors 1 and 2 (missing for 38 percent of people), and a probit for the missing indicator.
The education stage has 12 components: two advanced-course probits, a 4-choice track logit with a fixed
coefficient, an application probit, a 5-choice exploded nested logit with up to 5 ranks, a selection probit and a
test score, an enrollment probit, a 4-choice switching logit and three graduation probits. The labor stage has two
wage equations. Each person has three covariates plus the endogenous regressors HJV uses. N = 2,000 unless noted.

## Results

| Step | Question | Result |
|---|---|---|
| T1 | Same likelihood, gradient and Hessian at the true parameters? | Yes. Log-likelihood to 1e-10; every gradient and Hessian element to about 1e-13 relative, all 290 free parameters |
| T2 | Same Stage 1 estimates? | Yes. 96 parameters; estimates to 1e-8, SEs to 4e-7, log-likelihood to 9e-12 |
| T3 | Same education and labor estimates with Stage 1 fixed? | Yes. 176 and 18 parameters; same precision as T2 |
| T4 | Same factor scores? | No, by design. Legacy: maximum likelihood with a shrunk SE. R: posterior mode. RMSE against the true factors 0.65, 0.58, 0.92 (legacy) against 0.39, 0.47, 0.57 (R) |
| T5 | Same adaptive integrator given the same scores? | Identical with legacy scores. With R scores, estimates differ by up to 0.42 SE because the fallback rules differ (factor level in R, person level in legacy) |
| T6 | Which combination is least biased? | Posterior-mode scores with R's fallback rule: mean absolute bias 0.31 SE and RMSE 1.12 against 0.23 and 1.10 for full quadrature. The legacy production combination (ML scores) gives 0.55 and 1.38, with factor loadings biased by up to 6.7 SE |
| T7 | Speed | Near parity per evaluation (R 0.7 to 1.5 times legacy's time on likelihood, gradient and Hessian) and R 2 to 3 times faster to convergence, because `nlminb` needs a third of Ipopt's iterations. The first version of this row (R 3 to 11 times slower) timed an unoptimized build, see "Speed" |
| T8 | Which adaptive rule gives the best accuracy per unit of time? | Posterior-mode scores with the posterior SE and threshold 0.5: 8 points per person, 4.5 percent of the full-quadrature time, 0.14 SE mean deviation from full quadrature. Threshold 0.3 (35 points) matches full quadrature. Observation-only-SE rules are dominated |
| T9 | Do the codes agree with latent types in a joint final stage? | Yes, with factorana 1.8.0: 340 parameters to 7e-9 and SEs to 4e-7, up to a permutation of type labels; R converges faster. Versions before 1.8.0 mis-evaluate typed multinomial logits and mis-set a typed stage after an untyped one |

The log-likelihoods differ by a constant: legacy does not divide the Gauss-Hermite weights by sqrt(pi), so its
objective is lower by N times n_factors times 0.5 log(pi). Gradients and Hessians are unaffected.

## What the adaptive comparison shows

With the legacy production setting (ML scores, threshold 0.5, at most 8 points) the integration uses 3.0 points per
person on average, against 512 for full quadrature, and one point on factor 1 for 99.9 percent of people. The
education estimates then differ from the full-quadrature estimates by 0.49 SE on average and by up to 4.45 SE.
With R's posterior-mode scores the adaptive fits stay within 0.20 to 0.25 SE on average and 1.15 SE at most.

The reason is the legacy SE. It is the inverse-Hessian SE shrunk by factors of 0.8 until the likelihood ratio at
one SE is at most 1. For the 10 percent of people with neither cognitive nor leadership measures, legacy SEs
average 0.42, 0.62 and 0.51 where R's observation-only SEs average 1.14, 0.78 and 1.56. Small SEs give one
integration point, so the legacy pipeline is close to plug-in estimation for weakly measured people.

## Accuracy over 20 replications (T6)

N = 2,001 per replication, 176 education parameters, bias and RMSE in units of the full-quadrature SE. Stage 1 and
the full-quadrature reference were estimated by legacy (the codes agree to 1e-8); R's scores used those stage-1 values.

| Estimator | Points per person | Optimizer seconds | Mean bias | Mean absolute bias | RMSE | Mean deviation from full | Max deviation from full |
|---|---|---|---|---|---|---|---|
| Full quadrature (reference) | 512 | 38.6 legacy, 30.9 R | -0.04 | 0.23 | 1.10 | 0 | 0 |
| Legacy adaptive, legacy ML scores (legacy production) | 3.1 | 0.3 | -0.10 | 0.55 | 1.38 | 0.43 | 5.4 |
| Legacy adaptive, R posterior-mode scores | 17.0 | 1.5 | -0.08 | 0.35 | 1.14 | 0.22 | 1.8 |
| R adaptive, R posterior-mode scores (R pipeline, obs-only SE) | 34.1 | 3.3 | -0.07 | 0.31 | 1.12 | 0.18 | 1.5 |
| R adaptive, legacy ML scores | 3.1 | not timed | -0.10 | 0.55 | 1.38 | 0.43 | 5.4 |

Points per person follow from each rule applied to the factor-score SEs; seconds are the optimizer's own time for the
education stage (Ipopt for legacy, `nlminb` with 8 workers for R), averaged over the replications, from T8.

Three results. First, the score type decides the outcome: ML scores raise the mean absolute bias from 0.23 to 0.55 SE,
and the two ML-score rows are identical because the integrators coincide once no fallback fires. Second, with
posterior-mode scores R's per-factor fallback rule does slightly better than legacy's per-person rule (0.31 against
0.35, maximum deviation 1.5 against 1.8 SE). Third, the legacy production estimator attenuates factor-1 loadings of
the outcome equations: over 20 replications `swesat_loading_1` averages -6.7 SE from truth, `advmath_loading_1` -4.9,
`adveng_loading_1` -4.2, `HS_track_c3_loading_1` -2.7. The R pipeline keeps those at -2.9, -2.1, -1.3 and -1.1, and
the full-quadrature reference at -2.5, -2.2, -1.3 and -0.8, so part of the remaining gap is finite-sample bias of the
ML estimator itself. All 40 R adaptive fits converged.

## Cost against accuracy of the adaptive rules (T8)

The T6 rows differ in two things at once: which SE goes into the point-count formula, and what happens when a factor is
poorly determined. T8 separates them on the same 20 replications. Every estimator uses the R posterior-mode scores; the
rows differ in the SE (observation-only, as in the R pipeline so far, or the posterior SE from
`estimate_factorscores_rcpp(include_prior = TRUE)`), in the fallback (R's per-factor rule, legacy's per-person rule,
or a strict per-person rule that sends anyone with a factor SE above the prior SD to full quadrature), and in the
threshold t of 1 + 2 floor(SE / SD / t) points per factor, capped at 8. All 121 R fits converged.

| SE in the rule | Fallback | t | Points per person | Seconds | Time / full | Mean absolute bias | RMSE | Mean deviation from full | Max |
|---|---|---|---|---|---|---|---|---|---|
| (full quadrature, legacy / R) | | | 512 | 38.6 / 30.9 | 1 | 0.23 | 1.10 | 0 | 0 |
| (production: ML scores, legacy SE) | legacy | 0.5 | 3.1 | 0.3 | 0.008 | 0.55 | 1.38 | 0.43 | 5.4 |
| observation-only | legacy per-person (never fires) | 0.5 | 17.0 | 1.5 | 0.039 | 0.35 | 1.14 | 0.22 | 1.8 |
| observation-only | R per-factor | 0.5 | 34.1 | 3.3 | 0.107 | 0.31 | 1.12 | 0.18 | 1.5 |
| observation-only | R per-factor | 1.0 | 14.9 | 1.8 | 0.058 | 0.34 | 1.14 | 0.21 | 1.9 |
| observation-only | strict per-person (30 percent of persons at 512) | 0.5 | 156.7 | 11.8 | 0.382 | 0.28 | 1.11 | 0.14 | 1.1 |
| posterior | none can fire | 0.3 | 35.3 | 3.2 | 0.104 | 0.20 | 1.08 | 0.11 | 1.2 |
| **posterior** | **none can fire** | **0.5** | **8.1** | **1.4** | **0.045** | **0.28** | **1.10** | **0.14** | **1.2** |
| posterior | none can fire | 1.0 | 1.0 | 0.9 | 0.029 | 0.40 | 1.17 | 0.27 | 2.2 |

Three results. First, the posterior SE removes the need for a fallback rule. It is bounded by the prior SD, so a person
with no measures has score 0 and SE 1 and gets 1 + 2 floor(1/t) points per factor at the prior with unit spread, which
is standard quadrature with importance weight 1; everyone else gets points in proportion to posterior uncertainty.
Second, posterior SE with t = 0.5 costs 8 points per person, 4.5 percent of the full-quadrature time, and is as accurate
as the strict per-person fallback at 157 points (mean deviation 0.14 SE in both, maximum 1.2 against 1.1). With t = 0.3
it costs 35 points and is indistinguishable from full quadrature (0.20 / 1.08 against 0.23 / 1.10 is replication noise).
Third, every observation-only-SE rule is dominated: the T6 R pipeline spends 34 points for a mean deviation of 0.18.
One point per person (posterior SE, t = 1.0) is a Laplace-type plug-in that already beats the production setting at a
third of its bias and 3 percent of the full time, but 8 points at t = 0.5 is the better trade. Thresholds between the
two add little: t = 0.7 already puts 1.2 points on the average person and gives the same accuracy as t = 1.0 (mean
deviation 0.27 SE), because with the posterior SE a factor gets 3 points only when its SE exceeds t, and few SEs exceed
0.7. Raising t above 0.5 therefore saves at most a third of a time that is already 4.5 percent of full quadrature, and
doubles the deviation. Thresholds at or above 1 all give one point per person, the floor, since the posterior SE never
exceeds the prior SD. The floor sits at 2.9 percent of full time in R against 0.8 percent in legacy not because the
integration simplifies less in R (single-core, one Hessian evaluation falls from 22.9 s to 0.11 s in R and from 50.5 s
to 0.20 s in legacy, so R is faster on every row) but because each optimizer iteration pays about 0.12 s of doParallel
round-trip and serialization cost, which equals the whole compute per iteration at this size; eight workers do not
speed up the one-point fit at all (1.01 s on one core, 0.83 s on eight). That cost is fixed per iteration, so at HJV's
sample size the adaptive stage runs at the compute ratio.

## Latent types and the joint final stage (T9)

Production estimates education and wage models jointly in one group with 8 latent types: a type-probability logit
with an intercept and free factor loadings, type-specific intercepts in every outcome model (one per alternative for a
logit), none in the measurement system. T9 gives the mini model 3 types with the same layout, estimates both codes on
the same simulated data (N = 2,000, seed 1), and times 8 types.

| Check | Result |
|---|---|
| Stage 1 (typed factor model in R, untyped model on the typed data in legacy) | 96 parameters to 5e-9, SEs to 4e-7, log-likelihood to 2e-12 |
| T1 at the truth, 344 parameters, 3 types | Log-likelihood to 1.4e-10. Gradient over the 244 outcome parameters to 2e-13. Hessian: 59,524 of 59,536 outcome elements to 1e-6; the 12 that differ are legacy's (see below) |
| Joint final stage, 3 types, R truth start vs legacy truth start | 340 parameters to 7e-9, SEs to 4e-7: PASS |
| Joint final stage, R default start vs legacy default start | Same log-likelihood (1e-7). Legacy relabels types 2 and 3; after the swap every parameter agrees to 0.0e+00 |
| Time to convergence, 3 types, 8 cores | Legacy 561 s (28 Ipopt iterations), R 430 s (15 iterations) |
| One evaluation, 3 types, 1 core (likelihood / gradient / Hessian) | Legacy 5.9 / 9.1 / 204 s, R 3.2 / 6.7 / 256 s |
| One evaluation, 8 types, 465 parameters, 1 core | Legacy 19.6 / 29.8 / 1,146 s, R 10.9 / 23.4 / 1,383 s |

T9 took three attempts, and that is the main finding: the typed, two-stage path of the R package had never been run
against an external reference and carried four defects, all fixed in 1.8.0 and now covered by tests.

1. **Multinomial logit with types used one type intercept for every alternative.** The parameter layout defines
   `<comp>_c<k>_type_<t>_intercept` per alternative, but the C++ added a single value per type to every alternative's
   index and differentiated the block as if it held `n_types - 1` entries. On the mini model the typed log-likelihood
   at the truth was off by 108, entirely from the three logits (HS track, application, major switching); binary and
   linear components were right. HJV's final stage has three such logits.
2. **A factor model with `n_types > 1` and no typed component crashed** ("Fixed values size mismatch"), which is
   exactly the stage-1 call of a typed pipeline.
3. **Typed components after an untyped `previous_stage` were mis-set.** The stage-1 vector was reused without the
   type-probability block, the C++ fixed that block as if it belonged to stage 1, and a naming difference between two
   internal paths dropped the multinomial logit values on the way. The fit converged with code 0 to a point 32.6
   log-likelihood units below the optimum, with a residual standard deviation on its bound.
4. **`fix_type_intercepts(choice = )`** fixed the wrong slot for multinomial logit, and `initialize_parameters()`
   dropped a fixed per-alternative intercept from the vector.

Two legacy findings came out as well. First, the legacy Hessian is wrong on the cross terms between a type intercept
and a regressor in typed linear outcomes (12 elements here; a finite difference of the shared gradient confirms the R
values to all printed digits). The error is proportional to the score, so it vanishes at the optimum: evaluated at the
converged estimates the two Hessians agree to 3e-5 and the reported SEs to 3e-7, and both codes converge to the same
point. Away from the optimum the curvature is off (at the truth, -20 against -145 on one element), which can slow or
misdirect Ipopt's Newton steps. That is consistent with the convergence trouble legacy showed in production and with
the help that scaling the start values by 0.001 gave, although the mini model cannot prove the link. Second, `Est_outcomes` multiplies every value read from `init_par.txt` by `loadingMultiplier`
(0.001), and stage 1 cannot run with `ntypes > 1` when the type model sits in a skipped group.

## Recommendation

1. Move the re-estimation to the R package. The likelihood code is the same to machine precision and every converged
   estimate matches legacy to 1e-8. At equal cores the R package reaches convergence in about half the legacy wall time
   on the measurement system and in about the same time on the education stage. Build the package with optimization
   (`R CMD INSTALL`, or `pkgbuild::compile_dll(debug = FALSE)`); a `load_all()` build runs 3 to 11 times slower.
2. Use R's posterior-mode factor scores with the posterior SE: `estimate_factorscores_rcpp(..., include_prior = TRUE)`
   (`estimate_hjv.R` line 397 currently takes the default, observation-only SE). With `adapt_int_thresh = 0.5` and 8
   points per factor this uses 8 integration points per person on the mini model, 2.6 times the production count and
   a twentieth of full quadrature, and no fallback rule is involved. The legacy ML scores with their shrunk SEs are the
   source of the largest biases, so no ML-score option is needed in R.
3. Expect those adaptive estimates to sit 0.14 SE from full quadrature on average, with individual coefficients up to
   1.2 SE away. If the final stage can afford four times that integration cost, `adapt_int_thresh = 0.3` (35 points
   per person) reproduces full quadrature within replication noise.
4. Earlier legacy results that used adaptive integration carry the ML-score bias. On this mini model it reached 6.7 SE
   on a factor loading. The size on the real HJV estimates is unknown until they are re-estimated.

## Speed (T7)

Same model, same parameters, same quadrature (8 points per factor), both codes compiled at `-O2`, machine otherwise
idle. Legacy times are wall time of `mpirun` (Ipopt's own time in parentheses); R times are wall time of
`estimate_model_rcpp()`, which includes starting the workers.

| Measurement | Legacy | R | Ratio R / legacy |
|---|---|---|---|
| One likelihood evaluation, N = 5,000, 1 core | 4.0 s | 3.0 s | 0.74 |
| Likelihood + gradient | 6.1 s | 8.6 s | 1.4 |
| Likelihood + gradient + Hessian | 123 s | 160 s | 1.3 |
| Same three at N = 2,000 (medians of 3) | 1.59, 2.43, 49.2 s | 1.49, 3.03, 72.9 s | 0.94, 1.25, 1.48 |
| Stage 1 to convergence, N = 5,000, 8 cores | 131 s (Ipopt 123 s, 38 iterations) | 56.5 s (12 nlminb iterations) | 0.43 |
| Education stage to convergence, N = 5,000, 8 cores | Ipopt 79 s (118 s wall with the wage group) | 76 s | 0.96 |
| Stage 1, N = 2,000: 1, 4, 8 cores | 394, 96, 56 s (Ipopt 374, 89, 52) | 123, 38, 20 s | 0.31, 0.39, 0.36 |
| Adaptive education stage, N = 2,000, 8 cores | 3 s (Ipopt 1.1 s) | 3.9 s (optimizer 1.8 s) | 1.3 |

The two codes are within a factor 1.5 of each other per evaluation: R is faster on the likelihood alone and slower on
the gradient and Hessian, with the Hessian the largest gap (R's Hessian call costs 24 times its gradient call, legacy's
20 times). To convergence R is faster, since `nlminb` takes a third of Ipopt's iterations on stage 1 and each code
computes a Hessian per iteration. Per-evaluation times scale linearly in N on both sides. Parallel speedups from 1 to
8 cores are 7.0 for legacy (MPI) and 6.1 for R (doParallel). Peak memory could not be compared: `/usr/bin/time -l`
reports the launcher process only.

**Correction.** The first version of this table (R 3.2, 4.8 and 10.7 times slower per evaluation, 2.8 to 5.7 times
to convergence) timed the shared object that `pkgload::load_all()` builds. pkgload compiles through
`pkgbuild::compile_dll(debug = TRUE)`, whose flags are `-UNDEBUG -Wall -pedantic -g -O0`, and that build is 4.4, 5.0 and
11.3 times slower than `-O2` on the likelihood, gradient and Hessian. The February package (v1.0.2) and v1.4.0 installed
at `-O2` time the same as the current source, so there is no regression in the C++ code. The T2 to T6 results are
unaffected (estimates do not depend on the build); only their quoted run times were.

## Findings that need a decision

1. **Labor-stage likelihood.** Legacy estimates each outcome group with the measurement models only, ignoring earlier
   outcome groups. `estimate_hjv.R` Stage 4 keeps the education components in the likelihood. Both are consistent
   under the model; the R version uses more information and gives different wage estimates. The comparison above
   used the legacy convention on both sides.
2. **R package bug, fixed.** `define_model_system(previous_stage = )` stored `variance_fixed` and `variance_value` in
   the factor model, so a later `previous_stage` call failed the identity check; `estimate_hjv.R` Stage 4 hit this.
   The check now compares the specification only, `variance_value` covers every factor, and a regression test
   (`test-previous-stage-chain.R`) passes. Not yet committed.
3. **Factor scores and fallback rule: settled by T6 and T8.** Keep R's posterior-mode scores, and feed the adaptive
   integration the posterior SE (`include_prior = TRUE`) so that no fallback rule is needed.
4. **Latent types and the joint final stage: settled by T9, after four package fixes.** Use factorana 1.8.0 or later
   for any typed estimation; every earlier version mis-evaluates a typed multinomial logit with more than two
   alternatives and mis-sets a typed stage that follows an untyped stage. With 1.8.0 the typed joint final stage
   matches legacy to 7e-9 and converges faster. Type labels can permute between starts (legacy's default start
   relabeled types 2 and 3 here); compare typed fits up to that permutation.

## Legacy pitfalls met on the way

- With `stage = 1`, `AddModel` skips outcome models, but `LastModel_FixParamValue()` still binds to the last model
  added, which is then a measurement model. In the first run it silently fixed `stresstol_x2` at 0.
- Stage 0 (`TestCode`) evaluates measurement models with the likelihood-only flag whenever outcome models are
  present, so it does not return their derivatives.
- `PredictFactors` aborts under MPI with `MPI_ERR_TRUNCATE`; it runs serially (2 s for 2,000 people).

## Build on this Mac

Homebrew ROOT 6.40.04, Ipopt 3.14.20 with bundled MUMPS 5.6.2, Open MPI 5.0.11, Apple clang 17. Five patches,
listed in `legacy_build/PATCHES.md`: MUMPS in place of HSL ma57, two constraint callbacks returning `true` for a
problem with no constraints, C++20 and the ROOT 6 library list, and two diagnostic additions to stage 0.

## Files

`validation/hjv_mini/`: `generate_data.R`, `model_r.R`, `prob_hjv_mini.cc`, `run_r.R`, `run_legacy.sh`,
`compare_stage.R`, `t1_*.R`, `t4_*.R`, `t5_*.R`, `t6_*.R`, `run_t7.sh`, `NOTES.md` (all numbers), `results/`.
