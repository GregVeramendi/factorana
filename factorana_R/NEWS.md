# factorana 1.1.3

## New tests and documentation

* New `tests/testthat/test-dynamic-single-factor.R` covers the standard
  workflow for estimating an SE_linear dynamic structural equation on a
  single latent construct measured at two time points:
  - Stage 1 fits a 2-factor independent measurement model with
    `equality_constraints` tying factor loadings and residual sigmas
    across periods but leaving measurement intercepts period-specific.
  - A dummy 2-factor `previous_stage` object is built that carries the
    wave-1 intercepts into both factor slots (discarding the wave-2
    intercepts, which absorb the structural-equation mean shift).
  - Stage 2 fits `SE_linear` and recovers the structural intercept
    `se_intercept` (alpha), slope `se_linear_1` (beta), residual
    variance `se_residual_var`, and input factor variance
    `factor_var_1`.
* New vignette `vignettes/dynamic_structural.Rmd` walks through the
  same workflow with executable code and explains the motivation for
  using wave-1 intercepts in Stage 2 (naive pooling of measurement
  intercepts across periods biases them by `lambda_m * E[f_2] / 2`,
  which propagates into an under-estimate of alpha).
* Removed the parameter-recovery test from
  `test-two-stage-se-types.R` (its Stage-1-no-types ->
  Stage-2-with-types setup is not a canonical workflow); the shape
  test, FD gradient/Hessian test, and the skipped
  Stage-1-with-types known-issue placeholder remain in place.

# factorana 1.1.2

## Bug fixes

* Two-stage estimation with `previous_stage + SE_linear/SE_quadratic` and
  `n_types > 1` now correctly builds the Stage 2 parameter vector. Prior
  versions omitted the `typeprob_*` and `type_*_loading_*` slots, causing
  either a crash in `setup_parameter_constraints()` or silently mis-fixed
  parameters. The measurement-parameter filter in
  `initialize_parameters()` was also tightened so that factor-level type
  parameters from Stage 1 are no longer duplicated into the measurement
  block. Discovered during analysis of a structural model where types are
  introduced at Stage 2.

## Tests

* New `tests/testthat/test-two-stage-se-types.R` adds:
  - a shape test that verifies Stage 2 SE_linear + `n_types = 2` produces
    the expected parameter vector aligned with `build_parameter_metadata()`,
  - a finite-difference gradient and Hessian check at the DGP parameters,
  - a structural parameter recovery test (`se_linear_1`,
    `se_intercept_type_2`, `se_residual_var`, `factor_var_1`) with init
    `se_intercept = -0.5` (Stage 1 absorbs E[f2] into the measurement
    intercepts, so the MLE `se_intercept` is negative even if the DGP
    constant is 0), and
  - a skipped placeholder documenting a known Hessian-FD mismatch in the
    Stage-1-with-types -> Stage-2-SE_linear variant (not the common
    workflow; tracked for a future fix).

# factorana 1.1.1

## CRAN resubmission fixes
* Added method references (with DOIs) to the `DESCRIPTION` field: Heckman,
  Humphries & Veramendi (2016, 2018) and Humphries, Joensen & Veramendi (2024).
* Added `\value` (via `@return`) to all exported functions that were missing it,
  including `as_kv`, `estimate_and_write`, `write_model_config_csv`, the
  adaptive-quadrature and observation-weight setters, and every `print` method.
* Replaced `\dontrun{}` with runnable examples (`fix_coefficient`,
  `fix_type_intercepts`) or `\donttest{}` blocks (`results_table`,
  `results_to_latex`, `components_table`, `estimate_factorscores_rcpp`,
  `cleanup_parallel_workers`).
* `estimate_and_write()` and `write_model_config_csv()` no longer write to a
  default path; `results_dir` / `file` is now required (use `tempdir()` in
  examples and tests).
* Replaced the non-executing `introduction` vignette with two executable
  vignettes: `measurement_system` (two-factor CFA) and `roy_model` (sector
  choice with a latent ability factor).

## Bug fixes
* Fixed systematic-suite test that constructed an ordered-probit component
  with `intercept` in its covariates; this configuration is rejected (the
  intercept is absorbed into the cut points).

# factorana 1.0.2

## Improvements
* CRAN compliance fixes
* Added introductory vignette
* Improved test coverage for SE models, equality constraints, observation weights
* Documentation improvements: fixed adaptive integration formula in README
* Updated SE_linear example to use larger sample size for reliable convergence

# factorana 1.0.1

## Bug Fixes
* Fix binary logit initialization and add dedicated test

# factorana 1.0.0

## New Features
* Structural equation models (SE_linear, SE_quadratic) for causal factor relationships
* Mixture of normals factor distribution (n_mixtures = 1, 2, or 3)
* Equality constraints for measurement invariance via `equality_constraints` parameter
* Component-level type control via `use_types` parameter
* Observation weights for survey weights or importance sampling
* Checkpointing for long-running estimations via `checkpoint_file` parameter
* Exploded multinomial logit for ranked choice models
* Exploded nested logit with `exclude_chosen = FALSE`
* Rank-share corrections via `rankshare_var` parameter

## Core Features
* Multi-factor models with flexible loading normalization
* Linear, probit, ordered probit, and multinomial logit model components
* Analytical gradients and Hessians for fast convergence
* Parallel estimation via doParallel for large datasets
* Multi-stage (sequential) estimation with fixed early-stage parameters
* Adaptive integration for efficient two-stage estimation
* Factor interaction terms (quadratic and cross-product) via `factor_spec`
* Correlated two-factor models via `factor_structure = "correlation"`
