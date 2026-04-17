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
