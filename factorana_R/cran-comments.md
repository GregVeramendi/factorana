## Submission summary

This is a feature update (version 1.7.1). Since the version currently on CRAN
it adds, with documentation and tests:

* `simulate_factor_model()`: simulate data from a fitted or specified model.
* `vcov_factorana()` / `robust_se()`: robust and cluster-robust (sandwich)
  standard errors.
* A restartable, cluster-aware bootstrap
  (`generate_bootstrap_samples()`, `bootstrap_fit_sample()`,
  `collect_bootstrap()`, `bootstrap_factorana()`,
  `bootstrap_factorana_multistage()`).
* `SE_interactions` / `SE_full` structural-equation factor structures and
  `fix_factor_param()`.
* Clearer error messages for mis-coded ordered-probit / multinomial-logit
  outcomes.

## Test environments

* local: macOS (aarch64-apple-darwin20), R 4.5.2.
* win-builder: R-release and R-devel.

## R CMD check results

0 ERRORs, 0 WARNINGs.

There is 1 NOTE, only in the local check environment:

* "checking HTML version of manual ... NOTE": the local machine does not have a
  recent enough HTML Tidy and the V8 package is unavailable, so the HTML/maths
  validation is skipped. This NOTE does not appear on win-builder.

Finite-difference gradient/Hessian developer checks that compare analytical
derivatives to numerical ones are guarded with `skip_on_cran()`, because under
alternative BLAS implementations the finite-difference reference accumulates
floating-point round-off differently and can exceed the comparison tolerance on
a single near-zero matrix element. The analytical derivatives themselves are
correct; the package's user-facing functionality is fully tested.

## Reverse dependencies

There are no reverse dependencies.
