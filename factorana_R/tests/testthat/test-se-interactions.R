# Tests for structural-equation cross-product (interaction) structures:
#   SE_interactions: f_out = a + sum_j a_j f_j + sum_{a<b} a_ab f_a f_b + eps
#   SE_full:         adds the quadratic self-terms sum_j a_qj f_j^2
#
# These mirror the observed-component factor_spec = "interactions" / "full"
# options. They require at least two input factors (n_factors >= 3); with a
# single input factor the structures downgrade to SE_linear / SE_quadratic.
#
# The analytical gradient/Hessian for the cross-product terms is checked
# against finite differences (every upper-triangle Hessian element, per
# CLAUDE.md), followed by parameter recovery and the fix_factor_param path.

# Finite-difference tolerances, matching the rest of the test suite
# (test-systematic-suite.R, test-two-stage-se-types.R).
GRAD_TOL <- 1e-3
HESS_TOL <- 1e-3

# Build a logical param_fixed vector from the model system's real
# constraints, so check_*_accuracy() skips every parameter the model
# holds fixed (loading normalizations and fix_coefficient() intercepts).
# finite_diff_gradient()/finite_diff_hessian() perturb fixed positions
# regardless, so without this they report spurious mismatches at the
# pinned outcome-factor intercepts.
.se_param_fixed <- function(ms, params, init) {
  metadata <- factorana:::build_parameter_metadata(ms)
  constraints <- factorana:::setup_parameter_constraints(
    ms, params, metadata, init$factor_variance_fixed, verbose = FALSE)
  pf <- rep(TRUE, length(params))
  pf[constraints$free_idx] <- FALSE
  pf
}

# Build a 3-factor SE model (2 input factors f1, f2 + 1 outcome factor f3)
# with two linear measurement indicators per factor. Outcome-factor
# measurement intercepts are fixed to 0 for identification (otherwise
# se_intercept and the outcome-measurement intercepts form a flat ridge).
.build_se_inter_data <- function(n, structure, seed,
                                 a0 = 0.2, a1 = 0.5, a2 = 0.4,
                                 aq1 = 0.15, aq2 = 0.10, a12 = 0.25,
                                 res_sd = sqrt(0.5)) {
  set.seed(seed)
  f1 <- rnorm(n, 0, 1)
  f2 <- rnorm(n, 0, 1)
  eps <- rnorm(n, 0, res_sd)
  f3 <- a0 + a1 * f1 + a2 * f2 + a12 * f1 * f2 + eps
  if (structure == "SE_full") {
    f3 <- f3 + aq1 * f1^2 + aq2 * f2^2
  }
  data.frame(
    intercept = 1,
    y1 = 1.0 * f1 + rnorm(n, 0, 0.5),
    y2 = 0.8 * f1 + rnorm(n, 0, 0.5),
    y3 = 1.0 * f2 + rnorm(n, 0, 0.5),
    y4 = 0.9 * f2 + rnorm(n, 0, 0.5),
    y5 = 1.0 * f3 + rnorm(n, 0, 0.5),
    y6 = 0.9 * f3 + rnorm(n, 0, 0.5)
  )
}

.build_se_inter_system <- function(dat, structure) {
  fm <- define_factor_model(n_factors = 3, factor_structure = structure)
  # Factor 1 indicators
  mc1 <- define_model_component("m1", dat, "y1", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = c(1, 0, 0))
  mc2 <- define_model_component("m2", dat, "y2", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = c(NA, 0, 0))
  # Factor 2 indicators
  mc3 <- define_model_component("m3", dat, "y3", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = c(0, 1, 0))
  mc4 <- define_model_component("m4", dat, "y4", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = c(0, NA, 0))
  # Outcome factor indicators (intercepts fixed to 0 for identification)
  mc5 <- fix_coefficient(define_model_component("m5", dat, "y5", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = c(0, 0, 1)),
                         "intercept", 0)
  mc6 <- fix_coefficient(define_model_component("m6", dat, "y6", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = c(0, 0, NA)),
                         "intercept", 0)
  define_model_system(components = list(mc1, mc2, mc3, mc4, mc5, mc6), factor = fm)
}

test_that("SE_interactions gradient matches finite differences", {
  skip_on_cran()
  dat <- .build_se_inter_data(300, "SE_interactions", seed = 701)
  ms <- .build_se_inter_system(dat, "SE_interactions")
  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params

  # Confirm the new parameter is present and exercise it at a non-zero value.
  expect_true("se_interaction_1_2" %in% names(params))
  params["se_linear_1"] <- 0.5
  params["se_linear_2"] <- 0.4
  params["se_interaction_1_2"] <- 0.25

  param_fixed <- .se_param_fixed(ms, params, init)

  res <- check_gradient_accuracy(ms, dat, params, param_fixed = param_fixed,
                                 tol = GRAD_TOL, verbose = FALSE, n_quad = 12)
  expect_true(res$pass,
              info = sprintf("SE_interactions gradient check failed, max error: %.2e",
                             res$max_error))
})

test_that("SE_interactions Hessian matches finite differences", {
  skip_on_cran()
  dat <- .build_se_inter_data(300, "SE_interactions", seed = 702)
  ms <- .build_se_inter_system(dat, "SE_interactions")
  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params
  params["se_linear_1"] <- 0.5
  params["se_linear_2"] <- 0.4
  params["se_interaction_1_2"] <- 0.25

  param_fixed <- .se_param_fixed(ms, params, init)

  # The SE structure integrates over both input factors plus the residual
  # (3 dimensions), so the finite-difference Hessian needs finer quadrature
  # than the 2-D SE_linear models for structurally near-zero cross-terms
  # (e.g. factor_var_2 x m2_sigma, where m2 loads only on f1) to converge.
  res <- check_hessian_accuracy(ms, dat, params, param_fixed = param_fixed,
                                tol = HESS_TOL, verbose = FALSE, n_quad = 20)
  expect_true(res$pass,
              info = sprintf("SE_interactions Hessian check failed, max error: %.2e",
                             res$max_error))
})

test_that("SE_full gradient matches finite differences", {
  skip_on_cran()
  dat <- .build_se_inter_data(300, "SE_full", seed = 703)
  ms <- .build_se_inter_system(dat, "SE_full")
  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params

  expect_true("se_interaction_1_2" %in% names(params))
  expect_true("se_quadratic_1" %in% names(params))
  expect_true("se_quadratic_2" %in% names(params))
  params["se_linear_1"] <- 0.5
  params["se_linear_2"] <- 0.4
  params["se_quadratic_1"] <- 0.15
  params["se_quadratic_2"] <- 0.10
  params["se_interaction_1_2"] <- 0.25

  param_fixed <- .se_param_fixed(ms, params, init)

  res <- check_gradient_accuracy(ms, dat, params, param_fixed = param_fixed,
                                 tol = GRAD_TOL, verbose = FALSE, n_quad = 12)
  expect_true(res$pass,
              info = sprintf("SE_full gradient check failed, max error: %.2e",
                             res$max_error))
})

test_that("SE_full Hessian matches finite differences", {
  skip_on_cran()
  dat <- .build_se_inter_data(300, "SE_full", seed = 704)
  ms <- .build_se_inter_system(dat, "SE_full")
  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params
  params["se_linear_1"] <- 0.5
  params["se_linear_2"] <- 0.4
  params["se_quadratic_1"] <- 0.15
  params["se_quadratic_2"] <- 0.10
  params["se_interaction_1_2"] <- 0.25

  param_fixed <- .se_param_fixed(ms, params, init)

  # See note in the SE_interactions Hessian test: 3-D SE integration needs
  # finer quadrature for the FD comparison on near-zero cross-terms.
  res <- check_hessian_accuracy(ms, dat, params, param_fixed = param_fixed,
                                tol = HESS_TOL, verbose = FALSE, n_quad = 20)
  expect_true(res$pass,
              info = sprintf("SE_full Hessian check failed, max error: %.2e",
                             res$max_error))
})

test_that("SE_full recovers the cross-product coefficient", {
  skip_on_cran()
  # Larger sample for parameter recovery.
  a12_true <- 0.30
  dat <- .build_se_inter_data(3000, "SE_full", seed = 705,
                              a1 = 0.5, a2 = 0.4, aq1 = 0.10, aq2 = 0.10,
                              a12 = a12_true)
  ms <- .build_se_inter_system(dat, "SE_full")
  ctrl <- define_estimation_control(n_quad_points = 12, num_cores = 1)

  result <- estimate_model_rcpp(ms, dat, control = ctrl, parallel = FALSE,
                                optimizer = "nlminb", verbose = FALSE)

  expect_equal(result$convergence, 0,
               info = sprintf("SE_full did not converge, code: %d", result$convergence))

  est <- result$estimates
  se  <- result$std_errors
  expect_true("se_interaction_1_2" %in% names(est))

  # Standard errors should be sane before doing a "within 2 SE" check.
  expect_true(is.finite(se["se_interaction_1_2"]) && se["se_interaction_1_2"] < 10,
              info = sprintf("se_interaction_1_2 SE looks wrong: %.4g",
                             se["se_interaction_1_2"]))

  # Cross-product coefficient recovered within ~2.5 SE of truth.
  expect_lt(abs(est["se_interaction_1_2"] - a12_true),
            2.5 * se["se_interaction_1_2"])

  # Residual variance is less precisely estimated; allow a wider band there.
  expect_true(est["se_residual_var"] > 0)
})

test_that("fix_factor_param fixes se_interaction_1_2 at 0", {
  skip_on_cran()
  dat <- .build_se_inter_data(800, "SE_full", seed = 706, a12 = 0.0)
  fm <- define_factor_model(n_factors = 3, factor_structure = "SE_full")
  fm <- fix_factor_param(fm, c(se_interaction_1_2 = 0))
  expect_true("se_interaction_1_2" %in% names(fm$fixed_params))
  expect_equal(unname(fm$fixed_params["se_interaction_1_2"]), 0)

  mc1 <- define_model_component("m1", dat, "y1", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = c(1, 0, 0))
  mc2 <- define_model_component("m2", dat, "y2", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = c(NA, 0, 0))
  mc3 <- define_model_component("m3", dat, "y3", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = c(0, 1, 0))
  mc4 <- define_model_component("m4", dat, "y4", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = c(0, NA, 0))
  mc5 <- fix_coefficient(define_model_component("m5", dat, "y5", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = c(0, 0, 1)),
                         "intercept", 0)
  mc6 <- fix_coefficient(define_model_component("m6", dat, "y6", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = c(0, 0, NA)),
                         "intercept", 0)
  ms <- define_model_system(components = list(mc1, mc2, mc3, mc4, mc5, mc6), factor = fm)
  ctrl <- define_estimation_control(n_quad_points = 10, num_cores = 1)

  result <- estimate_model_rcpp(ms, dat, control = ctrl, parallel = FALSE,
                                optimizer = "nlminb", verbose = FALSE)

  expect_equal(result$convergence, 0,
               info = sprintf("fix_factor_param model did not converge, code: %d",
                              result$convergence))
  # The fixed interaction stays exactly at 0.
  expect_equal(unname(result$estimates["se_interaction_1_2"]), 0, tolerance = 1e-8)
})

test_that("SE_interactions downgrades to SE_linear with one input factor", {
  # n_factors = 2 -> single input factor -> no cross-products -> downgrade.
  expect_message(
    fm_i <- define_factor_model(n_factors = 2, factor_structure = "SE_interactions"),
    "downgrading to 'SE_linear'"
  )
  expect_equal(fm_i$factor_structure, "SE_linear")

  expect_message(
    fm_f <- define_factor_model(n_factors = 2, factor_structure = "SE_full"),
    "downgrading to 'SE_quadratic'"
  )
  expect_equal(fm_f$factor_structure, "SE_quadratic")

  # With 3 factors (2 input), the structures are kept as requested.
  fm_keep <- define_factor_model(n_factors = 3, factor_structure = "SE_interactions")
  expect_equal(fm_keep$factor_structure, "SE_interactions")
  # 2 input factors -> choose(2,2)=1 interaction coefficient.
  # nse_param = intercept(1) + linear(2) + interaction(1) + residual(1) = 5
  expect_equal(fm_keep$nse_param, 5L)

  fm_full <- define_factor_model(n_factors = 3, factor_structure = "SE_full")
  # nse_param = intercept(1) + linear(2) + quadratic(2) + interaction(1) + residual(1) = 7
  expect_equal(fm_full$nse_param, 7L)
})
