# Tests for structural equation models (SE_linear, SE_quadratic)

test_that("SE_linear gradient is accurate vs finite differences", {
  skip_on_cran()
  set.seed(123)
  n <- 300

  # Generate data with SE_linear structure: f2 = alpha + alpha1*f1 + epsilon
  f1 <- rnorm(n, 0, 1)
  epsilon <- rnorm(n, 0, sqrt(0.5))  # Residual variance = 0.5
  f2 <- 0.2 + 0.8 * f1 + epsilon  # SE structure

  # Measurement equations
  y1 <- 1.0 * f1 + rnorm(n, 0, 0.5)  # f1 measure 1 (loading = 1, fixed)
  y2 <- 0.8 * f1 + rnorm(n, 0, 0.5)  # f1 measure 2
  y3 <- 1.0 * f2 + rnorm(n, 0, 0.5)  # f2 measure 1 (loading = 1, fixed)
  y4 <- 0.9 * f2 + rnorm(n, 0, 0.5)  # f2 measure 2

  dat <- data.frame(y1 = y1, y2 = y2, y3 = y3, y4 = y4, intercept = 1)

  # Define SE_linear model
  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_linear")

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(1, 0))
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(NA, 0))
  mc3 <- define_model_component("m3", dat, "y3", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, 1))
  mc4 <- define_model_component("m4", dat, "y4", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, NA))

  ms <- define_model_system(components = list(mc1, mc2, mc3, mc4), factor = fm)
  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params

  # Determine which parameters are fixed
  param_fixed <- rep(FALSE, length(params))
  param_fixed[1] <- TRUE  # factor_var_1 is fixed when loading is fixed

  # Check gradient accuracy - use looser tolerance for SE models due to numerical precision
  result <- check_gradient_accuracy(ms, dat, params, param_fixed = param_fixed,
                                     tol = 5e-2, verbose = FALSE, n_quad = 8)
  expect_true(result$pass, info = sprintf("SE_linear gradient check failed, max error: %.2e", result$max_error))
})

test_that("SE_linear Hessian is accurate vs finite differences", {
  skip_on_cran()
  set.seed(124)
  n <- 300

  # Generate data with SE_linear structure
  f1 <- rnorm(n, 0, 1)
  epsilon <- rnorm(n, 0, sqrt(0.5))
  f2 <- 0.2 + 0.8 * f1 + epsilon

  y1 <- 1.0 * f1 + rnorm(n, 0, 0.5)
  y2 <- 0.8 * f1 + rnorm(n, 0, 0.5)
  y3 <- 1.0 * f2 + rnorm(n, 0, 0.5)
  y4 <- 0.9 * f2 + rnorm(n, 0, 0.5)

  dat <- data.frame(y1 = y1, y2 = y2, y3 = y3, y4 = y4, intercept = 1)

  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_linear")

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(1, 0))
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(NA, 0))
  mc3 <- define_model_component("m3", dat, "y3", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, 1))
  mc4 <- define_model_component("m4", dat, "y4", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, NA))

  ms <- define_model_system(components = list(mc1, mc2, mc3, mc4), factor = fm)
  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params

  param_fixed <- rep(FALSE, length(params))
  param_fixed[1] <- TRUE

  # Check Hessian accuracy - use looser tolerance for SE models
  result <- check_hessian_accuracy(ms, dat, params, param_fixed = param_fixed,
                                    tol = 5e-2, verbose = FALSE, n_quad = 8)
  expect_true(result$pass, info = sprintf("SE_linear Hessian check failed, max error: %.2e", result$max_error))
})

test_that("SE_quadratic gradient is accurate vs finite differences", {
  skip_on_cran()
  set.seed(125)
  n <- 300

  # Generate data with SE_quadratic structure: f2 = alpha + alpha1*f1 + alpha2*f1^2 + epsilon
  f1 <- rnorm(n, 0, 1)
  epsilon <- rnorm(n, 0, sqrt(0.5))
  f2 <- 0.2 + 0.6 * f1 + 0.1 * f1^2 + epsilon

  y1 <- 1.0 * f1 + rnorm(n, 0, 0.5)
  y2 <- 0.8 * f1 + rnorm(n, 0, 0.5)
  y3 <- 1.0 * f2 + rnorm(n, 0, 0.5)
  y4 <- 0.9 * f2 + rnorm(n, 0, 0.5)

  dat <- data.frame(y1 = y1, y2 = y2, y3 = y3, y4 = y4, intercept = 1)

  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_quadratic")

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(1, 0))
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(NA, 0))
  mc3 <- define_model_component("m3", dat, "y3", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, 1))
  mc4 <- define_model_component("m4", dat, "y4", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, NA))

  ms <- define_model_system(components = list(mc1, mc2, mc3, mc4), factor = fm)
  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params

  param_fixed <- rep(FALSE, length(params))
  param_fixed[1] <- TRUE

  result <- check_gradient_accuracy(ms, dat, params, param_fixed = param_fixed,
                                     tol = 5e-2, verbose = FALSE, n_quad = 8)
  expect_true(result$pass, info = sprintf("SE_quadratic gradient check failed, max error: %.2e", result$max_error))
})

test_that("SE_quadratic Hessian is accurate vs finite differences", {
  skip_on_cran()
  set.seed(126)
  n <- 300

  f1 <- rnorm(n, 0, 1)
  epsilon <- rnorm(n, 0, sqrt(0.5))
  f2 <- 0.2 + 0.6 * f1 + 0.1 * f1^2 + epsilon

  y1 <- 1.0 * f1 + rnorm(n, 0, 0.5)
  y2 <- 0.8 * f1 + rnorm(n, 0, 0.5)
  y3 <- 1.0 * f2 + rnorm(n, 0, 0.5)
  y4 <- 0.9 * f2 + rnorm(n, 0, 0.5)

  dat <- data.frame(y1 = y1, y2 = y2, y3 = y3, y4 = y4, intercept = 1)

  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_quadratic")

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(1, 0))
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(NA, 0))
  mc3 <- define_model_component("m3", dat, "y3", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, 1))
  mc4 <- define_model_component("m4", dat, "y4", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, NA))

  ms <- define_model_system(components = list(mc1, mc2, mc3, mc4), factor = fm)
  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params

  param_fixed <- rep(FALSE, length(params))
  param_fixed[1] <- TRUE

  result <- check_hessian_accuracy(ms, dat, params, param_fixed = param_fixed,
                                    tol = 5e-2, verbose = FALSE, n_quad = 8)
  expect_true(result$pass, info = sprintf("SE_quadratic Hessian check failed, max error: %.2e", result$max_error))
})

test_that("SE_linear model converges", {
  skip_on_cran()
  set.seed(128)
  n <- 500

  # Generate data with SE_linear structure
  f1 <- rnorm(n, 0, 1)
  epsilon <- rnorm(n, 0, sqrt(0.5))
  f2 <- 0.7 * f1 + epsilon

  y1 <- 1.0 * f1 + rnorm(n, 0, 0.5)
  y2 <- 0.8 * f1 + rnorm(n, 0, 0.5)
  y3 <- 1.0 * f2 + rnorm(n, 0, 0.5)
  y4 <- 0.9 * f2 + rnorm(n, 0, 0.5)

  dat <- data.frame(y1 = y1, y2 = y2, y3 = y3, y4 = y4, intercept = 1)

  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_linear")

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(1, 0))
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(NA, 0))
  mc3 <- define_model_component("m3", dat, "y3", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, 1))
  mc4 <- define_model_component("m4", dat, "y4", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, NA))

  ms <- define_model_system(components = list(mc1, mc2, mc3, mc4), factor = fm)
  ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)

  result <- estimate_model_rcpp(ms, dat, control = ctrl, parallel = FALSE,
                                 optimizer = "nlminb", verbose = FALSE)

  # Check convergence (0 = success, 1 = relative convergence also acceptable)
  expect_true(result$convergence %in% c(0, 1),
               info = sprintf("SE_linear model did not converge, code: %d", result$convergence))

  # Check SE parameters exist and are finite
  expect_true(is.finite(result$estimates["se_intercept"]))
  expect_true(is.finite(result$estimates["se_linear_1"]))
  expect_true(is.finite(result$estimates["se_residual_var"]))
  expect_true(result$estimates["se_residual_var"] > 0,
              info = "SE residual variance should be positive")
})
