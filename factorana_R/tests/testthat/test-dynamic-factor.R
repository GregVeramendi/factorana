# Tests for Dynamic Factor Model Components
#
# This test suite validates:
# 1. Dynamic factor model component creation
# 2. Parameter structure with outcome factor
# 3. Analytical gradients vs finite differences
# 4. Analytical Hessians vs finite differences
# 5. Model estimation with dynamic structural equations

# Test configuration
VERBOSE <- Sys.getenv("FACTORANA_TEST_VERBOSE", "FALSE") == "TRUE"
GRAD_TOL <- 1e-3
HESS_TOL <- 1e-3

# ==============================================================================
# Test 1: Basic dynamic model component structure
# ==============================================================================

test_that("define_dyn_model_component creates correct structure", {
  skip_on_cran()

  set.seed(301)
  n <- 100
  dat <- data.frame(
    intercept = 1,
    x = rnorm(n),
    Y = rnorm(n),
    eval = 1
  )

  # Create 2-factor model
  fm <- define_factor_model(n_factors = 2, n_types = 1)

  # Create dynamic component: f1 = intercept + lambda2*f2 + epsilon
  dyn_comp <- define_dyn_model_component(
    name = "structural",
    data = dat,
    outcome_factor = 1,
    factor = fm,
    covariates = NULL,
    intercept = TRUE,
    evaluation_indicator = "eval"
  )

  # Check structure
  expect_true(inherits(dyn_comp, "dyn_model_component"))
  expect_true(inherits(dyn_comp, "model_component"))
  expect_true(dyn_comp$is_dynamic)
  expect_equal(dyn_comp$outcome_factor, 1)
  expect_equal(dyn_comp$dyn_type, "linear")
  expect_equal(dyn_comp$k, 2)

  # Loading normalization: -1 for outcome factor, NA for others
  expect_equal(dyn_comp$loading_normalization[1], -1)
  expect_true(is.na(dyn_comp$loading_normalization[2]))

  # Check covariates include intercept
  expect_true("intercept" %in% dyn_comp$covariates)
})

# ==============================================================================
# Test 2: Dynamic model without intercept
# ==============================================================================

test_that("define_dyn_model_component works without intercept", {
  skip_on_cran()

  set.seed(302)
  n <- 100
  dat <- data.frame(
    x = rnorm(n),
    eval = 1
  )

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  # Create dynamic component without intercept: f1 = lambda2*f2 + epsilon
  dyn_comp <- define_dyn_model_component(
    name = "struct_no_int",
    data = dat,
    outcome_factor = 1,
    factor = fm,
    covariates = NULL,
    intercept = FALSE,
    evaluation_indicator = "eval"
  )

  # No covariates when intercept = FALSE and covariates = NULL
  expect_equal(length(dyn_comp$covariates), 0)
  expect_true(dyn_comp$is_dynamic)
})

# ==============================================================================
# Test 3: Dynamic model with covariates
# ==============================================================================

test_that("define_dyn_model_component works with covariates", {
  skip_on_cran()

  set.seed(303)
  n <- 100
  dat <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    eval = 1
  )

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  # Create dynamic component: f1 = intercept + x1*beta1 + x2*beta2 + lambda2*f2 + epsilon
  dyn_comp <- define_dyn_model_component(
    name = "struct_covs",
    data = dat,
    outcome_factor = 1,
    factor = fm,
    covariates = c("x1", "x2"),
    intercept = TRUE,
    evaluation_indicator = "eval"
  )

  # Check covariates: intercept + x1 + x2
  expect_equal(dyn_comp$covariates, c("intercept", "x1", "x2"))
  expect_true(dyn_comp$is_dynamic)
})

# ==============================================================================
# Test 4: Dynamic model with different outcome factor
# ==============================================================================

test_that("outcome_factor can be any valid factor index", {
  skip_on_cran()

  set.seed(304)
  n <- 100
  dat <- data.frame(
    intercept = 1,
    eval = 1
  )

  fm <- define_factor_model(n_factors = 3, n_types = 1)

  # Test outcome_factor = 2
  dyn_comp <- define_dyn_model_component(
    name = "struct_f2",
    data = dat,
    outcome_factor = 2,
    factor = fm,
    intercept = TRUE,
    evaluation_indicator = "eval"
  )

  expect_equal(dyn_comp$outcome_factor, 2)
  expect_equal(dyn_comp$loading_normalization[1], NA_real_)  # f1 loading is free
  expect_equal(dyn_comp$loading_normalization[2], -1)        # f2 is outcome
  expect_equal(dyn_comp$loading_normalization[3], NA_real_)  # f3 loading is free
})

# ==============================================================================
# Test 5: Validation - invalid outcome_factor
# ==============================================================================

test_that("invalid outcome_factor is rejected", {
  skip_on_cran()

  set.seed(305)
  n <- 100
  dat <- data.frame(intercept = 1, eval = 1)

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  # outcome_factor = 0 (invalid)
  expect_error(
    define_dyn_model_component("test", dat, outcome_factor = 0, factor = fm),
    "between 1 and"
  )

  # outcome_factor = 3 (exceeds n_factors)
  expect_error(
    define_dyn_model_component("test", dat, outcome_factor = 3, factor = fm),
    "between 1 and"
  )
})

# ==============================================================================
# Test 6: Dynamic model requires n_factors >= 2
# ==============================================================================

test_that("dynamic model requires n_factors >= 2", {
  skip_on_cran()

  set.seed(306)
  n <- 100
  dat <- data.frame(intercept = 1, eval = 1)

  fm_1fac <- define_factor_model(n_factors = 1, n_types = 1)

  # Should error: need at least 2 factors for dynamic model
  expect_error(
    define_dyn_model_component("test", dat, outcome_factor = 1, factor = fm_1fac),
    "n_factors >= 2"
  )
})

# ==============================================================================
# Test 7: Dynamic component in model system
# ==============================================================================

test_that("dynamic component works in model system with measurement equations", {
  skip_on_cran()

  set.seed(307)
  n <- 200

  # Simulate data for 2-factor model
  # f1 = true latent factor 1
  # f2 = true latent factor 2
  # Y1 measures f1
  # Y2 measures f2
  # Structural equation: f1 = alpha*f2 + epsilon

  true_alpha <- 0.7
  true_sigma_struct <- 0.5

  f2 <- rnorm(n)
  f1 <- true_alpha * f2 + rnorm(n, 0, true_sigma_struct)

  Y1 <- f1 + rnorm(n, 0, 0.3)  # Measure f1
  Y2 <- f2 + rnorm(n, 0, 0.3)  # Measure f2

  dat <- data.frame(
    intercept = 1,
    Y1 = Y1,
    Y2 = Y2,
    eval = 1
  )

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  # Measurement equation for Y1 (loads on factor 1)
  mc1 <- define_model_component(
    name = "Y1", data = dat, outcome = "Y1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1, 0),  # Fix loading on f1=1, f2=0
    evaluation_indicator = "eval"
  )

  # Measurement equation for Y2 (loads on factor 2)
  mc2 <- define_model_component(
    name = "Y2", data = dat, outcome = "Y2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, 1),  # Fix loading on f1=0, f2=1
    evaluation_indicator = "eval"
  )

  # Structural equation: f1 = beta0 + lambda2*f2 + epsilon
  dyn <- define_dyn_model_component(
    name = "structural",
    data = dat,
    outcome_factor = 1,
    factor = fm,
    intercept = TRUE,
    evaluation_indicator = "eval"
  )

  # Build model system
  ms <- define_model_system(components = list(mc1, mc2, dyn), factor = fm)

  # Check system structure
  expect_equal(length(ms$components), 3)
  expect_true(ms$components[[3]]$is_dynamic)

  # Initialize parameters
  init <- initialize_parameters(ms, dat)

  # Check that we have the right parameters
  expect_true("structural_loading_2" %in% init$param_names)  # Loading on f2
  expect_true("structural_sigma" %in% init$param_names)
  expect_false("structural_loading_1" %in% init$param_names)  # Outcome factor loading not estimated
})

# ==============================================================================
# Test 8: Parameter initialization for dynamic model
# ==============================================================================

test_that("parameter initialization works for dynamic model", {
  skip_on_cran()

  set.seed(308)
  n <- 100
  dat <- data.frame(
    intercept = 1,
    x = rnorm(n),
    eval = 1
  )

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  dyn <- define_dyn_model_component(
    name = "dyn",
    data = dat,
    outcome_factor = 1,
    factor = fm,
    covariates = "x",
    intercept = TRUE,
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(dyn), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Check initialization values are reasonable
  # For dynamic model, sigma should be initialized to 1.0 (not 0)
  sigma_idx <- which(init$param_names == "dyn_sigma")
  expect_equal(unname(init$init_params[sigma_idx]), 1.0)

  # Coefficients should be initialized to small values
  intercept_idx <- which(init$param_names == "dyn_intercept")
  x_idx <- which(init$param_names == "dyn_x")
  expect_true(abs(init$init_params[intercept_idx]) < 1)
  expect_true(abs(init$init_params[x_idx]) < 1)
})

# ==============================================================================
# Test 9: Dynamic model with factor_spec = "quadratic"
# ==============================================================================

test_that("dynamic model with factor_spec='quadratic' skips outcome factor", {
  skip_on_cran()

  set.seed(309)
  n <- 100
  dat <- data.frame(
    intercept = 1,
    eval = 1
  )

  # 3-factor model for quadratic terms
  fm <- define_factor_model(n_factors = 3, n_types = 1)

  # Dynamic: f1 = intercept + lambda2*f2 + lambda3*f3 + lambda_q2*f2^2 + lambda_q3*f3^2 + epsilon
  # (No f1^2 term since f1 is the outcome)
  dyn <- define_dyn_model_component(
    name = "dyn",
    data = dat,
    outcome_factor = 1,
    factor = fm,
    intercept = TRUE,
    factor_spec = "quadratic",
    evaluation_indicator = "eval"
  )

  # Should have 2 quadratic loadings (for f2 and f3, not f1)
  expect_equal(dyn$n_quadratic_loadings, 2)

  # Initialize and check parameter names
  ms <- define_model_system(components = list(dyn), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Should have quad loadings for f2 and f3, not f1
  expect_true("dyn_loading_quad_2" %in% init$param_names)
  expect_true("dyn_loading_quad_3" %in% init$param_names)
  expect_false("dyn_loading_quad_1" %in% init$param_names)
})

# ==============================================================================
# Test 10: Dynamic model with factor_spec = "interactions"
# ==============================================================================

test_that("dynamic model with factor_spec='interactions' skips outcome factor", {
  skip_on_cran()

  set.seed(310)
  n <- 100
  dat <- data.frame(
    intercept = 1,
    eval = 1
  )

  # 3-factor model for interaction terms
  fm <- define_factor_model(n_factors = 3, n_types = 1)

  # Dynamic: f1 = intercept + lambda2*f2 + lambda3*f3 + lambda_inter_23*f2*f3 + epsilon
  # (No f1*f2 or f1*f3 terms since f1 is the outcome)
  dyn <- define_dyn_model_component(
    name = "dyn",
    data = dat,
    outcome_factor = 1,
    factor = fm,
    intercept = TRUE,
    factor_spec = "interactions",
    evaluation_indicator = "eval"
  )

  # Should have 1 interaction loading (f2*f3 only, not f1*f2 or f1*f3)
  expect_equal(dyn$n_interaction_loadings, 1)

  # Initialize and check parameter names
  ms <- define_model_system(components = list(dyn), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Should have interaction loading for f2*f3 only
  expect_true("dyn_loading_inter_2_3" %in% init$param_names)
  expect_false("dyn_loading_inter_1_2" %in% init$param_names)
  expect_false("dyn_loading_inter_1_3" %in% init$param_names)
})

# ==============================================================================
# Test 11: Print method for dynamic model component
# ==============================================================================

test_that("print method works for dynamic model component", {
  skip_on_cran()

  set.seed(311)
  n <- 50
  dat <- data.frame(intercept = 1, eval = 1)

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  dyn <- define_dyn_model_component(
    name = "dyn",
    data = dat,
    outcome_factor = 1,
    factor = fm,
    intercept = TRUE,
    evaluation_indicator = "eval"
  )

  # Print should work without error
  output <- capture.output(print(dyn))
  expect_true(length(output) > 0)
  expect_true(any(grepl("Dynamic", output)))
  expect_true(any(grepl("Outcome Factor", output)))
})

# ==============================================================================
# Test 12: Gradient accuracy for basic dynamic model
# ==============================================================================

test_that("Gradient accuracy for basic dynamic model", {
  skip_on_cran()

  set.seed(312)
  n <- 200

  # Generate data for 2-factor dynamic model
  f2 <- rnorm(n)
  f1 <- 0.5 + 0.7 * f2 + rnorm(n, 0, 0.5)

  # Measurement equations
  Y1 <- f1 + rnorm(n, 0, 0.3)
  Y2 <- f2 + rnorm(n, 0, 0.3)

  dat <- data.frame(
    intercept = 1,
    Y1 = Y1,
    Y2 = Y2,
    eval = 1
  )

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  # Measurement equations
  mc1 <- define_model_component(
    name = "Y1", data = dat, outcome = "Y1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1, 0),
    evaluation_indicator = "eval"
  )

  mc2 <- define_model_component(
    name = "Y2", data = dat, outcome = "Y2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, 1),
    evaluation_indicator = "eval"
  )

  # Dynamic structural equation
  dyn <- define_dyn_model_component(
    name = "dyn",
    data = dat,
    outcome_factor = 1,
    factor = fm,
    intercept = TRUE,
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc1, mc2, dyn), factor = fm)

  # Initialize and check gradient
  init <- initialize_parameters(ms, dat)
  test_params <- init$init_params

  grad_check <- check_gradient_accuracy(
    model_system = ms,
    data = dat,
    params = test_params,
    tol = GRAD_TOL,
    verbose = VERBOSE,
    n_quad = 5
  )

  # Check gradient accuracy
  expect_true(grad_check$pass)
})
