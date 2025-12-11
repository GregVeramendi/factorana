# Tests for fix_type_intercepts function

test_that("fix_type_intercepts validates inputs correctly", {
  set.seed(123)
  n <- 100
  dat <- data.frame(intercept = 1, x1 = rnorm(n), Y = rnorm(n), eval = 1)

  # n_types = 1 should error (no type intercepts exist)
  fm1 <- define_factor_model(n_factors = 1, n_types = 1)
  mc1 <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm1,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )

  expect_error(
    fix_type_intercepts(mc1),
    "n_types < 2"
  )

  # n_types = 2 should work
  fm2 <- define_factor_model(n_factors = 1, n_types = 2)
  mc2 <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm2,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )

  mc2_fixed <- fix_type_intercepts(mc2)
  expect_s3_class(mc2_fixed, "model_component")
  expect_equal(length(mc2_fixed$fixed_type_intercepts), 1)
  expect_equal(mc2_fixed$fixed_type_intercepts[[1]]$type, 2)
  expect_equal(mc2_fixed$fixed_type_intercepts[[1]]$value, 0)

  # Error on non-model_component
  expect_error(
    fix_type_intercepts(list(a = 1)),
    "must be an object of class"
  )

  # Error on invalid type number
  expect_error(
    fix_type_intercepts(mc2, types = 1),  # Type 1 is reference
    "must be integers between 2"
  )
  expect_error(
    fix_type_intercepts(mc2, types = 5),  # Type 5 doesn't exist (n_types=2)
    "must be integers between 2"
  )
})

test_that("fix_type_intercepts correctly reduces parameter count", {
  set.seed(456)
  n <- 100
  dat <- data.frame(intercept = 1, x1 = rnorm(n), Y = rnorm(n), eval = 1)

  # 2-type model
  fm <- define_factor_model(n_factors = 1, n_types = 2)
  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval"
  )

  nparam_before <- mc$nparam_model
  mc_fixed <- fix_type_intercepts(mc)
  nparam_after <- mc_fixed$nparam_model

  # Should reduce by 1 (one type intercept for type 2)
  expect_equal(nparam_after, nparam_before - 1)
})

test_that("fix_type_intercepts excludes parameters from initialization", {
  set.seed(789)
  n <- 100
  dat <- data.frame(intercept = 1, x1 = rnorm(n), Y = rnorm(n), eval = 1)

  fm <- define_factor_model(n_factors = 1, n_types = 2)

  # Without fix_type_intercepts
  mc1 <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval"
  )
  ms1 <- define_model_system(components = list(mc1), factor = fm)
  init1 <- initialize_parameters(ms1, dat, verbose = FALSE)

  # With fix_type_intercepts
  mc2 <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval"
  )
  mc2 <- fix_type_intercepts(mc2)
  ms2 <- define_model_system(components = list(mc2), factor = fm)
  init2 <- initialize_parameters(ms2, dat, verbose = FALSE)

  # Type intercept should be in init1 but not in init2
  expect_true("Y_type_2_intercept" %in% init1$param_names)
  expect_false("Y_type_2_intercept" %in% init2$param_names)

  # Total parameters should differ by 1
  expect_equal(length(init2$init_params), length(init1$init_params) - 1)
})

test_that("fix_type_intercepts works with 3 types and partial fixing", {
  set.seed(101)
  n <- 100
  dat <- data.frame(intercept = 1, x1 = rnorm(n), Y = rnorm(n), eval = 1)

  fm <- define_factor_model(n_factors = 1, n_types = 3)
  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval"
  )

  # Fix only type 2 (leave type 3 free)
  mc_partial <- fix_type_intercepts(mc, types = 2)

  ms <- define_model_system(components = list(mc_partial), factor = fm)
  init <- initialize_parameters(ms, dat, verbose = FALSE)

  # Type 2 should be excluded, type 3 should be included
  expect_false("Y_type_2_intercept" %in% init$param_names)
  expect_true("Y_type_3_intercept" %in% init$param_names)

  # Fix all types (default behavior)
  mc_all <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval"
  )
  mc_all <- fix_type_intercepts(mc_all)  # Fixes types 2 and 3

  ms_all <- define_model_system(components = list(mc_all), factor = fm)
  init_all <- initialize_parameters(ms_all, dat, verbose = FALSE)

  # Both should be excluded
  expect_false("Y_type_2_intercept" %in% init_all$param_names)
  expect_false("Y_type_3_intercept" %in% init_all$param_names)
})

test_that("is_type_intercept_fixed helper works correctly", {
  set.seed(202)
  n <- 100
  dat <- data.frame(intercept = 1, x1 = rnorm(n), Y = rnorm(n), eval = 1)

  fm <- define_factor_model(n_factors = 1, n_types = 3)
  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval"
  )

  # Before fixing
  expect_false(is_type_intercept_fixed(mc, 2))
  expect_false(is_type_intercept_fixed(mc, 3))

  # Fix type 2 only
  mc <- fix_type_intercepts(mc, types = 2)

  expect_true(is_type_intercept_fixed(mc, 2))
  expect_false(is_type_intercept_fixed(mc, 3))
})
