#make_toycreates a Y that's continuous, not binary. the test should throw an error
#saying that it's supposed to be either 0 or 1. this enforces the model type constraint

test_that("probit requires 0/1 outcome", {
  dat <- make_toy()
  fm  <- define_factor_model(1, 1, 8)
  expect_error(
    define_model_component("sel", dat, "Y", fm,
                           covariates = "Z1", model_type = "probit"),
    regexp = "0/1"
  )
})

test_that("ordered probit returns J-1 ordered cutpoints", {
  skip_if_not_installed("MASS")
  dat <- make_toy()
  fm  <- define_factor_model(1, 1, 8)

  # make an ordinal outcome with 4 categories
  dat$Yord <- cut(dat$Y,
                  breaks = quantile(dat$Y, probs = seq(0, 1, by = 0.25)),
                  include.lowest = TRUE, ordered_result = TRUE)

  mc <- define_model_component("Y_ord", dat, "Yord", fm,
                               evaluation_indicator = "eval_y1",
                               covariates = "X1",
                               model_type = "oprobit")
  ini <- initialize_parameters(mc)

  expect_true(is.numeric(ini$thresholds))
  expect_length(ini$thresholds, 3)          # 4 cats => 3 thresholds
  expect_true(all(diff(ini$thresholds) > 0))# strictly increasing
})

test_that("oprobit works with multi-factor loading normalization", {
  set.seed(42)
  dat <- make_toy()

  # create an ordered outcome with 4 categories
  z <- dat$Y
  dat$Yord <- cut(z,
                  breaks = quantile(z, seq(0, 1, 0.25), na.rm = TRUE),
                  include.lowest = TRUE, ordered_result = TRUE)

  # two factors; fix the 2nd loading to 1, leave the 1st free (NA)
  fm <- define_factor_model(2, 1, 8, loading_normalization = c(NA, 1))

  mc <- define_model_component(
    "Y_ord2", dat, "Yord", fm,
    evaluation_indicator = "eval_y1",
    covariates = "X1",
    model_type = "oprobit"
  )
  ini <- initialize_parameters(mc)

  # outcome coerced to ordered, thresholds valid
  expect_true(is.ordered(mc$data[[mc$outcome]]))
  expect_true(is.numeric(ini$thresholds))
  expect_length(ini$thresholds, 3)
  expect_true(all(diff(ini$thresholds) > 0))

  # loadings: length k and the 2nd fixed to 1
  expect_length(ini$loading, 2)
  expect_true(is.numeric(ini$loading))
  expect_equal(ini$loading[2], 1)
})
