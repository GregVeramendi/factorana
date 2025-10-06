#assert:intercept is a single number, betas have length 1 (only passed 1 covar)
#loading is one number
#no cutpoints for non ordered probit

test_that("initialize_parameters returns expected shapes", {
  dat <- make_toy()
  fm  <- define_factor_model(1, 1, 8)

  mc <- define_model_component("Y1", dat, "Y", fm,
                               evaluation_indicator = "eval_y1",
                               covariates = c("X1"),
                               model_type = "linear")
  ini <- initialize_parameters(mc)

  expect_true(is.numeric(ini$intercept) && length(ini$intercept) == 1)
  expect_length(ini$betas, 1)
  expect_true(is.numeric(ini$loading) && length(ini$loading) == 1)
  # non-oprobit -> cutpoints absent or empty
  expect_true(is.null(ini$cutpoints) || length(ini$cutpoints) == 0)
})
