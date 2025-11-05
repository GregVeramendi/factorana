#Validates the inputs: Valid outcome, covariate, and eval_indicator type

test_that("define_model_component validates inputs", {
  dat <- make_toy()
  fm  <- define_factor_model(1, 1)

  # good path : calls functino with valid arguments
  mc <- define_model_component("Y1", dat, "Y", fm,
                               evaluation_indicator = "eval_y1",
                               covariates = c("X1"),
                               model_type = "linear")
  expect_s3_class(mc, "model_component")

  # missing outcome ("NOPE" outcome name, doesn't exist)
  expect_error(
    define_model_component("bad", dat, "NOPE", fm,
                           covariates = "X1", model_type = "linear"),
    regexp = "outcome must be a column"
  )

  # missing covariate ("NOPE covariate doesn't exist)
  expect_error(
    define_model_component("bad", dat, "Y", fm,
                           covariates = "NOPE", model_type = "linear"),
    regexp = "Covariates not found"
  )

  # bad eval indicator type (eval_bad is a string instead of bool)
  dat2 <- dat; dat2$eval_bad <- "oops"
  expect_error(
    define_model_component("bad", dat2, "Y", fm,
                           evaluation_indicator = "eval_bad",
                           covariates = "X1", model_type = "linear"),
    regexp = "evaluation_indicator"
  )
})
