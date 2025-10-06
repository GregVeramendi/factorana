test_that("multi-factor loadings work with normalization (linear)", {
  set.seed(1)
  fm <- define_factor_model(3, 1, 8, loading_normalization = c(NA, 0, 1))
  dat <- data.frame(y = rnorm(120), x1 = rnorm(120), x2 = rnorm(120), eval = 1L)

  mc <- define_model_component("lin3", dat, "y", fm,
                               evaluation_indicator = "eval",
                               covariates = c("x1","x2"),
                               model_type = "linear")
  ini <- initialize_parameters(mc)

  expect_length(ini$loading, 3)
  expect_true(is.numeric(ini$loading))
  expect_equal(ini$loading[2], 0)
  expect_equal(ini$loading[3], 1)
})

test_that("multi-factor loadings work with normalization (probit)", {
  set.seed(2)
  fm <- define_factor_model(2, 1, 8, loading_normalization = c(NA, 1))
  # simple binary outcome
  x1 <- rnorm(150); x2 <- rnorm(150)
  p  <- plogis(0.5*x1 - 0.3*x2)
  y  <- as.integer(runif(150) < p)
  dat <- data.frame(y=y, x1=x1, x2=x2, eval=1L)

  mc <- define_model_component("prb2", dat, "y", fm,
                               evaluation_indicator = "eval",
                               covariates = c("x1","x2"),
                               model_type = "probit")
  ini <- initialize_parameters(mc)
  expect_length(ini$loading, 2)
  expect_equal(ini$loading[2], 1)
})
