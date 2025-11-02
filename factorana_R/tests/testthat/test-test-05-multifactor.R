test_that("multi-factor loadings work with normalization (linear)", {
  set.seed(1)
  fm <- define_factor_model(3, 1, 8)
  dat <- data.frame(y = rnorm(120), x1 = rnorm(120), x2 = rnorm(120), eval = 1L)

  # Specify loading normalization at component level: NA, 0, 1
  mc <- define_model_component("lin3", dat, "y", fm,
                               evaluation_indicator = "eval",
                               covariates = c("x1","x2"),
                               model_type = "linear",
                               loading_normalization = c(NA, 0, 1))

  # Create model system for initialization
  ms <- define_model_system(components = list(mc), factor = fm)
  ini <- initialize_parameters(ms, dat)

  # Check component loading normalization
  expect_length(mc$loading_normalization, 3)
  expect_true(is.na(mc$loading_normalization[1]))
  expect_equal(mc$loading_normalization[2], 0)
  expect_equal(mc$loading_normalization[3], 1)
})

test_that("multi-factor loadings work with normalization (probit)", {
  set.seed(2)
  fm <- define_factor_model(2, 1, 8)
  # simple binary outcome
  x1 <- rnorm(150); x2 <- rnorm(150)
  p  <- plogis(0.5*x1 - 0.3*x2)
  y  <- as.integer(runif(150) < p)
  dat <- data.frame(y=y, x1=x1, x2=x2, eval=1L)

  # Specify loading normalization at component level: NA, 1
  mc <- define_model_component("prb2", dat, "y", fm,
                               evaluation_indicator = "eval",
                               covariates = c("x1","x2"),
                               model_type = "probit",
                               loading_normalization = c(NA, 1))

  # Create model system for initialization
  ms <- define_model_system(components = list(mc), factor = fm)
  ini <- initialize_parameters(ms, dat)

  # Check component loading normalization
  expect_length(mc$loading_normalization, 2)
  expect_true(is.na(mc$loading_normalization[1]))
  expect_equal(mc$loading_normalization[2], 1)
})
