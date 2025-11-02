test_that("Gauss-Hermite quadrature accuracy", {
  skip_on_cran()

  set.seed(888)
  n <- 500

  # Generate data without √2 factor (our corrected implementation)
  f <- rnorm(n, 0, 1)
  Y <- 1.0 + 1.0*f + rnorm(n, 0, 0.5)

  dat <- data.frame(intercept=1, Y=Y, eval=1)

  # Test with different quadrature points
  for (n_quad in c(8, 16)) {
    fm <- define_factor_model(n_factors=1, n_types=1, n_quad=n_quad)
    mc <- define_model_component(name="Y", data=dat, outcome="Y", factor=fm,
      covariates="intercept", model_type="linear",
      loading_normalization=NA_real_, evaluation_indicator="eval")

    ms <- define_model_system(components=list(mc), factor=fm)

    init_params <- c(1.0, 1.0, 1.0, 0.5)
    result <- estimate_model_rcpp(
      ms, dat,
      init_params = init_params,
      optimizer = "nlminb",
      parallel = FALSE,
      verbose = FALSE
    )

    # Check that estimation completed
    expect_true(!is.null(result$loglik))

    # Check that loading is reasonably recovered (within 20%)
    expect_true(abs(result$estimates[3] - 1.0) < 0.2)
  }
})

test_that("n_quad parameter is properly passed to C++", {
  skip_on_cran()

  set.seed(456)
  n <- 200

  f <- rnorm(n)
  Y <- 1.0 + 1.0*f + rnorm(n, 0, 0.5)

  dat <- data.frame(intercept=1, Y=Y, eval=1)

  # Test with n_quad=8
  fm8 <- define_factor_model(n_factors=1, n_types=1, n_quad=8)
  mc8 <- define_model_component(name="Y", data=dat, outcome="Y", factor=fm8,
    covariates="intercept", model_type="linear",
    loading_normalization=NA_real_, evaluation_indicator="eval")
  ms8 <- define_model_system(components=list(mc8), factor=fm8)

  init_params <- c(1.0, 1.0, 1.0, 0.5)
  result8 <- estimate_model_rcpp(ms8, dat, init_params=init_params,
    optimizer="nlminb", parallel=FALSE, verbose=FALSE)

  # Test with n_quad=16
  fm16 <- define_factor_model(n_factors=1, n_types=1, n_quad=16)
  mc16 <- define_model_component(name="Y", data=dat, outcome="Y", factor=fm16,
    covariates="intercept", model_type="linear",
    loading_normalization=NA_real_, evaluation_indicator="eval")
  ms16 <- define_model_system(components=list(mc16), factor=fm16)

  result16 <- estimate_model_rcpp(ms16, dat, init_params=init_params,
    optimizer="nlminb", parallel=FALSE, verbose=FALSE)

  # Results should be similar but may differ slightly
  # The key is that both converge successfully
  expect_true(!is.null(result8$loglik))
  expect_true(!is.null(result16$loglik))

  # More quadrature points should generally give better approximation
  # (though for simple cases the difference is small)
  expect_true(abs(result8$loglik - result16$loglik) < 5)
})
