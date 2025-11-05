test_that("multi-factor loadings work with normalization (linear)", {
  set.seed(1)
  fm <- define_factor_model(3, 1)
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
  fm <- define_factor_model(2, 1)
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

test_that("two-factor CFA with ordered probit converges and recovers parameters", {
  skip_if_not_installed("MASS")

  # Generate data with known two-factor structure
  set.seed(123)
  n <- 500

  # Two independent factors
  f1 <- rnorm(n, 0, 1)
  f2 <- rnorm(n, 0, 1)

  # True loadings for factor 1 measures (m1, m2, m3)
  lambda1_1 <- 1.0   # Fixed for identification
  lambda1_2 <- 0.8
  lambda1_3 <- 1.2

  # True loadings for factor 2 measures (m4, m5, m6)
  lambda2_4 <- 1.0   # Fixed for identification
  lambda2_5 <- 0.9
  lambda2_6 <- 1.1

  # Generate latent continuous variables
  y1_star <- lambda1_1 * f1 + rnorm(n, 0, 0.5)
  y2_star <- lambda1_2 * f1 + rnorm(n, 0, 0.5)
  y3_star <- lambda1_3 * f1 + rnorm(n, 0, 0.5)
  y4_star <- lambda2_4 * f2 + rnorm(n, 0, 0.5)
  y5_star <- lambda2_5 * f2 + rnorm(n, 0, 0.5)
  y6_star <- lambda2_6 * f2 + rnorm(n, 0, 0.5)

  # Convert to ordered categories (5 categories each)
  # Use integers 1-5, define_model_component will convert to ordered factor
  make_ordered <- function(y_star) {
    as.integer(cut(y_star,
                   breaks = c(-Inf, -1, -0.3, 0.3, 1, Inf),
                   labels = FALSE))
  }

  dat <- data.frame(
    m1 = make_ordered(y1_star),
    m2 = make_ordered(y2_star),
    m3 = make_ordered(y3_star),
    m4 = make_ordered(y4_star),
    m5 = make_ordered(y5_star),
    m6 = make_ordered(y6_star)
  )

  # Define 2-factor model
  fm <- define_factor_model(n_factors = 2, n_types = 1)

  # Factor 1 measures: loading on f1, zero on f2
  # m1: fixed loading = 1.0 on f1 for identification
  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "m1", factor = fm,
    covariates = NULL,
    model_type = "oprobit",
    loading_normalization = c(1.0, 0),  # f1=1.0 (fixed), f2=0 (zero)
    num_choices = 5
  )

  # m2: free loading on f1
  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "m2", factor = fm,
    covariates = NULL,
    model_type = "oprobit",
    loading_normalization = c(NA, 0),  # f1=free, f2=0 (zero)
    num_choices = 5
  )

  # m3: free loading on f1
  mc3 <- define_model_component(
    name = "m3", data = dat, outcome = "m3", factor = fm,
    covariates = NULL,
    model_type = "oprobit",
    loading_normalization = c(NA, 0),  # f1=free, f2=0 (zero)
    num_choices = 5
  )

  # Factor 2 measures: zero on f1, loading on f2
  # m4: fixed loading = 1.0 on f2 for identification
  mc4 <- define_model_component(
    name = "m4", data = dat, outcome = "m4", factor = fm,
    covariates = NULL,
    model_type = "oprobit",
    loading_normalization = c(0, 1.0),  # f1=0 (zero), f2=1.0 (fixed)
    num_choices = 5
  )

  # m5: free loading on f2
  mc5 <- define_model_component(
    name = "m5", data = dat, outcome = "m5", factor = fm,
    covariates = NULL,
    model_type = "oprobit",
    loading_normalization = c(0, NA),  # f1=0 (zero), f2=free
    num_choices = 5
  )

  # m6: free loading on f2
  mc6 <- define_model_component(
    name = "m6", data = dat, outcome = "m6", factor = fm,
    covariates = NULL,
    model_type = "oprobit",
    loading_normalization = c(0, NA),  # f1=0 (zero), f2=free
    num_choices = 5
  )

  # Create model system
  ms <- define_model_system(
    components = list(mc1, mc2, mc3, mc4, mc5, mc6),
    factor = fm
  )

  # Initialize parameters
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)

  # Estimate the model
  ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)
  result <- estimate_model_rcpp(
    model_system = ms,
    data = dat,
    init_params = init_result$init_params,
    control = ctrl,
    parallel = FALSE,
    optimizer = "nlminb",
    verbose = FALSE
  )

  # Check convergence
  expect_equal(result$convergence, 0,
               info = "Model should converge successfully")

  # Check that we have reasonable estimates
  expect_true(is.numeric(result$estimates))
  expect_true(all(is.finite(result$estimates)))

  # Check standard errors are computed
  expect_true(is.numeric(result$std_errors))
  expect_true(all(is.finite(result$std_errors)))
  expect_true(all(result$std_errors > 0))

  # Extract estimated loadings
  # Note: Parameter order depends on implementation, but we can check structure
  param_names <- names(result$estimates)

  # Check that factor variances are estimated
  expect_true(any(grepl("sigma_f", param_names, ignore.case = TRUE)),
              info = "Should have factor variance parameters")

  # Check that loadings are estimated for the correct measures
  expect_true(any(grepl("m2.*Factor", param_names)),
              info = "m2 should have a free loading on f1")
  expect_true(any(grepl("m3.*Factor", param_names)),
              info = "m3 should have a free loading on f1")
  expect_true(any(grepl("m5.*Factor", param_names)),
              info = "m5 should have a free loading on f2")
  expect_true(any(grepl("m6.*Factor", param_names)),
              info = "m6 should have a free loading on f2")

  # Print results for manual inspection
  cat("\n=== Two-Factor CFA Results ===\n")
  cat("Convergence:", result$convergence, "\n")
  cat("Log-likelihood:", result$loglik, "\n")
  cat("Number of parameters:", length(result$estimates), "\n\n")

  if (!is.null(param_names) && length(param_names) > 0) {
    # Find loading parameters
    loading_idx <- grep("Factor", param_names)
    if (length(loading_idx) > 0) {
      cat("Estimated Loadings:\n")
      for (i in loading_idx) {
        cat(sprintf("  %s: %.4f (SE: %.4f)\n",
                    param_names[i],
                    result$estimates[i],
                    result$std_errors[i]))
      }
    }
  }

  # Parameter recovery check: loadings should be reasonably close to true values
  # We'll do a rough check on the free loadings
  # (exact recovery depends on discretization and identification)
  if (!is.null(param_names)) {
    # Check m2 loading (true = 0.8)
    m2_load_idx <- grep("m2.*Factor", param_names)
    if (length(m2_load_idx) > 0) {
      m2_est <- result$estimates[m2_load_idx[1]]
      expect_true(abs(m2_est - 0.8) < 0.3,
                  info = sprintf("m2 loading should be near 0.8, got %.3f", m2_est))
    }

    # Check m3 loading (true = 1.2)
    m3_load_idx <- grep("m3.*Factor", param_names)
    if (length(m3_load_idx) > 0) {
      m3_est <- result$estimates[m3_load_idx[1]]
      expect_true(abs(m3_est - 1.2) < 0.3,
                  info = sprintf("m3 loading should be near 1.2, got %.3f", m3_est))
    }

    # Check m5 loading (true = 0.9)
    m5_load_idx <- grep("m5.*Factor", param_names)
    if (length(m5_load_idx) > 0) {
      m5_est <- result$estimates[m5_load_idx[1]]
      expect_true(abs(m5_est - 0.9) < 0.3,
                  info = sprintf("m5 loading should be near 0.9, got %.3f", m5_est))
    }

    # Check m6 loading (true = 1.1)
    m6_load_idx <- grep("m6.*Factor", param_names)
    if (length(m6_load_idx) > 0) {
      m6_est <- result$estimates[m6_load_idx[1]]
      expect_true(abs(m6_est - 1.1) < 0.3,
                  info = sprintf("m6 loading should be near 1.1, got %.3f", m6_est))
    }
  }
})
