# Parallelization Test Suite for Factorana Package
#
# This test suite verifies that parallel estimation produces correct results
# and achieves meaningful speedup on larger datasets.

# Test configuration
VERBOSE <- Sys.getenv("FACTORANA_TEST_VERBOSE", "FALSE") == "TRUE"
SAVE_LOGS <- Sys.getenv("FACTORANA_TEST_SAVE_LOGS", "TRUE") == "TRUE"

# ==============================================================================
# Test: Roy Model Parallelization (n=10,000)
# ==============================================================================

test_that("Parallelization: Roy model with 1, 2, and 4 cores produces identical results", {
  skip_on_cran()

  set.seed(108)

  # Simulate Roy model data with n=10,000 for meaningful parallelization
  n <- 10000
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  f <- rnorm(n)  # Latent factor (ability)

  # True parameters
  true_params <- c(
    1.0,  # Factor variance
    # T1: int, sigma (loading FIXED to 1.0)
    2.0, 0.5,
    # T2: int, lambda, sigma
    1.5, 1.2, 0.6,
    # T3: int, lambda, sigma
    1.0, 0.8, 0.4,
    # Wage0: int, beta1, beta2, sigma (loading FIXED to 0.0)
    2.0, 0.5, 0.3, 0.6,
    # Wage1: int, beta1, lambda, sigma
    2.5, 0.6, 1.0, 0.7,
    # Sector: int, beta1, loading (probit)
    0.0, 0.4, 0.8
  )

  # Generate test scores
  T1 <- true_params[2] + 1.0*f + rnorm(n, 0, true_params[3])
  T2 <- true_params[4] + true_params[5]*f + rnorm(n, 0, true_params[6])
  T3 <- true_params[7] + true_params[8]*f + rnorm(n, 0, true_params[9])

  # Generate potential wages
  wage0 <- true_params[10] + true_params[11]*x1 + true_params[12]*x2 +
           rnorm(n, 0, true_params[13])
  wage1 <- true_params[14] + true_params[15]*x1 + true_params[16]*f +
           rnorm(n, 0, true_params[17])

  # Sector choice
  z_sector <- true_params[18] + true_params[19]*x2 + true_params[20]*f
  sector <- as.numeric(runif(n) < pnorm(z_sector))

  # Observed wage
  wage <- ifelse(sector == 1, wage1, wage0)

  dat <- data.frame(
    intercept = 1,
    x1 = x1, x2 = x2,
    T1 = T1, T2 = T2, T3 = T3,
    wage = wage,
    sector = sector,
    eval_tests = 1,
    eval_wage0 = 1 - sector,
    eval_wage1 = sector,
    eval_sector = 1
  )

  # Create model system
  fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 16)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval_tests"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval_tests"
  )
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval_tests"
  )
  mc_wage0 <- define_model_component(
    name = "wage0", data = dat, outcome = "wage", factor = fm,
    covariates = c("intercept", "x1", "x2"), model_type = "linear",
    loading_normalization = 0.0, evaluation_indicator = "eval_wage0"
  )
  mc_wage1 <- define_model_component(
    name = "wage1", data = dat, outcome = "wage", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval_wage1"
  )
  mc_sector <- define_model_component(
    name = "sector", data = dat, outcome = "sector", factor = fm,
    covariates = c("intercept", "x2"), model_type = "probit",
    loading_normalization = NA_real_, evaluation_indicator = "eval_sector"
  )

  ms <- define_model_system(
    components = list(mc_T1, mc_T2, mc_T3, mc_wage0, mc_wage1, mc_sector),
    factor = fm
  )

  if (VERBOSE) {
    cat("\n========================================\n")
    cat("Parallelization Test: Roy Model\n")
    cat("========================================\n")
    cat("Sample size: n =", n, "\n")
    cat("Parameters:", length(true_params), "\n\n")
  }

  # Test with 1, 2, and 4 cores
  test_cores <- c(1, 2, 4)
  results <- list()
  timings <- numeric(length(test_cores))

  for (i in seq_along(test_cores)) {
    nc <- test_cores[i]

    if (VERBOSE) cat(sprintf("Testing with %d core(s)...\n", nc))

    ctrl <- define_estimation_control(num_cores = nc)

    timings[i] <- system.time({
      results[[i]] <- estimate_model_rcpp(
        model_system = ms,
        data = dat,
        init_params = NULL,  # Use automatic initialization
        control = ctrl,
        parallel = (nc > 1),
        optimizer = "nloptr",
        verbose = FALSE
      )
    })["elapsed"]

    if (VERBOSE) {
      cat(sprintf("  Time: %.2f seconds\n", timings[i]))
      cat(sprintf("  Log-likelihood: %.4f\n", results[[i]]$loglik))
      cat(sprintf("  Convergence: %d\n\n", results[[i]]$convergence))
    }
  }

  # Verify all results match
  loglik_1core <- results[[1]]$loglik
  loglik_2core <- results[[2]]$loglik
  loglik_4core <- results[[3]]$loglik

  diff_2core <- abs(loglik_1core - loglik_2core)
  diff_4core <- abs(loglik_1core - loglik_4core)

  if (VERBOSE) {
    cat("========================================\n")
    cat("Results Verification\n")
    cat("========================================\n")
    cat(sprintf("1 core log-likelihood: %.6f\n", loglik_1core))
    cat(sprintf("2 core log-likelihood: %.6f\n", loglik_2core))
    cat(sprintf("4 core log-likelihood: %.6f\n", loglik_4core))
    cat(sprintf("Difference (1 vs 2 cores): %.2e\n", diff_2core))
    cat(sprintf("Difference (1 vs 4 cores): %.2e\n\n", diff_4core))
  }

  # Performance metrics
  speedup_2core <- timings[1] / timings[2]
  speedup_4core <- timings[1] / timings[3]
  efficiency_2core <- 100 * speedup_2core / 2
  efficiency_4core <- 100 * speedup_4core / 4

  if (VERBOSE) {
    cat("========================================\n")
    cat("Performance Metrics\n")
    cat("========================================\n")
    cat(sprintf("Time (1 core):  %.2f seconds\n", timings[1]))
    cat(sprintf("Time (2 cores): %.2f seconds (%.2fx speedup, %.1f%% efficiency)\n",
                timings[2], speedup_2core, efficiency_2core))
    cat(sprintf("Time (4 cores): %.2f seconds (%.2fx speedup, %.1f%% efficiency)\n\n",
                timings[3], speedup_4core, efficiency_4core))
  }

  # Collect diagnostics
  diagnostics <- list(
    sample_size = n,
    test_cores = test_cores,
    timings = timings,
    speedup_2core = speedup_2core,
    speedup_4core = speedup_4core,
    efficiency_2core = efficiency_2core,
    efficiency_4core = efficiency_4core,
    loglik_1core = loglik_1core,
    loglik_2core = loglik_2core,
    loglik_4core = loglik_4core,
    diff_2core = diff_2core,
    diff_4core = diff_4core
  )

  # Save log
  if (SAVE_LOGS) {
    log_file <- save_diagnostics_to_log("test_parallelization_roy", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(diff_2core < 1e-6,
              info = sprintf("2-core result differs from 1-core: %.2e", diff_2core))
  expect_true(diff_4core < 1e-6,
              info = sprintf("4-core result differs from 1-core: %.2e", diff_4core))

  # Expect at least some speedup with 4 cores (conservative threshold)
  # On some systems parallelization overhead may be high, so we only require >1.2x
  expect_true(speedup_4core > 1.2,
              info = sprintf("Expected speedup >1.2x with 4 cores, got %.2fx", speedup_4core))

  if (VERBOSE) cat("✓ All parallelization tests passed!\n\n")
})
