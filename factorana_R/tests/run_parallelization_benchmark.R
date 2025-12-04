#!/usr/bin/env Rscript
# Parallelization Benchmark
#
# Tests parallelization speedup using the INSTALLED package.
# Install the package first: R CMD INSTALL --no-multiarch .
#
# Run with: Rscript tests/run_parallelization_benchmark.R

cat("=== Parallelization Benchmark ===\n\n")

# Load installed package
if (!requireNamespace("factorana", quietly = TRUE)) {
  stop("factorana package not installed. Run: R CMD INSTALL --no-multiarch .")
}
library(factorana)

cat("Using factorana from:",
    dirname(system.file("DESCRIPTION", package = "factorana")), "\n\n")

set.seed(108)

# Simulate Roy model data
n <- 10000
x1 <- rnorm(n)
x2 <- rnorm(n)
f <- rnorm(n)

true_params <- c(
  1.0, 2.0, 0.5, 1.5, 1.2, 0.6, 1.0, 0.8, 0.4,
  2.0, 0.5, 0.3, 0.6, 2.5, 0.6, 1.0, 0.7, 0.0, 0.4, 0.8
)

T1 <- true_params[2] + 1.0*f + rnorm(n, 0, true_params[3])
T2 <- true_params[4] + true_params[5]*f + rnorm(n, 0, true_params[6])
T3 <- true_params[7] + true_params[8]*f + rnorm(n, 0, true_params[9])
wage0 <- true_params[10] + true_params[11]*x1 + true_params[12]*x2 + rnorm(n, 0, true_params[13])
wage1 <- true_params[14] + true_params[15]*x1 + true_params[16]*f + rnorm(n, 0, true_params[17])
z_sector <- true_params[18] + true_params[19]*x2 + true_params[20]*f
sector <- as.numeric(runif(n) < pnorm(z_sector))
wage <- ifelse(sector == 1, wage1, wage0)

dat <- data.frame(
  intercept = 1, x1 = x1, x2 = x2,
  T1 = T1, T2 = T2, T3 = T3,
  wage = wage, sector = sector,
  eval_tests = 1, eval_wage0 = 1 - sector,
  eval_wage1 = sector, eval_sector = 1
)

# Define model
fm <- define_factor_model(n_factors = 1, n_types = 1)

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

cat("Model: Roy model with 6 components, 20 parameters\n")
cat("Sample size: n =", n, "\n")
cat("Available cores:", parallel::detectCores(), "\n\n")

# Warm-up
cat("Warm-up run...\n")
ctrl_warmup <- define_estimation_control(n_quad_points = 8, num_cores = 1)
invisible(capture.output(suppressMessages(estimate_model_rcpp(
  model_system = ms, data = dat, init_params = NULL,
  control = ctrl_warmup, parallel = FALSE, optimizer = "nlminb", verbose = FALSE
))))

# Benchmark
test_cores <- c(1, 2, 4)
test_cores <- test_cores[test_cores <= parallel::detectCores()]

results <- list()
timings <- numeric(length(test_cores))

cat("\n--- Running benchmarks ---\n\n")

for (i in seq_along(test_cores)) {
  nc <- test_cores[i]
  cat(sprintf("Testing with %d core(s)... ", nc))

  ctrl <- define_estimation_control(n_quad_points = 16, num_cores = nc)

  timings[i] <- system.time({
    invisible(capture.output(suppressMessages(
      results[[i]] <- estimate_model_rcpp(
        model_system = ms, data = dat, init_params = NULL,
        control = ctrl, parallel = (nc > 1),
        optimizer = "nlminb", verbose = FALSE
      )
    )))
  })[["elapsed"]]

  cat(sprintf("%.2f sec (loglik: %.4f)\n", timings[i], results[[i]]$loglik))
}

# Results
cat("\n========================================\n")
cat("RESULTS SUMMARY\n")
cat("========================================\n\n")

cat("Timings:\n")
for (i in seq_along(test_cores)) {
  cat(sprintf("  %d core(s): %.2f sec\n", test_cores[i], timings[i]))
}

if (length(test_cores) > 1) {
  cat("\nSpeedups (relative to 1-core):\n")
  for (i in 2:length(test_cores)) {
    speedup <- timings[1] / timings[i]
    efficiency <- 100 * speedup / test_cores[i]
    cat(sprintf("  %d cores: %.2fx speedup (%.1f%% efficiency)\n",
                test_cores[i], speedup, efficiency))
  }

  cat("\nLog-likelihood consistency:\n")
  for (i in 2:length(test_cores)) {
    diff <- abs(results[[1]]$loglik - results[[i]]$loglik)
    status <- if (diff < 1e-5) "✓ PASS" else "✗ FAIL"
    cat(sprintf("  1-core vs %d-core: diff = %.2e %s\n", test_cores[i], diff, status))
  }
}

cat("\n========================================\n")
