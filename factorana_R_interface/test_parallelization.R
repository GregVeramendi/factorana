#!/usr/bin/env Rscript
# Test parallelization in factorana
# Compares single-core vs multi-core performance

library(factorana)

cat("========================================\n")
cat("Parallelization Test\n")
cat("========================================\n\n")

# Generate data (measurement system + multinomial logit like Test D)
set.seed(107)
n <- 5000  # Larger sample for meaningful timing comparison
x1 <- rnorm(n)
f <- rnorm(n)

true_params <- c(1.0,  # Factor variance
                 2.0, 0.5,   # T1: int, sigma
                 1.5, 1.2, 0.6,   # T2: int, lambda, sigma
                 1.0, 0.8, 0.4,   # T3: int, lambda, sigma
                 0.5, 0.6, 0.7,   # Mlogit choice 1
                 1.0, -0.5, 0.9)  # Mlogit choice 2

T1 <- true_params[2] + 1.0*f + rnorm(n, 0, true_params[3])
T2 <- true_params[4] + true_params[5]*f + rnorm(n, 0, true_params[6])
T3 <- true_params[7] + true_params[8]*f + rnorm(n, 0, true_params[9])

z1 <- true_params[10] + true_params[11]*x1 + true_params[12]*f
z2 <- true_params[13] + true_params[14]*x1 + true_params[15]*f

exp_z0 <- 1
exp_z1 <- exp(z1)
exp_z2 <- exp(z2)
denom <- exp_z0 + exp_z1 + exp_z2

p0 <- exp_z0 / denom
p1 <- exp_z1 / denom
p2 <- exp_z2 / denom

y <- numeric(n)
for (i in seq_len(n)) {
  y[i] <- sample(1:3, 1, prob = c(p0[i], p1[i], p2[i]))
}

dat <- data.frame(intercept = 1, x1 = x1, T1 = T1, T2 = T2, T3 = T3, y = y, eval = 1)

# Create model system
fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8)

mc_T1 <- define_model_component(
  name = "T1", data = dat, outcome = "T1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = 1.0, evaluation_indicator = "eval"
)
mc_T2 <- define_model_component(
  name = "T2", data = dat, outcome = "T2", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval"
)
mc_T3 <- define_model_component(
  name = "T3", data = dat, outcome = "T3", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval"
)
mc_y <- define_model_component(
  name = "y", data = dat, outcome = "y", factor = fm,
  covariates = c("intercept", "x1"), model_type = "logit",
  num_choices = 3,
  loading_normalization = NA_real_,
  evaluation_indicator = "eval"
)

ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3, mc_y), factor = fm)

cat("Model: Measurement system (3 tests) + Multinomial logit\n")
cat("Sample size: n =", n, "\n")
cat("Parameters: 15\n\n")

# Test 1: Single core
cat("========================================\n")
cat("Test 1: Single core (num_cores = 1)\n")
cat("========================================\n\n")

ctrl_single <- define_estimation_control(num_cores = 1)

time_single <- system.time({
  result_single <- estimate_model_rcpp(
    model_system = ms,
    data = dat,
    init_params = true_params,
    control = ctrl_single,
    parallel = FALSE,
    optimizer = "nloptr",
    verbose = FALSE
  )
})

cat("Time elapsed:", round(time_single["elapsed"], 2), "seconds\n")
cat("Log-likelihood:", round(result_single$loglik, 4), "\n")
cat("Convergence:", result_single$convergence, "\n\n")

# Test 2: Multi-core
n_cores <- min(4, parallel::detectCores() - 1)  # Use up to 4 cores
if (n_cores < 2) n_cores <- 2

cat("========================================\n")
cat("Test 2: Multi-core (num_cores =", n_cores, ")\n")
cat("========================================\n\n")

ctrl_multi <- define_estimation_control(num_cores = n_cores)

time_multi <- system.time({
  result_multi <- estimate_model_rcpp(
    model_system = ms,
    data = dat,
    init_params = true_params,
    control = ctrl_multi,
    parallel = TRUE,
    optimizer = "nloptr",
    verbose = FALSE
  )
})

cat("Time elapsed:", round(time_multi["elapsed"], 2), "seconds\n")
cat("Log-likelihood:", round(result_multi$loglik, 4), "\n")
cat("Convergence:", result_multi$convergence, "\n\n")

# Comparison
cat("========================================\n")
cat("Comparison\n")
cat("========================================\n\n")

speedup <- time_single["elapsed"] / time_multi["elapsed"]
cat("Speedup:", round(speedup, 2), "x\n")
cat("Efficiency:", round(100 * speedup / n_cores, 1), "%\n\n")

# Check results match
loglik_diff <- abs(result_single$loglik - result_multi$loglik)
cat("Log-likelihood difference:", sprintf("%.2e", loglik_diff), "\n")

if (loglik_diff < 1e-6) {
  cat("✓ Results match (parallelization correct)\n")
} else {
  cat("✗ Results differ (parallelization bug!)\n")
}

cat("\n")
cat("Test completed successfully!\n")
