#!/usr/bin/env Rscript
# Debug gradient indexing issue for multinomial logit with factors

library(factorana)

set.seed(42)
n <- 150

# Simple data
f1 <- rnorm(n)
x1 <- rnorm(n)
intercept <- rep(1, n)

U1 <- 0.5*intercept + 0.3*x1 + 0.4*f1 + rnorm(n)
U2 <- -0.3*intercept + 0.5*x1 + 0.6*f1 + rnorm(n)

Y_choice <- rep(1, n)
Y_choice[U1 > 0 & U1 > U2] <- 2
Y_choice[U2 > 0 & U2 > U1] <- 3

dat <- data.frame(Y = Y_choice, intercept = intercept, x1 = x1, eval = 1)

# Define model with FREE factor loadings
fm <- define_factor_model(
  n_factors = 1,
  n_types = 1,
  n_quad = 8,
  loading_normalization = NA_real_  # FREE loading
)

mc <- define_model_component(
  name = "Y",
  data = dat,
  outcome = "Y",
  factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "logit",
  evaluation_indicator = "eval",
  num_choices = 3
)

ms <- define_model_system(factor = fm, components = list(mc))
fm_cpp <- initialize_factor_model_cpp(ms, dat, n_quad = 8)

param_info <- get_parameter_info_cpp(fm_cpp)
cat("Number of parameters:", param_info$n_param_free, "\n\n")

# Test parameters
test_params <- c(1.0, 0.5, 0.3, 0.4, -0.3, 0.5, 0.6)
param_names <- c("sigma_f^2", "beta0_c1", "beta_x1_c1", "lambda_c1",
                 "beta0_c2", "beta_x1_c2", "lambda_c2")

cat("=== ANALYTICAL GRADIENT ===\n")
result <- evaluate_likelihood_cpp(fm_cpp, test_params,
                                 compute_gradient = TRUE,
                                 compute_hessian = FALSE)
for (i in seq_along(result$gradient)) {
  cat(sprintf("[%d] %-12s: %10.6f\n", i, param_names[i], result$gradient[i]))
}

cat("\n=== INDIVIDUAL FINITE DIFFERENCES ===\n")
delta <- 1e-7

for (i in seq_along(test_params)) {
  h <- delta * (abs(test_params[i]) + 1.0)

  params_plus <- test_params
  params_plus[i] <- test_params[i] + h
  f_plus <- evaluate_loglik_only_cpp(fm_cpp, params_plus)

  params_minus <- test_params
  params_minus[i] <- test_params[i] - h
  f_minus <- evaluate_loglik_only_cpp(fm_cpp, params_minus)

  grad_fd <- (f_plus - f_minus) / (2.0 * h)

  cat(sprintf("[%d] %-12s: FD=%10.6f, Analytical=%10.6f, Diff=%10.6f\n",
              i, param_names[i], grad_fd, result$gradient[i],
              abs(grad_fd - result$gradient[i])))
}

cat("\n=== CHECKING PARAMETER PERTURBATIONS ===\n")
# Check if perturbing param[5] (beta0_c2) affects likelihood correctly
cat("Perturbing beta0_c2 (position 5):\n")
params_base <- test_params
ll_base <- evaluate_loglik_only_cpp(fm_cpp, params_base)

params_pert <- test_params
params_pert[5] <- params_pert[5] + 0.1  # Perturb beta0_c2
ll_pert <- evaluate_loglik_only_cpp(fm_cpp, params_pert)

cat(sprintf("  Base LL: %.6f\n", ll_base))
cat(sprintf("  Perturbed LL: %.6f\n", ll_pert))
cat(sprintf("  Change: %.6f\n", ll_pert - ll_base))
cat(sprintf("  Expected gradient direction: %s\n",
            ifelse(result$gradient[5] > 0, "positive", "negative")))
cat(sprintf("  Observed change direction: %s\n",
            ifelse(ll_pert - ll_base > 0, "positive", "negative")))
