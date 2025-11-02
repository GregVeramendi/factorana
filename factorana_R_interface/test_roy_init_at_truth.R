#!/usr/bin/env Rscript
# Test: Roy model starting optimization at TRUE parameters

library(factorana)

cat("=============================================================\n")
cat("Roy Model - Start Optimization at True Parameters\n")
cat("=============================================================\n\n")

set.seed(456)
n <- 2000

# TRUE PARAMETERS
true_sel_intercept <- 0.5
true_sel_beta_Z1 <- 1.0
true_sel_lambda_f <- 1.0

true_Y1_intercept <- 2.0
true_Y1_beta_X1 <- 1.0
true_Y1_lambda_f <- 2.0
true_Y1_sigma <- 0.5

true_Y0_intercept <- 1.0
true_Y0_beta_X1 <- 1.0
true_Y0_lambda_f <- 1.0
true_Y0_sigma <- 0.5

true_T1_intercept <- 1.0
true_T1_beta_Q1 <- 1.0
true_T1_lambda_f <- 1.0
true_T1_sigma <- 0.5

true_T2_intercept <- 0.5
true_T2_beta_Q2 <- 0.8
true_T2_lambda_f <- 1.2
true_T2_sigma <- 0.6

true_T3_intercept <- -0.5
true_T3_beta_Q3 <- 1.1
true_T3_lambda_f <- 0.9
true_T3_sigma <- 0.4

true_factor_var <- 1.0

# Generate data
f <- rnorm(n, 0, 1)
U1 <- rnorm(n, 0, 1) * true_Y1_sigma
U0 <- rnorm(n, 0, 1) * true_Y0_sigma
UV <- rnorm(n, 0, 1)
UT1 <- rnorm(n, 0, 1) * true_T1_sigma
UT2 <- rnorm(n, 0, 1) * true_T2_sigma
UT3 <- rnorm(n, 0, 1) * true_T3_sigma

intercept <- rep(1, n)
Z1 <- rnorm(n, 0, 1)
X1 <- rnorm(n, 0, 1)
Q1 <- rnorm(n, 0, 1)
Q2 <- rnorm(n, 0, 1)
Q3 <- rnorm(n, 0, 1)

Y1 <- true_Y1_intercept * intercept + true_Y1_beta_X1 * X1 + true_Y1_lambda_f * f + U1
Y0 <- true_Y0_intercept * intercept + true_Y0_beta_X1 * X1 + true_Y0_lambda_f * f + U0

I_latent <- true_sel_intercept * intercept + true_sel_beta_Z1 * Z1 + true_sel_lambda_f * f + UV
D <- as.integer(I_latent > 0)

T1_score <- true_T1_intercept * intercept + true_T1_beta_Q1 * Q1 + true_T1_lambda_f * f + UT1
T2_score <- true_T2_intercept * intercept + true_T2_beta_Q2 * Q2 + true_T2_lambda_f * f + UT2
T3_score <- true_T3_intercept * intercept + true_T3_beta_Q3 * Q3 + true_T3_lambda_f * f + UT3

dat <- data.frame(
  intercept = intercept, D = D, Z1 = Z1, X1 = X1,
  Y1 = Y1, Y0 = Y0,
  Q1 = Q1, Q2 = Q2, Q3 = Q3,
  T1 = T1_score, T2 = T2_score, T3 = T3_score,
  eval_Y1 = D, eval_Y0 = 1L - D, eval_all = 1L
)

cat(sprintf("  Sample size: %d\n", n))
cat(sprintf("  Selection rate (D=1): %.1f%%\n", 100 * mean(D)))
cat(sprintf("  Y1 observations: %d\n", sum(dat$eval_Y1)))
cat(sprintf("  Y0 observations: %d\n\n", sum(dat$eval_Y0)))

# Define model
fm <- define_factor_model(n_factors = 1, n_types = 1, n_mixtures = 1, n_quad = 16, correlation = FALSE)

mc_sel <- define_model_component(
  name = "Selection", data = dat, outcome = "D", factor = fm,
  evaluation_indicator = "eval_all", covariates = c("intercept", "Z1"),
  model_type = "probit", loading_normalization = NA_real_
)

mc_Y1 <- define_model_component(
  name = "Y1", data = dat, outcome = "Y1", factor = fm,
  evaluation_indicator = "eval_Y1", covariates = c("intercept", "X1"),
  model_type = "linear", loading_normalization = NA_real_
)

mc_Y0 <- define_model_component(
  name = "Y0", data = dat, outcome = "Y0", factor = fm,
  evaluation_indicator = "eval_Y0", covariates = c("intercept", "X1"),
  model_type = "linear", loading_normalization = NA_real_
)

mc_T1 <- define_model_component(
  name = "TestScore1", data = dat, outcome = "T1", factor = fm,
  evaluation_indicator = "eval_all", covariates = c("intercept", "Q1"),
  model_type = "linear", loading_normalization = NA_real_
)

mc_T2 <- define_model_component(
  name = "TestScore2", data = dat, outcome = "T2", factor = fm,
  evaluation_indicator = "eval_all", covariates = c("intercept", "Q2"),
  model_type = "linear", loading_normalization = NA_real_
)

mc_T3 <- define_model_component(
  name = "TestScore3", data = dat, outcome = "T3", factor = fm,
  evaluation_indicator = "eval_all", covariates = c("intercept", "Q3"),
  model_type = "linear", loading_normalization = NA_real_
)

ms <- define_model_system(
  components = list(mc_sel, mc_Y1, mc_Y0, mc_T1, mc_T2, mc_T3),
  factor = fm
)

fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 16)

# TRUE parameter vector (all free)
true_params <- c(
  true_factor_var,
  true_sel_intercept, true_sel_beta_Z1, true_sel_lambda_f,
  true_Y1_intercept, true_Y1_beta_X1, true_Y1_lambda_f, true_Y1_sigma,
  true_Y0_intercept, true_Y0_beta_X1, true_Y0_lambda_f, true_Y0_sigma,
  true_T1_intercept, true_T1_beta_Q1, true_T1_lambda_f, true_T1_sigma,
  true_T2_intercept, true_T2_beta_Q2, true_T2_lambda_f, true_T2_sigma,
  true_T3_intercept, true_T3_beta_Q3, true_T3_lambda_f, true_T3_sigma
)

param_names <- c(
  "factor_var",
  "Sel_int", "Sel_beta_Z1", "Sel_lambda",
  "Y1_int", "Y1_beta_X1", "Y1_lambda", "Y1_sigma",
  "Y0_int", "Y0_beta_X1", "Y0_lambda", "Y0_sigma",
  "T1_int", "T1_beta_Q1", "T1_lambda", "T1_sigma",
  "T2_int", "T2_beta_Q2", "T2_lambda", "T2_sigma",
  "T3_int", "T3_beta_Q3", "T3_lambda", "T3_sigma"
)

# Fix factor variance
param_fixed <- c(TRUE, rep(FALSE, 23))

# Evaluate at true parameters
ll_true <- evaluate_loglik_only_cpp(fm_cpp, true_params)
cat(sprintf("Log-likelihood at true parameters: %.6f\n\n", ll_true))

# Check gradient at true parameters
result_true <- evaluate_likelihood_cpp(fm_cpp, true_params, TRUE, FALSE)
grad_norm <- sqrt(sum(result_true$gradient[!param_fixed]^2))
cat(sprintf("Gradient L2 norm at true parameters: %.6e\n", grad_norm))
cat(sprintf("Max gradient component: %.6f\n\n", max(abs(result_true$gradient[!param_fixed]))))

cat("=== TEST 1: Start at True Parameters ===\n\n")

opt_from_truth <- nlminb(
  true_params[!param_fixed],
  function(p) {
    full_p <- true_params
    full_p[!param_fixed] <- p
    -evaluate_loglik_only_cpp(fm_cpp, full_p)
  },
  gradient = function(p) {
    full_p <- true_params
    full_p[!param_fixed] <- p
    grad_full <- evaluate_likelihood_cpp(fm_cpp, full_p, TRUE, FALSE)$gradient
    -grad_full[!param_fixed]
  },
  control = list(eval.max = 2000, iter.max = 2000, trace = 1)
)

final_from_truth <- true_params
final_from_truth[!param_fixed] <- opt_from_truth$par

cat(sprintf("\nConverged in %d iterations\n", opt_from_truth$iterations))
cat(sprintf("Final log-likelihood: %.6f\n", -opt_from_truth$objective))
cat(sprintf("Change from start: %.6f\n\n", -opt_from_truth$objective - ll_true))

cat("Parameter changes from truth:\n")
cat(sprintf("%-12s %10s %10s %10s\n", "Parameter", "True", "Final", "Change"))
cat(strrep("-", 50), "\n")

max_change <- 0
for (i in which(!param_fixed)) {
  change <- final_from_truth[i] - true_params[i]
  pct_change <- 100 * change / true_params[i]
  max_change <- max(max_change, abs(pct_change))
  if (abs(pct_change) > 1) {
    cat(sprintf("%-12s %10.4f %10.4f %9.1f%%\n",
                param_names[i], true_params[i], final_from_truth[i], pct_change))
  }
}

cat(sprintf("\nMax parameter change: %.1f%%\n", max_change))
cat(sprintf("Message: %s\n", opt_from_truth$message))
