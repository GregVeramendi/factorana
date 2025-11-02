#!/usr/bin/env Rscript
# Compare Roy model estimation from two different starting points

library(factorana)

set.seed(456)
n <- 2000

cat("=============================================================\n")
cat("Roy Model - Compare Starting Points\n")
cat("=============================================================\n\n")

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

cat(sprintf("Sample size: %d\n", n))
cat(sprintf("Selection rate (D=1): %.1f%%\n\n", 100 * mean(D)))

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

# Parameter setup
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

param_fixed <- c(TRUE, rep(FALSE, 23))

# Default initialization (from data)
init_default <- c(
  1.0,  # factor_var (fixed)
  0.0, 0.5, 0.8,  # Selection
  mean(dat$Y1[dat$eval_Y1==1], na.rm=TRUE), 0.8, 1.5, sd(dat$Y1[dat$eval_Y1==1], na.rm=TRUE),  # Y1
  mean(dat$Y0[dat$eval_Y0==1], na.rm=TRUE), 0.8, 0.8, sd(dat$Y0[dat$eval_Y0==1], na.rm=TRUE),  # Y0
  mean(dat$T1), 0.8, 0.8, sd(dat$T1),  # T1
  mean(dat$T2), 0.8, 1.0, sd(dat$T2),  # T2
  mean(dat$T3), 0.8, 0.8, sd(dat$T3)   # T3
)

# Objective and gradient functions
obj_fn <- function(p) {
  full_p <- init_default
  full_p[!param_fixed] <- p
  -evaluate_loglik_only_cpp(fm_cpp, full_p)
}

grad_fn <- function(p) {
  full_p <- init_default
  full_p[!param_fixed] <- p
  grad_full <- evaluate_likelihood_cpp(fm_cpp, full_p, TRUE, FALSE)$gradient
  -grad_full[!param_fixed]
}

cat("=== ESTIMATION 1: Start from Default Initialization ===\n\n")

opt1 <- nlminb(
  init_default[!param_fixed],
  obj_fn,
  gradient = grad_fn,
  control = list(eval.max = 2000, iter.max = 2000)
)

params1 <- init_default
params1[!param_fixed] <- opt1$par
ll1 <- -opt1$objective

cat(sprintf("Converged in %d iterations\n", opt1$iterations))
cat(sprintf("Log-likelihood: %.6f\n", ll1))
cat(sprintf("Message: %s\n\n", opt1$message))

cat("=== ESTIMATION 2: Start from TRUE Parameters ===\n\n")

# Update init to use true parameters as starting point
init_from_truth <- true_params

obj_fn2 <- function(p) {
  full_p <- init_from_truth
  full_p[!param_fixed] <- p
  -evaluate_loglik_only_cpp(fm_cpp, full_p)
}

grad_fn2 <- function(p) {
  full_p <- init_from_truth
  full_p[!param_fixed] <- p
  grad_full <- evaluate_likelihood_cpp(fm_cpp, full_p, TRUE, FALSE)$gradient
  -grad_full[!param_fixed]
}

opt2 <- nlminb(
  init_from_truth[!param_fixed],
  obj_fn2,
  gradient = grad_fn2,
  control = list(eval.max = 2000, iter.max = 2000)
)

params2 <- init_from_truth
params2[!param_fixed] <- opt2$par
ll2 <- -opt2$objective

cat(sprintf("Converged in %d iterations\n", opt2$iterations))
cat(sprintf("Log-likelihood: %.6f\n", ll2))
cat(sprintf("Message: %s\n\n", opt2$message))

cat("=== COMPARISON ===\n\n")

cat(sprintf("LL difference (Start2 - Start1): %.6f\n", ll2 - ll1))
cat(sprintf("Absolute difference: %.6f\n\n", abs(ll2 - ll1)))

if (abs(ll2 - ll1) < 0.01) {
  cat("✓ Both starting points converged to SAME optimum (within 0.01)\n\n")
  cat("This suggests a single global optimum (good!)\n")
} else {
  cat("⚠ Different optima found!\n\n")
  cat("This suggests multiple local optima.\n")
}

# Compare parameter estimates
cat("\nParameter comparison:\n")
cat(sprintf("%-12s %10s %10s %10s %10s\n",
            "Parameter", "True", "Start1", "Start2", "Diff(1-2)"))
cat(strrep("-", 65), "\n")

max_diff <- 0
for (i in seq_along(param_names)) {
  if (param_fixed[i]) {
    cat(sprintf("%-12s %10.4f %10s %10s %10s (fixed)\n",
                param_names[i], true_params[i], "-", "-", "-"))
  } else {
    diff <- params1[i] - params2[i]
    max_diff <- max(max_diff, abs(diff))
    cat(sprintf("%-12s %10.4f %10.4f %10.4f %10.4f\n",
                param_names[i], true_params[i], params1[i], params2[i], diff))
  }
}

cat(strrep("-", 65), "\n")
cat(sprintf("\nMax parameter difference: %.6f\n", max_diff))

if (max_diff < 0.001) {
  cat("\n✓ Parameters match (< 0.001 difference)\n")
  cat("Both starting points found the same MLE.\n")
} else if (max_diff < 0.01) {
  cat("\n≈ Parameters very close (< 0.01 difference)\n")
  cat("Likely the same optimum, minor numerical differences.\n")
} else {
  cat("\n⚠ Parameters differ significantly\n")
  cat("Multiple local optima exist.\n")
}
