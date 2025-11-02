#!/usr/bin/env Rscript
# Check if the "better" optimum found by optimizer has zero gradient

library(factorana)

set.seed(456)
n <- 2000

# Generate same data
f <- rnorm(n, 0, 1)
U1 <- rnorm(n, 0, 1) * 0.5
U0 <- rnorm(n, 0, 1) * 0.5
UV <- rnorm(n, 0, 1)
UT1 <- rnorm(n, 0, 1) * 0.5
UT2 <- rnorm(n, 0, 1) * 0.6
UT3 <- rnorm(n, 0, 1) * 0.4

intercept <- rep(1, n)
Z1 <- rnorm(n, 0, 1)
X1 <- rnorm(n, 0, 1)
Q1 <- rnorm(n, 0, 1)
Q2 <- rnorm(n, 0, 1)
Q3 <- rnorm(n, 0, 1)

Y1 <- 2.0 * intercept + 1.0 * X1 + 2.0 * f + U1
Y0 <- 1.0 * intercept + 1.0 * X1 + 1.0 * f + U0

I_latent <- 0.5 * intercept + 1.0 * Z1 + 1.0 * f + UV
D <- as.integer(I_latent > 0)

T1 <- 1.0 * intercept + 1.0 * Q1 + 1.0 * f + UT1
T2 <- 0.5 * intercept + 0.8 * Q2 + 1.2 * f + UT2
T3 <- -0.5 * intercept + 1.1 * Q3 + 0.9 * f + UT3

dat <- data.frame(
  intercept = intercept, D = D, Z1 = Z1, X1 = X1,
  Y1 = Y1, Y0 = Y0, Q1 = Q1, Q2 = Q2, Q3 = Q3,
  T1 = T1, T2 = T2, T3 = T3,
  eval_Y1 = D, eval_Y0 = 1L - D, eval_all = 1L
)

# Setup model
fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 16)

mc_sel <- define_model_component(name = "Selection", data = dat, outcome = "D", factor = fm,
  evaluation_indicator = "eval_all", covariates = c("intercept", "Z1"),
  model_type = "probit", loading_normalization = NA_real_)

mc_Y1 <- define_model_component(name = "Y1", data = dat, outcome = "Y1", factor = fm,
  evaluation_indicator = "eval_Y1", covariates = c("intercept", "X1"),
  model_type = "linear", loading_normalization = NA_real_)

mc_Y0 <- define_model_component(name = "Y0", data = dat, outcome = "Y0", factor = fm,
  evaluation_indicator = "eval_Y0", covariates = c("intercept", "X1"),
  model_type = "linear", loading_normalization = NA_real_)

mc_T1 <- define_model_component(name = "T1", data = dat, outcome = "T1", factor = fm,
  evaluation_indicator = "eval_all", covariates = c("intercept", "Q1"),
  model_type = "linear", loading_normalization = NA_real_)

mc_T2 <- define_model_component(name = "T2", data = dat, outcome = "T2", factor = fm,
  evaluation_indicator = "eval_all", covariates = c("intercept", "Q2"),
  model_type = "linear", loading_normalization = NA_real_)

mc_T3 <- define_model_component(name = "T3", data = dat, outcome = "T3", factor = fm,
  evaluation_indicator = "eval_all", covariates = c("intercept", "Q3"),
  model_type = "linear", loading_normalization = NA_real_)

ms <- define_model_system(components = list(mc_sel, mc_Y1, mc_Y0, mc_T1, mc_T2, mc_T3), factor = fm)
fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 16)

# True parameters
true_params <- c(1.0, 0.5, 1.0, 1.0, 2.0, 1.0, 2.0, 0.5, 1.0, 1.0, 1.0, 0.5,
                 1.0, 1.0, 1.0, 0.5, 0.5, 0.8, 1.2, 0.6, -0.5, 1.1, 0.9, 0.4)

# "Better" optimum found when starting from truth
better_params <- c(1.0, 0.543646, 0.954845, 0.869052, 2.07391, 0.987290, 1.77241, 0.604773,
                   1.05008, 0.986139, 0.917588, 0.512469, 1.02943, 1.00890, 0.895049, 0.512399,
                   0.541506, 0.786979, 1.06937, 0.614604, -0.477337, 1.09220, 0.805507, 0.414306)

param_names <- c("factor_var", "Sel_int", "Sel_beta_Z1", "Sel_lambda",
                 "Y1_int", "Y1_beta_X1", "Y1_lambda", "Y1_sigma",
                 "Y0_int", "Y0_beta_X1", "Y0_lambda", "Y0_sigma",
                 "T1_int", "T1_beta_Q1", "T1_lambda", "T1_sigma",
                 "T2_int", "T2_beta_Q2", "T2_lambda", "T2_sigma",
                 "T3_int", "T3_beta_Q3", "T3_lambda", "T3_sigma")

param_fixed <- c(TRUE, rep(FALSE, 23))

cat("=== Checking Gradient at Two Points ===\n\n")

# At truth
result_truth <- evaluate_likelihood_cpp(fm_cpp, true_params, TRUE, FALSE)
ll_truth <- result_truth$loglik
grad_truth <- result_truth$gradient[!param_fixed]
grad_norm_truth <- sqrt(sum(grad_truth^2))

cat(sprintf("TRUE PARAMETERS:\n"))
cat(sprintf("  Log-likelihood: %.6f\n", ll_truth))
cat(sprintf("  Gradient L2 norm: %.6e\n", grad_norm_truth))
cat(sprintf("  Max |gradient|: %.6f\n\n", max(abs(grad_truth))))

# At "better" optimum
result_better <- evaluate_likelihood_cpp(fm_cpp, better_params, TRUE, FALSE)
ll_better <- result_better$loglik
grad_better <- result_better$gradient[!param_fixed]
grad_norm_better <- sqrt(sum(grad_better^2))

cat(sprintf('"BETTER" OPTIMUM (found from truth):\n'))
cat(sprintf("  Log-likelihood: %.6f\n", ll_better))
cat(sprintf("  Gradient L2 norm: %.6e\n", grad_norm_better))
cat(sprintf("  Max |gradient|: %.6f\n\n", max(abs(grad_better))))

cat(sprintf("Likelihood improvement: %.2f\n", ll_better - ll_truth))
cat(sprintf("Gradient is %s at 'better' point\n",
    ifelse(grad_norm_better < 1e-3, "NEAR ZERO (true optimum)", "NONZERO (not optimum)")))

# Check largest gradient components at both points
cat("\n=== Largest Gradient Components ===\n\n")
cat("At TRUE parameters:\n")
idx_truth <- order(abs(grad_truth), decreasing = TRUE)[1:5]
for (i in idx_truth) {
  cat(sprintf("  %-12s: %10.4f\n", param_names[!param_fixed][i], grad_truth[i]))
}

cat("\nAt 'BETTER' optimum:\n")
idx_better <- order(abs(grad_better), decreasing = TRUE)[1:5]
for (i in idx_better) {
  cat(sprintf("  %-12s: %10.4f\n", param_names[!param_fixed][i], grad_better[i]))
}
