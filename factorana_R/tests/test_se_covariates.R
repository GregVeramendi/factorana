# Test script for SE covariates feature
# Tests SE_linear with covariates: f_2 = intercept + alpha*f_1 + beta*X + epsilon
#
# Model structure:
#   Input factor: f_1 ~ N(0, sigma^2_1)
#   Outcome factor: f_2 = alpha_0 + alpha_1 * f_1 + beta * X + epsilon
#   Measurement equations:
#     Y1 = 1.0 * f_1 + e1   (loading fixed for identification)
#     Y2 = lambda_2 * f_1 + e2
#     Y3 = 1.0 * f_2 + e3   (loading fixed for identification)
#     Y4 = lambda_4 * f_2 + e4
#
# Tests:
# 1. Gradient finite difference check
# 2. Hessian finite difference check
# 3. Parameter recovery

library(factorana)
set.seed(54321)

cat("========================================\n")
cat("Test: SE covariates\n")
cat("========================================\n")

# Generate data
n_obs <- 500

# True parameters
true_factor_var_1 <- 1.0      # Input factor variance
true_se_intercept <- 0.5      # SE intercept
true_se_linear <- 0.6         # SE linear coefficient on f_1
true_se_beta <- 0.4           # SE covariate coefficient
true_se_residual_var <- 0.8   # SE residual variance

true_loading_Y2 <- 0.8        # Loading of Y2 on f_1
true_loading_Y4 <- 0.7        # Loading of Y4 on f_2

true_sigma <- c(0.5, 0.5, 0.5, 0.5)  # Measurement error SDs

# Generate covariate
X <- rnorm(n_obs, mean = 2, sd = 1)
X_demeaned <- X - mean(X)

# Generate latent factors
f1 <- rnorm(n_obs, mean = 0, sd = sqrt(true_factor_var_1))
eps <- rnorm(n_obs, mean = 0, sd = sqrt(true_se_residual_var))
f2 <- true_se_intercept + true_se_linear * f1 + true_se_beta * X_demeaned + eps

# Generate outcomes
Y1 <- 1.0 * f1 + rnorm(n_obs, sd = sqrt(true_sigma[1]))
Y2 <- true_loading_Y2 * f1 + rnorm(n_obs, sd = sqrt(true_sigma[2]))
Y3 <- 1.0 * f2 + rnorm(n_obs, sd = sqrt(true_sigma[3]))
Y4 <- true_loading_Y4 * f2 + rnorm(n_obs, sd = sqrt(true_sigma[4]))

dat <- data.frame(Y1 = Y1, Y2 = Y2, Y3 = Y3, Y4 = Y4, X = X)

cat("Data generated:\n")
cat("  n_obs:", n_obs, "\n")
cat("  True SE intercept:", true_se_intercept, "\n")
cat("  True SE linear (alpha_1):", true_se_linear, "\n")
cat("  True SE covariate (beta):", true_se_beta, "\n")
cat("  True SE residual var:", true_se_residual_var, "\n")
cat("\n")

# Define factor model with SE structure and covariates
Factor <- define_factor_model(
  n_factors = 2,
  factor_structure = "SE_linear",
  se_covariates = c("X")
)

cat("Factor model:\n")
print(Factor)
cat("\n")

# Define model components
# Y1 and Y2 load on f_1 (input factor)
# Y3 and Y4 load on f_2 (outcome factor)
comp1 <- define_model_component(
  name = "Y1",
  data = dat,
  outcome = "Y1",
  factor = Factor,
  covariates = NULL,
  model_type = "linear",
  intercept = FALSE,
  loading_normalization = c(1, NA_real_)  # Fix loading on f_1 to 1, NA for f_2
)

comp2 <- define_model_component(
  name = "Y2",
  data = dat,
  outcome = "Y2",
  factor = Factor,
  covariates = NULL,
  model_type = "linear",
  intercept = FALSE,
  loading_normalization = c(NA_real_, NA_real_)  # Free loading on f_1
)

comp3 <- define_model_component(
  name = "Y3",
  data = dat,
  outcome = "Y3",
  factor = Factor,
  covariates = NULL,
  model_type = "linear",
  intercept = FALSE,
  loading_normalization = c(NA_real_, 1)  # Fix loading on f_2 to 1
)

comp4 <- define_model_component(
  name = "Y4",
  data = dat,
  outcome = "Y4",
  factor = Factor,
  covariates = NULL,
  model_type = "linear",
  intercept = FALSE,
  loading_normalization = NULL  # All free
)

model_system <- define_model_system(
  components = list(comp1, comp2, comp3, comp4),
  factor = Factor
)

cat("Model system created\n\n")

# Initialize parameters
init_result <- initialize_parameters(model_system, dat)
init_params <- init_result$init_params

cat("Initial parameters:\n")
print(init_params)
cat("\n")

# Create FactorModel object
fm_ptr <- factorana:::initialize_factor_model_cpp(
  model_system = model_system,
  data = dat,
  n_quad = 8
)

param_info <- factorana:::get_parameter_info_cpp(fm_ptr)
cat("C++ parameter info:\n")
cat("  n_param:", param_info$n_param, "\n")
cat("  n_param_free:", param_info$n_param_free, "\n")
cat("  R init_params length:", length(init_params), "\n\n")

# ========================================
# Test 1: Gradient finite difference check
# ========================================
cat("========================================\n")
cat("Test 1: Gradient finite difference check\n")
cat("========================================\n")

result <- factorana:::evaluate_likelihood_cpp(
  fm_ptr = fm_ptr,
  params = init_params,
  compute_gradient = TRUE,
  compute_hessian = FALSE
)

cat("Log-likelihood:", result$logLikelihood, "\n")
analytic_grad <- result$gradient
n_params <- length(init_params)

# Compute numerical gradient
eps <- 1e-5
numerical_grad <- rep(0, n_params)

for (i in 1:n_params) {
  params_plus <- init_params
  params_minus <- init_params
  params_plus[i] <- params_plus[i] + eps
  params_minus[i] <- params_minus[i] - eps

  result_plus <- factorana:::evaluate_loglik_only_cpp(fm_ptr = fm_ptr, params = params_plus)
  result_minus <- factorana:::evaluate_loglik_only_cpp(fm_ptr = fm_ptr, params = params_minus)

  numerical_grad[i] <- (result_plus - result_minus) / (2 * eps)
}

# Compare
grad_diff <- abs(analytic_grad - numerical_grad)
max_grad_diff <- max(grad_diff)
rel_grad_diff <- grad_diff / (abs(numerical_grad) + 1e-8)
max_rel_grad_diff <- max(rel_grad_diff)

cat("Gradient comparison:\n")
cat(sprintf("  %-25s %12s %12s %12s\n", "Parameter", "Analytic", "Numerical", "Rel Diff"))
for (i in 1:n_params) {
  cat(sprintf("  %-25s %12.6f %12.6f %12.6e\n",
              names(init_params)[i], analytic_grad[i], numerical_grad[i], rel_grad_diff[i]))
}
cat("\n")
cat("Max absolute difference:", max_grad_diff, "\n")
cat("Max relative difference:", max_rel_grad_diff, "\n")

if (max_rel_grad_diff < 1e-4) {
  cat("GRADIENT CHECK: PASSED\n\n")
} else {
  cat("GRADIENT CHECK: FAILED\n\n")
}

# ========================================
# Test 2: Hessian finite difference check
# ========================================
cat("========================================\n")
cat("Test 2: Hessian finite difference check\n")
cat("========================================\n")

result_hess <- factorana:::evaluate_likelihood_cpp(
  fm_ptr = fm_ptr,
  params = init_params,
  compute_gradient = TRUE,
  compute_hessian = TRUE
)

# Expand upper triangular Hessian to full matrix
hess_upper <- result_hess$hessian
analytic_hess <- matrix(0, n_params, n_params)
idx <- 1
for (i in 1:n_params) {
  for (j in i:n_params) {
    analytic_hess[i, j] <- hess_upper[idx]
    analytic_hess[j, i] <- hess_upper[idx]
    idx <- idx + 1
  }
}

# Compute numerical Hessian (from gradient)
numerical_hess <- matrix(0, n_params, n_params)
eps_hess <- 1e-4

for (i in 1:n_params) {
  params_plus <- init_params
  params_minus <- init_params
  params_plus[i] <- params_plus[i] + eps_hess
  params_minus[i] <- params_minus[i] - eps_hess

  result_plus <- factorana:::evaluate_likelihood_cpp(
    fm_ptr = fm_ptr,
    params = params_plus,
    compute_gradient = TRUE,
    compute_hessian = FALSE
  )

  result_minus <- factorana:::evaluate_likelihood_cpp(
    fm_ptr = fm_ptr,
    params = params_minus,
    compute_gradient = TRUE,
    compute_hessian = FALSE
  )

  numerical_hess[, i] <- (result_plus$gradient - result_minus$gradient) / (2 * eps_hess)
}

# Symmetrize
numerical_hess <- (numerical_hess + t(numerical_hess)) / 2

# Compare
hess_diff <- abs(analytic_hess - numerical_hess)
max_hess_diff <- max(hess_diff)

hess_scale <- pmax(abs(numerical_hess), 1e-6)
rel_hess_diff <- hess_diff / hess_scale
max_rel_hess_diff <- max(rel_hess_diff)

cat("Hessian comparison (diagonal elements):\n")
cat(sprintf("  %-25s %12s %12s %12s\n", "Parameter", "Analytic", "Numerical", "Rel Diff"))
for (i in 1:n_params) {
  cat(sprintf("  %-25s %12.4f %12.4f %12.6e\n",
              names(init_params)[i], analytic_hess[i, i], numerical_hess[i, i],
              rel_hess_diff[i, i]))
}
cat("\n")

# Show SE covariate row specifically
se_cov_idx <- grep("se_cov", names(init_params))
if (length(se_cov_idx) > 0) {
  cat("Hessian row for SE covariate parameter:\n")
  cat(sprintf("  %-25s %12s %12s %12s\n", "Cross with", "Analytic", "Numerical", "Rel Diff"))
  for (j in 1:n_params) {
    cat(sprintf("  %-25s %12.4f %12.4f %12.6e\n",
                names(init_params)[j],
                analytic_hess[se_cov_idx, j],
                numerical_hess[se_cov_idx, j],
                rel_hess_diff[se_cov_idx, j]))
  }
  cat("\n")
}

cat("Max absolute difference:", max_hess_diff, "\n")
cat("Max relative difference:", max_rel_hess_diff, "\n")

if (max_rel_hess_diff < 0.05) {
  cat("HESSIAN CHECK: PASSED\n\n")
} else {
  cat("HESSIAN CHECK: FAILED\n\n")
}

# ========================================
# Test 3: Parameter recovery
# ========================================
cat("========================================\n")
cat("Test 3: Parameter recovery\n")
cat("========================================\n")

result_est <- estimate_model_rcpp(
  model_system = model_system,
  data = dat,
  control = define_estimation_control(n_quad_points = 12),
  init_params = init_params
)

cat("\nEstimation results:\n")
cat("Convergence:", result_est$convergence, "\n")
cat("Log-likelihood:", result_est$logLikelihood, "\n\n")

cat("Parameter recovery:\n")
cat(sprintf("  %-25s %12s %12s %12s\n", "Parameter", "True", "Estimated", "Std Error"))

# Input factor variance
est_var <- result_est$estimates["factor_var_1"]
se_var <- result_est$std_errors["factor_var_1"]
cat(sprintf("  %-25s %12.4f %12.4f %12.4f\n", "factor_var_1", true_factor_var_1, est_var, se_var))

# SE intercept
est_int <- result_est$estimates["se_intercept"]
se_int <- result_est$std_errors["se_intercept"]
cat(sprintf("  %-25s %12.4f %12.4f %12.4f\n", "se_intercept", true_se_intercept, est_int, se_int))

# SE linear
est_lin <- result_est$estimates["se_linear_1"]
se_lin <- result_est$std_errors["se_linear_1"]
cat(sprintf("  %-25s %12.4f %12.4f %12.4f\n", "se_linear_1", true_se_linear, est_lin, se_lin))

# SE residual var
est_res <- result_est$estimates["se_residual_var"]
se_res <- result_est$std_errors["se_residual_var"]
cat(sprintf("  %-25s %12.4f %12.4f %12.4f\n", "se_residual_var", true_se_residual_var, est_res, se_res))

# SE covariate (the key test!)
est_beta <- result_est$estimates["se_cov_X"]
se_beta <- result_est$std_errors["se_cov_X"]
cat(sprintf("  %-25s %12.4f %12.4f %12.4f\n", "se_cov_X (beta)", true_se_beta, est_beta, se_beta))

# Loadings
est_l2 <- result_est$estimates["Y2_loading_1"]
se_l2 <- result_est$std_errors["Y2_loading_1"]
cat(sprintf("  %-25s %12.4f %12.4f %12.4f\n", "Y2_loading_1", true_loading_Y2, est_l2, se_l2))

est_l4 <- result_est$estimates["Y4_loading_2"]
se_l4 <- result_est$std_errors["Y4_loading_2"]
cat(sprintf("  %-25s %12.4f %12.4f %12.4f\n", "Y4_loading_2", true_loading_Y4, est_l4, se_l4))

cat("\n")

# Check if key parameters are within 2 SEs of true values
beta_in_ci <- abs(est_beta - true_se_beta) < 2 * se_beta
linear_in_ci <- abs(est_lin - true_se_linear) < 2 * se_lin
cat("SE covariate (beta) within 2 SEs of true:", beta_in_ci, "\n")
cat("SE linear (alpha_1) within 2 SEs of true:", linear_in_ci, "\n")

if (beta_in_ci && linear_in_ci) {
  cat("PARAMETER RECOVERY: PASSED\n")
} else {
  cat("PARAMETER RECOVERY: CHECK MANUALLY\n")
}

cat("\n========================================\n")
cat("All tests completed\n")
cat("========================================\n")
