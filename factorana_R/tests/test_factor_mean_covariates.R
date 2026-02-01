# Test script for factor mean covariates feature
# Tests:
# 1. Basic functionality
# 2. Gradient finite difference check
# 3. Hessian finite difference check
# 4. Parameter recovery (simple simulation)
#
# Model structure:
#   Factor distribution: f_i ~ N(gamma * (age_i - mean(age)), sigma^2_f)
#   Measurement equations:
#     Y1 = 1.0 * f_i + beta_age * age_i + epsilon_1  (age in both factor AND outcome)
#     Y2 = lambda_2 * f_i + epsilon_2
#     Y3 = lambda_3 * f_i + epsilon_3
#
# This tests cross-derivatives between gamma (factor mean) and beta_age (model covariate)

library(factorana)
set.seed(12345)

cat("========================================\n")
cat("Test 1: Basic functionality\n")
cat("========================================\n")

# Create simple model with factor mean covariates
n_obs <- 500
n_factors <- 1

# Generate data: one factor measured by 3 continuous outcomes
# Factor mean depends on a covariate (age)
# Y1 also has a direct effect of age (same covariate in factor AND model)
true_loading <- c(1.0, 0.8, 0.6)  # First loading fixed to 1 for identification
true_sigma <- c(0.5, 0.5, 0.5)
true_factor_var <- 1.0
true_gamma <- 0.3      # Coefficient on factor mean covariate
true_beta_age <- 0.05  # Direct effect of age on Y1

# Generate covariate (will be demeaned internally for factor mean)
age <- rnorm(n_obs, mean = 40, sd = 10)

# Generate factor values with covariate-dependent mean
# f_i ~ N(gamma * (age_i - mean(age)), factor_var)
age_demeaned <- age - mean(age)
f <- rnorm(n_obs, mean = true_gamma * age_demeaned, sd = sqrt(true_factor_var))

# Generate outcomes
# Y1 has both factor loading AND direct age effect
Y1 <- true_loading[1] * f + true_beta_age * age + rnorm(n_obs, sd = sqrt(true_sigma[1]))
Y2 <- true_loading[2] * f + rnorm(n_obs, sd = sqrt(true_sigma[2]))
Y3 <- true_loading[3] * f + rnorm(n_obs, sd = sqrt(true_sigma[3]))

# Create data frame
dat <- data.frame(Y1 = Y1, Y2 = Y2, Y3 = Y3, age = age)

cat("Data generated:\n")
cat("  n_obs:", n_obs, "\n")
cat("  True gamma (factor mean coef):", true_gamma, "\n")
cat("  True beta_age (direct effect on Y1):", true_beta_age, "\n")
cat("  True factor variance:", true_factor_var, "\n")
cat("  Note: age appears in BOTH factor distribution AND Y1 measurement equation\n")
cat("\n")

# Define factor model with covariates
Factor <- define_factor_model(
  n_factors = 1,
  factor_covariates = c("age")
)

cat("Factor model created:\n")
print(Factor)
cat("\n")

# Define model components (using correct API)
# Y1 includes age as a covariate - same variable used in factor distribution
# This tests cross-derivatives between gamma (factor mean) and beta (model covariate)
comp1 <- define_model_component(
  name = "Y1",
  data = dat,
  outcome = "Y1",
  factor = Factor,
  covariates = c("age"),
  model_type = "linear",
  intercept = FALSE,
  loading_normalization = c(1)  # Fix loading to 1 for identification
)

comp2 <- define_model_component(
  name = "Y2",
  data = dat,
  outcome = "Y2",
  factor = Factor,
  covariates = NULL,
  model_type = "linear",
  intercept = FALSE
)

comp3 <- define_model_component(
  name = "Y3",
  data = dat,
  outcome = "Y3",
  factor = Factor,
  covariates = NULL,
  model_type = "linear",
  intercept = FALSE
)

# Define model system
model_system <- define_model_system(
  components = list(comp1, comp2, comp3),
  factor = Factor
)

cat("Model system created\n\n")

cat("========================================\n")
cat("Test 2: Gradient finite difference check\n")
cat("========================================\n")

# Initialize parameters
init_result <- initialize_parameters(model_system, dat)
init_params <- init_result$init_params
cat("Initial parameters:\n")
print(init_params)
cat("\n")

# Set up estimation control
control <- define_estimation_control(n_quad_points = 8)

# Create FactorModel object using internal function
fm_ptr <- factorana:::initialize_factor_model_cpp(
  model_system = model_system,
  data = dat,
  n_quad = 8
)

# Debug: check parameter info
param_info <- factorana:::get_parameter_info_cpp(fm_ptr)
cat("C++ parameter info:\n")
cat("  n_param:", param_info$n_param, "\n")
cat("  n_param_free:", param_info$n_param_free, "\n")
cat("  R init_params length:", length(init_params), "\n")

# Evaluate at initial params
result <- factorana:::evaluate_likelihood_cpp(
  fm_ptr = fm_ptr,
  params = init_params,
  compute_gradient = TRUE,
  compute_hessian = FALSE
)

cat("Debug: Log-likelihood =", result$logLikelihood, "\n")
cat("Debug: Gradient =", result$gradient, "\n")

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


cat("========================================\n")
cat("Test 3: Hessian finite difference check\n")
cat("========================================\n")

# Get analytic Hessian
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

# Relative difference (careful with near-zero entries)
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

# Show factor_mean row/column specifically
factor_mean_idx <- grep("factor_mean", names(init_params))
if (length(factor_mean_idx) > 0) {
  cat("Hessian row for factor_mean parameter:\n")
  cat(sprintf("  %-25s %12s %12s %12s\n", "Cross with", "Analytic", "Numerical", "Rel Diff"))
  for (j in 1:n_params) {
    cat(sprintf("  %-25s %12.4f %12.4f %12.6e\n",
                names(init_params)[j],
                analytic_hess[factor_mean_idx, j],
                numerical_hess[factor_mean_idx, j],
                rel_hess_diff[factor_mean_idx, j]))
  }
  cat("\n")
}

# Show cross-derivative between factor_mean_1_age and Y1_age specifically
# This is the key test: same covariate in both factor distribution AND model
Y1_age_idx <- grep("Y1_age", names(init_params))
if (length(factor_mean_idx) > 0 && length(Y1_age_idx) > 0) {
  cat("KEY TEST: Cross-derivative between factor_mean_1_age and Y1_age:\n")
  cat("  (Same covariate appears in both factor distribution and measurement equation)\n")
  cat(sprintf("  Analytic: %12.6f\n", analytic_hess[factor_mean_idx, Y1_age_idx]))
  cat(sprintf("  Numerical: %12.6f\n", numerical_hess[factor_mean_idx, Y1_age_idx]))
  cat(sprintf("  Rel Diff: %12.6e\n", rel_hess_diff[factor_mean_idx, Y1_age_idx]))
  cat("\n")
}

cat("Max absolute difference:", max_hess_diff, "\n")
cat("Max relative difference:", max_rel_hess_diff, "\n")

if (max_rel_hess_diff < 0.05) {
  cat("HESSIAN CHECK: PASSED\n\n")
} else {
  cat("HESSIAN CHECK: FAILED (may need looser tolerance for Hessian)\n\n")
}


cat("========================================\n")
cat("Test 4: Parameter recovery\n")
cat("========================================\n")

# Run estimation
result_est <- estimate_model_rcpp(
  model_system = model_system,
  data = dat,
  control = define_estimation_control(n_quad_points = 12),
  init_params = init_params
)

cat("\nEstimation results:\n")
cat("Convergence:", result_est$convergence, "\n")
cat("Log-likelihood:", result_est$logLikelihood, "\n\n")

# Compare estimates to true values
cat("Parameter recovery:\n")
cat(sprintf("  %-25s %12s %12s %12s\n", "Parameter", "True", "Estimated", "Std Error"))

# Factor variance
est_factor_var <- result_est$estimates["factor_var_1"]
se_factor_var <- result_est$std_errors["factor_var_1"]
cat(sprintf("  %-25s %12.4f %12.4f %12.4f\n", "factor_var_1", true_factor_var, est_factor_var, se_factor_var))

# Factor mean coefficient (gamma)
est_gamma <- result_est$estimates["factor_mean_1_age"]
se_gamma <- result_est$std_errors["factor_mean_1_age"]
cat(sprintf("  %-25s %12.4f %12.4f %12.4f\n", "factor_mean_1_age", true_gamma, est_gamma, se_gamma))

# Direct effect of age on Y1 (beta_age) - same covariate as in factor distribution
est_beta_age <- result_est$estimates["Y1_age"]
se_beta_age <- result_est$std_errors["Y1_age"]
cat(sprintf("  %-25s %12.4f %12.4f %12.4f\n", "Y1_age (direct effect)", true_beta_age, est_beta_age, se_beta_age))

# Loadings
est_loading2 <- result_est$estimates["Y2_loading_1"]
se_loading2 <- result_est$std_errors["Y2_loading_1"]
cat(sprintf("  %-25s %12.4f %12.4f %12.4f\n", "Y2_loading_1", true_loading[2], est_loading2, se_loading2))

est_loading3 <- result_est$estimates["Y3_loading_1"]
se_loading3 <- result_est$std_errors["Y3_loading_1"]
cat(sprintf("  %-25s %12.4f %12.4f %12.4f\n", "Y3_loading_1", true_loading[3], est_loading3, se_loading3))

# Sigmas
est_sigma1 <- result_est$estimates["Y1_sigma"]
se_sigma1 <- result_est$std_errors["Y1_sigma"]
cat(sprintf("  %-25s %12.4f %12.4f %12.4f\n", "Y1_sigma", sqrt(true_sigma[1]), est_sigma1, se_sigma1))

cat("\n")

# Check if gamma and beta_age are within 2 SEs of true values
gamma_in_ci <- abs(est_gamma - true_gamma) < 2 * se_gamma
beta_in_ci <- abs(est_beta_age - true_beta_age) < 2 * se_beta_age
cat("Factor mean coefficient (gamma) within 2 SEs of true:", gamma_in_ci, "\n")
cat("Direct effect (beta_age) within 2 SEs of true:", beta_in_ci, "\n")

if (gamma_in_ci && beta_in_ci) {
  cat("PARAMETER RECOVERY: PASSED\n")
} else {
  cat("PARAMETER RECOVERY: CHECK MANUALLY\n")
}

cat("\n========================================\n")
cat("All tests completed\n")
cat("========================================\n")
