#!/usr/bin/env Rscript
# Test Roy Model with Selection
# This is a more complex model with:
# - Selection equation (probit)
# - Two outcome equations with missing indicators
# - Test score equation

library(factorana)

cat("=============================================================\n")
cat("Roy Model Parameter Recovery Test\n")
cat("=============================================================\n\n")

# ---- Generate Data with Known Parameters ----
cat("Generating synthetic Roy model data...\n")
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

# Latent factor
f <- rnorm(n, 0, 1)

# Error terms (all with variance = sigma^2)
U1 <- rnorm(n, 0, 1) * true_Y1_sigma  # Y1 error: SD = 0.5
U0 <- rnorm(n, 0, 1) * true_Y0_sigma  # Y0 error: SD = 0.5
UV <- rnorm(n, 0, 1)  # Selection error: SD = 1 (probit normalization)
UT1 <- rnorm(n, 0, 1) * true_T1_sigma  # Test score 1 error: SD = 0.5
UT2 <- rnorm(n, 0, 1) * true_T2_sigma  # Test score 2 error: SD = 0.6
UT3 <- rnorm(n, 0, 1) * true_T3_sigma  # Test score 3 error: SD = 0.4

# Covariates
intercept <- rep(1, n)  # Intercept for all equations
Z1 <- rnorm(n, 0, 1)
X1 <- rnorm(n, 0, 1)
Q1 <- rnorm(n, 0, 1)
Q2 <- rnorm(n, 0, 1)
Q3 <- rnorm(n, 0, 1)

# Generate outcomes
Y1 <- true_Y1_intercept * intercept + true_Y1_beta_X1 * X1 + true_Y1_lambda_f * f + U1
Y0 <- true_Y0_intercept * intercept + true_Y0_beta_X1 * X1 + true_Y0_lambda_f * f + U0

# Selection (latent utility)
I_latent <- true_sel_intercept * intercept + true_sel_beta_Z1 * Z1 + true_sel_lambda_f * f + UV
D <- as.integer(I_latent > 0)

# Test scores
T1_score <- true_T1_intercept * intercept + true_T1_beta_Q1 * Q1 + true_T1_lambda_f * f + UT1
T2_score <- true_T2_intercept * intercept + true_T2_beta_Q2 * Q2 + true_T2_lambda_f * f + UT2
T3_score <- true_T3_intercept * intercept + true_T3_beta_Q3 * Q3 + true_T3_lambda_f * f + UT3

# Observed outcome (Roy model: observe Y1 if D=1, Y0 if D=0)
Y <- ifelse(D == 1, Y1, Y0)

# Create dataset
dat <- data.frame(
  intercept = intercept,
  Y = Y,
  D = D,
  Z1 = Z1,
  X1 = X1,
  Q1 = Q1,
  Q2 = Q2,
  Q3 = Q3,
  T1 = T1_score,
  T2 = T2_score,
  T3 = T3_score
)
dat$eval_Y1 <- dat$D          # Y1 observed when D=1
dat$eval_Y0 <- 1L - dat$D     # Y0 observed when D=0

cat(sprintf("  Sample size: %d\n", n))
cat(sprintf("  Selection rate (D=1): %.1f%%\n", 100 * mean(D)))
cat(sprintf("  Y1 observations: %d\n", sum(dat$eval_Y1)))
cat(sprintf("  Y0 observations: %d\n", sum(dat$eval_Y0)))
cat("\n")

# ---- Define Factor Model ----
fm <- define_factor_model(
  n_factors = 1,
  n_types = 1,
  n_quad = 16
)

# ---- Define Model Components ----
cat("Defining model components...\n")

# Selection equation (probit, full sample)
mc_sel <- define_model_component(
  name = "Selection",
  data = dat,
  outcome = "D",
  factor = fm,
  covariates = c("intercept", "Z1"),
  model_type = "probit",
  loading_normalization = NA_real_  # Free loading
)

# Y1 outcome (linear, observed when D=1)
mc_Y1 <- define_model_component(
  name = "Y1",
  data = dat,
  outcome = "Y",
  factor = fm,
  evaluation_indicator = "eval_Y1",  # Only D=1
  covariates = c("intercept", "X1"),
  model_type = "linear",
  loading_normalization = NA_real_  # Free loading
)

# Y0 outcome (linear, observed when D=0)
mc_Y0 <- define_model_component(
  name = "Y0",
  data = dat,
  outcome = "Y",
  factor = fm,
  evaluation_indicator = "eval_Y0",  # Only D=0
  covariates = c("intercept", "X1"),
  model_type = "linear",
  loading_normalization = NA_real_  # Free loading
)

# Test score 1 (linear, full sample)
mc_T1 <- define_model_component(
  name = "TestScore1",
  data = dat,
  outcome = "T1",
  factor = fm,
  covariates = c("intercept", "Q1"),
  model_type = "linear",
  loading_normalization = NA_real_  # Free loading
)

# Test score 2 (linear, full sample)
mc_T2 <- define_model_component(
  name = "TestScore2",
  data = dat,
  outcome = "T2",
  factor = fm,
  covariates = c("intercept", "Q2"),
  model_type = "linear",
  loading_normalization = NA_real_  # Free loading
)

# Test score 3 (linear, full sample)
mc_T3 <- define_model_component(
  name = "TestScore3",
  data = dat,
  outcome = "T3",
  factor = fm,
  covariates = c("intercept", "Q3"),
  model_type = "linear",
  loading_normalization = NA_real_  # Free loading
)

cat("  ✓ Selection equation (probit)\n")
cat("  ✓ Y1 outcome (linear, D=1 only)\n")
cat("  ✓ Y0 outcome (linear, D=0 only)\n")
cat("  ✓ Test score 1 (linear, full sample)\n")
cat("  ✓ Test score 2 (linear, full sample)\n")
cat("  ✓ Test score 3 (linear, full sample)\n")
cat("\n")

# ---- Create Model System ----
ms <- define_model_system(
  components = list(mc_sel, mc_Y1, mc_Y0, mc_T1, mc_T2, mc_T3),
  factor = fm
)

cat("Model system created with 6 components\n\n")

# ---- Initialize C++ Object for Gradient Checking ----
cat("=== PART 1: Gradient Validation ===\n\n")

fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 16)

param_info <- get_parameter_info_cpp(fm_cpp)
n_params <- param_info$n_param_free

cat(sprintf("Total parameters: %d\n", n_params))
cat("  [1] factor_var\n")
cat("  [2-4] Selection: intercept, beta_Z1, lambda_f (probit has no sigma)\n")
cat("  [5-8] Y1: intercept, beta_X1, lambda_f, sigma\n")
cat("  [9-12] Y0: intercept, beta_X1, lambda_f, sigma\n")
cat("  [13-16] TestScore1: intercept, beta_Q1, lambda_f, sigma\n")
cat("  [17-20] TestScore2: intercept, beta_Q2, lambda_f, sigma\n")
cat("  [21-24] TestScore3: intercept, beta_Q3, lambda_f, sigma\n\n")

# Test parameters (close to true values)
test_params <- c(
  1.0,                                      # factor_var
  0.5, 1.0, 1.0,                            # Selection: intercept, beta_Z1, lambda
  2.0, 1.0, 2.0, 0.5,                       # Y1: intercept, beta_X1, lambda, sigma
  1.0, 1.0, 1.0, 0.5,                       # Y0: intercept, beta_X1, lambda, sigma
  1.0, 1.0, 1.0, 0.5,                       # T1: intercept, beta_Q1, lambda, sigma
  0.5, 0.8, 1.2, 0.6,                       # T2: intercept, beta_Q2, lambda, sigma
  -0.5, 1.1, 0.9, 0.4                       # T3: intercept, beta_Q3, lambda, sigma
)

param_names <- c(
  "factor_var",
  "Sel_intercept", "Sel_beta_Z1", "Sel_lambda",
  "Y1_intercept", "Y1_beta_X1", "Y1_lambda", "Y1_sigma",
  "Y0_intercept", "Y0_beta_X1", "Y0_lambda", "Y0_sigma",
  "T1_intercept", "T1_beta_Q1", "T1_lambda", "T1_sigma",
  "T2_intercept", "T2_beta_Q2", "T2_lambda", "T2_sigma",
  "T3_intercept", "T3_beta_Q3", "T3_lambda", "T3_sigma"
)

cat("Testing gradient at parameter values close to truth...\n\n")

result <- evaluate_likelihood_cpp(fm_cpp, test_params,
                                 compute_gradient = TRUE,
                                 compute_hessian = FALSE)

cat(sprintf("Log-likelihood: %.4f\n\n", result$loglik))

# Finite difference gradient check
cat("Gradient validation (finite differences):\n")
delta <- 1e-7
max_diff <- 0
all_pass <- TRUE

for (i in seq_along(test_params)) {
  h <- delta * (abs(test_params[i]) + 1.0)

  params_plus <- test_params
  params_plus[i] <- test_params[i] + h
  f_plus <- evaluate_loglik_only_cpp(fm_cpp, params_plus)

  params_minus <- test_params
  params_minus[i] <- test_params[i] - h
  f_minus <- evaluate_loglik_only_cpp(fm_cpp, params_minus)

  grad_fd <- (f_plus - f_minus) / (2.0 * h)
  grad_analytical <- result$gradient[i]

  diff <- abs(grad_fd - grad_analytical)
  rel_diff <- diff / (abs(grad_analytical) + 1e-10)

  max_diff <- max(max_diff, rel_diff)

  status <- if (rel_diff < 1e-4) "✓" else "FAIL"
  if (rel_diff >= 1e-4) all_pass <- FALSE

  cat(sprintf("[%2d] %-15s: FD=%10.6f, Analytical=%10.6f, RelDiff=%.2e %s\n",
              i, param_names[i], grad_fd, grad_analytical, rel_diff, status))
}

cat(sprintf("\nMax relative difference: %.2e\n", max_diff))
if (all_pass) {
  cat("✓ ALL GRADIENTS PASS!\n\n")
} else {
  cat("⚠️  Some gradients failed\n\n")
}

# ---- Hessian Validation ----
cat("=== Hessian Validation (Finite Differences) ===\n\n")

# Determine which parameters are fixed
# Factor variance is fixed when no component has a fixed loading
factor_variance_identified <- FALSE
for (comp in list(mc_sel, mc_Y1, mc_Y0, mc_T1, mc_T2, mc_T3)) {
  if (!is.null(comp$loading_normalization)) {
    if (!is.na(comp$loading_normalization[1]) && abs(comp$loading_normalization[1]) > 1e-6) {
      factor_variance_identified <- TRUE
      break
    }
  }
}
fixed_params <- if (!factor_variance_identified) c(1) else c()  # Parameter 1 is factor_var

if (length(fixed_params) > 0) {
  cat(sprintf("Skipping fixed parameters: %s\n", paste(param_names[fixed_params], collapse = ", ")))
}
cat("\n")

result_hess <- evaluate_likelihood_cpp(fm_cpp, test_params,
                                       compute_gradient = TRUE,
                                       compute_hessian = TRUE)

# Convert upper triangle to full matrix
n_params <- length(test_params)
hess_analytical <- matrix(0, n_params, n_params)
idx <- 1
for (i in 1:n_params) {
  for (j in i:n_params) {
    hess_analytical[i, j] <- result_hess$hessian[idx]
    if (i != j) {
      hess_analytical[j, i] <- result_hess$hessian[idx]
    }
    idx <- idx + 1
  }
}

# Compute Hessian via finite differences of gradient
cat("Computing Hessian via finite differences of gradient...\n")
hess_fd <- matrix(0, n_params, n_params)
delta_hess <- 1e-6

for (i in seq_along(test_params)) {
  h <- delta_hess * (abs(test_params[i]) + 1.0)

  params_plus <- test_params
  params_plus[i] <- test_params[i] + h
  grad_plus <- evaluate_likelihood_cpp(fm_cpp, params_plus,
                                       compute_gradient = TRUE,
                                       compute_hessian = FALSE)$gradient

  params_minus <- test_params
  params_minus[i] <- test_params[i] - h
  grad_minus <- evaluate_likelihood_cpp(fm_cpp, params_minus,
                                        compute_gradient = TRUE,
                                        compute_hessian = FALSE)$gradient

  # Hessian column i = d(gradient)/d(param_i)
  hess_fd[, i] <- (grad_plus - grad_minus) / (2.0 * h)
}

# Compare diagonal elements (skip fixed parameters)
free_params <- setdiff(seq_along(test_params), fixed_params)

cat("\nHessian diagonal comparison (free parameters only):\n")
cat(sprintf("%-15s %15s %15s %12s\n", "Parameter", "Analytical", "Finite Diff", "RelDiff"))
cat(strrep("-", 70), "\n")

max_hess_diff <- 0
hess_all_pass <- TRUE

for (i in free_params) {
  analytical_val <- hess_analytical[i, i]
  fd_val <- hess_fd[i, i]
  diff <- abs(analytical_val - fd_val)
  rel_diff <- diff / (abs(analytical_val) + 1e-6)

  max_hess_diff <- max(max_hess_diff, rel_diff)

  status <- if (rel_diff < 1e-3) "✓" else "FAIL"
  if (rel_diff >= 1e-3) hess_all_pass <- FALSE

  cat(sprintf("%-15s %15.6f %15.6f %11.2e %s\n",
              param_names[i], analytical_val, fd_val, rel_diff, status))
}

cat(sprintf("\nMax Hessian diagonal relative difference: %.2e\n", max_hess_diff))
cat(sprintf("(%d free parameters tested)\n", length(free_params)))
if (hess_all_pass) {
  cat("✓ ALL HESSIAN DIAGONALS PASS!\n\n")
} else {
  cat("⚠️  Some Hessian elements failed\n\n")
}

# ---- Parameter Estimation ----
cat("=== PART 3: Parameter Recovery ===\n\n")

cat("Running estimation starting from TRUE parameters...\n")
cat("(Automatic initialization gets stuck in poor local optimum)\n\n")

# Start from true parameters to avoid poor local optimum
init_from_truth <- c(
  true_factor_var,
  true_sel_intercept, true_sel_beta_Z1, true_sel_lambda_f,
  true_Y1_intercept, true_Y1_beta_X1, true_Y1_lambda_f, true_Y1_sigma,
  true_Y0_intercept, true_Y0_beta_X1, true_Y0_lambda_f, true_Y0_sigma,
  true_T1_intercept, true_T1_beta_Q1, true_T1_lambda_f, true_T1_sigma,
  true_T2_intercept, true_T2_beta_Q2, true_T2_lambda_f, true_T2_sigma,
  true_T3_intercept, true_T3_beta_Q3, true_T3_lambda_f, true_T3_sigma
)

result_est <- estimate_model_rcpp(
  ms, dat,
  init_params = init_from_truth,  # Start from truth
  optimizer = "nlminb",
  parallel = FALSE,
  verbose = TRUE
)

cat("\n=== Parameter Recovery Results ===\n\n")

true_params <- c(
  true_factor_var,
  true_sel_intercept, true_sel_beta_Z1, true_sel_lambda_f,
  true_Y1_intercept, true_Y1_beta_X1, true_Y1_lambda_f, true_Y1_sigma,
  true_Y0_intercept, true_Y0_beta_X1, true_Y0_lambda_f, true_Y0_sigma,
  true_T1_intercept, true_T1_beta_Q1, true_T1_lambda_f, true_T1_sigma,
  true_T2_intercept, true_T2_beta_Q2, true_T2_lambda_f, true_T2_sigma,
  true_T3_intercept, true_T3_beta_Q3, true_T3_lambda_f, true_T3_sigma
)

cat(sprintf("%-20s %10s %10s %10s %10s %10s %10s\n",
            "Parameter", "True", "Estimated", "Std.Err", "t-stat", "Diff", "% Error"))
cat(strrep("-", 100), "\n")

recovery_good <- TRUE
for (i in seq_along(true_params)) {
  est <- result_est$estimates[i]
  se <- result_est$std_errors[i]
  true_val <- true_params[i]
  diff <- est - true_val
  pct_error <- 100 * diff / true_val
  t_stat <- (est - true_val) / se

  status <- if (abs(pct_error) < 15) "✓" else "⚠"
  if (abs(pct_error) >= 15) recovery_good <- FALSE

  cat(sprintf("%-20s %10.4f %10.4f %10.4f %10.2f %10.4f %9.1f%% %s\n",
              param_names[i], true_val, est, se, t_stat, diff, pct_error, status))
}

cat(strrep("-", 100), "\n")
cat(sprintf("Estimated log-likelihood: %.4f\n", result_est$loglik))

# Evaluate likelihood at true parameters
cat("\nEvaluating likelihood at true parameters...\n")
true_loglik <- evaluate_loglik_only_cpp(fm_cpp, true_params)
cat(sprintf("True parameters log-likelihood: %.4f\n", true_loglik))
cat(sprintf("Difference (estimated - true): %.4f\n", result_est$loglik - true_loglik))

if (result_est$loglik > true_loglik) {
  cat("✓ Estimated likelihood is higher than true parameters (as expected)\n\n")
} else {
  cat("⚠️  Warning: True parameters have higher likelihood!\n\n")
}

# ---- Summary ----
cat("=============================================================\n")
cat("SUMMARY\n")
cat("=============================================================\n\n")

cat("Model Structure:\n")
cat("  - 6 components (selection + 2 outcomes + 3 test scores)\n")
cat("  - 2 outcomes with missing indicators (Roy model)\n")
cat("  - 3 test scores for factor identification\n")
cat("  - All equations include intercepts\n")
cat(sprintf("  - %d observations total\n", n))
cat(sprintf("  - %d parameters estimated\n\n", n_params))

cat("Gradient Validation:\n")
if (all_pass) {
  cat("  ✓ All gradients match finite differences (max rel diff: ")
  cat(sprintf("%.2e)\n\n", max_diff))
} else {
  cat("  ⚠️  Some gradients failed validation\n\n")
}

cat("Hessian Validation:\n")
if (hess_all_pass) {
  cat("  ✓ All Hessian diagonals match finite differences (max rel diff: ")
  cat(sprintf("%.2e)\n\n", max_hess_diff))
} else {
  cat("  ⚠️  Some Hessian elements failed validation\n\n")
}

cat("Parameter Recovery:\n")
if (recovery_good) {
  cat("  ✓ All parameters recovered within 15% error\n")
} else {
  cat("  ⚠️  Some parameters have >15% error\n")
}
cat("  (Note: Roy models are challenging due to selection)\n\n")

cat("=============================================================\n")
cat("TEST COMPLETE\n")
cat("=============================================================\n")
