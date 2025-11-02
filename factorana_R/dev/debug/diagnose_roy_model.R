#!/usr/bin/env Rscript
# Diagnostic: Check gradient and Hessian at TRUE parameters

library(factorana)

cat("=============================================================\n")
cat("Roy Model Diagnostic: Gradient & Hessian at True Parameters\n")
cat("=============================================================\n\n")

# ---- Generate Data with Known Parameters ----
set.seed(456)
n <- 2000

# Latent factor
f <- rnorm(n, 0, 1)

# Error terms
U1 <- rnorm(n, 0, 0.5)  # Y1 error
U0 <- rnorm(n, 0, 0.5)  # Y0 error
UV <- rnorm(n, 0, 0.5)  # Selection error
UT1 <- rnorm(n, 0, 0.5)
UT2 <- rnorm(n, 0, 0.5)
UT3 <- rnorm(n, 0, 0.5)

# Covariates
Z0 <- rep(1, n)
Z1 <- rnorm(n, 0, 1)
X0 <- rep(1, n)
X1 <- rnorm(n, 0, 1)
Q0 <- rep(1, n)
Q1 <- rnorm(n, 0, 1)
Q2 <- rnorm(n, 0, 1)
Q3 <- rnorm(n, 0, 1)

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

# Generate outcomes
Y1 <- true_Y1_intercept * X0 + true_Y1_beta_X1 * X1 + true_Y1_lambda_f * f + U1
Y0 <- true_Y0_intercept * X0 + true_Y0_beta_X1 * X1 + true_Y0_lambda_f * f + U0

# Selection (latent utility)
I_latent <- true_sel_intercept * Z0 + true_sel_beta_Z1 * Z1 + true_sel_lambda_f * f + UV
D <- as.integer(I_latent > 0)

# Test scores
T1_score <- true_T1_intercept * Q0 + true_T1_beta_Q1 * Q1 + true_T1_lambda_f * f + UT1
T2_score <- true_T2_intercept * Q0 + true_T2_beta_Q2 * Q2 + true_T2_lambda_f * f + UT2
T3_score <- true_T3_intercept * Q0 + true_T3_beta_Q3 * Q3 + true_T3_lambda_f * f + UT3

# Observed outcome (Roy model)
Y <- ifelse(D == 1, Y1, Y0)

# Create dataset
dat <- data.frame(
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
dat$eval_Y1 <- dat$D
dat$eval_Y0 <- 1L - dat$D

cat(sprintf("Sample size: %d\n", n))
cat(sprintf("Selection rate (D=1): %.1f%%\n\n", 100 * mean(D)))

# ---- Define Factor Model ----
fm <- define_factor_model(
  n_factors = 1,
  n_types = 1,
  n_quad = 8
)

# ---- Define Model Components ----
mc_sel <- define_model_component(
  name = "Selection",
  data = dat,
  outcome = "D",
  factor = fm,
  covariates = "Z1",
  model_type = "probit",
  loading_normalization = NA_real_
)

mc_Y1 <- define_model_component(
  name = "Y1",
  data = dat,
  outcome = "Y",
  factor = fm,
  evaluation_indicator = "eval_Y1",
  covariates = "X1",
  model_type = "linear",
  loading_normalization = NA_real_
)

mc_Y0 <- define_model_component(
  name = "Y0",
  data = dat,
  outcome = "Y",
  factor = fm,
  evaluation_indicator = "eval_Y0",
  covariates = "X1",
  model_type = "linear",
  loading_normalization = NA_real_
)

mc_T1 <- define_model_component(
  name = "TestScore1",
  data = dat,
  outcome = "T1",
  factor = fm,
  covariates = "Q1",
  model_type = "linear",
  loading_normalization = NA_real_
)

mc_T2 <- define_model_component(
  name = "TestScore2",
  data = dat,
  outcome = "T2",
  factor = fm,
  covariates = "Q2",
  model_type = "linear",
  loading_normalization = NA_real_
)

mc_T3 <- define_model_component(
  name = "TestScore3",
  data = dat,
  outcome = "T3",
  factor = fm,
  covariates = "Q3",
  model_type = "linear",
  loading_normalization = NA_real_
)

# ---- Create Model System ----
ms <- define_model_system(
  components = list(mc_sel, mc_Y1, mc_Y0, mc_T1, mc_T2, mc_T3),
  factor = fm
)

# ---- Initialize C++ Object ----
fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8)

# ---- TRUE PARAMETERS (those we used to generate data) ----
# Order: factor_var, Sel params, Y1 params, Y0 params, T1 params, T2 params, T3 params
true_params <- c(
  true_factor_var,
  true_sel_beta_Z1, true_sel_lambda_f,
  true_Y1_beta_X1, true_Y1_lambda_f, true_Y1_sigma,
  true_Y0_beta_X1, true_Y0_lambda_f, true_Y0_sigma,
  true_T1_beta_Q1, true_T1_lambda_f, true_T1_sigma,
  true_T2_beta_Q2, true_T2_lambda_f, true_T2_sigma,
  true_T3_beta_Q3, true_T3_lambda_f, true_T3_sigma
)

param_names <- c(
  "factor_var",
  "Sel_beta_Z1", "Sel_lambda",
  "Y1_beta_X1", "Y1_lambda", "Y1_sigma",
  "Y0_beta_X1", "Y0_lambda", "Y0_sigma",
  "T1_beta_Q1", "T1_lambda", "T1_sigma",
  "T2_beta_Q2", "T2_lambda", "T2_sigma",
  "T3_beta_Q3", "T3_lambda", "T3_sigma"
)

cat("=== Evaluating Gradient and Hessian at TRUE Parameters ===\n\n")

# Evaluate likelihood, gradient, and Hessian at true parameters
result <- evaluate_likelihood_cpp(fm_cpp, true_params,
                                 compute_gradient = TRUE,
                                 compute_hessian = TRUE)

cat(sprintf("Log-likelihood at true parameters: %.4f\n\n", result$loglik))

# ---- Check Gradient ----
cat("Gradient at true parameters:\n")
cat("(Should be close to zero if true params are a local optimum)\n\n")

max_grad <- max(abs(result$gradient))
cat(sprintf("%-20s %15s\n", "Parameter", "Gradient"))
cat(strrep("-", 40), "\n")

for (i in seq_along(true_params)) {
  status <- if (abs(result$gradient[i]) < 1e-2) "✓" else "⚠"
  cat(sprintf("%-20s %15.6f %s\n", param_names[i], result$gradient[i], status))
}

cat(strrep("-", 40), "\n")
cat(sprintf("Max absolute gradient: %.6f\n\n", max_grad))

if (max_grad < 1e-2) {
  cat("✓ Gradients are very close to zero (true params are near a stationary point)\n\n")
} else {
  cat("⚠️  Gradients are NOT close to zero!\n")
  cat("This suggests the true parameters are NOT a local maximum.\n")
  cat("Possible issues:\n")
  cat("  1. Model specification doesn't match data generation\n")
  cat("  2. Missing intercepts in some equations\n")
  cat("  3. Wrong parameter ordering\n\n")
}

# ---- Check Hessian ----
cat("=== Hessian Analysis at True Parameters ===\n\n")

# Convert upper triangle to full matrix
n_params <- length(true_params)
hess_matrix <- matrix(0, n_params, n_params)
idx <- 1
for (i in 1:n_params) {
  for (j in i:n_params) {
    hess_matrix[i, j] <- result$hessian[idx]
    if (i != j) {
      hess_matrix[j, i] <- result$hessian[idx]
    }
    idx <- idx + 1
  }
}

# Negate because C++ returns Hessian of loglik, but we minimize -loglik
hess_matrix <- -hess_matrix

# Eigenanalysis
eig_result <- eigen(hess_matrix, symmetric = TRUE)
eig_vals <- eig_result$values

cat("Hessian eigenvalues (should all be positive for local minimum):\n")
cat(sprintf("  Min eigenvalue: %.4e\n", min(eig_vals)))
cat(sprintf("  Max eigenvalue: %.4e\n", max(eig_vals)))
cat(sprintf("  Condition number: %.4e\n\n", max(abs(eig_vals)) / min(abs(eig_vals))))

n_negative <- sum(eig_vals < 0)
n_zero <- sum(abs(eig_vals) < 1e-6)
n_positive <- sum(eig_vals > 1e-6)

cat(sprintf("Eigenvalue counts: %d negative, %d ~zero, %d positive\n\n",
            n_negative, n_zero, n_positive))

if (n_negative > 0) {
  cat("⚠️  Hessian has NEGATIVE eigenvalues!\n")
  cat("This means true parameters are NOT at a local minimum.\n")
  cat("The negative log-likelihood surface is not concave at this point.\n\n")

  cat("Negative eigenvalues:\n")
  neg_eigs <- eig_vals[eig_vals < 0]
  for (i in seq_along(neg_eigs)) {
    cat(sprintf("  [%d] %.4e\n", i, neg_eigs[i]))
  }
  cat("\n")
}

# Show Hessian diagonal
cat("Hessian diagonal elements:\n")
cat(sprintf("%-20s %15s\n", "Parameter", "Hess[i,i]"))
cat(strrep("-", 40), "\n")
for (i in seq_along(true_params)) {
  status <- if (hess_matrix[i, i] > 0) "✓" else "⚠"
  cat(sprintf("%-20s %15.6f %s\n", param_names[i], hess_matrix[i, i], status))
}

cat("\n=== DIAGNOSIS ===\n\n")

if (max_grad > 1e-2) {
  cat("❌ PROBLEM IDENTIFIED: Gradient not zero at true parameters\n")
  cat("The model likelihood is NOT maximized at the true parameters.\n")
  cat("This indicates a MISMATCH between:\n")
  cat("  - How the data was generated\n")
  cat("  - How the model is specified\n\n")
  cat("Likely causes:\n")
  cat("  1. MISSING INTERCEPTS: Check if all equations include intercepts\n")
  cat("  2. Parameter ordering mismatch\n")
  cat("  3. Different normalization assumptions\n\n")
} else if (n_negative > 0) {
  cat("⚠️  Gradient is zero but Hessian has negative eigenvalues\n")
  cat("True parameters are at a saddle point or local maximum (not minimum).\n\n")
} else {
  cat("✓ True parameters appear to be a valid local minimum\n")
  cat("The optimization should converge near these values.\n\n")
}

cat("=============================================================\n")
