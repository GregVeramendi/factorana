#!/usr/bin/env Rscript
# Debug: what is TModel returning in modHess for sigma?

library(factorana)

# Create minimal test: 1 test score, 1 factor, fixed loading
set.seed(123)
n <- 100

f <- rnorm(n, 0, 1)
T1 <- 1.0 + 1.0 * f + rnorm(n, 0, 0.5)  # intercept=1, lambda=1 (fixed), sigma=0.5

dat <- data.frame(intercept = 1, T1 = T1, eval = 1)

fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 1)  # nquad=1 to simplify

mc_T1 <- define_model_component(
  name = "T1", data = dat, outcome = "T1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = 1.0, evaluation_indicator = "eval"
)

ms <- define_model_system(components = list(mc_T1), factor = fm)
fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 1)

# Parameters: factor_var, intercept, sigma (lambda is fixed)
params <- c(1.0, 1.0, 0.5)

result <- evaluate_likelihood_cpp(fm_cpp, params, TRUE, TRUE)

cat("Parameters: factor_var=1.0, intercept=1.0, sigma=0.5\n")
cat("nquad = 1 (single integration point at mean)\n\n")

cat("Gradient:\n")
print(result$gradient)

cat("\nHessian (upper triangle):\n")
print(result$hessian)

cat("\nHessian as 3x3 matrix:\n")
hess <- matrix(0, 3, 3)
idx <- 1
for (i in 1:3) {
  for (j in i:3) {
    hess[i,j] <- result$hessian[idx]
    hess[j,i] <- result$hessian[idx]
    idx <- idx + 1
  }
}
print(hess)

cat("\nSigma diagonal (should be sum over obs of: -3*Z^2/sigma^4 + 1/sigma^2):\n")
cat("where Z = Y - intercept - lambda*f\n")
cat("At f=0 (single quad point), Z = Y - 1.0\n\n")

Y <- dat$T1
Z <- Y - 1.0
sigma <- 0.5
manual_sigma_hess <- sum(-3*Z^2/sigma^4 + 1/sigma^2)

cat(sprintf("Manual calculation: %.6f\n", manual_sigma_hess))
cat(sprintf("C++ returned:       %.6f\n", hess[3,3]))
cat(sprintf("Match: %s\n", abs(manual_sigma_hess - hess[3,3]) < 1e-6))

# Now finite difference
delta <- 1e-6
params_plus <- params
params_plus[3] <- params[3] + delta
grad_plus <- evaluate_likelihood_cpp(fm_cpp, params_plus, TRUE, FALSE)$gradient[3]

params_minus <- params
params_minus[3] <- params[3] - delta  
grad_minus <- evaluate_likelihood_cpp(fm_cpp, params_minus, TRUE, FALSE)$gradient[3]

hess_fd <- (grad_plus - grad_minus) / (2 * delta)

cat(sprintf("\nFinite difference: %.6f\n", hess_fd))
cat(sprintf("Ratio (C++/FD):    %.6f\n", hess[3,3] / hess_fd))
