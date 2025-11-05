# Debug script for two-factor model convergence issues
library(factorana)
library(MASS)

set.seed(123)
n <- 500

# Two independent factors
f1 <- rnorm(n, 0, 1)
f2 <- rnorm(n, 0, 1)

# True loadings for factor 1 measures (m1, m2, m3)
lambda1_1 <- 1.0   # Fixed for identification
lambda1_2 <- 0.8
lambda1_3 <- 1.2

# True loadings for factor 2 measures (m4, m5, m6)
lambda2_4 <- 1.0   # Fixed for identification
lambda2_5 <- 0.9
lambda2_6 <- 1.1

# Generate latent continuous variables
y1_star <- lambda1_1 * f1 + rnorm(n, 0, 0.5)
y2_star <- lambda1_2 * f1 + rnorm(n, 0, 0.5)
y3_star <- lambda1_3 * f1 + rnorm(n, 0, 0.5)
y4_star <- lambda2_4 * f2 + rnorm(n, 0, 0.5)
y5_star <- lambda2_5 * f2 + rnorm(n, 0, 0.5)
y6_star <- lambda2_6 * f2 + rnorm(n, 0, 0.5)

# Convert to ordered categories (5 categories each)
make_ordered <- function(y_star) {
  as.integer(cut(y_star,
                 breaks = c(-Inf, -1, -0.3, 0.3, 1, Inf),
                 labels = FALSE))
}

dat <- data.frame(
  m1 = make_ordered(y1_star),
  m2 = make_ordered(y2_star),
  m3 = make_ordered(y3_star),
  m4 = make_ordered(y4_star),
  m5 = make_ordered(y5_star),
  m6 = make_ordered(y6_star)
)

cat("Data summary:\n")
print(summary(dat))
cat("\n\n")

# Define 2-factor model
fm <- define_factor_model(n_factors = 2, n_types = 1, n_quad_points = 16)

# Factor 1 measures: loading on f1, zero on f2
cat("Defining model components...\n")

mc1 <- define_model_component(
  name = "m1", data = dat, outcome = "m1", factor = fm,
  covariates = NULL,
  model_type = "oprobit",
  loading_normalization = c(1.0, 0),  # f1=1.0 (fixed), f2=0 (zero)
  num_choices = 5
)

mc2 <- define_model_component(
  name = "m2", data = dat, outcome = "m2", factor = fm,
  covariates = NULL,
  model_type = "oprobit",
  loading_normalization = c(NA, 0),  # f1=free, f2=0 (zero)
  num_choices = 5
)

mc3 <- define_model_component(
  name = "m3", data = dat, outcome = "m3", factor = fm,
  covariates = NULL,
  model_type = "oprobit",
  loading_normalization = c(NA, 0),  # f1=free, f2=0 (zero)
  num_choices = 5
)

mc4 <- define_model_component(
  name = "m4", data = dat, outcome = "m4", factor = fm,
  covariates = NULL,
  model_type = "oprobit",
  loading_normalization = c(0, 1.0),  # f1=0 (zero), f2=1.0 (fixed)
  num_choices = 5
)

mc5 <- define_model_component(
  name = "m5", data = dat, outcome = "m5", factor = fm,
  covariates = NULL,
  model_type = "oprobit",
  loading_normalization = c(0, NA),  # f1=0 (zero), f2=free
  num_choices = 5
)

mc6 <- define_model_component(
  name = "m6", data = dat, outcome = "m6", factor = fm,
  covariates = NULL,
  model_type = "oprobit",
  loading_normalization = c(0, NA),  # f1=0 (zero), f2=free
  num_choices = 5
)

# Create model system
cat("Creating model system...\n")
ms <- define_model_system(
  components = list(mc1, mc2, mc3, mc4, mc5, mc6),
  factor = fm
)

# Initialize parameters
cat("\nInitializing parameters...\n")
init_result <- initialize_parameters(ms, dat, verbose = TRUE)
cat("\nInitial parameter vector:\n")
print(init_result$init_params)
cat("\nLength of initial parameters:", length(init_result$init_params), "\n")
cat("Factor variance fixed:", init_result$factor_variance_fixed, "\n\n")

# Estimate the model
cat("\nStarting estimation...\n")
ctrl <- define_estimation_control(num_cores = 1)
result <- estimate_model_rcpp(
  model_system = ms,
  data = dat,
  init_params = init_result$init_params,
  control = ctrl,
  parallel = FALSE,
  optimizer = "nlminb",
  verbose = TRUE
)

cat("\n\n=== RESULTS ===\n")
cat("Convergence code:", result$convergence, "\n")
cat("Log-likelihood:", result$loglik, "\n")
cat("Number of parameters:", length(result$estimates), "\n")
cat("\nParameter names:\n")
print(names(result$estimates))
cat("\nEstimates:\n")
print(result$estimates)
cat("\nStandard errors:\n")
print(result$std_errors)
