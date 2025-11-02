#!/usr/bin/env Rscript

# ---- Roy one-factor pipeline: simple smoke script (no testthat) ----

# tiny assertion helper
assert_true <- function(cond, msg = "assertion failed") {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

# Source your package functions if they exist (don’t source tests)
to_source <- c(
  "define_model_component.R",
  "define_model_system.R",
  "define_estimation_control.R",
  "define_factor_model.R",
  "estimate_model.R",
  "initialize_parameters.R"
)
for (f in to_source) {
  p <- file.path("R", f)
  if (file.exists(p)) source(p, chdir = TRUE)
}

# If define_factor_model() isn’t available yet, provide a minimal one
if (!exists("define_factor_model", mode = "function")) {
  define_factor_model <- function(nfactors = 1L, n_types = 1L, gh_nodes = 8L) {
    structure(list(
      nfactors = as.integer(nfactors),
      n_types  = as.integer(n_types),
      gh_nodes = as.integer(gh_nodes)
    ), class = "factor_model")
  }
}

set.seed(123)

# ---- Generate Fake Data (homework D) ----
n  <- 5000
f  <- rnorm(n, 0, 1)
U1 <- rnorm(n, 0, 1)
U0 <- rnorm(n, 0, 1)
UV <- rnorm(n, 0, 1)
UT <- rnorm(n, 0, 1)

Z0 <- X0 <- Q0 <- 1   # intercepts
Z1 <- rnorm(n, 0, 2)
X1 <- rnorm(n, 0, 2)
Q1 <- rnorm(n, 0, 1)

Y1 <- 2*Z0 + X1 + 2*f + U1
Y0 <- 1*Z0 + X1 + 1*f + U0
I  <- 0.5*Z0 + Z1 + f + UV
D  <- as.integer(I > 0)
Tscore <- Q0 + Q1 + f + UT
Y  <- ifelse(D == 1L, Y1, Y0)

dat <- data.frame(Y, D, Z1, X1, Q1, T = Tscore)
dat$Z0 <- 1; dat$X0 <- 1; dat$Q0 <- 1
dat$eval_y0 <- 1L - dat$D

# ---- Factor model ----
factor_1 <- define_factor_model(1, 1, 8)
assert_true(is.list(factor_1), "factor_1 must be a list")
assert_true(all(c("nfactors","n_types") %in% names(factor_1)),
            "factor_1 must include nfactors and n_types")

# ---- Components ----
mc_sel <- define_model_component("selection", dat, "D", factor_1,
                                 covariates = c("Z0","Z1"), model_type = "probit")
mc_y1  <- define_model_component("Y1", dat, "Y", factor_1,
                                 evaluation_indicator = "D",
                                 covariates = c("X0","X1"), model_type = "linear")
mc_y0  <- define_model_component("Y0", dat, "Y", factor_1,
                                 evaluation_indicator = "eval_y0",
                                 covariates = c("X0","X1"), model_type = "linear")
mc_T   <- define_model_component("T", dat, "T", factor_1,
                                 covariates = c("Q0","Q1"), model_type = "linear")

# ---- System & control ----
ms   <- define_model_system(factor = factor_1, components = list(mc_sel, mc_y1, mc_y0, mc_T))
ctrl <- define_estimation_control(num_cores = 1)

# ---- Initialize + basic checks ----
inits <- lapply(ms$components, initialize_parameters)
assert_true(length(inits) == 4, "need 4 components/inits")
for (ini in inits) {
  must <- c("intercept","betas","loading","factor_var","factor_cor")
  assert_true(all(must %in% names(ini)), "initializer missing fields")
  assert_true(is.finite(ini$intercept) && is.finite(ini$loading), "bad init values")
}

# ---- Estimate model & write CSV (if your estimator doesn’t) ----
out <- estimate_model(ms, ctrl)
assert_true(is.data.frame(out), "estimate_model() must return a data.frame")

if (!dir.exists("results")) dir.create("results", recursive = TRUE)
csv_path <- file.path("results", paste0("system_inits_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
try(utils::write.csv(out, csv_path, row.names = FALSE), silent = TRUE)

if (file.exists(csv_path)) {
  cat("✅ CSV written:", csv_path, "\n")
} else {
  cat("ℹ️ No CSV written by script; estimator may handle CSV itself.\n")
}

# ---- Sanity stats ----
m1 <- mean(dat$Y[dat$D == 1]); m0 <- mean(dat$Y[dat$D == 0])
sel <- mean(dat$D)
cat(sprintf("Y mean (treated): %.3f\n", m1))
cat(sprintf("Y mean (control): %.3f\n", m0))
cat(sprintf("Selection rate  : %.3f\n", sel))
cat("✅ Smoke run finished without errors.\n")
