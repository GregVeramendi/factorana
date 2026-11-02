# Multinomial logit with latent types: one type intercept per non-reference alternative.
# Before 1.8.0 Model::Eval() received a single type intercept per (model, type) and added it to every
# alternative, so the typed likelihood of a logit with more than two alternatives (and its gradient and
# Hessian) used the last alternative's intercept everywhere. Found against the legacy C++ code (T9).
#
# Design: one factor measured by three linear indicators without types; two types whose probability
# depends on the factor; type intercepts in two linear outcomes, one probit and one 3-alternative logit.
# Several typed components are needed for the type model to be identified.

.mnl_types_data <- function(n = 1500, seed = 5) {
  set.seed(seed)
  f <- rnorm(n); x1 <- rnorm(n)
  type <- 1L + rbinom(n, 1, plogis(0.3 + 0.8 * f)); t2 <- as.numeric(type == 2)
  V2 <- 0.2 + 0.3 * x1 + 0.7 * f + 0.9 * t2
  V3 <- -0.3 + 0.5 * x1 + 1.1 * f - 0.7 * t2
  U <- cbind(0, V2, V3) - log(-log(matrix(runif(3 * n), n)))
  data.frame(intercept = 1, x1 = x1,
             y1 = 0.1 + 1.0 * f + rnorm(n, 0, 0.6),
             y2 = -0.2 + 0.8 * f + rnorm(n, 0, 0.7),
             y3 = 0.2 + 0.6 * f + rnorm(n, 0, 0.5),
             w1 = 0.5 + 0.2 * x1 + 0.6 * f + 1.0 * t2 + rnorm(n, 0, 0.7),
             w2 = -0.4 - 0.3 * x1 + 0.5 * f - 0.8 * t2 + rnorm(n, 0, 0.6),
             d = as.integer(-0.2 + 0.4 * x1 + 0.8 * f + 0.7 * t2 + rnorm(n) > 0),
             out = max.col(U, ties.method = "first"))
}

.mnl_types_truth <- c(factor_var_1 = 1, typeprob_2_intercept = 0.3, type_2_loading_1 = 0.8,
                      y1_intercept = 0.1, y1_sigma = 0.6,          # y1 loading fixed at 1, not a parameter
                      y2_intercept = -0.2, y2_loading_1 = 0.8, y2_sigma = 0.7,
                      y3_intercept = 0.2, y3_loading_1 = 0.6, y3_sigma = 0.5,
                      w1_intercept = 0.5, w1_x1 = 0.2, w1_loading_1 = 0.6, w1_sigma = 0.7, w1_type_2_intercept = 1.0,
                      w2_intercept = -0.4, w2_x1 = -0.3, w2_loading_1 = 0.5, w2_sigma = 0.6, w2_type_2_intercept = -0.8,
                      d_intercept = -0.2, d_x1 = 0.4, d_loading_1 = 0.8, d_type_2_intercept = 0.7,
                      out_c1_intercept = 0.2, out_c1_x1 = 0.3, out_c1_loading_1 = 0.7,
                      out_c2_intercept = -0.3, out_c2_x1 = 0.5, out_c2_loading_1 = 1.1,
                      out_c1_type_2_intercept = 0.9, out_c2_type_2_intercept = -0.7)

.mnl_types_system <- function(dat, fix_choice = NULL) {
  fm <- define_factor_model(n_factors = 1, n_types = 2)
  ind <- function(nm, norm) define_model_component(nm, dat, nm, fm, covariates = "intercept",
                                                   model_type = "linear", loading_normalization = norm)
  typed <- function(nm, type) define_model_component(nm, dat, nm, fm, covariates = c("intercept", "x1"),
                                                     model_type = type, loading_normalization = NA_real_, use_types = TRUE)
  out <- define_model_component("out", dat, "out", fm, covariates = c("intercept", "x1"), model_type = "logit",
                                num_choices = 3, loading_normalization = NA_real_, use_types = TRUE)
  if (!is.null(fix_choice)) out <- fix_type_intercepts(out, types = 2, choice = fix_choice)
  define_model_system(components = list(ind("y1", 1), ind("y2", NA_real_), ind("y3", NA_real_),
                                        typed("w1", "linear"), typed("w2", "linear"), typed("d", "probit"), out), factor = fm)
}

.mnl_types_setup <- function(ms, dat, p) {
  init <- initialize_parameters(ms, dat, verbose = FALSE)
  stopifnot(setequal(names(p), init$param_names))
  p <- p[init$param_names]
  meta <- factorana:::build_parameter_metadata(ms)
  cons <- factorana:::setup_parameter_constraints(ms, p, meta, init$factor_variance_fixed, verbose = FALSE)
  list(p = p, free_idx = cons$free_idx, fixed = !(seq_along(p) %in% cons$free_idx))
}

.mnl_types_loglik <- function(ms, dat, p) {
  s <- .mnl_types_setup(ms, dat, p)
  fm_ptr <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8L, init_params = s$p)
  evaluate_loglik_only_cpp(fm_ptr, s$p[s$free_idx])
}

test_that("each alternative's type intercept moves the typed logit likelihood on its own", {
  dat <- .mnl_types_data(); ms <- .mnl_types_system(dat)
  p0 <- .mnl_types_truth; p0[c("out_c1_type_2_intercept", "out_c2_type_2_intercept")] <- 0
  pa <- p0; pa["out_c1_type_2_intercept"] <- 0.9
  pb <- p0; pb["out_c2_type_2_intercept"] <- -0.7
  ll0 <- .mnl_types_loglik(ms, dat, p0); lla <- .mnl_types_loglik(ms, dat, pa); llb <- .mnl_types_loglik(ms, dat, pb)
  expect_gt(abs(lla - ll0), 1e-3)
  expect_gt(abs(llb - ll0), 1e-3)
  expect_gt(abs(lla - llb), 1e-3)
})

test_that("typed multinomial logit: analytical gradient and Hessian match finite differences", {
  dat <- .mnl_types_data(); ms <- .mnl_types_system(dat)
  s <- .mnl_types_setup(ms, dat, .mnl_types_truth)
  # Gradient: central differences with step 1e-5. The helper's default step, sqrt(eps) * max(|p|, 1), is
  # roundoff-dominated when |logL| is about 8,000 (errors near 1e-4 for every parameter, typed or not),
  # so the reference is computed here with a better-conditioned step; the 1e-5 criterion is unchanged.
  g_an <- evaluate_likelihood_rcpp(ms, dat, s$p, compute_gradient = TRUE, n_quad = 8)$gradient
  g_fd <- finite_diff_gradient(ms, dat, s$p, h = 1e-5, n_quad = 8)
  rel <- abs(g_an - g_fd) / pmax(abs(g_fd), 1e-10)
  expect_true(all(rel[s$free_idx] < 1e-5), info = sprintf("max gradient relative error %.2e", max(rel[s$free_idx])))
  h <- check_hessian_accuracy(ms, dat, s$p, param_fixed = s$fixed, tol = 1e-3, n_quad = 8)
  expect_true(h$pass, info = sprintf("max Hessian error %.2e", h$max_error))
})

test_that("typed multinomial logit recovers alternative-specific type intercepts", {
  skip_on_cran()
  dat <- .mnl_types_data(); ms <- .mnl_types_system(dat)
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  fit <- estimate_model_rcpp(ms, dat, control = ctrl, optimizer = "nlminb", parallel = FALSE, verbose = FALSE)
  expect_equal(fit$convergence, 0)
  typed <- grep("type", names(fit$estimates), value = TRUE)
  expect_true(all(fit$std_errors[typed] < 10))
  z <- (fit$estimates[typed] - .mnl_types_truth[typed]) / fit$std_errors[typed]
  expect_true(all(abs(z) < 3.5), info = paste(sprintf("%s z=%.2f", typed, z), collapse = "; "))
})

test_that("fix_type_intercepts(choice = 1) fixes only that alternative's type intercept", {
  skip_on_cran()
  dat <- .mnl_types_data(); ms <- .mnl_types_system(dat, fix_choice = 1)
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  fit <- estimate_model_rcpp(ms, dat, control = ctrl, optimizer = "nlminb", parallel = FALSE, verbose = FALSE)
  expect_equal(fit$convergence, 0)
  expect_equal(unname(fit$estimates["out_c1_type_2_intercept"]), 0)
  expect_gt(abs(fit$estimates["out_c2_type_2_intercept"]), 0.05)
})
