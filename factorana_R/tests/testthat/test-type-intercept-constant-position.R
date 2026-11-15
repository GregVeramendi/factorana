# The constant of a typed component can sit anywhere among its covariates.
#
# A type intercept shifts the linear predictor as the coefficient of the component's constant does,
# and FactorModel.cpp takes its gradient and Hessian from the slots of that coefficient. Up to 1.8.2
# it took them from the first covariate. Every typed test listed "intercept" first, so nothing showed
# that a model with another covariate first had a right likelihood and wrong derivatives for every
# type intercept. nlminb then stops without converging. The HJV model lists its constant twelfth.
#
# Design: one factor measured by three untyped linear indicators, two types whose probability depends
# on the factor, and a typed linear, probit and three-alternative logit component. The typed
# components take their covariates in the order given, c("x1", "intercept") or c("intercept", "x1").

.tic_data <- function(n = 1200, seed = 5) {
  set.seed(seed)
  f <- rnorm(n); x1 <- rnorm(n)
  t2 <- rbinom(n, 1, plogis(0.3 + 0.8 * f))
  V2 <- 0.2 + 0.3 * x1 + 0.7 * f + 0.9 * t2
  V3 <- -0.3 + 0.5 * x1 + 1.1 * f - 0.7 * t2
  U <- cbind(0, V2, V3) - log(-log(matrix(runif(3 * n), n)))
  data.frame(intercept = 1, x1 = x1,
             y1 = 0.1 + 1.0 * f + rnorm(n, 0, 0.6),
             y2 = -0.2 + 0.8 * f + rnorm(n, 0, 0.7),
             y3 = 0.2 + 0.6 * f + rnorm(n, 0, 0.5),
             w = 0.5 + 0.2 * x1 + 0.6 * f + 1.0 * t2 + rnorm(n, 0, 0.7),
             d = as.integer(-0.2 + 0.4 * x1 + 0.8 * f + 0.7 * t2 + rnorm(n) > 0),
             out = max.col(U, ties.method = "first"))
}

.tic_system <- function(dat, covs) {
  fm <- define_factor_model(n_factors = 1, n_types = 2)
  ind <- function(nm, norm) define_model_component(nm, dat, nm, fm, covariates = "intercept",
                                                   model_type = "linear", loading_normalization = norm)
  typed <- function(nm, type, ...) define_model_component(nm, dat, nm, fm, covariates = covs, model_type = type,
                                                          loading_normalization = NA_real_, use_types = TRUE, ...)
  define_model_system(components = list(ind("y1", 1), ind("y2", NA_real_), ind("y3", NA_real_),
                                        typed("w", "linear"), typed("d", "probit"),
                                        typed("out", "logit", num_choices = 3)), factor = fm)
}

# A point away from any optimum, the same for both covariate orders (values are set by name)
.tic_point <- function(ms, dat) {
  init <- initialize_parameters(ms, dat, verbose = FALSE)
  nm <- init$param_names
  p <- setNames(0.3 * sin(seq_along(sort(nm))), sort(nm))[nm]
  p[grepl("_sigma$", nm)] <- 0.7
  p["factor_var_1"] <- 1
  meta <- factorana:::build_parameter_metadata(ms)
  cons <- factorana:::setup_parameter_constraints(ms, p, meta, init$factor_variance_fixed, verbose = FALSE)
  list(p = p, free = nm[cons$free_idx], fixed = !(seq_along(p) %in% cons$free_idx))
}

.tic_unpack <- function(v, nm) {
  n <- length(nm)
  H <- matrix(0, n, n, dimnames = list(nm, nm))
  H[lower.tri(H, diag = TRUE)] <- v
  H <- t(H)
  H[lower.tri(H)] <- t(H)[lower.tri(H)]
  H
}

test_that("a typed component names its constant, wherever it is listed", {
  dat <- .tic_data(200)
  fm <- define_factor_model(n_factors = 1, n_types = 2)
  mc <- function(covs, typed) define_model_component("w", dat, "w", fm, covariates = covs, model_type = "linear",
                                                     loading_normalization = NA_real_, use_types = typed)
  expect_identical(mc(c("x1", "intercept"), TRUE)$type_intercept_covariate, "intercept")
  expect_identical(mc(c("intercept", "x1"), TRUE)$type_intercept_covariate, "intercept")
  expect_null(mc(c("x1", "intercept"), FALSE)$type_intercept_covariate)
})

test_that("a typed component without a constant is an error at definition", {
  dat <- .tic_data(200)
  fm <- define_factor_model(n_factors = 1, n_types = 2)
  expect_error(
    suppressWarnings(define_model_component("w", dat, "w", fm, covariates = "x1", model_type = "linear",
                                            loading_normalization = NA_real_, use_types = TRUE)),
    "uses types but has no constant")
})

test_that("constant listed second: analytical gradient and Hessian match finite differences", {
  skip_on_cran()
  dat <- .tic_data()
  ms <- .tic_system(dat, c("x1", "intercept"))
  s <- .tic_point(ms, dat)
  # Central differences with step 1e-5 (see test-types-mnl.R on the helper's default step)
  g_an <- setNames(evaluate_likelihood_rcpp(ms, dat, s$p, compute_gradient = TRUE, n_quad = 8)$gradient, names(s$p))
  g_fd <- setNames(finite_diff_gradient(ms, dat, s$p, h = 1e-5, n_quad = 8), names(s$p))
  err <- abs(g_an - g_fd) / pmax(abs(g_fd), 1)
  type_int <- grep("_type_2_intercept$", s$free, value = TRUE)
  expect_length(type_int, 4)                      # w, d, out_c1, out_c2
  expect_gt(min(abs(g_fd[type_int])), 0.1)        # the type intercepts are not at a stationary point
  expect_true(all(err[type_int] < 1e-5), info = sprintf("max type-intercept gradient error %.2e", max(err[type_int])))
  expect_true(all(err[s$free] < 1e-5), info = sprintf("max gradient error %.2e", max(err[s$free])))
  h <- check_hessian_accuracy(ms, dat, s$p, param_fixed = s$fixed, tol = 1e-3, n_quad = 8)
  expect_true(h$pass, info = sprintf("max Hessian error %.2e", h$max_error))
})

test_that("the order of the covariates does not change the typed likelihood or its derivatives", {
  skip_on_cran()
  dat <- .tic_data()
  msA <- .tic_system(dat, c("intercept", "x1"))
  msB <- .tic_system(dat, c("x1", "intercept"))
  sA <- .tic_point(msA, dat)
  sB <- .tic_point(msB, dat)
  expect_setequal(names(sA$p), names(sB$p))
  expect_equal(sB$p[names(sA$p)], sA$p)
  eA <- evaluate_likelihood_rcpp(msA, dat, sA$p, compute_gradient = TRUE, compute_hessian = TRUE, n_quad = 8)
  eB <- evaluate_likelihood_rcpp(msB, dat, sB$p, compute_gradient = TRUE, compute_hessian = TRUE, n_quad = 8)
  expect_equal(eB$loglik, eA$loglik, tolerance = 1e-12)
  gA <- setNames(eA$gradient, names(sA$p))
  gB <- setNames(eB$gradient, names(sB$p))
  expect_equal(gB[sA$free], gA[sA$free], tolerance = 1e-10)
  HA <- .tic_unpack(eA$hessian, names(sA$p))
  HB <- .tic_unpack(eB$hessian, names(sB$p))
  expect_equal(HB[sA$free, sA$free], HA[sA$free, sA$free], tolerance = 1e-10)
})
