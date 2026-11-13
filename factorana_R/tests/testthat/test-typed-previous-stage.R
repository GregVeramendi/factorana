# A typed stage used as previous_stage.
#
# Workflow being covered (the legacy HJV sequence): stage 1 estimates an untyped
# measurement system, stage A adds the type-probability model and typed outcome
# components (a probit and a multinomial logit), and stage B estimates further typed
# components with stage A held fixed. Stage B therefore fixes the type-probability
# parameters and stage A's type intercepts, and it fixes one of its own type intercepts
# at zero with fix_type_intercepts(), the normalization for a type that is rare in a cell.
#
# The other typed two-stage tests enter types in the last stage only
# (test-two-stage-se-types.R, validation/hjv_mini T9). Here the types are estimated in
# the middle stage and frozen in the next one. The reference is the joint model, which
# has stage A's and stage B's components in one system on top of stage 1: at the same
# parameter values the staged model must return the same log-likelihood, and the same
# gradient and Hessian for stage B's parameters.

.typed_prev_dgp <- function(n = 1500, seed = 7) {
  set.seed(seed)
  f <- rnorm(n)
  x <- rnorm(n)
  # Three types whose probabilities depend on the factor
  u <- cbind(0, 0.3 + 0.6 * f, -0.2 - 0.5 * f)
  p <- exp(u) / rowSums(exp(u))
  ty <- apply(p, 1, function(pp) sample(3, 1, prob = pp))
  t2 <- as.numeric(ty == 2)
  t3 <- as.numeric(ty == 3)
  d <- data.frame(
    intercept = 1, x = x,
    m1 = 0.5 + 1.0 * f + rnorm(n, 0, 0.6),
    m2 = -0.2 + 0.8 * f + rnorm(n, 0, 0.7),
    m3 = 0.1 + 1.2 * f + rnorm(n, 0, 0.5)
  )
  # Stage A outcomes: a probit and a three-alternative logit, both with type intercepts
  d$D <- as.integer(0.2 + 0.5 * x + 0.7 * f + 0.9 * t2 - 0.8 * t3 + rnorm(n) > 0)
  v <- cbind(0,
             0.1 + 0.4 * x + 0.5 * f + 1.0 * t2 - 0.6 * t3,
             -0.3 - 0.3 * x - 0.4 * f - 0.7 * t2 + 1.1 * t3)
  pc <- exp(v) / rowSums(exp(v))
  d$C <- apply(pc, 1, function(pp) sample(3, 1, prob = pp))
  # Stage B outcomes: a linear model and a probit. The probit's type-3 intercept is zero.
  d$Y <- 1.0 + 0.5 * x + 0.6 * f + 0.8 * t2 - 1.0 * t3 + rnorm(n, 0, 0.8)
  d$G <- as.integer(-0.1 + 0.3 * x + 0.5 * f + 0.7 * t2 + rnorm(n) > 0)
  d
}

.typed_prev_components <- function(dat, fm, which) {
  mk <- function(name, type, covariates = c("intercept", "x"), norm = NA_real_, typed = TRUE, ...) {
    define_model_component(name = name, data = dat, outcome = name, factor = fm,
                           covariates = covariates, model_type = type,
                           loading_normalization = norm, use_types = typed, ...)
  }
  switch(which,
    meas = list(m1 = mk("m1", "linear", "intercept", 1, FALSE),
                m2 = mk("m2", "linear", "intercept", NA_real_, FALSE),
                m3 = mk("m3", "linear", "intercept", NA_real_, FALSE)),
    A = list(D = mk("D", "probit"), C = mk("C", "logit", num_choices = 3)),
    B = list(Y = mk("Y", "linear"),
             G = fix_type_intercepts(mk("G", "probit"), types = 3))
  )
}

# Packed upper triangle (as evaluate_likelihood_rcpp() returns it) to a named matrix
.typed_prev_unpack <- function(v, nm) {
  n <- length(nm)
  H <- matrix(0, n, n, dimnames = list(nm, nm))
  k <- 1
  for (i in seq_len(n)) for (j in i:n) {
    H[i, j] <- v[k]
    H[j, i] <- v[k]
    k <- k + 1
  }
  H
}

# One fit of the three stages, shared by the tests below
.typed_prev_fits <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    dat <- .typed_prev_dgp()
    fm <- define_factor_model(n_factors = 1, n_types = 3)
    ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
    est <- function(ms, ...) estimate_model_rcpp(ms, dat, control = ctrl, optimizer = "nlminb",
                                                 parallel = FALSE, verbose = FALSE, ...)
    fit1 <- est(define_model_system(.typed_prev_components(dat, fm, "meas"), fm))
    fm1 <- fit1$model_system$factor
    fitA <- est(define_model_system(.typed_prev_components(dat, fm1, "A"), fm1, previous_stage = fit1))
    fmA <- fitA$model_system$factor
    msB <- define_model_system(.typed_prev_components(dat, fmA, "B"), fmA, previous_stage = fitA)
    fitB <- est(msB)
    msJ <- define_model_system(c(.typed_prev_components(dat, fm1, "A"), .typed_prev_components(dat, fm1, "B")),
                               fm1, previous_stage = fit1)
    cache <<- list(dat = dat, ctrl = ctrl, est = est, fit1 = fit1, fitA = fitA, fitB = fitB, msB = msB, msJ = msJ)
    cache
  }
})

test_that("a typed stage chains as previous_stage and keeps its type parameters fixed", {
  skip_on_cran()
  s <- .typed_prev_fits()
  expect_equal(s$fit1$convergence, 0)
  expect_equal(s$fitA$convergence, 0)
  expect_equal(s$fitB$convergence, 0)

  fixedB <- s$msB$previous_stage_info$fixed_param_names
  type_block <- c("typeprob_2_intercept", "typeprob_3_intercept", "type_2_loading_1", "type_3_loading_1")
  typed_A <- c("D_type_2_intercept", "D_type_3_intercept",
               "C_c1_type_2_intercept", "C_c1_type_3_intercept",
               "C_c2_type_2_intercept", "C_c2_type_3_intercept")
  expect_true(all(type_block %in% fixedB))
  expect_true(all(typed_A %in% fixedB))
  expect_setequal(fixedB, names(s$fitA$estimates))

  # Stage B leaves stage A's values untouched and estimates only its own parameters
  expect_identical(s$fitB$estimates[names(s$fitA$estimates)], s$fitA$estimates)
  pt <- s$fitB$param_table
  free_names <- pt$name[!pt$fixed]
  expect_setequal(free_names, c("Y_intercept", "Y_x", "Y_loading_1", "Y_sigma",
                                "Y_type_2_intercept", "Y_type_3_intercept",
                                "G_intercept", "G_x", "G_loading_1", "G_type_2_intercept"))
  # The normalized type intercept is in the vector, fixed at zero
  expect_true(pt$fixed[pt$name == "G_type_3_intercept"])
  expect_equal(unname(s$fitB$estimates["G_type_3_intercept"]), 0)
})

test_that("the staged typed model equals the joint model at the same parameters", {
  skip_on_cran()
  s <- .typed_prev_fits()
  invisible(capture.output({
    pJ <- initialize_parameters(s$msJ, s$dat, verbose = FALSE)$init_params
    pB <- initialize_parameters(s$msB, s$dat, verbose = FALSE)$init_params
  }))
  expect_identical(names(pJ), names(pB))

  # Stage 1 and stage A at their estimates, stage B moved off its optimum so that the
  # gradient is not zero
  b_names <- setdiff(names(s$fitB$estimates), names(s$fitA$estimates))
  b_free <- setdiff(b_names, "G_type_3_intercept")
  theta <- s$fitB$estimates
  theta[b_free] <- theta[b_free] + c(0.10, -0.05, 0.08, 0.03, 0.20, -0.15, -0.10, 0.05, 0.07, 0.12)
  pJ[names(theta)] <- theta
  pB[names(theta)] <- theta

  invisible(capture.output({
    eJ <- evaluate_likelihood_rcpp(s$msJ, s$dat, pJ, compute_gradient = TRUE, compute_hessian = TRUE, n_quad = 8)
    eB <- evaluate_likelihood_rcpp(s$msB, s$dat, pB, compute_gradient = TRUE, compute_hessian = TRUE, n_quad = 8)
  }))

  expect_equal(eB$loglik, eJ$loglik, tolerance = 1e-12)

  gJ <- setNames(eJ$gradient, names(pJ))
  gB <- setNames(eB$gradient, names(pB))
  expect_gt(max(abs(gJ[b_free])), 1)                       # the comparison is not at a stationary point
  expect_equal(gB[b_free], gJ[b_free], tolerance = 1e-10)
  expect_true(all(gB[setdiff(names(pB), b_free)] == 0))    # nothing outside stage B is free

  HJ <- .typed_prev_unpack(eJ$hessian, names(pJ))
  HB <- .typed_prev_unpack(eB$hessian, names(pB))
  expect_equal(HB[b_free, b_free], HJ[b_free, b_free], tolerance = 1e-10)
})

test_that("stage B after a typed stage: analytical gradient and Hessian match finite differences", {
  skip_on_cran()
  s <- .typed_prev_fits()
  params <- s$fitB$estimates
  b_free <- s$fitB$param_table$name[!s$fitB$param_table$fixed]
  params[b_free] <- params[b_free] + 0.05
  param_fixed <- !(names(params) %in% b_free)

  invisible(capture.output({
    grad_check <- check_gradient_accuracy(s$msB, s$dat, params, param_fixed = param_fixed,
                                          tol = 1e-4, verbose = FALSE, n_quad = 8)
    hess_check <- check_hessian_accuracy(s$msB, s$dat, params, param_fixed = param_fixed,
                                         tol = 1e-3, verbose = FALSE, n_quad = 8)
  }))
  expect_true(grad_check$pass, info = sprintf("gradient FD failed (max err %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian FD failed (max err %.2e)", hess_check$max_error))
})

test_that("the joint model started from the staged fits does not lose likelihood", {
  skip_on_cran()
  s <- .typed_prev_fits()
  invisible(capture.output(init <- initialize_parameters(s$msJ, s$dat, verbose = FALSE)$init_params))
  init[names(s$fitB$estimates)] <- s$fitB$estimates
  fitJ <- s$est(s$msJ, init_params = init)
  expect_equal(fitJ$convergence, 0)
  # The staged fit is a feasible point of the joint problem with the same likelihood
  # (previous test), so the joint optimum cannot lie below it
  expect_gte(fitJ$loglik, s$fitB$loglik - 1e-6)
  # The normalization carries into the joint model
  expect_equal(unname(fitJ$estimates["G_type_3_intercept"]), 0)
  expect_true(fitJ$param_table$fixed[fitJ$param_table$name == "G_type_3_intercept"])
})
