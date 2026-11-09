# cmp_eval.R <lib> <label> <mode>    mode: correctness | timing | peakmem
# Evaluates likelihood, gradient and Hessian on fixed, seeded model systems so that two
# builds of factorana can be compared element by element and timed on identical work.
args <- commandArgs(trailingOnly = TRUE)
lib <- normalizePath(args[1]); label <- args[2]; mode <- args[3]
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(factorana))
pkg_dir <- normalizePath(system.file(package = "factorana"))
cat(sprintf("factorana %s loaded from %s\n", packageVersion("factorana"), pkg_dir))
stopifnot(startsWith(pkg_dir, lib))
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])))
out_file <- file.path(script_dir, sprintf("res_%s_%s.rds", mode, label))

# ---------------------------------------------------------------- data ----
gen_data <- function(n, n_ind, seed = 1) {
  set.seed(seed)
  f <- matrix(rnorm(n * 3), n, 3)
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
  dat <- data.frame(id = seq_len(n), intercept = 1, x1 = x1, x2 = x2, eval = 1L)
  for (j in seq_len(n_ind)) {
    pat <- (j - 1) %% 3
    lam <- c(0, 0, 0)
    if (pat == 0) lam[1] <- 0.6 + 0.4 * runif(1)
    else if (pat == 1) lam[2] <- 0.6 + 0.4 * runif(1)
    else { lam[1] <- 0.3 + 0.3 * runif(1); lam[3] <- 0.6 + 0.4 * runif(1) }
    dat[[paste0("m", j)]] <- as.numeric(0.2 + 0.3 * x1 - 0.2 * x2 + f %*% lam + rnorm(n, 0, 0.6))
  }
  gum <- function() -log(-log(runif(n)))
  dat$pr <- as.integer(0.1 + 0.4 * x1 + 0.8 * f[, 1] + rnorm(n) > 0)
  u <- cbind(gum(), 0.2 + 0.3 * x1 + 0.7 * f[, 1] + 0.3 * f[, 3] + gum(),
             -0.1 - 0.2 * x2 + 0.4 * f[, 1] + 0.8 * f[, 3] + gum())
  dat$lg <- max.col(u)
  dat$y1 <- 0.5 + 0.4 * x1 + 0.5 * x2 + 0.6 * f[, 1] + 0.3 * f[, 2] - 0.4 * f[, 3] + rnorm(n, 0, 0.7)
  dat$y2 <- as.integer(-0.2 + 0.3 * x1 + 0.5 * f[, 1] + 0.4 * f[, 2] + rnorm(n) > 0)
  lat3 <- 0.3 * x1 + 0.5 * f[, 1] + 0.4 * f[, 3] + rnorm(n)
  dat$y3 <- as.integer(cut(lat3, breaks = c(-Inf, -0.5, 0.6, Inf), labels = FALSE))
  list(dat = dat, f = f)
}

# -------------------------------------------------------------- models ----
fm3 <- function(n_types = 1)
  fix_factor_param(define_factor_model(n_factors = 3, n_types = n_types),
                   c(factor_var_1 = 1, factor_var_2 = 1, factor_var_3 = 1))

mk_meas <- function(dat, fm, n_ind, extra = TRUE) {
  comps <- list()
  for (j in seq_len(n_ind)) {
    pat <- (j - 1) %% 3
    norm <- switch(pat + 1, c(NA_real_, 0, 0), c(0, NA_real_, 0), c(NA_real_, 0, NA_real_))
    comps[[paste0("m", j)]] <- define_model_component(paste0("m", j), dat, paste0("m", j), fm,
      covariates = c("intercept", "x1", "x2"), model_type = "linear",
      loading_normalization = norm, evaluation_indicator = "eval")
  }
  if (extra) {
    comps$pr <- define_model_component("pr", dat, "pr", fm, covariates = c("intercept", "x1"),
      model_type = "probit", loading_normalization = c(NA_real_, 0, 0), evaluation_indicator = "eval")
    lg <- define_model_component("lg", dat, "lg", fm, covariates = c("intercept", "x1", "x2"),
      model_type = "logit", num_choices = 3, loading_normalization = c(NA_real_, 0, NA_real_),
      evaluation_indicator = "eval")
    comps$lg <- fix_coefficient(lg, "x2", 0, choice = 1)
  }
  comps
}

mk_out <- function(dat, fm, use_types = FALSE) {
  y1 <- define_model_component("y1", dat, "y1", fm, covariates = c("intercept", "x1", "x2"),
    model_type = "linear", loading_normalization = c(NA_real_, NA_real_, NA_real_),
    evaluation_indicator = "eval", use_types = use_types)
  y1 <- fix_coefficient(y1, "x2", 0.5)
  y2 <- define_model_component("y2", dat, "y2", fm, covariates = c("intercept", "x1"),
    model_type = "probit", loading_normalization = c(NA_real_, NA_real_, 0),
    evaluation_indicator = "eval", use_types = use_types)
  y3 <- define_model_component("y3", dat, "y3", fm, covariates = "x1", model_type = "oprobit",
    num_choices = 3, loading_normalization = c(NA_real_, 0, NA_real_),
    evaluation_indicator = "eval", use_types = use_types)
  list(y1 = y1, y2 = y2, y3 = y3)
}

fake_fit <- function(ms1, dat) {
  est <- suppressWarnings(initialize_parameters(ms1, dat, verbose = FALSE)$init_params)
  structure(list(model_system = ms1, estimates = est,
                 std_errors = setNames(rep(0.1, length(est)), names(est))),
            class = "factorana_result")
}

adaptive_inputs <- function(f, seed = 7) {
  set.seed(seed)
  n <- nrow(f)
  list(scores = f + matrix(rnorm(n * 3, 0, 0.3), n, 3),
       ses = matrix(runif(n * 3, 0.15, 1.3), n, 3), vars = c(1, 1, 1), threshold = 0.5)
}

make_ptr <- function(ms, dat, n_quad, adaptive = NULL) {
  init <- suppressWarnings(initialize_parameters(ms, dat, verbose = FALSE)$init_params)
  fm_ptr <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad, init)
  if (!is.null(adaptive))
    set_adaptive_quadrature_cpp(fm_ptr, adaptive$scores, adaptive$ses, adaptive$vars,
                                threshold = adaptive$threshold, max_quad = n_quad, verbose = FALSE)
  free <- extract_free_params_cpp(fm_ptr, init)
  list(ptr = fm_ptr, init = init, free = free)
}

eval_all <- function(p) {
  l1 <- evaluate_likelihood_cpp(p$ptr, p$free, FALSE, FALSE)
  g2 <- evaluate_likelihood_cpp(p$ptr, p$free, TRUE, FALSE)
  h3 <- evaluate_likelihood_cpp(p$ptr, p$free, TRUE, TRUE)
  list(nparam = length(p$init), nfree = length(p$free),
       loglik_1 = l1$logLikelihood, loglik_2 = g2$logLikelihood, loglik_3 = h3$logLikelihood,
       gradient_2 = g2$gradient, gradient_3 = h3$gradient, hessian = h3$hessian)
}

scenario <- function(name, expr) {
  cat(sprintf("  %-28s ", name))
  res <- tryCatch(expr, error = function(e) list(error = conditionMessage(e)))
  if (!is.null(res$error)) cat("ERROR:", res$error, "\n")
  else cat(sprintf("nparam %4d  nfree %4d  loglik %.6f\n", res$nparam, res$nfree, res$loglik_3))
  res
}

# ---------------------------------------------------------- correctness ----
if (mode == "correctness") {
  n <- 300; n_quad <- 4
  g <- gen_data(n, n_ind = 9, seed = 11); dat <- g$dat
  res <- list()

  fm <- fm3()
  meas <- mk_meas(dat, fm, 9)
  ms1 <- define_model_system(components = meas, factor = fm)
  fit1 <- fake_fit(ms1, dat)
  outs <- mk_out(dat, fm)

  res$fixed_heavy <- scenario("fixed_heavy", {
    ms2 <- define_model_system(components = outs, factor = fm, previous_stage = fit1)
    eval_all(make_ptr(ms2, dat, n_quad))
  })
  res$all_free <- scenario("all_free", {
    ms <- define_model_system(components = c(meas, outs), factor = fm)
    eval_all(make_ptr(ms, dat, n_quad))
  })
  res$equality_all_free <- scenario("equality_all_free", {
    ms <- define_model_system(components = c(meas, outs), factor = fm,
      equality_constraints = list(c("m1_sigma", "m4_sigma"), c("m2_loading_2", "m5_loading_2")))
    eval_all(make_ptr(ms, dat, n_quad))
  })
  res$adaptive_fixed_heavy <- scenario("adaptive_fixed_heavy", {
    ms2 <- define_model_system(components = outs, factor = fm, previous_stage = fit1)
    eval_all(make_ptr(ms2, dat, n_quad, adaptive = adaptive_inputs(g$f)))
  })
  res$adaptive_all_free <- scenario("adaptive_all_free", {
    ms <- define_model_system(components = c(meas, outs), factor = fm)
    eval_all(make_ptr(ms, dat, n_quad, adaptive = adaptive_inputs(g$f)))
  })

  fmt <- fm3(n_types = 2)
  meas_t <- mk_meas(dat, fmt, 9)
  outs_t <- mk_out(dat, fmt, use_types = TRUE)
  res$typed_all_free <- scenario("typed_all_free", {
    ms <- define_model_system(components = c(meas_t, outs_t), factor = fmt)
    eval_all(make_ptr(ms, dat, n_quad))
  })
  res$typed_fixed_heavy <- scenario("typed_fixed_heavy", {
    ms1t <- define_model_system(components = meas_t, factor = fmt)
    ms2t <- define_model_system(components = outs_t, factor = fmt, previous_stage = fake_fit(ms1t, dat))
    eval_all(make_ptr(ms2t, dat, n_quad))
  })
  res$adaptive_typed_fixed_heavy <- scenario("adaptive_typed_fixed_heavy", {
    ms1t <- define_model_system(components = meas_t, factor = fmt)
    ms2t <- define_model_system(components = outs_t, factor = fmt, previous_stage = fake_fit(ms1t, dat))
    eval_all(make_ptr(ms2t, dat, n_quad, adaptive = adaptive_inputs(g$f)))
  })

  res$mixture <- scenario("mixture_2fac_2mix", {
    fmm <- define_factor_model(n_factors = 2, n_mixtures = 2)
    mk <- function(nm, norm) define_model_component(nm, dat, nm, fmm, covariates = c("intercept", "x1"),
      model_type = "linear", loading_normalization = norm, evaluation_indicator = "eval")
    m4 <- fix_coefficient(mk("m4", c(NA_real_, 0)), "x1", 0.3)
    ms <- define_model_system(components = list(mk("m1", c(1, 0)), m4, mk("m2", c(0, 1)),
                                                mk("m5", c(0, NA_real_))), factor = fmm)
    eval_all(make_ptr(ms, dat, n_quad))
  })
  res$mixture_probit <- scenario("mixture_1fac_2mix_probit", {
    fmm <- define_factor_model(n_factors = 1, n_mixtures = 2)
    mk <- function(nm, type, norm) define_model_component(nm, dat, nm, fmm, covariates = c("intercept", "x1"),
      model_type = type, loading_normalization = norm, evaluation_indicator = "eval")
    ms <- define_model_system(components = list(mk("m1", "linear", 1), mk("m4", "linear", NA_real_),
                                                mk("pr", "probit", NA_real_)), factor = fmm)
    eval_all(make_ptr(ms, dat, n_quad))
  })
  saveRDS(res, out_file)
  cat("saved", out_file, "\n")
}

# --------------------------------------------------------------- timing ----
time_one <- function(p, grad, hess, reps) {
  t <- numeric(reps)
  for (r in seq_len(reps)) t[r] <- system.time(evaluate_likelihood_cpp(p$ptr, p$free, grad, hess))[["elapsed"]]
  t
}
timing_block <- function(name, p, reps_lgh = c(3, 3, 2)) {
  tl <- time_one(p, FALSE, FALSE, reps_lgh[1]); tg <- time_one(p, TRUE, FALSE, reps_lgh[2]); th <- time_one(p, TRUE, TRUE, reps_lgh[3])
  out <- list(name = name, nparam = length(p$init), nfree = length(p$free),
              lik = tl, grad = tg, hess = th)
  cat(sprintf("  %-22s nparam %4d nfree %4d | L %7.3f s | G %7.3f s | H %8.3f s | H/L %7.1f\n", name, out$nparam, out$nfree,
              median(tl), median(tg), median(th), median(th) / median(tl)))
  out
}
build_T <- function(n, n_ind, n_quad, fixed_heavy) {
  g <- gen_data(n, n_ind, seed = 21); dat <- g$dat
  fm <- fm3(); meas <- mk_meas(dat, fm, n_ind); outs <- mk_out(dat, fm)
  if (fixed_heavy) {
    ms1 <- define_model_system(components = meas, factor = fm)
    ms <- define_model_system(components = outs, factor = fm, previous_stage = fake_fit(ms1, dat))
  } else ms <- define_model_system(components = c(meas, outs), factor = fm)
  make_ptr(ms, dat, n_quad)
}
if (mode == "timing") {
  res <- list()
  res$T1 <- timing_block("T1 fixed_heavy N5000", build_T(5000, 120, 4, TRUE))
  res$T3 <- timing_block("T3 fixed_heavy N2000", build_T(2000, 40, 4, TRUE))
  res$T2 <- timing_block("T2 all_free N2000", build_T(2000, 40, 4, FALSE))
  saveRDS(res, out_file)
  cat("saved", out_file, "\n")
}
if (mode == "peakmem") {
  p <- build_T(5000, 120, 4, TRUE)
  r <- evaluate_likelihood_cpp(p$ptr, p$free, TRUE, TRUE)
  cat(sprintf("peakmem run: nparam %d nfree %d loglik %.6f\n", length(p$init), length(p$free), r$logLikelihood))
}

# ------------------------------------------------- big fixed-heavy model ----
# 60 linear indicators with n_xtra extra covariates each, so the total parameter count reaches
# the HJV Stage 4 scale (about 3,550) without the likelihood underflow that hundreds of
# indicators would cause. Only the outcome components are free.
gen_data_x <- function(n, n_ind, n_xtra, seed = 31) {
  g <- gen_data(n, n_ind, seed = seed); dat <- g$dat
  set.seed(seed + 1000)
  for (k in seq_len(n_xtra)) {
    z <- rnorm(n); dat[[paste0("z", k)]] <- z
    for (j in seq_len(n_ind)) dat[[paste0("m", j)]] <- dat[[paste0("m", j)]] + 0.02 * rnorm(1) * z
  }
  list(dat = dat, f = g$f)
}
mk_meas_x <- function(dat, fm, n_ind, n_xtra) {
  covs <- c("intercept", "x1", "x2", paste0("z", seq_len(n_xtra)))
  comps <- list()
  for (j in seq_len(n_ind)) {
    pat <- (j - 1) %% 3
    norm <- switch(pat + 1, c(NA_real_, 0, 0), c(0, NA_real_, 0), c(NA_real_, 0, NA_real_))
    comps[[paste0("m", j)]] <- define_model_component(paste0("m", j), dat, paste0("m", j), fm,
      covariates = covs, model_type = "linear", loading_normalization = norm, evaluation_indicator = "eval")
  }
  comps
}
build_T4 <- function(n, n_ind, n_xtra, n_quad) {
  g <- gen_data_x(n, n_ind, n_xtra); dat <- g$dat
  fm <- fm3(); meas <- mk_meas_x(dat, fm, n_ind, n_xtra); outs <- mk_out(dat, fm)
  ms1 <- define_model_system(components = meas, factor = fm)
  ms <- define_model_system(components = outs, factor = fm, previous_stage = fake_fit(ms1, dat))
  make_ptr(ms, dat, n_quad)
}
if (mode == "timing_big") {
  t0 <- proc.time()[["elapsed"]]
  p <- build_T4(1000, 60, 53, 4)
  cat(sprintf("  build time %.1f s\n", proc.time()[["elapsed"]] - t0))
  res <- list(T4 = timing_block("T4 fixed_heavy N1000", p, reps_lgh = c(3, 3, 2)))
  cat(sprintf("  loglik %.6f\n", evaluate_likelihood_cpp(p$ptr, p$free, FALSE, FALSE)$logLikelihood))
  saveRDS(res, out_file)
  cat("saved", out_file, "\n")
}
