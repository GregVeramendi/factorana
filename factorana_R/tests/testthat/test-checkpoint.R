# The iteration log of nlminb, the checkpoint file, and a restart from it.
#
# estimate_model_rcpp() prints one line per iteration of nlminb from its gradient function, which
# nlminb calls at its starting point and at each point it moves to. With checkpoint_file it saves
# the best parameters so far each time the Hessian is evaluated, and read_checkpoint() returns
# them as starting values.
#
# Design: one factor measured by three linear indicators (stage 1), and a linear and a probit
# outcome that load on the factor (stage 2, with the stage 1 parameters held fixed).

.ckpt_data <- function(n = 600, seed = 11) {
  set.seed(seed)
  f <- rnorm(n); x <- rnorm(n)
  data.frame(intercept = 1, x = x,
             y1 = 0.2 + 1.0 * f + rnorm(n, 0, 0.6),
             y2 = -0.1 + 0.8 * f + rnorm(n, 0, 0.7),
             y3 = 0.3 + 0.6 * f + rnorm(n, 0, 0.5),
             w = 0.5 + 0.4 * x + 0.7 * f + rnorm(n, 0, 0.8),
             d = as.integer(-0.2 + 0.5 * x + 0.9 * f + rnorm(n) > 0))
}

.ckpt_stage1 <- function(dat) {
  fm <- define_factor_model(n_factors = 1)
  ind <- function(nm, norm) define_model_component(nm, dat, nm, fm, covariates = "intercept",
                                                   model_type = "linear", loading_normalization = norm)
  define_model_system(components = list(ind("y1", 1), ind("y2", NA_real_), ind("y3", NA_real_)), factor = fm)
}

.ckpt_stage2 <- function(dat, stage1_fit) {
  fm <- stage1_fit$model_system$factor
  out <- function(nm, type) define_model_component(nm, dat, nm, fm, covariates = c("intercept", "x"),
                                                   model_type = type, loading_normalization = NA_real_)
  define_model_system(components = list(out("w", "linear"), out("d", "probit")), factor = fm,
                      previous_stage = stage1_fit)
}

# Run an estimation and keep what it prints: the messages (the log) and standard output
.ckpt_run <- function(...) {
  log <- character(0)
  out <- capture.output(
    fit <- withCallingHandlers(estimate_model_rcpp(...),
                               message = function(m) {
                                 log <<- c(log, conditionMessage(m))
                                 invokeRestart("muffleMessage")
                               }))
  list(fit = fit, log = sub("\n$", "", log), out = out)
}

.ckpt_iterations <- function(log) {
  lines <- grep("^  Iteration [0-9]+: loglik = ", log, value = TRUE)
  data.frame(iteration = as.integer(sub("^  Iteration ([0-9]+):.*$", "\\1", lines)),
             loglik = as.numeric(sub("^.*loglik = (-?[0-9.]+),.*$", "\\1", lines)))
}

test_that("nlminb logs one line per iteration with the log-likelihood", {
  skip_on_cran()
  dat <- .ckpt_data()
  ms <- .ckpt_stage1(dat)
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  run <- .ckpt_run(ms, dat, control = ctrl, optimizer = "nlminb", parallel = FALSE, verbose = TRUE)
  it <- .ckpt_iterations(run$log)
  expect_equal(run$fit$convergence, 0)
  # one line per gradient evaluation, numbered as nlminb numbers its iterations
  expect_identical(it$iteration, seq_len(nrow(it)) - 1L)
  expect_equal(nrow(it), unname(run$fit$evaluations[["gradient"]]))
  expect_equal(max(it$iteration), run$fit$iterations)
  # the log-likelihood rises from line to line and ends at the estimate
  expect_true(all(diff(it$loglik) > -1e-6))
  expect_equal(it$loglik[nrow(it)], run$fit$loglik, tolerance = 1e-6)
  # nlminb's own trace, which prints every parameter, is off
  expect_false(any(grepl("^ *[0-9]+: +[0-9.]+:", run$out)))
  # no iteration lines without verbose
  quiet <- .ckpt_run(ms, dat, control = ctrl, optimizer = "nlminb", parallel = FALSE, verbose = FALSE)
  expect_length(quiet$log, 0)
  expect_length(quiet$out, 0)
})

test_that("the parallel gradient reports the same log-likelihood", {
  skip_on_cran()
  dat <- .ckpt_data()
  ms <- .ckpt_stage1(dat)
  serial <- .ckpt_run(ms, dat, control = define_estimation_control(n_quad_points = 8, num_cores = 1),
                      optimizer = "nlminb", parallel = FALSE, verbose = TRUE)
  par <- .ckpt_run(ms, dat, control = define_estimation_control(n_quad_points = 8, num_cores = 2),
                   optimizer = "nlminb", parallel = TRUE, verbose = TRUE)
  expect_equal(par$fit$loglik, serial$fit$loglik, tolerance = 1e-8)
  expect_equal(.ckpt_iterations(par$log), .ckpt_iterations(serial$log), tolerance = 1e-6)
})

test_that("a checkpoint holds the best parameters and restarts the estimation", {
  skip_on_cran()
  dat <- .ckpt_data()
  ms <- .ckpt_stage1(dat)
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  ckpt <- tempfile(fileext = ".csv")
  on.exit(unlink(c(ckpt, paste0(ckpt, ".tmp"))), add = TRUE)
  run <- .ckpt_run(ms, dat, control = ctrl, optimizer = "nlminb", parallel = FALSE, verbose = TRUE,
                   checkpoint_file = ckpt)
  fit <- run$fit
  expect_true(file.exists(ckpt))
  expect_false(file.exists(paste0(ckpt, ".tmp")))
  # one checkpoint per iteration: nlminb evaluates the Hessian at each point it moves to
  expect_equal(sum(grepl("Checkpoint [0-9]+ saved", run$log)), nrow(.ckpt_iterations(run$log)))

  saved <- read_checkpoint(ckpt, ms)
  expect_named(saved, c("params", "loglik", "saved"))
  expect_identical(names(saved$params), fit$param_names)
  expect_length(saved$params, length(fit$estimates))
  expect_equal(unname(saved$params), unname(fit$estimates), tolerance = 1e-4)
  expect_equal(saved$loglik, fit$loglik, tolerance = 1e-6)
  expect_false(is.na(saved$saved))

  # the log-likelihood in the file is that of its parameters
  free <- saved$params[fit$free_idx]
  fm_ptr <- initialize_factor_model_cpp(ms, as.matrix(dat), 8, saved$params)
  expect_equal(evaluate_loglik_only_cpp(fm_ptr, free), saved$loglik, tolerance = 1e-8)

  # a run that starts from the checkpoint is at the optimum: same estimates, few iterations
  refit <- estimate_model_rcpp(ms, dat, init_params = saved$params, control = ctrl, optimizer = "nlminb",
                               parallel = FALSE, verbose = FALSE)
  expect_equal(refit$convergence, 0)
  expect_equal(refit$loglik, fit$loglik, tolerance = 1e-8)
  expect_equal(unname(refit$estimates), unname(fit$estimates), tolerance = 1e-4)
  expect_lte(refit$iterations, 2)
  expect_lt(refit$iterations, fit$iterations)
})

test_that("a checkpoint written partway restarts the run at that point", {
  skip_on_cran()
  dat <- .ckpt_data()
  ms <- .ckpt_stage1(dat)
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  ckpt <- tempfile(fileext = ".csv")
  kept <- tempfile(fileext = ".csv")
  on.exit(unlink(c(ckpt, kept)), add = TRUE)
  # keep the file as it stands after the second iteration, as a stopped run would leave it
  full <- withCallingHandlers(
    estimate_model_rcpp(ms, dat, control = ctrl, optimizer = "nlminb", parallel = FALSE, verbose = TRUE,
                        checkpoint_file = ckpt),
    message = function(m) {
      if (grepl("Checkpoint 3 saved", conditionMessage(m))) file.copy(ckpt, kept, overwrite = TRUE)
      invokeRestart("muffleMessage")
    })
  expect_gt(full$iterations, 3)
  partway <- read_checkpoint(kept, ms)
  expect_lt(partway$loglik, full$loglik - 1e-3)
  resumed <- .ckpt_run(ms, dat, init_params = partway$params, control = ctrl, optimizer = "nlminb",
                       parallel = FALSE, verbose = TRUE)
  # the resumed run starts at the log-likelihood of the checkpoint and reaches the same optimum
  expect_equal(.ckpt_iterations(resumed$log)$loglik[1], partway$loglik, tolerance = 1e-6)
  expect_equal(resumed$fit$convergence, 0)
  expect_equal(resumed$fit$loglik, full$loglik, tolerance = 1e-6)
  expect_equal(unname(resumed$fit$estimates), unname(full$estimates), tolerance = 1e-3)
})

test_that("read_checkpoint stops on a file that does not fit the model", {
  skip_on_cran()
  dat <- .ckpt_data()
  ms1 <- .ckpt_stage1(dat)
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  fit1 <- estimate_model_rcpp(ms1, dat, control = ctrl, optimizer = "nlminb", parallel = FALSE, verbose = FALSE)
  ms2 <- .ckpt_stage2(dat, fit1)
  ckpt1 <- tempfile(fileext = ".csv")
  ckpt2 <- tempfile(fileext = ".csv")
  on.exit(unlink(c(ckpt1, ckpt2)), add = TRUE)
  estimate_model_rcpp(ms1, dat, control = ctrl, optimizer = "nlminb", parallel = FALSE, verbose = FALSE,
                      checkpoint_file = ckpt1)
  fit2 <- estimate_model_rcpp(ms2, dat, control = ctrl, optimizer = "nlminb", parallel = FALSE, verbose = FALSE,
                              checkpoint_file = ckpt2)

  expect_error(read_checkpoint(tempfile()), "not found")
  # the stage 1 file has other parameters than the stage 2 model
  expect_error(read_checkpoint(ckpt1, ms2), "does not fit the model")
  # without a model the file is read as it is
  expect_length(read_checkpoint(ckpt1)$params, length(fit1$estimates))

  # the stage 2 file fits the stage 2 model, and restarts it with the stage 1 values fixed
  saved <- read_checkpoint(ckpt2, ms2)
  fixed <- ms2$previous_stage_info$fixed_param_names
  expect_gt(length(fixed), 0)
  expect_equal(unname(saved$params[fixed]), unname(fit1$estimates[fixed]), tolerance = 1e-12)
  refit <- estimate_model_rcpp(ms2, dat, init_params = saved$params, control = ctrl, optimizer = "nlminb",
                               parallel = FALSE, verbose = FALSE)
  expect_equal(refit$loglik, fit2$loglik, tolerance = 1e-8)
  expect_equal(unname(refit$estimates[fixed]), unname(fit1$estimates[fixed]), tolerance = 1e-12)

  # against other stage 1 results the file is refused
  other1 <- fit1
  other1$estimates["y2_sigma"] <- other1$estimates["y2_sigma"] * 1.01
  expect_error(read_checkpoint(ckpt2, .ckpt_stage2(dat, other1)), "other previous-stage results")

  # a truncated file is refused
  writeLines(readLines(ckpt2)[1:3], ckpt2)
  expect_error(read_checkpoint(ckpt2, ms2), "Not a complete checkpoint")
})
