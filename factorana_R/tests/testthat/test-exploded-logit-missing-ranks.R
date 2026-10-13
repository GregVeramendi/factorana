# A rank that holds no valid choice is skipped (src/Model.cpp, the rank loop).
# Coding an unused rank as NA, -9999, 0 or an out-of-range number must give the
# same likelihood and gradient.

test_that("unused ranks give the same likelihood whatever their code", {
  set.seed(42)
  n <- 200
  K <- 4
  dat <- data.frame(intercept = 1, x1 = rnorm(n), f = rnorm(n))
  for (r in 1:3) dat[[paste0("rank", r)]] <- sample(1:K, n, replace = TRUE)
  dat$T1 <- dat$f + rnorm(n, 0, 0.5)
  dat$T2 <- 0.8 * dat$f + rnorm(n, 0, 0.6)
  dat$eval <- 1
  u <- runif(n)
  drop1 <- u < 0.2                # first rank unused, later ranks used
  drop3 <- u >= 0.2 & u < 0.4     # last rank unused
  dropall <- u >= 0.4 & u < 0.45  # no rank used

  build <- function(code) {
    d <- dat
    d$rank1[drop1] <- code
    d$rank3[drop3] <- code
    d[dropall, c("rank1", "rank2", "rank3")] <- code
    fm <- define_factor_model(n_factors = 1)
    comps <- list(
      define_model_component(name = "T1", data = d, outcome = "T1", factor = fm,
        covariates = NULL, model_type = "linear", loading_normalization = 1,
        evaluation_indicator = "eval", intercept = FALSE),
      define_model_component(name = "T2", data = d, outcome = "T2", factor = fm,
        covariates = NULL, model_type = "linear", loading_normalization = NA_real_,
        evaluation_indicator = "eval", intercept = FALSE),
      define_model_component(name = "choice", data = d,
        outcome = c("rank1", "rank2", "rank3"), factor = fm,
        covariates = c("intercept", "x1"), model_type = "logit", num_choices = K,
        exclude_chosen = FALSE, loading_normalization = NA_real_,
        evaluation_indicator = "eval"))
    list(ms = define_model_system(components = comps, factor = fm), d = d)
  }

  b0 <- build(0)
  params <- initialize_parameters(b0$ms, b0$d, verbose = FALSE)$init_params
  eval_at <- function(b) {
    ptr <- initialize_factor_model_cpp(b$ms, b$d, 8)
    evaluate_likelihood_cpp(ptr, params, compute_gradient = TRUE)
  }
  ref <- eval_at(b0)
  expect_true(is.finite(ref$logLikelihood))
  for (code in list(NA_real_, -9999, 1e10)) {
    res <- eval_at(build(code))
    expect_identical(res$logLikelihood, ref$logLikelihood)
    expect_identical(res$gradient, ref$gradient)
  }
})
