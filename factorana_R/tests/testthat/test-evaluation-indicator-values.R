# define_model_component() requires every evaluation indicator value to be 0 or 1.

ev_component <- function(dat, fm) {
  define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = NULL, model_type = "linear",
    loading_normalization = 1, evaluation_indicator = "ev",
    intercept = FALSE
  )
}

test_that("0/1 numeric and FALSE/TRUE indicators are accepted", {
  set.seed(7)
  n <- 30
  fm <- define_factor_model(n_factors = 1)
  dat <- data.frame(T1 = rnorm(n), ev = rep(c(0, 1, 1), length.out = n))
  expect_s3_class(ev_component(dat, fm), "model_component")
  dat$ev <- dat$ev == 1
  expect_s3_class(ev_component(dat, fm), "model_component")
})

test_that("an indicator value other than 0 or 1 stops the definition", {
  set.seed(7)
  n <- 30
  fm <- define_factor_model(n_factors = 1)
  dat <- data.frame(T1 = rnorm(n), ev = 1)
  dat$ev[3] <- 2
  expect_error(ev_component(dat, fm), "must be 0 or 1 in every row. 1 row\\(s\\) hold another value \\(2\\)")
  dat$ev[3] <- NA
  expect_error(ev_component(dat, fm), "must be 0 or 1 in every row")
  dat$ev <- TRUE
  dat$ev[3] <- NA
  expect_error(ev_component(dat, fm), "must be 0 or 1 in every row")
})
