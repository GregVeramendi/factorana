# Test: Three tests (no missing data), start from truth
library(factorana)

set.seed(888)
n <- 1000

f <- rnorm(n)
T1 <- 1.0 + 1.0*f + rnorm(n, 0, 0.5)
T2 <- 0.5 + 1.2*f + rnorm(n, 0, 0.6)
T3 <- -0.3 + 0.8*f + rnorm(n, 0, 0.4)

dat <- data.frame(intercept=1, T1=T1, T2=T2, T3=T3, eval=1)

fm <- define_factor_model(n_factors=1, n_types=1, n_quad=16)
mc_T1 <- define_model_component(name="T1", data=dat, outcome="T1", factor=fm,
  covariates="intercept", model_type="linear",
  loading_normalization=NA_real_, evaluation_indicator="eval")
mc_T2 <- define_model_component(name="T2", data=dat, outcome="T2", factor=fm,
  covariates="intercept", model_type="linear",
  loading_normalization=NA_real_, evaluation_indicator="eval")
mc_T3 <- define_model_component(name="T3", data=dat, outcome="T3", factor=fm,
  covariates="intercept", model_type="linear",
  loading_normalization=NA_real_, evaluation_indicator="eval")

ms <- define_model_system(components=list(mc_T1,mc_T2,mc_T3), factor=fm)
fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad=16)

true_p <- c(1, 1,1,0.5, 0.5,1.2,0.6, -0.3,0.8,0.4)
param_fixed <- c(TRUE, rep(FALSE,9))

ll_true <- evaluate_loglik_only_cpp(fm_cpp, true_p)
grad_true <- evaluate_likelihood_cpp(fm_cpp, true_p, TRUE, FALSE)$gradient[!param_fixed]

cat(sprintf("At TRUE parameters (no missing data):\n"))
cat(sprintf("  LL: %.4f\n", ll_true))
cat(sprintf("  Grad L2 norm: %.2e\n", sqrt(sum(grad_true^2))))
cat(sprintf("  Max |grad|: %.2f\n\n", max(abs(grad_true))))

opt <- nlminb(true_p[!param_fixed],
  function(p) { full_p <- true_p; full_p[!param_fixed] <- p; -evaluate_loglik_only_cpp(fm_cpp, full_p) },
  gradient = function(p) { full_p <- true_p; full_p[!param_fixed] <- p;
    -evaluate_likelihood_cpp(fm_cpp, full_p, TRUE, FALSE)$gradient[!param_fixed] })

ll_final <- -opt$objective

cat(sprintf("After optimization from truth:\n"))
cat(sprintf("  LL: %.4f (change: %+.4f)\n", ll_final, ll_final - ll_true))
cat(sprintf("  Converged: %d iterations\n", opt$iterations))

if (abs(ll_final - ll_true) < 0.01) {
  cat("\n✓ Optimizer stayed at truth (no missing data case works!)\n")
} else {
  cat("\n⚠ Optimizer moved away (problem exists even without missing data)\n")
}
