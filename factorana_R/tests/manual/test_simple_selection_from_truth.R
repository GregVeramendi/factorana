# Test: Simple selection model, start from truth
library(factorana)

set.seed(777)
n <- 1000

f <- rnorm(n); Z <- rnorm(n); X <- rnorm(n)
D_latent <- 0.5 + 1.0*Z + 1.0*f + rnorm(n)
D <- as.integer(D_latent > 0)
Y1 <- 2.0 + 1.0*X + 2.0*f + rnorm(n, 0, 0.5)
T1 <- 1.0 + 1.0*f + rnorm(n, 0, 0.5)
T2 <- 0.5 + 1.2*f + rnorm(n, 0, 0.6)

dat <- data.frame(intercept=1, D=D, Z=Z, X=X, Y1=Y1, T1=T1, T2=T2,
                  eval_Y1=D, eval_all=1)

fm <- define_factor_model(n_factors=1, n_types=1, n_quad=16)
mc_D <- define_model_component(name="D", data=dat, outcome="D", factor=fm,
  covariates=c("intercept","Z"), model_type="probit",
  loading_normalization=NA_real_, evaluation_indicator="eval_all")
mc_Y1 <- define_model_component(name="Y1", data=dat, outcome="Y1", factor=fm,
  covariates=c("intercept","X"), model_type="linear",
  loading_normalization=NA_real_, evaluation_indicator="eval_Y1")
mc_T1 <- define_model_component(name="T1", data=dat, outcome="T1", factor=fm,
  covariates="intercept", model_type="linear",
  loading_normalization=NA_real_, evaluation_indicator="eval_all")
mc_T2 <- define_model_component(name="T2", data=dat, outcome="T2", factor=fm,
  covariates="intercept", model_type="linear",
  loading_normalization=NA_real_, evaluation_indicator="eval_all")

ms <- define_model_system(components=list(mc_D, mc_Y1, mc_T1, mc_T2), factor=fm)
fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad=16)

true_p <- c(1, 0.5,1,1, 2,1,2,0.5, 1,1,0.5, 0.5,1.2,0.6)
param_fixed <- c(TRUE, rep(FALSE,13))

ll_true <- evaluate_loglik_only_cpp(fm_cpp, true_p)
grad_true <- evaluate_likelihood_cpp(fm_cpp, true_p, TRUE, FALSE)$gradient[!param_fixed]

cat(sprintf("At TRUE parameters:\n"))
cat(sprintf("  LL: %.4f\n", ll_true))
cat(sprintf("  Grad L2 norm: %.2e\n", sqrt(sum(grad_true^2))))
cat(sprintf("  Max |grad|: %.2f\n\n", max(abs(grad_true))))

opt <- nlminb(true_p[!param_fixed],
  function(p) { full_p <- true_p; full_p[!param_fixed] <- p; -evaluate_loglik_only_cpp(fm_cpp, full_p) },
  gradient = function(p) { full_p <- true_p; full_p[!param_fixed] <- p;
    -evaluate_likelihood_cpp(fm_cpp, full_p, TRUE, FALSE)$gradient[!param_fixed] })

final_p <- true_p; final_p[!param_fixed] <- opt$par
ll_final <- -opt$objective

cat(sprintf("After optimization from truth:\n"))
cat(sprintf("  LL: %.4f (change: %+.4f)\n", ll_final, ll_final - ll_true))
cat(sprintf("  Converged: %d iterations\n", opt$iterations))

param_names <- c("fvar","D_int","D_Z","D_lam","Y1_int","Y1_X","Y1_lam","Y1_sig",
                 "T1_int","T1_lam","T1_sig","T2_int","T2_lam","T2_sig")
cat("\nParameter changes:\n")
for (i in which(!param_fixed)) {
  change <- 100*(final_p[i]-true_p[i])/true_p[i]
  if (abs(change) > 5) cat(sprintf("  %s: %+.1f%%\n", param_names[i], change))
}
