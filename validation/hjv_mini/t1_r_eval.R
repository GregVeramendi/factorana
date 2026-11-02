# T1, R side: log-likelihood, gradient and Hessian at the true parameters, full system, full quadrature.
# Usage: Rscript t1_r_eval.R <N> [full|meas] [ntypes]. Writes results/t1_r<tag>_<N>.rds and ..._params.csv (R name order).
# tag: "" full untyped system, "m" measurement only, "t<T>" full system with T latent types (T9; data suffix _t<T>).
args <- commandArgs(trailingOnly = TRUE); N <- as.integer(args[1]); sys <- if (length(args) >= 2) args[2] else "full"
ntypes <- if (length(args) >= 3) as.integer(args[3]) else 1L
tag <- if (sys != "full") "m" else if (ntypes > 1) sprintf("t%d", ntypes) else ""
dtag <- if (ntypes > 1) sprintf("_t%d", ntypes) else ""
suppressMessages(pkgload::load_all("../../factorana_R", compile = FALSE, quiet = TRUE))
source("model_r.R"); source("truth_params.R")
dat <- read.csv(sprintf("data/mini_%d%s.csv", N, dtag)); tr <- readRDS(sprintf("data/truth_%d%s.rds", N, dtag))$truth
fm <- mini_factor_model(ntypes)
ms <- if (sys == "full") define_model_system(components = c(build_mini_system(dat, 1, ntypes = ntypes)$components,
                                                            build_mini_system(dat, 5, ntypes = ntypes)$components), factor = fm) else build_mini_system(dat, 1)
init <- factorana:::initialize_parameters(ms, dat)
p <- truth_r_params(tr, init$param_names, subset = TRUE)
# evaluate_likelihood_cpp() takes the FREE parameters only; fixed values (factor variances, fixed coefficient)
# come from init_params given to initialize_factor_model_cpp(). free_idx follows the package's own constraint setup.
meta <- factorana:::build_parameter_metadata(ms)
cons <- factorana:::setup_parameter_constraints(ms, p, meta, init$factor_variance_fixed, verbose = FALSE)
free_idx <- cons$free_idx
cat("free parameters (R):", length(free_idx), " fixed:", length(p) - length(free_idx), "\n")
fm_ptr <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8L, init_params = p)
t <- system.time(res <- evaluate_likelihood_cpp(fm_ptr, p[free_idx], compute_gradient = TRUE, compute_hessian = TRUE))
cat(sprintf("R T1 N=%d: loglik = %.10f  (%.1f s for lkhd+grad+hess)\n", N, res$logLikelihood, t[["elapsed"]]))
cat("gradient length:", length(res$gradient), " hessian length:", length(res$hessian), "\n")
saveRDS(list(params = p, loglik = res$logLikelihood, gradient = res$gradient, hessian = res$hessian,
             param_names = init$param_names, free_names = init$param_names[free_idx], secs = t[["elapsed"]]), sprintf("results/t1_r%s_%d.rds", tag, N))
write.csv(data.frame(idx = seq_along(p) - 1, param = names(p), value = unname(p)), sprintf("results/t1_r%s_%d_params.csv", tag, N), row.names = FALSE)
