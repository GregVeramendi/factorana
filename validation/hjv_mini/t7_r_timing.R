# T7, R side: time one likelihood, one gradient and one Hessian evaluation at the truth (full system, 8 quadrature points,
# single core, as in legacy TestCode). Usage: Rscript t7_r_timing.R <N> <repeats>
args <- commandArgs(trailingOnly = TRUE); N <- as.integer(args[1]); reps <- if (length(args) >= 2) as.integer(args[2]) else 5L
suppressMessages(pkgload::load_all("../../factorana_R", compile = FALSE, quiet = TRUE)); source("model_r.R"); source("truth_params.R")
dat <- read.csv(sprintf("data/mini_%d.csv", N)); tr <- readRDS(sprintf("data/truth_%d.rds", N))$truth
fm <- mini_factor_model()
ms <- define_model_system(components = c(build_mini_system(dat, 1)$components, build_mini_system(dat, 3)$components, build_mini_system(dat, 4)$components), factor = fm)
init <- factorana:::initialize_parameters(ms, dat); p <- truth_r_params(tr, init$param_names)
meta <- factorana:::build_parameter_metadata(ms); cons <- factorana:::setup_parameter_constraints(ms, p, meta, init$factor_variance_fixed, verbose = FALSE)
fm_ptr <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8L, init_params = p); pf <- p[cons$free_idx]
# time a function, not an expression: a promise would be evaluated once and time zero on later repeats
tm <- function(f) median(vapply(seq_len(reps), function(i) system.time(f())[["elapsed"]], numeric(1)))
t_l <- tm(function() evaluate_loglik_only_cpp(fm_ptr, pf))
t_g <- tm(function() evaluate_likelihood_cpp(fm_ptr, pf, TRUE, FALSE))
t_h <- tm(function() evaluate_likelihood_cpp(fm_ptr, pf, TRUE, TRUE))
cat(sprintf("R N=%d (median of %d): loglik %.3f s | loglik+grad %.3f s | loglik+grad+hess %.3f s | free params %d\n", N, reps, t_l, t_g, t_h, length(pf)))
write.csv(data.frame(code = "R", N = N, loglik = t_l, grad = t_g, hess = t_h, reps = reps), sprintf("results/t7_r_eval_%d.csv", N), row.names = FALSE)
