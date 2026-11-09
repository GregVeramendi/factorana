# T7c. Usage: Rscript t7c_version_timing.R <lib> <label>   (single core, N = 20000 mini stage 1, same parameter point for every version)
args <- commandArgs(TRUE); lib <- args[1]; label <- args[2]; reps <- 3L
.libPaths(c(lib, .libPaths())); suppressMessages(library(factorana)); ver <- as.character(packageVersion("factorana"))
# run from validation/hjv_mini; the library given as <lib> is searched first
if (!"use_types" %in% names(formals(factorana::define_model_component)))
  define_model_component <- function(..., use_types = NULL) factorana::define_model_component(...)
if (!exists("fix_factor_param", asNamespace("factorana"))) fix_factor_param <- function(fm, ...) fm   # all loadings free: variances fixed at 1 by the identification rule
source("model_r.R"); dat <- read.csv("data/mini_20000.csv"); ms <- build_mini_system(dat, 1)
init <- factorana:::initialize_parameters(ms, dat, verbose = FALSE); p <- init$init_params
est <- readRDS("results/r_s1_20000.rds")$fit$estimates; shared <- intersect(names(p), names(est)); p[shared] <- est[shared]
fm_ptr <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8L, init_params = p); pf <- extract_free_params_cpp(fm_ptr, p)
tm <- function(f, r = reps) median(vapply(seq_len(r), function(i) system.time(f())[["elapsed"]], numeric(1)))
cat(sprintf("%s (%s): nparam %d, free %d, %d of %d set to the converged estimates, loglik %.8f\n", label, ver, length(p), length(pf), length(shared), length(p), evaluate_loglik_only_cpp(fm_ptr, pf)))
cat(sprintf("%s (%s): loglik %.3f s   gradient %.3f s   hessian %.3f s   (medians of %d, Hessian once)\n", label, ver,
            tm(function() evaluate_loglik_only_cpp(fm_ptr, pf)), tm(function() evaluate_likelihood_cpp(fm_ptr, pf, TRUE, FALSE)),
            tm(function() evaluate_likelihood_cpp(fm_ptr, pf, TRUE, TRUE), 1L), reps))
