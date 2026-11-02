# T7b, R side: per-evaluation timing for a chosen stage system, single core, at the truth, with a chosen package
# library (current working copy, or an installed older version). Comparable to legacy TestCode with MEAS_ONLY (stage 1).
# Usage: Rscript t7b_r_timing.R <N> <stage 1|all> <repeats> [lib_dir or ""] [ntypes]
# ntypes > 1 (T9): full system = measurement + joint education/labor stage with latent types, data suffix _t<T>.
args <- commandArgs(trailingOnly = TRUE); N <- as.integer(args[1]); stage <- args[2]
reps <- if (length(args) >= 3) as.integer(args[3]) else 3L; lib <- if (length(args) >= 4) args[4] else ""
ntypes <- if (length(args) >= 5) as.integer(args[5]) else 1L; dtag <- if (ntypes > 1) sprintf("_t%d", ntypes) else ""
if (nzchar(lib)) { suppressMessages(library(factorana, lib.loc = lib)); ver <- as.character(packageVersion("factorana", lib.loc = lib)) } else {
  suppressMessages(pkgload::load_all("../../factorana_R", compile = FALSE, quiet = TRUE)); ver <- paste0("working copy ", read.dcf("../../factorana_R/DESCRIPTION", "Version")) }
source("model_r.R"); source("truth_params.R")
# Packages before 1.3.0 have no fix_factor_param(); the variances are then free parameters held at 1 in the
# evaluated vector, so the objective is the same and the free-parameter count is 3 higher.
if (!exists("fix_factor_param")) mini_factor_model <- function(ntypes = 1) define_factor_model(n_factors = 3, n_types = ntypes)
dat <- read.csv(sprintf("data/mini_%d%s.csv", N, dtag)); tr <- readRDS(sprintf("data/truth_%d%s.rds", N, dtag))$truth
fm <- mini_factor_model(ntypes)
comps <- if (stage == "1") build_mini_system(dat, 1, ntypes = ntypes)$components else
  c(build_mini_system(dat, 1, ntypes = ntypes)$components, build_mini_system(dat, 5, ntypes = ntypes)$components)
ms <- define_model_system(components = comps, factor = fm)
init <- factorana:::initialize_parameters(ms, dat); p <- truth_r_params(tr, init$param_names, subset = TRUE)
meta <- factorana:::build_parameter_metadata(ms); cons <- factorana:::setup_parameter_constraints(ms, p, meta, init$factor_variance_fixed, verbose = FALSE)
fm_ptr <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8L, init_params = p); pf <- p[cons$free_idx]
tm <- function(f) median(vapply(seq_len(reps), function(i) system.time(f())[["elapsed"]], numeric(1)))
t_l <- tm(function() evaluate_loglik_only_cpp(fm_ptr, pf))
t_g <- tm(function() evaluate_likelihood_cpp(fm_ptr, pf, TRUE, FALSE))
t_h <- tm(function() evaluate_likelihood_cpp(fm_ptr, pf, TRUE, TRUE))
ll <- evaluate_loglik_only_cpp(fm_ptr, pf)
cat(sprintf("R [%s] N=%d stage=%s types=%d (median of %d): loglik %.3f s | +grad %.3f s | +hess %.3f s | free params %d | loglik %.6f\n",
            ver, N, stage, ntypes, reps, t_l, t_g, t_h, length(pf), ll))
