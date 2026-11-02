# T9: Hessians and SEs at the converged typed optimum. Writes legacy run B's estimates as init_par.txt for the typed
# TestCode, runs it (full dump), evaluates R at the same point, and compares the outcome-block Hessians and the SEs
# from their inverses with each other and with the SEs both codes reported. Usage: Rscript t9_hess_at_opt.R <N> <T>
args <- commandArgs(trailingOnly = TRUE); N <- as.integer(args[1]); T <- as.integer(args[2])
suppressMessages(pkgload::load_all("../../factorana_R", compile = FALSE, quiet = TRUE)); source("model_r.R"); source("legacy_names.R")
dtag <- sprintf("_t%d", T); d0 <- sprintf("results/legacy_s0t%d_%d", T, N); dB <- sprintf("results/legacy_s2t%dB_%d", T, N)
dat <- read.csv(sprintf("data/mini_%d%s.csv", N, dtag)); fm <- mini_factor_model(T)
ms <- define_model_system(components = c(build_mini_system(dat, 1, ntypes = T)$components, build_mini_system(dat, 5, ntypes = T)$components), factor = fm)
init <- factorana:::initialize_parameters(ms, dat, verbose = FALSE)
tab <- legacy_param_table(file.path(dB, "run.log")); par <- read.table(file.path(dB, "fullmodel_par.txt"), col.names = c("idx", "est", "se"))
leg <- merge(tab[, c("idx", "param", "fixed")], par, by = "idx")
p <- setNames(rep(NA_real_, length(init$param_names)), init$param_names); p[leg$param] <- leg$est; p[grep("^factor_var", names(p))] <- 1
stopifnot(!anyNA(p))
# legacy TestCode at this point
map <- read.csv(file.path(d0, "legacy_index_map.csv")); file.copy(file.path(d0, "init_par.txt"), file.path(d0, "init_par_truth.txt"), overwrite = TRUE)
write.table(data.frame(map$idx, sprintf("%.12f", p[map$param]), 0.1), file.path(d0, "init_par.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
rc <- system(sprintf("./prob_t9_s0t%d > %s/run_at_opt.log 2>&1", T, d0)); stopifnot(rc == 0)
dump <- scan(file.path(d0, "t1_dump.txt"), quiet = TRUE); free_l <- scan(file.path(d0, "t1_freeparlist.txt"), quiet = TRUE); nf <- length(free_l)
names_l <- map$param[match(free_l, map$idx)]; hess_l <- dump[-(1:(1 + nf))]
file.copy(file.path(d0, "init_par_truth.txt"), file.path(d0, "init_par.txt"), overwrite = TRUE)
# R at the same point
meta <- factorana:::build_parameter_metadata(ms); cons <- factorana:::setup_parameter_constraints(ms, p, meta, init$factor_variance_fixed, verbose = FALSE)
fm_ptr <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8L, init_params = p); fi <- cons$free_idx; fn <- names(p)[fi]
res <- evaluate_likelihood_cpp(fm_ptr, p[fi], TRUE, TRUE)
unpack <- function(h, n) { H <- matrix(0, n, n); k <- 1; for (i in 1:n) for (j in i:n) { H[i, j] <- H[j, i] <- h[k]; k <- k + 1 }; H }
H_r <- unpack(res$hessian, length(fi)); dimnames(H_r) <- list(fn, fn)
H_l <- -unpack(hess_l, nf); dimnames(H_l) <- list(names_l, names_l); H_l <- H_l[fn, fn]      # legacy: Hessian of -logL
out <- !grepl("^(g_|cog|misslead|lead|stresstol|factor_var)", fn)
D <- H_r[out, out] - H_l[out, out]; big <- abs(H_r[out, out]) > 1e-6 * max(abs(H_r[out, out]))
cat(sprintf("gradient at optimum: max |R| = %.2e (outcome params)\n", max(abs(res$gradient[out]))))
cat(sprintf("Hessian at optimum, %d outcome parameters: elements with |R - legacy| > 1e-6: %d of %d; max |diff| = %.3e; max relative gap on non-tiny = %.2e\n",
            sum(out), sum(abs(D) > 1e-6), length(D), max(abs(D)), max(abs(D[big]) / abs(H_r[out, out][big]))))
idx <- which(abs(D) > 1e-4, arr.ind = TRUE)
if (nrow(idx)) { o <- order(-abs(D[idx]))[seq_len(min(12, nrow(idx)))]
  print(data.frame(row = rownames(D)[idx[o, 1]], col = colnames(D)[idx[o, 2]], H_R = round(H_r[out, out][idx[o, , drop = FALSE]], 4), H_legacy = round(H_l[out, out][idx[o, , drop = FALSE]], 4)), row.names = FALSE) }
# SEs from the two Hessians (outcome block; measurement parameters are fixed in the final stage)
se_r <- sqrt(diag(solve(-H_r[out, out]))); se_l <- sqrt(diag(solve(-H_l[out, out])))
rep_l <- setNames(leg$se, leg$param)[fn[out]]; fitR <- readRDS(sprintf("results/r_s5t%dT_%d.rds", T, N))$fit; rep_r <- fitR$std_errors[fn[out]]
cat(sprintf("SEs from inverse Hessians at the optimum: max |SE_legacyHess / SE_RHess - 1| = %.2e (%s)\n", max(abs(se_l / se_r - 1)), fn[out][which.max(abs(se_l / se_r - 1))]))
cat(sprintf("legacy's REPORTED SEs vs SEs from legacy's own Hessian at the optimum: max |ratio - 1| = %.2e\n", max(abs(rep_l / se_l - 1))))
cat(sprintf("legacy's reported SEs vs R's reported SEs: max |ratio - 1| = %.2e; R reported vs R Hessian here: %.2e\n", max(abs(rep_l / rep_r - 1)), max(abs(rep_r / se_r - 1))))
aff <- c("wage1_type_3_intercept", "wage1_x1", "swesat_type_3_intercept", "swesat_x1", "wage2_type_3_intercept", "wage2_x1", "wage2_type_2_intercept", "wage2_intercept", "swesat_type_2_intercept", "swesat_intercept", "wage1_type_2_intercept", "wage1_intercept")
print(data.frame(param = aff, se_legacy_reported = signif(rep_l[aff], 5), se_R = signif(rep_r[aff], 5), ratio = round(rep_l[aff] / rep_r[aff], 6)), row.names = FALSE)
