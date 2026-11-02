# T6 summary: bias and RMSE of the education-stage estimators across replications, in units of the full-quadrature SE.
# Estimators: full (legacy, full quadrature), legacyAdapt_MLscores (legacy production), legacyAdapt_Rscores,
#             Radapt_Rscores (R pipeline), Radapt_MLscores. Usage: Rscript t6_summary.R [N]
args <- commandArgs(trailingOnly = TRUE); N <- if (length(args) >= 1) as.integer(args[1]) else 2001L
source("legacy_names.R"); source("truth_params.R")
reps <- sort(as.integer(sub("rep", "", list.files("results/t6", "^rep[0-9]+$"))))
read_leg <- function(dir) { tab <- legacy_param_table(file.path(dir, "run.log")); par <- read.table(file.path(dir, "fullmodel_par.txt"), col.names = c("idx","est","se"))
  m <- merge(tab[!tab$fixed, c("idx","param")], par, by = "idx"); setNames(m$est, m$param) }
read_leg_se <- function(dir) { tab <- legacy_param_table(file.path(dir, "run.log")); par <- read.table(file.path(dir, "fullmodel_par.txt"), col.names = c("idx","est","se"))
  m <- merge(tab[!tab$fixed, c("idx","param")], par, by = "idx"); setNames(m$se, m$param) }
read_r <- function(f) { d <- read.csv(f); setNames(d$est, d$param) }
rows <- list(); nonconv <- c(R_L = 0, R_R = 0)
for (r in reps) {
  d <- sprintf("results/t6/rep%d", r); if (!file.exists(file.path(d, "r_adapt_R.csv"))) next
  tr <- readRDS(file.path(d, sprintf("truth_%d.rds", N)))$truth
  full <- read_leg(file.path(d, sprintf("legacy_s2_%d", N))); se_full <- read_leg_se(file.path(d, sprintf("legacy_s2_%d", N)))
  ed <- names(full)[!grepl("^(g_|cog|misslead|lead|stresstol|wage|factor_var)", names(full))]
  truth <- truth_r_params(tr, ed, subset = TRUE)
  est <- list(full = full[ed], legacyAdapt_MLscores = read_leg(file.path(d, sprintf("legacy_s2aL_%d", N)))[ed],
              legacyAdapt_Rscores = read_leg(file.path(d, sprintf("legacy_s2aR_%d", N)))[ed],
              Radapt_Rscores = read_r(file.path(d, "r_adapt_R.csv"))[ed], Radapt_MLscores = read_r(file.path(d, "r_adapt_L.csv"))[ed])
  cv <- read.csv(file.path(d, "r_adapt_L.csv"))$conv[1]; if (cv != 0) nonconv["R_L"] <- nonconv["R_L"] + 1
  cv <- read.csv(file.path(d, "r_adapt_R.csv"))$conv[1]; if (cv != 0) nonconv["R_R"] <- nonconv["R_R"] + 1
  for (nm in names(est)) rows[[length(rows) + 1]] <- data.frame(rep = r, estimator = nm, param = ed, z = (est[[nm]] - truth) / se_full[ed],
                                                                 dev_full = (est[[nm]] - full[ed]) / se_full[ed])
}
d <- do.call(rbind, rows); cat("replications used:", length(unique(d$rep)), " education parameters:", length(unique(d$param)), "\n")
cat("R adaptive non-convergence counts:", paste(names(nonconv), nonconv, collapse = ", "), "\n")
by_est <- function(f) do.call(rbind, lapply(split(d, d$estimator), function(x) { pm <- tapply(x$z, x$param, mean); data.frame(estimator = x$estimator[1],
  mean_bias_SE = round(mean(pm), 3), mean_abs_bias_SE = round(mean(abs(pm)), 3), RMSE_SE = round(sqrt(mean(x$z^2)), 3),
  mean_abs_dev_from_full_SE = round(mean(abs(x$dev_full)), 3), max_abs_dev_from_full_SE = round(max(abs(x$dev_full)), 3)) }))
res <- by_est(); res <- res[order(match(res$estimator, c("full","legacyAdapt_MLscores","legacyAdapt_Rscores","Radapt_Rscores","Radapt_MLscores"))), ]
print(res, row.names = FALSE)
cat("\nBias per parameter (mean z over replications) with |bias| > 0.5 SE for any adaptive estimator:\n")
pm <- tapply(d$z, list(d$param, d$estimator), mean); flag <- apply(abs(pm[, colnames(pm) != "full", drop = FALSE]), 1, max) > 0.5
if (any(flag)) print(round(pm[flag, , drop = FALSE], 2)) else cat("(none)\n")
write.csv(d, sprintf("results/t6_long_%d.csv", N), row.names = FALSE); write.csv(res, sprintf("results/t6_summary_%d.csv", N), row.names = FALSE)
