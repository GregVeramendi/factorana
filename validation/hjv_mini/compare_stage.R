# Usage: Rscript compare_stage.R <r_stage 1|3|4> <legacy_stage 1|2> <N>
# Joins R estimates to legacy estimates by exact parameter name and reports the gaps.
args <- commandArgs(trailingOnly = TRUE)
rs <- args[1]; ls_ <- args[2]; N <- as.integer(args[3])   # stage numbers or tags such as s3aL / s2aL
group <- if (length(args) >= 4) as.integer(args[4]) else NA   # legacy estimation group whose objective to use
rtag <- if (grepl("^[0-9]+$", rs)) paste0("s", rs) else rs; ltag <- if (grepl("^[0-9]+$", ls_)) paste0("s", ls_) else ls_
r <- read.csv(sprintf("results/r_%s_%d.csv", rtag, N))
rr <- readRDS(sprintf("results/r_%s_%d.rds", rtag, N))
ldir <- sprintf("results/legacy_%s_%d", ltag, N)
log <- readLines(file.path(ldir, "run.log"))
par <- read.table(file.path(ldir, if (grepl("^s1", ltag)) "meas_par.txt" else "fullmodel_par.txt"), col.names = c("idx", "est", "se"))

source("legacy_names.R")
leg <- legacy_param_table(file.path(ldir, "run.log"))
leg <- merge(leg[, c("idx", "param", "fixed")], par, by = "idx")   # est/se from the parameter file, full precision
leg <- leg[!leg$fixed, ]
j <- merge(r, leg[, c("param", "est", "se")], by = "param", suffixes = c("_r", "_l"))
r_free <- r$param[!is.na(r$se) & r$se > 0]
unmatched_r <- setdiff(r_free, j$param); unmatched_l <- setdiff(leg$param, j$param)
j$d_est <- j$est_r - j$est_l; j$z <- j$d_est / j$se_l; j$se_ratio <- j$se_r / j$se_l
cat(sprintf("matched %d parameters; unmatched R: %d, unmatched legacy: %d\n", nrow(j), length(unmatched_r), length(unmatched_l)))
if (length(unmatched_r)) cat("  R only:", paste(head(unmatched_r, 12), collapse = ", "), "\n")
if (length(unmatched_l)) cat("  legacy only:", paste(head(unmatched_l, 12), collapse = ", "), "\n")
free <- j[j$se_l > 0, ]
cat(sprintf("estimates: max |diff| = %.2e (%s); max |diff|/SE = %.4f\n", max(abs(free$d_est)), free$param[which.max(abs(free$d_est))], max(abs(free$z))))
cat(sprintf("SEs:       max |ratio-1| = %.2e (%s)\n", max(abs(free$se_ratio - 1)), free$param[which.max(abs(free$se_ratio - 1))]))
# Ipopt prints the objective to 17 digits; "Final Likelihood=" has only one decimal.
obj <- grep("^Objective\\.+:", log, value = TRUE)
pick <- if (is.na(group)) length(obj) else group
lk_l <- as.numeric(strsplit(trimws(sub("^Objective\\.+:", "", obj[pick])), " +")[[1]][1])
nfac <- 3; const <- N * nfac * 0.5 * log(pi)
cat(sprintf("loglik: R = %.6f, legacy objective = %.6f, R(-loglik) - legacy - N*nfac*0.5*log(pi) = %.3e\n",
            rr$fit$loglik, lk_l, -rr$fit$loglik - lk_l - const))
cat(sprintf("pass (|diff|/SE <= 0.01, SE ratio within 1e-3, loglik gap <= 1e-4): %s\n",
            all(abs(free$z) <= 0.01) && all(abs(free$se_ratio - 1) <= 1e-3) && abs(-rr$fit$loglik - lk_l - const) <= 1e-4))
write.csv(j[order(-abs(j$z)), ], sprintf("results/compare_%s_vs_%s_%d.csv", rtag, ltag, N), row.names = FALSE)
print(head(j[order(-abs(j$z)), c("param", "est_r", "est_l", "se_r", "se_l", "z")], 8), row.names = FALSE, digits = 7)
