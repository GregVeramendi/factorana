# T5: adaptive estimates against the full-quadrature estimates (R stage 3 and legacy group 1), education parameters only.
args <- commandArgs(trailingOnly = TRUE); N <- as.integer(args[1]); src <- args[2]
source("legacy_names.R")
full_r <- read.csv(sprintf("results/r_s3_%d.csv", N)); full_r <- full_r[!is.na(full_r$se) & full_r$se > 0, ]
ad_r <- read.csv(sprintf("results/r_s3a%s_%d.csv", src, N))
lt <- legacy_param_table(sprintf("results/legacy_s2a%s_%d/run.log", src, N))
lp <- read.table(sprintf("results/legacy_s2a%s_%d/fullmodel_par.txt", src, N), col.names = c("idx","est","se"))
ad_l <- merge(lt[!lt$fixed, c("idx","param")], lp, by = "idx")
meas <- grepl("^(g_|cog|misslead|lead|stresstol|wage)", full_r$param)
ed <- full_r[!meas, ]
j <- merge(ed, ad_r, by = "param", suffixes = c("_full", "_ar")); j <- merge(j, ad_l[, c("param","est")], by = "param"); names(j)[ncol(j)] <- "est_al"
cat(sprintf("education parameters: %d\n", nrow(j)))
cat(sprintf("R adaptive  vs full: max |diff|/SE_full = %.3f (%s); mean |diff|/SE = %.3f\n", max(abs(j$est_ar - j$est_full) / j$se_full), j$param[which.max(abs(j$est_ar - j$est_full) / j$se_full)], mean(abs(j$est_ar - j$est_full) / j$se_full)))
cat(sprintf("legacy adaptive vs full: max |diff|/SE_full = %.3f (%s); mean |diff|/SE = %.3f\n", max(abs(j$est_al - j$est_full) / j$se_full), j$param[which.max(abs(j$est_al - j$est_full) / j$se_full)], mean(abs(j$est_al - j$est_full) / j$se_full)))
cat(sprintf("R adaptive vs legacy adaptive (same scores %s): max |diff|/SE_full = %.4f\n", src, max(abs(j$est_ar - j$est_al) / j$se_full)))
write.csv(j, sprintf("results/t5_vs_full_%s_%d.csv", src, N), row.names = FALSE)
