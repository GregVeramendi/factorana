# Write R factor scores in the legacy factor_predictions.txt format (id, then score/SE per factor), for the T5 cross run.
args <- commandArgs(trailingOnly = TRUE); N <- as.integer(args[1]); out <- args[2]
fs <- readRDS(sprintf("results/r_scores_%d.rds", N))$obs_se$scores
write.table(data.frame(fs$id, fs$factor_1, fs$se_factor_1, fs$factor_2, fs$se_factor_2, fs$factor_3, fs$se_factor_3),
            out, row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("wrote", out, nrow(fs), "rows\n")
