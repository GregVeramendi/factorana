# Write init_par.txt (idx value err) for the legacy stage-0 (TestCode) run at the TRUE parameters.
# Usage: Rscript make_init_par.R <N> [logfile] [outdir] [expected rows, 0 = no check] [data tag, e.g. _t3].
# Names come from the legacy stage-2 log; values from truth_params.R.
args <- commandArgs(trailingOnly = TRUE); N <- as.integer(args[1])
logfile <- if (length(args) >= 2) args[2] else sprintf("results/legacy_s2_%d/run.log", N)
out <- if (length(args) >= 3) args[3] else sprintf("results/legacy_s0_%d", N)
expected <- if (length(args) >= 4) as.integer(args[4]) else 294L
dtag <- if (length(args) >= 5) args[5] else ""
suppressMessages(pkgload::load_all("../../factorana_R", compile = FALSE, quiet = TRUE))
source("legacy_names.R"); source("truth_params.R")
tab <- legacy_param_table(logfile)
tr <- readRDS(sprintf("data/truth_%d%s.rds", N, dtag))$truth
p <- truth_r_params(tr, tab$param, subset = TRUE)   # stops if any legacy name is not a model parameter
stopifnot(expected <= 0 || nrow(tab) == expected, !any(duplicated(tab$idx)), all(tab$idx == seq_len(nrow(tab)) - 1))
dir.create(out, showWarnings = FALSE)
write.table(data.frame(tab$idx, sprintf("%.12f", p[tab$param]), 0.1), file.path(out, "init_par.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("wrote", file.path(out, "init_par.txt"), "with", nrow(tab), "rows;", sum(tab$fixed), "fixed in legacy\n")
write.csv(tab[, c("idx", "param", "fixed")], file.path(out, "legacy_index_map.csv"), row.names = FALSE)
