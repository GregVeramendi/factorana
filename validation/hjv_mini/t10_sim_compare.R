# T10: legacy C++ simulation_data.csv against simulate_factor_model(layout = "legacy") on the mini model with types.
# Usage: Rscript t10_sim_compare.R <N> <T> <NSIM>, after run_t10.sh has run the C++ simulation. Checks: every legacy
# parameter maps to an R parameter through legacy_param_table(), the two files have the same header and the same
# field formats, and the column means, standard deviations and missing shares agree within sampling error.
args <- commandArgs(trailingOnly = TRUE); N <- as.integer(args[1]); T <- as.integer(args[2]); NSIM <- as.integer(args[3])
suppressMessages(pkgload::load_all("../../factorana_R", compile = FALSE, quiet = TRUE))
source("model_r.R")
W <- sprintf("results/legacy_s4t%d_%d", T, N)
dat <- read.csv(sprintf("data/mini_%d_t%d.csv", N, T))

# ---- model: measurement + education + labor market in one system, 3 types ----
ms1 <- build_mini_system(dat, 1, ntypes = T)
p1 <- initialize_parameters(ms1, dat, verbose = FALSE)$init_params
fake <- structure(list(model_system = ms1, estimates = p1, std_errors = p1 * 0), class = "factorana_result")
ms <- suppressMessages(build_mini_system(dat, 5, previous = fake, ntypes = T))
cn <- vapply(ms$components, function(c) c$name, "")
cat("components:", length(cn), "\n")
endo <- list(g_swe10 = "HS_track", g_sports10 = "HS_track", g_gpahs = "HS_track",
             HS_track = c("advmath", "adveng"), wage1 = "HS_track", wage2 = "HS_track")
for (nm in c("apply", "application", "Dswesat", "swesat", "enroll", "switch_major", "grad1", "grad2", "grad3"))
  endo[[nm]] <- c("advmath", "adveng", "HS_track")
for (nm in names(endo)) ms$components[[match(nm, cn)]]$endogenous <- endo[[nm]]

# ---- parameters: legacy fullmodel_par.txt mapped to R names by (model, variable) ----
vl <- read.table(file.path(W, "fullmodel_par_varlist.txt"), col.names = c("idx", "model", "variable"))
fp <- read.table(file.path(W, "fullmodel_par.txt"), col.names = c("idx", "estimate", "se"))
leg <- merge(vl, fp, by = "idx")
p0 <- initialize_parameters(ms, dat, verbose = FALSE)$init_params
tab <- as.data.frame(legacy_param_table(ms, p0, type_intercept_name = "intercept"))
key <- function(d) { k <- paste(d$model, d$variable); paste(k, ave(seq_along(k), k, FUN = seq_along)) }  # k-th occurrence
m <- match(key(tab), key(leg))
cat("R parameters:", nrow(tab), " legacy parameters:", nrow(leg), " R rows without a legacy match:", sum(is.na(m)),
    " legacy rows without an R match:", sum(!key(leg) %in% key(tab)), "\n")
if (any(is.na(m))) print(tab[is.na(m), ])
if (any(!key(leg) %in% key(tab))) print(leg[!key(leg) %in% key(tab), ])
params <- p0
val <- leg$estimate[m]
is_fv <- grepl("^factor_var_", tab$parameter)
val[is_fv] <- val[is_fv]^2                      # legacy sigma -> variance
ok <- !is.na(m) & tab$parameter %in% names(params)
params[tab$parameter[ok]] <- val[ok]

# ---- R simulation in the legacy layout ----
csv <- "results/t10/r_simulation_data.csv"
t0 <- Sys.time()
L <- simulate_factor_model(ms, dat, n = NSIM, params = params, seed = 20261110, layout = "legacy",
                           index_var = "id", component_order = c("advmath", "adveng", "HS_track"),
                           csv_file = csv)
cat(sprintf("R simulation + CSV: %.1f s, %d x %d, csv %.0f MB\n", as.numeric(difftime(Sys.time(), t0, units = "secs")),
            nrow(L), ncol(L), file.size(csv) / 1e6))

# ---- 1. header ----
hc <- readLines(file.path(W, "simulation_data.csv"), n = 1); hr <- readLines(csv, n = 1)
cat("header identical:", identical(hc, hr), "\n")
if (!identical(hc, hr)) {
  a <- strsplit(hc, ", ")[[1]]; b <- strsplit(hr, ", ")[[1]]
  cat(" C++ only:", head(setdiff(a, b), 20), "\n R only:", head(setdiff(b, a), 20), "\n same set:", setequal(a, b), "\n")
  k <- which(a[seq_len(min(length(a), length(b)))] != b[seq_len(min(length(a), length(b)))])[1]
  cat(" first position that differs:", k, a[k], b[k], "\n")
}
# ---- 2. line format ----
lc <- readLines(file.path(W, "simulation_data.csv"), n = 2001)[-1]; lr <- readLines(csv, n = 2001)[-1]
shape <- function(x) gsub("[0-9]", "9", x)
fmt_cols <- function(lines) {     # per column: field widths and whether values are written as integers or with 5 decimals
  f <- do.call(rbind, strsplit(sub(" $", "", lines), ", ", fixed = TRUE))
  apply(f, 2, function(v) paste(paste(sort(unique(nchar(v))), collapse = "/"),
                                paste(sort(unique(sub("^ *-?[0-9]+", "I", sub("\\.[0-9]+$", ".D", v)))), collapse = "|")))
}
fc <- fmt_cols(lc); fr <- fmt_cols(lr)
cat("field formats in the C++ file:"); print(table(fc))
cat("line widths C++:", paste(range(nchar(lc)), collapse = "-"), " R:", paste(range(nchar(lr)), collapse = "-"),
    "  trailing space C++/R:", all(grepl(" $", lc)), all(grepl(" $", lr)), "\n")
cat("columns whose field shapes differ in the first 2000 lines:", sum(fc != fr), "of", length(fc), "\n")
if (any(fc != fr)) { h <- trimws(strsplit(hc, ",")[[1]]); print(head(data.frame(col = h[fc != fr], cpp = fc[fc != fr], r = fr[fc != fr]), 12)) }
# ---- 3. moments ----
Cx <- data.table::fread(file.path(W, "simulation_data.csv"), header = TRUE, sep = ",", data.table = FALSE, showProgress = FALSE)
names(Cx) <- trimws(names(Cx)); Rx <- as.data.frame(L)
stopifnot(setequal(names(Cx), names(Rx))); Rx <- Rx[names(Cx)]
Cx[Cx == -9999] <- NA
res <- do.call(rbind, lapply(names(Cx), function(v) {
  a <- Cx[[v]]; b <- Rx[[v]]
  na_a <- mean(is.na(a)); na_b <- mean(is.na(b)); a <- a[!is.na(a)]; b <- b[!is.na(b)]
  se <- sqrt(var(a) / length(a) + var(b) / length(b))
  data.frame(col = v, na_cpp = na_a, na_r = na_b, mean_cpp = mean(a), mean_r = mean(b), sd_cpp = sd(a), sd_r = sd(b),
             z_mean = if (se > 0) (mean(a) - mean(b)) / se else 0,
             z_na = { p <- (na_a + na_b) / 2; s <- sqrt(2 * p * (1 - p) / NSIM); if (s > 0) (na_a - na_b) / s else 0 },
             # Kolmogorov-Smirnov on the values rounded as the files hold them, continuous columns only
             ks = if (length(unique(round(a, 5))) > 1000) suppressWarnings(ks.test(round(a, 5), round(b, 5))$p.value) else NA)
}))
saveRDS(res, "results/t10/moments.rds")
cat("columns compared:", nrow(res), "\n")
cat("max |z| of mean difference:", round(max(abs(res$z_mean)), 2), "  columns with |z| > 3:", sum(abs(res$z_mean) > 3),
    " (expected by chance at 0.27%:", round(0.0027 * nrow(res), 1), ")\n")
cat("max |z| of missing-share difference:", round(max(abs(res$z_na)), 2), "\n")
cat("SD ratio R/C++ range (columns with sd > 0):", paste(round(range((res$sd_r / res$sd_cpp)[res$sd_cpp > 0]), 4), collapse = " to "), "\n")
cat("KS p-values (continuous columns): n =", sum(!is.na(res$ks)), " min =", signif(min(res$ks, na.rm = TRUE), 3),
    " share below 0.01 =", round(mean(res$ks < 0.01, na.rm = TRUE), 4), "\n")
print(head(res[order(-abs(res$z_mean)), c("col", "mean_cpp", "mean_r", "sd_cpp", "sd_r", "z_mean")], 12), row.names = FALSE, digits = 4)
miss <- grep("_miss$", names(Cx), value = TRUE)
mt <- t(sapply(miss, function(v) c(cpp = tabulate(Cx[[v]] + 1L, 3), r = tabulate(Rx[[v]] + 1L, 3))))
colnames(mt) <- c("cpp_0", "cpp_1", "cpp_2", "r_0", "r_1", "r_2"); print(mt[rowSums(mt[, c(2, 3, 5, 6)]) > 0, , drop = FALSE])
