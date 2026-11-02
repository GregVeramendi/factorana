# Tiny 1-factor, 3-measure problem used only to check that the legacy build runs.
set.seed(7); n <- 2000
f <- rnorm(n); x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4)
truth <- list(m1 = c(b0 = 0.5, b1 = 0.3, b2 = -0.2, lam = 1.0, sig = 0.6),
              m2 = c(b0 = -0.1, b1 = 0.2, b2 = 0.4, lam = 0.8, sig = 0.7),
              m3 = c(b0 = 0.2, b1 = -0.3, b2 = 0.1, lam = 1.2, sig = 0.5))
mk <- function(p) p["b0"] + p["b1"] * x1 + p["b2"] * x2 + p["lam"] * f + rnorm(n, 0, p["sig"])
d <- data.frame(M1 = mk(truth$m1), M2 = mk(truth$m2), M3 = mk(truth$m3), x1 = x1, x2 = x2, const = 1)
write.table(d, "tiny.raw", row.names = FALSE, col.names = FALSE)
writeLines(c(ncol(d), names(d)), "tiny_varlist.txt")
write.csv(d, "tiny.csv", row.names = FALSE)
saveRDS(list(truth = truth, f = f), "tiny_truth.rds")
suppressMessages(pkgload::load_all("../../../factorana_R", compile = FALSE, quiet = TRUE))
fm <- define_factor_model(n_factors = 1)
mc <- lapply(c("M1", "M2", "M3"), function(y) define_model_component(
  tolower(y), d, y, fm, covariates = c("const", "x1", "x2"), model_type = "linear",
  loading_normalization = if (y == "M1") 1 else NA_real_))
ms <- define_model_system(components = mc, factor = fm)
ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
t0 <- Sys.time()
fit <- estimate_model_rcpp(ms, d, control = ctrl, optimizer = "nlminb", parallel = FALSE, verbose = FALSE)
cat(sprintf("R: convergence=%d loglik=%.6f secs=%.1f\n", fit$convergence, fit$loglik, as.numeric(Sys.time() - t0, units = "secs")))
out <- data.frame(param = names(fit$estimates), est = round(fit$estimates, 5), se = round(fit$std_errors, 5))
write.csv(out, "results/r_estimates.csv", row.names = FALSE); print(out, row.names = FALSE)
