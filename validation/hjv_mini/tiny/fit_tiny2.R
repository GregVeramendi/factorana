suppressMessages(pkgload::load_all("../../../factorana_R", compile = FALSE, quiet = TRUE))
d <- read.csv("tiny.csv")
fm <- fix_factor_param(define_factor_model(n_factors = 1), c(factor_var_1 = 1))
mc <- lapply(c("M1","M2","M3"), function(y) suppressWarnings(define_model_component(
  tolower(y), d, y, fm, covariates = c("const","x1","x2"), model_type = "linear", loading_normalization = NA_real_)))
ms <- define_model_system(components = mc, factor = fm)
fit <- estimate_model_rcpp(ms, d, control = define_estimation_control(n_quad_points = 8, num_cores = 1),
                           optimizer = "nlminb", parallel = FALSE, verbose = FALSE)
cat(sprintf("R: convergence=%d loglik=%.6f\n", fit$convergence, fit$loglik))
out <- data.frame(param = names(fit$estimates), est = round(fit$estimates, 6), se = round(fit$std_errors, 6))
write.csv(out, "results2/r_estimates.csv", row.names = FALSE); print(out, row.names = FALSE)
