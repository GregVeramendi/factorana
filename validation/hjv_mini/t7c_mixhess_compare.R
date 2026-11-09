args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
b <- readRDS(file.path(dir, "res_correctness_before.rds")); a <- readRDS(file.path(dir, "res_correctness_after.rds"))
cat(sprintf("%-28s %6s %6s %12s %12s %12s %9s\n", "scenario", "nparam", "nfree", "d_loglik", "max|d_grad|", "max|d_hess|", "identical"))
for (nm in names(b)) {
  x <- b[[nm]]; y <- a[[nm]]
  if (!is.null(x$error) || !is.null(y$error)) { cat(sprintf("%-28s ERROR before: %s | after: %s\n", nm, x$error, y$error)); next }
  dl <- max(abs(c(x$loglik_1 - y$loglik_1, x$loglik_2 - y$loglik_2, x$loglik_3 - y$loglik_3)))
  dg <- max(abs(c(x$gradient_2 - y$gradient_2, x$gradient_3 - y$gradient_3)))
  dh <- max(abs(x$hessian - y$hessian))
  cat(sprintf("%-28s %6d %6d %12.3e %12.3e %12.3e %9s\n", nm, x$nparam, x$nfree, dl, dg, dh,
              identical(x[c("loglik_1","loglik_2","loglik_3","gradient_2","gradient_3","hessian")],
                        y[c("loglik_1","loglik_2","loglik_3","gradient_2","gradient_3","hessian")])))
}
