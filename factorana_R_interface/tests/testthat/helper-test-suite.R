# Helper functions for systematic test suite
# These functions support gradient checking, Hessian checking, and estimation comparison

#' Wrapper for evaluating likelihood with model_system
#'
#' @param model_system Model system object
#' @param data Data frame
#' @param params Parameter vector
#' @param compute_gradient Logical
#' @param compute_hessian Logical
#' @return List with loglik, gradient (if requested), hessian (if requested)
evaluate_likelihood_rcpp <- function(model_system, data, params,
                                     compute_gradient = FALSE, compute_hessian = FALSE) {
  # Convert data to matrix
  data_mat <- as.matrix(data)
  n_quad <- model_system$factor$n_quad

  # Initialize C++ model
  fm_ptr <- initialize_factor_model_cpp(model_system, data_mat, n_quad)

  # Evaluate
  result <- evaluate_likelihood_cpp(fm_ptr, params,
                                   compute_gradient = compute_gradient,
                                   compute_hessian = compute_hessian)

  # Standardize field names (C++ returns 'logLikelihood', we want 'loglik')
  if (!is.null(result$logLikelihood)) {
    result$loglik <- result$logLikelihood
  }

  return(result)
}

#' Compute finite difference gradient
#'
#' @param model_system Model system object
#' @param data Data frame
#' @param params Parameter vector
#' @param h Step size (if NULL, uses adaptive step size)
#' @return Gradient vector computed via finite differences
finite_diff_gradient <- function(model_system, data, params, h = NULL) {
  n_params <- length(params)
  grad <- numeric(n_params)

  for (i in seq_len(n_params)) {
    # Adaptive step size: sqrt(machine epsilon) * max(|param|, 1)
    if (is.null(h)) {
      h_i <- sqrt(.Machine$double.eps) * max(abs(params[i]), 1)
    } else {
      h_i <- h
    }

    # Forward and backward steps
    params_fwd <- params
    params_bwd <- params
    params_fwd[i] <- params[i] + h_i
    params_bwd[i] <- params[i] - h_i

    # Evaluate likelihood at both points
    result_fwd <- tryCatch({
      evaluate_likelihood_rcpp(model_system, data, params_fwd,
                               compute_gradient = FALSE, compute_hessian = FALSE)
    }, error = function(e) list(loglik = NA))

    result_bwd <- tryCatch({
      evaluate_likelihood_rcpp(model_system, data, params_bwd,
                               compute_gradient = FALSE, compute_hessian = FALSE)
    }, error = function(e) list(loglik = NA))

    # Central difference
    if (!is.na(result_fwd$loglik) && !is.na(result_bwd$loglik)) {
      grad[i] <- (result_fwd$loglik - result_bwd$loglik) / (2 * h_i)
    } else {
      grad[i] <- NA
    }
  }

  return(grad)
}

#' Compute finite difference Hessian
#'
#' @param model_system Model system object
#' @param data Data frame
#' @param params Parameter vector
#' @param h Step size (if NULL, uses adaptive step size)
#' @return Hessian matrix computed via finite differences
finite_diff_hessian <- function(model_system, data, params, h = NULL) {
  n_params <- length(params)
  hess <- matrix(0, n_params, n_params)

  # Get function value at base point
  result_base <- tryCatch({
    evaluate_likelihood_rcpp(model_system, data, params,
                             compute_gradient = TRUE, compute_hessian = FALSE)
  }, error = function(e) list(loglik = NA, gradient = rep(NA, n_params)))

  if (is.na(result_base$loglik)) {
    return(matrix(NA, n_params, n_params))
  }

  f_base <- result_base$loglik
  grad_base <- result_base$gradient

  # Compute Hessian via finite differences of gradient
  for (i in seq_len(n_params)) {
    # Adaptive step size
    if (is.null(h)) {
      h_i <- sqrt(.Machine$double.eps) * max(abs(params[i]), 1)
    } else {
      h_i <- h
    }

    # Forward step
    params_fwd <- params
    params_fwd[i] <- params[i] + h_i

    result_fwd <- tryCatch({
      evaluate_likelihood_rcpp(model_system, data, params_fwd,
                               compute_gradient = TRUE, compute_hessian = FALSE)
    }, error = function(e) list(gradient = rep(NA, n_params)))

    if (any(is.na(result_fwd$gradient))) {
      hess[i, ] <- NA
      hess[, i] <- NA
    } else {
      # Finite difference of gradient
      hess[i, ] <- (result_fwd$gradient - grad_base) / h_i
    }
  }

  # Symmetrize
  hess <- (hess + t(hess)) / 2

  return(hess)
}

#' Check gradient accuracy against finite differences
#'
#' @param model_system Model system object
#' @param data Data frame
#' @param params Parameter vector to check at
#' @param tol Tolerance for relative error
#' @param verbose Print detailed diagnostics
#' @return List with pass (logical), max_error (numeric), diagnostics (data.frame)
check_gradient_accuracy <- function(model_system, data, params, tol = 1e-5, verbose = FALSE) {
  # Get analytical gradient
  result <- tryCatch({
    evaluate_likelihood_rcpp(model_system, data, params,
                             compute_gradient = TRUE, compute_hessian = FALSE)
  }, error = function(e) {
    if (verbose) cat("Error computing analytical gradient:", e$message, "\n")
    return(list(gradient = rep(NA, length(params))))
  })

  grad_analytical <- result$gradient

  # Get finite difference gradient
  grad_fd <- finite_diff_gradient(model_system, data, params)

  # Compute errors
  abs_error <- abs(grad_analytical - grad_fd)
  rel_error <- abs_error / pmax(abs(grad_fd), 1e-10)

  # Create diagnostics table
  diagnostics <- data.frame(
    param_idx = seq_along(params),
    param_value = params,
    analytical = grad_analytical,
    finite_diff = grad_fd,
    abs_error = abs_error,
    rel_error = rel_error,
    pass = rel_error < tol
  )

  # Check if all pass
  all_pass <- all(diagnostics$pass, na.rm = TRUE) && !any(is.na(diagnostics$pass))
  max_error <- if (all(is.na(rel_error))) NA else max(rel_error, na.rm = TRUE)

  if (verbose) {
    cat("\n=== Gradient Check ===\n")
    cat(sprintf("Max relative error: %.2e (tolerance: %.2e)\n", max_error, tol))
    cat(sprintf("Result: %s\n\n", ifelse(all_pass, "PASS", "FAIL")))
    print(diagnostics, row.names = FALSE)
    cat("\n")
  }

  return(list(
    pass = all_pass,
    max_error = max_error,
    diagnostics = diagnostics
  ))
}

#' Check Hessian accuracy against finite differences
#'
#' @param model_system Model system object
#' @param data Data frame
#' @param params Parameter vector to check at
#' @param tol Tolerance for relative error
#' @param verbose Print detailed diagnostics
#' @return List with pass (logical), max_error (numeric), diagnostics (data.frame)
check_hessian_accuracy <- function(model_system, data, params, tol = 1e-4, verbose = FALSE) {
  # Get analytical Hessian
  result <- tryCatch({
    evaluate_likelihood_rcpp(model_system, data, params,
                             compute_gradient = TRUE, compute_hessian = TRUE)
  }, error = function(e) {
    if (verbose) cat("Error computing analytical Hessian:", e$message, "\n")
    return(list(hessian = matrix(NA, length(params), length(params))))
  })

  hess_analytical <- result$hessian

  # Get finite difference Hessian
  hess_fd <- finite_diff_hessian(model_system, data, params)

  # Compute errors (element-wise)
  abs_error <- abs(hess_analytical - hess_fd)
  rel_error <- abs_error / pmax(abs(hess_fd), 1e-10)

  # Create diagnostics table (for upper triangle only)
  n_params <- length(params)
  diag_list <- list()
  idx <- 1

  for (i in seq_len(n_params)) {
    for (j in i:n_params) {
      diag_list[[idx]] <- data.frame(
        row = i,
        col = j,
        analytical = hess_analytical[i, j],
        finite_diff = hess_fd[i, j],
        abs_error = abs_error[i, j],
        rel_error = rel_error[i, j],
        pass = rel_error[i, j] < tol
      )
      idx <- idx + 1
    }
  }

  diagnostics <- do.call(rbind, diag_list)

  # Check if all pass
  all_pass <- all(diagnostics$pass, na.rm = TRUE) && !any(is.na(diagnostics$pass))
  max_error <- if (all(is.na(rel_error))) NA else max(rel_error, na.rm = TRUE)

  if (verbose) {
    cat("\n=== Hessian Check ===\n")
    cat(sprintf("Max relative error: %.2e (tolerance: %.2e)\n", max_error, tol))
    cat(sprintf("Result: %s\n\n", ifelse(all_pass, "PASS", "FAIL")))

    # Print summary statistics
    cat("Summary of relative errors:\n")
    print(summary(diagnostics$rel_error))
    cat("\n")

    # Print worst cases
    cat("Worst 10 elements:\n")
    worst <- diagnostics[order(diagnostics$rel_error, decreasing = TRUE)[1:min(10, nrow(diagnostics))], ]
    print(worst, row.names = FALSE)
    cat("\n")
  }

  return(list(
    pass = all_pass,
    max_error = max_error,
    diagnostics = diagnostics
  ))
}

#' Run estimation comparison test
#'
#' @param model_system Model system object
#' @param data Data frame
#' @param true_params True parameter vector (used for simulation)
#' @param param_names Optional parameter names
#' @param verbose Print detailed diagnostics
#' @return List with results and diagnostics
run_estimation_comparison <- function(model_system, data, true_params,
                                      param_names = NULL, verbose = FALSE) {

  if (is.null(param_names)) {
    param_names <- paste0("param_", seq_along(true_params))
  }

  # Get default initial parameters
  default_init <- initialize_parameters(model_system, data, verbose = FALSE)

  if (verbose) {
    cat("\n=== Estimation Comparison ===\n\n")
  }

  # Estimate with default initialization
  if (verbose) cat("Estimating with default initialization...\n")
  result_default <- tryCatch({
    estimate_model_rcpp(model_system, data, init_params = default_init,
                       optimizer = "nlminb", parallel = FALSE, verbose = FALSE)
  }, error = function(e) {
    if (verbose) cat("Error with default init:", e$message, "\n")
    return(list(estimates = rep(NA, length(true_params)),
                std_errors = rep(NA, length(true_params)),
                loglik = NA, convergence = 999))
  })

  # Estimate with true parameters as initialization
  if (verbose) cat("Estimating with true parameter initialization...\n")
  result_true_init <- tryCatch({
    estimate_model_rcpp(model_system, data, init_params = true_params,
                       optimizer = "nlminb", parallel = FALSE, verbose = FALSE)
  }, error = function(e) {
    if (verbose) cat("Error with true init:", e$message, "\n")
    return(list(estimates = rep(NA, length(true_params)),
                std_errors = rep(NA, length(true_params)),
                loglik = NA, convergence = 999))
  })

  # Evaluate likelihood at true parameters
  if (verbose) cat("Evaluating likelihood at true parameters...\n")
  result_true_params <- tryCatch({
    evaluate_likelihood_rcpp(model_system, data, true_params,
                            compute_gradient = FALSE, compute_hessian = FALSE)
  }, error = function(e) {
    if (verbose) cat("Error evaluating at true params:", e$message, "\n")
    return(list(loglik = NA))
  })

  # Handle case where true_params and default_init include fixed parameters
  # If lengths don't match, they likely include fixed params
  n_estimates <- length(result_default$estimates)

  if (length(true_params) > n_estimates) {
    if (verbose) {
      cat(sprintf("Note: true_params has %d values but only %d are estimated (some fixed)\n",
                 length(true_params), n_estimates))
    }
    # For now, just use NA for true values - tests will need to provide correct subset
    true_values_for_comparison <- rep(NA, n_estimates)
  } else {
    true_values_for_comparison <- true_params
  }

  if (length(default_init) > n_estimates) {
    if (verbose) {
      cat(sprintf("Note: default_init has %d values but only %d are estimated (some fixed)\n",
                 length(default_init), n_estimates))
    }
    # For now, just use NA for default init - ideally we'd extract the free params
    default_init_for_comparison <- rep(NA, n_estimates)
  } else {
    default_init_for_comparison <- default_init
  }

  # Adjust param_names length if needed
  if (length(param_names) != n_estimates) {
    if (verbose) {
      cat(sprintf("Warning: param_names has %d entries but %d parameters were estimated\n",
                 length(param_names), n_estimates))
    }
    # Pad or truncate param_names to match
    if (length(param_names) < n_estimates) {
      param_names <- c(param_names, paste0("param_", (length(param_names)+1):n_estimates))
    } else {
      param_names <- param_names[1:n_estimates]
    }
  }

  # Create comparison table
  comparison <- data.frame(
    parameter = param_names,
    true_value = true_values_for_comparison,
    default_init = default_init_for_comparison,
    est_default = result_default$estimates,
    se_default = result_default$std_errors,
    est_true_init = result_true_init$estimates,
    se_true_init = result_true_init$std_errors
  )

  # Add z-statistics (comparing estimates to true values)
  comparison$z_default <- (comparison$est_default - comparison$true_value) / comparison$se_default
  comparison$z_true_init <- (comparison$est_true_init - comparison$true_value) / comparison$se_true_init

  # Likelihood comparison
  loglik_comparison <- data.frame(
    source = c("True parameters", "Est (default init)", "Est (true init)"),
    loglik = c(result_true_params$loglik, result_default$loglik, result_true_init$loglik),
    convergence = c(NA, result_default$convergence, result_true_init$convergence)
  )

  if (verbose) {
    cat("\n--- Parameter Estimates ---\n")
    print(comparison, row.names = FALSE, digits = 4)

    cat("\n--- Likelihood Comparison ---\n")
    print(loglik_comparison, row.names = FALSE, digits = 6)
    cat("\n")
  }

  # Check if estimation was successful
  default_converged <- !is.na(result_default$loglik) && result_default$convergence == 0
  true_init_converged <- !is.na(result_true_init$loglik) && result_true_init$convergence == 0

  # Check if estimates are reasonable (within 3 SEs of true value for most params)
  if (default_converged) {
    reasonable_default <- mean(abs(comparison$z_default) < 3, na.rm = TRUE) > 0.7
  } else {
    reasonable_default <- FALSE
  }

  if (true_init_converged) {
    reasonable_true_init <- mean(abs(comparison$z_true_init) < 3, na.rm = TRUE) > 0.7
  } else {
    reasonable_true_init <- FALSE
  }

  return(list(
    comparison = comparison,
    loglik_comparison = loglik_comparison,
    default_converged = default_converged,
    true_init_converged = true_init_converged,
    reasonable_default = reasonable_default,
    reasonable_true_init = reasonable_true_init
  ))
}

#' Save diagnostics to log file
#'
#' @param test_name Name of the test
#' @param diagnostics List of diagnostic information
#' @param log_dir Directory to save log files
save_diagnostics_to_log <- function(test_name, diagnostics, log_dir = "test_logs") {
  # Create log directory if it doesn't exist
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
  }

  # Create log file name with timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  log_file <- file.path(log_dir, sprintf("%s_%s.txt", test_name, timestamp))

  # Open connection
  con <- file(log_file, open = "w")

  tryCatch({
    # Write header
    writeLines(sprintf("=== Test: %s ===", test_name), con)
    writeLines(sprintf("Timestamp: %s", Sys.time()), con)
    writeLines("", con)

    # Write gradient check results
    if (!is.null(diagnostics$gradient_check)) {
      writeLines("--- Gradient Check ---", con)
      writeLines(sprintf("Pass: %s", diagnostics$gradient_check$pass), con)
      writeLines(sprintf("Max error: %.2e", diagnostics$gradient_check$max_error), con)
      writeLines("", con)
      write.table(diagnostics$gradient_check$diagnostics, con,
                  row.names = FALSE, quote = FALSE)
      writeLines("", con)
    }

    # Write Hessian check results
    if (!is.null(diagnostics$hessian_check)) {
      writeLines("--- Hessian Check ---", con)
      writeLines(sprintf("Pass: %s", diagnostics$hessian_check$pass), con)
      writeLines(sprintf("Max error: %.2e", diagnostics$hessian_check$max_error), con)
      writeLines("", con)
      writeLines("Summary of relative errors:", con)
      capture.output(summary(diagnostics$hessian_check$diagnostics$rel_error),
                    file = con)
      writeLines("", con)
    }

    # Write estimation comparison
    if (!is.null(diagnostics$estimation)) {
      writeLines("--- Estimation Comparison ---", con)
      writeLines(sprintf("Default init converged: %s",
                        diagnostics$estimation$default_converged), con)
      writeLines(sprintf("True init converged: %s",
                        diagnostics$estimation$true_init_converged), con)
      writeLines(sprintf("Default estimates reasonable: %s",
                        diagnostics$estimation$reasonable_default), con)
      writeLines(sprintf("True init estimates reasonable: %s",
                        diagnostics$estimation$reasonable_true_init), con)
      writeLines("", con)
      writeLines("Parameter estimates:", con)
      write.table(diagnostics$estimation$comparison, con,
                  row.names = FALSE, quote = FALSE)
      writeLines("", con)
      writeLines("Likelihood comparison:", con)
      write.table(diagnostics$estimation$loglik_comparison, con,
                  row.names = FALSE, quote = FALSE)
      writeLines("", con)
    }

    # Write overall result
    writeLines("--- Overall Result ---", con)
    writeLines(sprintf("Test passed: %s", diagnostics$overall_pass), con)

  }, finally = {
    close(con)
  })

  return(log_file)
}
