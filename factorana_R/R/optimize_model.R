# ---- Internal helper functions ----

# Build parameter metadata from model system
build_parameter_metadata <- function(model_system) {
  n_factors <- model_system$factor$n_factors

  # Initialize vectors
  param_names <- character(0)
  param_types <- character(0)  # "factor_var", "intercept", "beta", "loading", "sigma", "cutpoint"
  component_id <- integer(0)

  # Add factor variances
  for (k in seq_len(n_factors)) {
    param_names <- c(param_names, sprintf("factor_var_%d", k))
    param_types <- c(param_types, "factor_var")
    component_id <- c(component_id, 0)  # 0 = factor model
  }

  # Add component parameters
  for (i in seq_along(model_system$components)) {
    comp <- model_system$components[[i]]
    comp_name <- comp$name

    # Special handling for multinomial logit (mlogit) with >2 choices
    # Parameters are organized per choice (except reference choice 0)
    if (comp$model_type == "logit" && !is.null(comp$num_choices) && comp$num_choices > 2) {
      n_alternatives <- comp$num_choices - 1  # Exclude reference category

      for (alt in seq_len(n_alternatives)) {
        # Intercept for this alternative
        if (comp$intercept) {
          param_names <- c(param_names, sprintf("%s_intercept_alt%d", comp_name, alt))
          param_types <- c(param_types, "intercept")
          component_id <- c(component_id, i)
        }

        # Covariate coefficients for this alternative
        if (!is.null(comp$covariates) && length(comp$covariates) > 0) {
          for (cov in comp$covariates) {
            if (cov != "intercept") {
              param_names <- c(param_names, sprintf("%s_beta_%s_alt%d", comp_name, cov, alt))
              param_types <- c(param_types, "beta")
              component_id <- c(component_id, i)
            }
          }
        }

        # Factor loadings for this alternative
        if (is.null(comp$loading_normalization)) {
          for (k in seq_len(n_factors)) {
            param_names <- c(param_names, sprintf("%s_loading_%d_alt%d", comp_name, k, alt))
            param_types <- c(param_types, "loading")
            component_id <- c(component_id, i)
          }
        } else {
          for (k in seq_len(n_factors)) {
            if (is.na(comp$loading_normalization[k]) ||
                abs(comp$loading_normalization[k]) < 1e-10) {
              param_names <- c(param_names, sprintf("%s_loading_%d_alt%d", comp_name, k, alt))
              param_types <- c(param_types, "loading")
              component_id <- c(component_id, i)
            }
          }
        }
      }
    } else {
      # Standard handling for all other model types
      # Intercept
      if (comp$intercept) {
        if (comp$model_type == "oprobit" && comp$num_choices > 2) {
          n_intercepts <- comp$num_choices - 1
          for (j in seq_len(n_intercepts)) {
            param_names <- c(param_names, sprintf("%s_intercept_%d", comp_name, j))
            param_types <- c(param_types, "intercept")
            component_id <- c(component_id, i)
          }
        } else {
          param_names <- c(param_names, sprintf("%s_intercept", comp_name))
          param_types <- c(param_types, "intercept")
          component_id <- c(component_id, i)
        }
      }

      # Covariate coefficients
      if (!is.null(comp$covariates) && length(comp$covariates) > 0) {
        for (cov in comp$covariates) {
          if (cov != "intercept") {
            param_names <- c(param_names, sprintf("%s_beta_%s", comp_name, cov))
            param_types <- c(param_types, "beta")
            component_id <- c(component_id, i)
          }
        }
      }

      # Factor loadings
      if (is.null(comp$loading_normalization)) {
        for (k in seq_len(n_factors)) {
          param_names <- c(param_names, sprintf("%s_loading_%d", comp_name, k))
          param_types <- c(param_types, "loading")
          component_id <- c(component_id, i)
        }
      } else {
        for (k in seq_len(n_factors)) {
          if (is.na(comp$loading_normalization[k])) {
            param_names <- c(param_names, sprintf("%s_loading_%d", comp_name, k))
            param_types <- c(param_types, "loading")
            component_id <- c(component_id, i)
          }
        }
      }
    }

    # Residual variance
    if (comp$model_type %in% c("linear", "oprobit")) {
      param_names <- c(param_names, sprintf("%s_sigma", comp_name))
      param_types <- c(param_types, "sigma")
      component_id <- c(component_id, i)
    }

    # Cutpoints
    if (comp$model_type == "oprobit" && !is.null(comp$n_cutpoints)) {
      for (j in seq_len(comp$n_cutpoints)) {
        param_names <- c(param_names, sprintf("%s_cutpoint_%d", comp_name, j))
        param_types <- c(param_types, "cutpoint")
        component_id <- c(component_id, i)
      }
    }
  }

  return(list(
    names = param_names,
    types = param_types,
    component_id = component_id,
    n_params = length(param_names)
  ))
}

# Setup parameter bounds and identify fixed parameters
setup_parameter_constraints <- function(model_system, init_params, param_metadata) {
  n_params <- length(init_params)
  n_factors <- model_system$factor$n_factors

  lower_bounds <- rep(-Inf, n_params)
  upper_bounds <- rep(Inf, n_params)

  # Identify which factor variances are identified
  # A factor is identified if at least one component has a fixed non-zero loading on that factor
  factor_variance_identified <- rep(FALSE, n_factors)
  for (comp in model_system$components) {
    if (!is.null(comp$loading_normalization)) {
      for (k in seq_len(n_factors)) {
        if (!is.na(comp$loading_normalization[k]) &&
            abs(comp$loading_normalization[k]) > 1e-6) {
          factor_variance_identified[k] <- TRUE
        }
      }
    }
  }

  param_fixed <- rep(FALSE, n_params)

  # Handle previous_stage: mark all previous-stage parameters as fixed
  if (!is.null(model_system$previous_stage_info)) {
    n_fixed <- model_system$previous_stage_info$n_params_fixed
    if (n_fixed > 0 && n_fixed <= n_params) {
      # Fix first n_fixed parameters (from previous stage)
      param_fixed[1:n_fixed] <- TRUE
      # Set bounds to fixed values
      fixed_values <- model_system$previous_stage_info$fixed_param_values
      lower_bounds[1:n_fixed] <- fixed_values
      upper_bounds[1:n_fixed] <- fixed_values
    }
  }

  # Track cutpoint indices per component to identify incremental thresholds
  cutpoint_counter <- list()

  for (i in seq_len(n_params)) {
    param_type <- param_metadata$types[i]
    comp_id <- param_metadata$component_id[i]

    # Fix non-identified factor variances
    if (param_type == "factor_var") {
      factor_idx <- i
      if (!factor_variance_identified[factor_idx]) {
        param_fixed[i] <- TRUE
        lower_bounds[i] <- init_params[i]
        upper_bounds[i] <- init_params[i]
      }
    }

    # Set lower bound for sigma parameters
    if (param_type == "sigma") {
      lower_bounds[i] <- 0.01
    }

    # Set lower bounds for ordered probit cutpoints (incremental parameterization)
    # The first cutpoint is unrestricted (can be any value)
    # Subsequent cutpoints are increments and must be positive to ensure ordering
    if (param_type == "cutpoint") {
      comp_key <- as.character(comp_id)
      if (is.null(cutpoint_counter[[comp_key]])) {
        cutpoint_counter[[comp_key]] <- 1
        # First cutpoint is unrestricted
      } else {
        cutpoint_counter[[comp_key]] <- cutpoint_counter[[comp_key]] + 1
        # Subsequent cutpoints are increments, must be positive
        lower_bounds[i] <- 0.01
      }
    }
  }

  free_idx <- which(!param_fixed)
  n_free <- length(free_idx)

  return(list(
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    param_fixed = param_fixed,
    free_idx = free_idx,
    n_free = n_free,
    lower_bounds_free = lower_bounds[free_idx],
    upper_bounds_free = upper_bounds[free_idx]
  ))
}

# ---- Main estimation function ----

#' Estimate factor model using R-based optimization
#'
#' This function estimates a factor model by optimizing the likelihood using
#' C++ for fast evaluation and R for optimization and parallelization.
#'
#' @param model_system A model_system object from define_model_system()
#' @param data Data frame containing all variables
#' @param init_params Initial parameter values (optional)
#' @param control Estimation control object from define_estimation_control()
#' @param optimizer Optimizer to use (default "nloptr"):
#'   \itemize{
#'     \item \code{"nloptr"} - L-BFGS via nloptr (gradient only, no Hessian)
#'     \item \code{"optim"} - L-BFGS-B via stats::optim (gradient only, no Hessian)
#'     \item \code{"nlminb"} - Uses analytical gradient AND Hessian (more efficient!)
#'     \item \code{"trust"} - Trust region method with analytical Hessian (requires trustOptim package)
#'   }
#' @param parallel Whether to use parallel computation (default TRUE)
#' @param verbose Whether to print progress (default TRUE)
#'
#' @details
#' For maximum efficiency, use \code{optimizer = "nlminb"} or \code{optimizer = "trust"}
#' which exploit the analytical Hessian computed in C++. The default L-BFGS methods
#' only use the gradient and approximate the Hessian from gradient history.
#'
#' @return List with parameter estimates, standard errors, log-likelihood, etc.
#' @export
estimate_model_rcpp <- function(model_system, data, init_params = NULL,
                                control = NULL, optimizer = "nlminb",
                                parallel = TRUE, verbose = TRUE) {

  # WORKAROUND: Deep copy model_system to avoid C++ reuse bug
  # Use serialize/unserialize for true deep copy
  model_system <- unserialize(serialize(model_system, NULL))

  # Default control if not provided
  if (is.null(control)) {
    control <- define_estimation_control()
  }

  # Convert data to matrix
  data_mat <- as.matrix(data)

  # Setup parallel cluster if requested
  cl <- NULL
  if (parallel && control$num_cores > 1) {
    # Split data across workers first to determine actual number of workers needed
    n_obs <- nrow(data_mat)
    n_per_worker <- ceiling(n_obs / control$num_cores)
    data_splits <- split(1:n_obs, ceiling(1:n_obs / n_per_worker))

    # Use actual number of data splits (may be less than requested cores)
    n_workers <- length(data_splits)

    if (verbose) {
      message(sprintf("Setting up parallel cluster with %d cores (requested %d)...",
                     n_workers, control$num_cores))
    }

    # Use doParallel for Windows compatibility
    cl <- parallel::makeCluster(n_workers)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)

  } else {
    # Single worker
    data_splits <- list(1:nrow(data_mat))
  }

  # Initialize factor models on each worker
  if (verbose) message("Initializing C++ likelihood evaluators...")

  # Get n_quad from control
  n_quad <- control$n_quad_points

  if (!is.null(cl)) {
    # Export necessary objects to workers
    parallel::clusterExport(cl, c("model_system", "data_mat", "n_quad", "data_splits"),
                           envir = environment())
    parallel::clusterEvalQ(cl, {
      library(factorana)
    })

    # Set worker IDs (1 to n_workers) in each worker's global environment
    for (i in seq_along(cl)) {
      parallel::clusterCall(cl[i], assign, ".self_id", i, envir = .GlobalEnv)
    }

    # Initialize FactorModel on each worker with its data subset
    # IMPORTANT: Store pointer in worker's global env to avoid serialization
    parallel::clusterEvalQ(cl, {
      worker_id <- .self_id
      idx <- data_splits[[worker_id]]
      data_subset <- data_mat[idx, , drop = FALSE]
      .fm_ptr <- initialize_factor_model_cpp(model_system, data_subset, n_quad)
    })

    fm_ptrs <- NULL  # Not used; pointers stored on workers
  } else {
    # Single worker initialization
    fm_ptrs <- list(initialize_factor_model_cpp(model_system, data_mat, n_quad))
  }

  # Get initial parameters if not provided
  if (is.null(init_params)) {
    # Use smart initialization based on separate component estimation
    init_result <- initialize_parameters(model_system, data, verbose = verbose)
    init_params <- init_result$init_params

    # Get parameter count from C++ object
    if (!is.null(cl)) {
      # For parallel case, get from first worker
      param_info <- parallel::clusterEvalQ(cl, get_parameter_info_cpp(.fm_ptr))[[1]]
    } else {
      param_info <- get_parameter_info_cpp(fm_ptrs[[1]])
    }
    n_params <- param_info$n_param_free

    if (length(init_params) != n_params) {
      stop(sprintf("Parameter initialization returned %d parameters but C++ expects %d",
                   length(init_params), n_params))
    }
  }

  # Build parameter metadata (names, types, etc.)
  param_metadata <- build_parameter_metadata(model_system)

  # Setup constraints and identify fixed/free parameters
  param_constraints <- setup_parameter_constraints(model_system, init_params, param_metadata)

  if (verbose) {
    message(sprintf("Total parameters: %d", length(init_params)))
    message(sprintf("Free parameters: %d", param_constraints$n_free))
    if (any(param_constraints$param_fixed)) {
      fixed_names <- param_metadata$names[param_constraints$param_fixed]
      message(sprintf("Fixed parameters (%d): %s",
                     sum(param_constraints$param_fixed),
                     paste(fixed_names, collapse = ", ")))
    }
    n_sigma <- sum(param_metadata$types == "sigma")
    if (n_sigma > 0) {
      message(sprintf("Sigma parameters (%d) have lower bound = 0.01", n_sigma))
    }
  }

  # Define objective function (operates on free parameters only)
  objective_fn <- function(params_free) {
    # Reconstruct full parameter vector
    params_full <- init_params
    params_full[param_constraints$free_idx] <- params_free

    if (!is.null(cl)) {
      # Parallel evaluation: aggregate across workers
      parallel::clusterExport(cl, "params_full", envir = environment())
      loglik_parts <- parallel::clusterEvalQ(cl, {
        evaluate_loglik_only_cpp(.fm_ptr, params_full)
      })
      loglik <- sum(unlist(loglik_parts))
    } else {
      # Single worker
      loglik <- evaluate_loglik_only_cpp(fm_ptrs[[1]], params_full)
    }
    return(-loglik)  # Negative for minimization
  }

  # Define gradient function (returns gradient for free parameters only)
  gradient_fn <- function(params_free) {
    # Reconstruct full parameter vector
    params_full <- init_params
    params_full[param_constraints$free_idx] <- params_free

    if (!is.null(cl)) {
      # Parallel evaluation: aggregate gradients
      parallel::clusterExport(cl, "params_full", envir = environment())
      grad_parts <- parallel::clusterEvalQ(cl, {
        result <- evaluate_likelihood_cpp(.fm_ptr, params_full,
                                         compute_gradient = TRUE,
                                         compute_hessian = FALSE)
        result$gradient
      })
      grad_full <- Reduce(`+`, grad_parts)
    } else {
      result <- evaluate_likelihood_cpp(fm_ptrs[[1]], params_full,
                                       compute_gradient = TRUE,
                                       compute_hessian = FALSE)
      grad_full <- result$gradient
    }

    # Return gradient for free parameters only
    return(-grad_full[param_constraints$free_idx])  # Negative for minimization
  }

  # Define Hessian function (returns Hessian for free parameters only)
  hessian_fn <- function(params_free) {
    # Reconstruct full parameter vector
    params_full <- init_params
    params_full[param_constraints$free_idx] <- params_free

    if (!is.null(cl)) {
      # Parallel evaluation: aggregate Hessians
      parallel::clusterExport(cl, "params_full", envir = environment())
      hess_parts <- parallel::clusterEvalQ(cl, {
        result <- evaluate_likelihood_cpp(.fm_ptr, params_full,
                                         compute_gradient = FALSE,
                                         compute_hessian = TRUE)
        result$hessian
      })
      # Aggregate Hessian matrices (they're upper triangles, addition is correct)
      hess <- Reduce(`+`, hess_parts)
    } else {
      result <- evaluate_likelihood_cpp(fm_ptrs[[1]], params_full,
                                       compute_gradient = FALSE,
                                       compute_hessian = TRUE)
      hess <- result$hessian
    }

    # Convert upper triangle vector to full symmetric matrix
    n_params_full <- length(params_full)
    hess_mat_full <- matrix(0, n_params_full, n_params_full)
    idx <- 1
    for (i in 1:n_params_full) {
      for (j in i:n_params_full) {
        hess_mat_full[i, j] <- hess[idx]
        hess_mat_full[j, i] <- hess[idx]  # Symmetrize
        idx <- idx + 1
      }
    }

    # Extract submatrix for free parameters only
    hess_mat_free <- hess_mat_full[param_constraints$free_idx, param_constraints$free_idx, drop = FALSE]

    return(-hess_mat_free)  # Negative for minimization
  }

  # Run optimization
  if (verbose) message(sprintf("Running optimization using %s...", optimizer))

  if (optimizer == "nloptr") {
    if (!requireNamespace("nloptr", quietly = TRUE)) {
      stop("nloptr package required but not installed")
    }

    result <- nloptr::nloptr(
      x0 = init_params[param_constraints$free_idx],
      eval_f = objective_fn,
      eval_grad_f = gradient_fn,
      lb = param_constraints$lower_bounds_free,
      ub = param_constraints$upper_bounds_free,
      opts = list(
        algorithm = "NLOPT_LD_LBFGS",  # L-BFGS with analytical gradient (no Hessian)
        maxeval = 1000,
        xtol_rel = 1e-6,
        print_level = if (verbose) 2 else 0
      )
    )

    # Reconstruct full parameter vector
    estimates <- init_params
    estimates[param_constraints$free_idx] <- result$solution
    loglik <- -result$objective
    convergence <- result$status

  } else if (optimizer == "optim") {
    result <- stats::optim(
      par = init_params[param_constraints$free_idx],
      fn = objective_fn,
      gr = gradient_fn,
      method = "L-BFGS-B",  # Quasi-Newton (no Hessian)
      lower = param_constraints$lower_bounds_free,
      upper = param_constraints$upper_bounds_free,
      control = list(trace = if (verbose) 1 else 0)
    )

    # Reconstruct full parameter vector
    estimates <- init_params
    estimates[param_constraints$free_idx] <- result$par
    loglik <- -result$value
    convergence <- result$convergence

  } else if (optimizer == "nlminb") {
    # nlminb can use analytical Hessian for more efficient optimization
    if (verbose) message("  Using analytical gradient and Hessian")

    result <- stats::nlminb(
      start = init_params[param_constraints$free_idx],
      objective = objective_fn,
      gradient = gradient_fn,
      hessian = hessian_fn,  # Uses analytical Hessian!
      lower = param_constraints$lower_bounds_free,
      upper = param_constraints$upper_bounds_free,
      control = list(
        trace = if (verbose) 1 else 0,
        eval.max = 5000,
        iter.max = 2000
      )
    )

    # Reconstruct full parameter vector
    estimates <- init_params
    estimates[param_constraints$free_idx] <- result$par
    loglik <- -result$objective
    convergence <- result$convergence

  } else if (optimizer == "trust") {
    # Trust region method - most efficient with analytical Hessian
    if (!requireNamespace("trustOptim", quietly = TRUE)) {
      stop("trustOptim package required but not installed. Install with: install.packages('trustOptim')")
    }

    # Warning: trust optimizer doesn't support parameter bounds
    if (any(param_constraints$param_fixed)) {
      warning("trust optimizer does not support fixed parameters. ",
              "Fixed parameters may drift from their initial values. ",
              "Consider using 'nlminb' or 'optim' instead.")
    }

    if (verbose) message("  Using trust region with analytical gradient and Hessian")

    result <- trustOptim::trust.optim(
      x = init_params[param_constraints$free_idx],
      fn = objective_fn,
      gr = gradient_fn,
      hs = hessian_fn,  # Uses analytical Hessian!
      method = "SR1",   # Symmetric rank-1 trust region
      control = list(
        report.level = if (verbose) 2 else 0,
        maxit = 500
      )
    )

    # Reconstruct full parameter vector
    estimates <- init_params
    estimates[param_constraints$free_idx] <- result$solution
    loglik <- -result$value
    convergence <- result$converged

  } else {
    stop("Unknown optimizer: ", optimizer, "\n",
         "Available: 'nloptr' (L-BFGS), 'optim' (L-BFGS-B), 'nlminb' (uses Hessian), 'trust' (uses Hessian)")
  }

  # Compute Hessian for standard errors
  if (verbose) message("Computing Hessian for standard errors...")

  if (!is.null(cl)) {
    parallel::clusterExport(cl, "estimates", envir = environment())
    hess_parts <- parallel::clusterEvalQ(cl, {
      result <- evaluate_likelihood_cpp(.fm_ptr, estimates,
                                       compute_gradient = TRUE,
                                       compute_hessian = TRUE)
      result$hessian
    })
    # Aggregate Hessian matrices (they're upper triangles, addition is correct)
    hessian <- Reduce(`+`, hess_parts)
  } else {
    result <- evaluate_likelihood_cpp(fm_ptrs[[1]], estimates,
                                     compute_gradient = TRUE,
                                     compute_hessian = TRUE)
    hessian <- result$hessian
  }

  # IMPORTANT: The C++ code returns the Hessian of the log-likelihood.
  # Since we're minimizing the NEGATIVE log-likelihood, we need to negate
  # the Hessian to get the correct second derivatives.
  hessian <- -hessian

  # Compute standard errors from Hessian
  # The Hessian is the matrix of second derivatives of the negative log-likelihood
  # Standard errors are sqrt(diag(inv(Hessian)))
  std_errors <- rep(NA, length(estimates))

  tryCatch({
    # Convert upper triangle to full symmetric matrix
    n_params <- length(estimates)
    hess_matrix <- matrix(0, n_params, n_params)

    # Fill in the upper triangle
    idx <- 1
    for (i in 1:n_params) {
      for (j in i:n_params) {
        hess_matrix[i, j] <- hessian[idx]
        if (i != j) {
          hess_matrix[j, i] <- hessian[idx]  # Symmetric
        }
        idx <- idx + 1
      }
    }

    # Verify the matrix is symmetric
    if (verbose) {
      max_asym <- max(abs(hess_matrix - t(hess_matrix)))
      if (max_asym > 1e-10) {
        warning(sprintf("Hessian matrix is not symmetric (max diff: %.2e)", max_asym))
      }
    }

    # Use parameter constraints to identify fixed vs free parameters
    fixed_params <- which(param_constraints$param_fixed)
    free_params <- param_constraints$free_idx

    if (length(free_params) > 0) {
      # Extract Hessian for free parameters only
      hess_free <- hess_matrix[free_params, free_params, drop = FALSE]

      # Check condition number of Hessian
      if (verbose) {
        eig_vals <- eigen(hess_free, only.values = TRUE)$values
        cond_num <- max(abs(eig_vals)) / min(abs(eig_vals))
        message(sprintf("  Hessian condition number: %.2e", cond_num))
        message(sprintf("  Eigenvalues: min=%.2e, max=%.2e", min(abs(eig_vals)), max(abs(eig_vals))))
        if (cond_num > 1e10) {
          message("  Warning: Hessian is poorly conditioned")
        }

        # Show Hessian diagonal elements
        message("\n  Hessian diagonal elements (2nd derivatives of -loglik):")
        hess_diag <- diag(hess_free)
        for (i in seq_along(free_params)) {
          param_idx <- free_params[i]
          message(sprintf("    [%2d] %.4e", param_idx, hess_diag[i]))
        }
      }

      # Invert the Hessian to get the covariance matrix for free parameters
      # Use SVD-based pseudoinverse for numerical stability
      if (verbose) {
        message("\n  Using SVD-based pseudoinverse for numerical stability")
      }

      svd_result <- svd(hess_free)

      # Set tolerance for singular values
      tol <- max(dim(hess_free)) * max(svd_result$d) * .Machine$double.eps * 100  # More conservative
      pos_idx <- svd_result$d > tol

      if (verbose) {
        message(sprintf("  SVD: %d/%d singular values kept (tol=%.2e)", sum(pos_idx), length(pos_idx), tol))
        message(sprintf("  Singular values range: [%.2e, %.2e]", min(svd_result$d), max(svd_result$d)))
      }

      # Compute pseudoinverse
      d_inv <- rep(0, length(svd_result$d))
      d_inv[pos_idx] <- 1 / svd_result$d[pos_idx]

      cov_free <- svd_result$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(svd_result$u)

      # Check if the covariance matrix is positive definite
      if (verbose) {
        cov_eig <- eigen(cov_free, symmetric = TRUE, only.values = TRUE)$values
        n_neg_eig <- sum(cov_eig < -1e-10)  # Allow small numerical errors
        if (n_neg_eig > 0) {
          message(sprintf("  Warning: Covariance matrix has %d negative eigenvalues", n_neg_eig))
          message(sprintf("  Eigenvalues range: [%.2e, %.2e]", min(cov_eig), max(cov_eig)))
        } else {
          message(sprintf("  Covariance matrix is positive semi-definite"))
        }
      }

      # Standard errors for free parameters
      cov_diag <- diag(cov_free)
      se_free <- sqrt(pmax(0, cov_diag))  # pmax ensures non-negative
      std_errors[free_params] <- se_free

      # Fixed parameters have zero standard error
      # Exception: previous_stage parameters should retain their SEs
      std_errors[fixed_params] <- 0.0

      # If we have previous_stage, use those standard errors
      if (!is.null(model_system$previous_stage_info)) {
        n_prev_stage <- model_system$previous_stage_info$n_params_fixed
        prev_se <- model_system$previous_stage_info$fixed_std_errors
        if (length(prev_se) == n_prev_stage) {
          std_errors[1:n_prev_stage] <- prev_se
        }
      }

      # Diagnostic: show covariance diagonal and standard errors
      if (verbose) {
        message("\n  Covariance matrix diagonal (variance estimates):")
        for (i in seq_along(free_params)) {
          param_idx <- free_params[i]
          message(sprintf("    [%2d] Var=%.4e, SE=%.4e", param_idx, cov_diag[i], se_free[i]))
        }

        n_small_se <- sum(se_free < 1e-6 & se_free > 0)
        n_zero_se <- sum(se_free == 0)
        if (n_small_se > 0) {
          message(sprintf("\n  Warning: %d parameters have very small SE (< 1e-6)", n_small_se))
        }
        if (n_zero_se > 0) {
          message(sprintf("  Warning: %d parameters have zero SE", n_zero_se))
        }
      }
    }

    if (verbose) {
      n_free <- length(free_params)
      n_fixed <- length(fixed_params)
      message(sprintf("Standard errors computed (%d free, %d fixed parameters)", n_free, n_fixed))
    }
  }, error = function(e) {
    if (verbose) {
      warning("Could not compute standard errors: ", e$message)
    }
  })

  if (verbose) {
    # Report convergence status based on optimizer
    if ((optimizer %in% c("optim", "nlminb") && convergence == 0) ||
        (optimizer %in% c("L-BFGS", "L-BFGS-B") && convergence == 0) ||
        (optimizer == "trust" && convergence == TRUE)) {
      message(sprintf("Converged successfully. Log-likelihood: %.4f", loglik))
    } else {
      # Convergence failed - determine reason
      if (optimizer == "nlminb") {
        reason <- switch(as.character(convergence),
                        "1" = "iteration limit reached",
                        "non-zero convergence code")
      } else if (optimizer %in% c("optim", "L-BFGS", "L-BFGS-B")) {
        reason <- switch(as.character(convergence),
                        "1" = "iteration limit reached",
                        "10" = "degeneracy in Nelder-Mead simplex",
                        "non-zero convergence code")
      } else {
        reason <- "failed to converge"
      }
      message(sprintf("Convergence FAILED (%s). Log-likelihood: %.4f", reason, loglik))
    }
  }

  # Return results
  list(
    estimates = estimates,
    std_errors = std_errors,
    loglik = loglik,
    convergence = convergence,
    model_system = model_system,
    optimizer = optimizer
  )
}


#' Helper function to convert model_system to format expected by C++
#'
#' @param model_system model_system object
#' @param data Data frame
#' @return List in format expected by initialize_factor_model_cpp
#' @keywords internal
prepare_model_system_for_cpp <- function(model_system, data) {
  # Extract variable names and create index mapping
  var_names <- names(data)
  var_index <- setNames(seq_along(var_names) - 1, var_names)  # 0-indexed

  # Prepare factor model info
  fm <- model_system$factor
  factor_info <- list(
    n_factors = fm$n_factors,
    n_types = fm$n_types,
    n_mixtures = fm$n_mixtures,
    correlation = fm$correlation,
    loading_normalization = fm$loading_normalization
  )

  # Prepare component info
  components_info <- lapply(model_system$components, function(comp) {
    # Get outcome index
    outcome_idx <- var_index[comp$outcome]

    # Get regressor indices
    regressor_idx <- var_index[comp$covariates]

    # Get missing indicator if present
    missing_idx <- if (!is.null(comp$evaluation_indicator)) {
      var_index[comp$evaluation_indicator]
    } else {
      -1L
    }

    list(
      name = comp$name,
      model_type = comp$model_type,
      outcome_idx = outcome_idx,
      missing_idx = missing_idx,
      regressor_idx = as.integer(regressor_idx),
      factor_normalizations = comp$loading_normalization,
      num_choices = if (comp$model_type == "oprobit") {
        length(unique(data[[comp$outcome]]))
      } else {
        2L
      }
    )
  })

  list(
    factor = factor_info,
    components = components_info
  )
}
