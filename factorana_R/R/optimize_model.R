# ---- Internal helper functions ----

# Build parameter metadata from model system
build_parameter_metadata <- function(model_system) {
  n_factors <- model_system$factor$n_factors

  # Initialize vectors
  param_names <- character(0)
  param_types <- character(0)  # "factor_var", "factor_corr", "intercept", "beta", "loading", "sigma", "cutpoint"
  component_id <- integer(0)

  # Add factor variances
  for (k in seq_len(n_factors)) {
    param_names <- c(param_names, sprintf("factor_var_%d", k))
    param_types <- c(param_types, "factor_var")
    component_id <- c(component_id, 0)  # 0 = factor model
  }

  # Add factor correlation if correlation = TRUE and n_factors = 2
  if (isTRUE(model_system$factor$correlation) && n_factors == 2) {
    param_names <- c(param_names, "factor_corr_1_2")
    param_types <- c(param_types, "factor_corr")
    component_id <- c(component_id, 0)  # 0 = factor model
  }

  # Add type model loading parameters if n_types > 1
  # Type model: log(P(type=t)/P(type=1)) = sum_k lambda_t_k * f_k
  # (n_types - 1) * n_factors parameters (type 1 is reference)
  n_types <- model_system$factor$n_types
  if (!is.null(n_types) && n_types > 1L) {
    for (t in 2:n_types) {
      for (k in seq_len(n_factors)) {
        param_names <- c(param_names, sprintf("type_%d_loading_%d", t, k))
        param_types <- c(param_types, "type_loading")
        component_id <- c(component_id, 0)  # 0 = factor model
      }
    }
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

        # Quadratic factor loadings for this alternative (if factor_spec includes quadratic)
        if (!is.null(comp$factor_spec) && comp$factor_spec %in% c("quadratic", "full")) {
          for (k in seq_len(n_factors)) {
            param_names <- c(param_names, sprintf("%s_loading_quad_%d_alt%d", comp_name, k, alt))
            param_types <- c(param_types, "loading_quad")
            component_id <- c(component_id, i)
          }
        }

        # Interaction factor loadings for this alternative (if factor_spec includes interactions)
        if (!is.null(comp$factor_spec) && comp$factor_spec %in% c("interactions", "full") && n_factors >= 2) {
          for (j in seq_len(n_factors - 1)) {
            for (kk in (j + 1):n_factors) {
              param_names <- c(param_names, sprintf("%s_loading_inter_%d_%d_alt%d", comp_name, j, kk, alt))
              param_types <- c(param_types, "loading_inter")
              component_id <- c(component_id, i)
            }
          }
        }
      }
    } else {
      # Standard handling for all other model types
      # Covariate coefficients - naming must match initialize_parameters.R: comp_name_covariate
      if (!is.null(comp$covariates) && length(comp$covariates) > 0) {
        for (cov in comp$covariates) {
          param_names <- c(param_names, sprintf("%s_%s", comp_name, cov))
          # Mark intercept type separately for constraint handling
          if (cov == "intercept") {
            param_types <- c(param_types, "intercept")
          } else {
            param_types <- c(param_types, "beta")
          }
          component_id <- c(component_id, i)
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

      # Quadratic factor loadings (if factor_spec includes quadratic)
      if (!is.null(comp$factor_spec) && comp$factor_spec %in% c("quadratic", "full")) {
        for (k in seq_len(n_factors)) {
          param_names <- c(param_names, sprintf("%s_loading_quad_%d", comp_name, k))
          param_types <- c(param_types, "loading_quad")
          component_id <- c(component_id, i)
        }
      }

      # Interaction factor loadings (if factor_spec includes interactions)
      if (!is.null(comp$factor_spec) && comp$factor_spec %in% c("interactions", "full") && n_factors >= 2) {
        for (j in seq_len(n_factors - 1)) {
          for (kk in (j + 1):n_factors) {
            param_names <- c(param_names, sprintf("%s_loading_inter_%d_%d", comp_name, j, kk))
            param_types <- c(param_types, "loading_inter")
            component_id <- c(component_id, i)
          }
        }
      }
    }

    # Residual variance (only for linear models - oprobit has normalized variance)
    if (comp$model_type == "linear") {
      param_names <- c(param_names, sprintf("%s_sigma", comp_name))
      param_types <- c(param_types, "sigma")
      component_id <- c(component_id, i)
    }

    # Thresholds for ordered probit (num_choices - 1 thresholds for identification)
    if (comp$model_type == "oprobit" && !is.null(comp$num_choices) && comp$num_choices > 1) {
      n_thresholds <- comp$num_choices - 1
      for (j in seq_len(n_thresholds)) {
        param_names <- c(param_names, sprintf("%s_thresh_%d", comp_name, j))
        param_types <- c(param_types, "cutpoint")
        component_id <- c(component_id, i)
      }
    }

    # Type-specific intercepts for this component (if n_types > 1)
    # Added after all other component parameters to match initialize_parameters.R ordering
    if (!is.null(n_types) && n_types > 1L) {
      for (t in 2:n_types) {
        param_names <- c(param_names, sprintf("%s_type_%d_intercept", comp_name, t))
        param_types <- c(param_types, "type_intercept")
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
setup_parameter_constraints <- function(model_system, init_params, param_metadata, factor_variance_fixed = NULL, verbose = FALSE) {
  n_params <- length(init_params)
  n_factors <- model_system$factor$n_factors

  lower_bounds <- rep(-Inf, n_params)
  upper_bounds <- rep(Inf, n_params)

  # Identify which factor variances are identified
  # If factor_variance_fixed was passed from initialize_parameters, use it
  # Otherwise, compute it by checking for fixed non-zero loadings
  if (!is.null(factor_variance_fixed)) {
    # Use the passed value: factor_variance_fixed=TRUE means variance is estimated (identified by loading)
    # So factor_variance_identified = factor_variance_fixed
    factor_variance_identified <- factor_variance_fixed
    if (verbose) {
      message(sprintf("Factor variance identification (from initialize_parameters): %s",
                      paste(factor_variance_identified, collapse=", ")))
    }
  } else {
    # Fallback: compute by checking for fixed non-zero loadings
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
    if (verbose) {
      message(sprintf("Factor variance identification (computed fallback): %s",
                      paste(factor_variance_identified, collapse=", ")))
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
      # For factor variances, check identification
      # Only the first n_factors parameters are factor variances
      factor_idx <- i
      if (factor_idx <= length(factor_variance_identified)) {
        if (verbose) {
          message(sprintf("  Constraint check: param %d (factor_var) -> factor_idx=%d, identified=%s",
                          i, factor_idx, factor_variance_identified[factor_idx]))
        }
        if (!factor_variance_identified[factor_idx]) {
          param_fixed[i] <- TRUE
          lower_bounds[i] <- init_params[i]
          upper_bounds[i] <- init_params[i]
          if (verbose) {
            message(sprintf("    -> FIXED at %.6f", init_params[i]))
          }
        } else {
          # Set lower bound for free factor variances to prevent numerical issues
          # (division by sqrt(factor_var) in gradient/Hessian chain rule)
          lower_bounds[i] <- 0.01
        }
      }
    }

    # Set bounds for factor correlation parameters
    if (param_type == "factor_corr") {
      # Correlation must be between -1 and 1
      lower_bounds[i] <- -0.99
      upper_bounds[i] <- 0.99
      if (verbose) {
        message(sprintf("  Constraint: param %d (factor_corr) -> bounds = [-0.99, 0.99]", i))
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

    # Fix oprobit intercepts at 0 (absorbed into thresholds)
    if (param_type == "intercept") {
      comp <- model_system$components[[comp_id]]
      if (!is.null(comp) && comp$model_type == "oprobit") {
        param_fixed[i] <- TRUE
        lower_bounds[i] <- 0.0
        upper_bounds[i] <- 0.0
        if (verbose) {
          message(sprintf("  Fixed oprobit intercept: param %d (%s) at 0",
                          i, param_metadata$names[i]))
        }
      }
    }

    # Fix type intercepts that were marked as fixed via fix_type_intercepts()
    if (param_type == "type_intercept") {
      comp <- model_system$components[[comp_id]]
      if (!is.null(comp) && !is.null(comp$fixed_type_intercepts) && length(comp$fixed_type_intercepts) > 0) {
        # Extract type number from parameter name (e.g., "Y_type_2_intercept" -> 2)
        param_name <- param_metadata$names[i]
        type_match <- regmatches(param_name, regexec("_type_([0-9]+)_intercept$", param_name))[[1]]
        if (length(type_match) >= 2) {
          type_num <- as.integer(type_match[2])
          if (is_type_intercept_fixed(comp, type_num, choice = NULL)) {
            param_fixed[i] <- TRUE
            lower_bounds[i] <- 0.0
            upper_bounds[i] <- 0.0
            if (verbose) {
              message(sprintf("  Fixed type intercept: param %d (%s) at 0",
                              i, param_name))
            }
          }
        }
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

# ---- Helper function for eigenvector-based saddle point escape ----

#' Escape saddle point using eigenvector-based restarts
#'
#' When optimization converges to a saddle point (detected by positive Hessian
#' eigenvalues), this function perturbs parameters in the direction of positive
#' eigenvectors to escape the saddle point.
#'
#' @param params Current parameter estimates at saddle point
#' @param hessian_fn Function to compute Hessian at given parameters
#' @param objective_fn Function to compute objective at given parameters
#' @param param_constraints List with lower/upper bounds and free_idx
#' @param verbose Whether to print progress
#' @return List with new_params and found_escape (TRUE if escape direction found)
#' @keywords internal
find_saddle_escape_direction <- function(params, hessian_fn, objective_fn,
                                          param_constraints, verbose = FALSE) {
  # Compute Hessian at current point
  hess_mat <- hessian_fn(params)

  # Get eigenvalue decomposition
  eig <- eigen(hess_mat, symmetric = TRUE)

  # For minimization: saddle point has at least one NEGATIVE eigenvalue

  # (positive eigenvalues are good - they indicate local minimum directions)
  neg_idx <- which(eig$values < -1e-6)

  if (length(neg_idx) == 0) {
    if (verbose) message("  No negative eigenvalues found - not a saddle point")
    return(list(new_params = params, found_escape = FALSE))
  }

  if (verbose) {
    message(sprintf("  Found %d negative eigenvalue(s): %s",
                   length(neg_idx),
                   paste(round(eig$values[neg_idx], 2), collapse = ", ")))
  }

  # Current objective value
  current_obj <- objective_fn(params)
  best_obj <- current_obj
  best_params <- params

  # Try moving in direction of negative eigenvector(s)
  for (idx in neg_idx) {
    eigvec <- eig$vectors[, idx]

    # Scale step: use a minimum step size to avoid getting stuck with tiny steps
    # Larger negative eigenvalues suggest steeper escape, but ensure at least 0.1
    base_step <- max(0.1, min(1.0, sqrt(abs(eig$values[idx]))))

    # Try both directions with different step sizes
    for (direction in c(-1, 1)) {
      for (step_mult in c(0.5, 1.0, 2.0, 5.0)) {
        step <- base_step * step_mult

        new_params <- params + direction * step * eigvec

        # Apply bounds
        new_params <- pmax(new_params, param_constraints$lower_bounds_free)
        new_params <- pmin(new_params, param_constraints$upper_bounds_free)

        # Evaluate objective at new point
        new_obj <- tryCatch(
          objective_fn(new_params),
          error = function(e) Inf
        )

        if (is.finite(new_obj) && new_obj < best_obj) {
          best_obj <- new_obj
          best_params <- new_params
          if (verbose) {
            message(sprintf("    Eigenvector %d, dir=%+d, step=%.2f: obj %.4f -> %.4f",
                           idx, direction, step, current_obj, new_obj))
          }
        }
      }
    }
  }

  found_escape <- best_obj < current_obj - 1e-6
  if (verbose && found_escape) {
    message(sprintf("  Found escape direction: obj improved by %.4f", current_obj - best_obj))
  }

  return(list(new_params = best_params, found_escape = found_escape))
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
#' @param max_restarts Maximum number of eigenvector-based restarts for escaping
#'   saddle points (default 5). Set to 0 to disable.
#'
#' @details
#' For maximum efficiency, use \code{optimizer = "nlminb"} or \code{optimizer = "trust"}
#' which exploit the analytical Hessian computed in C++. The default L-BFGS methods
#' only use the gradient and approximate the Hessian from gradient history.
#'
#' When optimization fails to converge (possibly at a saddle point), the function
#' will attempt to escape by moving in the direction of negative Hessian eigenvalues.
#' This is controlled by the \code{max_restarts} parameter.
#'
#' @return List with parameter estimates, standard errors, log-likelihood, etc.
#' @export
estimate_model_rcpp <- function(model_system, data, init_params = NULL,
                                control = NULL, optimizer = "nlminb",
                                parallel = TRUE, verbose = TRUE,
                                max_restarts = 5) {

  # WORKAROUND: Deep copy model_system to avoid C++ reuse bug
  # Use serialize/unserialize for true deep copy
  model_system <- unserialize(serialize(model_system, NULL))

  # Default control if not provided
  if (is.null(control)) {
    control <- define_estimation_control()
  }

  # Convert data to matrix
  data_mat <- as.matrix(data)

  # Validate observation weights if specified
  if (!is.null(model_system$weights)) {
    weights_var <- model_system$weights
    if (!weights_var %in% colnames(data_mat)) {
      stop(sprintf("Weights variable '%s' not found in data", weights_var))
    }
    weights_vec <- data_mat[, weights_var]
    if (any(is.na(weights_vec))) {
      stop("Observation weights contain NA values")
    }
    if (any(weights_vec <= 0)) {
      warning("Some observation weights are <= 0. This may cause issues.")
    }
  }

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

  # Get initial parameters FIRST (needed for fixed coefficient values in C++)
  factor_variance_fixed <- NULL
  full_init_params <- NULL  # Full parameter vector (including fixed)

  if (is.null(init_params)) {
    # Use smart initialization based on separate component estimation
    init_result <- initialize_parameters(model_system, data, verbose = verbose)
    full_init_params <- init_result$init_params
    factor_variance_fixed <- init_result$factor_variance_fixed
  } else {
    # User provided init_params (assumed to be full vector)
    full_init_params <- init_params
  }

  # Initialize factor models on each worker
  if (verbose) message("Initializing C++ likelihood evaluators...")

  # Get n_quad from control
  n_quad <- control$n_quad_points

  if (!is.null(cl)) {
    # Export library paths to workers (ensures workers can find factorana)
    current_lib_paths <- .libPaths()
    parallel::clusterExport(cl, "current_lib_paths", envir = environment())
    parallel::clusterEvalQ(cl, {
      .libPaths(current_lib_paths)
      library(factorana)
    })

    # Export necessary objects to workers
    # Note: model_system no longer contains data (stripped in define_model_component)
    parallel::clusterExport(cl, c("model_system", "data_mat", "n_quad", "data_splits", "full_init_params"),
                           envir = environment())

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
      .fm_ptr <- initialize_factor_model_cpp(model_system, data_subset, n_quad, full_init_params)
    })

    fm_ptrs <- NULL  # Not used; pointers stored on workers

    # Set observation weights if specified in model_system
    if (!is.null(model_system$weights)) {
      weights_var <- model_system$weights
      # Export weights variable name
      parallel::clusterExport(cl, "weights_var", envir = environment())
      # Each worker extracts and sets weights for its data subset
      parallel::clusterEvalQ(cl, {
        worker_id <- .self_id
        idx <- data_splits[[worker_id]]
        weights_subset <- data_mat[idx, weights_var]
        set_observation_weights_cpp(.fm_ptr, weights_subset)
      })
      if (verbose) {
        weights_vec <- data_mat[, weights_var]
        message(sprintf("Using observation weights from '%s' (range: %.3f to %.3f)",
                       weights_var, min(weights_vec), max(weights_vec)))
      }
    }
  } else {
    # Single worker initialization
    fm_ptrs <- list(initialize_factor_model_cpp(model_system, data_mat, n_quad, full_init_params))

    # Set observation weights if specified in model_system
    if (!is.null(model_system$weights)) {
      weights_var <- model_system$weights
      weights_vec <- data_mat[, weights_var]
      set_observation_weights_cpp(fm_ptrs[[1]], weights_vec)
      if (verbose) {
        message(sprintf("Using observation weights from '%s' (range: %.3f to %.3f)",
                       weights_var, min(weights_vec), max(weights_vec)))
      }
    }
  }

  # Get parameter count and extract free parameters
  if (!is.null(cl)) {
    # For parallel case, get from first worker
    param_info <- parallel::clusterEvalQ(cl, get_parameter_info_cpp(.fm_ptr))[[1]]
  } else {
    param_info <- get_parameter_info_cpp(fm_ptrs[[1]])
  }
  n_params_free <- param_info$n_param_free
  n_params_total <- param_info$n_param

  # Extract free parameters from full init_params for the optimizer
  if (n_params_free == n_params_total) {
    # No fixed params - use full vector
    init_params <- full_init_params
  } else {
    # Some params are fixed - extract only free ones
    if (!is.null(cl)) {
      # For parallel case, use first worker to extract
      parallel::clusterExport(cl, "full_init_params", envir = environment())
      init_params <- parallel::clusterEvalQ(cl, {
        extract_free_params_cpp(.fm_ptr, full_init_params)
      })[[1]]
    } else {
      init_params <- extract_free_params_cpp(fm_ptrs[[1]], full_init_params)
    }

    if (verbose) {
      n_fixed <- n_params_total - n_params_free
      message(sprintf("  %d parameters total, %d fixed, %d free", n_params_total, n_fixed, n_params_free))
    }
  }

  # Validate parameter count
  if (length(init_params) != n_params_free) {
    stop(sprintf("Parameter initialization returned %d free parameters but C++ expects %d",
                 length(init_params), n_params_free))
  }

  # Build parameter metadata (names, types, etc.)
  param_metadata <- build_parameter_metadata(model_system)

  # Setup constraints and identify fixed/free parameters
  param_constraints <- setup_parameter_constraints(model_system, init_params, param_metadata, factor_variance_fixed, verbose)

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

  # Benchmark likelihood/gradient/Hessian computation before optimization
  if (verbose) {
    message("\nBenchmarking computation times (single evaluation)...")

    # Benchmark log-likelihood only
    t_loglik <- system.time({
      if (!is.null(cl)) {
        parallel::clusterExport(cl, "full_init_params", envir = environment())
        loglik_parts <- parallel::clusterEvalQ(cl, {
          evaluate_loglik_only_cpp(.fm_ptr, full_init_params)
        })
        loglik_test <- sum(unlist(loglik_parts))
      } else {
        loglik_test <- evaluate_loglik_only_cpp(fm_ptrs[[1]], full_init_params)
      }
    })[3]

    # Benchmark gradient
    t_grad <- system.time({
      if (!is.null(cl)) {
        grad_parts <- parallel::clusterEvalQ(cl, {
          result <- evaluate_likelihood_cpp(.fm_ptr, full_init_params,
                                           compute_gradient = TRUE,
                                           compute_hessian = FALSE)
          result$gradient
        })
      } else {
        result <- evaluate_likelihood_cpp(fm_ptrs[[1]], full_init_params,
                                         compute_gradient = TRUE,
                                         compute_hessian = FALSE)
      }
    })[3]

    # Benchmark Hessian
    t_hess <- system.time({
      if (!is.null(cl)) {
        hess_parts <- parallel::clusterEvalQ(cl, {
          result <- evaluate_likelihood_cpp(.fm_ptr, full_init_params,
                                           compute_gradient = FALSE,
                                           compute_hessian = TRUE)
          result$hessian
        })
      } else {
        result <- evaluate_likelihood_cpp(fm_ptrs[[1]], full_init_params,
                                         compute_gradient = FALSE,
                                         compute_hessian = TRUE)
      }
    })[3]

    message(sprintf("  Log-likelihood:  %.3f sec", t_loglik))
    message(sprintf("  Gradient:        %.3f sec", t_grad))
    message(sprintf("  Hessian:         %.3f sec", t_hess))
    message(sprintf("  Initial loglik:  %.4f", loglik_test))
    message("")
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
    hess_vec_len <- length(hess)
    expected_len <- n_params_full * (n_params_full + 1) / 2

    if (verbose && hess_vec_len != expected_len) {
      message(sprintf("WARNING: Hessian vector length mismatch!"))
      message(sprintf("  Expected: %d (for %d params)", expected_len, n_params_full))
      message(sprintf("  Actual: %d", hess_vec_len))
      message(sprintf("  Difference: %d", expected_len - hess_vec_len))
    }

    hess_mat_full <- matrix(0, n_params_full, n_params_full)
    idx <- 1
    for (i in 1:n_params_full) {
      for (j in i:n_params_full) {
        if (idx > hess_vec_len) {
          if (verbose) {
            message(sprintf("ERROR: Trying to access hess[%d] but length is only %d", idx, hess_vec_len))
            message(sprintf("  At position i=%d, j=%d", i, j))
          }
          stop("Hessian vector length insufficient for reconstruction")
        }
        hess_mat_full[i, j] <- hess[idx]
        hess_mat_full[j, i] <- hess[idx]  # Symmetrize
        idx <- idx + 1
      }
    }

    # Extract submatrix for free parameters only
    hess_mat_free <- hess_mat_full[param_constraints$free_idx, param_constraints$free_idx, drop = FALSE]

    # Debug: Check for NA/NaN values
    if (any(is.na(hess_mat_free)) || any(is.infinite(hess_mat_free))) {
      if (verbose) {
        message("WARNING: Hessian contains NA/NaN/Inf values!")
        message(sprintf("  NA count: %d", sum(is.na(hess_mat_free))))
        message(sprintf("  Inf count: %d", sum(is.infinite(hess_mat_free))))
        message(sprintf("  params_full: %s", paste(head(params_full, 5), collapse=", ")))
        message(sprintf("  hess vector (first 10): %s", paste(head(hess, 10), collapse=", ")))
      }
    }

    return(-hess_mat_free)  # Negative for minimization
  }

  # Run optimization with eigenvector-based restarts
  if (verbose) message(sprintf("Running optimization using %s...", optimizer))

  # Track current starting point for restarts
  current_start <- init_params[param_constraints$free_idx]
  n_restarts_used <- 0

  # Helper function to run one optimization attempt
  run_one_optimization <- function(start_params, verbose_opt) {
    if (optimizer == "nloptr") {
      if (!requireNamespace("nloptr", quietly = TRUE)) {
        stop("nloptr package required but not installed")
      }

      opt_result <- nloptr::nloptr(
        x0 = start_params,
        eval_f = objective_fn,
        eval_grad_f = gradient_fn,
        lb = param_constraints$lower_bounds_free,
        ub = param_constraints$upper_bounds_free,
        opts = list(
          algorithm = "NLOPT_LD_LBFGS",
          maxeval = 1000,
          xtol_rel = 1e-6,
          print_level = if (verbose_opt) 2 else 0
        )
      )
      return(list(
        par = opt_result$solution,
        value = opt_result$objective,
        convergence = opt_result$status
      ))

    } else if (optimizer == "optim") {
      opt_result <- stats::optim(
        par = start_params,
        fn = objective_fn,
        gr = gradient_fn,
        method = "L-BFGS-B",
        lower = param_constraints$lower_bounds_free,
        upper = param_constraints$upper_bounds_free,
        control = list(trace = if (verbose_opt) 1 else 0)
      )
      return(list(
        par = opt_result$par,
        value = opt_result$value,
        convergence = opt_result$convergence
      ))

    } else if (optimizer == "nlminb") {
      if (verbose_opt) message("  Using analytical gradient and Hessian")

      opt_result <- stats::nlminb(
        start = start_params,
        objective = objective_fn,
        gradient = gradient_fn,
        hessian = hessian_fn,
        lower = param_constraints$lower_bounds_free,
        upper = param_constraints$upper_bounds_free,
        control = list(
          trace = if (verbose_opt) 1 else 0,
          eval.max = 5000,
          iter.max = 10000
        )
      )
      return(list(
        par = opt_result$par,
        value = opt_result$objective,
        convergence = opt_result$convergence,
        iterations = opt_result$iterations,
        evaluations = opt_result$evaluations
      ))

    } else if (optimizer == "trust") {
      if (!requireNamespace("trustOptim", quietly = TRUE)) {
        stop("trustOptim package required but not installed. Install with: install.packages('trustOptim')")
      }

      if (any(param_constraints$param_fixed)) {
        warning("trust optimizer does not support fixed parameters. ",
                "Fixed parameters may drift from their initial values. ",
                "Consider using 'nlminb' or 'optim' instead.")
      }

      if (verbose_opt) message("  Using trust region with analytical gradient and Hessian")

      opt_result <- trustOptim::trust.optim(
        x = start_params,
        fn = objective_fn,
        gr = gradient_fn,
        hs = hessian_fn,
        method = "SR1",
        control = list(
          report.level = if (verbose_opt) 2 else 0,
          maxit = 500
        )
      )
      return(list(
        par = opt_result$solution,
        value = opt_result$value,
        convergence = if (opt_result$converged) 0 else 1
      ))

    } else {
      stop("Unknown optimizer: ", optimizer, "\n",
           "Available: 'nloptr' (L-BFGS), 'optim' (L-BFGS-B), 'nlminb' (uses Hessian), 'trust' (uses Hessian)")
    }
  }

  # Determine convergence success based on optimizer
  is_converged <- function(conv_code) {
    if (optimizer %in% c("optim", "nlminb")) {
      return(conv_code == 0)
    } else if (optimizer == "nloptr") {
      return(conv_code >= 1 && conv_code <= 4)  # nloptr success codes
    } else if (optimizer == "trust") {
      return(conv_code == 0)  # Already converted to 0 = success
    }
    return(FALSE)
  }

  # Run initial optimization
  opt_result <- run_one_optimization(current_start, verbose)
  estimates_free <- opt_result$par
  loglik <- -opt_result$value
  convergence <- opt_result$convergence

  # Eigenvector-based restart loop
  if (max_restarts > 0 && !is_converged(convergence)) {
    for (restart in seq_len(max_restarts)) {
      if (verbose) {
        message(sprintf("  Optimization did not converge. Attempting eigenvector-based restart %d/%d...",
                       restart, max_restarts))
      }

      # Try to find escape direction using eigenvector analysis
      escape_result <- find_saddle_escape_direction(
        params = estimates_free,
        hessian_fn = hessian_fn,
        objective_fn = objective_fn,
        param_constraints = param_constraints,
        verbose = verbose
      )

      if (!escape_result$found_escape) {
        if (verbose) message("  No escape direction found. Stopping restarts.")
        break
      }

      n_restarts_used <- restart

      # Restart from new point
      current_start <- escape_result$new_params
      opt_result <- run_one_optimization(current_start, verbose)
      estimates_free <- opt_result$par
      loglik <- -opt_result$value
      convergence <- opt_result$convergence

      if (verbose) {
        message(sprintf("  Restart %d: loglik = %.4f, conv = %d",
                       restart, loglik, convergence))
      }

      if (is_converged(convergence)) {
        if (verbose) message(sprintf("  Converged after %d restart(s)!", restart))
        break
      }
    }
  }

  # Reconstruct full parameter vector
  estimates <- init_params
  estimates[param_constraints$free_idx] <- estimates_free

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

  # Return results with class for print/summary methods
  result <- list(
    estimates = estimates,
    std_errors = std_errors,
    param_names = param_metadata$names,
    loglik = loglik,
    convergence = convergence,
    n_restarts = n_restarts_used,
    iterations = opt_result$iterations,
    evaluations = opt_result$evaluations,
    model_system = model_system,
    optimizer = optimizer
  )
  class(result) <- "factorana_result"
  result
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
