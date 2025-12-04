#' Initialize parameters for factor model estimation
#'
#' Estimates each model component separately in R (ignoring factors) to obtain
#' good starting values. Also checks factor identification.
#'
#' @param model_system A model_system object from define_model_system()
#' @param data Data frame containing all variables
#' @param verbose Whether to print progress (default TRUE)
#'
#' @return List with:
#'   \itemize{
#'     \item \code{init_params} - Initial parameter values
#'     \item \code{factor_variance_fixed} - Logical vector indicating which factor variances must be fixed
#'   }
#'
#' @details
#' Factor identification: For each factor, if NO component has a non-zero fixed loading,
#' then the factor variance is not identified and must be fixed to 1.0.
#'
#' @export
initialize_parameters <- function(model_system, data, verbose = TRUE) {

  if (verbose) {
    message("Initializing parameters...")
  }

  # Get n_types from factor model (needed in all code paths)
  n_types <- model_system$factor$n_types
  if (is.null(n_types)) n_types <- 1L

  # ---- 1. Handle previous_stage if present ----
  if (!is.null(model_system$previous_stage_info)) {
    # Use previous-stage parameter values directly
    n_fixed_comps <- model_system$previous_stage_info$n_components
    init_params <- model_system$previous_stage_info$fixed_param_values
    # Inherit names from previous stage if available
    param_names <- if (!is.null(model_system$previous_stage_info$param_names)) {
      model_system$previous_stage_info$param_names
    } else {
      names(init_params)
    }
    if (is.null(param_names)) param_names <- character(0)

    if (verbose) {
      message(sprintf("Using %d fixed parameters from previous stage", length(init_params)))
      message(sprintf("  Previous stage had %d components", n_fixed_comps))
    }

    # Initialize only the new (second-stage) components
    start_comp_idx <- n_fixed_comps + 1
  } else {
    # Standard initialization: start from scratch
    n_factors <- model_system$factor$n_factors

    # Check factor identification
    factor_variance_fixed <- rep(FALSE, n_factors)

    for (comp in model_system$components) {
      if (!is.null(comp$loading_normalization)) {
        for (k in seq_len(n_factors)) {
          # Check if this loading is fixed to 1.0 (identification via unit loading)
          if (!is.na(comp$loading_normalization[k]) &&
              abs(comp$loading_normalization[k] - 1.0) < 1e-6) {
            factor_variance_fixed[k] <- TRUE
          }
        }
      }
    }

    if (verbose) {
      message("Factor identification:")
      for (k in seq_len(n_factors)) {
        status <- if (factor_variance_fixed[k]) "identified by fixed loading (variance will be estimated)" else "NOT identified (variance fixed to 1.0)"
        message(sprintf("  Factor %d: %s", k, status))
      }
    }

    # Initialize factor variances to 1.0
    init_params <- rep(1.0, n_factors)
    param_names <- paste0("factor_var_", seq_len(n_factors))

    # Add correlation parameters if correlation = TRUE
    if (isTRUE(model_system$factor$correlation) && n_factors == 2) {
      # For 2-factor correlated model, add one correlation parameter
      # Initialize to 0 (uncorrelated) as a neutral starting point
      init_params <- c(init_params, 0.0)
      param_names <- c(param_names, "factor_corr_1_2")
    }

    # Add type model parameters if n_types > 1
    # Type model: log(P(type=t)/P(type=1)) = sum_k lambda_t_k * f_k
    # (n_types - 1) * n_factors parameters (type 1 is reference)
    if (n_types > 1L) {
      type_loadings <- rep(0.0, (n_types - 1L) * n_factors)
      type_loading_names <- character(0)
      for (t in 2:n_types) {
        for (k in seq_len(n_factors)) {
          type_loading_names <- c(type_loading_names, paste0("type_", t, "_loading_", k))
        }
      }
      init_params <- c(init_params, type_loadings)
      param_names <- c(param_names, type_loading_names)
    }

    start_comp_idx <- 1
  }

  # ---- 2. Estimate each component separately ----
  # If previous_stage, only estimate new components
  for (i_comp in start_comp_idx:length(model_system$components)) {
    comp <- model_system$components[[i_comp]]

    if (verbose) {
      message(sprintf("\nEstimating component %d (%s)...", i_comp, comp$name))
    }

    # Apply evaluation indicator if present
    if (!is.null(comp$evaluation_indicator)) {
      idx <- data[[comp$evaluation_indicator]] == 1 & !is.na(data[[comp$evaluation_indicator]])
      comp_data <- data[idx, , drop = FALSE]
    } else {
      comp_data <- data
    }

    outcome <- comp_data[[comp$outcome]]
    X <- as.matrix(comp_data[, comp$covariates, drop = FALSE])

    # Count free factor loadings for this component
    n_free_loadings <- sum(is.na(comp$loading_normalization))

    # Estimate model depending on type
    comp_params <- NULL
    comp_param_names <- NULL

    # Get fixed coefficients list (may be empty)
    fixed_coefs <- if (!is.null(comp$fixed_coefficients)) comp$fixed_coefficients else list()

    if (comp$model_type == "linear") {
      # Linear regression
      fit <- lm(outcome ~ X - 1)  # -1 because intercept is already in covariates
      coefs <- coef(fit)
      sigma <- summary(fit)$sigma

      # Apply any fixed coefficient values
      coefs <- apply_fixed_coefficients(coefs, comp$covariates, fixed_coefs)

      if (verbose) {
        n_fixed <- length(fixed_coefs)
        msg <- sprintf("  Linear model: %d covariates, sigma = %.4f", length(coefs), sigma)
        if (n_fixed > 0) msg <- paste0(msg, sprintf(" (%d fixed)", n_fixed))
        message(msg)
      }

      # For multinomial logit with factors, we need parameters per choice
      comp_params <- c(coefs, rep(0.5, n_free_loadings), sigma)

      # Build parameter names
      coef_names <- paste0(comp$name, "_", comp$covariates)
      loading_names <- character(0)
      if (n_free_loadings > 0) {
        free_factor_idx <- which(is.na(comp$loading_normalization))
        loading_names <- paste0(comp$name, "_loading_", free_factor_idx)
      }
      comp_param_names <- c(coef_names, loading_names, paste0(comp$name, "_sigma"))

      # Add type-specific intercepts if n_types > 1
      # Initialize to small non-zero values to make types distinguishable
      # (avoids degenerate Hessian at initialization)
      if (n_types > 1L) {
        type_intercepts <- 0.1 * (seq_len(n_types - 1L))  # 0.1, 0.2, 0.3, ...
        type_intercept_names <- paste0(comp$name, "_type_", 2:n_types, "_intercept")
        comp_params <- c(comp_params, type_intercepts)
        comp_param_names <- c(comp_param_names, type_intercept_names)
      }

    } else if (comp$model_type == "probit") {
      # Binary probit
      fit <- glm(outcome ~ X - 1, family = binomial(link = "probit"))
      coefs <- coef(fit)

      # Apply any fixed coefficient values
      coefs <- apply_fixed_coefficients(coefs, comp$covariates, fixed_coefs)

      if (verbose) {
        n_fixed <- length(fixed_coefs)
        msg <- sprintf("  Probit model: %d covariates", length(coefs))
        if (n_fixed > 0) msg <- paste0(msg, sprintf(" (%d fixed)", n_fixed))
        message(msg)
      }

      comp_params <- c(coefs, rep(0.5, n_free_loadings))

      # Build parameter names
      coef_names <- paste0(comp$name, "_", comp$covariates)
      loading_names <- character(0)
      if (n_free_loadings > 0) {
        free_factor_idx <- which(is.na(comp$loading_normalization))
        loading_names <- paste0(comp$name, "_loading_", free_factor_idx)
      }
      comp_param_names <- c(coef_names, loading_names)

      # Add type-specific intercepts if n_types > 1
      # Initialize to small non-zero values to make types distinguishable
      if (n_types > 1L) {
        type_intercepts <- 0.1 * (seq_len(n_types - 1L))  # 0.1, 0.2, 0.3, ...
        type_intercept_names <- paste0(comp$name, "_type_", 2:n_types, "_intercept")
        comp_params <- c(comp_params, type_intercepts)
        comp_param_names <- c(comp_param_names, type_intercept_names)
      }

    } else if (comp$model_type == "logit") {
      if (comp$num_choices == 2) {
        # Binary logit
        fit <- glm(outcome ~ X - 1, family = binomial(link = "logit"))
        coefs <- coef(fit)

        # Apply any fixed coefficient values
        coefs <- apply_fixed_coefficients(coefs, comp$covariates, fixed_coefs)

        if (verbose) {
          n_fixed <- length(fixed_coefs)
          msg <- sprintf("  Binary logit: %d covariates", length(coefs))
          if (n_fixed > 0) msg <- paste0(msg, sprintf(" (%d fixed)", n_fixed))
          message(msg)
        }

        comp_params <- c(coefs, rep(0.5, n_free_loadings))

        # Build parameter names
        coef_names <- paste0(comp$name, "_", comp$covariates)
        loading_names <- character(0)
        if (n_free_loadings > 0) {
          free_factor_idx <- which(is.na(comp$loading_normalization))
          loading_names <- paste0(comp$name, "_loading_", free_factor_idx)
        }
        comp_param_names <- c(coef_names, loading_names)

        # Add type-specific intercepts if n_types > 1
        # Initialize to small non-zero values to make types distinguishable
        if (n_types > 1L) {
          type_intercepts <- 0.1 * (seq_len(n_types - 1L))  # 0.1, 0.2, 0.3, ...
          type_intercept_names <- paste0(comp$name, "_type_", 2:n_types, "_intercept")
          comp_params <- c(comp_params, type_intercepts)
          comp_param_names <- c(comp_param_names, type_intercept_names)
        }

      } else {
        # Multinomial logit
        if (!requireNamespace("nnet", quietly = TRUE)) {
          stop("nnet package required for multinomial logit initialization. Install with: install.packages('nnet')")
        }

        fit <- nnet::multinom(outcome ~ X - 1, trace = FALSE)
        coefs_mat <- coef(fit)

        if (verbose) {
          n_fixed <- length(fixed_coefs)
          msg <- sprintf("  Multinomial logit: %d choices, %d covariates", comp$num_choices, ncol(X))
          if (n_fixed > 0) msg <- paste0(msg, sprintf(" (%d fixed)", n_fixed))
          message(msg)
        }

        # Flatten parameters: for each choice (except reference), add covariates + loadings
        comp_params <- c()
        comp_param_names <- c()
        for (choice in seq_len(comp$num_choices - 1)) {
          # Get coefficients for this choice
          if (comp$num_choices == 2) {
            choice_coefs <- coefs_mat
          } else {
            choice_coefs <- coefs_mat[choice, ]
          }

          # Apply any fixed coefficient values for this choice
          choice_coefs <- apply_fixed_coefficients(choice_coefs, comp$covariates, fixed_coefs, choice = choice)

          # Add covariates and loadings
          comp_params <- c(comp_params, choice_coefs, rep(0.5, n_free_loadings))

          # Build parameter names for this choice
          coef_names <- paste0(comp$name, "_c", choice, "_", comp$covariates)
          loading_names <- character(0)
          if (n_free_loadings > 0) {
            free_factor_idx <- which(is.na(comp$loading_normalization))
            loading_names <- paste0(comp$name, "_c", choice, "_loading_", free_factor_idx)
          }
          comp_param_names <- c(comp_param_names, coef_names, loading_names)
        }

        # Add type-specific intercepts if n_types > 1
        # For multinomial logit, each non-reference choice gets type-specific intercepts
        # Initialize to small non-zero values to make types distinguishable
        if (n_types > 1L) {
          for (choice in seq_len(comp$num_choices - 1)) {
            type_intercepts <- 0.1 * (seq_len(n_types - 1L))  # 0.1, 0.2, 0.3, ...
            type_intercept_names <- paste0(comp$name, "_c", choice, "_type_", 2:n_types, "_intercept")
            comp_params <- c(comp_params, type_intercepts)
            comp_param_names <- c(comp_param_names, type_intercept_names)
          }
        }
      }

    } else if (comp$model_type == "oprobit") {
      # Ordered probit
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("MASS package required for ordered probit initialization. Install with: install.packages('MASS')")
      }

      # For ordered probit, intercept is not identified (absorbed into thresholds)
      # Remove intercept column if present
      X_no_int <- X
      intercept_col <- which(tolower(comp$covariates) == "intercept")
      if (length(intercept_col) > 0) {
        X_no_int <- X[, -intercept_col, drop = FALSE]
      }

      # Fit model without intercept
      # Suppress MASS::polr warning about intercept (it's expected for ordered models)
      if (ncol(X_no_int) > 0) {
        fit <- suppressWarnings(
          MASS::polr(as.ordered(outcome) ~ X_no_int - 1, method = "probit")
        )
        coefs_no_int <- coef(fit)
      } else {
        # No covariates: fit intercept-only model to get thresholds
        fit <- suppressWarnings(
          MASS::polr(as.ordered(outcome) ~ 1, method = "probit")
        )
        coefs_no_int <- numeric(0)
      }
      thresholds_abs <- fit$zeta

      # Convert absolute thresholds to incremental parameterization
      # thresh1, thresh2, thresh3 -> thresh1, (thresh2-thresh1), (thresh3-thresh2)
      thresholds <- c(thresholds_abs[1])
      if (length(thresholds_abs) > 1) {
        for (i in 2:length(thresholds_abs)) {
          thresholds <- c(thresholds, thresholds_abs[i] - thresholds_abs[i-1])
        }
      }

      # Scale thresholds to match factor model variance at initialization
      # MASS::polr assumes Var(y*) = 1, but the factor model has:
      # Var(y*) = sum(lambda_k^2 * sigma_k^2) + 1
      # where lambda_k are loadings and sigma_k^2 are factor variances (all = 1.0 at init)

      # Calculate expected variance contribution from factors at initialization
      variance_from_factors <- 0.0

      if (!is.null(comp$loading_normalization)) {
        for (k in seq_along(comp$loading_normalization)) {
          loading_value <- comp$loading_normalization[k]

          if (is.na(loading_value)) {
            # Free loading: will be initialized to 0.5
            variance_from_factors <- variance_from_factors + (0.5^2) * 1.0
          } else if (abs(loading_value) > 1e-6) {
            # Fixed non-zero loading (e.g., 1.0)
            variance_from_factors <- variance_from_factors + (loading_value^2) * 1.0
          }
          # Fixed zero loadings contribute 0
        }
      }

      # Total variance at initialization: factor contributions + residual (1.0)
      total_variance <- variance_from_factors + 1.0

      # Scale thresholds if total variance differs from 1.0 (polr's assumption)
      if (abs(total_variance - 1.0) > 0.01) {
        scale_factor <- sqrt(total_variance)
        thresholds <- thresholds * scale_factor

        if (verbose) {
          message(sprintf("  Scaling thresholds by %.3f (factor model variance = %.2f at initialization)",
                          scale_factor, total_variance))
        }
      }

      # Build parameter vector: put intercept (0.0) first if it was in covariates
      # Note: For oprobit, intercept is theoretically absorbed into thresholds,
      # but we keep it in the parameter vector (as 0.0) to match C++ expectations
      if (length(intercept_col) > 0) {
        # Reconstruct in original order: intercept gets 0.0, others get fitted values
        coefs <- rep(0.0, ncol(X))
        coefs[-intercept_col] <- coefs_no_int
      } else {
        coefs <- coefs_no_int
      }

      # Apply any fixed coefficient values
      coefs <- apply_fixed_coefficients(coefs, comp$covariates, fixed_coefs)

      if (verbose) {
        n_fixed <- length(fixed_coefs)
        msg <- sprintf("  Ordered probit: %d covariates, %d thresholds (incremental form)", length(coefs), length(thresholds))
        if (n_fixed > 0) msg <- paste0(msg, sprintf(" (%d fixed)", n_fixed))
        message(msg)
      }

      comp_params <- c(coefs, rep(0.5, n_free_loadings), thresholds)

      # Build parameter names
      coef_names <- character(0)
      if (length(coefs) > 0 && length(comp$covariates) > 0) {
        coef_names <- paste0(comp$name, "_", comp$covariates)
      }
      loading_names <- character(0)
      if (n_free_loadings > 0) {
        free_factor_idx <- which(is.na(comp$loading_normalization))
        loading_names <- paste0(comp$name, "_loading_", free_factor_idx)
      }
      n_thresholds <- length(thresholds)
      threshold_names <- paste0(comp$name, "_thresh_", seq_len(n_thresholds))
      comp_param_names <- c(coef_names, loading_names, threshold_names)

      # Add type-specific intercepts if n_types > 1
      # For oprobit, type-specific intercepts shift all thresholds by a constant
      # Initialize to small non-zero values to make types distinguishable
      if (n_types > 1L) {
        type_intercepts <- 0.1 * (seq_len(n_types - 1L))  # 0.1, 0.2, 0.3, ...
        type_intercept_names <- paste0(comp$name, "_type_", 2:n_types, "_intercept")
        comp_params <- c(comp_params, type_intercepts)
        comp_param_names <- c(comp_param_names, type_intercept_names)
      }
    }

    # Add component parameters to overall parameter vector
    init_params <- c(init_params, comp_params)
    param_names <- c(param_names, comp_param_names)
  }

  if (verbose) {
    message(sprintf("\nInitialized %d parameters total", length(init_params)))
  }

  # Determine factor_variance_fixed status
  if (!is.null(model_system$previous_stage_info)) {
    # For previous_stage, factor variance is always fixed
    factor_variance_fixed_status <- TRUE
  } else {
    # Standard case: use the computed status (vector for multifactor models)
    factor_variance_fixed_status <- factor_variance_fixed
  }

  # Assign names to parameter vector
  names(init_params) <- param_names

  list(
    init_params = init_params,
    param_names = param_names,
    factor_variance_fixed = factor_variance_fixed_status
  )
}


#' Apply fixed coefficient values to estimated coefficients
#'
#' Helper function to replace estimated coefficients with fixed values
#' where specified in the component's fixed_coefficients list.
#'
#' @param coefs Named numeric vector of estimated coefficients
#' @param covariates Character vector of covariate names (in order)
#' @param fixed_coefficients List of fixed coefficient constraints from component
#' @param choice Integer or NULL. For multinomial logit, which choice these coefs belong to.
#' @return Numeric vector with fixed values applied
#' @keywords internal
apply_fixed_coefficients <- function(coefs, covariates, fixed_coefficients, choice = NULL) {
  if (length(fixed_coefficients) == 0) {
    return(coefs)
  }

  for (fc in fixed_coefficients) {
    # Check if this constraint applies (matching choice for mlogit)
    if (!identical(fc$choice, choice)) {
      next
    }

    # Find the position of this covariate
    pos <- match(fc$covariate, covariates)
    if (!is.na(pos)) {
      coefs[pos] <- fc$value
    }
  }

  return(coefs)
}
