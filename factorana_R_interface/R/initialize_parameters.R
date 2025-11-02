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

  # ---- 1. Check factor identification ----
  # A factor variance is identified if ANY component has a non-NA, non-zero normalization
  n_factors <- model_system$factor$n_factors
  n_types <- model_system$factor$n_types

  factor_variance_identified <- rep(FALSE, n_factors)

  for (comp in model_system$components) {
    if (!is.null(comp$loading_normalization)) {
      for (k in seq_len(n_factors)) {
        # Check if this loading is fixed to a non-zero value
        if (!is.na(comp$loading_normalization[k]) &&
            abs(comp$loading_normalization[k]) > 1e-6) {
          factor_variance_identified[k] <- TRUE
        }
      }
    }
  }

  if (verbose) {
    message("Factor identification:")
    for (k in seq_len(n_factors)) {
      status <- if (factor_variance_identified[k]) "identified (variance will be estimated)" else "NOT identified (variance fixed to 1.0)"
      message(sprintf("  Factor %d: %s", k, status))
    }
  }

  # ---- 2. Initialize factor variances to 1.0 ----
  init_params <- rep(1.0, n_factors)

  # ---- 3. Estimate each component separately ----
  for (i_comp in seq_along(model_system$components)) {
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

    if (comp$model_type == "linear") {
      # Linear regression
      fit <- lm(outcome ~ X - 1)  # -1 because intercept is already in covariates
      coefs <- coef(fit)
      sigma <- summary(fit)$sigma

      if (verbose) {
        message(sprintf("  Linear model: %d covariates, sigma = %.4f", length(coefs), sigma))
      }

      # For multinomial logit with factors, we need parameters per choice
      comp_params <- c(coefs, rep(0.01, n_free_loadings), sigma)

    } else if (comp$model_type == "probit") {
      # Binary probit
      fit <- glm(outcome ~ X - 1, family = binomial(link = "probit"))
      coefs <- coef(fit)

      if (verbose) {
        message(sprintf("  Probit model: %d covariates", length(coefs)))
      }

      comp_params <- c(coefs, rep(0.01, n_free_loadings))

    } else if (comp$model_type == "logit") {
      if (comp$num_choices == 2) {
        # Binary logit
        fit <- glm(outcome ~ X - 1, family = binomial(link = "logit"))
        coefs <- coef(fit)

        if (verbose) {
          message(sprintf("  Binary logit: %d covariates", length(coefs)))
        }

        comp_params <- c(coefs, rep(0.01, n_free_loadings))

      } else {
        # Multinomial logit
        if (!requireNamespace("nnet", quietly = TRUE)) {
          stop("nnet package required for multinomial logit initialization. Install with: install.packages('nnet')")
        }

        fit <- nnet::multinom(outcome ~ X - 1, trace = FALSE)
        coefs_mat <- coef(fit)

        if (verbose) {
          message(sprintf("  Multinomial logit: %d choices, %d covariates", comp$num_choices, ncol(X)))
        }

        # Flatten parameters: for each choice (except reference), add covariates + loadings
        comp_params <- c()
        for (choice in seq_len(comp$num_choices - 1)) {
          # Get coefficients for this choice
          if (comp$num_choices == 2) {
            choice_coefs <- coefs_mat
          } else {
            choice_coefs <- coefs_mat[choice, ]
          }

          # Add covariates and loadings
          comp_params <- c(comp_params, choice_coefs, rep(0.01, n_free_loadings))
        }
      }

    } else if (comp$model_type == "oprobit") {
      # Ordered probit
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("MASS package required for ordered probit initialization. Install with: install.packages('MASS')")
      }

      fit <- MASS::polr(as.ordered(outcome) ~ X - 1, method = "probit")
      coefs <- coef(fit)
      thresholds <- fit$zeta

      if (verbose) {
        message(sprintf("  Ordered probit: %d covariates, %d thresholds", length(coefs), length(thresholds)))
      }

      comp_params <- c(coefs, rep(0.01, n_free_loadings), thresholds)
    }

    # Add component parameters to overall parameter vector
    init_params <- c(init_params, comp_params)
  }

  if (verbose) {
    message(sprintf("\nInitialized %d parameters total", length(init_params)))
  }

  list(
    init_params = init_params,
    factor_variance_fixed = !factor_variance_identified
  )
}
