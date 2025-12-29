#' Define a dynamic factor model component
#'
#' Creates a structural equation between latent factors, where one factor
#' is the "outcome" that depends on other factors and optional covariates.
#' This is useful for specifying dynamic relationships in factor models.
#'
#' The model form is:
#' \deqn{f_{outcome} = X'\beta + \sum_{k \neq outcome} \lambda_k f_k + \epsilon}
#'
#' With `factor_spec`, higher-order terms can be added (quadratic and interactions),
#' but only for non-outcome factors.
#'
#' @param name Character. Name of the dynamic model component.
#' @param data data.frame for validation and covariates.
#' @param outcome_factor Integer. 1-based index of the factor that is the outcome
#'   (left-hand side of the structural equation).
#' @param factor Object of class `factor_model`. Must have `n_factors >= 2`.
#' @param covariates Character vector (optional). Names of covariates from `data`.
#'   Can be NULL for a pure factor-to-factor relationship.
#' @param dyn_type Character. Type of dynamic specification. Currently only
#'   `"linear"` is supported. Future options: `"cobb_douglas"`, `"ces"`.
#' @param factor_spec Character. Specification for non-outcome factor terms:
#'   - `"linear"` (default): Only linear factor terms (\eqn{\lambda_k f_k})
#'   - `"quadratic"`: Linear + quadratic terms (\eqn{\lambda_k f_k + \lambda_{qk} f_k^2})
#'   - `"interactions"`: Linear + interaction terms (\eqn{\lambda_k f_k + \lambda_{jk} f_j f_k})
#'   - `"full"`: Linear + quadratic + interaction terms
#'   Note: Interaction terms require `n_factors >= 3` (at least 2 non-outcome factors).
#' @param intercept Logical. Whether to include an intercept (default = TRUE).
#' @param evaluation_indicator Character (optional). Variable name for subsample
#'   evaluation (1 = include, 0 = exclude).
#'
#' @return An object of class `c("dyn_model_component", "model_component")`.
#'
#' @details
#' The dynamic factor model component wraps the linear model internally by:
#' \itemize{
#'   \item Creating a dummy outcome variable of zeros
#'   \item Setting the loading on the outcome factor to -1 (fixed)
#'   \item Estimating loadings on all other factors
#' }
#'
#' The equation \eqn{f_{out} = X'\beta + \sum \lambda_k f_k + \epsilon} is
#' equivalent to \eqn{0 = X'\beta + \sum \lambda_k f_k - f_{out} + \epsilon}.
#'
#' @examples
#' \dontrun{
#' # Define a 3-factor model
#' fm <- define_factor_model(n_factors = 3, n_types = 1)
#'
#' # Define measurement equations for factors 2 and 3
#' mc1 <- define_model_component("Y1", data, "y1", fm, covariates = "intercept",
#'                               model_type = "linear", loading_normalization = c(0, 1, 0))
#' mc2 <- define_model_component("Y2", data, "y2", fm, covariates = "intercept",
#'                               model_type = "linear", loading_normalization = c(0, 0, 1))
#'
#' # Define structural equation: f1 = beta0 + lambda2*f2 + lambda3*f3 + epsilon
#' dyn <- define_dyn_model_component("structural", data, outcome_factor = 1, fm,
#'                                    intercept = TRUE)
#'
#' # Build model system
#' ms <- define_model_system(factor = fm, components = list(mc1, mc2, dyn))
#' }
#'
#' @export
define_dyn_model_component <- function(name,
                                        data,
                                        outcome_factor,
                                        factor,
                                        covariates = NULL,
                                        dyn_type = c("linear"),
                                        factor_spec = c("linear", "quadratic", "interactions", "full"),
                                        intercept = TRUE,
                                        evaluation_indicator = NULL) {

  # ---- 1. Basic argument validation ----

  # Check data
 if (!is.data.frame(data)) stop("`data` must be a data.frame (or data.table/tibble).")
  data <- as.data.frame(data)

  # Check factor model
  if (!inherits(factor, "factor_model")) {
    stop("`factor` must be an object of class 'factor_model'.")
  }

  k <- as.integer(factor$n_factors)
  if (is.na(k) || k < 2L) {
    stop("Dynamic factor model requires n_factors >= 2.")
  }

  # Check outcome_factor
  if (!is.numeric(outcome_factor) || length(outcome_factor) != 1L) {
    stop("`outcome_factor` must be a single integer.")
  }
  outcome_factor <- as.integer(outcome_factor)
  if (outcome_factor < 1L || outcome_factor > k) {
    stop("`outcome_factor` must be between 1 and n_factors (", k, ").")
  }

  # Match arguments
  dyn_type <- match.arg(dyn_type)
  factor_spec <- match.arg(factor_spec)

  # ---- 2. Validate factor_spec for non-outcome factors ----
  # For interactions, we need at least 2 non-outcome factors (k-1 >= 2, i.e., k >= 3)

  n_non_outcome <- k - 1L
  if (factor_spec %in% c("interactions", "full") && n_non_outcome < 2L) {
    warning("factor_spec='", factor_spec, "' requires at least 2 non-outcome factors (n_factors >= 3). ",
            "Downgrading to '", if (factor_spec == "full") "quadratic" else "linear", "'.")
    factor_spec <- if (factor_spec == "full") "quadratic" else "linear"
  }

  # ---- 3. Handle intercept and covariates ----

  if (!is.logical(intercept) || length(intercept) != 1L) {
    stop("`intercept` must be a single logical value.")
  }

  # Process covariates
  if (is.null(covariates)) {
    covariates <- character(0)
  } else if (!is.character(covariates)) {
    stop("`covariates` must be NULL or a character vector.")
  }

  # Add intercept if requested and not already present
  if (intercept) {
    if (!("intercept" %in% names(data))) {
      data[["intercept"]] <- 1
    }
    if (!("intercept" %in% covariates)) {
      covariates <- c("intercept", covariates)
    }
  }

  # Validate covariates exist in data
  if (length(covariates) > 0) {
    missing_covs <- setdiff(covariates, names(data))
    if (length(missing_covs) > 0) {
      stop("Covariates not found in `data`: ", paste(missing_covs, collapse = ", "))
    }
  }

  # ---- 4. Validate evaluation indicator ----

  if (!is.null(evaluation_indicator)) {
    if (!is.character(evaluation_indicator) || length(evaluation_indicator) != 1L) {
      stop("`evaluation_indicator` must be a single column name (character).")
    }
    if (!(evaluation_indicator %in% names(data))) {
      stop("`evaluation_indicator` not found in `data`.")
    }
  }

  # ---- 5. Create dummy outcome column ----
  # The dynamic model uses Y=0, with loading on outcome factor = -1
  # This gives: 0 = Xβ + Σλ_k*f_k + (-1)*f_out + ε
  # Which rearranges to: f_out = Xβ + Σλ_k*f_k + ε

  dummy_outcome_name <- paste0(".__dyn_", name, "__")
  data[[dummy_outcome_name]] <- 0

  # ---- 6. Build loading_normalization ----
  # Outcome factor gets -1 (fixed), all others are free (NA)

  loading_norm <- rep(NA_real_, k)
  loading_norm[outcome_factor] <- -1

  # ---- 7. Call define_model_component ----
  # Use the linear model with our special setup

  component <- define_model_component(
    name = name,
    data = data,
    outcome = dummy_outcome_name,
    factor = factor,
    evaluation_indicator = evaluation_indicator,
    covariates = if (length(covariates) > 0) covariates else NULL,
    model_type = "linear",
    intercept = FALSE,  # Already handled above
    loading_normalization = loading_norm,
    factor_spec = factor_spec
  )

  # ---- 8. Adjust parameter counts for dynamic model ----
  # Quadratic and interaction terms only apply to non-outcome factors

  # Recalculate counts based on non-outcome factors only
  if (factor_spec %in% c("quadratic", "full")) {
    # k-1 quadratic terms (one per non-outcome factor)
    component$n_quadratic_loadings <- n_non_outcome
  }

  if (factor_spec %in% c("interactions", "full") && n_non_outcome >= 2L) {
    # (k-1)*(k-2)/2 interaction terms among non-outcome factors
    component$n_interaction_loadings <- as.integer(n_non_outcome * (n_non_outcome - 1L) / 2L)
  }

  # Recalculate nparam_model
  # For linear dynamic model: betas + free_loadings + quadratic + interactions + sigma
  n_free_loadings <- sum(is.na(loading_norm))  # k-1 (all except outcome factor)
  n_second_order <- component$n_quadratic_loadings + component$n_interaction_loadings
  component$nparam_model <- length(covariates) + n_free_loadings + n_second_order + 1L  # +1 for sigma

  # ---- 9. Add dynamic-specific attributes ----

  component$is_dynamic <- TRUE
  component$outcome_factor <- outcome_factor
  component$dyn_type <- dyn_type

  # ---- 10. Set class and return ----

  class(component) <- c("dyn_model_component", "model_component")
  return(component)
}


#' Print method for dynamic model component
#'
#' @param x A dyn_model_component object
#' @param ... Additional arguments (ignored)
#' @export
print.dyn_model_component <- function(x, ...) {
  cat("Dynamic Factor Model Component\n")
  cat("==============================\n")
  cat("Name:", x$name, "\n")
  cat("Type:", x$dyn_type, "\n")
  cat("Outcome Factor:", x$outcome_factor, "of", x$k, "factors\n")
  cat("Factor Spec:", x$factor_spec, "\n")

  if (length(x$covariates) > 0) {
    cat("Covariates:", paste(x$covariates, collapse = ", "), "\n")
  } else {
    cat("Covariates: (none - pure factor-to-factor)\n")
  }

  cat("Observations:", x$n_obs, "\n")
  cat("Parameters:", x$nparam_model, "\n")

  if (!is.null(x$evaluation_indicator)) {
    cat("Evaluation indicator:", x$evaluation_indicator, "\n")
  }

  # Show loading structure
  cat("\nFactor Loadings:\n")
  for (i in seq_len(x$k)) {
    if (i == x$outcome_factor) {
      cat(sprintf("  Factor %d: -1 (outcome, fixed)\n", i))
    } else {
      cat(sprintf("  Factor %d: free\n", i))
    }
  }

  if (x$n_quadratic_loadings > 0) {
    cat(sprintf("\nQuadratic terms: %d (non-outcome factors)\n", x$n_quadratic_loadings))
  }
  if (x$n_interaction_loadings > 0) {
    cat(sprintf("Interaction terms: %d (among non-outcome factors)\n", x$n_interaction_loadings))
  }

  invisible(x)
}
