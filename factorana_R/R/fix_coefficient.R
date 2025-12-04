#' Fix a coefficient in a model component
#'
#' Constrains a regression coefficient (beta) to a fixed value during estimation.
#' The coefficient will not be optimized - it remains at the specified value.
#'
#' @param component A model_component object from define_model_component()
#' @param covariate Character. Name of the covariate whose coefficient to fix.
#'   Must be one of the covariates specified when creating the component.
#' @param value Numeric. The value to fix the coefficient to.
#' @param choice Integer (optional). For multinomial logit models with num_choices > 2,
#'   specifies which choice's coefficient to fix. Must be between 1 and (num_choices - 1).
#'   Choice 0 (reference category) has no parameters.
#'   For other model types, this should be NULL (default).
#'
#' @return The modified model_component object with the fixed coefficient constraint added.
#'
#' @details
#' Fixed coefficients are stored in the component and used during:
#' \itemize{
#'   \item Parameter initialization: fixed values are used directly (no estimation)
#'   \item Optimization: fixed parameters are excluded from the optimization
#'   \item Likelihood evaluation: fixed values are inserted into the parameter vector
#' }
#'
#' Note: This function only fixes regression coefficients (betas), not:
#' \itemize{
#'   \item Factor loadings (use loading_normalization in define_model_component)
#'   \item Sigma (for linear models)
#'   \item Thresholds (for ordered probit)
#' }
#'
#' @examples
#' \dontrun{
#' # Fix intercept to 0
#' mc <- define_model_component(name = "wage", ...)
#' mc <- fix_coefficient(mc, covariate = "intercept", value = 0.0)
#'
#' # Fix a coefficient to a specific value
#' mc <- fix_coefficient(mc, covariate = "education", value = 0.1)
#'
#' # For multinomial logit: fix coefficient for choice 2
#' mc <- define_model_component(name = "occupation", ..., num_choices = 4)
#' mc <- fix_coefficient(mc, covariate = "age", value = 0.05, choice = 2)
#' }
#'
#' @export
fix_coefficient <- function(component, covariate, value, choice = NULL) {


  # ---- 1. Validate component ----
  if (!inherits(component, "model_component")) {
    stop("`component` must be an object of class 'model_component'.")
  }

  # ---- 2. Validate covariate ----
  if (!is.character(covariate) || length(covariate) != 1L) {
    stop("`covariate` must be a single character string.")
  }

  if (!(covariate %in% component$covariates)) {
    stop("Covariate '", covariate, "' not found in component. ",
         "Available covariates: ", paste(component$covariates, collapse = ", "))
  }

  # ---- 3. Validate value ----
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop("`value` must be a single finite numeric value.")
  }

  # ---- 4. Validate choice (for multinomial logit) ----
  if (component$model_type == "logit" && component$num_choices > 2) {
    # Multinomial logit: choice is required or defaults to applying to all choices
    if (!is.null(choice)) {
      choice <- as.integer(choice)
      if (!is.finite(choice) || choice < 1L || choice > (component$num_choices - 1L)) {
        stop("`choice` must be between 1 and ", component$num_choices - 1L,
             " (choice 0 is the reference category with no parameters).")
      }
    }
  } else {
    # Non-multinomial models: choice must be NULL
    if (!is.null(choice)) {
      stop("`choice` argument is only valid for multinomial logit models with num_choices > 2.")
    }
  }

  # ---- 5. Check for duplicate constraints ----
  for (fc in component$fixed_coefficients) {
    if (fc$covariate == covariate && identical(fc$choice, choice)) {
      stop("Coefficient for covariate '", covariate, "'",
           if (!is.null(choice)) paste0(" (choice ", choice, ")") else "",
           " is already fixed. Remove and recreate component to change.")
    }
  }

  # ---- 6. Add the constraint ----
  new_constraint <- list(
    covariate = covariate,
    value = value,
    choice = choice
  )

  component$fixed_coefficients <- c(component$fixed_coefficients, list(new_constraint))

  # ---- 7. Update parameter count ----
  # Reduce nparam_model by 1 for each fixed coefficient
  component$nparam_model <- component$nparam_model - 1L

  return(component)
}


#' Get the number of fixed coefficients in a model component
#'
#' @param component A model_component object
#' @return Integer count of fixed coefficients
#' @keywords internal
n_fixed_coefficients <- function(component) {
  length(component$fixed_coefficients)
}


#' Check if a covariate's coefficient is fixed
#'
#' @param component A model_component object
#' @param covariate Character. Name of the covariate to check.
#' @param choice Integer (optional). For multinomial logit, which choice to check.
#' @return Logical TRUE if fixed, FALSE otherwise
#' @keywords internal
is_coefficient_fixed <- function(component, covariate, choice = NULL) {
  for (fc in component$fixed_coefficients) {
    if (fc$covariate == covariate && identical(fc$choice, choice)) {
      return(TRUE)
    }
  }
  return(FALSE)
}


#' Get the fixed value for a covariate's coefficient
#'
#' @param component A model_component object
#' @param covariate Character. Name of the covariate.
#' @param choice Integer (optional). For multinomial logit, which choice.
#' @return The fixed value, or NULL if not fixed
#' @keywords internal
get_fixed_coefficient_value <- function(component, covariate, choice = NULL) {
  for (fc in component$fixed_coefficients) {
    if (fc$covariate == covariate && identical(fc$choice, choice)) {
      return(fc$value)
    }
  }
  return(NULL)
}
