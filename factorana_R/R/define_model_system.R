#' Define a model system
#'
#' @param components A named list of model_component objects
#' @param factor A factor_model object
#' @param previous_stage Optional result from a previous estimation stage.
#'   If provided, the previous stage components and parameters will be fixed
#'   and prepended to the new components. This enables sequential/multi-stage
#'   estimation where early stages are held fixed while later stages are optimized.
#' @param weights Optional. Name of a variable in the data containing observation weights.
#'   When specified, each observation's contribution to the log-likelihood is multiplied
#'   by its weight. Useful for survey weights, importance sampling, or giving different
#'   observations different influence on the estimation. Weights should be positive.
#'   The variable is extracted from the data passed to \code{estimate_model_rcpp()}.
#'
#' @return An object of class "model_system". A list of model_component objects and one factor_model object.
#' @export
define_model_system <- function(components, factor, previous_stage = NULL, weights = NULL) {
  # Validate the inputs:

  if (!is.list(components) || !all(sapply(components, inherits, "model_component"))) {
    stop("Input must be a list of model_component objects.")
  }

  if (is.null(names(components))){
    stop("model_component objects must be named")
  }

  if (!inherits(factor, "factor_model")) {
    stop("`factor` must be of class 'factor_model'")
  }


  if (!all(sapply(components, function(x) inherits(x, "model_component")))) {
    stop("All elements must be of class 'model_component'")
  }

  #check if each component has the same factor model as factor

 if (!all(vapply(components, function(comp) identical(get_factor(comp), factor), logical(1)))) {
   stop("All model_components must have Factor as the factor_model")
 }

  # Validate dynamic components
  for (comp in components) {
    if (isTRUE(comp$is_dynamic)) {
      if (is.null(comp$outcome_factor) ||
          comp$outcome_factor < 1L ||
          comp$outcome_factor > factor$n_factors) {
        stop(sprintf("Dynamic component '%s' has invalid outcome_factor (%s). Must be between 1 and %d.",
                     comp$name, comp$outcome_factor, factor$n_factors))
      }
    }
  }

  # Validate weights parameter
  if (!is.null(weights)) {
    if (!is.character(weights) || length(weights) != 1) {
      stop("`weights` must be a single character string (variable name in data)")
    }
  }

  # Handle previous_stage if provided
  previous_stage_info <- NULL
  if (!is.null(previous_stage)) {
    # Validate previous_stage
    if (!is.list(previous_stage) ||
        !all(c("model_system", "estimates", "std_errors") %in% names(previous_stage))) {
      stop("`previous_stage` must be a result object from estimate_model_rcpp() with model_system, estimates, and std_errors")
    }

    prev_ms <- previous_stage$model_system
    if (!inherits(prev_ms, "model_system")) {
      stop("`previous_stage$model_system` must be a model_system object")
    }

    # Check factor models match
    if (!identical(prev_ms$factor, factor)) {
      stop("Factor model in previous_stage must be identical to current factor model")
    }

    # Mark all previous stage components as having fixed parameters
    prev_components <- prev_ms$components
    for (i in seq_along(prev_components)) {
      prev_components[[i]]$all_params_fixed <- TRUE
    }

    # Prepend previous stage components to new components
    components <- c(prev_components, components)

    # Store metadata about previous stage
    previous_stage_info <- list(
      n_components = length(prev_components),
      fixed_param_values = previous_stage$estimates,
      fixed_param_names = names(previous_stage$estimates),
      fixed_std_errors = previous_stage$std_errors,
      n_params_fixed = length(previous_stage$estimates)
    )

    # Mark factor variance as fixed
    factor$variance_fixed <- TRUE
    factor$variance_value <- previous_stage$estimates[1]  # First param is always factor variance
  }

  out <- list(
    components = components,
    factor = factor,
    previous_stage_info = previous_stage_info,
    weights = weights)
  class(out) <- "model_system"
  return(out)
}

#' @export
print.model_system <- function(x, ...) {

  stopifnot(inherits(x, "model_system")) #check that x is a class model_system
  comps <- x$components
  n <- length(comps)
  nms <- names(comps) #captures the list names (check if there are errors in doing it this way)

  cat("Model System\n")
  cat("------------\n")
  cat("Components:", n, "\n")
  if (!is.null(x$weights)) {
    cat("Observation weights:", x$weights, "\n")
  }

  if (!n) return(invisible(x))

  for (i in seq_along(comps)) { #loops over indices, seq_along safer for empty vectors
    # header for this component
    label <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else paste0("<unnamed-", i, ">") #get label for component, otherwise call it unnamed - number
    cat("\n[", i, "] ", label, "\n", sep = "") #print section header
    cat(strrep("-", 2 + nchar(label) + nchar(as.character(i))), "\n", sep = "") #print underline

    # delegate to the component's own print method
    print(comps[[i]], ...)   # passes ... along to the model_component's own print method
  }

  invisible(x) #if 0 components, stop printing and return x invisibly
}

