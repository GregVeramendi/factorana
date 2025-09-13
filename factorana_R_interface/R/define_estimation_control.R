#' Define estimation control settings
#'
#' @param model_system Object of class model_system.
#' @param num_cores Integer. Number of processes to use for parallel estimation (default = 1)
#'
#' @return An object of class estimation_control containing control settings
#' @export
define_estimation_control <- function(num_cores = 1) {
  if (!is.numeric(num_cores) || num_cores < 1) {
    stop("num_cores must be a positive integer.")
  }

  out <- list(
    num_cores = as.integer(num_cores)
  )

  class(out) <- "estimation_control"
  return(out)
}


#' ###should be a function that calculates and outputs how many parameters
#' #placeholder below
#'
#' #' @keywords internal
#' component_param_count <- function(mc) {
#'   # If you already computed/stored it, trust that:
#'   if (!is.null(mc$nparam_model)) return(as.integer(mc$nparam_model))
#'
#'   # Fallback: derive something sensible (adjust if your spec changes)
#'   k_cov  <- if (!is.null(mc$covariates)) length(mc$covariates) else 0L
#'   k_int  <- if (isTRUE(mc$intercept)) 1L else 0L
#'   k_fac  <- if (!is.null(mc$factor$n_factors)) as.integer(mc$factor$n_factors) else 0L
#'   k_types_minus_1 <- if (!is.null(mc$factor$n_types)) max(as.integer(mc$factor$n_types) - 1L, 0L) else 0L
#'   as.integer(k_cov + k_int + k_fac + k_types_minus_1)
#' }
#'
#' #' @keywords internal
#' param_counts <- function(model_system, factor_model) {
#'   stopifnot(inherits(model_system, "model_system"),
#'             inherits(factor_model, "factor_model"))
#'
#'   comps <- model_system$components
#'   per_component <- vapply(comps, component_param_count, integer(1))
#'   names(per_component) <- names(comps)
#'
#'   total_component_params <- sum(per_component)
#'   factor_dist_params     <- as.integer(factor_model$nfac_param)
#'   total_params           <- as.integer(total_component_params + factor_dist_params)
#'
#'   list(
#'     per_component            = per_component,                 # named int vector
#'     total_component_params   = total_component_params,        # int
#'     factor_dist_params       = factor_dist_params,            # int
#'     total_params             = total_params                   # int
#'   )
#' }





#' @export
print.estimation_control <- function(x, ...) {
  cat("Estimation Control\n")
  cat("------------------\n")
#  cat("Number of unobserved types:", x$num_types, "\n")
  cat("Number of cores for parallelization:", x$num_cores, "\n")
}
