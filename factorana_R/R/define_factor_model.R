#' Define latent factor model structure
#'
#' Creates an object of class `"factor_model"` that specifies the structure of the
#' unobserved latent factors. This includes the number of factors, mixture components,
#' and quadrature points for numerical integration. Loading constraints are now
#' specified at the component level via `define_model_component()`.
#'
#' @param n_factors Integer. Number of latent factors (>=0). Use 0 for models without latent factors.
#' @param n_types Integer. Number of types (>=1)
#' @param n_quad_points Integer. Number of Gauss-Hermite quadrature points (>=1)
#' @param correlation Logical. Whether to allow correlation between two factors (default = FALSE)
#' @param n_mixtures Integer. Number of discrete mixtures (default = 1, allowed: 1-3)
#'
#' @return An object of class "factor_model"
#' @export
define_factor_model <- function(n_factors,
                                n_types,
                                n_quad_points,
                                correlation = FALSE,
                                n_mixtures = 1) {

  # ---- 1. Input validation ----
  # Check all arguments are the correct type and within supported range.

  if (!is.numeric(n_factors) || n_factors < 0) stop("n_factors must be a non-negative integer.")
  if (!is.numeric(n_types) || n_types < 1) stop("n_types must be a positive integer")
  if (!is.numeric(n_quad_points) || n_quad_points < 1) stop("n_quad_points must be a positive integer.")
  if (!is.logical(correlation)) stop("correlation must be either TRUE or FALSE")
  if (!n_mixtures %in% 1:3) stop("n_mixtures should be between 1-3") #currently this is the case, might change later

  # ---- 2. Compute number of variance/covariance parameters ----
  # If correlation = TRUE, include all unique covariance terms among k factors.

  f_nvariance <- if (isTRUE(correlation)) n_factors * (n_factors + 1L) / 2L else n_factors

  # ---- 3. Compute total number of parameters for factor distribution ----
  # Combines variances/covariances with mixture-related parameters.

  nfac_param <- as.integer(f_nvariance * n_mixtures +
                             (n_mixtures - 1L) * n_factors +
                             (n_mixtures - 1L))

  # ---- 4. Construct the factor_model object ----
  # Bundle all parameters and metadata into a list, and assign class.

  out <- list(
    n_factors = as.integer(n_factors),
    n_types = as.integer(n_types),
    n_quad_points = as.integer(n_quad_points),
    correlation = correlation,
    n_mixtures = as.integer(n_mixtures),
    nfac_param = nfac_param,
    params = rep(0.0, nfac_param)
  )

  class(out) <- "factor_model"
  return(out)
}


#' @export
print.factor_model <- function(x, ...) {
  cat("Factor Model\n")
  cat("------------\n")
  cat("Number of latent factors:", x$n_factors, "\n")
  cat("Number of types:", x$n_types, "\n")
  cat("Number of quadrature points:", x$n_quad_points, "\n")
  cat("Correlation allowed?:", x$correlation, "\n")
  cat("Number of mixtures:", x$n_mixtures, "\n")
  cat("Number of parameters in latent factor distribution:", x$nfac_param, "\n")
  invisible(x)
}

