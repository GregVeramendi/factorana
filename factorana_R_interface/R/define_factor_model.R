#' Define latent factor model structure
#'
#' Creates an object of class `"factor_model"` that specifies the structure of the
#' unobserved latent factors. This includes the number of factors, mixture components,
#' quadrature points for numerical integration, and any fixed loading constraints.
#'
#' @param n_factors Integer. Number of latent factors (>=1)
#' @param n_types Integer. Number of types (>=1)
#' @param n_quad_points Integer. Number of Gauss-Hermite quadrature points (>=1)
#' @param correlation Logical. Whether to allow correlation between two factors (default = FALSE)
#' @param n_mixtures Integer. Number of discrete mixtures (default = 1, allowed: 1-3)
#' @param loading_normalization Numeric vector of length `n_factors`.
#'   Specifies which loadings are fixed or free:
#'   - `NA` → loading is free (estimated later).
#'   - numeric value → loading is fixed at that value (e.g. `1` to normalize a factor).
#'
#' @return An object of class "factor_model"
#' @export
define_factor_model <- function(n_factors,
                                n_types,
                                n_quad_points,
                                correlation = FALSE,
                                n_mixtures = 1,
                                loading_normalization = NULL) {

  # ---- 1. Input validation ----
  # Check all arguments are the correct type and within supported range.

  if (!is.numeric(n_factors) || n_factors < 1) stop("n_factors must be a positive integer.")
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

  k <- as.integer(n_factors)

  # ---- 4. Handle loading normalization vector ----
  # loading_normalization controls which loadings are fixed/free.
  # Default: all NA = free parameters to be estimated.

  if (is.null(loading_normalization)) {
    loading_normalization <- rep(NA_real_, k)
  } else {
    stopifnot(is.numeric(loading_normalization),
              length(loading_normalization) == k)
  }

  # ---- 5. Construct the factor_model object ----
  # Bundle all parameters and metadata into a list, and assign class.

  out <- list(
    n_factors = as.integer(n_factors),
    n_types = as.integer(n_types),
    n_quad_points = as.integer(n_quad_points),
    correlation = correlation,
    n_mixtures = as.integer(n_mixtures),
    nfac_param = nfac_param, #replaced with lines 36-38
    params = rep(0.0, nfac_param),
    loading_normalization = as.numeric(loading_normalization)
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
  cat("Number of quadrature points:        ", x$n_quad_points, "\n")
  cat("Correlation allowed?:     ", x$correlation, "\n")
  cat("Number of mixtures:       ", x$n_mixtures, "\n")
  cat("Number of parameters in latent factor distribution:", x$nfac_param, "\n")
  cat("Loading Normalization:", x$loading_normalization)
}

