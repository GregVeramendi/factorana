#' Define latent factor model structure
#'
#' @param n_factors Integer. Number of latent factors (>=1)
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

  if (!is.numeric(n_factors) || n_factors < 1) stop("n_factors must be a positive integer.")
  if (!is.numeric(n_types) || n_types < 1) stop("n_types must be a positive integer")
  if (!is.numeric(n_quad_points) || n_quad_points < 1) stop("n_quad_points must be a positive integer.")
  if (!is.logical(correlation)) stop("correlation must be either TRUE or FALSE")
  if (!n_mixtures %in% 1:3) stop("n_mixtures should be between 1-3") #currently this is the case, might change later?

#instead of is.numeric should i use is.integer,,, why do i have to convert to integer after?

  #part of code in C++ to copy in R:
  ###  f_nvariance = nfac;
  ###  if (fac_corr!=0) f_nvariance = (nfac*(nfac+1)/2);      correlation -> n_factors*(n_factors+1)/2, no-correlation: variances only -> n_factors
  ###  nfac_param = f_nvariance*fac_nmix + (fac_nmix-1)*nfac + (fac_nmix-1);

  # f_nvariance = number of unique variance-covariance params per mixture:
  f_nvariance <- if (isTRUE(correlation)) n_factors * (n_factors + 1L) / 2L else n_factors

  nfac_param <- as.integer(f_nvariance * n_mixtures +
                             (n_mixtures - 1L) * n_factors +
                             (n_mixtures - 1L))

  out <- list(
    n_factors = as.integer(n_factors),
    n_types = as.integer(n_types),
    n_quad_points = as.integer(n_quad_points),
    correlation = correlation,
    n_mixtures = as.integer(n_mixtures),
    nfac_param = nfac_param, #replaced with lines 36-38
    params = rep(0.0, nfac_param)
  )

  class(out) <- "factor_model"
  return(out)
}

#calculate nfac_param in C++
#// nmix = 1,2,3 implemented (nmix=1 --> no mixture)
#// fac_nmix var-covar matrices, (nmix-1)*nfac means, nmix-1 weights
#if ( (fac_nmix==1) || (fac_nmix==2) || (fac_nmix==3)) {
#  // nmix variance, nmix-1 mean, and nmix-1 weight parameters
###  f_nvariance = nfac;
###  if (fac_corr!=0) f_nvariance = (nfac*(nfac+1)/2);
###  nfac_param = f_nvariance*fac_nmix + (fac_nmix-1)*nfac + (fac_nmix-1);
#  nparam = nfac_param;
#  parconstrained.assign(nfac_param,-1);
#}


#' @export
print.factor_model <- function(x, ...) {
  cat("Factor Model\n")
  cat("------------\n")
  cat("Number of latent factors:", x$n_factors, "\n")
  cat("Number of types:", x$n_types, "\n")
  cat("Number of quadrature points:        ", x$n_quad_points, "\n")
  cat("Correlation allowed?:     ", x$correlation, "\n")
  cat("Number of mixtures:       ", x$n_mixtures, "\n")
  cat("Number of parameters in latent factor distribution:", x$nparam)
}

