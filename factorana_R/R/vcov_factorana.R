#' Robust and cluster-robust covariance for a fitted factorana model
#'
#' Computes a sandwich (Huber-White) or cluster-robust covariance matrix for a
#' model fitted with [estimate_model_rcpp()], as a post-hoc alternative to the
#' default inverse-Hessian (model-based) standard errors. The estimator is
#' \deqn{V = H^{-1} B H^{-1},}
#' where the "bread" \eqn{H^{-1}} is the inverse observed information already
#' computed at fit time and the "meat" \eqn{B} is built from the
#' per-observation scores of the marginal (factor-integrated) log-likelihood:
#' \eqn{B = \sum_i g_i g_i'} (robust) or \eqn{B = \sum_c s_c s_c'} with
#' \eqn{s_c = \sum_{i \in c} g_i} (cluster).
#'
#' @details
#' One factorana observation is one row of `data`, with its measurements
#' integrated jointly over the latent factors. The latent factor therefore
#' already models the within-row dependence. Clustering is meaningful only when
#' a cluster groups several rows (for example several individuals in the same
#' household or school, or several stacked person-wave rows for one person in a
#' long-format second stage). With one row per cluster, the cluster estimator
#' reduces to the ordinary robust estimator.
#'
#' For two-stage estimates these standard errors treat the first-stage
#' (measurement) parameters and any plug-in factor scores as known, so they
#' understate uncertainty (the generated-regressor problem). Use a bootstrap
#' that re-runs both stages for honest two-stage inference.
#'
#' @param object A `factorana_result` from [estimate_model_rcpp()].
#' @param data The data frame used to fit `object` (same rows, same order).
#' @param type Covariance type: `"hessian"` (inverse information, the default,
#'   reproduces the fit's standard errors), `"robust"` (Huber-White), or
#'   `"cluster"` (cluster-robust).
#' @param cluster For `type = "cluster"`: either the name of a column in `data`
#'   holding the cluster id, or a vector of length `nrow(data)`.
#' @param n_quad Number of quadrature points for recomputing the per-observation
#'   scores. Defaults to the value used at estimation (stored on `object`).
#' @param finite_sample Logical; apply the usual finite-sample scaling
#'   (`n/(n-k)` for robust; `G/(G-1) * (N-1)/(N-k)` for cluster). Default `TRUE`.
#'
#' @return A `nparam x nparam` covariance matrix over the full parameter vector,
#'   with row/column names equal to the parameter names. Fixed and
#'   equality-tied parameters have zero rows and columns. Standard errors are
#'   `sqrt(diag(.))`.
#'
#' @seealso [estimate_model_rcpp()]
#' @export
vcov_factorana <- function(object, data,
                           type = c("hessian", "robust", "cluster"),
                           cluster = NULL, n_quad = NULL,
                           finite_sample = TRUE) {
  type <- match.arg(type)
  if (!inherits(object, "factorana_result")) {
    stop("`object` must be a 'factorana_result' from estimate_model_rcpp().")
  }

  estimates <- object$estimates
  free_idx  <- object$free_idx
  bread     <- object$cov_free
  n_full    <- length(estimates)
  n_free    <- length(free_idx)

  if (is.null(bread) || is.null(free_idx)) {
    stop("This fit does not carry the inverse-Hessian covariance needed for a ",
         "sandwich estimator (standard errors may have failed at fit time, or ",
         "the model was fit by an older version). Refit with the current ",
         "version of factorana.")
  }

  # Place a free-parameter covariance into the full parameter space.
  to_full <- function(V_free) {
    V <- matrix(0.0, n_full, n_full,
                dimnames = list(names(estimates), names(estimates)))
    V[free_idx, free_idx] <- V_free
    V
  }

  if (type == "hessian") {
    return(to_full(bread))
  }

  # ---- robust / cluster: need the per-observation score matrix at the MLE ----
  if (is.null(n_quad)) n_quad <- object$n_quad
  if (is.null(n_quad)) {
    stop("Could not determine `n_quad`; pass the number of quadrature points ",
         "used at estimation.")
  }

  data_mat <- as.matrix(data)
  fm_ptr <- initialize_factor_model_cpp(object$model_system, data_mat,
                                        as.integer(n_quad), estimates)
  pf <- extract_free_params_cpp(fm_ptr, estimates)
  if (length(pf) != n_free) {
    stop("Free-parameter count from the rebuilt model (", length(pf),
         ") does not match the fit (", n_free, "). Make sure `data` is the ",
         "same data frame used to fit `object`.")
  }

  G <- evaluate_obs_scores_cpp(fm_ptr, pf)  # nobs x n_free, marginal scores
  nobs <- nrow(G)

  if (type == "robust") {
    meat <- crossprod(G)  # sum_i g_i g_i'
    if (finite_sample) meat <- meat * (nobs / (nobs - n_free))
  } else {  # cluster
    if (is.null(cluster)) {
      stop("type = 'cluster' requires `cluster`: a column name in `data` or a ",
           "vector of length nrow(data).")
    }
    if (length(cluster) == 1L && is.character(cluster)) {
      if (!cluster %in% colnames(data)) {
        stop("cluster column '", cluster, "' not found in `data`.")
      }
      cl <- data[[cluster]]
    } else {
      cl <- cluster
    }
    if (length(cl) != nobs) {
      stop("cluster length (", length(cl), ") must equal the number of ",
           "observations (", nobs, ").")
    }
    if (anyNA(cl)) stop("`cluster` must not contain NA.")

    S <- rowsum(G, group = cl)  # n_clusters x n_free
    n_clusters <- nrow(S)
    if (n_clusters < 2L) {
      stop("Cluster-robust covariance needs at least 2 clusters; found ",
           n_clusters, ".")
    }
    meat <- crossprod(S)  # sum_c s_c s_c'
    if (finite_sample) {
      meat <- meat * (n_clusters / (n_clusters - 1)) *
        ((nobs - 1) / (nobs - n_free))
    }
    if (n_clusters < 30L) {
      warning("Only ", n_clusters, " clusters: cluster-robust standard errors ",
              "can be unreliable with few clusters; consider a (wild) cluster ",
              "bootstrap.", call. = FALSE)
    }
  }

  V_free <- bread %*% meat %*% bread
  V_free <- (V_free + t(V_free)) / 2  # enforce symmetry against round-off
  to_full(V_free)
}

#' Robust / cluster-robust standard errors for a fitted factorana model
#'
#' Convenience wrapper returning `sqrt(diag(vcov_factorana(...)))` as a named
#' vector over the full parameter vector. See [vcov_factorana()] for details
#' and caveats.
#'
#' @inheritParams vcov_factorana
#' @return Named numeric vector of standard errors (fixed/tied parameters are 0).
#' @seealso [vcov_factorana()]
#' @export
robust_se <- function(object, data,
                      type = c("robust", "cluster"),
                      cluster = NULL, n_quad = NULL, finite_sample = TRUE) {
  type <- match.arg(type)
  V <- vcov_factorana(object, data, type = type, cluster = cluster,
                      n_quad = n_quad, finite_sample = finite_sample)
  sqrt(pmax(0, diag(V)))
}
