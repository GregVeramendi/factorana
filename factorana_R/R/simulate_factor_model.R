# Simulate data from a fitted (or fully specified) factorana model. Follows the
# legacy TMinLkhd::Simulate / TModel::Sim algorithm:
#   1. resample a base row from `data` (with weights) to supply the covariates X,
#   2. draw the latent factor(s) from the unconditional mixture-of-normals,
#   3. draw the latent type from the type model,
#   4. for each component, form V = X*beta + loadings*f (+ type effects) and draw
#      the outcome by model type (linear / probit / logit / ordered probit).
#
# Output is "full detail" like the legacy simulation_data.csv: the resampled
# covariates, the drawn factors and type, and per component the linear predictor,
# the shock, the choice probabilities, and the simulated outcome.

# ---- small helpers --------------------------------------------------------

# Loading on factor k for a single-outcome component block: a fixed value from
# loading_normalization, otherwise the free parameter <comp>_loading_<k>.
.sim_loading_value <- function(params, comp, k) {
  ln <- comp$loading_normalization[k]
  if (!is.na(ln)) return(ln)
  .sim_param(params, sprintf("%s_loading_%d", comp$name, k))
}

.sim_param <- function(params, nm) {
  if (!nm %in% names(params)) {
    stop("simulate_factor_model: parameter '", nm, "' not found in `params`.")
  }
  unname(params[nm])
}

# Beta for one covariate. Falls back to a fixed coefficient stored on the
# component (fix_coefficient()) when the parameter is absent from `params`.
# `choice` is the multinomial-logit alternative index (NULL for single-outcome).
.sim_beta_value <- function(params, comp, cov, choice = NULL) {
  nm <- if (is.null(choice)) sprintf("%s_%s", comp$name, cov)
        else sprintf("%s_c%d_%s", comp$name, choice, cov)
  if (nm %in% names(params)) return(unname(params[nm]))
  for (fc in comp$fixed_coefficients) {
    if (identical(fc$covariate, cov) &&
        ((is.null(choice) && is.null(fc$choice)) || identical(fc$choice, choice))) {
      return(fc$value)
    }
  }
  stop("simulate_factor_model: parameter '", nm, "' not found in `params` and no ",
       "fixed coefficient for covariate '", cov, "' on component '", comp$name, "'.")
}

# Linear predictor X*beta for a single-outcome component (intercept named
# <comp>_intercept, other covariates <comp>_<cov>).
.sim_xbeta_single <- function(params, comp, Xrow) {
  v <- 0.0
  for (cov in comp$covariates) {
    v <- v + .sim_beta_value(params, comp, cov) * Xrow[[cov]]
  }
  v
}

# Factor contribution sum_k loading_k * f_k for a single-outcome component.
.sim_facsum_single <- function(params, comp, fac) {
  s <- 0.0
  for (k in seq_along(fac)) {
    s <- s + .sim_loading_value(params, comp, k) * fac[k]
  }
  s
}

# Draw the latent factors for all n units. Returns an n x nfac matrix. Handles
# independent, correlated (2-factor), mixture-of-normals, and structural-
# equation (SE_*) factor structures.
.sim_draw_factors <- function(fm, params, n) {
  nfac <- fm$n_factors
  if (nfac == 0) return(matrix(0.0, n, 0))
  structure <- fm$factor_structure
  nmix <- if (is.null(fm$n_mixtures)) 1L else fm$n_mixtures

  # ---- structural-equation structures: draw input factors + residual ----
  if (structure %in% c("SE_linear", "SE_quadratic", "SE_interactions", "SE_full")) {
    n_input <- nfac - 1L
    fmat <- matrix(0.0, n, nfac)
    for (k in seq_len(n_input)) {
      sd_k <- sqrt(.sim_param(params, sprintf("factor_var_%d", k)))
      fmat[, k] <- sd_k * stats::rnorm(n)
    }
    se_int <- .sim_param(params, "se_intercept")
    res_sd <- sqrt(.sim_param(params, "se_residual_var"))
    fout <- se_int + res_sd * stats::rnorm(n)
    for (j in seq_len(n_input)) {
      fout <- fout + .sim_param(params, sprintf("se_linear_%d", j)) * fmat[, j]
    }
    if (structure %in% c("SE_quadratic", "SE_full")) {
      for (j in seq_len(n_input)) {
        fout <- fout + .sim_param(params, sprintf("se_quadratic_%d", j)) * fmat[, j]^2
      }
    }
    if (structure %in% c("SE_interactions", "SE_full") && n_input >= 2L) {
      for (a in seq_len(n_input - 1L)) for (b in (a + 1L):n_input) {
        fout <- fout + .sim_param(params, sprintf("se_interaction_%d_%d", a, b)) *
          fmat[, a] * fmat[, b]
      }
    }
    fmat[, nfac] <- fout
    return(fmat)
  }

  # ---- correlated 2-factor structure (Cholesky) ----
  if (structure == "correlation") {
    if (nfac != 2L) stop("simulate_factor_model: 'correlation' supports exactly 2 factors.")
    sd1 <- sqrt(.sim_param(params, "factor_var_1"))
    sd2 <- sqrt(.sim_param(params, "factor_var_2"))
    rho <- .sim_param(params, "factor_corr_1_2")
    z1 <- stats::rnorm(n); z2 <- stats::rnorm(n)
    fmat <- cbind(sd1 * z1, sd2 * (rho * z1 + sqrt(1 - rho^2) * z2))
    return(fmat)
  }

  # ---- independent factors, possibly a mixture of normals ----
  if (nmix > 1L) {
    # mixture weights via the legacy logit scheme; last mixture is the reference
    lw <- vapply(seq_len(nmix - 1L), function(m)
      .sim_param(params, sprintf("mix%d_logweight", m)), numeric(1))
    den <- 1 + sum(exp(lw))
    w <- c(exp(lw) / den, 1 / den)
    # per-mixture means; the last mixture's mean keeps the overall mean at 0
    means <- matrix(0.0, nmix, nfac)
    for (m in seq_len(nmix - 1L)) for (k in seq_len(nfac)) {
      means[m, k] <- .sim_param(params, sprintf("mix%d_factor_mean_%d", m, k))
    }
    for (k in seq_len(nfac)) {
      means[nmix, k] <- -sum(exp(lw) * means[seq_len(nmix - 1L), k])
    }
    sds <- matrix(0.0, nmix, nfac)
    for (m in seq_len(nmix)) for (k in seq_len(nfac)) {
      sds[m, k] <- sqrt(.sim_param(params, sprintf("mix%d_factor_var_%d", m, k)))
    }
    comp_draw <- sample.int(nmix, n, replace = TRUE, prob = w)
    fmat <- matrix(0.0, n, nfac)
    for (k in seq_len(nfac)) {
      fmat[, k] <- means[comp_draw, k] + sds[comp_draw, k] * stats::rnorm(n)
    }
    attr(fmat, "mixture") <- comp_draw
    return(fmat)
  }

  # plain independent normals
  fmat <- matrix(0.0, n, nfac)
  for (k in seq_len(nfac)) {
    sd_k <- sqrt(.sim_param(params, sprintf("factor_var_%d", k)))
    fmat[, k] <- sd_k * stats::rnorm(n)
  }
  fmat
}

# Draw a latent type (1..n_types) per unit from the type model: a multinomial
# logit on the factors with utility U_t = typeprob_t_intercept +
# sum_k type_t_loading_k * f_k for t = 2..n_types (type 1 is the reference).
.sim_draw_types <- function(fm, params, fmat, n) {
  ntyp <- if (is.null(fm$n_types)) 1L else fm$n_types
  if (ntyp <= 1L) return(rep(1L, n))
  nfac <- ncol(fmat)
  U <- matrix(0.0, n, ntyp)  # column t = type t; type 1 reference stays 0
  for (t in 2:ntyp) {
    u <- rep(.sim_param(params, sprintf("typeprob_%d_intercept", t)), n)
    for (k in seq_len(nfac)) {
      u <- u + .sim_param(params, sprintf("type_%d_loading_%d", t, k)) * fmat[, k]
    }
    U[, t] <- u
  }
  ex <- exp(U - apply(U, 1, max))
  prob <- ex / rowSums(ex)
  apply(prob, 1, function(p) sample.int(ntyp, 1, prob = p))
}

#' Simulate data from a factorana model
#'
#' Generates a simulated data set from a fitted or fully specified factorana
#' model, following the legacy simulation algorithm: for each simulated unit a
#' base row is drawn from `data` (with weights) to supply the covariates, the
#' latent factor(s) are drawn from the model's unconditional distribution, and
#' each component's outcome is drawn from its model (linear, probit, logit, or
#' ordered probit) given the covariates and factors.
#'
#' Supports all four model types (linear, probit, logit including multinomial,
#' ordered probit), latent types, mixtures of normals, and the independent,
#' correlated, and structural-equation (`SE_*`) factor structures. The simulated
#' outcome column is named after each component's outcome variable, so the
#' returned data frame can be re-estimated with the same model system.
#'
#' @param object A `factorana_result` from [estimate_model_rcpp()], or a
#'   `model_system`. If a `model_system`, supply `params`.
#' @param data Data frame whose rows are resampled to supply the covariates.
#' @param n Number of units to simulate (default: `nrow(data)`).
#' @param params Named parameter vector. Defaults to `object$estimates` when
#'   `object` is a fitted result.
#' @param weights Optional resampling weights: a column name in `data` or a
#'   numeric vector of length `nrow(data)`.
#' @param seed Optional integer seed.
#' @param detail If `TRUE` (default), include per-component linear predictor
#'   (`Vobs`), factor contributions (`Vfac<k>`), shock (`eps`), and choice
#'   probabilities, like the legacy detailed output.
#' @return A data frame with one row per simulated unit: an `xid`, the resampled
#'   covariates, the drawn latent factor(s) (`factor_<k>`), and per component the
#'   simulated outcome (named after the component) plus, when `detail = TRUE`,
#'   the diagnostic columns.
#' @seealso [estimate_model_rcpp()]
#' @export
simulate_factor_model <- function(object, data, n = NULL, params = NULL,
                                  weights = NULL, seed = NULL, detail = TRUE) {
  if (inherits(object, "factorana_result")) {
    ms <- object$model_system
    if (is.null(params)) params <- object$estimates
  } else if (inherits(object, "model_system")) {
    ms <- object
    if (is.null(params)) stop("`params` is required when `object` is a model_system.")
  } else {
    stop("`object` must be a 'factorana_result' or a 'model_system'.")
  }
  if (!is.null(seed)) set.seed(seed)
  if (is.null(n)) n <- nrow(data)

  fm <- ms$factor
  nfac <- fm$n_factors

  # ---- resample base rows for covariates ----
  if (is.null(weights)) {
    prob <- NULL
  } else {
    wv <- if (length(weights) == 1L && is.character(weights)) data[[weights]] else weights
    if (length(wv) != nrow(data)) stop("`weights` must have length nrow(data).")
    if (anyNA(wv) || any(wv < 0)) stop("`weights` must be non-negative and non-missing.")
    prob <- wv / sum(wv)
  }
  base_idx <- sample.int(nrow(data), n, replace = TRUE, prob = prob)
  base <- data[base_idx, , drop = FALSE]

  # ---- draw factors and latent type ----
  fmat <- .sim_draw_factors(fm, params, n)
  types <- .sim_draw_types(fm, params, fmat, n)

  out <- data.frame(xid = base_idx)
  # carry covariates actually used by any component
  used_cov <- unique(unlist(lapply(ms$components, function(c) c$covariates)))
  for (cov in used_cov) out[[cov]] <- base[[cov]]
  if (nfac > 0) for (k in seq_len(nfac)) out[[sprintf("factor_%d", k)]] <- fmat[, k]
  if (!is.null(fm$n_types) && fm$n_types > 1) out[["type"]] <- types
  if (!is.null(attr(fmat, "mixture"))) out[["mixture"]] <- attr(fmat, "mixture")

  # ---- simulate each component ----
  for (comp in ms$components) {
    sim <- .simulate_component(comp, params, base, fmat, types, detail)
    for (nm in names(sim)) out[[nm]] <- sim[[nm]]
  }
  out
}

# Type-specific intercept per unit for a component that uses types (0 for the
# reference type 1 or when the component does not use types). For multinomial
# logit, pass the per-alternative prefix (e.g. "ml_c1") as `name_pref`.
.sim_type_intercept <- function(params, comp, types, name_pref = comp$name) {
  ti <- numeric(length(types))
  if (!isTRUE(comp$use_types)) return(ti)
  for (t in unique(types[types >= 2L])) {
    ti[types == t] <- .sim_param(params, sprintf("%s_type_%d_intercept", name_pref, t))
  }
  ti
}

# Simulate one component for all n units. Returns a named list of columns.
.simulate_component <- function(comp, params, base, fmat, types, detail) {
  n <- nrow(base)
  nfac <- ncol(fmat)
  mt <- comp$model_type
  cn <- comp$name
  oc <- comp$outcome
  if (length(oc) != 1L) {
    stop("simulate_factor_model: component '", cn, "' has a multi-column outcome ",
         "(exploded/ranked logit), which is not yet supported.")
  }
  res <- list()

  # Missing indicator from the evaluation_indicator (rows where the component
  # is not evaluated get NA outcomes).
  miss <- rep(0L, n)
  if (!is.null(comp$evaluation_indicator) && comp$evaluation_indicator %in% colnames(base)) {
    ev <- base[[comp$evaluation_indicator]]
    miss <- as.integer(!(!is.na(ev) & (ev == 1L | ev == TRUE)))
  }
  res[[paste0(cn, "_miss")]] <- miss

  # type-specific intercept per unit (0 unless the component uses types)
  ti <- .sim_type_intercept(params, comp, types)

  if (mt == "linear") {
    V <- vapply(seq_len(n), function(i)
      .sim_xbeta_single(params, comp, base[i, , drop = FALSE]) +
      .sim_facsum_single(params, comp, fmat[i, ]), numeric(1)) + ti
    sigma <- abs(.sim_param(params, sprintf("%s_sigma", cn)))
    eps <- stats::rnorm(n, 0, sigma)
    y <- V + eps
    if (detail) { res[[paste0(cn, "_Vobs")]] <- V; res[[paste0(cn, "_eps")]] <- eps }
    res[[oc]] <- y

  } else if (mt == "probit") {
    V <- vapply(seq_len(n), function(i)
      .sim_xbeta_single(params, comp, base[i, , drop = FALSE]) +
      .sim_facsum_single(params, comp, fmat[i, ]), numeric(1)) + ti
    eps <- stats::rnorm(n)
    prob <- stats::pnorm(V)
    y <- as.integer((V + eps) > 0)
    if (detail) { res[[paste0(cn, "_Vobs")]] <- V; res[[paste0(cn, "_eps")]] <- eps
                  res[[paste0(cn, "_prob")]] <- prob }
    res[[oc]] <- y

  } else if (mt == "oprobit") {
    V <- vapply(seq_len(n), function(i)
      .sim_xbeta_single(params, comp, base[i, , drop = FALSE]) +
      .sim_facsum_single(params, comp, fmat[i, ]), numeric(1)) + ti
    K <- comp$num_choices
    # cumulative thresholds tau_1..tau_{K-1} (incremental, abs after the first)
    incr <- vapply(seq_len(K - 1L), function(j)
      .sim_param(params, sprintf("%s_thresh_%d", cn, j)), numeric(1))
    tau <- cumsum(c(incr[1], abs(incr[-1])))
    eps <- stats::rnorm(n)
    latent <- V + eps
    y <- findInterval(latent, tau) + 1L  # category in 1..K
    if (detail) { res[[paste0(cn, "_Vobs")]] <- V; res[[paste0(cn, "_eps")]] <- eps }
    res[[oc]] <- y

  } else if (mt == "logit") {
    K <- comp$num_choices
    # choice-specific linear predictors V_a for a = 2..K (reference choice 1 = 0).
    # Parameter names follow the fitted-estimate convention: binary logit (K = 2)
    # uses the single-block form <comp>_<cov> / <comp>_loading_<k>; multinomial
    # logit uses a per-alternative prefix <comp>_c<a>_... for a = 1..K-1.
    Vc <- matrix(0.0, n, K)  # column j = choice j; column 1 (reference) stays 0
    for (a in seq_len(K - 1L)) {
      pref <- if (K == 2L) cn else sprintf("%s_c%d", cn, a)
      ch <- if (K == 2L) NULL else a
      beta_v <- numeric(n)
      for (cov in comp$covariates) {
        beta_v <- beta_v + .sim_beta_value(params, comp, cov, choice = ch) * base[[cov]]
      }
      facsum <- numeric(n)
      for (k in seq_len(nfac)) {
        ln <- comp$loading_normalization[k]
        lv <- if (!is.na(ln)) ln else .sim_param(params, sprintf("%s_loading_%d", pref, k))
        facsum <- facsum + lv * fmat[, k]
      }
      Vc[, a + 1L] <- beta_v + facsum + .sim_type_intercept(params, comp, types, pref)
    }
    # softmax probabilities and a Gumbel-max draw
    ex <- exp(Vc - apply(Vc, 1, max))
    prob <- ex / rowSums(ex)
    g <- -log(-log(matrix(stats::runif(n * K), n, K)))  # standard Gumbel
    y <- max.col(Vc + g, ties.method = "first")
    if (detail) {
      for (a in seq_len(K - 1L)) res[[sprintf("%s_Vobs_alt%d", cn, a)]] <- Vc[, a + 1L]
      for (j in seq_len(K)) res[[sprintf("%s_prob_C%d", cn, j)]] <- prob[, j]
    }
    res[[oc]] <- y

  } else {
    stop("simulate_factor_model: unsupported model_type '", mt, "'.")
  }

  # apply missingness
  if (any(miss != 0L)) res[[oc]][miss != 0L] <- NA
  res
}
