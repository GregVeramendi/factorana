# Simulate data from a fitted (or fully specified) factorana model. Follows the
# legacy TMinLkhd::Simulate / TModel::Sim algorithm:
#   1. resample a base row from `data` (with weights) to supply the covariates X,
#   2. draw the latent factor(s) from the unconditional mixture-of-normals,
#   3. draw the latent type from the type model,
#   4. for each component, form V = X*beta + loadings*f (+ type effects) and draw
#      the outcome by model type (linear / probit / logit / ordered probit). A
#      ranked (exploded) logit component draws one alternative per rank column,
#      mirroring the rank loop of Model::EvalLogit() in src/Model.cpp.
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

# ---- endogenous-regressor resolution -------------------------------------

# Resolve the declared endogenous-regressor links among the model components.
# A component's `endogenous` field names OTHER components' outcomes whose
# simulated values feed this component's covariates during simulation. Returns:
#   deps   : list; deps[[i]] = integer vector of producer indices component i needs
#   writes : named list keyed by PRODUCER OUTCOME name; each element is the
#            character vector of columns to (re)write into the data after that
#            producer runs -- the raw outcome name and/or its one-hot dummies that
#            are actually used as covariates downstream.
#   needs  : list; needs[[i]] = the covariates of component i that its declared
#            producers write. Only these read simulated values, as in the legacy
#            TModel::SetEndogenousRegs(), where a regressor is endogenous only in
#            the models that declare it.
#   miss_cols : list; miss_cols[[i]] = covariates of component i named
#            `<producer outcome>_miss` for a declared producer. The legacy
#            simulation sets their contribution to zero.
#   producer : named list; producer[[outcome]] = index of the component that
#            produces that outcome.
# `bins` is the `endogenous_bins` argument of simulate_factor_model(): thresholds
# per linear producer outcome, whose consumers read bin dummies.
# Errors on a declared name that matches no component outcome, a self-loop, or a
# declared name that drives none of the consumer's covariates.
.sim_resolve_endogenous <- function(components, bins = NULL) {
  n <- length(components)
  deps <- vector("list", n)
  needs <- vector("list", n)      # per consumer: covariates replaced by simulated values
  miss_cols <- vector("list", n)  # per consumer: `<producer outcome>_miss` covariates

  # map each single-column outcome name -> producing component index, and each
  # rank column of a ranked (exploded) logit -> its component index
  outcome_to_index <- list()
  ranked_to_index <- list()
  for (i in seq_len(n)) {
    oc <- components[[i]]$outcome
    # when several components share an outcome variable, the first one listed is
    # the producer, as in the legacy TModel::SetEndogenousRegs()
    if (length(oc) == 1L) { if (is.null(outcome_to_index[[oc]])) outcome_to_index[[oc]] <- i }
    else for (col in oc) ranked_to_index[[col]] <- i
  }

  writes <- list()
  for (i in seq_len(n)) {
    endo <- components[[i]]$endogenous
    if (is.null(endo) || length(endo) == 0L) next
    cov_i <- components[[i]]$covariates
    for (nm in unique(endo)) {
      p <- outcome_to_index[[nm]]
      if (is.null(p)) {
        rp <- ranked_to_index[[nm]]
        if (!is.null(rp)) {
          stop(sprintf(paste0("simulate_factor_model: component '%s' declares endogenous ",
                              "regressor '%s', which is a rank column of the ranked (exploded) ",
                              "logit component '%s'. A ranked outcome cannot be an endogenous ",
                              "regressor in simulation."),
                       components[[i]]$name, nm, components[[rp]]$name))
        }
        stop(sprintf(paste0("simulate_factor_model: component '%s' declares endogenous ",
                            "regressor '%s', but no component produces an outcome named '%s'."),
                     components[[i]]$name, nm, nm))
      }
      if (p == i) {
        stop(sprintf("simulate_factor_model: component '%s' cannot be endogenous on its own outcome '%s'.",
                     components[[i]]$name, nm))
      }
      P <- components[[p]]
      if (P$model_type == "linear" && !is.null(bins[[nm]])) {
        # binned linear producer: bin dummies nm1..nmK (K = thresholds + 1), or raw `nm`
        K <- length(bins[[nm]]) + 1L
        needed <- intersect(c(nm, paste0(nm, seq_len(K))), cov_i)
      } else if (P$model_type %in% c("linear", "probit")) {
        # continuous / binary producer: substitute the raw value into covariate `nm`
        needed <- intersect(nm, cov_i)
      } else {
        # discrete producer (logit / oprobit): one-hot dummies nm1..nmK, or raw `nm`
        K <- P$num_choices
        needed <- intersect(c(nm, paste0(nm, seq_len(K))), cov_i)
      }
      if (length(needed) == 0L) {
        stop(sprintf(paste0("simulate_factor_model: component '%s' declares endogenous ",
                            "regressor '%s' (produced by '%s'), but none of its covariates ",
                            "match '%s' or its dummies '%s1', '%s2', ..."),
                     components[[i]]$name, nm, P$name, nm, nm, nm))
      }
      deps[[i]] <- unique(c(deps[[i]], p))
      writes[[nm]] <- union(writes[[nm]], needed)
      needs[[i]] <- union(needs[[i]], needed)
      mc <- paste0(nm, "_miss")
      if (mc %in% cov_i) miss_cols[[i]] <- union(miss_cols[[i]], mc)
    }
  }
  list(deps = deps, writes = writes, needs = needs, miss_cols = miss_cols,
       producer = outcome_to_index)
}

# Topologically order component indices so every producer precedes its consumers.
# `deps[[i]]` lists the producer indices that component i depends on (edges p -> i).
# Returns an integer permutation of 1:n; when there are no endogenous links this is
# the identity 1:n (preserving the current behavior). Errors on a cyclic dependency.
.sim_topo_order <- function(n, deps) {
  if (n == 0L) return(integer(0))
  indeg <- vapply(seq_len(n), function(i) length(deps[[i]]), integer(1))
  done <- logical(n)
  ord <- integer(0)
  repeat {
    ready <- which(indeg == 0L & !done)
    if (length(ready) == 0L) break
    v <- ready[1L]            # smallest index -> stable, respects listing order
    done[v] <- TRUE
    ord <- c(ord, v)
    for (i in seq_len(n)) {
      if (!done[i] && (v %in% deps[[i]])) indeg[i] <- indeg[i] - 1L
    }
  }
  if (length(ord) < n) {
    cyc <- which(!done)
    stop("simulate_factor_model: endogenous regressors form a cyclic dependency among ",
         "components ", paste(cyc, collapse = ", "),
         ". A component's covariate cannot (directly or transitively) depend on its own outcome.")
  }
  ord
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
#' Supports all four model types (linear, probit, logit including multinomial
#' and ranked (exploded) logit, ordered probit), latent types, mixtures of
#' normals, and the independent, correlated, and structural-equation (`SE_*`)
#' factor structures. The simulated outcome column is named after each
#' component's outcome variable, so the returned data frame can be re-estimated
#' with the same model system.
#'
#' @section Ranked (exploded) logit:
#' A component with a multi-column `outcome` is simulated rank by rank, as the
#' likelihood evaluates it. The utility of a non-reference alternative is
#' `X beta_a + lambda_a f` plus any type intercept, and at rank `r` the
#' rank-share correction for that rank and alternative is added when
#' `rankshare_var` is set. One alternative is drawn from those available and
#' written to the rank's outcome column. With `exclude_chosen = TRUE` the drawn
#' alternative leaves the choice set for later ranks, and with
#' `exclude_chosen = FALSE` every alternative stays available. The rank-share
#' block, the `(num_choices - 1) * nrank` columns that follow `rankshare_var`
#' in `data`, is read by column position, as the likelihood reads it, with
#' `NA` treated as 0. The block is copied into the
#' simulated data in the same order, so the data can be re-estimated with the
#' same model system. A rank column cannot be an endogenous regressor of
#' another component.
#'
#' `fill_ranks = "first"` reproduces the legacy `TModel::Sim()`, which has no
#' rank loop: one multinomial draw over all alternatives, written to the first
#' outcome column, with the other rank columns left `NA` so that a
#' re-estimation uses rank 1 only. The utilities follow the likelihood at rank
#' 1, so a rank-share correction for rank 1 is applied when `rankshare_var` is
#' set, whereas the legacy simulation ignores it.
#'
#' @section Endogenous regressors:
#' A covariate reads a simulated value only in the components that declare its
#' producer through the `endogenous` argument of [define_model_component()].
#' Every other component reads the resampled data value, as in the legacy
#' `TModel::SetEndogenousRegs()`. The returned standard-layout data frame holds
#' the simulated value in the column.
#'
#' @section Legacy layout:
#' `layout = "legacy"` returns the table that the legacy C++
#' `TMinLkhd::Simulate()` wrote to `simulation_data.csv`, so that code written
#' for that file runs on the returned object. The table is a `data.table` when
#' the data.table package is installed and a data frame otherwise. Missing
#' values are `NA`, in memory and in the file written by `csv_file`, where the
#' C++ file has `-9999`.
#' Columns, in order:
#' \itemize{
#'   \item the index: `index_var` when given, otherwise `xid`, the row number of
#'     the resampled base row.
#'   \item when the model has latent types, the type-probability model under
#'     the name `type_model_name`, laid out as a multinomial logit component.
#'   \item one block per component, in simulation order (see
#'     `component_order`). A block holds
#'     `sim_<name>_miss`, the detail columns when `detail` covers the component,
#'     the choice probabilities of a discrete component, and the simulated
#'     outcome `sim_<name>`.
#'   \item the drawn factors `sim_fac0`, `sim_fac1`, ..., the type indicators
#'     `sim_type1`, ... for the non-reference types, and `sim_flag` (always 1).
#' }
#' The detail columns split the index of the outcome equation as the C++ code
#' does: `_Vobs` (covariates that are not endogenous), `_Vend` (endogenous
#' covariates, present only when the component declares some), `_Vfac<k>` (the
#' loading times factor `k`, numbered from 0), `_Vtyp<t>` (the type intercept
#' of non-reference type `t + 2` times its indicator, numbered from 0, present
#' whenever the model has types), and `_eps` (the shock). The outcome can be
#' rebuilt from these pieces. A multinomial logit with more than two
#' alternatives repeats each piece with the suffix `_C<j>` for alternatives
#' `j = 2, ..., K`, and its `_eps_C<j>` is the difference between the
#' extreme-value shocks of alternative `j` and alternative 1. Three conventions
#' of the C++ file are kept: `_prob` of a probit is the probability of 1,
#' `_prob_C<j>` of a logit is the probability of alternative `j`, and the
#' `K - 1` columns `_prob_C2`, ..., `_prob_CK` of an ordered probit hold the
#' probabilities of categories `1, ..., K - 1`.
#'
#' The legacy layout simulates every component for every unit
#' (`simulate_unevaluated = TRUE`), as the C++ code does. `sim_<name>_miss` is
#' 0 for a unit whose base row evaluates the component, 1 for a unit whose base
#' row does not (the outcome is simulated all the same), and 2 for a unit with
#' a missing covariate (`NA`) or a missing endogenous producer, whose block is
#' missing. The package reads only `NA` as missing: a covariate that still holds
#' the code -9999 of the legacy data files stops the simulation with a request
#' to recode it. A ranked logit component is simulated as
#' the C++ code simulates it, with one draw over all alternatives
#' (`fill_ranks = "first"`). A covariate named `<outcome>_miss` of a declared
#' producer is set to 0, as in the C++ simulation.
#'
#' With the same `seed`, the legacy layout and the standard layout hold the
#' same draws. The extreme-value shocks of the type-probability model are drawn
#' last, conditional on the drawn type.
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
#' @param detail `TRUE` (default), `FALSE`, or a character vector of component
#'   names. Covered components include their detail columns. In the standard
#'   layout these are the linear predictor (`<comp>_Vobs`, covariates plus
#'   factor and type terms), the shock (`<comp>_eps`), and the choice
#'   probabilities. A ranked logit component adds per rank `r` the columns
#'   `<comp>_r<r>_Vobs_alt<a>` (the utility of non-reference alternative `a` at
#'   that rank, rank-share correction included) and `<comp>_r<r>_prob_C<j>`
#'   (the probability of alternative `j` over the alternatives available at
#'   that rank, 0 for an alternative excluded by an earlier rank). In the
#'   legacy layout the detail columns are those of the section on that layout,
#'   and `type_model_name` may be listed to cover the type-probability model.
#' @param fill_ranks For ranked (exploded) logit components, which rank columns
#'   a simulated unit fills. `"base"` (default) fills the ranks in which the
#'   resampled base row holds a valid choice (a value in `1..num_choices`), so
#'   the simulated data keep the observed distribution of the number of ranked
#'   alternatives. `"all"` fills every rank. `"first"` fills rank 1 only, as the
#'   legacy `TModel::Sim()` does. An unfilled rank is `NA`, which the
#'   likelihood skips. The legacy layout uses `"first"`.
#' @param layout `"standard"` (default) returns the resampled covariates and
#'   the simulated outcomes under the outcome variable names, ready for
#'   re-estimation. `"legacy"` returns the column layout of the legacy C++
#'   `simulation_data.csv`.
#' @param simulate_unevaluated If `FALSE`, the outcome of a component is `NA`
#'   for a unit whose base row does not evaluate it (its
#'   `evaluation_indicator` is off). If `TRUE`, the outcome is simulated for
#'   every unit with complete covariates, as the legacy C++ code does, and
#'   `<comp>_miss` takes the values 0, 1 and 2 described in the section on the
#'   legacy layout. Defaults to `FALSE` in the standard layout and `TRUE` in
#'   the legacy layout.
#' @param component_order Optional character vector of component names. The
#'   named components are simulated first, in that order, and the others
#'   follow in their listing order, as with the legacy `SetSimModelOrder()`.
#'   A component that is not named still waits for the producers of its
#'   endogenous regressors, and a named component cannot precede one of its
#'   producers. Without this argument the components run in listing order,
#'   with producers moved ahead of their consumers. The order is also the
#'   order of the component blocks in the legacy layout.
#' @param endogenous_bins Optional named list of increasing threshold vectors,
#'   named by the outcome of a linear producer component. Consumers that
#'   declare that outcome as endogenous then read the bin dummies
#'   `<outcome>1`, ..., `<outcome>K` for `K = length(thresholds) + 1` bins,
#'   where the bin of a unit is 1 plus the number of thresholds strictly below
#'   its simulated outcome. This is the legacy `LastModel_SetXtileThresh()`. A
#'   covariate named `<outcome>` still reads the simulated value itself.
#' @param index_var Legacy layout only: name of a column of `data` that
#'   identifies the base row (the C++ index variable). It becomes the first
#'   column of the table, under its own name.
#' @param type_model_name Legacy layout only: name under which the
#'   type-probability model is written. The C++ code names it `"TypeProb"`.
#' @param csv_file Legacy layout only: path of a CSV file to write in the
#'   format of the C++ `simulation_data.csv`. See
#'   [write_legacy_simulation_csv()].
#' @param rds_file Path of a compressed R file ([saveRDS()]) holding the
#'   returned table, in either layout. Read it back with [readRDS()].
#' @return Standard layout: a data frame with one row per simulated unit: an
#'   `xid`, the resampled covariates (and, for a ranked logit component with
#'   `rankshare_var`, its rank-share block), the drawn latent factor(s)
#'   (`factor_<k>`), and per component the simulated outcome (named after the
#'   component's outcome variable) plus the detail columns. Legacy layout: the
#'   table described in the section on that layout.
#' @seealso [estimate_model_rcpp()], [write_legacy_simulation_csv()],
#'   [legacy_param_table()]
#' @export
simulate_factor_model <- function(object, data, n = NULL, params = NULL,
                                  weights = NULL, seed = NULL, detail = TRUE,
                                  fill_ranks = c("base", "all", "first"),
                                  layout = c("standard", "legacy"),
                                  simulate_unevaluated = NULL, component_order = NULL,
                                  endogenous_bins = NULL,
                                  index_var = NULL, type_model_name = "TypeProb",
                                  csv_file = NULL, rds_file = NULL) {
  layout <- match.arg(layout)
  legacy <- layout == "legacy"
  if (legacy) {
    if (!missing(fill_ranks) && !identical(fill_ranks, "first")) {
      stop("simulate_factor_model: layout = \"legacy\" writes one simulated choice per ",
           "logit component, as the legacy TModel::Sim() does. Leave `fill_ranks` unset ",
           "or use fill_ranks = \"first\".")
    }
    fill_ranks <- "first"
  } else {
    fill_ranks <- match.arg(fill_ranks)
    if (!is.null(csv_file)) {
      stop("simulate_factor_model: `csv_file` writes the legacy C++ file format and ",
           "needs layout = \"legacy\". Use `rds_file` to save the standard layout.")
    }
  }
  if (is.null(simulate_unevaluated)) simulate_unevaluated <- legacy
  if (!isTRUE(simulate_unevaluated) && !isFALSE(simulate_unevaluated)) {
    stop("simulate_factor_model: `simulate_unevaluated` must be TRUE or FALSE.")
  }

  if (inherits(object, "factorana_result")) {
    ms <- object$model_system
    if (is.null(params)) params <- object$estimates
  } else if (inherits(object, "model_system")) {
    ms <- object
    if (is.null(params)) stop("`params` is required when `object` is a model_system.")
  } else {
    stop("`object` must be a 'factorana_result' or a 'model_system'.")
  }

  comp_names <- vapply(ms$components, function(c) c$name, character(1))
  if (is.character(detail)) {
    unknown <- setdiff(detail, c(comp_names, type_model_name))
    if (length(unknown) > 0L) {
      stop("simulate_factor_model: `detail` names no component: ",
           paste(sQuote(unknown), collapse = ", "), ".")
    }
    detail_on <- stats::setNames(comp_names %in% detail, comp_names)
    detail_types <- type_model_name %in% detail
  } else {
    detail_on <- stats::setNames(rep(isTRUE(detail), length(comp_names)), comp_names)
    detail_types <- isTRUE(detail)
  }
  if (!is.null(index_var) && !(length(index_var) == 1L && index_var %in% names(data))) {
    stop("simulate_factor_model: `index_var` must name one column of `data`.")
  }

  if (!is.null(seed)) set.seed(seed)
  if (is.null(n)) n <- nrow(data)

  fm <- ms$factor
  nfac <- fm$n_factors
  ntyp <- if (is.null(fm$n_types)) 1L else as.integer(fm$n_types)
  # A model declared with types that no component uses carries no type
  # parameters (the estimation collapses it to one type), so it has no types to draw.
  if (ntyp > 1L &&
      !any(vapply(ms$components, function(c) isTRUE(c$use_types), logical(1))) &&
      !("typeprob_2_intercept" %in% names(params))) {
    ntyp <- 1L
  }

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
  # Only the columns the simulation reads are resampled, in the order of `data`,
  # so the rank-share block of a ranked logit component stays contiguous.
  keep <- .sim_needed_columns(ms, data, index_var)
  base <- lapply(as.list(data)[keep], function(col) col[base_idx])
  attr(base, "row.names") <- .set_row_names(n)
  class(base) <- "data.frame"

  # ---- draw factors and latent type ----
  fmat <- .sim_draw_factors(fm, params, n)
  types <- if (ntyp > 1L) .sim_draw_types(fm, params, fmat, n) else rep(1L, n)

  # ---- resolve endogenous-regressor links and order producers before consumers ----
  if (!is.null(endogenous_bins)) {
    lin <- unlist(lapply(ms$components, function(c)
      if (c$model_type == "linear" && length(c$outcome) == 1L) c$outcome))
    ok <- is.list(endogenous_bins) && !is.null(names(endogenous_bins)) &&
      all(names(endogenous_bins) %in% lin) &&
      all(vapply(endogenous_bins, function(b)
        is.numeric(b) && length(b) >= 1L && !anyNA(b) && !is.unsorted(b, strictly = TRUE), logical(1)))
    if (!ok) {
      stop("simulate_factor_model: `endogenous_bins` must be a list of increasing numeric ",
           "threshold vectors named by the outcomes of linear components.")
    }
  }
  endo <- .sim_resolve_endogenous(ms$components, endogenous_bins)
  sim_order <- .sim_topo_order(length(ms$components), endo$deps)
  if (!is.null(component_order)) {
    sim_order <- .sim_explicit_order(component_order, comp_names, endo$deps)
  }

  if (!legacy) {
    out <- data.frame(xid = base_idx)
    # carry covariates actually used by any component (endogenous dummy columns may
    # be absent from `data` -- they are materialized when their producer is simulated)
    used_cov <- unique(unlist(lapply(ms$components, function(c) c$covariates)))
    for (cov in used_cov) if (cov %in% names(base)) out[[cov]] <- base[[cov]]
    # carry the rank-share block of every ranked logit component, in data order and
    # contiguous, because the likelihood addresses the block by column position
    for (comp in ms$components) {
      if (length(comp$outcome) > 1L && !is.null(comp$rankshare_var)) {
        rs_cols <- .sim_rankshare_block(comp, base)
        for (col in rs_cols) out[[col]] <- base[[col]]
        pos <- match(rs_cols, names(out))
        if (any(diff(pos) != 1L)) {
          stop("simulate_factor_model: the rank-share columns of component '", comp$name,
               "' are not contiguous in the simulated data (one of them is also a ",
               "covariate or an earlier column). Rename the block so it stands alone.")
        }
      }
    }
    if (nfac > 0) for (k in seq_len(nfac)) out[[sprintf("factor_%d", k)]] <- fmat[, k]
    if (!is.null(fm$n_types) && fm$n_types > 1 && ntyp > 1L) out[["type"]] <- types
    if (!is.null(attr(fmat, "mixture"))) out[["mixture"]] <- attr(fmat, "mixture")
  }

  # ---- simulate each component (producers before consumers) ----
  simcols <- list()   # simulated values of endogenous covariates, by column name
  blocks <- list()    # legacy layout: the columns of each component
  for (ci in sim_order) {
    comp <- ms$components[[ci]]
    # The covariate frame of this component: the resampled data, with the
    # covariates its declared producers write replaced by their simulated values.
    Xi <- base
    for (col in endo$needs[[ci]]) Xi[[col]] <- simcols[[col]]
    if (simulate_unevaluated) for (col in endo$miss_cols[[ci]]) Xi[[col]] <- rep(0, n)
    endo_cols <- as.character(c(endo$needs[[ci]], endo$miss_cols[[ci]]))

    sim <- .simulate_component(
      comp, params, Xi, fmat, types,
      detail = !legacy && detail_on[[comp$name]],
      fill_ranks = fill_ranks, simulate_unevaluated = simulate_unevaluated,
      endo_cols = endo_cols, ntyp = ntyp,
      pieces = if (!legacy) "none" else if (detail_on[[comp$name]]) "full" else "basic")
    if (legacy) {
      # built here so that the decomposition of one component at a time is in memory
      blocks[[comp$name]] <- .sim_legacy_block(comp$name, attr(sim, "pieces"),
                                               detail_on[[comp$name]], nfac, ntyp)
      attr(sim, "pieces") <- NULL
    } else {
      for (nm in names(sim)) out[[nm]] <- sim[[nm]]
    }

    # Keep this component's simulated outcome for the consumers that declared it
    # endogenous, and write it into the returned standard-layout frame.
    oc <- comp$outcome
    if (length(oc) == 1L && !is.null(endo$writes[[oc]]) && identical(endo$producer[[oc]], ci)) {
      y <- sim[[oc]]
      # category of the outcome: the outcome itself for a discrete producer, the
      # threshold bin for a linear producer listed in `endogenous_bins`
      ycat <- if (is.null(endogenous_bins[[oc]])) y
              else findInterval(y, endogenous_bins[[oc]], left.open = TRUE) + 1L
      for (col in endo$writes[[oc]]) {
        if (identical(col, oc)) {
          val <- y                                   # raw value (continuous/binary/category)
        } else {
          j <- as.integer(substring(col, nchar(oc) + 1L))
          val <- as.integer(ycat == j)               # one-hot dummy for category j
        }
        simcols[[col]] <- val
        if (!legacy) out[[col]] <- val
      }
    }
  }

  if (legacy) {
    index <- if (is.null(index_var)) base_idx else base[[index_var]]
    out <- .sim_legacy_table(
      blocks = blocks[comp_names[sim_order]], fmat = fmat, types = types, ntyp = ntyp,
      index = index, index_name = if (is.null(index_var)) "xid" else index_var,
      type_block = if (ntyp > 1L) {
        .sim_legacy_type_block(type_model_name, params, fmat, types, ntyp, detail_types)
      } else NULL)
    rm(blocks)
    if (requireNamespace("data.table", quietly = TRUE)) data.table::setDT(out)
    if (!is.null(csv_file)) write_legacy_simulation_csv(out, csv_file)
  }
  if (!is.null(rds_file)) saveRDS(out, rds_file)
  out
}

# Names of the columns of `data` that the simulation reads, in the order of
# `data`: covariates, evaluation indicators, the rank columns and the rank-share
# block of ranked logit components, and the index variable.
.sim_needed_columns <- function(ms, data, index_var = NULL) {
  need <- index_var
  for (comp in ms$components) {
    need <- c(need, comp$covariates, comp$evaluation_indicator)
    if (length(comp$outcome) > 1L) {
      need <- c(need, comp$outcome)
      if (!is.null(comp$rankshare_var)) need <- c(need, .sim_rankshare_block(comp, data))
    }
  }
  names(data)[names(data) %in% need]
}

# Simulation order given by the user, as the legacy SetSimModelOrder(): the
# components named in `component_order` come first, in that order, and the
# others follow in listing order, a consumer waiting for its producers. A named
# component cannot precede one of its producers.
.sim_explicit_order <- function(component_order, comp_names, deps) {
  unknown <- setdiff(component_order, comp_names)
  if (length(unknown) > 0L || anyDuplicated(component_order) > 0L) {
    stop("simulate_factor_model: `component_order` must list distinct component names",
         if (length(unknown) > 0L) paste0(". Unknown: ", paste(sQuote(unknown), collapse = ", ")),
         ".")
  }
  n <- length(comp_names)
  done <- rep(FALSE, n)
  ord <- integer(0)
  for (i in match(component_order, comp_names)) {
    late <- deps[[i]][!done[deps[[i]]]]
    if (length(late) > 0L) {
      stop("simulate_factor_model: `component_order` simulates component '", comp_names[i],
           "' before its endogenous producer(s) ", paste(sQuote(comp_names[late]), collapse = ", "),
           ". List the producers first.")
    }
    done[i] <- TRUE
    ord <- c(ord, i)
  }
  repeat {
    ready <- which(!done & vapply(seq_len(n), function(i) all(done[deps[[i]]]), logical(1)))
    if (length(ready) == 0L) break
    done[ready[1L]] <- TRUE          # smallest index: the listing order
    ord <- c(ord, ready[1L])
  }
  ord
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

# The type intercept split by type: a list with one vector per non-reference
# type t = 2..ntyp, holding the intercept of type t for the units of that type
# and 0 elsewhere (all zeros when the component does not use types). The
# vectors add up to .sim_type_intercept().
.sim_type_parts <- function(params, comp, types, ntyp, name_pref = comp$name) {
  parts <- vector("list", max(ntyp - 1L, 0L))
  for (t in seq_len(ntyp)[-1L]) {
    v <- numeric(length(types))
    if (isTRUE(comp$use_types)) {
      sel <- types == t
      if (any(sel)) v[sel] <- .sim_param(params, sprintf("%s_type_%d_intercept", name_pref, t))
    }
    parts[[t - 1L]] <- v
  }
  parts
}

# Simulate one component for all n units. Returns a named list of columns.
#
# `base` is the covariate frame of this component (endogenous covariates already
# hold simulated values) and `endo_cols` names those covariates. With
# simulate_unevaluated = FALSE the outcome is NA for a unit whose evaluation
# indicator is off. With simulate_unevaluated = TRUE the outcome is simulated
# for every unit, and `<comp>_miss` is 2, and the outcome NA, for a unit with a
# missing covariate, as the legacy TModel::Sim() does for its -9999 code. Only NA
# is missing here: a covariate at or below -9998 is rejected, since it would enter
# the index as a number.
#
# `pieces` asks for the decomposition the legacy layout writes, returned in the
# attribute "pieces": "basic" keeps the missing code, the legacy probabilities
# and the outcome, and "full" adds the index split into covariates (Vx), the
# endogenous covariates among them (Vend), each factor term (Vfac), each type
# term (Vtyp), and the shock (eps), one column per non-reference alternative.
.simulate_component <- function(comp, params, base, fmat, types, detail,
                                fill_ranks = "base", simulate_unevaluated = FALSE,
                                endo_cols = character(0), ntyp = 1L, pieces = "none") {
  n <- nrow(base)
  nfac <- ncol(fmat)
  mt <- comp$model_type
  cn <- comp$name
  oc <- comp$outcome
  res <- list()

  absent <- setdiff(comp$covariates, names(base))
  if (length(absent) > 0L) {
    stop("simulate_factor_model: covariate(s) ", paste(sQuote(absent), collapse = ", "),
         " of component '", cn, "' not found in `data`.")
  }

  # Missing indicator from the evaluation_indicator.
  miss <- rep(0L, n)
  if (!is.null(comp$evaluation_indicator) && comp$evaluation_indicator %in% colnames(base)) {
    ev <- base[[comp$evaluation_indicator]]
    miss <- as.integer(!(!is.na(ev) & (ev == 1L | ev == TRUE)))
  }
  # Units that cannot be simulated because a covariate is missing.
  blocked <- rep(FALSE, n)
  if (simulate_unevaluated) {
    for (cov in comp$covariates) {
      x <- base[[cov]]
      if (is.numeric(x) && any(x <= -9998, na.rm = TRUE)) {
        stop("simulate_factor_model: covariate '", cov, "' of component '", cn,
             "' holds values at or below -9998. A missing value must be NA: the package ",
             "does not read the -9999 code of the legacy data files as missing. ",
             "Recode those values to NA before simulating.")
      }
      blocked <- blocked | is.na(x)
    }
    miss[blocked] <- 2L
  }
  res[[paste0(cn, "_miss")]] <- miss
  # units whose outcome is not reported
  drop <- if (simulate_unevaluated) blocked else miss != 0L
  has_endog <- length(endo_cols) > 0L

  # ranked (exploded) logit: one draw per rank column
  if (length(oc) > 1L) {
    if (mt != "logit") {
      stop("simulate_factor_model: component '", cn, "' has a multi-column outcome ",
           "but model_type '", mt, "'. Only logit supports ranked outcomes.")
    }
    rk <- .sim_ranked_logit(comp, params, base, fmat, types, detail, drop, fill_ranks,
                            endo_cols, ntyp, pieces)
    pc <- attr(rk, "pieces")
    res <- c(res, rk)
    if (!is.null(pc)) {
      pc$miss <- miss
      pc$has_endog <- has_endog
      attr(res, "pieces") <- pc
    }
    return(res)
  }

  pc <- NULL
  if (mt %in% c("linear", "probit", "oprobit")) {
    # index V = X*beta + sum_k loading_k * f_k + type intercept
    xb <- numeric(n)
    for (cov in comp$covariates) {
      xb <- xb + .sim_beta_value(params, comp, cov) * base[[cov]]
    }
    vfac <- vector("list", nfac)
    facsum <- numeric(n)
    for (k in seq_len(nfac)) {
      vfac[[k]] <- .sim_loading_value(params, comp, k) * fmat[, k]
      facsum <- facsum + vfac[[k]]
    }
    V <- (xb + facsum) + .sim_type_intercept(params, comp, types)
    K <- 2L
    prob_legacy <- NULL

    if (mt == "linear") {
      sigma <- abs(.sim_param(params, sprintf("%s_sigma", cn)))
      eps <- stats::rnorm(n, 0, sigma)
      y <- V + eps
      if (detail) { res[[paste0(cn, "_Vobs")]] <- V; res[[paste0(cn, "_eps")]] <- eps }

    } else if (mt == "probit") {
      eps <- stats::rnorm(n)
      prob <- stats::pnorm(V)
      y <- as.integer((V + eps) > 0)
      if (detail) { res[[paste0(cn, "_Vobs")]] <- V; res[[paste0(cn, "_eps")]] <- eps
                    res[[paste0(cn, "_prob")]] <- prob }
      prob_legacy <- matrix(prob, n, 1L)

    } else {
      K <- comp$num_choices
      # cumulative thresholds tau_1..tau_{K-1} (incremental, abs after the first)
      incr <- vapply(seq_len(K - 1L), function(j)
        .sim_param(params, sprintf("%s_thresh_%d", cn, j)), numeric(1))
      tau <- cumsum(c(incr[1], abs(incr[-1])))
      eps <- stats::rnorm(n)
      latent <- V + eps
      y <- findInterval(latent, tau) + 1L  # category in 1..K
      if (detail) { res[[paste0(cn, "_Vobs")]] <- V; res[[paste0(cn, "_eps")]] <- eps }
      if (pieces != "none") {
        # probabilities of categories 1..K-1, the K-1 values the legacy file holds
        cdf <- matrix(vapply(tau, function(tj) stats::pnorm(tj - V), numeric(n)), n, K - 1L)
        prob_legacy <- cdf
        if (K > 2L) prob_legacy[, -1L] <- cdf[, -1L] - cdf[, -(K - 1L)]
      }
    }
    res[[oc]] <- y

    if (pieces != "none") {
      pc <- list(type = mt, K = K, prob = prob_legacy)
      if (pieces == "full") {
        vend <- numeric(n)
        for (cov in intersect(comp$covariates, endo_cols)) {
          vend <- vend + .sim_beta_value(params, comp, cov) * base[[cov]]
        }
        pc$Vx <- matrix(xb, n, 1L)
        pc$Vend <- matrix(vend, n, 1L)
        pc$Vfac <- lapply(vfac, function(v) matrix(v, n, 1L))
        pc$Vtyp <- lapply(.sim_type_parts(params, comp, types, ntyp),
                          function(v) matrix(v, n, 1L))
        pc$eps <- matrix(eps, n, 1L)
      }
    }

  } else if (mt == "logit") {
    K <- comp$num_choices
    # choice-specific linear predictors V_a for a = 2..K (reference choice 1 = 0)
    Vc <- .sim_logit_utilities(params, comp, base, fmat, types, endo_cols, ntyp,
                               parts = pieces == "full")
    parts <- attr(Vc, "parts")
    attr(Vc, "parts") <- NULL
    # softmax probabilities and a Gumbel-max draw
    ex <- exp(Vc - .sim_rowmax(Vc))
    prob <- ex / rowSums(ex)
    g <- -log(-log(matrix(stats::runif(n * K), n, K)))  # standard Gumbel
    y <- max.col(Vc + g, ties.method = "first")
    if (detail) {
      for (a in seq_len(K - 1L)) res[[sprintf("%s_Vobs_alt%d", cn, a)]] <- Vc[, a + 1L]
      for (j in seq_len(K)) res[[sprintf("%s_prob_C%d", cn, j)]] <- prob[, j]
    }
    res[[oc]] <- y
    if (pieces != "none") {
      pc <- list(type = "logit", K = K, prob = prob[, -1L, drop = FALSE])
      if (pieces == "full") {
        pc <- c(pc, parts)
        pc$eps <- g[, -1L, drop = FALSE] - g[, 1L]
      }
    }

  } else {
    stop("simulate_factor_model: unsupported model_type '", mt, "'.")
  }

  # apply missingness
  if (any(drop)) res[[oc]][drop] <- NA
  if (!is.null(pc)) {
    pc$miss <- miss
    pc$has_endog <- has_endog
    pc$y <- res[[oc]]
    attr(res, "pieces") <- pc
  }
  res
}

# Row-wise maximum of a numeric matrix without a per-row function call.
.sim_rowmax <- function(M) {
  m <- M[, 1L]
  for (j in seq_len(ncol(M))[-1L]) m <- pmax(m, M[, j])
  m
}

# Rank-invariant multinomial logit utilities: an n x K matrix whose column a
# holds V_a = X*beta_a + sum_k lambda_{a,k} f_k + type intercept for the
# non-reference alternatives a = 2..K, and 0 for the reference alternative in
# column 1. Parameter names follow the fitted-estimate convention: binary logit
# (K = 2) uses the single-block form <comp>_<cov> / <comp>_loading_<k>, and
# multinomial logit uses a per-alternative prefix <comp>_c<a>_... for a = 1..K-1.
# With parts = TRUE the attribute "parts" holds the same utilities split into
# n x (K-1) matrices: Vx (covariates), Vend (the covariates in `endo_cols`),
# Vfac (a list with one matrix per factor) and Vtyp (one per non-reference type).
.sim_logit_utilities <- function(params, comp, base, fmat, types,
                                 endo_cols = character(0), ntyp = 1L, parts = FALSE) {
  n <- nrow(base)
  nfac <- ncol(fmat)
  K <- comp$num_choices
  cn <- comp$name
  Vc <- matrix(0.0, n, K)
  if (parts) {
    Vx <- matrix(0.0, n, K - 1L)
    Vend <- matrix(0.0, n, K - 1L)
    Vfac <- lapply(seq_len(nfac), function(k) matrix(0.0, n, K - 1L))
    Vtyp <- lapply(seq_len(max(ntyp - 1L, 0L)), function(t) matrix(0.0, n, K - 1L))
  }
  for (a in seq_len(K - 1L)) {
    pref <- if (K == 2L) cn else sprintf("%s_c%d", cn, a)
    ch <- if (K == 2L) NULL else a
    beta_v <- numeric(n)
    vend <- numeric(n)
    for (cov in comp$covariates) {
      term <- .sim_beta_value(params, comp, cov, choice = ch) * base[[cov]]
      beta_v <- beta_v + term
      if (parts && cov %in% endo_cols) vend <- vend + term
    }
    facsum <- numeric(n)
    for (k in seq_len(nfac)) {
      ln <- comp$loading_normalization[k]
      lv <- if (!is.na(ln)) ln else .sim_param(params, sprintf("%s_loading_%d", pref, k))
      fk <- lv * fmat[, k]
      facsum <- facsum + fk
      if (parts) Vfac[[k]][, a] <- fk
    }
    Vc[, a + 1L] <- beta_v + facsum + .sim_type_intercept(params, comp, types, pref)
    if (parts) {
      Vx[, a] <- beta_v
      Vend[, a] <- vend
      tp <- .sim_type_parts(params, comp, types, ntyp, pref)
      for (t in seq_along(tp)) Vtyp[[t]][, a] <- tp[[t]]
    }
  }
  if (parts) attr(Vc, "parts") <- list(Vx = Vx, Vend = Vend, Vfac = Vfac, Vtyp = Vtyp)
  Vc
}

# Column names of the rank-share block of a ranked logit component: the
# (K-1) * nrank columns of `base` that start at `rankshare_var`, in data order.
# The likelihood addresses the block by column position (rcpp_interface.cpp finds
# `rankshare_var` by name and Model::EvalLogit() reads position + (K-1)*irank +
# icat), so the simulation reads the same positions. Returns character(0) when
# the component has no rank-share variable.
.sim_rankshare_block <- function(comp, base) {
  rs <- comp$rankshare_var
  if (is.null(rs)) return(character(0))
  K <- comp$num_choices
  nrank <- length(comp$outcome)
  j0 <- match(rs, names(base))
  if (is.na(j0)) {
    stop("simulate_factor_model: rankshare_var '", rs, "' of component '", comp$name,
         "' not found in `data`.")
  }
  nblock <- (K - 1L) * nrank
  if (j0 + nblock - 1L > ncol(base)) {
    stop(sprintf(paste0("simulate_factor_model: component '%s' needs %d rank-share columns ",
                        "starting at '%s' ((num_choices - 1) * nrank = %d * %d), but only ",
                        "%d columns of `data` start there."),
                 comp$name, nblock, rs, K - 1L, nrank, ncol(base) - j0 + 1L))
  }
  names(base)[j0 + seq_len(nblock) - 1L]
}

# Simulate a ranked (exploded) multinomial logit component, mirroring the rank
# loop of Model::EvalLogit() in src/Model.cpp. At rank r the utility of
# alternative a is V_a + corr_{r,a}: V_a is the rank-invariant index from
# .sim_logit_utilities() and corr_{r,a} the rank-share correction read from the
# data block (0 for the reference alternative and for NA). One alternative is drawn from those available with a Gumbel-max draw.
# With exclude_chosen = TRUE an alternative drawn at an earlier rank leaves the
# choice set, and with exclude_chosen = FALSE every alternative stays available.
# A unit fills rank r when the likelihood would evaluate that rank: under
# fill_ranks = "base" when the resampled base row holds a valid choice (1..K) in
# that rank column, under fill_ranks = "all" at every rank, and under
# fill_ranks = "first" at rank 1 only (the legacy TModel::Sim() draw: one
# multinomial choice over all alternatives, no later ranks). A unit in `drop`
# (logical) fills no rank. An unfilled rank is NA, which the likelihood skips,
# and it leaves the choice set unchanged. Detail columns per drawn rank:
# <comp>_r<r>_Vobs_alt<a> for the non-reference alternatives and
# <comp>_r<r>_prob_C<j> over the available set (0 for an excluded alternative,
# NA for a unit with no alternative left). `pieces` returns the decomposition of
# the rank-1 draw for the legacy layout (see .simulate_component()).
.sim_ranked_logit <- function(comp, params, base, fmat, types, detail, drop, fill_ranks,
                              endo_cols = character(0), ntyp = 1L, pieces = "none") {
  n <- nrow(base)
  cn <- comp$name
  oc <- comp$outcome
  K <- comp$num_choices
  nrank <- length(oc)
  excl <- isTRUE(comp$exclude_chosen)
  res <- list()

  # which ranks each unit fills
  if (fill_ranks == "all") {
    fill <- matrix(TRUE, n, nrank)
  } else if (fill_ranks == "first") {
    fill <- matrix(FALSE, n, nrank)
    fill[, 1L] <- TRUE
  } else {
    absent <- setdiff(oc, names(base))
    if (length(absent) > 0L) {
      stop("simulate_factor_model: fill_ranks = \"base\" needs the rank columns of ",
           "component '", cn, "' in `data`, but ", paste(sQuote(absent), collapse = ", "),
           " are absent. Add them or use fill_ranks = \"all\" or \"first\".")
    }
    fill <- matrix(FALSE, n, nrank)
    for (r in seq_len(nrank)) {
      v <- suppressWarnings(as.numeric(base[[oc[r]]]))
      fill[, r] <- !is.na(v) & v >= 1 & v < K + 1
    }
  }
  fill[drop, ] <- FALSE

  Vc <- .sim_logit_utilities(params, comp, base, fmat, types, endo_cols, ntyp,
                             parts = pieces == "full")
  parts <- attr(Vc, "parts")
  attr(Vc, "parts") <- NULL
  rs_cols <- .sim_rankshare_block(comp, base)
  pc <- NULL

  last_rank <- if (fill_ranks == "first") 1L else nrank
  avail <- matrix(TRUE, n, K)
  for (r in seq_len(last_rank)) {
    V <- Vc
    if (length(rs_cols) > 0L) {
      for (a in 2:K) {
        corr <- as.numeric(base[[rs_cols[(K - 1L) * (r - 1L) + (a - 1L)]]])
        corr[is.na(corr)] <- 0
        V[, a] <- V[, a] + corr
        # the legacy decomposition counts the correction with the covariates
        if (r == 1L && !is.null(parts)) parts$Vx[, a - 1L] <- parts$Vx[, a - 1L] + corr
      }
    }
    Vav <- V
    Vav[!avail] <- -Inf
    has_choice <- rowSums(avail) > 0L
    draw <- fill[, r] & has_choice
    g <- -log(-log(matrix(stats::runif(n * K), n, K)))  # standard Gumbel
    y <- max.col(Vav + g, ties.method = "first")
    y[!draw] <- NA_integer_
    want_pieces <- pieces != "none" && r == 1L
    if (detail || want_pieces) {
      m <- .sim_rowmax(Vav)
      ex <- exp(Vav - m)                       # 0 for excluded alternatives
      prob <- ex / rowSums(ex)
      prob[!has_choice, ] <- NA_real_
    }
    if (detail) {
      for (a in seq_len(K - 1L)) {
        res[[sprintf("%s_r%d_Vobs_alt%d", cn, r, a)]] <- V[, a + 1L]
      }
      for (j in seq_len(K)) res[[sprintf("%s_r%d_prob_C%d", cn, r, j)]] <- prob[, j]
    }
    if (want_pieces) {
      pc <- list(type = "logit", K = K, prob = prob[, -1L, drop = FALSE], y = y)
      if (pieces == "full") {
        pc <- c(pc, parts)
        pc$eps <- g[, -1L, drop = FALSE] - g[, 1L]
      }
    }
    if (excl && any(draw)) avail[cbind(which(draw), y[draw])] <- FALSE
    res[[oc[r]]] <- y
  }
  for (r in seq_len(nrank)[-seq_len(last_rank)]) res[[oc[r]]] <- rep(NA_integer_, n)
  if (!is.null(pc)) attr(res, "pieces") <- pc
  res
}
