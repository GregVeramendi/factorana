# The legacy layout of simulate_factor_model(): the table the C++
# TMinLkhd::Simulate() / TModel::Sim() wrote to simulation_data.csv
# (libraries/TMinLkhd.cc and libraries/TModel.cc), its file writer, and the
# parameter table that the files fullmodel_par_varlist.txt and fullmodel_par.txt
# gave when merged.

# ---- assembling the table --------------------------------------------------

# The columns of one component, in the order TMinLkhd::Simulate() writes its
# header and TModel::Sim() writes its values:
#   sim_<name>_miss
#   detail only, per non-reference alternative: _Vobs, _Vend (if the component
#     has endogenous covariates), _Vfac0.., _Vtyp0.. (if the model has types and
#     this is not the type model)
#   detail only, per non-reference alternative: _eps
#   discrete components: the probability columns
#   sim_<name>
# A multinomial logit with K > 2 alternatives carries the suffix _C<j>,
# j = 2..K, on every per-alternative column. `p` is the "pieces" list of
# .simulate_component(). A unit with miss == 2 has a missing block.
.sim_legacy_block <- function(name, p, detail, nfac, ntyp, typed = TRUE) {
  pre <- paste0("sim_", name)
  K <- p$K
  multi <- p$type == "logit" && K > 2L
  A <- if (multi) K - 1L else 1L
  sfx <- if (multi) sprintf("_C%d", 2:K) else ""
  blocked <- p$miss == 2L
  blank <- function(x) {
    if (any(blocked)) x[blocked] <- NA
    x
  }
  cols <- list()
  cols[[paste0(pre, "_miss")]] <- p$miss
  if (detail) {
    for (a in seq_len(A)) {
      vobs <- p$Vx[, a]
      if (p$has_endog) vobs <- vobs - p$Vend[, a]
      cols[[paste0(pre, "_Vobs", sfx[a])]] <- blank(vobs)
      if (p$has_endog) cols[[paste0(pre, "_Vend", sfx[a])]] <- blank(p$Vend[, a])
      for (k in seq_len(nfac)) {
        cols[[sprintf("%s_Vfac%d%s", pre, k - 1L, sfx[a])]] <- blank(p$Vfac[[k]][, a])
      }
      if (typed) {
        for (t in seq_len(ntyp - 1L)) {
          cols[[sprintf("%s_Vtyp%d%s", pre, t - 1L, sfx[a])]] <- blank(p$Vtyp[[t]][, a])
        }
      }
    }
    for (a in seq_len(A)) cols[[paste0(pre, "_eps", sfx[a])]] <- blank(p$eps[, a])
  }
  if (!is.null(p$prob)) {
    psfx <- if (K > 2L) sprintf("_C%d", 2:K) else ""
    for (j in seq_len(ncol(p$prob))) {
      cols[[paste0(pre, "_prob", psfx[j])]] <- blank(p$prob[, j])
    }
  }
  cols[[pre]] <- blank(p$y)
  cols
}

# Extreme-value shocks of a multinomial logit conditional on the chosen
# alternative. `U` is the n x J matrix of utilities and `y` the chosen column.
# The maximum of U_j + eps_j is Gumbel with location log(sum_j exp(U_j)) and is
# independent of which alternative attains it. Given the maximum M, the chosen
# alternative has U_y + eps_y = M and every other alternative has a Gumbel(U_j)
# value truncated above at M. Returns the n x J matrix of shocks eps.
.sim_conditional_gumbel <- function(U, y) {
  n <- nrow(U)
  J <- ncol(U)
  m <- .sim_rowmax(U)
  lse <- m + log(rowSums(exp(U - m)))
  M <- lse - log(-log(stats::runif(n)))
  G <- U - log(-log(matrix(stats::runif(n * J), n, J)))
  # -log(exp(-M) + exp(-G)), computed without overflow
  W <- pmin(M, G) - log1p(exp(-abs(M - G)))
  W[cbind(seq_len(n), y)] <- M
  W - U
}

# The block of the type-probability model, which the legacy code holds as a
# multinomial logit named "TypeProb" with a constant and the factors: utility
# of type t = typeprob_t_intercept + sum_k type_t_loading_k * f_k for
# t = 2..ntyp, type 1 the reference. Its shocks are drawn conditional on the
# type already drawn, so asking for this block does not change any other draw.
.sim_legacy_type_block <- function(name, params, fmat, types, ntyp, detail) {
  n <- nrow(fmat)
  nfac <- ncol(fmat)
  U <- matrix(0.0, n, ntyp)
  Vx <- matrix(0.0, n, ntyp - 1L)
  Vfac <- lapply(seq_len(nfac), function(k) matrix(0.0, n, ntyp - 1L))
  for (t in 2:ntyp) {
    b0 <- .sim_param(params, sprintf("typeprob_%d_intercept", t))
    u <- rep(b0, n)
    Vx[, t - 1L] <- b0
    for (k in seq_len(nfac)) {
      fk <- .sim_param(params, sprintf("type_%d_loading_%d", t, k)) * fmat[, k]
      u <- u + fk
      Vfac[[k]][, t - 1L] <- fk
    }
    U[, t] <- u
  }
  ex <- exp(U - .sim_rowmax(U))
  prob <- ex / rowSums(ex)
  p <- list(type = "logit", K = ntyp, miss = rep(0L, n), has_endog = FALSE,
            prob = prob[, -1L, drop = FALSE], y = as.integer(types))
  if (detail) {
    e <- .sim_conditional_gumbel(U, types)
    p$Vx <- Vx
    p$Vfac <- Vfac
    p$eps <- e[, -1L, drop = FALSE] - e[, 1L]
  }
  .sim_legacy_block(name, p, detail, nfac, ntyp, typed = FALSE)
}

# The whole table: index, type model, the component blocks (a list of
# .sim_legacy_block() results, in simulation order), the factors (numbered from
# 0), the indicators of the non-reference types (numbered from 1), and sim_flag.
.sim_legacy_table <- function(blocks, fmat, types, ntyp, index, index_name,
                              type_block = NULL) {
  n <- nrow(fmat)
  nfac <- ncol(fmat)
  index <- as.vector(unclass(index))
  if (is.numeric(index) && !is.integer(index) && !anyNA(index) &&
      all(index == round(index)) && all(abs(index) < .Machine$integer.max)) {
    index <- as.integer(index)
  }
  cols <- list()
  cols[[index_name]] <- index
  if (!is.null(type_block)) cols <- c(cols, type_block)
  cols <- c(cols, unlist(unname(blocks), recursive = FALSE))
  for (k in seq_len(nfac)) cols[[sprintf("sim_fac%d", k - 1L)]] <- fmat[, k]
  for (t in seq_len(ntyp)[-1L]) cols[[sprintf("sim_type%d", t - 1L)]] <- as.numeric(types == t)
  cols[["sim_flag"]] <- rep(1L, n)
  dup <- unique(names(cols)[duplicated(names(cols))])
  if (length(dup) > 0L) {
    stop("simulate_factor_model: the legacy layout has duplicated column names: ",
         paste(sQuote(dup), collapse = ", "), ". Rename the components or `index_var`.")
  }
  attr(cols, "row.names") <- .set_row_names(n)
  class(cols) <- "data.frame"
  cols
}

# ---- writing the file ------------------------------------------------------

#' Write a simulation in the format of the legacy C++ `simulation_data.csv`
#'
#' Writes a table in the legacy layout of [simulate_factor_model()] with the
#' formatting of the C++ `TMinLkhd::Simulate()`: a header of column names
#' separated by a comma and a space, and one line per simulated unit in which
#' integer columns are written with the C format `%10d` and numeric columns with
#' `%10.5f`. A missing value is written as `NA`, right-aligned in its field,
#' where the C++ file has `-9999`. As in the C++ file, each data line ends with
#' a space. The storage type of a column decides its
#' format, so the index, the `_miss` columns, the outcomes of discrete
#' components and `sim_flag` are integers and everything else has five
#' decimals.
#'
#' The column names read from the file carry the leading space of the
#' separator, as they do in the C++ file, so code that trims the names of the
#' C++ file needs no change.
#'
#' @param x A table in the legacy layout, as returned by
#'   `simulate_factor_model(layout = "legacy")`.
#' @param file Path of the file to write.
#' @param chunk_rows Number of rows formatted at a time, which bounds the
#'   memory the character conversion needs.
#' @return `file`, invisibly.
#' @seealso [simulate_factor_model()]
#' @export
write_legacy_simulation_csv <- function(x, file, chunk_rows = 50000L) {
  if (!is.data.frame(x)) stop("`x` must be a data frame in the legacy simulation layout.")
  cols <- as.list(x)
  p <- length(cols)
  n <- nrow(x)
  if (p == 0L) stop("`x` has no columns.")
  ok <- vapply(cols, function(v) is.numeric(v) || is.logical(v), logical(1))
  if (!all(ok)) {
    stop("write_legacy_simulation_csv: non-numeric column(s): ",
         paste(sQuote(names(cols)[!ok]), collapse = ", "), ".")
  }
  is_int <- vapply(cols, function(v) is.integer(v) || is.logical(v), logical(1))
  missing_code <- sprintf("%10s", "NA")

  con <- file(file, open = "wb")
  on.exit(close(con))
  writeLines(paste(names(cols), collapse = ", "), con, sep = "\n")
  chunk_rows <- max(1L, as.integer(chunk_rows))
  for (start in seq.int(1L, max(n, 1L), by = chunk_rows)) {
    if (n == 0L) break
    idx <- start:min(start + chunk_rows - 1L, n)
    txt <- vector("list", p)
    for (j in seq_len(p)) {
      v <- cols[[j]][idx]
      s <- if (is_int[j]) sprintf("%10d", as.integer(v)) else sprintf("%10.5f", as.numeric(v))
      na <- is.na(v)
      if (any(na)) s[na] <- missing_code
      txt[[j]] <- s
    }
    lines <- do.call(paste, c(txt, sep = ", "))
    writeLines(paste0(lines, " "), con, sep = "\n")
  }
  invisible(file)
}

# ---- the parameter table ---------------------------------------------------

#' Parameter table in the naming of the legacy C++ code
#'
#' Returns the parameters of a model with the model and variable names that
#' the legacy C++ code wrote to `fullmodel_par_varlist.txt`. It is the table
#' obtained by merging that file with `fullmodel_par.txt`, so code that looks
#' up a coefficient as `params[model == "m" & variable == "x", estimate]` runs
#' on it.
#'
#' The naming follows `TMinLkhd::PrintParam_Varlist()`:
#' \itemize{
#'   \item the standard deviation of factor `k` is model `factor<k-1>_mix0`,
#'     variable `sigma`. factorana estimates the variance, so the estimate is
#'     its square root and the standard error follows from the delta method.
#'   \item a regression coefficient is model `<component>`, variable
#'     `<covariate>`. A multinomial logit with more than three alternatives
#'     uses model `<component>_ch<j>` for alternative `j = 2, ..., K`. With
#'     exactly three alternatives the C++ list repeats the plain model name for
#'     both non-reference alternatives, in the order of the alternatives, and
#'     so does this table.
#'   \item a free loading on factor `k` is variable `factor<k>`, and the
#'     intercept of non-reference type `t` is variable `type<t-1>`.
#'   \item the residual standard deviation of a linear component is variable
#'     `sigma`, and threshold `j` of an ordered probit is `threshold<j>`.
#'   \item the type-probability model is model `type_model_name` (with the
#'     suffix `_ch<t>` when there are more than three types), with variable
#'     `type_intercept_name` for the intercept and `factor<k>` for the
#'     loadings.
#' }
#' Coefficients fixed with [fix_coefficient()] are listed with their fixed
#' value. Parameters without a legacy name (for example mixture or
#' structural-equation parameters) are listed last with a missing `model` and
#' the factorana name as `variable`.
#'
#' @param object A `factorana_result` from [estimate_model_rcpp()], or a
#'   `model_system` together with `params`.
#' @param params Named parameter vector. Defaults to `object$estimates`.
#' @param type_model_name Model name of the type-probability model.
#' @param type_intercept_name Variable name of the intercept of the
#'   type-probability model. The legacy code uses the name of the constant
#'   regressor passed to `AddTypesModel()`.
#' @return A table (a `data.table` when the data.table package is installed,
#'   otherwise a data frame) with columns `V1` (a running index from 0),
#'   `model`, `variable`, `estimate`, `se` (when `object` is a fitted result
#'   with standard errors, otherwise `NA`) and `parameter` (the factorana
#'   parameter name).
#' @seealso [simulate_factor_model()]
#' @export
legacy_param_table <- function(object, params = NULL, type_model_name = "TypeProb",
                               type_intercept_name = "constant") {
  se_vec <- NULL
  if (inherits(object, "factorana_result")) {
    ms <- object$model_system
    if (is.null(params)) params <- object$estimates
    se_vec <- object$std_errors
    if (!is.null(se_vec) && is.null(names(se_vec)) && length(se_vec) == length(object$estimates)) {
      names(se_vec) <- names(object$estimates)
    }
  } else if (inherits(object, "model_system")) {
    ms <- object
    if (is.null(params)) stop("`params` is required when `object` is a model_system.")
  } else {
    stop("`object` must be a 'factorana_result' or a 'model_system'.")
  }
  get_se <- function(nm) {
    if (is.null(se_vec) || is.null(names(se_vec)) || !(nm %in% names(se_vec))) return(NA_real_)
    unname(se_vec[nm])
  }

  model <- character(0); variable <- character(0)
  estimate <- numeric(0); se <- numeric(0); parameter <- character(0)
  add <- function(m, v, est, s, rname) {
    model <<- c(model, m); variable <<- c(variable, v)
    estimate <<- c(estimate, est); se <<- c(se, s); parameter <<- c(parameter, rname)
  }
  add_named <- function(m, v, rname) {
    if (rname %in% names(params)) add(m, v, unname(params[rname]), get_se(rname), rname)
  }

  fm <- ms$factor
  nfac <- fm$n_factors
  ntyp <- if (is.null(fm$n_types)) 1L else as.integer(fm$n_types)
  if (ntyp > 1L && !("typeprob_2_intercept" %in% names(params))) ntyp <- 1L

  # factor standard deviations
  for (k in seq_len(nfac)) {
    rname <- sprintf("factor_var_%d", k)
    if (rname %in% names(params)) {
      sd_k <- sqrt(unname(params[rname]))
      add(sprintf("factor%d_mix0", k - 1L), "sigma", sd_k, get_se(rname) / (2 * sd_k), rname)
    }
  }

  # type-probability model
  for (t in seq_len(ntyp)[-1L]) {
    m <- if (ntyp > 3L) sprintf("%s_ch%d", type_model_name, t) else type_model_name
    add_named(m, type_intercept_name, sprintf("typeprob_%d_intercept", t))
    for (k in seq_len(nfac)) add_named(m, sprintf("factor%d", k), sprintf("type_%d_loading_%d", t, k))
  }

  # components, in listing order
  for (comp in ms$components) {
    cn <- comp$name
    K <- comp$num_choices
    multi <- comp$model_type == "logit" && !is.null(K) && K > 2L
    for (a in if (multi) seq_len(K - 1L) else 1L) {
      m <- if (multi && K > 3L) sprintf("%s_ch%d", cn, a + 1L) else cn
      pref <- if (multi) sprintf("%s_c%d", cn, a) else cn
      for (cov in comp$covariates) {
        rname <- sprintf("%s_%s", pref, cov)
        add(m, cov, .sim_beta_value(params, comp, cov, choice = if (multi) a else NULL),
            get_se(rname), rname)
      }
      for (k in seq_len(nfac)) {
        if (is.na(comp$loading_normalization[k])) {
          add_named(m, sprintf("factor%d", k), sprintf("%s_loading_%d", pref, k))
        }
      }
      if (isTRUE(comp$use_types)) {
        for (t in seq_len(ntyp)[-1L]) {
          add_named(m, sprintf("type%d", t - 1L), sprintf("%s_type_%d_intercept", pref, t))
        }
      }
    }
    if (comp$model_type == "linear") add_named(cn, "sigma", sprintf("%s_sigma", cn))
    if (comp$model_type == "oprobit") {
      for (j in seq_len(K - 1L)) add_named(cn, sprintf("threshold%d", j), sprintf("%s_thresh_%d", cn, j))
    }
  }

  # parameters without a legacy name
  for (rname in setdiff(names(params), parameter)) {
    add(NA_character_, rname, unname(params[rname]), get_se(rname), rname)
  }

  out <- data.frame(V1 = seq_along(model) - 1L, model = model, variable = variable,
                    estimate = estimate, se = se, parameter = parameter,
                    stringsAsFactors = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) data.table::setDT(out)
  out
}
