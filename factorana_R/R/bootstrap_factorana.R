# Bootstrap inference for factorana, designed for the "estimate one stage at a
# time, across many bootstrap samples, on a compute cluster" workflow used by
# the legacy backend. The pieces are deliberately small and restartable:
#
#   generate_bootstrap_samples()  -- draw and persist the resampling weights once
#   bootstrap_fit_sample()        -- estimate ONE stage for ONE sample (idempotent)
#   collect_bootstrap()           -- read finished samples -> SE + percentile CIs
#   bootstrap_factorana()         -- single-node convenience driver over the above
#
# Resampling uses integer frequency weights (cluster or row multiplicities)
# rather than physically duplicating rows. Estimating with these weights is
# numerically identical to estimating on the expanded data, and it reuses the
# existing per-observation weight machinery (model_system$weights). A cluster
# drawn k times simply gets weight k on its rows, which sidesteps any
# duplicate-id bookkeeping.

#' Generate and persist bootstrap resampling weights
#'
#' Draws `R` bootstrap replicates by resampling clusters (or rows) with
#' replacement and records, for each replicate, the integer frequency weight of
#' every row (the number of times its cluster was drawn). The weights are saved
#' once so that the same replicate is reproducible across nodes and across
#' restarts.
#'
#' @param data Data frame used for estimation.
#' @param R Number of bootstrap replicates.
#' @param cluster Optional cluster id: a column name in `data` or a vector of
#'   length `nrow(data)`. If `NULL`, an ordinary row (nonparametric) bootstrap
#'   is used (each row is its own cluster).
#' @param dir Optional directory to write the samples object to (as
#'   `bootstrap_samples.rds`). Created if it does not exist.
#' @param seed Optional integer seed for reproducibility.
#' @return A `factorana_boot_samples` list with `weights` (an integer
#'   `nrow(data)` x `R` matrix), and metadata (`R`, `nobs`, `n_clusters`).
#'   Invisibly returns the same object when `dir` is supplied.
#' @seealso [bootstrap_fit_sample()], [collect_bootstrap()], [bootstrap_factorana()]
#' @export
generate_bootstrap_samples <- function(data, R, cluster = NULL, dir = NULL,
                                       seed = NULL) {
  nobs <- nrow(data)
  if (!is.numeric(R) || length(R) != 1L || R < 1) stop("`R` must be a positive integer.")
  R <- as.integer(R)
  if (!is.null(seed)) set.seed(seed)

  if (is.null(cluster)) {
    cl <- seq_len(nobs)  # row bootstrap: every row its own cluster
  } else if (length(cluster) == 1L && is.character(cluster)) {
    if (!cluster %in% colnames(data)) stop("cluster column '", cluster, "' not found in `data`.")
    cl <- data[[cluster]]
  } else {
    if (length(cluster) != nobs) stop("`cluster` length must equal nrow(data).")
    cl <- cluster
  }
  if (anyNA(cl)) stop("`cluster` must not contain NA.")

  uniq <- unique(cl)
  G <- length(uniq)
  if (G < 2L) stop("Need at least 2 clusters to bootstrap; found ", G, ".")
  cl_idx <- match(cl, uniq)
  rows_by_cluster <- split(seq_len(nobs), cl_idx)

  W <- matrix(0L, nobs, R)
  for (b in seq_len(R)) {
    mult <- tabulate(sample.int(G, G, replace = TRUE), nbins = G)
    w <- integer(nobs)
    for (g in which(mult > 0L)) w[rows_by_cluster[[g]]] <- mult[g]
    W[, b] <- w
  }

  samples <- list(weights = W, R = R, nobs = nobs, n_clusters = G,
                  cluster = if (is.null(cluster)) NULL else cluster, seed = seed)
  class(samples) <- "factorana_boot_samples"

  if (!is.null(dir)) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    saveRDS(samples, file.path(dir, "bootstrap_samples.rds"))
    return(invisible(samples))
  }
  samples
}

#' Estimate one stage for one bootstrap sample (restartable)
#'
#' Fits `model_system` to the data reweighted by bootstrap sample `sample_id`
#' and writes the result to `stage_dir/sample_<id>.rds`. If that file already
#' exists it returns immediately, so an interrupted run (or a rebooted node)
#' simply re-runs the missing samples. This is the unit of work to scatter
#' across a compute cluster (one array-job task per `sample_id`).
#'
#' Rows whose bootstrap weight is zero are dropped, and the remaining rows are
#' given their integer frequency weight via `model_system$weights`. For a
#' multi-stage workflow, set `model_system$previous_stage` (built from this
#' sample's previous-stage fit) before calling, so each replicate chains on its
#' own earlier stage.
#'
#' @param model_system A `model_system` (already carrying `previous_stage` for
#'   this sample if estimating a later stage).
#' @param data The full data frame used to generate the samples.
#' @param samples A `factorana_boot_samples` object or a path to a
#'   `bootstrap_samples.rds` file.
#' @param sample_id Integer replicate index (1..R).
#' @param stage_dir Directory for this stage's per-sample result files.
#' @param control Estimation control (see [define_estimation_control()]).
#' @param ... Passed to [estimate_model_rcpp()] (e.g. `optimizer`, `parallel`).
#' @param weight_col Name of the temporary weight column added to the data.
#' @param overwrite If `TRUE`, re-estimate even if the output file exists.
#' @return Invisibly, the path to the written result file.
#' @seealso [generate_bootstrap_samples()], [collect_bootstrap()]
#' @export
bootstrap_fit_sample <- function(model_system, data, samples, sample_id,
                                 stage_dir, control = NULL, ...,
                                 weight_col = ".bootw", overwrite = FALSE) {
  if (is.character(samples)) samples <- readRDS(samples)
  if (!inherits(samples, "factorana_boot_samples")) {
    stop("`samples` must be a factorana_boot_samples object or a path to one.")
  }
  if (!dir.exists(stage_dir)) dir.create(stage_dir, recursive = TRUE)
  out_file <- file.path(stage_dir, sprintf("sample_%d.rds", as.integer(sample_id)))
  if (file.exists(out_file) && !overwrite) return(invisible(out_file))

  sample_id <- as.integer(sample_id)
  if (sample_id < 1L || sample_id > samples$R) {
    stop("`sample_id` (", sample_id, ") out of range 1..", samples$R, ".")
  }
  if (nrow(data) != samples$nobs) {
    stop("nrow(data) (", nrow(data), ") does not match the sampled nobs (",
         samples$nobs, ").")
  }
  if (weight_col %in% colnames(data)) {
    stop("`weight_col` ('", weight_col, "') already exists in `data`; choose another.")
  }

  w <- samples$weights[, sample_id]
  keep <- w > 0L
  data_b <- data[keep, , drop = FALSE]
  data_b[[weight_col]] <- w[keep]

  ms_b <- model_system
  ms_b$weights <- weight_col

  fit <- estimate_model_rcpp(ms_b, data_b, control = control, ...)
  fit$bootstrap_sample_id <- sample_id
  saveRDS(fit, out_file)
  invisible(out_file)
}

#' Collect finished bootstrap samples into standard errors and intervals
#'
#' Reads every `sample_<id>.rds` in `stage_dir`, stacks the estimates, and
#' returns the bootstrap covariance, bootstrap standard errors, and percentile
#' confidence intervals. Non-converged replicates are dropped by default.
#'
#' @param stage_dir Directory of per-sample result files.
#' @param conf_level Confidence level for percentile intervals.
#' @param require_convergence Drop replicates with `convergence != 0`.
#' @return A `factorana_bootstrap` list with `boot_se`, `boot_cov`, `ci`
#'   (a data frame of percentile intervals and bootstrap SEs), the stacked
#'   `estimates`, and counts (`n_total`, `n_converged`, `n_used`).
#' @seealso [generate_bootstrap_samples()], [bootstrap_fit_sample()]
#' @export
collect_bootstrap <- function(stage_dir, conf_level = 0.95,
                              require_convergence = TRUE) {
  files <- list.files(stage_dir, pattern = "^sample_[0-9]+\\.rds$", full.names = TRUE)
  if (length(files) == 0L) stop("No sample_*.rds files found in '", stage_dir, "'.")
  ids <- as.integer(sub("^sample_([0-9]+)\\.rds$", "\\1", basename(files)))
  files <- files[order(ids)]; ids <- sort(ids)

  fits <- lapply(files, readRDS)
  conv <- vapply(fits, function(f) isTRUE(f$convergence == 0L), logical(1))
  pnames <- names(fits[[1]]$estimates)
  E <- do.call(rbind, lapply(fits, function(f) f$estimates[pnames]))
  rownames(E) <- ids
  colnames(E) <- pnames

  use <- if (require_convergence) conv else rep(TRUE, length(fits))
  if (sum(use) < 2L) {
    stop("Fewer than 2 usable bootstrap replicates (",
         sum(use), " of ", length(fits),
         if (require_convergence) " converged" else "", ").")
  }
  Ek <- E[use, , drop = FALSE]

  boot_cov <- stats::cov(Ek)
  boot_se  <- sqrt(diag(boot_cov))
  alpha <- 1 - conf_level
  qs <- t(apply(Ek, 2L, stats::quantile, probs = c(alpha / 2, 1 - alpha / 2),
                na.rm = TRUE))

  out <- list(
    boot_se = boot_se,
    boot_cov = boot_cov,
    ci = data.frame(parameter = pnames,
                    lower = unname(qs[, 1]), upper = unname(qs[, 2]),
                    boot_se = unname(boot_se), row.names = NULL),
    estimates = E,
    converged = setNames(conv, ids),
    n_total = length(fits),
    n_converged = sum(conv),
    n_used = nrow(Ek),
    conf_level = conf_level
  )
  class(out) <- "factorana_bootstrap"
  out
}

#' Single-node bootstrap driver (convenience over the primitives)
#'
#' Generates `R` bootstrap samples, fits every replicate in a local loop (each
#' fit is restartable via `stage_dir/sample_<id>.rds`), and collects the
#' results. For a single estimation stage on one machine. For multi-stage,
#' distributed, or restart-heavy workflows, use the primitives
#' ([generate_bootstrap_samples()], [bootstrap_fit_sample()],
#' [collect_bootstrap()]) directly so each stage and sample is an independent,
#' resumable job.
#'
#' @inheritParams generate_bootstrap_samples
#' @param model_system A `model_system`.
#' @param control Estimation control.
#' @param stage_dir Directory for per-sample result files (and the samples object).
#' @param conf_level Confidence level for percentile intervals.
#' @param ... Passed to [estimate_model_rcpp()].
#' @return A `factorana_bootstrap` object (see [collect_bootstrap()]).
#' @seealso [generate_bootstrap_samples()], [bootstrap_fit_sample()], [collect_bootstrap()]
#' @export
bootstrap_factorana <- function(model_system, data, R, cluster = NULL,
                                stage_dir, control = NULL, seed = NULL,
                                conf_level = 0.95, ...) {
  samples <- generate_bootstrap_samples(data, R = R, cluster = cluster,
                                         dir = stage_dir, seed = seed)
  for (b in seq_len(samples$R)) {
    bootstrap_fit_sample(model_system, data, samples, sample_id = b,
                         stage_dir = stage_dir, control = control, ...)
  }
  collect_bootstrap(stage_dir, conf_level = conf_level)
}

#' Multi-stage single-node bootstrap driver
#'
#' Runs a multi-stage bootstrap on one machine: for each replicate, fits each
#' stage in turn, chaining every stage on that replicate's own earlier-stage
#' fits. Each (stage, sample) fit is written to `dir/stage_<k>/sample_<id>.rds`
#' and skipped if it already exists, so an interrupted run resumes where it
#' stopped. If a stage fails to converge for a replicate, its later stages are
#' skipped (so a bad earlier fit is not chained forward).
#'
#' For distributed (multi-node) runs, drive the same per-stage, per-sample work
#' with [bootstrap_fit_sample()] directly (one array-job task per sample), which
#' is what this driver does internally.
#'
#' @param stage_builders A list of functions, one per stage, each with the
#'   signature `function(prev_fits, data)` returning the `model_system` for that
#'   stage. For the first stage `prev_fits` is an empty list; for later stages
#'   `prev_fits[[k]]` is the fitted result of stage `k` for the current
#'   replicate (use it to construct the stage's `previous_stage`).
#' @param data Data frame.
#' @param R Number of bootstrap replicates.
#' @param cluster Optional cluster id (column name in `data` or a vector); see
#'   [generate_bootstrap_samples()].
#' @param dir Directory holding the samples object and the per-stage
#'   subdirectories (`stage_1/`, `stage_2/`, ...).
#' @param control Estimation control.
#' @param seed Optional seed for sample generation.
#' @param conf_level Confidence level for percentile intervals.
#' @param ... Passed to [estimate_model_rcpp()] for every stage fit.
#' @return A `factorana_bootstrap_multistage` list with a per-stage
#'   [collect_bootstrap()] summary in `stages` and `final` pointing at the last
#'   stage.
#' @seealso [bootstrap_factorana()], [bootstrap_fit_sample()], [collect_bootstrap()]
#' @export
bootstrap_factorana_multistage <- function(stage_builders, data, R, cluster = NULL,
                                           dir, control = NULL, seed = NULL,
                                           conf_level = 0.95, ...) {
  if (!is.list(stage_builders) || length(stage_builders) < 1L ||
      !all(vapply(stage_builders, is.function, logical(1)))) {
    stop("`stage_builders` must be a non-empty list of functions ",
         "with signature function(prev_fits, data).")
  }
  n_stages <- length(stage_builders)
  samples <- generate_bootstrap_samples(data, R = R, cluster = cluster,
                                        dir = dir, seed = seed)
  stage_dirs <- file.path(dir, sprintf("stage_%d", seq_len(n_stages)))

  for (b in seq_len(samples$R)) {
    prev_fits <- list()
    for (k in seq_len(n_stages)) {
      if (k > 1L && !isTRUE(prev_fits[[k - 1L]]$convergence == 0L)) {
        message(sprintf("Sample %d: stage %d did not converge; skipping later stages.",
                        b, k - 1L))
        break
      }
      ms_k <- stage_builders[[k]](prev_fits, data)
      f <- bootstrap_fit_sample(ms_k, data, samples, sample_id = b,
                                stage_dir = stage_dirs[k], control = control, ...)
      prev_fits[[k]] <- readRDS(f)
    }
  }

  stages <- lapply(stage_dirs, function(d) {
    if (length(list.files(d, pattern = "^sample_[0-9]+\\.rds$")) > 0L) {
      collect_bootstrap(d, conf_level = conf_level, require_convergence = TRUE)
    } else {
      NULL
    }
  })
  names(stages) <- sprintf("stage_%d", seq_len(n_stages))

  out <- list(stages = stages, final = stages[[n_stages]],
              n_stages = n_stages, n_replicates = samples$R)
  class(out) <- "factorana_bootstrap_multistage"
  out
}

#' @export
print.factorana_bootstrap_multistage <- function(x, ...) {
  cat("factorana multi-stage bootstrap:", x$n_stages, "stages,",
      x$n_replicates, "replicates\n\n")
  for (k in seq_len(x$n_stages)) {
    st <- x$stages[[k]]
    if (is.null(st)) {
      cat(sprintf("Stage %d: (no completed samples)\n", k))
    } else {
      cat(sprintf("Stage %d: %d of %d replicates used (%d converged)\n",
                  k, st$n_used, st$n_total, st$n_converged))
    }
  }
  cat("\nFinal-stage percentile intervals:\n")
  if (!is.null(x$final)) print(utils::head(x$final$ci, 20), row.names = FALSE)
  invisible(x)
}

#' @export
print.factorana_bootstrap <- function(x, ...) {
  cat("factorana bootstrap:", x$n_used, "of", x$n_total,
      "replicates used (", x$n_converged, "converged)\n")
  cat("Confidence level:", x$conf_level, "\n\n")
  print(utils::head(x$ci, 20), row.names = FALSE)
  if (nrow(x$ci) > 20) cat("... (", nrow(x$ci) - 20, "more parameters)\n")
  invisible(x)
}
