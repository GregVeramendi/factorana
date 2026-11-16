#' Read the checkpoint of an estimation
#'
#' Reads the file that [estimate_model_rcpp()] writes to its `checkpoint_file`, so that an
#' estimation that was stopped can start again from the parameters it had reached.
#'
#' @param file Path to the checkpoint file.
#' @param model_system Optional model system from [define_model_system()], the one the
#'   estimation is restarted with. The function then stops unless the checkpoint fits it:
#'   the parameter names must be those of the model, in its order, and the parameters the
#'   model holds fixed from a previous stage must have the values of that stage. Fixed
#'   parameters take their values from `init_params`, so a checkpoint written against
#'   other previous-stage results would carry those results into the new run.
#'
#' @return A list with
#'   \item{params}{named numeric vector of all parameters (free and fixed), to pass as
#'     `init_params`}
#'   \item{loglik}{log-likelihood at these parameters}
#'   \item{saved}{time the checkpoint was written, as text}
#' @examples
#' \donttest{
#' set.seed(1); n <- 200
#' f <- rnorm(n)
#' dat <- data.frame(intercept = 1,
#'   y1 = 1.0 * f + rnorm(n, 0, 0.5),
#'   y2 = 0.8 * f + rnorm(n, 0, 0.5),
#'   y3 = 0.6 * f + rnorm(n, 0, 0.5))
#' fm <- define_factor_model(n_factors = 1)
#' mc <- function(y, norm) define_model_component(y, dat, y, fm,
#'   covariates = "intercept", model_type = "linear", loading_normalization = norm)
#' ms <- define_model_system(
#'   components = list(mc("y1", 1), mc("y2", NA_real_), mc("y3", NA_real_)), factor = fm)
#' ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
#' ckpt <- tempfile(fileext = ".csv")
#' fit <- estimate_model_rcpp(ms, dat, control = ctrl, parallel = FALSE, verbose = FALSE,
#'   checkpoint_file = ckpt)
#' # A later run starts from the saved parameters
#' saved <- read_checkpoint(ckpt, ms)
#' refit <- estimate_model_rcpp(ms, dat, init_params = saved$params, control = ctrl,
#'   parallel = FALSE, verbose = FALSE)
#' }
#' @seealso [estimate_model_rcpp()]
#' @export
read_checkpoint <- function(file, model_system = NULL) {
  if (!file.exists(file)) stop("Checkpoint file not found: ", file)
  lines <- readLines(file, warn = FALSE)
  header_field <- function(label) {
    prefix <- paste0("# ", label, ": ")
    hit <- lines[startsWith(lines, prefix)]
    if (length(hit)) substring(hit[1], nchar(prefix) + 1) else NA_character_
  }
  tab <- tryCatch(
    utils::read.csv(text = lines, comment.char = "#", stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(tab) || !identical(names(tab), c("parameter", "value")) || nrow(tab) == 0 ||
      !is.numeric(tab$value) || anyNA(tab$value)) {
    stop("Not a complete checkpoint file: ", file)
  }
  params <- stats::setNames(as.numeric(tab$value), tab$parameter)

  if (!is.null(model_system)) {
    model_names <- build_parameter_metadata(model_system)$names
    if (!identical(names(params), model_names)) {
      stop(sprintf(paste0("The checkpoint %s does not fit the model: it has %d parameters and the ",
                          "model has %d, of which %d have the same name in the same position. ",
                          "Delete the file to start from the default values."),
                   file, length(params), length(model_names),
                   sum(names(params)[seq_len(min(length(params), length(model_names)))] ==
                         model_names[seq_len(min(length(params), length(model_names)))])))
    }
    previous <- model_system$previous_stage_info
    if (!is.null(previous) && length(previous$fixed_param_names) > 0) {
      # The parameters held fixed from the previous stage, located as the estimator locates them
      idx <- match(previous$fixed_param_names, model_names)
      found <- !is.na(idx)
      stage_values <- as.numeric(previous$fixed_param_values[found])
      differs <- which(abs(params[idx[found]] - stage_values) > 1e-8 * pmax(1, abs(stage_values)))
      if (length(differs) > 0) {
        stop(sprintf(paste0("The checkpoint %s was written against other previous-stage results: %d of ",
                            "the %d parameters held fixed from the previous stage have other values ",
                            "(first: %s). Delete the file to start from the default values."),
                     file, length(differs), sum(found), names(params)[idx[found]][differs[1]]))
      }
    }
  }

  list(params = params,
       loglik = as.numeric(header_field("Log-likelihood")),
       saved = header_field("Checkpoint saved"))
}
