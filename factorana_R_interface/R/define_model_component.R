#' Define a model component
#'
#' @param name name of model component ####Long name and s name (like in the C++)
#' @param data data.frame for validation
#' @param outcome Character. Name of the outcome variable.
#' @param factor model. Object of Class factor_model.
#' @param evaluation_indicator Character (optional). Variable used for evaluation subsample.
#' @param covariates List of Character vectors. Names of covariates.
#' @param model_type Character. Type of model (e.g., "linear", "logit", "probit").
#' @param intercept Logical. Whether to include an intercept (default = TRUE).
#' @param factor_normalization Numeric. Factor loading normalization constant (default = 1).
#' @param num_choices Integer. Number of choices (for multinomial models).
#' @param nrank Integer (optional). Rank for exploded multinomial logit.
#'
#' @return An object of class "model_component". A list representing the model component
#' @export

define_model_component <- function(name,
                                   data,
                                   outcome,
                                   factor,
                                   evaluation_indicator = NULL,
                                   covariates,
                                   model_type = c("linear", "logit", "probit", "oprobit"),
                                   intercept = TRUE,
                                   factor_normalization = 1,
                                   num_choices = 2,
                                   nrank = NULL) {
  #_____Basic checks____
  #is the data a dataframe?
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")

  #Does 'outcome' exist in data.frame
  if (!(outcome %in% names(data))) stop("outcome must be a column in data.frame")

  #is the factor model a list?
  if (!is.list(factor)) stop("`factor` must be a list.") #should it be a list or object of S3 class???

  #evaluation indicator?
  if (!is.null(evaluation_indicator)) {
    if (!is.character(evaluation_indicator) || length(evaluation_indicator) != 1L)
      stop("`evaluation_indicator` must be a single column name (character).")
    if (!(evaluation_indicator %in% names(data)))
      stop("`evaluation_indicator` not found in `data`.")
  }

  #covariates existence:

  # covariates must be a non-empty character vector
  if (!is.character(covariates) || length(covariates) < 1L) {
    stop("`covariates` must be a non-empty character vector of column names.")
  }

  missing <- setdiff(covariates, names(data))
  if (length(missing)) {
    stop("Covariates not found in `data`: ", paste(missing, collapse = ", "))
  }

  #do the model types match? e.g. not "linnnear"
  model_type <- match.arg(model_type)

  #intercept: is it a boolean?
  if (!is.logical(intercept) || length(intercept) != 1L) {
    stop("`intercept` must be a single TRUE/FALSE.")
  }

  #factor_normalization: is it a constant?
  if (!is.numeric(factor_normalization) || length(factor_normalization) != 1L) {
    stop("`factor_normalization` must be a single numeric value.")
  }

  #num_choices: is it within X? (speciy X later)
  num_choices <- as.integer(num_choices)
  if (model_type %in% c("logit", "probit") && num_choices != 2L) {
    stop("For binary ", model_type, ", `num_choices` must be 2.")
  }

  #n_rank: is it either NULL or > 1
  if (!is.null(nrank)) {
    nrank <- as.integer(nrank)
    if (!is.finite(nrank) || nrank < 1L) stop("`nrank` must be a positive integer when provided.")
  }

  ## ---- Data validity checks ----

  # outcome must not have missing when evaluation_indicator is 1
  idx <- rep(TRUE, nrow(data))  # default: check all rows
  if (!is.null(evaluation_indicator)) {
    ei <- data[[evaluation_indicator]]
    if (is.logical(ei)) {
      idx <- !is.na(ei) & ei
    } else if (is.numeric(ei) || is.integer(ei)) {
      idx <- !is.na(ei) & (ei == 1L)
    } else {
      stop("`evaluation_indicator` must be logical or 0/1 numeric.")
    }
  }

  # *** FIX: subset the ENTIRE data frame here ***
  data <- data[idx, , drop = FALSE]
  rownames(data) <- NULL


  if (nrow(data) == 0L) stop("Evaluation subset has zero rows")

  idx <- rep(TRUE, nrow(data))

  # 2) If ordered probit, coerce OUTCOME to ordered factor *now* (post-subset)
  if (model_type == "oprobit") {
    y_sub <- data[[outcome]]

    if (is.factor(y_sub)) {
      # Make sure it is 'ordered'
      if (!is.ordered(y_sub)) y_sub <- ordered(y_sub)
    } else {
      # Accept integer-like labels (1..J or 0..J-1), otherwise error
      if (!is.numeric(y_sub) && !is.integer(y_sub))
        stop("oprobit outcome must be integer-like or an ordered factor.")
      # map to contiguous 1..J and mark ordered
      u   <- sort(unique(na.omit(as.integer(y_sub))))
      map <- match(as.integer(y_sub), u)
      y_sub <- ordered(map)
    }

    if (nlevels(y_sub) < 3L)
      stop("Ordered probit requires an outcome with ≥ 3 ordered categories.")

    data[[outcome]] <- y_sub
  }

  # Check outcome missingness on the subset
  if (anyNA(data[[outcome]][idx])) {
    stop("Missing values in outcome variable within evaluation subset.")
  }

  # covariates must not have missing (on the same subset)
  for (cov in covariates) {
    if (anyNA(data[[cov]][idx])) {
      stop("Missing values found in covariate: ", cov)
    }
  }


  # outcome must not have missing when evaluation_indicator is 1
#
#   if (any(is.na(data[outcome]))){
#     stop("Missing values in outcome variable")
#   }
#
#   # covariates must not have missing
#   for (cov in covariates) {
#     if (any(is.na(data[cov]))){
#       stop("Missing values found in covariate:", cov)
#     }
#   }

  # outcome must match model_type (checks for probit and logit)
  y <- data[[outcome]]
  if (model_type == "probit" && !all(y %in% c(0, 1))) {
    stop("Outcome for probit must be coded 0/1.")
  }

  if (model_type == "logit" && !all(y %in% 0:(max(y)))) {
    stop("Outcome for logit must be integers 0,...,K.")
  }


  ###ordered probit checks
  n_cats <- NULL
  if (model_type == "oprobit") {
    y_sub <- data[[outcome]][idx]
    # Accept integers 1..J, 0..J-1, or an ordered factor; coerce to ordered factor
    if (is.factor(y_sub)) {
      if (!is.ordered(y_sub)) y_sub <- ordered(y_sub) #coerce to an ordered factor here.
    } else {
      if (!is.numeric(y_sub) && !is.integer(y_sub))
        stop("oprobit outcome must be integer-like or ordered factor.")

      u <- sort(unique(na.omit(as.integer(y_sub))))
      # make contiguous levels starting at 1
      map <- match(as.integer(y_sub), u)
      y_sub <- ordered(map)
    }
    n_cats <- length(levels(y_sub))
    if (n_cats < 2L) stop("oprobit needs at least 2 ordered categories.")
    # (Optional) replace in data to keep consistency downstream:
    data[[outcome]] <- y_sub
  }


  #___output___

  out <- list(
    name = name,
    data = as.data.frame(data),
    outcome = outcome,
    factor = factor,
    evaluation_indicator = evaluation_indicator,
    covariates = covariates,
    model_type = model_type,
    intercept = intercept,
    factor_normalization = factor_normalization,
    num_choices = num_choices,
    nrank = nrank,
    nparam_model = length(covariates) + factor$n_factors + factor$n_types - 1 #this function needs to know the number of factors and number of types, need to pass from the factor model (?)
  )

  class(out) <- "model_component"
  return(out)
}

#  Getter for name(S3 generic + method)
#' @export
get_component_name <- function(x, ...) UseMethod("get_component_name") #this is the S3 generic

#' @export
get_component_name.model_component <- function(x, ...) {
  nm <- x$name
  if (is.null(nm) || is.na(nm) || !nzchar(nm)) NA_character_ else nm
}

# Getter for factor (S3 generic + method)
#' @export
get_factor <- function(x, ...) UseMethod("get_factor") #S3 generic

#' @export
get_factor.model_component <- function(x, ...) {
  return(x$factor)
}


#' Print method for model_component objects
#'
#' @param x An object of class "model_component".
#' @param ... Not used.
#' @export
print.model_component <- function(x, ...) {
  cat("Model Component\n")
  cat("------------------------------\n")
  cat("Model:                   ", x$name, "\n")
  cat("Outcome variable:        ", x$outcome, "\n")
  cat("Model type:              ", x$model_type, "\n")
  cat("Intercept:               ", ifelse(x$intercept, "Yes", "No"), "\n")
  cat("Factor normalization:    ", x$factor_normalization, "\n")
  cat("Number of choices:       ", x$num_choices, "\n")
  if (!is.null(x$nrank)) {
    cat("Rank (nrank):            ", x$nrank, "\n")
  }
  if (!is.null(x$evaluation_indicator)) {
    cat("Evaluation indicator:    ", x$evaluation_indicator, "\n")
  }
  cat("Covariates:              ", paste(x$covariates, collapse = ", "), "\n")

  cat("Total parameters:        ", x$nparam_model, "\n")
  invisible(x)
}

