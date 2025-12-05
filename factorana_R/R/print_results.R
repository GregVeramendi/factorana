#' Print and Summary Methods for Factor Model Results
#'
#' Functions for displaying and exporting factor model estimation results
#' in formatted tables, including LaTeX output.

#' Print method for factorana_result objects
#'
#' @param x A factorana_result object from estimate_model_rcpp()
#' @param digits Number of decimal places (default 4)
#' @param ... Additional arguments (ignored)
#'
#' @export
print.factorana_result <- function(x, digits = 4, ...) {
  cat("\n")
  cat("Factor Model Estimation Results\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")

  # Convergence status
  converged <- (x$convergence == 0)
  cat("Convergence: ", if (converged) "Yes" else "No", "\n")
  cat("Optimizer:   ", x$optimizer, "\n")
  cat("Log-lik:     ", format(round(x$loglik, digits), nsmall = digits), "\n")
  cat("Parameters:  ", length(x$estimates), "\n\n")

  # Build parameter table
  param_df <- .build_param_table(x)

  # Print by component
  components <- unique(param_df$component)

  for (comp in components) {
    comp_df <- param_df[param_df$component == comp, ]
    cat(paste0("--- ", comp, " ---\n"))

    for (i in seq_len(nrow(comp_df))) {
      est_str <- format(round(comp_df$estimate[i], digits), nsmall = digits, width = 10)
      se_str <- format(round(comp_df$se[i], digits), nsmall = digits, width = 10)
      cat(sprintf("  %-25s %s (%s)\n",
                  comp_df$param_label[i], est_str, se_str))
    }
    cat("\n")
  }

  invisible(x)
}


#' Summary method for factorana_result objects
#'
#' @param object A factorana_result object from estimate_model_rcpp()
#' @param ... Additional arguments (ignored)
#'
#' @return A summary object with formatted tables
#' @export
summary.factorana_result <- function(object, ...) {
  param_df <- .build_param_table(object)

  # Add z-values and p-values
  param_df$z_value <- param_df$estimate / param_df$se
  param_df$p_value <- 2 * (1 - pnorm(abs(param_df$z_value)))

  # Add significance stars
  param_df$signif <- ""
  param_df$signif[param_df$p_value < 0.1] <- "."

  param_df$signif[param_df$p_value < 0.05] <- "*"
  param_df$signif[param_df$p_value < 0.01] <- "**"
  param_df$signif[param_df$p_value < 0.001] <- "***"

  result <- list(
    coefficients = param_df,
    loglik = object$loglik,
    convergence = object$convergence,
    optimizer = object$optimizer,
    n_params = length(object$estimates)
  )
  class(result) <- "summary.factorana_result"
  result
}


#' Print method for summary.factorana_result
#'
#' @param x A summary.factorana_result object
#' @param digits Number of decimal places (default 4)
#' @param ... Additional arguments (ignored)
#'
#' @export
print.summary.factorana_result <- function(x, digits = 4, ...) {
  cat("\n")
  cat("Factor Model Estimation Results\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  cat("Convergence: ", if (x$convergence == 0) "Yes" else "No", "\n")
  cat("Optimizer:   ", x$optimizer, "\n")
  cat("Log-lik:     ", format(round(x$loglik, digits), nsmall = digits), "\n")
  cat("Parameters:  ", x$n_params, "\n\n")

  # Print coefficient table
  coef_df <- x$coefficients

  # Header
  cat(sprintf("%-25s %10s %10s %8s %8s %s\n",
              "Parameter", "Estimate", "Std.Error", "z value", "Pr(>|z|)", ""))
  cat(paste(rep("-", 70), collapse = ""), "\n")

  components <- unique(coef_df$component)

  for (comp in components) {
    cat(paste0("\n", comp, ":\n"))
    comp_df <- coef_df[coef_df$component == comp, ]

    for (i in seq_len(nrow(comp_df))) {
      row <- comp_df[i, ]
      p_str <- if (row$p_value < 0.0001) "< 0.0001" else format(round(row$p_value, 4), nsmall = 4)
      cat(sprintf("  %-23s %10.4f %10.4f %8.3f %8s %s\n",
                  row$param_label,
                  row$estimate,
                  row$se,
                  row$z_value,
                  p_str,
                  row$signif))
    }
  }

  cat("\n---\n")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

  invisible(x)
}


#' Create a formatted results table for multiple models
#'
#' Creates a table comparing estimates across multiple models, with standard
#' errors in parentheses below each estimate.
#'
#' @param ... One or more factorana_result objects, or a list of them
#' @param model_names Optional character vector of model names for column headers
#' @param digits Number of decimal places (default 3)
#' @param se_format How to display standard errors: "parentheses" (default),
#'   "brackets", or "below"
#' @param include_loglik Include log-likelihood row (default TRUE)
#' @param include_n Include number of observations (default TRUE)
#' @param stars Add significance stars (default TRUE)
#'
#' @return A data frame with the formatted table
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare two models
#' results_table(model1, model2, model_names = c("Baseline", "Extended"))
#'
#' # With custom formatting
#' results_table(model1, model2, model3,
#'               digits = 4,
#'               se_format = "brackets",
#'               stars = FALSE)
#' }
results_table <- function(..., model_names = NULL, digits = 3,
                          se_format = c("parentheses", "brackets", "below"),
                          include_loglik = TRUE, include_n = TRUE,
                          stars = TRUE) {

  se_format <- match.arg(se_format)

  # Collect models
  args <- list(...)
  if (length(args) == 1 && is.list(args[[1]]) && !inherits(args[[1]], "factorana_result")) {
    models <- args[[1]]
  } else {
    models <- args
  }

  n_models <- length(models)
  if (n_models == 0) {
    stop("At least one model result is required")
  }
  if (n_models > 5) {
    warning("More than 5 models may not display well. Consider splitting into multiple tables.")
  }

  # Set model names
  if (is.null(model_names)) {
    model_names <- paste0("(", seq_len(n_models), ")")
  } else if (length(model_names) != n_models) {
    stop("Length of model_names must match number of models")
  }

  # Collect all parameter names across all models
  all_params <- list()
  for (m in seq_len(n_models)) {
    param_df <- .build_param_table(models[[m]])
    for (i in seq_len(nrow(param_df))) {
      key <- paste0(param_df$component[i], "___", param_df$param_label[i])
      if (!(key %in% names(all_params))) {
        all_params[[key]] <- list(
          component = param_df$component[i],
          param_label = param_df$param_label[i]
        )
      }
    }
  }

  # Build results matrix
  param_keys <- names(all_params)
  n_params <- length(param_keys)

  # Create data frame structure
  result_df <- data.frame(
    Parameter = character(n_params),
    stringsAsFactors = FALSE
  )

  # Add model columns
  for (m in seq_len(n_models)) {
    result_df[[model_names[m]]] <- character(n_params)
  }

  # Fill in parameter names and values
  for (i in seq_along(param_keys)) {
    key <- param_keys[i]
    result_df$Parameter[i] <- all_params[[key]]$param_label
  }

  # Fill in estimates for each model
  for (m in seq_len(n_models)) {
    param_df <- .build_param_table(models[[m]])

    # Add significance if requested
    if (stars) {
      param_df$z_value <- param_df$estimate / param_df$se
      param_df$p_value <- 2 * (1 - pnorm(abs(param_df$z_value)))
      param_df$star <- ""
      param_df$star[param_df$p_value < 0.1] <- "."
      param_df$star[param_df$p_value < 0.05] <- "*"
      param_df$star[param_df$p_value < 0.01] <- "**"
      param_df$star[param_df$p_value < 0.001] <- "***"
    }

    for (i in seq_along(param_keys)) {
      key <- param_keys[i]
      info <- all_params[[key]]

      # Find matching row in this model's parameters
      match_idx <- which(param_df$component == info$component &
                           param_df$param_label == info$param_label)

      if (length(match_idx) == 1) {
        est <- round(param_df$estimate[match_idx], digits)
        se <- round(param_df$se[match_idx], digits)
        star_str <- if (stars) param_df$star[match_idx] else ""

        # Format based on se_format
        if (se_format == "parentheses") {
          result_df[[model_names[m]]][i] <- sprintf("%.*f%s\n(%.*f)",
                                                    digits, est, star_str, digits, se)
        } else if (se_format == "brackets") {
          result_df[[model_names[m]]][i] <- sprintf("%.*f%s\n[%.*f]",
                                                    digits, est, star_str, digits, se)
        } else {  # "below"
          result_df[[model_names[m]]][i] <- sprintf("%.*f%s\n%.*f",
                                                    digits, est, star_str, digits, se)
        }
      } else {
        result_df[[model_names[m]]][i] <- ""
      }
    }
  }

  # Add log-likelihood row
  if (include_loglik) {
    loglik_row <- data.frame(Parameter = "Log-likelihood", stringsAsFactors = FALSE)
    for (m in seq_len(n_models)) {
      loglik_row[[model_names[m]]] <- format(round(models[[m]]$loglik, 2), nsmall = 2)
    }
    result_df <- rbind(result_df, loglik_row)
  }

  # Add component grouping as attribute for nicer display
  components <- sapply(param_keys, function(k) all_params[[k]]$component)
  attr(result_df, "components") <- components
  attr(result_df, "model_names") <- model_names
  attr(result_df, "se_format") <- se_format

  class(result_df) <- c("factorana_table", "data.frame")
  result_df
}


#' Print method for factorana_table
#'
#' @param x A factorana_table object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.factorana_table <- function(x, ...) {
  # Get components for grouping
  components <- attr(x, "components")

  cat("\n")
  cat("Factor Model Comparison\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  # Print header
  model_names <- attr(x, "model_names")
  header_fmt <- paste0("%-25s", paste(rep(" %12s", length(model_names)), collapse = ""), "\n")
  cat(do.call(sprintf, c(list(header_fmt, "Parameter"), as.list(model_names))))
  cat(paste(rep("-", 60), collapse = ""), "\n")

  # Print by component
  if (!is.null(components)) {
    unique_comps <- unique(components)
    for (comp in unique_comps) {
      cat(paste0("\n", comp, ":\n"))
      comp_idx <- which(components == comp)

      for (i in comp_idx) {
        # Split estimate and SE
        vals <- strsplit(x[[2]][i], "\n")[[1]]
        if (length(vals) >= 1) {
          row_fmt <- paste0("  %-23s", paste(rep(" %12s", length(model_names)), collapse = ""), "\n")

          # Get values for all models
          model_vals <- character(length(model_names))
          for (m in seq_along(model_names)) {
            cell <- x[[m + 1]][i]
            parts <- strsplit(cell, "\n")[[1]]
            model_vals[m] <- if (length(parts) >= 1) parts[1] else ""
          }

          cat(do.call(sprintf, c(list(row_fmt, x$Parameter[i]), as.list(model_vals))))

          # Print SE row
          se_vals <- character(length(model_names))
          for (m in seq_along(model_names)) {
            cell <- x[[m + 1]][i]
            parts <- strsplit(cell, "\n")[[1]]
            se_vals[m] <- if (length(parts) >= 2) parts[2] else ""
          }
          cat(do.call(sprintf, c(list(row_fmt, ""), as.list(se_vals))))
        }
      }
    }
  }

  # Print log-likelihood row if present
  loglik_idx <- which(x$Parameter == "Log-likelihood")
  if (length(loglik_idx) > 0) {
    cat("\n")
    cat(paste(rep("-", 60), collapse = ""), "\n")
    row_fmt <- paste0("%-25s", paste(rep(" %12s", length(model_names)), collapse = ""), "\n")
    model_vals <- sapply(seq_along(model_names), function(m) x[[m + 1]][loglik_idx])
    cat(do.call(sprintf, c(list(row_fmt, "Log-likelihood"), as.list(model_vals))))
  }

  cat("\n")
  cat("Standard errors in ", attr(x, "se_format"), "\n", sep = "")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1\n")

  invisible(x)
}


#' Export results table to LaTeX
#'
#' Creates a LaTeX table from factorana results, suitable for academic papers.
#'
#' @param ... One or more factorana_result objects, or a list of them
#' @param model_names Optional character vector of model names for column headers
#' @param digits Number of decimal places (default 3)
#' @param file Optional file path to write the LaTeX table
#' @param caption Table caption (optional)
#' @param label Table label for cross-referencing (optional)
#' @param stars Add significance stars (default TRUE)
#' @param note Footnote text (optional, default includes significance codes)
#' @param booktabs Use booktabs package formatting (default TRUE)
#' @param include_loglik Include log-likelihood row (default TRUE)
#' @param param_labels Optional named list mapping parameter names to display labels
#'
#' @return A character string containing the LaTeX table code
#' @export
#'
#' @examples
#' \dontrun{
#' # Export to file
#' results_to_latex(model1, model2,
#'                  model_names = c("Baseline", "Extended"),
#'                  file = "results_table.tex",
#'                  caption = "Factor Model Estimates",
#'                  label = "tab:results")
#'
#' # Custom parameter labels
#' results_to_latex(model1,
#'                  param_labels = list(
#'                    "intercept" = "Constant",
#'                    "loading_1" = "$\\lambda_1$"
#'                  ))
#' }
results_to_latex <- function(..., model_names = NULL, digits = 3,
                              file = NULL, caption = NULL, label = NULL,
                              stars = TRUE, note = NULL, booktabs = TRUE,
                              include_loglik = TRUE, param_labels = NULL) {

  # Collect models
  args <- list(...)
  if (length(args) == 1 && is.list(args[[1]]) && !inherits(args[[1]], "factorana_result")) {
    models <- args[[1]]
  } else {
    models <- args
  }

  n_models <- length(models)

  # Set model names
  if (is.null(model_names)) {
    model_names <- paste0("(", seq_len(n_models), ")")
  }

  # Build the results table data
  tbl <- results_table(models, model_names = model_names, digits = digits,
                       se_format = "parentheses", include_loglik = include_loglik,
                       stars = stars)
  components <- attr(tbl, "components")

  # Start building LaTeX
  lines <- character()

  # Table preamble
  lines <- c(lines, "\\begin{table}[htbp]")
  lines <- c(lines, "\\centering")

  if (!is.null(caption)) {
    lines <- c(lines, sprintf("\\caption{%s}", caption))
  }
  if (!is.null(label)) {
    lines <- c(lines, sprintf("\\label{%s}", label))
  }

  # Column specification
  col_spec <- paste0("l", paste(rep("c", n_models), collapse = ""))
  lines <- c(lines, sprintf("\\begin{tabular}{%s}", col_spec))

  if (booktabs) {
    lines <- c(lines, "\\toprule")
  } else {
    lines <- c(lines, "\\hline")
  }

  # Header row
  header <- paste(c("", model_names), collapse = " & ")
  lines <- c(lines, paste0(header, " \\\\"))

  if (booktabs) {
    lines <- c(lines, "\\midrule")
  } else {
    lines <- c(lines, "\\hline")
  }

  # Parameter rows by component
  if (!is.null(components)) {
    unique_comps <- unique(components)

    for (comp_idx in seq_along(unique_comps)) {
      comp <- unique_comps[comp_idx]

      # Component header
      lines <- c(lines, sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", n_models + 1, comp))

      param_indices <- which(components == comp)

      for (i in param_indices) {
        # Get parameter label
        param_name <- tbl$Parameter[i]
        if (!is.null(param_labels) && param_name %in% names(param_labels)) {
          param_name <- param_labels[[param_name]]
        }

        # Get estimates row
        est_vals <- character(n_models)
        se_vals <- character(n_models)

        for (m in seq_len(n_models)) {
          cell <- tbl[[m + 1]][i]
          parts <- strsplit(cell, "\n")[[1]]
          est_vals[m] <- if (length(parts) >= 1 && parts[1] != "") parts[1] else ""
          se_vals[m] <- if (length(parts) >= 2) parts[2] else ""
        }

        # Estimate row
        est_row <- paste(c(paste0("\\quad ", param_name), est_vals), collapse = " & ")
        lines <- c(lines, paste0(est_row, " \\\\"))

        # SE row
        se_row <- paste(c("", se_vals), collapse = " & ")
        lines <- c(lines, paste0(se_row, " \\\\"))
      }

      # Add small space between components
      if (comp_idx < length(unique_comps)) {
        lines <- c(lines, "[0.5em]")
      }
    }
  }

  # Log-likelihood row
  loglik_idx <- which(tbl$Parameter == "Log-likelihood")
  if (length(loglik_idx) > 0) {
    if (booktabs) {
      lines <- c(lines, "\\midrule")
    } else {
      lines <- c(lines, "\\hline")
    }

    loglik_vals <- sapply(seq_len(n_models), function(m) tbl[[m + 1]][loglik_idx])
    loglik_row <- paste(c("Log-likelihood", loglik_vals), collapse = " & ")
    lines <- c(lines, paste0(loglik_row, " \\\\"))
  }

  # Bottom rule
  if (booktabs) {
    lines <- c(lines, "\\bottomrule")
  } else {
    lines <- c(lines, "\\hline")
  }

  lines <- c(lines, "\\end{tabular}")

  # Footnote
  if (is.null(note) && stars) {
    note <- "Standard errors in parentheses. $^{***}p<0.001$, $^{**}p<0.01$, $^*p<0.05$, $^{.}p<0.1$"
  }
  if (!is.null(note)) {
    lines <- c(lines, sprintf("\\\\\\footnotesize{%s}", note))
  }

  lines <- c(lines, "\\end{table}")

  # Combine lines
  latex_code <- paste(lines, collapse = "\n")

  # Write to file if specified
  if (!is.null(file)) {
    writeLines(latex_code, file)
    message(sprintf("LaTeX table written to: %s", file))
  }

  invisible(latex_code)
}


#' Create a components table for a single model
#'
#' Creates a table with model components as columns and parameter types as rows.
#' This is similar to how SEM/CFA software displays results.
#'
#' @param result A factorana_result object from estimate_model_rcpp()
#' @param digits Number of decimal places (default 3)
#' @param show_se Show standard errors in parentheses (default TRUE)
#' @param stars Add significance stars (default TRUE)
#'
#' @return A data frame with the formatted table
#' @export
#'
#' @examples
#' \dontrun{
#' # Display components table for a Roy model
#' components_table(result)
#' }
components_table <- function(result, digits = 3, show_se = TRUE, stars = TRUE) {
  if (!inherits(result, "factorana_result")) {
    stop("result must be a factorana_result object")
  }

  # Get the parameter table
  param_df <- .build_param_table(result)

  # Add z-values and p-values if showing stars
  if (stars) {
    param_df$z_value <- param_df$estimate / param_df$se
    param_df$p_value <- 2 * (1 - pnorm(abs(param_df$z_value)))
    param_df$star <- ""
    param_df$star[param_df$p_value < 0.1] <- "."
    param_df$star[param_df$p_value < 0.05] <- "*"
    param_df$star[param_df$p_value < 0.01] <- "**"
    param_df$star[param_df$p_value < 0.001] <- "***"
  }

  # Get unique components (excluding "Factor Model")
  components <- unique(param_df$component)
  components <- components[components != "Factor Model"]

  # Standardize parameter types for alignment across components
  # Map specific parameter names to generic types
  param_df$generic_type <- param_df$param_type

  # For betas, extract the variable name
  for (i in seq_len(nrow(param_df))) {
    if (param_df$param_type[i] == "beta") {
      # Use the actual covariate name as the generic type
      param_df$generic_type[i] <- paste0("beta_", param_df$param_label[i])
    } else if (param_df$param_type[i] == "intercept") {
      param_df$generic_type[i] <- "Intercept"
    } else if (param_df$param_type[i] == "loading") {
      param_df$generic_type[i] <- param_df$param_label[i]
    } else if (param_df$param_type[i] == "sigma") {
      param_df$generic_type[i] <- "Sigma"
    } else if (param_df$param_type[i] == "cutpoint") {
      param_df$generic_type[i] <- param_df$param_label[i]
    } else if (param_df$param_type[i] == "factor_var") {
      param_df$generic_type[i] <- param_df$param_label[i]
    }
  }

  # Collect all unique parameter types across all components
  all_param_types <- unique(param_df$generic_type[param_df$component != "Factor Model"])

  # Order parameter types logically
  type_order <- c("Intercept")
  # Add betas
  beta_types <- all_param_types[grepl("^beta_", all_param_types)]
  type_order <- c(type_order, sort(beta_types))
  # Add loadings
  loading_types <- all_param_types[grepl("^Loading", all_param_types)]
  type_order <- c(type_order, sort(loading_types))
  # Add sigma
  if ("Sigma" %in% all_param_types) {
    type_order <- c(type_order, "Sigma")
  }
  # Add thresholds
  threshold_types <- all_param_types[grepl("^Threshold", all_param_types)]
  type_order <- c(type_order, sort(threshold_types))
  # Add any remaining types
  remaining <- setdiff(all_param_types, type_order)
  type_order <- c(type_order, remaining)

  # Remove types not in data
  type_order <- type_order[type_order %in% all_param_types]

  # Build the result matrix
  result_df <- data.frame(
    Parameter = type_order,
    stringsAsFactors = FALSE
  )

  # Add column for each component
  for (comp in components) {
    result_df[[comp]] <- character(length(type_order))

    comp_params <- param_df[param_df$component == comp, ]

    for (i in seq_along(type_order)) {
      ptype <- type_order[i]
      match_idx <- which(comp_params$generic_type == ptype)

      if (length(match_idx) == 1) {
        est <- round(comp_params$estimate[match_idx], digits)
        se <- round(comp_params$se[match_idx], digits)
        star_str <- if (stars) comp_params$star[match_idx] else ""

        if (show_se) {
          result_df[[comp]][i] <- sprintf("%.*f%s\n(%.*f)", digits, est, star_str, digits, se)
        } else {
          result_df[[comp]][i] <- sprintf("%.*f%s", digits, est, star_str)
        }
      } else {
        result_df[[comp]][i] <- ""
      }
    }
  }

  # Collect n_obs for each component from model_system
  n_obs_vec <- integer(length(components))
  names(n_obs_vec) <- components
  for (i in seq_along(components)) {
    comp_name <- components[i]
    # Find matching component in model_system
    for (comp in result$model_system$components) {
      if (comp$name == comp_name) {
        n_obs_vec[i] <- comp$n_obs
        break
      }
    }
  }

  # Add factor variance info as separate section
  factor_params <- param_df[param_df$component == "Factor Model", ]
  attr(result_df, "factor_params") <- factor_params
  attr(result_df, "components") <- components
  attr(result_df, "n_obs") <- n_obs_vec
  attr(result_df, "show_se") <- show_se
  attr(result_df, "loglik") <- result$loglik
  attr(result_df, "digits") <- digits

  class(result_df) <- c("components_table", "data.frame")
  result_df
}


#' Print method for components_table
#'
#' @param x A components_table object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.components_table <- function(x, ...) {
  components <- attr(x, "components")
  factor_params <- attr(x, "factor_params")
  n_obs <- attr(x, "n_obs")
  show_se <- attr(x, "show_se")
  digits <- attr(x, "digits")
  loglik <- attr(x, "loglik")

  n_cols <- length(components)

  cat("\n")
  cat("Factor Model Results by Component\n")
  cat(paste(rep("=", 15 + n_cols * 15), collapse = ""), "\n\n")

  # Determine column width based on component name lengths
  col_width <- max(12, max(nchar(components)) + 2)

  # Print header
  header_fmt <- paste0("%-15s", paste(rep(sprintf(" %%%ds", col_width), n_cols), collapse = ""), "\n")
  cat(do.call(sprintf, c(list(header_fmt, "Parameter"), as.list(components))))
  cat(paste(rep("-", 15 + n_cols * (col_width + 1)), collapse = ""), "\n")

  # Print each parameter row
  for (i in seq_len(nrow(x))) {
    # Print estimate row
    vals <- character(n_cols)
    for (j in seq_along(components)) {
      cell <- x[[components[j]]][i]
      parts <- strsplit(cell, "\n")[[1]]
      vals[j] <- if (length(parts) >= 1) parts[1] else ""
    }
    row_fmt <- paste0("%-15s", paste(rep(sprintf(" %%%ds", col_width), n_cols), collapse = ""), "\n")
    cat(do.call(sprintf, c(list(row_fmt, x$Parameter[i]), as.list(vals))))

    # Print SE row if applicable
    if (show_se) {
      se_vals <- character(n_cols)
      for (j in seq_along(components)) {
        cell <- x[[components[j]]][i]
        parts <- strsplit(cell, "\n")[[1]]
        se_vals[j] <- if (length(parts) >= 2) parts[2] else ""
      }
      cat(do.call(sprintf, c(list(row_fmt, ""), as.list(se_vals))))
    }
  }

  # Print factor parameters
  if (nrow(factor_params) > 0) {
    cat("\n")
    cat(paste(rep("-", 15 + n_cols * (col_width + 1)), collapse = ""), "\n")
    cat("Factor Parameters:\n")
    for (i in seq_len(nrow(factor_params))) {
      est <- round(factor_params$estimate[i], digits)
      se <- round(factor_params$se[i], digits)
      if (show_se) {
        cat(sprintf("  %-20s %.*f (%.*f)\n", factor_params$param_label[i], digits, est, digits, se))
      } else {
        cat(sprintf("  %-20s %.*f\n", factor_params$param_label[i], digits, est))
      }
    }
  }

  # Print n_obs row
  cat("\n")
  cat(paste(rep("-", 15 + n_cols * (col_width + 1)), collapse = ""), "\n")
  row_fmt <- paste0("%-15s", paste(rep(sprintf(" %%%ds", col_width), n_cols), collapse = ""), "\n")
  cat(do.call(sprintf, c(list(row_fmt, "N"), as.list(as.character(n_obs)))))

  # Print log-likelihood
  cat(sprintf("Log-likelihood: %.2f\n", loglik))
  cat("\n")
  if (show_se) {
    cat("Standard errors in parentheses\n")
  }
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1\n")

  invisible(x)
}


#' Export components table to LaTeX
#'
#' Creates a LaTeX table from a factorana_result with components as columns.
#'
#' @param result A factorana_result object from estimate_model_rcpp()
#' @param digits Number of decimal places (default 3)
#' @param file Optional file path to write the LaTeX table
#' @param caption Table caption (optional)
#' @param label Table label for cross-referencing (optional)
#' @param stars Add significance stars (default TRUE)
#' @param note Footnote text (optional, default includes significance codes)
#' @param booktabs Use booktabs package formatting (default TRUE)
#'
#' @return A character string containing the LaTeX table code
#' @export
components_to_latex <- function(result, digits = 3, file = NULL, caption = NULL,
                                 label = NULL, stars = TRUE, note = NULL,
                                 booktabs = TRUE) {
  # Get the components table
  tbl <- components_table(result, digits = digits, show_se = TRUE, stars = stars)
  components <- attr(tbl, "components")
  factor_params <- attr(tbl, "factor_params")
  n_obs <- attr(tbl, "n_obs")
  loglik <- attr(tbl, "loglik")
  n_cols <- length(components)

  # Start building LaTeX
  lines <- character()

  # Table preamble
  lines <- c(lines, "\\begin{table}[htbp]")
  lines <- c(lines, "\\centering")

  if (!is.null(caption)) {
    lines <- c(lines, sprintf("\\caption{%s}", caption))
  }
  if (!is.null(label)) {
    lines <- c(lines, sprintf("\\label{%s}", label))
  }

  # Column specification
  col_spec <- paste0("l", paste(rep("c", n_cols), collapse = ""))
  lines <- c(lines, sprintf("\\begin{tabular}{%s}", col_spec))

  if (booktabs) {
    lines <- c(lines, "\\toprule")
  } else {
    lines <- c(lines, "\\hline")
  }

  # Header row
  header <- paste(c("", components), collapse = " & ")
  lines <- c(lines, paste0(header, " \\\\"))

  if (booktabs) {
    lines <- c(lines, "\\midrule")
  } else {
    lines <- c(lines, "\\hline")
  }

  # Parameter rows
  for (i in seq_len(nrow(tbl))) {
    # Estimate row
    est_vals <- character(n_cols)
    se_vals <- character(n_cols)

    for (j in seq_along(components)) {
      cell <- tbl[[components[j]]][i]
      parts <- strsplit(cell, "\n")[[1]]
      est_vals[j] <- if (length(parts) >= 1 && parts[1] != "") parts[1] else ""
      se_vals[j] <- if (length(parts) >= 2) parts[2] else ""
    }

    est_row <- paste(c(tbl$Parameter[i], est_vals), collapse = " & ")
    lines <- c(lines, paste0(est_row, " \\\\"))

    se_row <- paste(c("", se_vals), collapse = " & ")
    lines <- c(lines, paste0(se_row, " \\\\"))
  }

  # Factor parameters section
  if (nrow(factor_params) > 0) {
    if (booktabs) {
      lines <- c(lines, "\\midrule")
    } else {
      lines <- c(lines, "\\hline")
    }

    lines <- c(lines, sprintf("\\multicolumn{%d}{l}{\\textit{Factor Parameters}} \\\\", n_cols + 1))

    for (i in seq_len(nrow(factor_params))) {
      est <- round(factor_params$estimate[i], digits)
      se <- round(factor_params$se[i], digits)
      # Factor variance spans all columns
      fac_row <- paste(c(factor_params$param_label[i], sprintf("%.*f", digits, est),
                         rep("", n_cols - 1)), collapse = " & ")
      lines <- c(lines, paste0(fac_row, " \\\\"))
      se_row <- paste(c("", sprintf("(%.*f)", digits, se), rep("", n_cols - 1)), collapse = " & ")
      lines <- c(lines, paste0(se_row, " \\\\"))
    }
  }

  # N and Log-likelihood rows
  if (booktabs) {
    lines <- c(lines, "\\midrule")
  } else {
    lines <- c(lines, "\\hline")
  }

  # N row (observations per component)
  n_row <- paste(c("N", as.character(n_obs)), collapse = " & ")
  lines <- c(lines, paste0(n_row, " \\\\"))

  loglik_row <- paste(c("Log-likelihood", sprintf("%.2f", loglik),
                        rep("", n_cols - 1)), collapse = " & ")
  lines <- c(lines, paste0(loglik_row, " \\\\"))

  # Bottom rule
  if (booktabs) {
    lines <- c(lines, "\\bottomrule")
  } else {
    lines <- c(lines, "\\hline")
  }

  lines <- c(lines, "\\end{tabular}")

  # Footnote
  if (is.null(note) && stars) {
    note <- "Standard errors in parentheses. $^{***}p<0.001$, $^{**}p<0.01$, $^*p<0.05$, $^{.}p<0.1$"
  }
  if (!is.null(note)) {
    lines <- c(lines, sprintf("\\\\\\footnotesize{%s}", note))
  }

  lines <- c(lines, "\\end{table}")

  # Combine lines
  latex_code <- paste(lines, collapse = "\n")

  # Write to file if specified
  if (!is.null(file)) {
    writeLines(latex_code, file)
    message(sprintf("LaTeX table written to: %s", file))
  }

  invisible(latex_code)
}


# ---- Internal helper functions ----

#' Build parameter table from result object
#' @keywords internal
.build_param_table <- function(result) {
  estimates <- result$estimates
  std_errors <- result$std_errors
  model_system <- result$model_system

  n_factors <- model_system$factor$n_factors

  # Build the same metadata as initialize_parameters.R
  param_names <- character(0)
  param_types <- character(0)
  components <- character(0)

  # Factor variances (always included, value may be fixed to 1.0 but still in param vector)
  for (k in seq_len(n_factors)) {
    param_names <- c(param_names, sprintf("Var(Factor %d)", k))
    param_types <- c(param_types, "factor_var")
    components <- c(components, "Factor Model")
  }

  # Component parameters - matching initialize_parameters.R logic
  # Key insight: intercept is NOT a separate parameter, it's handled by C++ internally
  for (comp in model_system$components) {
    comp_name <- comp$name

    # Handle multinomial logit with >2 choices
    if (comp$model_type == "logit" && !is.null(comp$num_choices) && comp$num_choices > 2) {
      n_alts <- comp$num_choices - 1

      for (alt in seq_len(n_alts)) {
        # Covariates (no explicit intercept parameter)
        if (!is.null(comp$covariates)) {
          for (cov in comp$covariates) {
            param_names <- c(param_names, sprintf("%s (alt %d)", cov, alt))
            param_types <- c(param_types, "beta")
            components <- c(components, comp_name)
          }
        }

        # Loadings
        if (is.null(comp$loading_normalization)) {
          for (k in seq_len(n_factors)) {
            param_names <- c(param_names, sprintf("Loading %d (alt %d)", k, alt))
            param_types <- c(param_types, "loading")
            components <- c(components, comp_name)
          }
        } else {
          for (k in seq_len(n_factors)) {
            if (is.na(comp$loading_normalization[k]) ||
                abs(comp$loading_normalization[k]) < 1e-10) {
              param_names <- c(param_names, sprintf("Loading %d (alt %d)", k, alt))
              param_types <- c(param_types, "loading")
              components <- c(components, comp_name)
            }
          }
        }
      }
    } else {
      # Standard model types - covariates (NO intercept as separate parameter)
      if (!is.null(comp$covariates)) {
        for (cov in comp$covariates) {
          param_names <- c(param_names, cov)
          param_types <- c(param_types, "beta")
          components <- c(components, comp_name)
        }
      }

      # Loadings (free loadings only, where loading_normalization is NA)
      if (is.null(comp$loading_normalization)) {
        for (k in seq_len(n_factors)) {
          param_names <- c(param_names, sprintf("Loading %d", k))
          param_types <- c(param_types, "loading")
          components <- c(components, comp_name)
        }
      } else {
        for (k in seq_len(n_factors)) {
          if (is.na(comp$loading_normalization[k])) {
            param_names <- c(param_names, sprintf("Loading %d", k))
            param_types <- c(param_types, "loading")
            components <- c(components, comp_name)
          }
        }
      }

      # Sigma for linear models
      if (comp$model_type == "linear") {
        param_names <- c(param_names, "Sigma")
        param_types <- c(param_types, "sigma")
        components <- c(components, comp_name)
      }

      # Thresholds for ordered probit
      if (comp$model_type == "oprobit" && !is.null(comp$num_choices)) {
        n_thresh <- comp$num_choices - 1
        for (j in seq_len(n_thresh)) {
          param_names <- c(param_names, sprintf("Threshold %d", j))
          param_types <- c(param_types, "cutpoint")
          components <- c(components, comp_name)
        }
      }
    }
  }

  # Handle case where names don't match (e.g., if result has different structure)
  n_expected <- length(param_names)
  n_actual <- length(estimates)

  if (n_expected != n_actual) {
    # Fall back to simple numbering if structure doesn't match
    warning(sprintf("Parameter count mismatch (expected %d, got %d). Using simple labels.",
                    n_expected, n_actual))
    param_names <- paste0("param_", seq_len(n_actual))
    param_types <- rep("unknown", n_actual)
    components <- rep("Model", n_actual)
  }

  data.frame(
    component = components,
    param_label = param_names,
    param_type = param_types,
    estimate = estimates,
    se = std_errors,
    stringsAsFactors = FALSE
  )
}
