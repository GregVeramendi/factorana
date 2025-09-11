#' Define a model system
#'
#' @param components A named list of model_component objects
#'
#' @return An object of class "model_system". A list of model_component objects and one factor_model object.
#' @export
define_model_system <- function(components, factor) {
  # Validate the inputs:

  if (!is.list(components) || !all(sapply(components, inherits, "model_component"))) {
    stop("Input must be a list of model_component objects.")
  }

  if (!is.null(names(components))){
    stop("model_component objects must be named")
  }

  if (!inherits(factor, "factor_model")) {
    stop("`factor` must be of class 'factor_model'")
  }


  if (!all(sapply(components, function(x) inherits(x, "model_component")))) {
    stop("All elements must be of class 'model_component'")
  }

  #check if each component has the same factor model as factor

 if (!all(vapply(components, function(comp) identical(get_factor(comp), factor), logical(1)))) {
   stop("All model_components must have Factor as the factor_model")
 }


  out <- list(
    components = components,
    factor = factor)
  class(out) <- "model_system"
  return(out)
}

#' @export
print.model_system <- function(x, ...) {
  stopifnot(inherits(x, "model_system")) #check that x is a class model_system
  comps <- x$components
  n <- length(comps)
  nms <- names(comps) #captures the list names (check if there are errors in doing it this way)

  cat("Model System\n")
  cat("------------\n")
  cat("Components:", n, "\n")

  if (!n) return(invisible(x))

  for (i in seq_along(comps)) { #loops over indices, seq_along safer for empty vectors
    # header for this component
    label <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else paste0("<unnamed-", i, ">") #get label for component, otherwise call it unnamed - number
    cat("\n[", i, "] ", label, "\n", sep = "") #print section header
    cat(strrep("-", 2 + nchar(label) + nchar(as.character(i))), "\n", sep = "") #print underline

    # delegate to the component's own print method
    print(comps[[i]], ...)   # passes ... along to the model_component's own print method
  }

  invisible(x) #if 0 components, stop printing and return x invisibly
}

