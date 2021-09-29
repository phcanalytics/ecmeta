#' Create a log hazard ratio data object
#' 
#' Create a data object that stores estimates of log hazard ratios from a survival
#' model or covert an existing object with `as_loghr_data()`.
#' 
#' @param estimate The point estimate of the log hazard ratio.
#' @param standard_error The standard error of the estimate of the 
#' log hazard ratio.
#' 
#' @return An object of class `loghr_data`, which inherits from `data.frame` and 
#' contains one column named `estimate` and a second column named `standard_error`.
#' @seealso [`loghr_data`] objects can also be created from existing `R` objects 
#' with [`as_loghr_data()`].
#' @export
loghr_data <- function(estimate, standard_error) {
  
  if (length(estimate) != length(standard_error)) {
    stop("'estimate' and 'standard_error' must be the same length.")
  }
  
  obj <- data.frame(estimate = estimate,
                    standard_error = standard_error)
  class(obj) <- c("loghr_data", class(obj))
  obj
}

#' Convert to `loghr_data`
#' 
#' Methods to convert an object to a [`loghr_data`] object.
#' @param x An `R` object.
#' @export
as_loghr_data <- function(x, ...) {
  UseMethod("as_loghr_data")
}

#' @param estimate Name of the column containing the point estimate of 
#' the log hazard ratio.
#' @param standard_error Name of the column containing the standard error of 
#' the estimate of the log hazard ratio.
#' @rdname as_loghr_data
#' @export
as_loghr_data.data.frame <- function(x, estimate = "estimate", 
                         standard_error = "standard_error", ...) {
  loghr_data(estimate = x[[estimate]],
             standard_error = x[[standard_error]])
}

