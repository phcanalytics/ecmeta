as_array_ecmeta <- function(x, name = "mcmc", ...) {
  y <-  as.array(x[[name]], drop = FALSE, ...)
  if (!is.matrix(y)) {
    y <- aperm(y, perm = c(1, 3, 2))
  }
  y
}

#' Convert object to array
#' 
#' This is a helper function to convert a object to a 3-D
#' array as required by [`bayesplot`] in order to diagnosis the 
#' Markov Chain Monte Carlo (MCMC) sampling. See [`bayesplot::MCMC`] for more
#' information about the format of the array. 
#' 
#' @param x An object of the appropriate class.
#' @param ... Additional arguments to pass to [`coda::as.array.mcmc.list()`] other
#' than `drop`, which is set to `FALSE`.
#' 
#' @return A 3-D array in the format described in [`bayesplot::MCMC`].
#' @importFrom coda as.array.mcmc.list
#' @name as.array.ecmeta
#' @export
as.array.ecmeta_jags <- function(x, ...) {
  as_array_ecmeta(x, ...)
}

#' @rdname as.array.ecmeta
#' @export
as.array.ecmeta_jags_prediction <- function(x, ...) {
  as_array_ecmeta(x, name = "loghr", ...)
}

#' @rdname as.array.ecmeta
#' @export
as.array.ecmeta_ml_prediction <- function(x, ...) {
  as_array_ecmeta(x, name = "loghr", ...)
}