#' Prior distributions 
#' 
#' A collection of functions used to specify priors to pass to Bayesian models.
#' @param mean,location Prior for the location parameter. This typically the mean,
#' but for a Cauchy distribution (which is a `student_t` with `df = 1`),
#' it is the median and the mean does not exist.
#' @param sd,scale Prior for the scale parameter, which determines the
#' dispersion of the distribution.
#' @param shape,rate Shape and rate parameters for the inverse gamma distribution.
#' @param min,max Lower and upper limits of the distribution. Must be finite.
#' @name priors
#' @return A named list to be used inside `ecmeta()` and [`predict.ecmeta`]
#' functions.
NULL


#' @rdname priors
#' @export
normal <- function(mean = 0, sd = 10) {
  if (sd < 0) stop("'sd' cannot be negative.")
  list(mean = mean, sd = sd, dist = "normal")
}

#' @rdname priors
#' @export
invgamma <- function(shape = 0.001, rate = 0.001) {
  if (shape <= 0 || rate <= 0) stop("'shape' and 'rate' must be strictly positive.")
  list(shape = shape, rate = rate, dist = "invgamma")
}

#' @rdname priors
#' @export
student_t <- function(location = 0, scale = 25, df = 1) {
  if (location < 0) stop("'location' must be non-negative.")
  if ( scale <= 0 || df <= 0) stop("'scale' and 'df' must be strictly positive.")
  list(location = location, scale = scale, df = df, dist = "student_t")
}

#' @rdname priors
#' @export
uniform <- function(min = 0, max = 1) {
  list(min = min, max = max, dist = "uniform")
}