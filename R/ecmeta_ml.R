# Optimize log-likelihood ------------------------------------------------------
optimize_loglik <- function(estimate, standard_error, hessian = TRUE) {
  
  v <- standard_error^2
  
  loglik <- function(estimate, mu, sigma, v) {
    0.5 * sum(log(sigma * sigma + v)) + 0.5 * sum((estimate - mu)^2/(sigma * sigma + v))
  }
  
  optim_fun <- function(par) {
    loglik(estimate = estimate, mu = par[1], sigma = exp(par[2]), v = v)
  }
  
  inits <- c(mu = mean(estimate), logsigma = log(sd(estimate)))
  
  params <- stats::optim(par = inits, fn = optim_fun, hessian = hessian)
  params
} 

# Main ecmeta_ml() function ----------------------------------------------------
#' Meta-analysis of external control studies using maximum likelihood
#' 
#' Perform a meta-analysis of studies in which external control arms were 
#' compared to internal control arms by optimizing the log-likelihood function. 
#' You can view the [ecmeta()] documentation for more general details. 
#' 
#' @inheritParams ecmeta
#' @param hessian Logical passed to [stats::optim()] indicating whether a Hessian
#' matrix should be computed. This is required for computation of standard errors 
#' for `mu` and `sigma`, but not needed to simulate predicted log hazard ratios
#' with [predict.ecmeta_ml()].
#' 
#' 
#' @return An `ecmeta_ml` object which is a list containing the following elements:
#' \describe{
#' \item{estimates}{Point estimates of model parameters `mu` and `logsigma` (the
#' log of `sigma`).}
#' \item{vcov}{The variance covariance matrix of the model parameters. Only included
#' if `hessian = TRUE`.}
#' \item{data}{The `data` object.}
#' }
#' 
#' @seealso [ecmeta()]
#' @aliases print.ecmeta_ml
#' @export
ecmeta_ml <- function(data, hessian = TRUE) {
  
  opt <- optimize_loglik(data$estimate, data$standard_error, hessian = hessian)
  if (hessian) vcov <- solve(opt$hessian) else vcov <- NULL
  
  obj <- list(
    estimates = opt$par,
    vcov = vcov,
    data = data
  )
  class(obj) <- "ecmeta_ml"
  obj
  
}