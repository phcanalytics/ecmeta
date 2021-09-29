#' R-hat convergence diagnostic
#' 
#' This is a generic method for computing the R-hat convergence diagnostic 
#' from multiple chains using [coda::gelman.diag()]. Values substantially above 1
#' indicate a lack of convergence. 
#' 
#' @param An [`coda::mcmc.list`] object with more than one chain as required
#' by [coda::gelman.diag()].
#' 
#' @return A vector containing R-hat for each parameter.
rhat_jags <- function(x) {
  if (coda::nchain(x) > 1) {
    gelman_diag <- coda::gelman.diag(x)
    return(gelman_diag$psrf[, "Point est."])
  } else {
    stop("Requires more than one chain.", call. = FALSE)
  }
}