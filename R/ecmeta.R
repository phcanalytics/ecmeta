# ecmeta() ---------------------------------------------------------------------
#' Meta-analysis of external control studies
#' 
#' Perform a meta-analysis of studies in which external control arms were
#' compared to internal control arms. For more information on the methodology
#' please see `vignette("methodology")`. 
#' 
#' @param data A [`loghr_data`] object that stores estimates of the log hazard ratios
#' of the internal control arms relative to the external control arms.
#' @param method Method/software used to estimate the meta-analytic model. `jags`
#' uses Markov Chain Monte Carlo to sample from the posterior distribution of
#' the parameters while `ml` using maximum likelihood to fit the model. 
#' @param ... Additional arguments to pass to the underlying meta-analytic estimation
#' function. If `method == "jags"`, these are arguments to pass to [ecmeta_jags()].
#' Unused if `method == "ml`.
#' @seealso See `vignette("guide")` for an example.
#' 
#' @return Returns an object with class based on `method`. If `method == "jags"`,
#' then an [`ecmeta_jags`] object is returned; if `method == "ml`, then
#' an [`ecmeta_ml`] object is returned. 
#' @export
ecmeta <- function(data, method = c("jags", "ml"), ...) {
  
  if (!inherits(data, "loghr_data") ){
    stop("'data' must be an object of class 'loghr_data'.")
  }
  
  # Implementation according to selected method
  method <- match.arg(method)
  if (method == "jags") {
    return(ecmeta_jags(data = data, ...))
  } else {
    return(ecmeta_ml(data = data, ...))
  }
}