# Generic documentation for summary methods for ecmeta predictions -------------
#' Summary of predicted log hazard ratios
#' 
#' Summary statistics for true log hazard ratios (HRs) predicted using 
#' `ecmeta` [prediction functions][predict.ecmeta].
#' 
#' @param object A object of the appropriate class.
#' @param  exponentiate	Logical indicating whether or not to exponentiate the 
#' the log HR. Defaults to `FALSE`; if `TRUE`, then HRs
#' are summarized.
#' @param ... Additional arguments to pass to [coda::summary.mcmc.list()].
#' 
#' @return A `data.frame` with one row for each log HR comparison:
#' treatment vs. internal control (`loghr_trt_ic`),
#' treatment vs. external control (`loghr_trt_ec`), and
#' internal control vs. external control (`loghr_ic_ec`).
#' 
#' The columns contain summary statistics and quantiles. 
#' Statistics reported are the mean (`mean`), standard deviation (`sd`), 
#' and naive (ignoring autocorrelation across iterations) standard error
#' of the mean. Additional columns contain quantiles of the distribution as 
#' specified by `quantiles`.
#' 
#' 
#' If a Bayesian model was used to draw from the posterior distribution of 
#' `loghr_trt_ec` and multiple chains were used, then the column 
#' `rhat` is included as well. This is the R-hat convergence diagnostic with 
#' values substantially above 1 indicating a lack of convergence. Note that 
#' `rhat` is only provided for `loghr_trt_ec`. Diagnostics for `loghr_ic_ec` should
#'  be based on `mu` and `sigma` from [ecmeta()] and `loghr_trt_ic` is
#'  deterministically related to `loghr_trt_ec` and `loghr_ic_ec`.
#'  
#' @name summary.ecmeta_prediction
NULL

# Summary methods for ecmeta_jags predictions ----------------------------------
#' @rdname summary.ecmeta_prediction
#' @export
summary.ecmeta_jags_prediction <- function(object, exponentiate = FALSE,
                                           quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                           ...) {
  
  if (exponentiate) {
    fun <- exp  
    colnames(object$loghr) <- paste0("hr_", colnames(object$loghr))
  } else{
    fun <- NULL
    colnames(object$loghr) <- paste0("loghr_", colnames(object$loghr))
  }
  
  s_df <- summary_mcmc(coda::as.mcmc(object$loghr),
                       quantiles = quantiles,
                       fun = fun)
  s_df <- s_df[, -which(colnames(s_df) == c("se_ts"))]
  
  if (coda::nchain(object$mcmc) > 1) {
    s_df$rhat <- c(NA, rhat_jags(object$mcmc), NA)
  }
  s_df
}

# Summary methods for ecmeta_ml predictions ------------------------------------
#' @rdname summary.ecmeta_prediction
#' @export
summary.ecmeta_ml_prediction <- function(object, exponentiate = FALSE,
                                         quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                         ...) {
  
  if (exponentiate) {
    fun <- exp  
    colnames(object$loghr) <- paste0("hr_", colnames(object$loghr))
  } else{
    fun <- NULL
    colnames(object$loghr) <- paste0("loghr_", colnames(object$loghr))
  }
  
  # Hack is to convert to MCMC object
  mcmc <- coda::as.mcmc(object$loghr)
  mcmc_summary <- summary_mcmc(mcmc, quantiles = quantiles, fun = fun, ...)
  mcmc_summary$se_ts <- NULL
  
  # Return
  mcmc_summary
}