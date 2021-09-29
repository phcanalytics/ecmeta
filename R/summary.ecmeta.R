# Helper functions -------------------------------------------------------------
ci_alpha <- function(level) {
  if (level > 1 | level < 0){
    stop("'level' must be in the interval (0,1)",
         call. = FALSE)
  }
  lower <- (1 - level)/2
  upper <- 1 - lower
  return(list(lower = lower, upper = upper))
}

transform_posterior <- function(object, ...) {
  UseMethod("transform_posterior")
}

transform_posterior.mcmc <- function(object, fun, ...) {
  fun(object)
}

summary_mcmc <- function(object, quantiles, fun = NULL, ...) {
  
  # Possibly transform the object
  if (!is.null(fun)) object <- transform_posterior(object, fun = fun)
  
  # First summarize mcmc.list object
  s <- summary(object, quantiles = quantiles, ...)
  
  # Handle cases where there is only one variable
  if (!is.matrix(s$statistics)) {
    varname <- colnames(object[[1]])
    stat <- matrix(s$statistics, nrow = 1)
    colnames(stat) <- names(s$statistics)
    rownames(stat) <- varname
    s$statistics <- stat
  }
  
  if (!is.matrix(s$quantiles)) {
    if (length(quantiles == 1)) {
      quant <- matrix(s$quantiles, ncol = 1)
    } else{
      quant <- matrix(s$quantiles, nrow = 1)
    }
    colnames(quant) <- paste0(quantiles * 100, "%")
    s$quantiles <- quant
  }
  
  
  # Create a data frame
  s_df <- data.frame(
    param = rownames(s$statistics),
    mean = s$statistics[, "Mean"],
    se_naive = s$statistics[, "Naive SE"],
    se_ts = s$statistics[, "Time-series SE"],
    sd = s$statistics[, "SD"]
  )
  s_df <- cbind(
    s_df,
    s$quantiles
  )
  
  # Add rhat
  if (coda::nchain(object) > 1) s_df$rhat <- rhat_jags(object)
  
  # Return
  rownames(s_df) <- NULL
  s_df
}


# Summary methods for ecmeta_jags ----------------------------------------------
#' Summary of external control meta-analysis using JAGS
#' 
#' Summary statistics for a meta-analytic model fit using [ecmeta()] with 
#' `method == "jags"`. 
#' 
#' @param object An [`ecmeta_jags`][ecmeta_jags()] object.
#' @param quantiles A vector of quantiles to evaluate for each variable.
#' @param ... Additional arguments to pass to [coda::summary.mcmc.list()].
#' 
#' @return A `data.frame` with one row for each parameter (`mu` and `sigma`).
#' The columns contain the same summary statistics and quantiles as in
#' [coda::summary.mcmc.list()]. Statistics reported are the mean (`mean`), 
#' standard deviation (`sd`), naive (ignoring autocorrelation of the chain) 
#' standard error of the mean (`se_mean`), and time-series standard error
#' based on an estimate of the spectral density at 0 (`se_ts`). Additional 
#' columns contain quantiles of the distribution as specified by `quantiles`.
#' Finally, if multiple chains were used, the column `rhat` is included, which
#' is the R-hat convergence diagnostic (values substantially above 1 indicate
#' a lack of convergence).
#' @export
summary.ecmeta_jags <- function(object, 
                                quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                ...) {
  summary_mcmc(object$mcmc, quantiles = quantiles, ...)
}

# Summary methods for ecmeta_ml ------------------------------------------------
#' Summary of external control meta-analysis using maximum likelihood
#' 
#' Summary statistics for a meta-analytic model fit using [ecmeta()] with 
#' `method == "ml"`. 
#' 
#' @param object An [`ecmeta_ml`][ecmeta_ml()] object.
#' @param level The confidence level. Default is `0.95` for a 95 percent 
#' confidence interval.
#' @param ... Currently unused.
#' 
#' @note `object` contains estimates of `logsigma` rather than `sigma`. Point
#' estimates and confidence intervals for `sigma` are computed by simply 
#' exponentiating `logsigma`; the standard error for `sigma` is estimated
#' using the delta method. 
#' 
#' @return A `data.frame` with one row for each parameter (`mu` and `sigma`).
#' There is always one column named `estimate` containing the point estimates. 
#' In addition, if a variance-covariance matrix was estimated, then there are 
#' additional columns for the standard errors as well as the lower and upper 
#' confidence limits.
#' @export
summary.ecmeta_ml <- function(object, level = 0.95, ...) {
  
  # Point estimates
  # Transform logsigma to sigma
  estimates <- object$estimates
  estimates["logsigma"] <- exp(estimates["logsigma"])
  names(estimates)[names(estimates) == "logsigma"] <- "sigma"
  
  # Store in data frame
  obj <- data.frame(param = names(estimates),
                    estimate = estimates)
  
  # Standard errors and confidence intervals
  if (!is.null(object$vcov)) {
    v <- diag(object$vcov) # variance
    
    ## Use delta method to get standard error for exp(logsigma)
    grad <- estimates["sigma"] # gradient is just exp(logsigma)
    v_sigma <- c(sigma = grad %*% v["logsigma"] %*% grad)
    obj$se <- c(sqrt(v["mu"]), sqrt(v_sigma))
    
    ## Compute confidence intervals
    if (level > 1 | level < 0){
      stop("'level' must be in the interval (0,1)",
           call. = FALSE)
    }
    alpha <- ci_alpha(level)
    obj$lower <- object$estimates + qnorm(alpha$lower) * sqrt(v)
    obj$upper <- object$estimates + qnorm(alpha$upper) * sqrt(v)
    
    ### Need to exponentiate CI for logsigma
    obj[obj$param == "sigma", "lower"] <- exp(obj[obj$param == "sigma", "lower"])
    obj[obj$param == "sigma", "upper"] <- exp(obj[obj$param == "sigma", "upper"])
    
    ### Rename lower/upper to make quantiles explicit
    colnames(obj)[colnames(obj) %in% c("lower", "upper")] <- 
      paste0(100 * unlist(alpha), "%")
  }
  
  rownames(obj) <- NULL
  obj
}