# Generic documentation for predictions ----------------------------------------
#' Predict log hazard ratio
#' 
#' Predict the true log hazard ratio in a new single-arm study by adjusting
#' for additional bias and variability caused by the non-randomized design. This
#' function first draws from the distribution of the true log hazard 
#' ratio comparing the treatment to the external control (`trt_ec`). It then
#' draws the true log hazard ratio comparing the internal control to the external
#' control (`ic_ec`) from its predictive distribution. Lastly, `trt_ec` and 
#' `ic_ec` are combined to produce simulated draws of the true log hazard ratio 
#' in a comparison of the treatment to the internal control (`trt_ic`).
#' 
#' 
#' @details
#' The implementation differs slightly between the Bayesian and maximum likelihood 
#' approaches. First, when using a Bayesian approach, MCMC is use to sample
#' `trt_ec` whereas `trt_ec` is sampled from a normal distribution based on
#' the log hazard ratio estimates and standard errors in the new study. In 
#' practice, these approaches will produce very similar results.
#' 
#' Second, when using the Bayesian approach, the posterior predictive 
#' distribution of `ic_ec` is drawn using the posterior samples of 
#' `mu` and `sigma` stored in `object`. In the maximum likelihood approach,
#' `ic_ec` is simulated by using the point estimates of `mu` and `sigma` and
#' drawing from a t distribution, which is the predictive distribution for a 
#' future observation.
#' 
#' @param object An object of the appropriate class.
#' @param newdata A `loghr_data` object that stores estimate of the log hazard 
#' ratio for a comparison of the treatment arm to the external control arm in
#' the new single-arm study.
#' @param n_sims Number of simulations to use. Only relevant when using a
#' maximum likelihood based approach.
#' @inheritParams ecmeta_jags
#' @param ... Currently unused.
#' 
#' @seealso See `vignette("methodology")` for a description of the method.
#' 
#' @return A list that may contain the following elements:
#' \describe{
#' \item{loghr}{A matrix with three columns containing draws of the true log 
#' hazard ratios. The columns are: `trt_ic` (a comparison of the treatment
#' to the internal control), `trt_ec` (a comparison of the treatment to the 
#' external control) and `ic_ec` (a comparison of the internal control to 
#' the external control).}
#' \item{mcmc}{When using a Bayesian approach, a [`coda::mcmc.list`] object is 
#' also included that contains posterior samples of the true log hazard ratio 
#' comparing the treatment to the external control. This is not included when 
#' using a maximum likelihood approach.}
#' }
#' @aliases print.ecmeta_jags_prediction 
#' @name predict.ecmeta
NULL

# Prediction method for ecmeta_jags object -------------------------------------
#' @rdname predict.ecmeta
#' @export
predict.ecmeta_jags <- function(object, newdata, 
                                n_burnin = 0, thin = 1, n_adapt = 1000, 
                                quiet = FALSE, ...) {
  
  if (!inherits(newdata, "loghr_data") ){
    stop("'newdata' must be an object of class 'loghr_data'.")
  }
  
  
  # Extract some settings
  n_iter <- coda::niter(object$mcmc) * thin
  n_chains <- coda::nchain(object$mcmc)
  
  # Posterior samples for true log HR of treatment vs. external control
  ## Initial values
  model_inits <- list(
    new_truth = newdata$estimate
  )
  
  ## Model description
  model_description <- "
    model {
    
    new_estimate ~ dnorm(new_truth, 1/new_standard_error^2) # data = {new_estimate, new_standard_error}
    
    new_truth ~ dnorm(0 , 0.0001)
    }
  "
  
  ## Model data
  model_data <- list(
    new_estimate = newdata$estimate,
    new_standard_error = newdata$standard_error
  )
  
  ## JAGS model object
  jags_model <- rjags::jags.model(
    file = textConnection(model_description),
    data = model_data,
    inits = model_inits,
    n.chains = n_chains,
    n.adapt = n_adapt,
    quiet = quiet
  )
  
  ## Burn-in
  update_jags(jags_model, n_burnin, quiet)
  
  ## Sample from the posterior
  model_samples <- rjags::coda.samples(
    jags_model,
    variable.names = c("new_truth"),
    n.iter = n_iter,
    thin = thin
  )
  coda::varnames(model_samples) <- "loghr_trt_ec"
  
  # Posterior samples for true log HR of treatment vs. internal control
  post_ic_ec <- as.matrix(object$mcmc)
  post_trt_ec <- as.matrix(model_samples)
  
  true_loghr_ic_ec_new <- rnorm(nrow(post_ic_ec), 
                                mean = post_ic_ec[, "mu"],
                                sd = post_ic_ec[, "sigma"]) 
  true_loghr_trt_ic_new <- c(post_trt_ec[, "loghr_trt_ec"]) - true_loghr_ic_ec_new
  true_loghr_new <- matrix(
    c(true_loghr_trt_ic_new, 
      post_trt_ec[, "loghr_trt_ec"],
      true_loghr_ic_ec_new), 
    nrow = n_iter * n_chains, ncol = 3
  )
  colnames(true_loghr_new) <- c("trt_ic", "trt_ec", "ic_ec")
  
  # Object
  obj <- list(
    mcmc = model_samples,
    loghr = true_loghr_new
  )
  class(obj) <- "ecmeta_jags_prediction"
  obj
}

# Prediction of ecmeta_ml object -----------------------------------------------
#' @rdname predict.ecmeta
#' @aliases print.ecmeta_ml_prediction 
#' @export
predict.ecmeta_ml <- function(object, newdata, n_sims = 1000, ...) {
  
  # (1) log HR  comparing internal vs. external control
  # The posterior predictive distribution is student t
  n_studies <- nrow(object$data)
  posterior_predictive_rng <- function(n_sims, n, mu, sd) {
    mu + sqrt(1 + 1/n) * sd * rt(n_sims, df = n - 1)
  }
  loghr_ic_ec <- posterior_predictive_rng(n_sims, n =  n_studies,
                                          mu = object$estimates["mu"],
                                          sd = exp(object$estimates["logsigma"]))
  
  
  # (2) log HR  comparing treatment to external control
  loghr_trt_ec <- rnorm(n_sims, newdata$estimate, newdata$standard_error)
  
  # (3) log HR comparing true log HR comparing treatment to internal control
  loghr_trt_ic <- loghr_trt_ec - loghr_ic_ec
  
  # Return
  loghr <- matrix(c(loghr_trt_ic, loghr_trt_ec, loghr_ic_ec),
                   ncol = 3)
  colnames(loghr) <- c("trt_ic", "trt_ec", "ic_ec")
  object <- list(loghr = loghr)
  class(object) <- "ecmeta_ml_prediction"
  object
}
