update_jags <- function(object, n_burnin, quiet) {
  if (quiet) progress_bar <- "text" else progress_bar <- "none"
  if (n_burnin > 0) {
    update(object, n.iter = n_burnin, progress.bar = progress_bar)
  }
}

# ecmeta_jags() main function --------------------------------------------------
#' Meta-analysis of external control studies using JAGS
#' 
#' Perform a meta-analysis of studies in which external control arms were 
#' compared to internal control arms using [`rjags`]. You can view the 
#' [ecmeta()] documentation for more general details. 
#' 
#' @inheritParams ecmeta
#' @param n_chains The number of parallel chains for the model. Passed to 
#' [rjags::jags.model()].
#' @param n_iter Number of iterations to monitor. Passed to [rjags::coda.samples()].
#' @param thin Thinning interval for monitors. Passed to [rjags::coda.samples()].
#' @param n_burnin Number of 
#' @param n_adapt The number of iterations for adaptation. Passed to [rjags::jags.model()].
#' @param quiet If `TRUE` then messages generated during compilation will be suppressed, 
#' as well as the progress bar during adaptation. Passed to [rjags::jags.model()].
#' @param prior_mu Prior for `mu`. Must be specified as a normal distribution with 
#' [normal()].
#' @param prior_scale Prior for the scale parameter, which is either `sigma` or 
#' `sigma^2` depending on the choice of prior. A half [student_t()] distribution
#' or uniform prior can be placed on `sigma` while an [invgamma()] prior can be
#' used for `sigma^2`.
#' 
#' @return An `ecmeta_jags` object, which is a list containing the following elements:
#' \describe{
#' \item{mcmc}{A [`coda::mcmc.list`] object containing samples of the model parameters.}
#' \item{data}{The `data` object.}
#' }
#' 
#' @aliases print.ecmeta_jags 
#' @export
ecmeta_jags <- function(data, n_chains = 1, n_iter = 1000,
                        thin = 1, n_burnin = 0, n_adapt = 1000, quiet = FALSE,
                        prior_mu = normal(0, 100),
                        prior_scale = invgamma(0.001, 0.001)) {
  
  if (!"rjags" %in% rownames(installed.packages())) {
    stop("You must install 'rjags'.")
  }
  
  # Initial values
  model_inits <- list(
    mu = mean(data$estimate)
  )
  if (prior_scale$dist == "invgamma") {
    model_inits$tau = 1/sd(data$estimate)
  } else if (prior_scale$dist == "cauchy") {
    model_inits$sigma = sd(data$estimate)
  }
  
  # Describe JAGS model
  ## Set priors
  prior_mu_jags <- paste0("mu ~ dnorm(", prior_mu$mean, ",", 
                          1/prior_mu$sd^2, ")")
  if (prior_scale$dist == "invgamma") {
    # tau ~ gamma(shape, rate)
    prior_tau_jags <- paste0("tau ~ dgamma(", prior_scale$shape, ",", 
                              prior_scale$rate, ")")
    prior_scale_jags <- paste0(
      prior_tau_jags,
      "sigma = 1.0/sqrt(tau)"
    )
  } else if (prior_scale$dist == "student_t")  {
    # sigma ~ dt(location, scale, 1)T(0, )
    prior_sigma_jags <- paste0(
      "sigma ~ dt(", 
      prior_scale$location, ",", 
      1/prior_scale$scale^2, ",",
      prior_scale$df, ")",
      "T(0,)"
    )
    prior_scale_jags <- paste0(
      prior_sigma_jags,
      "tau = 1.0/(sigma^2)"
    )
  } else if (prior_scale$dist == "uniform") {
    prior_sigma_jags <- paste0(
      "sigma ~ dunif(", 
      prior_scale$min, ",", 
      prior_scale$max,
      ")"
    )
    prior_scale_jags <- paste0(
      prior_sigma_jags,
      "tau = 1.0/(sigma^2)"
    )
  } else{
    stop("Distribution used for 'prior_scale' is not supported.")
  }

  ## Describe model
  model_description <- paste0(
    "model {
      for(i in 1:n_ref) { #data = {n_ref}
        truth[i] ~ dnorm(mu, tau)
        estimate[i] ~ dnorm(truth[i] , 1/standard_error[i]^2) #data = {estimate, standard_error}
      }",
    prior_mu_jags,
    prior_scale_jags,
   "}"
  )
  
  # Model data
  model_data <- list(
    n_ref = length(data$estimate),
    standard_error = data$standard_error,
    estimate = data$estimate
  )
  
  # JAGS model object
  jags_model <- rjags::jags.model(
    file = textConnection(model_description),
    data = model_data,
    inits = model_inits,
    n.chains = n_chains,
    n.adapt = n_adapt,
    quiet = quiet
  )
  
  # Burn-in
  update_jags(jags_model, n_burnin, quiet)
  
  # Sample from the posterior
  model_samples <- rjags::coda.samples(
    jags_model,
    variable.names = c("mu", "sigma"),
    n.iter = n_iter,
    thin = thin
  )
  
  # Return 
  obj <- list(
    mcmc = model_samples,
    data = data
  )
  class(obj) <- "ecmeta_jags"
  obj
}

