# Helpers ----------------------------------------------------------------------
cat_mcmc_list <- function(x) {
  mcpar <- coda::mcpar(x[[1]])
  start <- mcpar[1]; end <- mcpar[2]; thin <- mcpar[3]
  
  cat("Iterations = ", start, ":", end, "\n", sep = "")
  cat("Thinning interval =", thin, "\n")
  cat("Number of chains =", coda::nchain(x), "\n")
  cat("Sample size per chain =", (end - start)/thin + 1, "\n\n")
}

cat_predict_message <- function() {
  cat("The estimates of 'mu' and 'sigma' were used to draw the true log HR of the \n",
      "internal control relative to the external control ('loghr_ic_ec'), which was \n",
      "in turn, used to adjust 'log_trt_ec' to obtain the predicted log HR \n",
      "for the treatment relative to a hypothetical internal control ('loghr_trt_ic').",
      sep = "")
}

# Print method for ecmeta_jags -------------------------------------------------
#' @export
print.ecmeta_jags <- function(x, ...) {
  cat_mcmc_list(x$mcmc)
  print(summary(x))
  invisible(x)
}

# Print method for ecmeta_jags_prediction --------------------------------------
#' @export
print.ecmeta_jags_prediction <- function(x, ...) {
  s <- summary(x)
  
  cat("Posterior sampling was performed for the true log HR of the treatment \n",
      "relative to the external control ('loghr_trt_ec') with the following settings:",
      sep = "")
  cat("\n\n")
  cat_mcmc_list(x$mcmc)
  
  cat_predict_message()
  cat("\n\n")
  cat("Summaries of the posterior distributions are displayed below:")
  cat("\n")
  print(s)
  
  invisible(x)
}

# Print method for ecmeta_ml ---------------------------------------------------
#' @export
print.ecmeta_ml <- function(x, ...) {
  cat("Parameters:")
  cat("\n")
  print(summary(x))
  invisible(x)
}

# Print method for ecmeta_ml_prediction ----------------------------------------
#' @export
print.ecmeta_ml_prediction <- function(x, ...) {
  cat_predict_message()
  
  cat("\n\n")
  cat("Summaries of the distributions of the simulated parameters:")
  cat("\n")
  print(summary(x))
  invisible(x)
}