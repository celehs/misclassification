#' Simulation Studies
#' @param beta regression coefficients
#' @param lambda0 baseline hazards
#' @param theta sensitivities
#' @param phi specificities
#' @param n sample size
#' @param R number of repetitions
#' @param method optimization algorithm
simfun <- function(beta, lambda0, theta, phi, n, R, method = "BFGS") {
  truth <- c(beta, lambda0, theta, phi)
  start <- c(truth[1], qlogis(truth[-1]))
  M <- matrix(NA, R, length(truth))
  colnames(M) <- c("beta",
    paste0("lambda0_", 1:5),
    paste0("theta_", 0:4),
    paste0("phi_", 1:5))
  conv <- rep(NA, R) # convergence check
  for (r in 1:R) {
    # data generation
    data <- mkdata(
      beta = beta,
      lambda0 = lambda0,
      theta = theta,
      phi = phi,
      n = n,
      seed = r)
    t0 <- data[, 1]
    d0 <- data[, 2]
    X <- as.matrix(data[, -(1:2)])
    if (r == 1) {
      cat("\n", "outcome summary of 1st simulated data:", "\n")
      print(addmargins(table(d0, t0)))
      cat("\n")
    }
    # theta & phi unknown and varying
    fit <- optim(
      par = start,
      fn = function(par) {
        -loglik(
          beta = par[1],
          lambda0 = plogis(par[1 + 1:5]),
          theta = plogis(par[6 + 1:5]),
          phi = plogis(par[11 + 1:5]), 
          t0 = t0,
          d0 = d0,
          X = X)
      },
      method = method)
    # save simulation results
    cat(paste("simulation:", r), "\n")
    conv[r] <- fit$convergence
    M[r, ] <- c(fit$par[1], plogis(fit$par[-1]))
  }
  cat("\n")
  # convergence: An integer code.
  # 0: indicates successful completion
  # 1: indicates that the iteration limit maxit had been reached.
  cat("convergence of optimizations:", "\n")
  cat(paste0("# of successes: ",
             sum(conv == 0), "/", R, "\n\n"))
  M <- M[conv == 0, ]
  print(summary(M))
  DF <- data.frame(
    truth = truth,
    mean = colMeans(M),
    bias = colMeans(M) - truth,
    sd = apply(M, 2, sd))
  return(round(DF, 3))
}
