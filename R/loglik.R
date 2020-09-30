#' Log-likelihood
#' @param beta regression coefficients
#' @param lambda0 baseline hazards
#' @param theta sensitivities
#' @param phi specificities
#' @param t0 observed event time
#' @param d0 observed event indicator
#' @param X baseline covaraites
#' @export
loglik <- function(beta, lambda0, theta, phi, t0, d0, X) {
  n <- length(t0)
  eta <- c(as.matrix(X) %*% beta)
  gamma0 <- stats::qlogis(lambda0)
  M <- t(outer(1 + exp(gamma0), -exp(eta), "^"))
  lik <- rep(NA, n)
  for (i in 1:n) {
    phi0 <- phi[t0[i]]
    Gamma <- prod(phi[1:t0[i]]) *
      ((1 - phi0) / phi0)^d0[i]
    theta0 <- theta[t0[i]]
    Delta <- prod(1 - theta[1:t0[i]]) *
      (theta0 / (1 - theta0))^d0[i] # k = 1
    lik[i] <- prod(M[i, 1:t0[i]]) * Gamma
    lik[i] <- lik[i] + (1 - M[i, 1]) * Delta
    if (t0[i] > 1) {
      for (k in 2:t0[i]) {
        theta0 <- theta[t0[i] - k + 1]
        Delta <- prod(phi[1:(k - 1)]) *
          prod(1 - theta[0:(t0[i] - k) + 1]) *
          (theta0 / (1 - theta0))^d0[i]
        lik[i] <- lik[i] +
          prod(M[i, 1:(k - 1)]) *
          (1 - M[i, k]) * Delta
      }
    }
  }
  return(sum(log(lik)))
}
