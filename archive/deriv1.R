deriv1 <- function(beta, lambda0, theta, phi, t0, d0, X) {
  gamma0 <- stats::qlogis(lambda0)
  n <- length(t0)
  p <- length(beta)
  q <- length(gamma0)  
  eta <- c(X %*% beta)
  v <- 1 + exp(gamma0)
  M <- t(outer(v, -exp(eta), "^"))
  L <- rep(NA, n)
  S <- matrix(NA, n, p + q)
  colnames(S) <- c(paste0("beta_", 1:p), paste0("gamma0_", 1:q))
  for (i in 1:n) {
    Gamma <- phi^(t0[i] - d0[i]) * (1 - phi)^d0[i] # scalar
    Delta <- phi^(1:t0[i] - 1) * (1 - theta)^(
      t0[i] - 1:t0[i] + 1 - d0[i]) * theta^d0[i] # vector
    # likelihood
    L[i] <- Gamma * prod(M[i, 1:t0[i]]) 
    L[i] <- L[i] + Delta[1] * (1 - M[i, 1]) 
    if (t0[i] > 1) {
      for (j in 2:t0[i]) {
        L[i] <- L[i] + Delta[j] * (1 - M[i, j]) * prod(M[i, 1:(j - 1)])
      }    
    }
    # first derivative (beta)
    for (l in 1:p) {
      S[i, l] <- Gamma * X[i, l] * exp(eta[i]) * 
        prod(M[i, 1:t0[i]]) * sum(log(v[1:t0[i]]))
      S[i, l] <- S[i, l] + Delta[1] * M[i, 1] * log(v[1]) 
      if (t0[i] > 1) {
        for (j in 2:t0[i]) {
          S[i, l] <- S[i, l] + Delta[j] * X[i, l] * exp(eta[i]) * 
            prod(M[i, 1:(j - 1)]) * M[i, j] * sum(log(v[1:j])) - 
            sum(log(v[1:(j - 1)]))
        }
      }
    }
    # first derivative (gamma0)
    for (l in 1:q) { 
      S[i, p + l] <- 0
      if (l <= t0[i]) {
        tmp <- Gamma * prod(M[i, 1:t0[i]]) - Delta[l] * prod(M[i, 1:l])
        if (t0[i] > 1) {
          for (j in (l+1):t0[i]) {
            # tmp <- tmp + Delta[j] * (1 - M[i, j]) * prod(M[, 1:(j - 1)])
          }
        }
        S[i, p + l] <- -exp(eta[i]) * plogis(gamma0[l]) * tmp
      }      
    }    
  }
  return(colSums((1/L) * S))  
}

