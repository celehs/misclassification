VTM <- function(vc, dm) {
  matrix(vc, ncol = length(vc), nrow = dm, byrow = T)
}


resam1 <- function(I, initial, X, t0, d0) {

  #### parameter estimation
  ffun.re <- function(initial) {
    beta <- initial[1]
    lambda0 <- plogis(initial[2:6])
    eta <- c(X %*% beta)
    O <- t(outer(1 - lambda0, exp(eta), "^"))
    S <- t(apply(O, 1, cumprod))
    P <- cbind(1, S[, 1:4]) - S

    lambda0_s <- plogis(initial[7:11])
    A <- diag(lambda0_s[1], 5)
    for (k in 1:4) {
      for (kk in (k + 1):5) {
        A[kk, k] <- prod((1 - lambda0_s)[1:(kk - k)]) * lambda0_s[kk - k + 1]
      }
    }

    P0 <- A %*% t(P)
    S0 <- 1 - apply(P0, 2, cumsum)
    P0 <- t(P0)
    S0 <- t(S0)
    index <- cbind(1:n, t0)
    loglik <- I * (d0 * log(P0[index]) + (1 - d0) * log(S0[index]))
    -sum(loglik)
  }

  tryCatch(
    {
      temp <- optim(initial, ffun.re, method = "BFGS")
      code.re <- temp$convergence
      parhat.re <- temp$par

      #### accuracy measures
      betahat <- parhat.re[1]
      lambda0hat <- plogis(parhat.re[2:6])
      eta <- c(X %*% betahat)
      O <- t(outer(1 - lambda0hat, exp(eta), "^"))
      S <- t(apply(O, 1, cumprod))
      P <- cbind(1, S[, 1:4]) - S
      # range(apply(P,1,sum))

      F <- 1 - S
      zsort <- sort(z)
      auc.hat <- rep(0, 5)
      TPR.hat <- matrix(0, nrow = 5, ncol = 2)
      FPR.hat <- matrix(0, nrow = 5, ncol = 2)
      for (k in 1:5) {
        tpr <- apply(as.numeric(I) * F[, k] * (t(VTM(z, n)) >= VTM(zsort, n)), 2, sum) / sum(as.numeric(I) * F[, k])
        fpr <- apply(as.numeric(I) * S[, k] * (t(VTM(z, n)) >= VTM(zsort, n)), 2, sum) / sum(as.numeric(I) * S[, k])
        auc.hat[k] <- -sum(tpr[1:(n - 1)] * diff(fpr))
        TPR.hat[k, ] <- unique(tpr)
        FPR.hat[k, ] <- unique(fpr)
      }
      auchat.re <- auc.hat
      tpr.re <- TPR.hat[, 2]
      fpr.re <- FPR.hat[, 2]

      out <- c(code.re, parhat.re[1], plogis(parhat.re[-1]), tpr.re, fpr.re)

      out
    },
    error = function(e) {
      return(rep(NA, 17))
    }
  )
}


resam2 <- function(I, initial, X, t0, d0) {

  #### parameter estimation
  ffun.re <- function(initial) {
    beta <- initial[1]
    lambda0 <- plogis(initial[2:6])
    eta <- c(X %*% beta)
    O <- t(outer(1 - lambda0, exp(eta), "^"))
    S <- t(apply(O, 1, cumprod))
    P <- cbind(1, S[, 1:4]) - S

    lambda0_s <- plogis(initial[7:11])
    A <- diag(lambda0_s[1], 5)
    for (k in 1:4) {
      for (kk in (k + 1):5) {
        A[kk, k] <- prod((1 - lambda0_s)[1:(kk - k)]) * lambda0_s[kk - k + 1]
      }
    }

    P0 <- A %*% t(P)
    S0 <- 1 - apply(P0, 2, cumsum)
    P0 <- t(P0)
    S0 <- t(S0)
    index <- cbind(1:n, t0)
    loglik <- I * (d0 * log(P0[index]) + (1 - d0) * log(S0[index]))
    -sum(loglik)
  }

  tryCatch(
    {
      temp <- optim(initial, ffun.re, method = "BFGS")
      code.re <- temp$convergence
      parhat.re <- temp$par

      #### accuracy measures
      betahat <- parhat.re[1]
      lambda0hat <- plogis(parhat.re[2:6])
      eta <- c(X %*% betahat)
      O <- t(outer(1 - lambda0hat, exp(eta), "^"))
      S <- t(apply(O, 1, cumprod))
      P <- cbind(1, S[, 1:4]) - S
      # range(apply(P,1,sum))

      F <- 1 - S
      zsort <- sort(z)
      auc.hat <- rep(0, 5)
      for (k in 1:5) {
        tpr <- apply(as.numeric(I) * F[, k] * (t(VTM(z, n)) >= VTM(zsort, n)), 2, sum) / sum(as.numeric(I) * F[, k])
        fpr <- apply(as.numeric(I) * S[, k] * (t(VTM(z, n)) >= VTM(zsort, n)), 2, sum) / sum(as.numeric(I) * S[, k])
        auc.hat[k] <- -sum(tpr[1:(n - 1)] * diff(fpr))
      }
      auchat.re <- auc.hat

      out <- c(code.re, parhat.re[1], plogis(parhat.re[-1]), auchat.re)

      out
    },
    error = function(e) {
      return(rep(NA, 17))
    }
  )
}
