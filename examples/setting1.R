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

mkdata1 <- function(n, beta, lambda0, lambda0_s) {
  z <- 2 * rbinom(n, 1, 0.5) # set1 #runif(n, 5, 6)
  # lambda0 <-0.05 * (1:5)#rep(.01, 7)
  # z <- rnorm(n, 2, 1) # set2
  X <- matrix(z)
  eta <- c(X %*% beta)
  O <- t(outer(1 - lambda0, exp(eta), "^"))
  S <- t(apply(O, 1, cumprod))
  P <- cbind(1, S[, 1:4]) - S
  # range(apply(P,1,sum))
  
  A <- diag(lambda0_s[1], 5)
  for (k in 1:4) {
    for (kk in (k + 1):5) {
      A[kk, k] <- prod((1 - lambda0_s)[1:(kk - k)]) * lambda0_s[kk - k + 1]
    }
  }
  
  # ## true auc
  # F=1-S
  # zsort=sort(z)
  # auc=rep(0,5); TPR=matrix(0,nrow=5,ncol=2); FPR=matrix(0,nrow=5,ncol=2)
  # for (k in 1:5){
  # tpr=apply(F[,k]*(t(VTM(z,n))>=VTM(zsort,n)),2,sum)/sum(F[,k])
  # fpr=apply(S[,k]*(t(VTM(z,n))>=VTM(zsort,n)),2,sum)/sum(S[,k])
  # auc[k]=-sum(tpr[1:(n-1)]*diff(fpr))
  # TPR[k,]=unique(tpr)
  # FPR[k,]=unique(fpr)
  # }
  
  time0 <- rep(NA, n)
  for (k in 1:n) {
    p <- P[k, ]
    p0 <- A %*% matrix(p, ncol = 1)
    M <- expand.grid(l = 1:6, w = NA)
    M[, "w"] <- c(p0, 1 - sum(p0))
    idx <- sample(1:nrow(M), size = 1, prob = M[, "w"])
    time0[k] <- M[idx, "l"]
  }
  c <- 5 # rdunif(n, 3, 5) #runif(n, 3, 5) # censoring time
  d0 <- 1 * (time0 <= c)
  t0 <- apply(cbind(time0, c), 1, min)
  # print(addmargins(table(d0, t0)))
  
  dplyr::as_tibble(cbind(t0 = t0, d0 = d0, z = z))
}

ffun <- function(initial) {
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
  loglik <- d0 * log(P0[index]) + (1 - d0) * log(S0[index])
  -sum(loglik)
}

model_fit <- function(t0, d0, X, initial) {
  tryCatch(
    {
      temp <- optim(initial, ffun, method = "BFGS")
      code <- temp$convergence
      parhat <- temp$par
    },
    error = function(e) {}
  )
  return(temp)
}
