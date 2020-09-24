#' Data Generation
#' @param beta regression coefficients
#' @param lambda0 baseline hazards
#' @param theta sensitivities
#' @param phi specificities
#' @param n sample size
#' @param seed random seed 
mkdata <- function(beta, lambda0, theta, phi, n, seed = 1) {
  set.seed(seed)
  z <- rnorm(n, 3, 0.1) # rbinom(n, 1, 0.5)
  X <- matrix(z)
  eta <- c(X %*% beta)
  O <- t(outer(1 - lambda0, exp(eta), "^"))
  S <- t(apply(O, 1, cumprod))
  P <- cbind(1, S[, 1:4]) - S
  time <- time0 <- rep(NA, n)
  for (k in 1:n) {
    p <- P[k, ]
    M <- expand.grid(i = 1:5, j = 1:5, w = NA)
    for (r in 1:nrow(M)) {
      i <- M[r, "i"]
      j <- M[r, "j"]
      if (i <= j) {
        term1 <- if (i > 1) prod(phi[1:(i - 1)]) else 1
        term2 <- if (j > i) prod(1 - theta[1:(j - i)]) else 1
        M[r, "w"] <- term1 * term2 * theta[j - i + 1] * p[i]
      } else {
        term1 <- if (j > 1) prod(phi[1:(j - 1)]) else 1
        term2 <- if (i > j) prod(1 - phi[j:(i - 1)]) else 1
        M[r, "w"] <- term1 * term2 * theta[1] * p[i] 
      }
    }
    M <- rbind(M, c(6, 6, 1 - sum(M$w)))
    idx <- sample(1:nrow(M), size = 1, prob = M[, "w"])
    time[k] <- M[idx, "i"]
    time0[k] <- M[idx, "j"]
  }
  c <- 5 # sample(3:5, size = n, replace = TRUE)
  t0 <- pmin(time0, c)
  d0 <- 1 * (time0 <= c)
  data.frame(t0, d0, z)
}
