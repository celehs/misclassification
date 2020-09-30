#' Estimation of Accuracy Measures
#' @param Z baseline covariates
#' @param beta regression coefficients
#' @param lambda0 baseline hazards
#' @export
acc_est <- function(Z, beta, lambda0) {
  K <- length(lambda0)
  Lambda0 <- cumsum(lambda0)  
  Surv0 <- exp(-Lambda0)  
  score <- c(as.matrix(Z) %*% beta)
  Surv <- outer(Surv0, exp(score), "^")
  score_cut <- sort(unique(score), decreasing = TRUE)
  # score_cut <- seq(max(score), min(score), length.out = 100)
  I <- 1 * outer(score, score_cut, ">=")
  FPR <- TPR <- matrix(NA, length(score_cut), K)
  AUC <- rep(NA, K)
  for (k in 1:K) {
    w1 <- Surv[k, ] / sum(Surv[k, ])
    w2 <- (1 - Surv[k, ]) / sum(1 - Surv[k, ])
    FPR[, k] <- colSums(w1 * I) 
    TPR[, k] <- colSums(w2 * I)
    AUC[k] <- sum(diff(c(0, FPR[, k])) * TPR[, k])
  }
  list(AUC = AUC, TPR = TPR, FPR = FPR)
}
