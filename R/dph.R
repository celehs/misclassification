#' Discrete Proportional Hazards (DPH) Model
#' @param t discrete time (vector)
#' @param d event indicator (vector)
#' @param X covariates (matrix)
#' @export
dph_fit <- function(t, d, X) {
  X <- as.matrix(X)
  n <- nrow(X)
  stopifnot(t > 0)
  stopifnot(d %in% c(0, 1))
  stopifnot(length(t) == n)
  stopifnot(length(d) == n)  
  L <- vector("list", n)
  for (i in 1:n) {
    y <- rep(0, t[i])
    if (d[i] == 1) y[t[i]] <- 1
    L[[i]] <- cbind(y, z = 1:t[i])
  }
  DF <- cbind(do.call("rbind", L), 
              X[rep(seq(t), t), ])
  y <- DF[, 1]
  z <- factor(DF[, 2])
  Z <- DF[, -(1:2)]
  fit <- stats::glm(y ~ z + Z - 1, 
             family = stats::binomial("cloglog"))
  return(fit)
}
