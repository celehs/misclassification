options(warn = -1)
rm(list = ls())
timestart <- Sys.time()
# setwd("/Users/xuanwang/Desktop/misclassification")
setwd("/home/xw127/")
# library("fitur")
source("funs1009.R")

rep.proc <- 5
code <- rep(NA, rep.proc)
parhat <- matrix(0, ncol = rep.proc, nrow = 1 + 5 + 5)
auchat <- matrix(0, ncol = rep.proc, nrow = 5)
tpr.hat <- matrix(0, ncol = rep.proc, nrow = 5)
fpr.hat <- matrix(0, ncol = rep.proc, nrow = 5)
convergeno <- rep(NA, rep.proc)
parse <- matrix(0, ncol = rep.proc, nrow = 1 + 5 + 5 + 5 + 5)
code.m <- rep(NA, rep.proc)
parhat.m <- matrix(0, ncol = rep.proc, nrow = 1 + 5)

for (ll in 1:rep.proc) {
  # ll=1
  # set.seed(ll)
  n <- 1000
  beta <- 1
  lambda0 <- 0.1 * (1:5) # rep(.01, 7)
  z <- 2 * rbinom(n, 1, 0.5) # set1 #runif(n, 5, 6)
  # lambda0 <-0.05 * (1:5)#rep(.01, 7)
  # z <- rnorm(n, 2, 1) # set2
  X <- matrix(z)
  eta <- c(X %*% beta)
  O <- t(outer(1 - lambda0, exp(eta), "^"))
  S <- t(apply(O, 1, cumprod))
  P <- cbind(1, S[, 1:4]) - S
  # range(apply(P,1,sum))

  lambda0_s <- 0.1 * (8:4)
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

  #### parameter estimation
  truth <- c(beta, lambda0, lambda0_s)
  ini <- truth #* ( 1+rnorm(length(truth),0,0.2^2) )
  initial <- c(ini[1], qlogis(ini[-1]))
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
  tryCatch(
    {
      temp <- optim(initial, ffun, method = "BFGS")
      code[ll] <- temp$convergence
      parhat[, ll] <- temp$par
    },
    error = function(e) {}
  )

  #### accuracy measures
  betahat <- parhat[1, ll]
  lambda0hat <- plogis(parhat[2:6, ll])
  eta <- c(X %*% betahat)
  O <- t(outer(1 - lambda0hat, exp(eta), "^"))
  S <- t(apply(O, 1, cumprod))
  P <- cbind(1, S[, 1:4]) - S
  # range(apply(P,1,sum))

  # set1
  F <- 1 - S
  zsort <- sort(z)
  auc.hat <- rep(0, 5)
  TPR.hat <- matrix(0, nrow = 5, ncol = 2)
  FPR.hat <- matrix(0, nrow = 5, ncol = 2)
  for (k in 1:5) {
    tpr <- apply(F[, k] * (t(VTM(z, n)) >= VTM(zsort, n)), 2, sum) / sum(F[, k])
    fpr <- apply(S[, k] * (t(VTM(z, n)) >= VTM(zsort, n)), 2, sum) / sum(S[, k])
    auc.hat[k] <- -sum(tpr[1:(n - 1)] * diff(fpr))
    TPR.hat[k, ] <- unique(tpr)
    FPR.hat[k, ] <- unique(fpr)
  }
  auchat[, ll] <- auc.hat
  tpr.hat[, ll] <- TPR.hat[, 2]
  fpr.hat[, ll] <- FPR.hat[, 2]

  # # set2
  # F=1-S
  # zsort=sort(z)
  # auc.hat=rep(0,5); TPR.hat=matrix(0,nrow=5,ncol=2); FPR.hat=matrix(0,nrow=5,ncol=2)
  # for (k in 1:5){
  #   tpr=apply(F[,k]*(t(VTM(z,n))>=VTM(zsort,n)),2,sum)/sum(F[,k])
  #   fpr=apply(S[,k]*(t(VTM(z,n))>=VTM(zsort,n)),2,sum)/sum(S[,k])
  #   auc.hat[k]=-sum(tpr[1:(n-1)]*diff(fpr))
  # }
  # auchat[,ll]=auc.hat

  ######## variance estimation
  # set1
  resample <- 500

  I <- matrix(rexp(n * resample), nrow = n)

  parre <- apply(I, 2, resam1, initial, X, t0, d0)

  convergeno[ll] <- sum(parre[1, ] == 0, na.rm = T)
  parse[, ll] <- apply(parre[-1, parre[1, ] == 0], 1, sd, na.rm = T)

  # # set2
  # resample=500
  #
  # I=matrix(rexp(n*resample),nrow=n)
  #
  # parre=apply(I,2,resam1, initial, X, t0, d0)
  #
  # convergeno[ll]=sum(parre[1,]==0,na.rm=T)
  # parse[,ll]=apply(parre[-1,parre[1,]==0],1,sd,na.rm=T)

  ######## Meier
  theta <- rep(.8, 5) # c(0.9, 0.88, 0.85, 0.8, 0.75)
  phi <- rep(1, 5)
  truth.m <- c(beta, lambda0)
  initial.m <- c(truth.m[1], qlogis(truth.m[-1]))
  ffun <- function(initial.m) {
    beta <- initial.m[1]
    gamma0 <- initial.m[-1]
    eta <- c(X %*% beta)
    M <- t(outer(1 + exp(gamma0), -exp(eta), "^"))
    lik <- rep(NA, n)
    for (i in 1:n) {
      phi0 <- phi[t0[i]]
      Gamma <- prod(phi[1:t0[i]]) * ((1 - phi0) / phi0)^d0[i]
      theta0 <- theta[t0[i]]
      Delta <- prod(1 - theta[1:t0[i]]) * (theta0 / (1 - theta0))^d0[i] # k = 1
      lik[i] <- prod(M[i, 1:t0[i]]) * Gamma
      lik[i] <- lik[i] + (1 - M[i, 1]) * Delta
      if (t0[i] > 1) {
        for (k in 2:t0[i]) {
          theta0 <- theta[t0[i] - k + 1]
          Delta <- prod(phi[1:(k - 1)]) * prod(1 - theta[0:(t0[i] - k) + 1]) * (theta0 / (1 - theta0))^d0[i]
          lik[i] <- lik[i] + prod(M[i, 1:(k - 1)]) * (1 - M[i, k]) * Delta
        }
      }
    }
    (-sum(log(lik)))
  }
  temp <- optim(initial.m, ffun, method = "BFGS")
  code.m[ll] <- temp$convergence
  parhat.m[, ll] <- temp$par

  print(ll)
}

######## cluster
write.table(cbind(c(code)), paste(getwd(), "/", "code.csv", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
write.table(cbind(c(1111, parhat)), paste(getwd(), "/", "parhat.csv", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
# write.table(cbind(c(1111,auchat)), paste(getwd(),"/", "auchat.csv", sep =""),append=TRUE,col.names = FALSE, row.names= FALSE)
write.table(cbind(c(1111, tpr.hat)), paste(getwd(), "/", "tprhat.csv", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
write.table(cbind(c(1111, fpr.hat)), paste(getwd(), "/", "fprhat.csv", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
write.table(cbind(c(convergeno)), paste(getwd(), "/", "convergeno.csv", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
write.table(cbind(c(1111, parse)), paste(getwd(), "/", "parse.csv", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
write.table(cbind(c(code.m)), paste(getwd(), "/", "codem.csv", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
write.table(cbind(c(1111, parhat.m)), paste(getwd(), "/", "parhatm.csv", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)

######## results
rm(list = ls())
library(readr)
code <- read_csv("/Users/xuanwang/Desktop/O2_Cluster/code.csv")
code <- code[[1]]
parhat <- read_csv("/Users/xuanwang/Desktop/O2_Cluster/parhat.csv")
par <- parhat[[1]]
tprhat <- read_csv("/Users/xuanwang/Desktop/O2_Cluster/tprhat.csv")
tpr <- tprhat[[1]]
fprhat <- read_csv("/Users/xuanwang/Desktop/O2_Cluster/fprhat.csv")
fpr <- fprhat[[1]]
convergeno <- read_csv("/Users/xuanwang/Desktop/O2_Cluster/convergeno.csv")
convergeno <- convergeno[[1]]
parse <- read_csv("/Users/xuanwang/Desktop/O2_Cluster/parse.csv")
parse <- parse[[1]]
codem <- read_csv("/Users/xuanwang/Desktop/O2_Cluster/codem.csv")
code.m <- codem[[1]]
parhatm <- read_csv("/Users/xuanwang/Desktop/O2_Cluster/parhatm.csv")
par.m <- parhatm[[1]]
each <- 5

index <- which(par == 1111)
parhat <- NULL
for (j in 1:length(index)) {
  parhat <- cbind(parhat, matrix(par[(index[j] + 1):(index[j] + 11 * each)], ncol = each))
}
index <- which(tpr == 1111)
tprhat <- NULL
fprhat <- NULL
for (j in 1:length(index)) {
  tprhat <- cbind(tprhat, matrix(tpr[(index[j] + 1):(index[j] + 5 * each)], ncol = each))
  fprhat <- cbind(fprhat, matrix(fpr[(index[j] + 1):(index[j] + 5 * each)], ncol = each))
}
index <- which(parse == 1111)
parsehat <- NULL
for (j in 1:length(index)) {
  parsehat <- cbind(parsehat, matrix(parse[(index[j] + 1):(index[j] + 21 * each)], ncol = each))
}
index <- which(par.m == 1111)
parhat.m <- NULL
for (j in 1:length(index)) {
  parhat.m <- cbind(parhat.m, matrix(par.m[(index[j] + 1):(index[j] + 6 * each)], ncol = each))
}

#### Meier
print(initial.m)
print(code.m)
print(apply(parhat.m[, code.m == 0], 1, mean))
parhat.m[-1, ] <- plogis(parhat.m[-1, ])
true <- truth.m
result <- data.frame(truth = (true), mean = rowMeans(parhat.m), bias = rowMeans(parhat.m) - true, sd = apply(parhat.m, 1, sd))
print(round(result, 3))

#### proposed
print(initial)
print(code)
# index1=which(parhat[1,]>0)
# index2=which(code==0)
# index=intersect(index1,index2)
# parhat=parhat[,index]
print(apply(parhat, 1, mean))
# parsehat=parsehat[,index]
# auchat=auchat[,index]

parhat[-1, ] <- plogis(parhat[-1, ])
true <- truth
result <- data.frame(
  truth = (true), est = rowMeans(parhat), bias = rowMeans(parhat) - true, sd = apply(parhat, 1, sd),
  se = apply(parsehat[1:11, ], 1, mean),
  cp = apply((parhat + 1.96 * parsehat[1:11, ] > true) * (parhat - 1.96 * parsehat[1:11, ] < true), 1, mean)
)
rownames(result) <- c(
  "beta", "lambda0_1", "lambda0_2", "lambda0_3", "lambda0_4", "lambda0_5",
  "lambda0star_1", "lambda0star_2", "lambda0star_3", "lambda0star_4", "lambda0star_5"
)
print(round(result, 3))

## set1
tprtrue <- c(0.843, 0.764, 0.666, 0.587, 0.539)
fprtrue <- c(0.336, 0.109, 0.012, 0.000, 0.000)
# auctrue=c(0.664, 0.892, 0.988, 1.000, 1.000)
result.tpr <- data.frame(
  truth = (tprtrue), est = rowMeans(tprhat), bias = rowMeans(tprhat) - tprtrue, sd = apply(tprhat, 1, sd),
  se = apply(parsehat[12:16, ], 1, mean),
  cp = apply((tprhat + 1.96 * parsehat[12:16, ] > tprtrue) * (tprhat - 1.96 * parsehat[12:16, ] < tprtrue), 1, mean)
)
rownames(result.tpr) <- c("tpr_1", "tpr_2", "tpr_3", "tpr_4", "tpr_5")
print(round(result.tpr, 3))
result.fpr <- data.frame(
  truth = (fprtrue), est = rowMeans(fprhat), bias = rowMeans(fprhat) - fprtrue, sd = apply(fprhat, 1, sd),
  se = apply(parsehat[17:21, ], 1, mean),
  cp = apply((fprhat + 1.96 * parsehat[17:21, ] > fprtrue) * (fprhat - 1.96 * parsehat[17:21, ] < fprtrue), 1, mean)
)
rownames(result.fpr) <- c("fpr_1", "fpr_2", "fpr_3", "fpr_4", "fpr_5")
print(round(result.fpr, 3))

# # set2
# auctrue=c(0.790, 0.830, 0.869, 0.902, 0.929)
# result.auc=data.frame(truth = (auctrue), est = rowMeans(auchat),bias = rowMeans(auchat) - auctrue,sd = apply(auchat, 1, sd)#,
#                       # se=apply(parsehat[12:16,],1,mean),
#                       # cp=apply((auchat+1.96*parsehat[12:16,]>auctrue)*(auchat-1.96*parsehat[12:16,]<auctrue),1,mean)
#                       )
# rownames(result.auc)=c("auc1","auc2","auc3","auc4","auc5")
# print(round(result.auc,3))


all <- rbind(result, result.tpr, result.fpr)
library(xtable)
print(xtable(all, digits = 3), sanitize.text.function = function(x) {
  x
})


timeend <- Sys.time()
runningtime <- timeend - timestart
print(runningtime)
