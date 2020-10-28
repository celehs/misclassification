
# Evaluating Biomarkers using Panel Current Status Data with Misclassication Errors

**NOTE**: WORK IN PROGRESS (WIP)\!

**GitHub Repository**: <https://github.com/celehs/misclassification/>

![](flowchart/flowchart-misclassification.jpg)

## Installation

Install development version from GitHub.

``` r
# install.packages("devtools")
devtools::install_github("celehs/misclassification")
```

Load the package into R.

``` r
library(misclassification)
```

## NEW Examples (WIP)

  - [Example with Binary Biomarker](examples/setting1.md)

  - [Example with Continuous Biomarker](examples/setting2.md)

## OLD Example

``` r
source("examples/mkdata.R")
```

``` r
beta <- 1.5
lambda0 <- 0.01 * (1:5)
theta <- c(0.9, 0.88, 0.85, 0.8, 0.75)
phi <- rep(0.95, 5)
n <- 5000
data <- mkdata(
      beta = beta,
      lambda0 = lambda0,
      theta = theta,
      phi = phi,
      n = n,
      seed = 1)
data.table::data.table(data)
```

    ##       t0 d0        z
    ##    1:  1  1 2.937355
    ##    2:  2  1 3.018364
    ##    3:  2  1 2.916437
    ##    4:  1  1 3.159528
    ##    5:  1  1 3.032951
    ##   ---               
    ## 4996:  1  1 3.016255
    ## 4997:  2  1 3.098074
    ## 4998:  3  1 2.930786
    ## 4999:  2  1 2.999651
    ## 5000:  2  1 3.017158

``` r
with(data, addmargins(table(d0, t0)))
```

    ##      t0
    ## d0       1    2    3    4    5  Sum
    ##   0      0    0    0    0   25   25
    ##   1   2734 1727  442   67    5 4975
    ##   Sum 2734 1727  442   67   30 5000

### Model Fitting

``` r
# sensitivity unknown, specificity fixed
start <- c(beta, qlogis(lambda0), qlogis(theta)) # 1 + 5 + 5
fit <- optim(
   par = start,
   fn = function(par) {
      -loglik(
         beta = par[1],
         lambda0 = plogis(par[1 + 1:5]),
         theta = plogis(par[6 + 1:5]), # sensitivity
         phi = phi, # specificity
         t0 = data$t0,
         d0 = data$d0,
         X = data$z)
      },
      method = "BFGS")
# convergence (An integer code. 0 indicates successful completion)
fit$convergence
```

    ## [1] 0

``` r
data.frame(name = c("beta", 
                    paste0("lambda0_", 1:5), 
                    paste0("theta_", 1:5)),
           truth = c(start[1], plogis(start[-1])), 
           param = round(c(fit$par[1], plogis(fit$par[-1])), 3))
```

    ##         name truth param
    ## 1       beta  1.50 1.334
    ## 2  lambda0_1  0.01 0.014
    ## 3  lambda0_2  0.02 0.026
    ## 4  lambda0_3  0.03 0.035
    ## 5  lambda0_4  0.04 0.028
    ## 6  lambda0_5  0.05 0.000
    ## 7    theta_1  0.90 0.995
    ## 8    theta_2  0.88 0.204
    ## 9    theta_3  0.85 0.036
    ## 10   theta_4  0.80 0.639
    ## 11   theta_5  0.75 0.006

``` r
# cumulative baseline hazards
data.frame(name = paste0("Lambda0_", 1:5),
           truth = cumsum(lambda0),
           param = round(cumsum(plogis(fit$par[1 + 1:5])), 3))
```

    ##        name truth param
    ## 1 Lambda0_1  0.01 0.014
    ## 2 Lambda0_2  0.03 0.040
    ## 3 Lambda0_3  0.06 0.074
    ## 4 Lambda0_4  0.10 0.102
    ## 5 Lambda0_5  0.15 0.102

### Accuracy Measures

``` r
# AUC at estimated parameters
est <- acc_est(
   Z = data$z, 
   beta = fit$par[1], 
   lambda0 = plogis(fit$par[1 + 1:5]))
est$AUC
```

    ## [1] 0.5545099 0.5922320 0.6498856 0.6955478 0.6955681

``` r
# AUC at true parameters
acc_est(
   Z = data$z, 
   beta = beta, 
   lambda0 = lambda0)$AUC
```

    ## [1] 0.5654422 0.6203994 0.7086030 0.8068632 0.8903748

``` r
par(mfrow = c(3, 2))
for (m in 1:5) {
   x <- c(0, est$FPR[, m])
   y <- c(0, est$TPR[, m])
   plot(x, y, las = 1, 
        type = "s", col = "blue",
        xaxs = "i", yaxs = "i", 
        xlab = "FPR", ylab = "TPR",
        main = paste0("t = ", m, " (AUC = ", 
                     sprintf("%.3f", est$AUC[m]), ")"))
      abline(a = 0, b = 1, lty = 2)   
}
```

![](README_files/figure-gfm/ROC-1.png)<!-- -->

``` r
proc.time()
```

    ##    user  system elapsed 
    ##  31.655   2.456  34.109
