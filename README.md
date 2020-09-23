
# Evaluating Biomarkers using Panel Current Status Data with Misclassication Errors

**NOTE**: WORK IN PROGRESS\!

## Installation

Install development version from GitHub.

``` r
# install.packages("remotes")
remotes::install_github("celehs/misclassification")
```

Load the package into R.

``` r
library(misclassification)
```

## Simulated Example

``` r
source("examples/mkdata.R")
```

``` r
data <- mkdata(
      beta = 1,
      lambda0 = 0.01 * (1:5),
      theta = c(0.9, 0.88, 0.85, 0.8, 0.75),
      phi = rep(0.95, 5),
      n = 5000,
      seed = 1)
data.table::data.table(data)
```

    ##       t0 d0        z
    ##    1:  2  1 2.937355
    ##    2:  4  1 3.018364
    ##    3:  5  0 2.916437
    ##    4:  2  1 3.159528
    ##    5:  2  1 3.032951
    ##   ---               
    ## 4996:  2  1 3.016255
    ## 4997:  4  1 3.098074
    ## 4998:  4  1 2.930786
    ## 4999:  5  1 2.999651
    ## 5000:  5  0 3.017158

``` r
with(data, addmargins(table(d0, t0)))
```

    ##      t0
    ## d0       1    2    3    4    5  Sum
    ##   0      0    0    0    0  474  474
    ##   1    883 1296 1165  814  368 4526
    ##   Sum  883 1296 1165  814  842 5000

### Standard Discrete Proportional Hazards (DPH) Model

``` r
summary(with(data, dph_fit(t0, d0, z))) 
```

    ## 
    ## Call:
    ## stats::glm(formula = y ~ z + Z - 1, family = stats::binomial("cloglog"))
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.3383  -0.8929  -0.6319   1.2557   1.9741  
    ## 
    ## Coefficients:
    ##    Estimate Std. Error z value Pr(>|z|)    
    ## z1  -3.9075     0.4471  -8.740  < 2e-16 ***
    ## z2  -3.2401     0.4464  -7.258 3.92e-13 ***
    ## z3  -2.8943     0.4458  -6.492 8.45e-11 ***
    ## z4  -2.6522     0.4456  -5.951 2.66e-09 ***
    ## z5  -2.8114     0.4466  -6.295 3.07e-10 ***
    ## Z    0.7556     0.1483   5.096 3.46e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 23972  on 14436  degrees of freedom
    ## Residual deviance: 17039  on 14430  degrees of freedom
    ## AIC: 17051
    ## 
    ## Number of Fisher Scoring iterations: 5

### TODO: Adjusted Proportional Hazards (APH) Model

## GitHub Repository

<https://github.com/celehs/misclassification/>
