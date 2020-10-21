
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Example 2 (Continuous Biomarker)

``` r
source("R/setting2.R")
```

### Simulation Setting

``` r
n <- 1000
beta <- 1
lambda0 <- 0.05 * (1:5) # rep(.01, 7)
lambda0_s <- 0.1 * (8:4)
```

### Data Simulation

``` r
data2 <- mkdata2(
  n = n, 
  beta = beta, 
  lambda0 = lambda0, 
  lambda0_s = lambda0_s)
t0 <- data2$t0
d0 <- data2$d0
z <- data2$z
X <- matrix(z)
```

### Simulated Data

``` r
data2
#> # A tibble: 1,000 x 3
#>       t0    d0     z
#>    <dbl> <dbl> <dbl>
#>  1     2     1 2.25 
#>  2     1     1 2.06 
#>  3     1     1 3.14 
#>  4     5     0 0.590
#>  5     5     0 0.717
#>  6     2     1 1.40 
#>  7     1     1 3.29 
#>  8     2     1 0.883
#>  9     1     1 4.24 
#> 10     1     1 1.12 
#> # … with 990 more rows
```

### Model Fitting

``` r
ini <- c(beta, lambda0, lambda0_s)
model_fit(t0, d0, X, initial = c(ini[1], qlogis(ini[-1])))
#> $par
#>  [1]  0.9746244 -2.3605770 -2.0800356 -1.5747380 -1.0840034 -1.6044074
#>  [7]  0.6837044  1.1149218  0.2719195  0.2110195 -0.4963027
#> 
#> $value
#> [1] 1446.132
#> 
#> $counts
#> function gradient 
#>       47       21 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> NULL
```

``` r
proc.time()
#>    user  system elapsed 
#>   2.430   0.125   2.541
```
