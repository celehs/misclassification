
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Example 2 (Continuous Biomarker)

``` r
source("setting2.R")
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
#>  1     3     1 1.39 
#>  2     1     1 2.85 
#>  3     5     0 1.96 
#>  4     1     1 2.70 
#>  5     2     1 1.48 
#>  6     1     1 1.39 
#>  7     1     1 3.01 
#>  8     2     1 1.86 
#>  9     1     1 2.85 
#> 10     4     1 0.670
#> # … with 990 more rows
```

### Model Fitting

``` r
ini <- c(beta, lambda0, lambda0_s)
model_fit(t0, d0, X, initial = c(ini[1], qlogis(ini[-1])))
#> $par
#>  [1]  1.34227956 -3.29027933 -2.41588031 -2.03458207 -1.34731755 -1.23397348
#>  [7]  0.75866434  0.48520525  1.16137808 -0.03256791 -1.80012040
#> 
#> $value
#> [1] 1362.911
#> 
#> $counts
#> function gradient 
#>       51       19 
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
#>   2.472   0.133   2.581
```
