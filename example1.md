
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Example 1 (Binary Biomarker)

``` r
source("R/setting1.R")
```

### Simulation Setting

``` r
n <- 1000
beta <- 1
lambda0 <- 0.1 * (1:5) # rep(.01, 7)
lambda0_s <- 0.1 * (8:4)
```

### Data Simulation

``` r
data1 <- mkdata1(
  n = n, 
  beta = beta, 
  lambda0 = lambda0, 
  lambda0_s = lambda0_s)
t0 <- data1$t0
d0 <- data1$d0
z <- data1$z
X <- matrix(z)
```

### Simulated Data

``` r
data1
#> # A tibble: 1,000 x 3
#>       t0    d0     z
#>    <dbl> <dbl> <dbl>
#>  1     1     1     2
#>  2     5     0     0
#>  3     1     1     2
#>  4     2     1     2
#>  5     5     0     0
#>  6     2     1     0
#>  7     3     1     0
#>  8     1     1     2
#>  9     3     1     0
#> 10     5     1     0
#> # … with 990 more rows
```

### Model Fitting

``` r
ini <- c(beta, lambda0, lambda0_s)
model_fit(t0, d0, X, initial = c(ini[1], qlogis(ini[-1])))
#> $par
#>  [1]  1.1294923 -2.0921434 -1.0301946 -0.5281412 -0.9608987  0.3805482
#>  [7]  0.5557513  0.5879113  0.9084300 -1.6605184  1.7215740
#> 
#> $value
#> [1] 1501.566
#> 
#> $counts
#> function gradient 
#>       46       21 
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
#>   2.132   0.105   2.225
```
