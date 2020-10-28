
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Example 1 (Binary Biomarker)

``` r
source("setting1.R")
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
#>  1     3     1     2
#>  2     4     1     0
#>  3     5     0     0
#>  4     5     0     0
#>  5     3     1     2
#>  6     2     1     2
#>  7     1     1     2
#>  8     3     1     2
#>  9     2     1     2
#> 10     1     1     0
#> # … with 990 more rows
```

### Model Fitting

``` r
ini <- c(beta, lambda0, lambda0_s)
model_fit(t0, d0, X, initial = c(ini[1], qlogis(ini[-1])))
#> $par
#>  [1]  0.9153374 -2.2478150 -1.5981839 -1.0475622 -0.5560202 -0.1275710
#>  [7]  2.4446682  0.8254618 -1.0077586 -6.7846731  1.1032328
#> 
#> $value
#> [1] 1489.911
#> 
#> $counts
#> function gradient 
#>       48       24 
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
#>   2.513   0.086   2.591
```
