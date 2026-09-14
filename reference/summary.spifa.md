# Summarise spifa Posterior Samples

Computes posterior summary statistics for every parameter in a fitted
`spifa` object, via
[`summarise_draws`](https://mc-stan.org/posterior/reference/draws_summary.html):
mean, median, sd, mad, the 2.5%/10%/50%/90%/97.5% quantiles, effective
sample size (bulk and tail), and `rhat`. `rhat` is computed via
split-chain R-hat (Vehtari et al. 2021), which splits each chain in half
and compares the halves – so it remains a meaningful convergence
diagnostic even though
[`spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md) only
ever fits a single chain. Discrimination parameters (`A`) structurally
restricted to zero (via `constraints$discrimination`) are excluded,
since they are fixed by construction rather than estimated.

## Usage

``` r
# S3 method for class 'spifa'
summary(object, burnin = 0, thin = 1, select = NULL, ...)
```

## Arguments

- object:

  A fitted `spifa` object, as returned by
  [`spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md).

- burnin:

  Number of initial iterations to discard.

- thin:

  Thinning interval applied after discarding burn-in.

- select:

  Character vector of parameter blocks to summarise (defaults to all of
  them). An error if any requested block does not exist in the fitted
  model (e.g. `"T"` for a model with no spatial process).

- ...:

  Further arguments passed to methods (currently unused).

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
one row per parameter; see
[`summarise_draws`](https://mc-stan.org/posterior/reference/draws_summary.html)
for column details.

## Author

Erick A. Chacón-Montalván

## Examples

``` r
# \donttest{
data(ipixuna)
nitems <- ncol(ipixuna$items)
nfactors <- 3

# discrimination constraint: start with every item free to load on every
# factor, then restrict a few items per factor based on what each item is
# meant to measure (0 = no relationship, 1 = free parameter to estimate)
A <- matrix(1, nitems, nfactors)
A[c(4, 8), 1] <- 0
A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
A[c(5, 6), 3] <- 0
samples <- spifa(
  items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
  niter = 20, standardize = FALSE,
  constraints = list(discrimination = A))
summary(samples, select = "c")
#> # A tibble: 10 × 12
#>    variable    mean  median    sd    mad    q2.5     q10     q90   q97.5
#>    <chr>      <dbl>   <dbl> <dbl>  <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#>  1 c[1]     -0.398  -0.392  0.120 0.0902 -0.628  -0.560  -0.261  -0.217 
#>  2 c[2]     -0.467  -0.445  0.169 0.189  -0.744  -0.639  -0.270  -0.154 
#>  3 c[3]     -0.547  -0.642  0.292 0.212  -0.915  -0.810  -0.220   0.0702
#>  4 c[4]     -0.324  -0.309  0.124 0.131  -0.552  -0.508  -0.170  -0.158 
#>  5 c[5]     -0.0754 -0.0662 0.119 0.154  -0.245  -0.209   0.0480  0.140 
#>  6 c[6]      0.130   0.151  0.102 0.0734 -0.0782 -0.0313  0.225   0.260 
#>  7 c[7]      0.251   0.294  0.156 0.161  -0.0416  0.0655  0.440   0.465 
#>  8 c[8]      0.193   0.174  0.199 0.290  -0.0746 -0.0491  0.430   0.500 
#>  9 c[9]      0.551   0.547  0.137 0.142   0.339   0.373   0.719   0.810 
#> 10 c[10]     0.221   0.212  0.118 0.0849  0.0265  0.0920  0.396   0.441 
#> # ℹ 3 more variables: ess_bulk <dbl>, ess_tail <dbl>, rhat <dbl>
# }
```
