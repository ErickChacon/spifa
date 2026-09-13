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
ever fits a single chain.

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
#>    variable    mean  median    sd    mad    q2.5     q10      q90   q97.5
#>    <chr>      <dbl>   <dbl> <dbl>  <dbl>   <dbl>   <dbl>    <dbl>   <dbl>
#>  1 c[1]     -0.238  -0.323  0.342 0.234  -0.592  -0.517   0.0310   0.594 
#>  2 c[2]     -0.381  -0.372  0.131 0.0883 -0.638  -0.572  -0.222   -0.196 
#>  3 c[3]     -0.377  -0.407  0.291 0.271  -0.745  -0.701   0.0354   0.220 
#>  4 c[4]     -0.164  -0.155  0.138 0.133  -0.352  -0.324  -0.00761  0.107 
#>  5 c[5]     -0.239  -0.263  0.136 0.137  -0.507  -0.370  -0.0905  -0.0393
#>  6 c[6]      0.0349  0.0515 0.101 0.120  -0.139  -0.0959  0.145    0.176 
#>  7 c[7]      0.285   0.294  0.109 0.0874  0.111   0.142   0.440    0.471 
#>  8 c[8]      0.104   0.113  0.137 0.149  -0.0984 -0.0616  0.279    0.355 
#>  9 c[9]      0.585   0.625  0.205 0.143   0.104   0.356   0.755    0.866 
#> 10 c[10]     0.384   0.395  0.186 0.206   0.0865  0.101   0.590    0.647 
#> # ℹ 3 more variables: ess_bulk <dbl>, ess_tail <dbl>, rhat <dbl>
# }
```
