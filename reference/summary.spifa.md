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

  Character vector of parameter groups to summarise (defaults to all of
  them). An error if any requested group does not exist in the fitted
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
#>    variable    mean  median    sd    mad    q2.5     q10    q90   q97.5 ess_bulk
#>    <chr>      <dbl>   <dbl> <dbl>  <dbl>   <dbl>   <dbl>  <dbl>   <dbl>    <dbl>
#>  1 c[1]     -0.375  -0.403  0.244 0.221  -0.710  -0.639  -0.119  0.145      2.89
#>  2 c[2]     -0.681  -0.765  0.220 0.305  -0.980  -0.965  -0.426 -0.324      2.34
#>  3 c[3]     -0.581  -0.615  0.243 0.190  -0.857  -0.810  -0.296 -0.0222     8.73
#>  4 c[4]     -0.307  -0.316  0.161 0.109  -0.611  -0.423  -0.133  0.0176    13.2 
#>  5 c[5]     -0.229  -0.181  0.124 0.110  -0.470  -0.369  -0.121 -0.0424    17.7 
#>  6 c[6]      0.0696  0.0676 0.122 0.0830 -0.146  -0.0335  0.183  0.307      3.37
#>  7 c[7]      0.118   0.138  0.146 0.127  -0.151  -0.0922  0.300  0.335     18.8 
#>  8 c[8]      0.0848  0.0768 0.175 0.217  -0.211  -0.173   0.284  0.326      7.47
#>  9 c[9]      0.387   0.364  0.203 0.105   0.0569  0.173   0.631  0.811      6.55
#> 10 c[10]     0.142   0.174  0.165 0.194  -0.123  -0.0917  0.327  0.366      8.12
#> # ℹ 2 more variables: ess_tail <dbl>, rhat <dbl>
# }
```
