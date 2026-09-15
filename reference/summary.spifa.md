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
#>    variable     mean   median    sd   mad    q2.5     q10     q90     q97.5
#>    <chr>       <dbl>    <dbl> <dbl> <dbl>   <dbl>   <dbl>   <dbl>     <dbl>
#>  1 c[1]     -0.438   -0.472   0.198 0.228 -0.728  -0.643  -0.198  -0.100   
#>  2 c[2]     -0.531   -0.512   0.214 0.180 -0.930  -0.849  -0.299  -0.192   
#>  3 c[3]     -0.563   -0.552   0.330 0.361 -1.05   -0.994  -0.0560  0.000170
#>  4 c[4]     -0.466   -0.487   0.192 0.278 -0.741  -0.720  -0.235  -0.154   
#>  5 c[5]     -0.266   -0.257   0.113 0.126 -0.447  -0.439  -0.123  -0.0947  
#>  6 c[6]     -0.0101  -0.00559 0.116 0.111 -0.226  -0.143   0.122   0.169   
#>  7 c[7]      0.103    0.114   0.149 0.171 -0.149  -0.0822  0.304   0.331   
#>  8 c[8]     -0.00922 -0.0283  0.103 0.113 -0.146  -0.118   0.146   0.155   
#>  9 c[9]      0.503    0.491   0.168 0.191  0.284   0.317   0.725   0.838   
#> 10 c[10]     0.316    0.364   0.186 0.169 -0.0717  0.112   0.494   0.571   
#> # ℹ 3 more variables: ess_bulk <dbl>, ess_tail <dbl>, rhat <dbl>
# }
```
