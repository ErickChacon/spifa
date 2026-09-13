# Deviance Information Criterion for a spifa Model

Computes the Deviance Information Criterion (DIC) for a fitted `spifa`
model, useful for comparing candidate models (e.g. different numbers of
factors or different restriction structures) fitted to the same data.

## Usage

``` r
# S3 method for class 'spifa'
dic(x, burnin = 0, thin = 1, ...)
```

## Arguments

- x:

  A fitted `spifa` object, as returned by
  [`spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md).

- burnin:

  Number of initial iterations to discard.

- thin:

  Thinning interval applied after discarding burn-in.

- ...:

  Further arguments passed to methods (currently unused).

## Value

A one-row [`tibble`](https://tibble.tidyverse.org/reference/tibble.html)
with columns `mean_deviance` (posterior mean of the deviance), `p_eff`
(effective number of parameters), and `dic` – naming follows INLA's
`$dic` output (`mean.deviance`/`p.eff`/ `dic`).

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
dic(samples)
#> # A tibble: 1 × 3
#>   mean_deviance p_eff   dic
#>           <dbl> <dbl> <dbl>
#> 1          993.  135. 1128.
# }
```
