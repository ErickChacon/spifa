# Predict the Latent Factors of a spifa Model

Predicts the latent factors of (spatial) item factor analysis for new
subjects/locations and/or for new predictor values, using the posterior
samples from a fitted
[`spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md) model.

## Usage

``` r
# S3 method for class 'spifa'
predict(object, newdata = NULL, burnin = 0, thin = 1, joint = FALSE, ...)
```

## Arguments

- object:

  A fitted `spifa` object, as returned by
  [`spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md).

- newdata:

  New data to predict at, mirroring
  [`spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md)'s
  own `data` argument: an
  [`sf`](https://r-spatial.github.io/sf/reference/sf.html)/[`sfc`](https://r-spatial.github.io/sf/reference/sfc.html)
  object (for spatial, i.e. `spifa`/`spifa_pred`, models – its geometry
  gives the new locations, and its CRS is used directly, so it need not
  match the training data's CRS) or a plain data frame (for
  `cifa_pred`), with columns matching the predictors used on the
  right-hand side of `formula` when the model was fitted. Its design
  matrix is built the same way
  [`spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md)
  built the training one, using the same terms and factor levels. For a
  spatial model, its geometry may also be `POLYGON`/`MULTIPOLYGON` (e.g.
  a prediction grid): the centroid of each cell is used for the spatial
  kernel, but the original polygons are kept (see `Value`) so
  [`plot_predict`](https://ErickChacon.github.io/spifa/reference/plot_predict.md)
  can draw them as a filled map instead of points.

- burnin:

  Number of initial (post-fitting) iterations to discard before using
  the posterior samples for prediction.

- thin:

  Thinning interval applied to the posterior samples used for
  prediction.

- joint:

  Logical; for spatial models (`spifa`/`spifa_pred`), whether the
  posterior predictive draws should respect the full predictive
  covariance across new locations and factors (`TRUE`), or be drawn
  marginally/independently per location-factor combination (`FALSE`, the
  default, cheaper). Each draw still propagates posterior parameter
  uncertainty (one draw per retained MCMC iteration) either way – this
  only controls whether, within a single draw, the values across new
  locations/factors are jointly correlated as the model implies.
  Marginal draws are fine for per-location summaries (e.g. means,
  credible intervals computed independently per column); set `TRUE` when
  the samples themselves will be used as input to another model or
  computation that depends on their joint structure (e.g. a spatial
  contrast or aggregate across new locations). Ignored for `cifa_pred`,
  whose draws are already jointly correct across factors.

- ...:

  Further arguments (currently unused).

## Value

A
[`draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
of posterior predictive samples of the latent abilities (`theta`) for
the requested new locations and/or predictor values (or, if no
prediction was requested, for the originally observed subjects). For a
spatial model, the locations used (`newdata`, or the training locations
if `newdata` was omitted) are attached as the `"newdata"` attribute, so
[`plot_predict`](https://ErickChacon.github.io/spifa/reference/plot_predict.md)
can map the result without needing it supplied again.

## Details

If the fitted model has no spatial or predictor structure (`eifa` or
`cifa`), or if `newdata` is not supplied for a model that has one, there
is nothing to predict beyond the latent abilities' own posterior samples
already available from the fit, so those are returned directly (subject
to `burnin`/`thin`) – this case is fast, since no new sampling is
needed. Otherwise, new draws are sampled for the new locations and/or
predictor values, which takes longer.

If the fitted model has predictors (`cifa_pred`/`spifa_pred`), `newdata`
must include those predictor columns, the same as
[`predict.lm`](https://rdrr.io/r/stats/predict.lm.html) and similar
methods require – this is an error otherwise. There's no synthesized
reference-level fallback for missing predictors (e.g. an `sf`/`sfc`
object holding only new locations, with no predictor columns at all):
for a factor predictor under the no-intercept encoding
[`spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md) uses,
an automatically-filled all-zero row wouldn't correspond to any real
category, so any reference profile – including all zeros for numeric
predictors – must be supplied explicitly in `newdata`, matching the
original data's format.

## Author

Erick A. Chacón-Montalván

## Examples

``` r
# \donttest{
library(sf)
data(ipixuna)

nitems <- ncol(ipixuna$items)
nfactors <- 3
A <- matrix(1, nitems, nfactors)
A[c(4, 8), 1] <- 0
A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
A[c(5, 6), 3] <- 0

# Spifa model
samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, niter = 5,
  constraints = list(discrimination = A))
# latent abilities for observed locations
predict(samples)
#> # A draws_array: 5 iterations, 1 chains, and 300 variables
#> , , variable = Theta[1,1]
#> 
#>          chain
#> iteration     1
#>         1 -1.38
#>         2 -1.58
#>         3 -1.50
#>         4 -0.81
#>         5 -0.68
#> 
#> , , variable = Theta[2,1]
#> 
#>          chain
#> iteration     1
#>         1  1.64
#>         2  1.30
#>         3  0.92
#>         4  0.54
#>         5 -0.22
#> 
#> , , variable = Theta[3,1]
#> 
#>          chain
#> iteration     1
#>         1 -0.40
#>         2  1.11
#>         3 -0.13
#>         4  0.89
#>         5  0.88
#> 
#> , , variable = Theta[4,1]
#> 
#>          chain
#> iteration      1
#>         1  0.351
#>         2 -1.323
#>         3  0.018
#>         4 -0.514
#>         5  0.427
#> 
#> # ... with 296 more variables
# latent abilities for new locations
newdata <- st_make_grid(ipixuna, n = c(3, 2), what = "centers")
predict(samples, newdata = newdata)
#> # A draws_array: 5 iterations, 1 chains, and 18 variables
#> , , variable = Theta[1,1]
#> 
#>          chain
#> iteration      1
#>         1  0.572
#>         2  0.239
#>         3 -0.843
#>         4 -0.588
#>         5 -0.051
#> 
#> , , variable = Theta[2,1]
#> 
#>          chain
#> iteration      1
#>         1 -0.063
#>         2 -0.823
#>         3 -0.645
#>         4 -0.652
#>         5 -0.647
#> 
#> , , variable = Theta[3,1]
#> 
#>          chain
#> iteration     1
#>         1  0.69
#>         2  0.84
#>         3 -0.61
#>         4 -1.01
#>         5 -0.58
#> 
#> , , variable = Theta[4,1]
#> 
#>          chain
#> iteration      1
#>         1  0.636
#>         2  0.062
#>         3 -0.158
#>         4 -1.590
#>         5 -1.195
#> 
#> # ... with 14 more variables

# Spifa model with predictors
samples_pred <- spifa(items ~ wealth, data = ipixuna, nfactors = nfactors, niter = 5,
  constraints = list(discrimination = A))
newdata_pred <- st_sf(wealth = rnorm(6), geometry = newdata)
predict(samples_pred, newdata = newdata_pred)
#> # A draws_array: 5 iterations, 1 chains, and 18 variables
#> , , variable = Theta[1,1]
#> 
#>          chain
#> iteration      1
#>         1  0.387
#>         2 -0.019
#>         3  2.848
#>         4 -0.752
#>         5  0.905
#> 
#> , , variable = Theta[2,1]
#> 
#>          chain
#> iteration     1
#>         1  0.27
#>         2 -1.04
#>         3  0.38
#>         4  1.84
#>         5  0.29
#> 
#> , , variable = Theta[3,1]
#> 
#>          chain
#> iteration     1
#>         1 -1.36
#>         2 -0.48
#>         3  0.98
#>         4 -1.66
#>         5 -0.28
#> 
#> , , variable = Theta[4,1]
#> 
#>          chain
#> iteration     1
#>         1 -0.30
#>         2 -0.25
#>         3 -0.64
#>         4  2.64
#>         5  0.22
#> 
#> # ... with 14 more variables
# }
```
