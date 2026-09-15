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
to `burnin`/`thin`) instead of calling the `C++` sampler. Otherwise,
prediction for the new locations and/or predictor values is delegated to
the `C++` sampler.

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
#>         1 -2.39
#>         2 -1.72
#>         3 -0.63
#>         4 -0.40
#>         5 -0.93
#> 
#> , , variable = Theta[2,1]
#> 
#>          chain
#> iteration     1
#>         1 -1.02
#>         2  0.36
#>         3 -0.72
#>         4  0.55
#>         5  0.59
#> 
#> , , variable = Theta[3,1]
#> 
#>          chain
#> iteration     1
#>         1  0.15
#>         2  0.92
#>         3 -0.17
#>         4  1.02
#>         5  0.54
#> 
#> , , variable = Theta[4,1]
#> 
#>          chain
#> iteration     1
#>         1 -1.28
#>         2 -0.90
#>         3 -0.16
#>         4 -1.73
#>         5  0.43
#> 
#> # ... with 296 more variables
# latent abilities for new locations
newdata <- st_make_grid(ipixuna, n = c(3, 2), what = "centers")
predict(samples, newdata = newdata)
#> # A draws_array: 5 iterations, 1 chains, and 18 variables
#> , , variable = Theta[1,1]
#> 
#>          chain
#> iteration     1
#>         1 -0.69
#>         2  0.48
#>         3 -1.10
#>         4  0.41
#>         5  0.93
#> 
#> , , variable = Theta[2,1]
#> 
#>          chain
#> iteration     1
#>         1  0.48
#>         2 -0.14
#>         3  0.42
#>         4  0.90
#>         5  0.13
#> 
#> , , variable = Theta[3,1]
#> 
#>          chain
#> iteration     1
#>         1 -0.92
#>         2 -0.53
#>         3 -0.90
#>         4  1.37
#>         5  1.29
#> 
#> , , variable = Theta[4,1]
#> 
#>          chain
#> iteration      1
#>         1 -0.332
#>         2  0.629
#>         3 -0.002
#>         4  0.764
#>         5  1.349
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
#>         1  0.098
#>         2  0.670
#>         3  1.648
#>         4  0.162
#>         5 -0.798
#> 
#> , , variable = Theta[2,1]
#> 
#>          chain
#> iteration      1
#>         1 -0.018
#>         2  0.381
#>         3 -0.809
#>         4 -0.613
#>         5  0.233
#> 
#> , , variable = Theta[3,1]
#> 
#>          chain
#> iteration     1
#>         1 -0.66
#>         2 -0.20
#>         3 -0.90
#>         4 -0.42
#>         5 -0.77
#> 
#> , , variable = Theta[4,1]
#> 
#>          chain
#> iteration     1
#>         1 -2.09
#>         2  0.46
#>         3 -0.31
#>         4 -1.75
#>         5  1.79
#> 
#> # ... with 14 more variables
# }
```
