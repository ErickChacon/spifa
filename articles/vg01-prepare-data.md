# Preparing Your Data

## Introduction

[`spifa()`](https://ErickChacon.github.io/spifa/reference/spifa.md)
expects data in a specific shape: one row per respondent, and the binary
item responses stored as a single *matrix-valued column*
(`nobs x nitems`) rather than one column per item. This is the same
convention base R uses for multivariate
[`lm()`](https://rdrr.io/r/stats/lm.html). Survey data rarely arrives in
that shape – this vignette builds a small questionnaire dataset from
scratch, the way it typically comes out of a survey tool (one column per
item, `"Yes"`/ `"No"` responses, a household covariate, GPS
coordinates), and walks through getting it into the shape
[`spifa()`](https://ErickChacon.github.io/spifa/reference/spifa.md)
needs.

``` r

library(dplyr)
```

    #> 
    #> Attaching package: 'dplyr'

    #> The following objects are masked from 'package:stats':
    #> 
    #>     filter, lag

    #> The following objects are masked from 'package:base':
    #> 
    #>     intersect, setdiff, setequal, union

``` r

library(sf)
```

    #> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE

``` r

library(spifa)
```

## Simulating a raw questionnaire export

Five binary items, a household `wealth` covariate, and GPS coordinates,
for 20 households – exactly what a survey export typically looks like:

``` r

set.seed(1)
nobs <- 20
nitems <- 5

raw <- data.frame(
  id = seq_len(nobs),
  wealth = rnorm(nobs),
  lon = runif(nobs, -71.70, -71.68),
  lat = runif(nobs, -7.06, -7.04))

for (j in seq_len(nitems)) {
  raw[[paste0("item", j)]] <- factor(ifelse(rbinom(nobs, 1, 0.5) == 1, "Yes", "No"))
}

head(raw, 3)
```

    #>   id     wealth       lon       lat item1 item2 item3 item4 item5
    #> 1  1 -0.6264538 -71.68358 -7.041742    No   Yes   Yes   Yes    No
    #> 2  2  0.1836433 -71.68706 -7.054128   Yes    No    No   Yes   Yes
    #> 3  3 -0.8356286 -71.68434 -7.050819    No    No    No    No    No

## From one column per item to a matrix column

Recode each item to binary `0`/`1`
([`spifa()`](https://ErickChacon.github.io/spifa/reference/spifa.md)
also accepts `TRUE`/`FALSE`, coerced the same way), then assemble the
item columns into a single matrix column with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) + ordinary `$<-`
assignment (or [`I()`](https://rdrr.io/r/base/AsIs.html)), which
preserves it as a matrix rather than flattening it back into separate
columns:

``` r

data <- raw |>
  mutate(across(starts_with("item"), ~ as.integer(.x == "Yes")))

data$items <- as.matrix(select(data, starts_with("item")))
data <- select(data, id, wealth, lon, lat, items)

str(data)
```

    #> 'data.frame':    20 obs. of  5 variables:
    #>  $ id    : int  1 2 3 4 5 6 7 8 9 10 ...
    #>  $ wealth: num  -0.626 0.184 -0.836 1.595 0.33 ...
    #>  $ lon   : num  -71.7 -71.7 -71.7 -71.7 -71.7 ...
    #>  $ lat   : num  -7.04 -7.05 -7.05 -7.05 -7.05 ...
    #>  $ items : int [1:20, 1:5] 0 1 0 0 1 0 1 0 0 0 ...
    #>   ..- attr(*, "dimnames")=List of 2
    #>   .. ..$ : NULL
    #>   .. ..$ : chr [1:5] "item1" "item2" "item3" "item4" ...

Passing the original one-column-per-item `raw` data (before this step)
to [`spifa()`](https://ErickChacon.github.io/spifa/reference/spifa.md)
fails clearly rather than silently doing the wrong thing:

``` r

spifa(item1 ~ 1, data = raw, nfactors = 2, niter = 5)
```

    #> Error in `spifa()`:
    #> ! The left-hand side of 'formula' must be a matrix

## Missing item responses

Individual missing item responses (`NA`) don’t need to be dropped or
imputed –
[`spifa()`](https://ErickChacon.github.io/spifa/reference/spifa.md)
handles them natively, treating a missing response as an unobserved
auxiliary variable to be sampled along with everything else, rather than
requiring complete rows:

``` r

data_na <- data
data_na$items[1, 1] <- NA

fit <- spifa(items ~ 1, data = data_na, nfactors = 2, ngp = 0, niter = 5)
attr(fit, "fit_args")$nobs  # respondent with the NA item is still included
```

    #> [1] 20

## Missing predictor values

Missing *predictor* values behave differently, and are worth calling
out: a `NA` in a predictor column drops that respondent from the whole
fit, the same way [`lm()`](https://rdrr.io/r/stats/lm.html) and friends
do – unlike item responses, there’s no partial information to keep:

``` r

data_na <- data
data_na$wealth[1] <- NA

fit <- spifa(items ~ wealth, data = data_na, nfactors = 2, ngp = 0, niter = 5)
attr(fit, "fit_args")$nobs  # one respondent dropped
```

    #> [1] 19

## Adding spatial coordinates

Spatial models need `data` to be an `sf` object. If your coordinates are
plain numeric columns (e.g. from a GPS device), as here, convert with
[`sf::st_as_sf()`](https://r-spatial.github.io/sf/reference/st_as_sf.html):

``` r

data_sf <- st_as_sf(data, coords = c("lon", "lat"), crs = 4326)
```

[`spifa()`](https://ErickChacon.github.io/spifa/reference/spifa.md) then
adds the spatial Gaussian process automatically (see
[`vignette("spifa")`](https://ErickChacon.github.io/spifa/articles/spifa.md)):

``` r

attr(spifa(items ~ 1, data = data_sf, nfactors = 2, niter = 5),
  "fit_args")$model_type
```

    #> [1] "spifa"

## Putting it together

``` r

fit <- spifa(items ~ wealth, data = data_sf, nfactors = 2, niter = 5)
fit
```

    #> Item factor analysis model: spifa_pred
    #> Formula: items ~ wealth
    #> Dimensions: 20 respondents, 5 items, 2 latent factors, 2 spatial processes
    #> MCMC: 1 chain, iter = 5, thin = 1, samples = 5
    #> 
    #> Item model parameters:
    #>          mean median    sd     q10     q90 ess_bulk rhat
    #> c[1]   -0.175 -0.054 0.373 -0.5564  0.1199       NA  1.9
    #> c[2]   -0.085 -0.029 0.237 -0.3235  0.1017       NA  1.9
    #> c[3]    0.157  0.106 0.160  0.0093  0.3281       NA  1.9
    #> c[4]   -0.178 -0.137 0.071 -0.2560 -0.1201       NA  0.9
    #> c[5]    0.138  0.265 0.268 -0.1554  0.3583       NA  0.9
    #> A[1,1]  0.978  1.029 0.218  0.7471  1.1772       NA  0.9
    #> A[2,1]  0.361  0.149 0.441  0.0269  0.8478       NA  1.9
    #> A[3,1]  0.698  0.605 0.613  0.1768  1.2795       NA  1.9
    #> A[4,1]  0.185  0.157 0.273 -0.0901  0.4621       NA  1.9
    #> A[5,1]  0.275  0.234 0.150  0.1710  0.4234       NA  1.9
    #> A[2,2]  0.140  0.085 0.437 -0.2403  0.5987       NA  1.9
    #> A[3,2] -0.645 -0.891 0.680 -1.2609  0.0863       NA  1.9
    #> A[4,2] -0.294 -0.362 0.300 -0.5338  0.0069       NA  0.9
    #> A[5,2] -0.046 -0.019 0.094 -0.1464  0.0401       NA  1.9
    #> 
    #> Factor model parameters:
    #>              mean  median     sd    q10   q90 ess_bulk rhat
    #> B[1,1]     0.1159  0.1732 0.2755 -0.150 0.312       NA  0.9
    #> B[1,2]    -0.2285 -0.3519 0.2725 -0.429 0.068       NA  0.9
    #> T[1,1]     0.4830  0.4819 0.0200  0.463 0.502       NA  1.9
    #> T[2,2]     0.5197  0.5220 0.0129  0.506 0.530       NA  1.9
    #> phi[1]     0.0495  0.0499 0.0024  0.047 0.052       NA  1.9
    #> phi[2]     0.0498  0.0501 0.0015  0.048 0.051       NA  0.9
    #> Corr[2,1]  0.0049 -0.0045 0.0599 -0.052 0.069       NA  1.9
    #> 
    #> ess_bulk is the bulk effective sample size; rhat is the potential
    #> scale reduction factor on split chains (Rhat = 1 at convergence).
    #> Use summary() for the full set of statistics (incl. ess_tail).
