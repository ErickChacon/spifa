# Preparing data for spifa()

[`spifa()`](https://ErickChacon.github.io/spifa/reference/spifa.md)
expects data in a specific shape: one row per respondent, and the binary
item responses stored as a single *matrix-valued column*
(`nobs x nitems`) rather than one column per item. This is the same
convention base R uses for multivariate
[`lm()`](https://rdrr.io/r/stats/lm.html). This vignette builds a small
dataset acceptable by
[`spifa()`](https://ErickChacon.github.io/spifa/reference/spifa.md).

``` r

library(spifa)
library(sf)
```

    #> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE

## Simulating a raw questionnaire

Five binary items, a household `wealth` covariate, and GPS coordinates,
for 20 households:

``` r

set.seed(1)
nobs <- 20
nitems <- 5

raw <- data.frame(
  id = seq_len(nobs),
  wealth = rnorm(nobs),
  lon = runif(nobs, -71.70, -71.68),
  lat = runif(nobs, -7.06, -7.04),
  item1 = rbinom(nobs, 1, 0.1) == 1,
  item2 = rbinom(nobs, 1, 0.2) == 1,
  item3 = rbinom(nobs, 1, 0.3) == 1,
  item4 = rbinom(nobs, 1, 0.4) == 1,
  item5 = rbinom(nobs, 1, 0.5) == 1
)
raw
```

    #>    id      wealth       lon       lat item1 item2 item3 item4 item5
    #> 1   1 -0.62645381 -71.68358 -7.041742 FALSE FALSE  TRUE  TRUE FALSE
    #> 2   2  0.18364332 -71.68706 -7.054128 FALSE FALSE FALSE  TRUE  TRUE
    #> 3   3 -0.83562861 -71.68434 -7.050819 FALSE FALSE FALSE FALSE FALSE
    #> 4   4  1.59528080 -71.68894 -7.053352 FALSE  TRUE FALSE FALSE  TRUE
    #> 5   5  0.32950777 -71.68941 -7.046983 FALSE FALSE  TRUE  TRUE  TRUE
    #> 6   6 -0.82046838 -71.68421 -7.054840 FALSE FALSE FALSE FALSE FALSE
    #> 7   7  0.48742905 -71.69953 -7.050429 FALSE FALSE FALSE FALSE FALSE
    #> 8   8  0.73832471 -71.69046 -7.044674 FALSE FALSE FALSE  TRUE FALSE
    #> 9   9  0.57578135 -71.68535 -7.058315 FALSE  TRUE FALSE FALSE  TRUE
    #> 10 10 -0.30538839 -71.68615 -7.042494 FALSE FALSE FALSE  TRUE FALSE
    #> 11 11  1.51178117 -71.69045 -7.053219 FALSE  TRUE FALSE  TRUE  TRUE
    #> 12 12  0.38984324 -71.68278 -7.043211 FALSE FALSE FALSE FALSE  TRUE
    #> 13 13 -0.62124058 -71.69124 -7.053066 FALSE FALSE FALSE FALSE  TRUE
    #> 14 14 -2.21469989 -71.69510 -7.053325 FALSE FALSE FALSE FALSE FALSE
    #> 15 15  1.12493092 -71.69859 -7.050473 FALSE FALSE  TRUE FALSE FALSE
    #> 16 16 -0.04493361 -71.69801 -7.042156 FALSE FALSE FALSE FALSE  TRUE
    #> 17 17 -0.01619026 -71.69367 -7.042713 FALSE FALSE FALSE FALSE  TRUE
    #> 18 18  0.94383621 -71.68963 -7.052200 FALSE FALSE FALSE FALSE  TRUE
    #> 19 19  0.82122120 -71.68676 -7.044454 FALSE FALSE  TRUE FALSE  TRUE
    #> 20 20  0.59390132 -71.69186 -7.040788 FALSE FALSE FALSE FALSE  TRUE

## Data for `spifa()`

We extract the item as a matrix, and then provide it to a new
`data.frame` to later convert to `sf`:

``` r

# get items as a matrix
items <- subset(raw, select = grep("^item", names(raw))) |>
  as.matrix() |>
  unname()
# create data.frame and geo-reference it
data <- data.frame(items = I(items), wealth = raw$wealth, x = raw$lon, y = raw$lat) |>
    st_as_sf(coords = c("x", "y"), crs = 4326)
data
```

    #> Simple feature collection with 20 features and 2 fields
    #> Geometry type: POINT
    #> Dimension:     XY
    #> Bounding box:  xmin: -71.69953 ymin: -7.058315 xmax: -71.68278 ymax: -7.040788
    #> Geodetic CRS:  WGS 84
    #> First 10 features:
    #>    items.1 items.2 items.3 items.4 items.5     wealth                    geometry
    #> 1    FALSE   FALSE    TRUE    TRUE   FALSE -0.6264538 POINT (-71.68358 -7.041742)
    #> 2    FALSE   FALSE   FALSE    TRUE    TRUE  0.1836433 POINT (-71.68706 -7.054128)
    #> 3    FALSE   FALSE   FALSE   FALSE   FALSE -0.8356286 POINT (-71.68434 -7.050819)
    #> 4    FALSE    TRUE   FALSE   FALSE    TRUE  1.5952808 POINT (-71.68894 -7.053352)
    #> 5    FALSE   FALSE    TRUE    TRUE    TRUE  0.3295078 POINT (-71.68941 -7.046983)
    #> 6    FALSE   FALSE   FALSE   FALSE   FALSE -0.8204684  POINT (-71.68421 -7.05484)
    #> 7    FALSE   FALSE   FALSE   FALSE   FALSE  0.4874291 POINT (-71.69953 -7.050429)
    #> 8    FALSE   FALSE   FALSE    TRUE   FALSE  0.7383247 POINT (-71.69046 -7.044674)
    #> 9    FALSE    TRUE   FALSE   FALSE    TRUE  0.5757814 POINT (-71.68535 -7.058315)
    #> 10   FALSE   FALSE   FALSE    TRUE   FALSE -0.3053884 POINT (-71.68615 -7.042494)

Notice that we can easily access the response `items`:

``` r

data$items
```

    #>        [,1]  [,2]  [,3]  [,4]  [,5]
    #>  [1,] FALSE FALSE  TRUE  TRUE FALSE
    #>  [2,] FALSE FALSE FALSE  TRUE  TRUE
    #>  [3,] FALSE FALSE FALSE FALSE FALSE
    #>  [4,] FALSE  TRUE FALSE FALSE  TRUE
    #>  [5,] FALSE FALSE  TRUE  TRUE  TRUE
    #>  [6,] FALSE FALSE FALSE FALSE FALSE
    #>  [7,] FALSE FALSE FALSE FALSE FALSE
    #>  [8,] FALSE FALSE FALSE  TRUE FALSE
    #>  [9,] FALSE  TRUE FALSE FALSE  TRUE
    #> [10,] FALSE FALSE FALSE  TRUE FALSE
    #> [11,] FALSE  TRUE FALSE  TRUE  TRUE
    #> [12,] FALSE FALSE FALSE FALSE  TRUE
    #> [13,] FALSE FALSE FALSE FALSE  TRUE
    #> [14,] FALSE FALSE FALSE FALSE FALSE
    #> [15,] FALSE FALSE  TRUE FALSE FALSE
    #> [16,] FALSE FALSE FALSE FALSE  TRUE
    #> [17,] FALSE FALSE FALSE FALSE  TRUE
    #> [18,] FALSE FALSE FALSE FALSE  TRUE
    #> [19,] FALSE FALSE  TRUE FALSE  TRUE
    #> [20,] FALSE FALSE FALSE FALSE  TRUE

This dataset can directly be used to fit a model:

``` r

spifa(items ~ wealth, data = data, nfactors = 2, niter = 5)
```

    #> Item factor analysis model: spifa_pred
    #> Formula: items ~ wealth
    #> Dimensions: 20 respondents, 5 items, 2 latent factors, 2 spatial processes
    #> MCMC: 1 chain, iter = 5, thin = 1, samples = 5
    #> 
    #> Item model parameters:
    #>          mean median   sd   q10    q90 ess_bulk rhat
    #> c[1]   -1.478 -1.471 0.24 -1.71 -1.241       NA  0.9
    #> c[2]   -0.298 -0.353 0.20 -0.44 -0.091       NA  1.9
    #> c[3]   -0.651 -0.783 0.20 -0.79 -0.422       NA  0.9
    #> c[4]   -0.734 -0.705 0.30 -1.02 -0.488       NA  0.9
    #> c[5]   -0.039 -0.044 0.17 -0.19  0.111       NA  0.9
    #> A[1,1]  0.470  0.500 0.22  0.25  0.663       NA  1.9
    #> A[2,1]  0.770  0.753 0.16  0.63  0.917       NA  1.9
    #> A[3,1]  0.846  0.907 0.35  0.47  1.179       NA  0.9
    #> A[4,1] -0.031 -0.173 0.63 -0.55  0.625       NA  1.9
    #> A[5,1]  0.441  0.370 0.56 -0.11  1.010       NA  1.9
    #> A[2,2]  1.176  1.282 0.34  0.80  1.480       NA  1.9
    #> A[3,2]  0.035 -0.024 0.20 -0.15  0.252       NA  0.9
    #> A[4,2] -0.394 -0.343 0.38 -0.76 -0.070       NA  1.9
    #> A[5,2] -0.161 -0.254 0.27 -0.35  0.111       NA  1.9
    #> 
    #> Factor model parameters:
    #>              mean median      sd    q10   q90 ess_bulk rhat
    #> B[1,1]     0.1078 0.1789 0.33165 -0.250 0.402       NA  1.9
    #> B[1,2]     0.0023 0.1551 0.35755 -0.384 0.321       NA  1.9
    #> T[1,1]     0.4868 0.4828 0.01278  0.477 0.500       NA  0.9
    #> T[2,2]     0.6929 0.6908 0.00903  0.684 0.702       NA  1.9
    #> phi[1]     0.0471 0.0477 0.00097  0.046 0.048       NA  1.9
    #> phi[2]     0.0536 0.0539 0.00062  0.053 0.054       NA  0.9
    #> Corr[2,1] -0.0021 0.0023 0.02035 -0.023 0.018       NA  1.9
    #> 
    #> ess_bulk is the bulk effective sample size; rhat is the potential
    #> scale reduction factor on split chains (Rhat = 1 at convergence).
    #> Use summary() for the full set of statistics (incl. ess_tail).

Alternatively, you could also create the `sf` data first and then assign
the items:

``` r

data <- data.frame(wealth = raw$wealth, x = raw$lon, y = raw$lat) |>
    st_as_sf(coords = c("x", "y"), crs = 4326)
data$items <- items
data
```

    #> Simple feature collection with 20 features and 2 fields
    #> Geometry type: POINT
    #> Dimension:     XY
    #> Bounding box:  xmin: -71.69953 ymin: -7.058315 xmax: -71.68278 ymax: -7.040788
    #> Geodetic CRS:  WGS 84
    #> First 10 features:
    #>        wealth                    geometry items.1 items.2 items.3 items.4 items.5
    #> 1  -0.6264538 POINT (-71.68358 -7.041742)   FALSE   FALSE    TRUE    TRUE   FALSE
    #> 2   0.1836433 POINT (-71.68706 -7.054128)   FALSE   FALSE   FALSE    TRUE    TRUE
    #> 3  -0.8356286 POINT (-71.68434 -7.050819)   FALSE   FALSE   FALSE   FALSE   FALSE
    #> 4   1.5952808 POINT (-71.68894 -7.053352)   FALSE    TRUE   FALSE   FALSE    TRUE
    #> 5   0.3295078 POINT (-71.68941 -7.046983)   FALSE   FALSE    TRUE    TRUE    TRUE
    #> 6  -0.8204684  POINT (-71.68421 -7.05484)   FALSE   FALSE   FALSE   FALSE   FALSE
    #> 7   0.4874291 POINT (-71.69953 -7.050429)   FALSE   FALSE   FALSE   FALSE   FALSE
    #> 8   0.7383247 POINT (-71.69046 -7.044674)   FALSE   FALSE   FALSE    TRUE   FALSE
    #> 9   0.5757814 POINT (-71.68535 -7.058315)   FALSE    TRUE   FALSE   FALSE    TRUE
    #> 10 -0.3053884 POINT (-71.68615 -7.042494)   FALSE   FALSE   FALSE    TRUE   FALSE

``` r

spifa(items ~ wealth, data = data, nfactors = 2, niter = 5)
```

    #> Item factor analysis model: spifa_pred
    #> Formula: items ~ wealth
    #> Dimensions: 20 respondents, 5 items, 2 latent factors, 2 spatial processes
    #> MCMC: 1 chain, iter = 5, thin = 1, samples = 5
    #> 
    #> Item model parameters:
    #>          mean median   sd    q10     q90 ess_bulk rhat
    #> c[1]   -0.521 -0.620 1.15 -1.596  0.6837       NA  1.9
    #> c[2]   -0.395 -0.422 0.22 -0.586 -0.1786       NA  1.9
    #> c[3]   -0.680 -0.644 0.64 -1.319 -0.0379       NA  1.9
    #> c[4]   -0.269 -0.381 0.26 -0.506  0.0105       NA  0.9
    #> c[5]    0.343  0.310 0.14  0.240  0.4781       NA  0.9
    #> A[1,1]  0.719  0.623 0.26  0.570  0.9666       NA  1.9
    #> A[2,1]  0.733  0.877 0.24  0.464  0.9076       NA  1.9
    #> A[3,1] -0.418 -0.393 0.43 -0.853  0.0078       NA  1.9
    #> A[4,1]  0.072  0.084 0.13 -0.060  0.2057       NA  1.9
    #> A[5,1]  0.539  0.674 0.46  0.039  0.9192       NA  1.9
    #> A[2,2]  0.513  0.558 0.14  0.365  0.6409       NA  0.9
    #> A[3,2] -0.640 -0.627 0.30 -0.939 -0.3380       NA  1.9
    #> A[4,2]  0.338  0.604 0.59 -0.306  0.8515       NA  1.9
    #> A[5,2] -0.347 -0.385 0.23 -0.562 -0.1122       NA  1.9
    #> 
    #> Factor model parameters:
    #>             mean median      sd     q10   q90 ess_bulk rhat
    #> B[1,1]    0.0339 0.0249 0.36675 -0.3279 0.345       NA  0.9
    #> B[1,2]    0.0378 0.0188 0.19827 -0.1308 0.222       NA  1.9
    #> T[1,1]    0.5518 0.5554 0.01123  0.5399 0.560       NA  0.9
    #> T[2,2]    0.5037 0.5025 0.00435  0.5008 0.508       NA  0.9
    #> phi[1]    0.0487 0.0487 0.00079  0.0480 0.049       NA  1.9
    #> phi[2]    0.0512 0.0511 0.00082  0.0505 0.052       NA  1.9
    #> Corr[2,1] 0.0061 0.0043 0.01438 -0.0058 0.021       NA  1.9
    #> 
    #> ess_bulk is the bulk effective sample size; rhat is the potential
    #> scale reduction factor on split chains (Rhat = 1 at convergence).
    #> Use summary() for the full set of statistics (incl. ess_tail).

## Missing values

### Missing item responses

Individual missing item responses (`NA`) don’t need to be dropped or
imputed –
[`spifa()`](https://ErickChacon.github.io/spifa/reference/spifa.md)
handles them natively, treating a missing response as an unobserved
auxiliary variable to be sampled along with everything else, rather than
requiring complete rows:

``` r

data_na <- data
data_na$items[1:5, 1] <- NA

fit_na <- spifa(items ~ 1, data = data_na, nfactors = 2, ngp = 0, niter = 5)
fit_na
```

    #> Item factor analysis model: eifa
    #> Formula: items ~ 1
    #> Dimensions: 20 respondents, 5 items, 2 latent factors, 0 spatial processes
    #> MCMC: 1 chain, iter = 5, thin = 1, samples = 5
    #> 
    #> Item model parameters:
    #>         mean median   sd     q10   q90 ess_bulk rhat
    #> c[1]   -0.76 -0.759 0.41 -1.1637 -0.35       NA  1.9
    #> c[2]   -1.09 -1.114 0.34 -1.4281 -0.76       NA  1.9
    #> c[3]   -0.11 -0.083 0.29 -0.3795  0.18       NA  1.9
    #> c[4]   -0.51 -0.431 0.29 -0.8142 -0.24       NA  0.9
    #> c[5]    0.35  0.348 0.10  0.2492  0.45       NA  1.9
    #> A[1,1]  0.45  0.577 0.25  0.1766  0.67       NA  0.9
    #> A[2,1] -0.37 -0.461 0.25 -0.5951 -0.11       NA  1.9
    #> A[3,1]  0.25  0.284 0.27 -0.0021  0.48       NA  1.9
    #> A[4,1]  0.59  0.670 0.24  0.3540  0.74       NA  1.9
    #> A[5,1]  0.23  0.194 0.24  0.0111  0.47       NA  1.9
    #> A[2,2]  1.49  1.534 0.30  1.2155  1.72       NA  0.9
    #> A[3,2]  0.15  0.217 0.20 -0.0657  0.34       NA  0.9
    #> A[4,2]  0.28  0.413 0.38 -0.1132  0.54       NA  0.9
    #> A[5,2]  0.30  0.239 0.16  0.1605  0.47       NA  1.9
    #> 
    #> ess_bulk is the bulk effective sample size; rhat is the potential
    #> scale reduction factor on split chains (Rhat = 1 at convergence).
    #> Use summary() for the full set of statistics (incl. ess_tail).

As you can see in `Dimensions`, the missing responses are not removed
but rather modelled.

### Missing predictor values

Missing *predictor* values behave differently, and are worth calling
out: a `NA` in a predictor column drops that respondent from the whole
fit:

``` r

data_na <- data
data_na$wealth[1:5] <- NA

fit_na <- spifa(items ~ wealth, data = data_na, nfactors = 2, ngp = 0, niter = 5)
fit_na
```

    #> Item factor analysis model: cifa_pred
    #> Formula: items ~ wealth
    #> Dimensions: 15 respondents, 5 items, 2 latent factors, 0 spatial processes
    #> MCMC: 1 chain, iter = 5, thin = 1, samples = 5
    #> 
    #> Item model parameters:
    #>          mean median   sd   q10    q90 ess_bulk rhat
    #> c[1]   -0.286 -0.472 0.48 -0.71  0.249       NA  1.9
    #> c[2]   -0.724 -0.606 0.69 -1.39 -0.081       NA  1.9
    #> c[3]   -1.120 -1.276 0.40 -1.46 -0.673       NA  1.9
    #> c[4]   -0.614 -0.572 0.22 -0.84 -0.436       NA  1.9
    #> c[5]    0.214  0.299 0.14  0.06  0.311       NA  1.9
    #> A[1,1]  0.437  0.412 0.63 -0.19  1.082       NA  1.9
    #> A[2,1] -0.266 -0.233 0.18 -0.46 -0.111       NA  1.9
    #> A[3,1]  0.148 -0.040 0.38 -0.17  0.571       NA  1.9
    #> A[4,1]  0.099  0.151 0.24 -0.13  0.276       NA  1.9
    #> A[5,1] -0.601 -0.679 0.35 -0.90 -0.242       NA  1.9
    #> A[2,2]  0.547  0.559 0.31  0.25  0.815       NA  1.9
    #> A[3,2] -0.070 -0.022 0.14 -0.22  0.048       NA  0.9
    #> A[4,2] -0.224 -0.257 0.23 -0.44  0.019       NA  0.9
    #> A[5,2]  0.046 -0.070 0.27 -0.19  0.338       NA  1.9
    #> 
    #> Factor model parameters:
    #>             mean median    sd     q10   q90 ess_bulk rhat
    #> B[1,1]    -0.495  -0.45 0.329 -0.8420 -0.19       NA  1.9
    #> B[1,2]    -0.780  -0.91 0.247 -0.9809 -0.50       NA  0.9
    #> Corr[2,1]  0.016   0.01 0.022 -0.0031  0.04       NA  0.9
    #> 
    #> ess_bulk is the bulk effective sample size; rhat is the potential
    #> scale reduction factor on split chains (Rhat = 1 at convergence).
    #> Use summary() for the full set of statistics (incl. ess_tail).

As you can see in `Dimensions`, the 5 respondents with a missing
`wealth` value are dropped – 15 respondents remain, unlike the
item-response case above.
