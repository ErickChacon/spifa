# Model types

[`spifa()`](https://ErickChacon.github.io/spifa/reference/spifa.md) can
fit five model types that are determined automatically from `formula`
and the class of `data` (see
[`?spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md)):

- Non-spatial:
  - `eifa`: Exploratory IFA
  - `cifa`: Confirmatory IFA
  - `cifa_pred`: CIFA with predictors
- Spatial:
  - `spifa`: Spatial IFA
  - `spifa_pred`: SPIFA with predictors

The model type is determined based on the following:

1.  Exploratory structure (`eifa`) is used when no `constraints` are
    provided.
2.  Models with predictors (`cifa_pred`/`spifa_pred`) are activated if
    terms are included on the right-hand side of `formula`.
3.  Spatial models (`spifa`/`spifa_pred`) are only activated when the
    `data` is geo-referenced as an `sf` object, otherwise the
    non-spatial models are activated.
4.  Setting `ngp = 0` forces the model to be non-spatial
    (`eifa`/`cifa`/`cifa_pred`).

We will showcase the five model types using the `ipixuna` data.

``` r

library(spifa)

data(ipixuna)
nitems <- ncol(ipixuna$items)
```

## Exploratory item factor analysis (EIFA)

The simplest case is the `eifa` model, the user does not need to provide
`constraints` neither a spatial `data`. The constraints will be
internally defined as a lower-triangular matrix; and in case `data` is
an `sf` object, the spatial functionality can be dropped using
`ngp = 0`:

``` r

set.seed(123)
samples <- spifa(items ~ 1, data = ipixuna, nfactors = 2, ngp = 0, niter = 500)
samples
```

    #> Item factor analysis model: eifa
    #> Formula: items ~ 1
    #> Dimensions: 100 respondents, 10 items, 2 latent factors, 0 spatial processes
    #> MCMC: 1 chain, iter = 500, thin = 1, samples = 500
    #> 
    #> Item model parameters:
    #>            mean median   sd    q10     q90 ess_bulk rhat
    #> c[1]    -0.4102 -0.414 0.15 -0.590 -0.2106    121.8  1.0
    #> c[2]    -0.5215 -0.517 0.18 -0.754 -0.2975     57.4  1.0
    #> c[3]    -0.5718 -0.559 0.22 -0.851 -0.2868     61.6  1.0
    #> c[4]    -0.2970 -0.301 0.16 -0.499 -0.0884    118.0  1.0
    #> c[5]    -0.2675 -0.268 0.20 -0.514  0.0018     46.6  1.0
    #> c[6]     0.1068  0.088 0.23 -0.148  0.3979     50.6  1.0
    #> c[7]     0.2626  0.259 0.19  0.018  0.5172     45.4  1.0
    #> c[8]     0.2074  0.199 0.21 -0.035  0.4479     45.9  1.0
    #> c[9]     0.7201  0.684 0.26  0.427  1.0441     45.3  1.0
    #> c[10]    0.5153  0.519 0.23  0.240  0.7878     41.7  1.1
    #> A[1,1]   0.6013  0.589 0.23  0.326  0.9009      5.5  1.2
    #> A[2,1]   0.8696  0.839 0.34  0.475  1.2985     13.9  1.1
    #> A[3,1]   1.2459  1.219 0.39  0.787  1.7374     14.6  1.1
    #> A[4,1]   0.6308  0.619 0.25  0.361  0.9279     36.6  1.1
    #> A[5,1]  -0.0097  0.078 0.41 -0.654  0.4757     12.0  1.1
    #> A[6,1]  -0.2310 -0.207 0.35 -0.692  0.2276     19.6  1.0
    #> A[7,1]   0.5565  0.550 0.29  0.204  0.9338      3.8  1.2
    #> A[8,1]   1.0440  0.900 0.59  0.566  1.8901      2.8  1.3
    #> A[9,1]   1.3692  1.201 0.68  0.692  2.6265     11.3  1.1
    #> A[10,1]  1.0474  1.070 0.39  0.506  1.5035      4.2  1.2
    #> A[2,2]   0.4742  0.453 0.20  0.229  0.7346     38.6  1.0
    #> A[3,2]   0.3524  0.345 0.27  0.016  0.7145     10.6  1.1
    #> A[4,2]   0.0648  0.055 0.22 -0.203  0.3401     13.3  1.1
    #> A[5,2]   1.6031  1.531 0.46  1.102  2.2714     13.9  1.0
    #> A[6,2]   1.6781  1.601 0.63  0.923  2.5744     12.1  1.0
    #> A[7,2]   1.0608  1.050 0.29  0.687  1.4535     13.3  1.1
    #> A[8,2]   0.2790  0.273 0.27 -0.054  0.6172      7.2  1.1
    #> A[9,2]   1.1263  1.104 0.37  0.684  1.6746     27.9  1.0
    #> A[10,2]  1.5646  1.542 0.39  1.099  2.0739     17.6  1.0
    #> 
    #> ess_bulk is the bulk effective sample size; rhat is the potential
    #> scale reduction factor on split chains (Rhat = 1 at convergence).
    #> Use summary() for the full set of statistics (incl. ess_tail).

Notice how the constraints have been defined internally:

``` r

attr(samples, "fit_args")$constrain_L
```

    #>       [,1] [,2]
    #>  [1,]    1    0
    #>  [2,]    1    1
    #>  [3,]    1    1
    #>  [4,]    1    1
    #>  [5,]    1    1
    #>  [6,]    1    1
    #>  [7,]    1    1
    #>  [8,]    1    1
    #>  [9,]    1    1
    #> [10,]    1    1

EIFA also never estimates the correlation between the latent factors, so
only the easiness (c) and discrimination (A) parameters are estimated.

## Confirmatory item factor analysis (CIFA)

Supplying your own discrimination `constraints` leads to a `cifa` model:

``` r

nfactors <- 3
A <- matrix(1, nitems, nfactors)
A[c(4, 8), 1] <- 0
A[c(4, 5, 6, 7, 8, 10), 2] <- 0
A[c(5, 6), 3] <- 0

samples <- spifa(
  items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0, niter = 500,
  constraints = list(discrimination = A)
)
samples
```

    #> Item factor analysis model: cifa
    #> Formula: items ~ 1
    #> Dimensions: 100 respondents, 10 items, 3 latent factors, 0 spatial processes
    #> MCMC: 1 chain, iter = 500, thin = 1, samples = 500
    #> 
    #> Item model parameters:
    #>           mean median   sd     q10    q90 ess_bulk rhat
    #> c[1]    -0.639 -0.624 0.25 -0.9412 -0.338       49  1.0
    #> c[2]    -0.509 -0.494 0.18 -0.7397 -0.312       86  1.0
    #> c[3]    -0.712 -0.709 0.27 -1.0354 -0.385       48  1.0
    #> c[4]    -0.348 -0.347 0.15 -0.5406 -0.161       72  1.0
    #> c[5]    -0.381 -0.374 0.24 -0.6974 -0.099       51  1.0
    #> c[6]     0.026  0.025 0.17 -0.1868  0.242      102  1.0
    #> c[7]     0.218  0.209 0.18 -0.0057  0.482       50  1.0
    #> c[8]     0.143  0.148 0.20 -0.1105  0.397       48  1.0
    #> c[9]     0.802  0.747 0.32  0.4527  1.215       23  1.0
    #> c[10]    0.500  0.479 0.27  0.1580  0.838       47  1.0
    #> A[1,1]   0.474  0.504 0.30  0.0713  0.840       18  1.0
    #> A[2,1]   0.367  0.379 0.24  0.0567  0.668       33  1.0
    #> A[3,1]   0.217  0.211 0.28 -0.1100  0.590       36  1.0
    #> A[5,1]   1.832  1.808 0.45  1.3315  2.433       23  1.0
    #> A[6,1]   1.140  1.143 0.28  0.8016  1.500       13  1.1
    #> A[7,1]   1.061  1.045 0.32  0.6408  1.485       27  1.0
    #> A[9,1]   1.383  1.372 0.42  0.8790  1.928       27  1.1
    #> A[10,1]  1.577  1.588 0.44  1.0507  2.095       17  1.0
    #> A[1,2]  -1.642 -1.700 0.45 -2.1196 -1.124       32  1.0
    #> A[2,2]  -0.281 -0.296 0.25 -0.6011  0.044       30  1.0
    #> A[3,2]  -0.757 -0.756 0.35 -1.1897 -0.318       20  1.0
    #> A[9,2]  -1.252 -1.233 0.48 -1.8385 -0.682       18  1.1
    #> A[1,3]  -0.422 -0.463 0.45 -0.9941  0.166        9  1.1
    #> A[2,3]   0.484  0.489 0.29  0.1108  0.829       33  1.0
    #> A[3,3]   1.153  1.143 0.28  0.8127  1.518       19  1.1
    #> A[4,3]   0.763  0.753 0.24  0.4767  1.072       35  1.0
    #> A[7,3]   0.602  0.604 0.23  0.3013  0.907       36  1.0
    #> A[8,3]   1.305  1.340 0.43  0.7765  1.869       29  1.0
    #> A[9,3]   0.670  0.661 0.38  0.2168  1.213       14  1.0
    #> A[10,3]  1.467  1.436 0.45  0.8870  2.072       13  1.0
    #> 
    #> Factor model parameters:
    #>             mean  median   sd    q10   q90 ess_bulk rhat
    #> Corr[2,1]  0.015  0.0034 0.13 -0.182 0.186      9.6  1.1
    #> Corr[3,1]  0.147  0.1048 0.16 -0.016 0.425      7.6  1.1
    #> Corr[3,2] -0.233 -0.2421 0.18 -0.483 0.063      6.2  1.1
    #> 
    #> ess_bulk is the bulk effective sample size; rhat is the potential
    #> scale reduction factor on split chains (Rhat = 1 at convergence).
    #> Use summary() for the full set of statistics (incl. ess_tail).

Unlike `eifa`, the `cifa` model estimates the correlation (Corr) between
the latent factors.

## CIFA model with predictors

Adding predictor terms to the right-hand side of `formula` generates a
`cifa_pred`, which models the latent factors in terms of the predictor
terms:

``` r

samples <- spifa(
  items ~ poly(wealth, 2), data = ipixuna, nfactors = nfactors, ngp = 0, niter = 500,
  constraints = list(discrimination = A)
)
samples
```

    #> Item factor analysis model: cifa_pred
    #> Formula: items ~ poly(wealth, 2)
    #> Dimensions: 100 respondents, 10 items, 3 latent factors, 0 spatial processes
    #> MCMC: 1 chain, iter = 500, thin = 1, samples = 500
    #> 
    #> Item model parameters:
    #>           mean  median   sd    q10    q90 ess_bulk rhat
    #> c[1]    -0.419 -0.4127 0.15 -0.605 -0.243    127.1  1.0
    #> c[2]    -0.448 -0.4339 0.15 -0.653 -0.254    123.9  1.0
    #> c[3]    -0.412 -0.4142 0.14 -0.599 -0.222    131.9  1.0
    #> c[4]    -0.267 -0.2743 0.13 -0.435 -0.094    194.0  1.0
    #> c[5]    -0.249 -0.2525 0.16 -0.460 -0.041    105.8  1.0
    #> c[6]     0.046  0.0459 0.15 -0.142  0.226    115.1  1.0
    #> c[7]     0.255  0.2547 0.17  0.038  0.478    118.5  1.0
    #> c[8]     0.157  0.1496 0.14 -0.022  0.335    180.6  1.0
    #> c[9]     0.582  0.5792 0.18  0.354  0.814     36.8  1.0
    #> c[10]    0.428  0.4284 0.20  0.157  0.689     59.9  1.0
    #> A[1,1]   0.374  0.3834 0.18  0.147  0.586     47.9  1.0
    #> A[2,1]   0.428  0.4196 0.27  0.102  0.796     41.4  1.0
    #> A[3,1]   0.550  0.5758 0.32  0.140  0.915     16.5  1.1
    #> A[5,1]   1.037  1.0138 0.29  0.739  1.390     16.7  1.1
    #> A[6,1]   0.785  0.7982 0.23  0.544  1.059     55.6  1.0
    #> A[7,1]   0.784  0.7731 0.27  0.460  1.118     56.8  1.0
    #> A[9,1]   0.730  0.7208 0.28  0.377  1.069     56.5  1.0
    #> A[10,1]  0.912  0.9268 0.26  0.574  1.248     32.7  1.0
    #> A[1,2]   0.567  0.6136 0.41  0.079  1.025     13.6  1.1
    #> A[2,2]   0.395  0.3975 0.18  0.160  0.621      8.8  1.1
    #> A[3,2]   0.180  0.2980 0.51 -0.447  0.734      2.3  1.4
    #> A[9,2]   0.206  0.2353 0.41 -0.351  0.704      4.7  1.2
    #> A[1,3]  -0.206 -0.1908 0.28 -0.609  0.142     17.2  1.1
    #> A[2,3]  -0.379 -0.3777 0.25 -0.713 -0.067     35.6  1.0
    #> A[3,3]   0.019  0.0081 0.16 -0.165  0.237    133.0  1.0
    #> A[4,3]  -0.473 -0.4598 0.19 -0.729 -0.243     68.5  1.0
    #> A[7,3]  -0.444 -0.4405 0.22 -0.744 -0.153     69.0  1.0
    #> A[8,3]  -0.622 -0.6213 0.20 -0.881 -0.356     55.2  1.0
    #> A[9,3]  -0.623 -0.6187 0.27 -0.986 -0.281     36.5  1.0
    #> A[10,3] -0.599 -0.5902 0.27 -0.949 -0.262     55.9  1.0
    #> 
    #> Factor model parameters:
    #>            mean median   sd   q10   q90 ess_bulk rhat
    #> B[1,1]    -6.03  -6.04 1.56 -7.92 -4.20     15.7  1.1
    #> B[2,1]     0.46   0.51 1.14 -1.04  1.80     66.1  1.0
    #> B[1,2]     4.14   4.21 2.09  1.60  6.76     26.9  1.0
    #> B[2,2]     0.75   0.65 1.51 -1.11  2.75     39.6  1.0
    #> B[1,3]     3.16   2.98 1.88  0.80  5.48     25.4  1.0
    #> B[2,3]     0.77   0.76 1.36 -0.97  2.47     61.3  1.0
    #> Corr[2,1] -0.47  -0.49 0.18 -0.62 -0.19      5.2  1.1
    #> Corr[3,1] -0.35  -0.34 0.14 -0.53 -0.18     18.5  1.1
    #> Corr[3,2]  0.17   0.14 0.24 -0.10  0.51      2.9  1.3
    #> 
    #> ess_bulk is the bulk effective sample size; rhat is the potential
    #> scale reduction factor on split chains (Rhat = 1 at convergence).
    #> Use summary() for the full set of statistics (incl. ess_tail).

## Spatial item factor analysis (SPIFA)

If the `data` is an `sf` object and no predictor terms are included in
`formula`, the model type is `spifa` automatically. By default, it
includes a spatial term for each factor `ngp = nfactors`; but custom
number of `ngp` can be provided with a constraints `loading` matrix (T)
of size `nfactors` x `ngp`.

``` r

T_loading <- matrix(c(1, 1, 0, 0, 0, 1), nfactors, 2)
T_loading
```

    #>      [,1] [,2]
    #> [1,]    1    0
    #> [2,]    1    0
    #> [3,]    0    1

This structure means that 1 GP is shared between factor 1 and 2, and
another GP is introduced for factor 3:

``` r

samples <- spifa(
  items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 2, niter = 500,
  constraints = list(discrimination = A, loading = T_loading)
)
samples
```

    #> Item factor analysis model: spifa
    #> Formula: items ~ 1
    #> Dimensions: 100 respondents, 10 items, 3 latent factors, 2 spatial processes
    #> MCMC: 1 chain, iter = 500, thin = 1, samples = 500
    #> 
    #> Item model parameters:
    #>           mean median   sd    q10   q90 ess_bulk rhat
    #> c[1]    -0.414 -0.410 0.16 -0.640 -0.21    124.5  1.0
    #> c[2]    -0.508 -0.504 0.17 -0.726 -0.29     71.9  1.0
    #> c[3]    -0.519 -0.527 0.17 -0.740 -0.29     88.4  1.0
    #> c[4]    -0.323 -0.323 0.17 -0.539 -0.12     87.0  1.0
    #> c[5]    -0.265 -0.261 0.19 -0.512 -0.02     64.7  1.0
    #> c[6]     0.081  0.074 0.19 -0.159  0.34     69.7  1.0
    #> c[7]     0.260  0.258 0.18  0.037  0.49     10.5  1.1
    #> c[8]     0.177  0.154 0.20 -0.059  0.43     50.3  1.0
    #> c[9]     0.739  0.726 0.25  0.443  1.08      7.9  1.1
    #> c[10]    0.472  0.461 0.25  0.173  0.77     21.4  1.1
    #> A[1,1]   0.512  0.531 0.36  0.065  0.96    134.2  1.0
    #> A[2,1]  -0.188 -0.179 0.35 -0.650  0.27     64.7  1.0
    #> A[3,1]   0.108  0.106 0.56 -0.646  0.85     55.7  1.0
    #> A[5,1]   1.402  1.388 0.29  1.069  1.77     33.7  1.0
    #> A[6,1]   1.286  1.284 0.36  0.850  1.77     28.5  1.0
    #> A[7,1]   1.066  1.051 0.30  0.716  1.47     35.6  1.0
    #> A[9,1]   0.762  0.735 0.52  0.134  1.44     69.9  1.0
    #> A[10,1]  1.429  1.387 0.39  0.991  1.90     20.1  1.1
    #> A[1,2]  -0.346 -0.334 0.43 -0.904  0.17     67.7  1.0
    #> A[2,2]   0.682  0.676 0.33  0.241  1.12     81.0  1.0
    #> A[3,2]   0.173  0.204 0.57 -0.605  0.87     40.6  1.0
    #> A[9,2]   0.516  0.519 0.64 -0.328  1.25     35.5  1.1
    #> A[1,3]   0.401  0.396 0.20  0.155  0.65     56.1  1.0
    #> A[2,3]   0.650  0.647 0.23  0.369  0.93     46.0  1.0
    #> A[3,3]   1.058  1.036 0.24  0.753  1.38     31.2  1.0
    #> A[4,3]   0.763  0.709 0.31  0.417  1.16     25.6  1.0
    #> A[7,3]   0.546  0.533 0.22  0.287  0.82     34.2  1.0
    #> A[8,3]   1.100  1.008 0.45  0.655  1.67     12.2  1.1
    #> A[9,3]   1.264  1.227 0.41  0.802  1.78     19.0  1.0
    #> A[10,3]  1.085  1.071 0.31  0.709  1.46      3.8  1.2
    #> 
    #> Factor model parameters:
    #>            mean median     sd    q10   q90 ess_bulk rhat
    #> T[1,1]    0.919  0.915 0.0630  0.832 1.006     14.5  1.1
    #> T[2,1]    0.912  0.908 0.0667  0.814 1.011     12.9  1.1
    #> T[3,2]    0.963  0.951 0.0671  0.882 1.047     22.3  1.1
    #> phi[1]    0.049  0.049 0.0032  0.045 0.052      4.4  1.2
    #> phi[2]    0.055  0.055 0.0041  0.049 0.060     30.1  1.0
    #> Corr[2,1] 0.115  0.144 0.1652 -0.069 0.326     16.0  1.0
    #> Corr[3,1] 0.254  0.257 0.0768  0.134 0.336     24.3  1.1
    #> Corr[3,2] 0.213  0.218 0.1049  0.030 0.331     31.5  1.0
    #> 
    #> ess_bulk is the bulk effective sample size; rhat is the potential
    #> scale reduction factor on split chains (Rhat = 1 at convergence).
    #> Use summary() for the full set of statistics (incl. ess_tail).

Notice that the `loading` parameter is also available in the samples.

## SPIFA model with predictors

A `spifa_pred` model is simply defined by adding predictor terms to
`formula`.

``` r

samples <- spifa(
  items ~ wealth, data = ipixuna, nfactors = nfactors, niter = 500,
  constraints = list(discrimination = A)
)
samples
```

    #> Item factor analysis model: spifa_pred
    #> Formula: items ~ wealth
    #> Dimensions: 100 respondents, 10 items, 3 latent factors, 3 spatial processes
    #> MCMC: 1 chain, iter = 500, thin = 1, samples = 500
    #> 
    #> Item model parameters:
    #>             mean median   sd      q10   q90 ess_bulk rhat
    #> c[1]    -0.57260 -0.528 0.27 -0.93231 -0.29     33.7  1.0
    #> c[2]    -0.56493 -0.563 0.19 -0.80248 -0.32     79.9  1.0
    #> c[3]    -0.58779 -0.577 0.21 -0.87309 -0.32     65.5  1.0
    #> c[4]    -0.36991 -0.355 0.18 -0.62738 -0.15     48.5  1.0
    #> c[5]    -0.35614 -0.346 0.20 -0.61555 -0.11     45.6  1.0
    #> c[6]     0.00062  0.015 0.19 -0.22678  0.22     54.9  1.0
    #> c[7]     0.19205  0.197 0.22 -0.10039  0.47     35.1  1.0
    #> c[8]     0.19514  0.187 0.21 -0.06398  0.46     52.3  1.0
    #> c[9]     0.56597  0.559 0.24  0.26629  0.87     34.9  1.0
    #> c[10]    0.39335  0.403 0.26  0.04838  0.72     37.9  1.0
    #> A[1,1]   0.23607  0.234 0.19  0.00078  0.47     44.3  1.0
    #> A[2,1]   0.34453  0.329 0.19  0.12458  0.58     71.8  1.0
    #> A[3,1]   0.20850  0.212 0.20 -0.05182  0.46     25.2  1.1
    #> A[5,1]   1.24924  1.227 0.35  0.82162  1.66     14.7  1.1
    #> A[6,1]   1.07571  1.042 0.34  0.68615  1.48      5.8  1.1
    #> A[7,1]   0.95125  0.944 0.24  0.63978  1.27     27.1  1.1
    #> A[9,1]   0.80681  0.780 0.23  0.53003  1.13     29.0  1.0
    #> A[10,1]  1.10507  1.093 0.28  0.78731  1.46     36.6  1.0
    #> A[1,2]   1.11162  1.055 0.40  0.64652  1.68     13.9  1.2
    #> A[2,2]   0.66447  0.651 0.22  0.38657  0.94     50.6  1.0
    #> A[3,2]   0.48612  0.462 0.26  0.16289  0.84     38.4  1.0
    #> A[9,2]   0.71997  0.715 0.26  0.39855  1.02     44.1  1.0
    #> A[1,3]   0.18727  0.179 0.31 -0.19703  0.56      6.1  1.2
    #> A[2,3]   0.55237  0.569 0.26  0.21271  0.89     36.2  1.0
    #> A[3,3]   1.02677  1.023 0.26  0.70991  1.33     47.2  1.0
    #> A[4,3]   0.84582  0.839 0.30  0.45570  1.26     37.1  1.0
    #> A[7,3]   0.57838  0.580 0.21  0.32502  0.84     40.9  1.0
    #> A[8,3]   1.17996  1.136 0.34  0.76787  1.64     23.0  1.0
    #> A[9,3]   0.80017  0.767 0.33  0.38408  1.24     26.6  1.1
    #> A[10,3]  0.98692  0.981 0.33  0.62181  1.43     40.1  1.0
    #> 
    #> Factor model parameters:
    #>             mean median     sd    q10    q90 ess_bulk rhat
    #> B[1,1]    -0.489 -0.476 0.1617 -0.701 -0.312      7.5  1.1
    #> B[1,2]     0.109  0.117 0.1271 -0.068  0.265    258.6  1.0
    #> B[1,3]    -0.251 -0.248 0.1320 -0.423 -0.100    167.8  1.0
    #> T[1,1]     1.008  0.946 0.1597  0.862  1.228      1.4  2.0
    #> T[2,2]     0.920  0.925 0.0663  0.833  1.011     10.6  1.0
    #> T[3,3]     0.938  0.936 0.1048  0.823  1.069     16.3  1.0
    #> phi[1]     0.054  0.055 0.0052  0.046  0.060      9.8  1.1
    #> phi[2]     0.055  0.055 0.0102  0.040  0.069      2.4  1.4
    #> phi[3]     0.054  0.052 0.0078  0.044  0.063      9.3  1.1
    #> Corr[2,1]  0.055  0.043 0.1257 -0.097  0.220     12.2  1.0
    #> Corr[3,1]  0.039  0.028 0.1606 -0.160  0.275      2.0  1.6
    #> Corr[3,2]  0.025  0.080 0.3392 -0.484  0.469      1.8  1.6
    #> 
    #> ess_bulk is the bulk effective sample size; rhat is the potential
    #> scale reduction factor on split chains (Rhat = 1 at convergence).
    #> Use summary() for the full set of statistics (incl. ess_tail).

Notice that the predictor effects `B` is now included in the samples.
