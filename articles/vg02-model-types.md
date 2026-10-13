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
    #>            mean median   sd    q10    q90 ess_bulk rhat
    #> c[1]    -0.4327 -0.426 0.15 -0.632 -0.246    137.1  1.0
    #> c[2]    -0.5617 -0.536 0.22 -0.823 -0.304     37.2  1.1
    #> c[3]    -0.6383 -0.618 0.23 -0.984 -0.346     31.7  1.0
    #> c[4]    -0.3066 -0.305 0.14 -0.481 -0.128    118.5  1.0
    #> c[5]    -0.4407 -0.407 0.33 -0.918 -0.024      6.1  1.2
    #> c[6]    -0.0018  0.022 0.26 -0.353  0.330     23.6  1.1
    #> c[7]     0.1759  0.158 0.20 -0.059  0.425      9.3  1.1
    #> c[8]     0.1451  0.147 0.16 -0.048  0.340     70.7  1.0
    #> c[9]     0.5490  0.547 0.24  0.250  0.849     21.8  1.0
    #> c[10]    0.3631  0.362 0.24  0.035  0.676      4.6  1.2
    #> A[1,1]   0.6428  0.626 0.25  0.344  0.987     25.9  1.0
    #> A[2,1]   0.8755  0.831 0.32  0.523  1.332     13.1  1.1
    #> A[3,1]   1.3056  1.258 0.42  0.833  1.915     12.9  1.0
    #> A[4,1]   0.6098  0.612 0.27  0.297  0.930     29.0  1.0
    #> A[5,1]  -0.0376 -0.044 0.45 -0.616  0.570      7.8  1.2
    #> A[6,1]  -0.2582 -0.299 0.42 -0.777  0.305      7.4  1.2
    #> A[7,1]   0.5505  0.534 0.36  0.121  1.078      9.2  1.1
    #> A[8,1]   0.8162  0.820 0.29  0.462  1.166     21.2  1.0
    #> A[9,1]   1.1389  1.103 0.39  0.703  1.674      3.3  1.3
    #> A[10,1]  0.9873  0.929 0.51  0.384  1.811      2.6  1.3
    #> A[2,2]   0.4768  0.471 0.24  0.165  0.812     12.4  1.1
    #> A[3,2]   0.3165  0.345 0.33 -0.161  0.743      5.2  1.2
    #> A[4,2]   0.0630  0.055 0.24 -0.244  0.380      9.9  1.1
    #> A[5,2]   1.7342  1.611 0.61  1.088  2.546     11.0  1.1
    #> A[6,2]   1.7316  1.663 0.67  0.899  2.657      8.3  1.1
    #> A[7,2]   1.0515  1.021 0.29  0.700  1.434     21.3  1.0
    #> A[8,2]   0.2543  0.260 0.27 -0.111  0.584      5.2  1.2
    #> A[9,2]   0.9528  0.944 0.36  0.536  1.312     10.7  1.1
    #> A[10,2]  1.4111  1.347 0.44  0.896  2.006     12.9  1.0
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
    #>           mean median   sd     q10   q90 ess_bulk rhat
    #> c[1]    -0.578 -0.562 0.23 -0.8878 -0.31     36.9  1.0
    #> c[2]    -0.582 -0.584 0.21 -0.8350 -0.32     47.2  1.0
    #> c[3]    -0.719 -0.720 0.27 -1.0479 -0.35     26.5  1.1
    #> c[4]    -0.344 -0.335 0.17 -0.5775 -0.13     51.3  1.0
    #> c[5]    -0.285 -0.292 0.26 -0.6052  0.05     41.2  1.0
    #> c[6]     0.067  0.074 0.19 -0.1756  0.31     54.2  1.0
    #> c[7]     0.281  0.267 0.22  0.0055  0.58     40.5  1.0
    #> c[8]     0.211  0.199 0.21 -0.0516  0.49     41.4  1.0
    #> c[9]     0.703  0.666 0.27  0.4005  1.05     30.0  1.1
    #> c[10]    0.573  0.564 0.32  0.1819  1.00     29.5  1.1
    #> A[1,1]   0.341  0.313 0.31 -0.0170  0.79      2.2  1.5
    #> A[2,1]   0.434  0.422 0.24  0.1202  0.77      2.8  1.3
    #> A[3,1]   0.137  0.159 0.26 -0.1946  0.44     42.3  1.0
    #> A[5,1]   1.635  1.655 0.47  1.0603  2.22     19.4  1.0
    #> A[6,1]   1.336  1.260 0.45  0.8590  1.95      2.9  1.3
    #> A[7,1]   1.104  1.106 0.31  0.7405  1.44      6.8  1.1
    #> A[9,1]   1.034  1.001 0.32  0.6770  1.48      3.2  1.3
    #> A[10,1]  1.440  1.450 0.44  0.8891  1.97      3.8  1.2
    #> A[1,2]   1.408  1.320 0.73  0.6707  2.19     12.1  1.2
    #> A[2,2]   0.984  0.967 0.35  0.5439  1.43     12.1  1.0
    #> A[3,2]   0.818  0.786 0.38  0.3453  1.31     19.5  1.0
    #> A[9,2]   0.950  0.961 0.40  0.5284  1.43     11.6  1.0
    #> A[1,3]  -0.361 -0.227 0.53 -1.2113  0.15      9.8  1.1
    #> A[2,3]   0.314  0.350 0.29 -0.0640  0.64     31.4  1.0
    #> A[3,3]   1.134  1.099 0.35  0.7076  1.59     16.7  1.0
    #> A[4,3]   0.948  0.935 0.29  0.5885  1.29     35.8  1.1
    #> A[7,3]   0.647  0.662 0.23  0.3404  0.94     75.6  1.0
    #> A[8,3]   1.197  1.221 0.32  0.7962  1.58     41.1  1.0
    #> A[9,3]   0.571  0.602 0.43  0.0243  1.09      4.4  1.2
    #> A[10,3]  1.212  1.222 0.39  0.6675  1.73      5.7  1.1
    #> 
    #> Factor model parameters:
    #>            mean median   sd    q10  q90 ess_bulk rhat
    #> Corr[2,1] 0.043  0.038 0.14 -0.151 0.22      1.5  1.9
    #> Corr[3,1] 0.146  0.121 0.18 -0.084 0.40      9.0  1.1
    #> Corr[3,2] 0.321  0.259 0.18  0.135 0.69      4.9  1.3
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
    #>           mean median   sd    q10    q90 ess_bulk rhat
    #> c[1]    -0.436  -0.43 0.16 -0.612 -0.240     88.4  1.0
    #> c[2]    -0.467  -0.46 0.16 -0.674 -0.264    100.8  1.0
    #> c[3]    -0.456  -0.44 0.16 -0.653 -0.269    107.2  1.0
    #> c[4]    -0.283  -0.28 0.15 -0.479 -0.085    136.9  1.0
    #> c[5]    -0.274  -0.26 0.16 -0.472 -0.069    146.5  1.0
    #> c[6]     0.052   0.04 0.16 -0.133  0.268    129.1  1.0
    #> c[7]     0.232   0.24 0.17  0.023  0.439    114.2  1.0
    #> c[8]     0.139   0.13 0.14 -0.025  0.323    148.1  1.0
    #> c[9]     0.589   0.58 0.20  0.322  0.847     47.7  1.0
    #> c[10]    0.393   0.39 0.18  0.159  0.628    101.7  1.0
    #> A[1,1]   0.364   0.35 0.17  0.157  0.587     33.0  1.1
    #> A[2,1]  -0.046  -0.05 0.22 -0.331  0.236     53.9  1.0
    #> A[3,1]  -0.120  -0.14 0.22 -0.365  0.150     60.7  1.0
    #> A[5,1]  -1.166  -1.18 0.26 -1.475 -0.836     19.5  1.1
    #> A[6,1]  -0.949  -0.95 0.25 -1.268 -0.629     32.8  1.0
    #> A[7,1]  -0.872  -0.87 0.24 -1.175 -0.575     70.6  1.0
    #> A[9,1]  -0.491  -0.49 0.25 -0.777 -0.204     60.0  1.0
    #> A[10,1] -0.955  -0.94 0.29 -1.319 -0.634     45.0  1.0
    #> A[1,2]  -0.244  -0.24 0.37 -0.706  0.173      2.9  1.3
    #> A[2,2]   0.162   0.16 0.15 -0.028  0.350     23.8  1.1
    #> A[3,2]  -0.232  -0.24 0.27 -0.573  0.117     19.1  1.0
    #> A[9,2]  -0.280  -0.30 0.33 -0.663  0.114      5.3  1.2
    #> A[1,3]   0.464   0.47 0.25  0.146  0.784     60.5  1.0
    #> A[2,3]   0.699   0.71 0.23  0.414  0.982     56.5  1.0
    #> A[3,3]   0.556   0.56 0.14  0.369  0.722    164.1  1.0
    #> A[4,3]   0.477   0.47 0.19  0.251  0.707     92.8  1.0
    #> A[7,3]   0.446   0.44 0.20  0.209  0.710     75.3  1.0
    #> A[8,3]   0.642   0.63 0.18  0.419  0.871     74.0  1.0
    #> A[9,3]   0.711   0.70 0.26  0.390  1.033     56.0  1.0
    #> A[10,3]  0.666   0.66 0.22  0.408  0.973     66.2  1.0
    #> 
    #> Factor model parameters:
    #>             mean median   sd   q10   q90 ess_bulk rhat
    #> B[1,1]     6.218  6.276 1.14  4.92  7.59     10.1  1.1
    #> B[2,1]     0.091  0.072 1.05 -1.26  1.36     55.0  1.0
    #> B[1,2]     1.272  1.235 2.36 -1.65  4.56     16.9  1.2
    #> B[2,2]    -0.040  0.153 2.34 -3.01  2.71     20.2  1.1
    #> B[1,3]    -2.760 -2.705 1.52 -4.67 -0.84     40.2  1.0
    #> B[2,3]    -0.226 -0.279 1.28 -1.85  1.45     67.4  1.0
    #> Corr[2,1]  0.070  0.103 0.19 -0.15  0.26     11.1  1.1
    #> Corr[3,1] -0.246 -0.208 0.13 -0.39 -0.11     20.4  1.0
    #> Corr[3,2] -0.437 -0.447 0.16 -0.63 -0.16      7.8  1.1
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
    #>            mean  median   sd    q10     q90 ess_bulk rhat
    #> c[1]    -0.4108 -0.4102 0.16 -0.606 -0.2025     92.8  1.0
    #> c[2]    -0.5238 -0.5219 0.17 -0.739 -0.3159     83.5  1.0
    #> c[3]    -0.5969 -0.5802 0.19 -0.832 -0.3775     82.7  1.0
    #> c[4]    -0.3062 -0.3039 0.15 -0.509 -0.1080     90.1  1.0
    #> c[5]    -0.3176 -0.3310 0.20 -0.563 -0.0579     38.5  1.0
    #> c[6]     0.0527  0.0340 0.17 -0.144  0.2827     76.8  1.0
    #> c[7]     0.2360  0.2327 0.16  0.030  0.4546     57.9  1.0
    #> c[8]     0.1659  0.1621 0.18 -0.052  0.3853     93.7  1.0
    #> c[9]     0.6982  0.6597 0.27  0.384  1.0851     24.6  1.1
    #> c[10]    0.4167  0.3970 0.22  0.158  0.7417     42.2  1.0
    #> A[1,1]   0.6591  0.6656 0.27  0.324  1.0018    110.1  1.0
    #> A[2,1]  -0.6127 -0.6007 0.33 -1.054 -0.1898     24.4  1.1
    #> A[3,1]  -0.1575 -0.1580 0.41 -0.689  0.3714     77.2  1.0
    #> A[5,1]  -1.2883 -1.3432 0.41 -1.744 -0.8491      3.1  1.3
    #> A[6,1]  -1.1843 -1.1850 0.49 -1.782 -0.6149      2.2  1.4
    #> A[7,1]  -0.9580 -0.9238 0.39 -1.443 -0.6203      4.0  1.2
    #> A[9,1]  -0.3980 -0.4094 0.56 -1.080  0.3085     34.6  1.0
    #> A[10,1] -1.1507 -1.1825 0.39 -1.613 -0.6891      9.3  1.0
    #> A[1,2]  -0.7189 -0.7874 0.51 -1.247 -0.0017      5.1  1.1
    #> A[2,2]   0.4527  0.4491 0.29  0.066  0.8418     16.1  1.1
    #> A[3,2]  -0.0029  0.0097 0.48 -0.633  0.6063      8.3  1.1
    #> A[9,2]  -0.6908 -0.6900 0.58 -1.485  0.0369     19.9  1.1
    #> A[1,3]   0.4742  0.4706 0.23  0.172  0.7802     17.6  1.1
    #> A[2,3]   0.7622  0.7576 0.21  0.493  1.0092     27.3  1.1
    #> A[3,3]   1.1314  1.1347 0.20  0.876  1.3827     69.1  1.0
    #> A[4,3]   0.6599  0.6687 0.23  0.354  0.9354      7.2  1.1
    #> A[7,3]   0.5872  0.5858 0.21  0.334  0.8415     12.2  1.1
    #> A[8,3]   0.9465  0.9449 0.25  0.627  1.2529      5.2  1.2
    #> A[9,3]   1.3049  1.2550 0.41  0.824  1.8185     11.6  1.1
    #> A[10,3]  1.0439  1.0464 0.26  0.723  1.4011     12.7  1.1
    #> 
    #> Factor model parameters:
    #>             mean median     sd    q10    q90 ess_bulk rhat
    #> T[1,1]     0.980  0.943 0.1260  0.847 1.1832      3.7  1.2
    #> T[2,1]     0.931  0.935 0.0857  0.818 1.0295     20.3  1.1
    #> T[3,2]     0.967  0.976 0.0963  0.840 1.0878     12.1  1.0
    #> phi[1]     0.043  0.042 0.0044  0.036 0.0485     11.1  1.1
    #> phi[2]     0.047  0.047 0.0062  0.039 0.0558     17.2  1.1
    #> Corr[2,1] -0.143 -0.150 0.1176 -0.288 0.0063      6.5  1.1
    #> Corr[3,1]  0.106  0.040 0.2070 -0.117 0.4123     18.7  1.1
    #> Corr[3,2]  0.049  0.029 0.2443 -0.246 0.3741     20.2  1.3
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
    #>           mean median   sd     q10     q90 ess_bulk rhat
    #> c[1]    -0.631 -0.609 0.27 -0.9807 -0.3099      4.4  1.2
    #> c[2]    -0.534 -0.525 0.17 -0.7716 -0.3210     44.0  1.0
    #> c[3]    -0.561 -0.575 0.20 -0.8209 -0.3027     66.1  1.0
    #> c[4]    -0.344 -0.334 0.20 -0.5809 -0.1155     31.4  1.0
    #> c[5]    -0.282 -0.284 0.21 -0.5496 -0.0069     31.4  1.0
    #> c[6]     0.055  0.060 0.18 -0.1884  0.2854     92.3  1.0
    #> c[7]     0.248  0.248 0.20  0.0072  0.4947     50.5  1.0
    #> c[8]     0.189  0.173 0.19 -0.0481  0.4331     80.9  1.0
    #> c[9]     0.698  0.683 0.23  0.4161  0.9885     27.7  1.0
    #> c[10]    0.474  0.470 0.24  0.1668  0.7957     48.2  1.0
    #> A[1,1]   0.426  0.411 0.21  0.1733  0.6864     35.7  1.0
    #> A[2,1]   0.463  0.455 0.25  0.1540  0.7820      8.8  1.1
    #> A[3,1]   0.301  0.276 0.22  0.0537  0.5761     28.8  1.0
    #> A[5,1]   1.354  1.345 0.32  0.9339  1.7900     26.8  1.0
    #> A[6,1]   1.199  1.166 0.33  0.8004  1.6159     24.5  1.1
    #> A[7,1]   1.136  1.082 0.35  0.7420  1.6389     30.2  1.0
    #> A[9,1]   1.098  1.090 0.28  0.7554  1.4541     27.4  1.0
    #> A[10,1]  1.374  1.272 0.43  0.9224  1.9533      4.3  1.2
    #> A[1,2]   1.578  1.481 0.62  0.8921  2.5351     12.0  1.0
    #> A[2,2]   0.746  0.725 0.24  0.4674  1.0469     36.3  1.0
    #> A[3,2]   0.615  0.587 0.30  0.2416  1.0433     19.1  1.1
    #> A[9,2]   0.864  0.849 0.35  0.4535  1.3111     17.6  1.0
    #> A[1,3]   0.030  0.036 0.33 -0.3755  0.4294     14.8  1.1
    #> A[2,3]   0.546  0.539 0.29  0.2034  0.9263     42.4  1.1
    #> A[3,3]   0.975  0.992 0.28  0.6248  1.3002     39.7  1.0
    #> A[4,3]   0.995  0.915 0.40  0.5916  1.5380     17.1  1.0
    #> A[7,3]   0.637  0.638 0.25  0.3140  0.9524     64.3  1.0
    #> A[8,3]   1.195  1.195 0.38  0.7728  1.7007     25.6  1.1
    #> A[9,3]   0.796  0.809 0.33  0.3532  1.2036     24.4  1.0
    #> A[10,3]  1.053  1.057 0.43  0.5549  1.6315      8.4  1.1
    #> 
    #> Factor model parameters:
    #>             mean median     sd    q10    q90 ess_bulk rhat
    #> B[1,1]    -0.449 -0.442 0.1260 -0.612 -0.300    166.8  1.0
    #> B[1,2]     0.129  0.129 0.1257 -0.031  0.288    174.1  1.0
    #> B[1,3]    -0.240 -0.244 0.1261 -0.398 -0.073     35.3  1.0
    #> T[1,1]     0.947  0.944 0.0992  0.826  1.054     26.1  1.1
    #> T[2,2]     0.884  0.862 0.1069  0.763  1.027      2.9  1.3
    #> T[3,3]     0.886  0.884 0.1168  0.737  1.045      3.8  1.3
    #> phi[1]     0.049  0.049 0.0077  0.039  0.057     11.4  1.1
    #> phi[2]     0.043  0.042 0.0050  0.038  0.048      2.5  1.4
    #> phi[3]     0.050  0.050 0.0032  0.046  0.055     10.7  1.0
    #> Corr[2,1] -0.103 -0.081 0.0848 -0.219 -0.018     31.7  1.0
    #> Corr[3,1] -0.074 -0.046 0.1215 -0.237  0.063     21.7  1.1
    #> Corr[3,2] -0.081 -0.092 0.0877 -0.188  0.043      2.2  1.5
    #> 
    #> ess_bulk is the bulk effective sample size; rhat is the potential
    #> scale reduction factor on split chains (Rhat = 1 at convergence).
    #> Use summary() for the full set of statistics (incl. ess_tail).

Notice that the predictor effects `B` is now included in the samples.
