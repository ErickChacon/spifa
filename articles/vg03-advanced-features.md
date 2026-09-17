# Advanced features

This vignette goes deeper into two things
[`vignette("spifa")`](https://ErickChacon.github.io/spifa/articles/spifa.md)
and
[`vignette("vg02-model-types")`](https://ErickChacon.github.io/spifa/articles/vg02-model-types.md)
only touch briefly:

- how to control `constraints`/`priors` for every parameter group,
- how to continue sampling with
  [`update()`](https://rdrr.io/r/stats/update.html) when a fit hasn’t
  converged yet, and
- how to compare models using
  [`dic()`](https://ErickChacon.github.io/spifa/reference/dic.md).

``` r

library(spifa)

data(ipixuna)
nfactors <- 3
nitems <- ncol(ipixuna$items)
```

## Restrictions and priors

[`spifa()`](https://ErickChacon.github.io/spifa/reference/spifa.md)’s
parameter glossary (see
[`?spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md)) maps
six configurable parameter names – `easiness` (c), `discrimination` (A),
`effect` (B), `corr` (Corr), `loading` (T), `range` (phi). Two
[`spifa()`](https://ErickChacon.github.io/spifa/reference/spifa.md)
arguments use these names:

- `constraints`: fixes cells of a parameter to a value instead of
  estimating them. Only `discrimination`, `loading`, and `sd` support
  this (the model needs some fixed structure for identifiability).
- `priors`: sets the `initial` value and prior hyperparameters
  (`mean`/`sd`, or `eta` for `corr`) for a parameter that is estimated.
  Every one of the six names accepts this.

### `constraints`

- `discrimination` (nitems x nfactors matrix): which items inform which
  factors,
- `loading` (nfactors x ngp matrix): which Gaussian processes feed into
  which factors, and
- `sd`: fixed values for the latent factors’ residual standard deviation
  (length `nfactors`). Unlike `discrimination`/`loading`, this is not
  estimated, it is always fixed.

Define `discrimination` structure:

``` r

A <- matrix(1, nitems, nfactors)
A[c(4, 8), 1] <- 0
A[c(4, 5, 6, 7, 8, 10), 2] <- 0
A[c(5, 6), 3] <- 0
A
```

    #>       [,1] [,2] [,3]
    #>  [1,]    1    1    1
    #>  [2,]    1    1    1
    #>  [3,]    1    1    1
    #>  [4,]    0    0    1
    #>  [5,]    1    0    0
    #>  [6,]    1    0    0
    #>  [7,]    1    0    1
    #>  [8,]    0    0    1
    #>  [9,]    1    1    1
    #> [10,]    1    0    1

Define `loading` structure, factors 1 and 2 share one Gaussian process,
factor 3 gets its own:

``` r

ngp <- 2
T_loading <- matrix(c(1, 1, 0, 0, 0, 1), nfactors, ngp)
T_loading
```

    #>      [,1] [,2]
    #> [1,]    1    0
    #> [2,]    1    0
    #> [3,]    0    1

Fix the residual standard deviation of each factor to 0.2:

``` r

sd_fixed <- rep(0.2, nfactors)
sd_fixed
```

    #> [1] 0.2 0.2 0.2

### `priors`

Every element accepts `initial` (starting value for the sampler) and the
prior hyperparameters `mean`/`sd`. The parameter `corr` is the
exception: it takes `initial` and `eta` (the LKJ shape parameter)
instead.

By default,
[`spifa()`](https://ErickChacon.github.io/spifa/reference/spifa.md)
starts every parameter at a generic default and lets the data speak. If
you have domain knowledge, setting `priors` closer to a plausible value
can speed up convergence substantially. Here we set a plausible prior
for every parameter group:

- item easiness roughly increasing across items,
- a couple of items with a strong sign-fixing discrimination prior (as
  in
  [`vignette("spifa")`](https://ErickChacon.github.io/spifa/articles/spifa.md)),
- a wealth effect that differs in sign across factors,
- moderate loadings,
- factor correlations nudged away from independence, and
- a range parameter around `150` as in
  [`vignette("spifa")`](https://ErickChacon.github.io/spifa/articles/spifa.md):

``` r

# easiness hyperparameters
c_mean <- seq(-0.7, 0.9, length.out = nitems)
c_sd <- rep(0.3, nitems)

# discrimination hyperparameters
A_mean <- matrix(0, nitems, nfactors)
A_mean[c(5, 6), 1] <- 1
A_mean[c(4, 8), 3] <- -1
A_sd <- matrix(1, nitems, nfactors)
A_sd[A_mean != 0] <- 0.45

# predictor effects hyperparameters
B_mean <- matrix(c(-0.5, 0.5, 0.2), nrow = 1)
B_sd <- matrix(0.5, nrow = 1, ncol = nfactors)

# correlation initial value
corr_initial <- diag(nfactors)

# range parameter hyperparameters
phi_mean <- 150
phi_sd <- 0.4
```

### Fit with constraints and priors

In the following code, we explicitely defined the `constraints`, and
`priors` for all possible parameters:

``` r

set.seed(1)
fit_informed <- spifa(
  items ~ wealth, data = ipixuna, nfactors = nfactors, ngp = 3,
  burnin = 200, niter = 200, thin = 2,
  constraints = list(discrimination = A, loading = diag(nfactors), sd = sd_fixed),
  priors = list(
    easiness = list(initial = c_mean, mean = c_mean, sd = c_sd),
    discrimination = list(initial = A_mean, mean = A_mean, sd = A_sd),
    effect = list(initial = B_mean, mean = B_mean, sd = B_sd),
    corr = list(initial = corr_initial, eta = 2),
    loading = list(initial = 0.6, mean = 0.6, sd = 0.3),
    range = list(initial = phi_mean, mean = phi_mean, sd = phi_sd)))
fit_informed
```

    #> Item factor analysis model: spifa_pred
    #> Formula: items ~ wealth
    #> Dimensions: 100 respondents, 10 items, 3 latent factors, 3 spatial processes
    #> MCMC: 1 chain, iter = 199, thin = 2, samples = 100
    #> 
    #> Item model parameters:
    #>           mean median   sd    q10      q90 ess_bulk rhat
    #> c[1]    -0.289 -0.286 0.19 -0.539 -0.06609     29.6 1.02
    #> c[2]    -0.299 -0.295 0.19 -0.535 -0.04725     17.3 1.05
    #> c[3]    -0.361 -0.335 0.21 -0.635 -0.06040     32.1 0.99
    #> c[4]    -0.220 -0.229 0.17 -0.410  0.00021     47.4 1.07
    #> c[5]    -0.259 -0.258 0.17 -0.489 -0.05588     51.8 1.00
    #> c[6]     0.014  0.012 0.12 -0.137  0.17382     61.7 1.02
    #> c[7]     0.165  0.160 0.14 -0.019  0.32253     67.8 0.99
    #> c[8]     0.171  0.152 0.17 -0.039  0.40831     56.8 1.03
    #> c[9]     0.508  0.492 0.16  0.283  0.72629     55.9 0.99
    #> c[10]    0.339  0.340 0.21  0.048  0.59850     33.4 1.00
    #> A[1,1]   0.496  0.462 0.25  0.184  0.85301      5.1 1.16
    #> A[2,1]   0.358  0.353 0.29 -0.039  0.75439      4.6 1.19
    #> A[3,1]   0.144  0.185 0.28 -0.217  0.47117      3.2 1.32
    #> A[5,1]   1.060  1.040 0.20  0.788  1.30934     20.6 1.07
    #> A[6,1]   1.017  1.018 0.20  0.792  1.26751     21.6 1.00
    #> A[7,1]   1.077  1.099 0.27  0.686  1.42900      4.7 1.16
    #> A[9,1]   0.966  0.947 0.31  0.607  1.35897     11.2 1.03
    #> A[10,1]  1.167  1.189 0.30  0.773  1.53663      2.9 1.29
    #> A[1,2]   1.142  1.133 0.38  0.769  1.67707     10.1 1.20
    #> A[2,2]   0.545  0.482 0.26  0.266  0.92465     17.3 1.00
    #> A[3,2]   0.352  0.328 0.25  0.035  0.67004      4.0 1.26
    #> A[9,2]   0.490  0.475 0.24  0.205  0.78129     12.2 1.08
    #> A[1,3]  -0.619 -0.626 0.26 -0.962 -0.25927      2.0 1.58
    #> A[2,3]  -0.829 -0.760 0.29 -1.191 -0.52113      4.6 1.19
    #> A[3,3]  -1.452 -1.412 0.45 -1.997 -0.88985      9.6 1.02
    #> A[4,3]  -0.717 -0.718 0.18 -0.930 -0.50310      4.9 1.16
    #> A[7,3]  -0.558 -0.543 0.20 -0.820 -0.30168     16.6 1.11
    #> A[8,3]  -0.859 -0.833 0.21 -1.163 -0.60641     11.6 1.04
    #> A[9,3]  -0.936 -0.852 0.34 -1.380 -0.55862      2.9 1.30
    #> A[10,3] -0.978 -0.947 0.30 -1.377 -0.62697      8.8 1.13
    #> 
    #> Factor model parameters:
    #>              mean  median     sd     q10    q90 ess_bulk rhat
    #> B[1,1]     -0.730  -0.738  0.080  -0.836  -0.63     66.9 1.00
    #> B[1,2]      0.690   0.687  0.080   0.599   0.78    101.6 0.99
    #> B[1,3]      0.333   0.331  0.092   0.219   0.44     64.4 1.02
    #> T[1,1]      0.752   0.762  0.058   0.665   0.82      7.2 1.16
    #> T[2,2]      0.712   0.732  0.063   0.617   0.77      7.5 1.13
    #> T[3,3]      0.915   0.926  0.062   0.831   0.98      3.6 1.20
    #> phi[1]    171.230 167.900 31.368 136.123 223.36      3.1 1.29
    #> phi[2]    151.843 150.887 24.972 120.307 187.44      3.7 1.27
    #> phi[3]    141.064 125.641 29.877 100.789 177.57      2.4 1.43
    #> Corr[2,1]   0.110   0.123  0.106   0.010   0.20      2.8 1.34
    #> Corr[3,1]   0.060   0.078  0.088  -0.038   0.19      2.7 1.39
    #> Corr[3,2]   0.098   0.110  0.082  -0.012   0.21      6.0 1.13
    #> 
    #> ess_bulk is the bulk effective sample size; rhat is the potential
    #> scale reduction factor on split chains (Rhat = 1 at convergence).
    #> Use summary() for the full set of statistics (incl. ess_tail).

## Continue sampling

Fits don’t always converge on the first try. If you take a look at
`rhat` in the previous section, you will notice that some values are
much higher that `1.05` for example. To continue sampling from the
current state you can use
[`update()`](https://rdrr.io/r/stats/update.html):

``` r

fit1 <- update(fit_informed, niter = 1000)
fit1
```

    #> Item factor analysis model: spifa_pred
    #> Formula: items ~ wealth
    #> Dimensions: 100 respondents, 10 items, 3 latent factors, 3 spatial processes
    #> MCMC: 1 chain, iter = 1000, thin = 1, samples = 1000
    #> 
    #> Item model parameters:
    #>           mean median   sd     q10    q90 ess_bulk rhat
    #> c[1]    -0.332 -0.331 0.19 -0.5868 -0.089    122.1  1.0
    #> c[2]    -0.404 -0.403 0.16 -0.6104 -0.184     92.1  1.0
    #> c[3]    -0.417 -0.419 0.17 -0.6301 -0.204     94.5  1.0
    #> c[4]    -0.306 -0.305 0.16 -0.5253 -0.091     77.0  1.0
    #> c[5]    -0.286 -0.281 0.18 -0.5141 -0.061    116.1  1.0
    #> c[6]    -0.011 -0.015 0.17 -0.2295  0.207    101.0  1.0
    #> c[7]     0.116  0.125 0.17 -0.1065  0.339     90.7  1.0
    #> c[8]     0.127  0.124 0.17 -0.0862  0.329    119.5  1.0
    #> c[9]     0.474  0.468 0.20  0.2336  0.738    125.0  1.0
    #> c[10]    0.244  0.239 0.20 -0.0086  0.505     68.3  1.0
    #> A[1,1]   0.530  0.520 0.28  0.1679  0.884      4.3  1.2
    #> A[2,1]   0.440  0.423 0.27  0.1075  0.784     46.9  1.0
    #> A[3,1]   0.189  0.195 0.29 -0.1809  0.545      5.6  1.2
    #> A[5,1]   1.357  1.328 0.28  1.0224  1.728     11.7  1.1
    #> A[6,1]   1.250  1.240 0.27  0.9049  1.625     43.6  1.0
    #> A[7,1]   1.074  1.029 0.33  0.7052  1.560     35.4  1.1
    #> A[9,1]   1.090  1.066 0.36  0.6769  1.533     56.0  1.0
    #> A[10,1]  1.227  1.182 0.37  0.8122  1.738     40.9  1.0
    #> A[1,2]   1.306  1.264 0.40  0.8601  1.755     14.8  1.1
    #> A[2,2]   0.765  0.745 0.28  0.4248  1.150     59.2  1.0
    #> A[3,2]   0.490  0.483 0.29  0.1131  0.858     21.5  1.1
    #> A[9,2]   0.706  0.697 0.34  0.2854  1.130     39.5  1.0
    #> A[1,3]  -0.524 -0.487 0.35 -0.9917 -0.107     16.5  1.1
    #> A[2,3]  -0.840 -0.790 0.31 -1.2653 -0.496     16.3  1.1
    #> A[3,3]  -1.171 -1.114 0.40 -1.7303 -0.677      7.6  1.1
    #> A[4,3]  -0.830 -0.812 0.23 -1.1244 -0.560     60.7  1.0
    #> A[7,3]  -0.607 -0.588 0.25 -0.9326 -0.301     83.4  1.0
    #> A[8,3]  -0.926 -0.907 0.24 -1.2475 -0.629     84.5  1.0
    #> A[9,3]  -1.034 -1.019 0.37 -1.5087 -0.579      9.0  1.1
    #> A[10,3] -1.065 -1.054 0.35 -1.5468 -0.627     12.5  1.1
    #> 
    #> Factor model parameters:
    #>               mean  median     sd    q10    q90 ess_bulk rhat
    #> B[1,1]     -0.6795  -0.678  0.088  -0.80  -0.57     71.4  1.0
    #> B[1,2]      0.6892   0.690  0.079   0.59   0.78    448.3  1.0
    #> B[1,3]      0.5137   0.508  0.092   0.40   0.63    356.2  1.0
    #> T[1,1]      0.7792   0.787  0.090   0.64   0.89      5.0  1.2
    #> T[2,2]      0.6811   0.679  0.082   0.57   0.78     14.5  1.2
    #> T[3,3]      0.8676   0.856  0.134   0.70   1.06     17.0  1.1
    #> phi[1]    125.2429 123.358 29.167  83.26 164.11     29.8  1.0
    #> phi[2]    163.7712 157.699 32.435 126.64 217.32      5.4  1.2
    #> phi[3]    151.8532 149.498 50.787  94.53 226.39      4.6  1.2
    #> Corr[2,1]   0.2094   0.224  0.105   0.08   0.34     20.0  1.2
    #> Corr[3,1]   0.0016   0.027  0.156  -0.21   0.16      4.2  1.2
    #> Corr[3,2]   0.0623   0.081  0.215  -0.24   0.32     22.4  1.0
    #> 
    #> ess_bulk is the bulk effective sample size; rhat is the potential
    #> scale reduction factor on split chains (Rhat = 1 at convergence).
    #> Use summary() for the full set of statistics (incl. ess_tail).

`rhat` has moved closer to 1 with the extra iterations. In practice, you
have to run many more iterations.

## Compare models with DIC

[`dic()`](https://ErickChacon.github.io/spifa/reference/dic.md) computes
the Deviance Information Criterion for a fitted model, useful for
comparing candidate models fitted to the same data. You should compare
them once convergence is ensured. For example, let’s compare with a
`spifa` with 2 latent factors and default `constraints` and `priors`:

``` r

fit2 <- spifa(items ~ wealth, data = ipixuna, nfactors = 2, niter = 1000, burnin = 1000)
fit2
```

    #> Item factor analysis model: spifa_pred
    #> Formula: items ~ wealth
    #> Dimensions: 100 respondents, 10 items, 2 latent factors, 2 spatial processes
    #> MCMC: 1 chain, iter = 1000, thin = 1, samples = 1000
    #> 
    #> Item model parameters:
    #>            mean median   sd    q10    q90 ess_bulk rhat
    #> c[1]    -0.3654 -0.363 0.13 -0.531 -0.201    443.3  1.0
    #> c[2]    -0.5366 -0.528 0.18 -0.775 -0.320     74.2  1.0
    #> c[3]    -0.6185 -0.606 0.25 -0.940 -0.321     69.9  1.0
    #> c[4]    -0.3479 -0.351 0.16 -0.565 -0.133    147.9  1.0
    #> c[5]    -0.3818 -0.368 0.30 -0.769 -0.042     33.0  1.0
    #> c[6]     0.0011  0.012 0.24 -0.280  0.265     84.3  1.0
    #> c[7]     0.2067  0.208 0.20 -0.051  0.466     65.3  1.0
    #> c[8]     0.1822  0.171 0.20 -0.072  0.435     91.9  1.0
    #> c[9]     0.6242  0.612 0.25  0.318  0.930     59.6  1.0
    #> c[10]    0.5143  0.483 0.33  0.107  0.970     41.5  1.1
    #> A[1,1]  -0.2321 -0.222 0.17 -0.453 -0.014     75.3  1.0
    #> A[2,1]  -0.4842 -0.462 0.28 -0.867 -0.142     28.9  1.0
    #> A[3,1]  -1.0417 -1.042 0.37 -1.488 -0.585     23.8  1.0
    #> A[4,1]  -0.8171 -0.768 0.36 -1.273 -0.409     29.2  1.0
    #> A[5,1]   0.9144  0.944 0.46  0.280  1.511     16.1  1.1
    #> A[6,1]   0.8598  0.829 0.41  0.352  1.393     22.7  1.0
    #> A[7,1]   0.0206  0.089 0.39 -0.558  0.486     16.9  1.0
    #> A[8,1]  -0.9578 -0.861 0.48 -1.519 -0.492      9.4  1.1
    #> A[9,1]  -0.4763 -0.427 0.44 -1.107  0.042     13.8  1.0
    #> A[10,1] -0.3895 -0.346 0.63 -1.250  0.374     12.9  1.1
    #> A[2,2]   0.5878  0.604 0.25  0.254  0.893     16.5  1.0
    #> A[3,2]   0.5240  0.562 0.44 -0.083  1.041     11.8  1.1
    #> A[4,2]   0.0576  0.116 0.37 -0.469  0.459      9.7  1.1
    #> A[5,2]   1.9341  1.873 0.66  1.081  2.922     13.5  1.1
    #> A[6,2]   1.4581  1.401 0.41  0.977  2.035     33.8  1.0
    #> A[7,2]   1.2592  1.240 0.33  0.850  1.681     43.0  1.0
    #> A[8,2]   0.3153  0.402 0.50 -0.270  0.833      9.7  1.1
    #> A[9,2]   1.2482  1.197 0.51  0.662  1.928     13.2  1.0
    #> A[10,2]  2.1616  2.118 0.56  1.436  2.945     18.7  1.0
    #> 
    #> Factor model parameters:
    #>             mean median     sd    q10    q90 ess_bulk rhat
    #> B[1,1]     0.095  0.090 0.1880 -0.144  0.348     11.4  1.1
    #> B[1,2]    -0.483 -0.483 0.1220 -0.646 -0.329    142.7  1.0
    #> T[1,1]     0.503  0.477 0.0996  0.383  0.669     29.1  1.0
    #> T[2,2]     0.559  0.550 0.1093  0.420  0.700      2.9  1.4
    #> phi[1]     0.051  0.051 0.0096  0.037  0.064     40.1  1.0
    #> phi[2]     0.048  0.047 0.0089  0.037  0.059     29.3  1.0
    #> Corr[2,1] -0.487 -0.494 0.2954 -0.858 -0.037      3.0  1.3
    #> 
    #> ess_bulk is the bulk effective sample size; rhat is the potential
    #> scale reduction factor on split chains (Rhat = 1 at convergence).
    #> Use summary() for the full set of statistics (incl. ess_tail).

Computing
[`dic()`](https://ErickChacon.github.io/spifa/reference/dic.md) for each
model:

``` r

# M1: 3 factors
dic(fit1, burnin = 100, thin = 2)
```

    #> # A tibble: 1 × 3
    #>   mean_deviance p_eff   dic
    #>           <dbl> <dbl> <dbl>
    #> 1          848.  146.  995.

``` r

# M2: 2 factors
dic(fit2, burnin = 100, thin = 2)
```

    #> # A tibble: 1 × 3
    #>   mean_deviance p_eff   dic
    #>           <dbl> <dbl> <dbl>
    #> 1          878.  146. 1024.

Assuming both models have converged, the 3-factor model has the lower
DIC, favouring the discrimination structure used above over an 2-factor
alternative.
