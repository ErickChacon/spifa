# Advanced features

## Introduction

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
    #>            mean median   sd     q10      q90 ess_bulk rhat
    #> c[1]    -0.2728 -0.272 0.20 -0.5391 -0.01118     34.2 1.02
    #> c[2]    -0.2918 -0.295 0.17 -0.5119 -0.04725     26.1 1.05
    #> c[3]    -0.3660 -0.357 0.20 -0.6340 -0.12828     37.2 1.00
    #> c[4]    -0.1988 -0.208 0.16 -0.3885  0.00021     30.0 1.06
    #> c[5]    -0.2691 -0.283 0.17 -0.4889 -0.06159     35.6 1.01
    #> c[6]     0.0089  0.014 0.14 -0.1669  0.18128     58.3 1.02
    #> c[7]     0.1594  0.160 0.16 -0.0282  0.34168     45.2 1.00
    #> c[8]     0.1802  0.159 0.17 -0.0370  0.41219     48.1 1.01
    #> c[9]     0.4917  0.477 0.15  0.2841  0.68948     58.9 0.99
    #> c[10]    0.3049  0.338 0.20  0.0481  0.56284     43.2 1.01
    #> A[1,1]   0.5144  0.526 0.28  0.1393  0.87232      8.1 1.10
    #> A[2,1]   0.5026  0.516 0.22  0.2411  0.79028     14.1 1.05
    #> A[3,1]   0.2330  0.226 0.19  0.0073  0.49672     12.3 1.10
    #> A[5,1]   1.1038  1.086 0.20  0.8728  1.35704     30.3 1.02
    #> A[6,1]   1.0077  0.999 0.21  0.7668  1.27589     26.9 1.00
    #> A[7,1]   1.1291  1.139 0.27  0.7777  1.47842     13.6 1.06
    #> A[9,1]   0.9366  0.933 0.22  0.6453  1.20621     17.7 1.02
    #> A[10,1]  1.3098  1.310 0.25  1.0038  1.56808     39.0 1.00
    #> A[1,2]   1.1351  1.140 0.42  0.7044  1.69366      9.3 1.21
    #> A[2,2]   0.6120  0.557 0.27  0.2842  1.02290     15.8 1.03
    #> A[3,2]   0.3943  0.392 0.27  0.0350  0.70584     11.7 1.09
    #> A[9,2]   0.4232  0.429 0.22  0.1256  0.73274     26.5 1.01
    #> A[1,3]  -0.5373 -0.550 0.24 -0.8625 -0.24338      4.4 1.24
    #> A[2,3]  -0.6796 -0.666 0.19 -0.9349 -0.48988     38.6 1.00
    #> A[3,3]  -1.3625 -1.301 0.50 -2.0296 -0.77968      9.7 1.02
    #> A[4,3]  -0.7173 -0.715 0.17 -0.9417 -0.51119     23.6 1.02
    #> A[7,3]  -0.5553 -0.546 0.18 -0.7772 -0.31375     19.0 1.09
    #> A[8,3]  -0.8525 -0.833 0.18 -1.1117 -0.63895     26.6 0.99
    #> A[9,3]  -0.7693 -0.767 0.18 -0.9619 -0.52137     20.6 1.06
    #> A[10,3] -0.8936 -0.867 0.29 -1.3075 -0.55830      3.8 1.23
    #> 
    #> Factor model parameters:
    #>              mean   median     sd     q10    q90 ess_bulk rhat
    #> B[1,1]     -0.728  -0.7332  0.074  -0.821  -0.63     77.7 1.00
    #> B[1,2]      0.693   0.6969  0.080   0.606   0.79    103.4 0.99
    #> B[1,3]      0.333   0.3314  0.097   0.210   0.46     82.2 1.00
    #> T[1,1]      0.748   0.7569  0.060   0.691   0.81      7.8 1.08
    #> T[2,2]      0.712   0.7250  0.062   0.656   0.77     10.5 1.05
    #> T[3,3]      0.925   0.9242  0.057   0.856   1.00     12.8 1.05
    #> phi[1]    164.243 156.7435 29.294 135.845 219.54     10.1 1.09
    #> phi[2]    160.257 163.7669 30.937 120.307 193.31      2.1 1.49
    #> phi[3]    152.577 155.7586 36.778 100.789 199.03      1.5 1.99
    #> Corr[2,1]   0.135   0.1544  0.138   0.010   0.29      1.9 1.55
    #> Corr[3,1]   0.033   0.0087  0.103  -0.088   0.19      1.7 1.70
    #> Corr[3,2]   0.066   0.0308  0.089  -0.025   0.19      2.1 1.54
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
    #>            mean  median   sd    q10    q90 ess_bulk rhat
    #> c[1]    -0.3841 -0.3763 0.19 -0.622 -0.152    149.0  1.0
    #> c[2]    -0.4209 -0.4187 0.16 -0.626 -0.209    170.2  1.0
    #> c[3]    -0.4393 -0.4323 0.17 -0.651 -0.234    154.2  1.0
    #> c[4]    -0.3062 -0.3055 0.15 -0.507 -0.111    195.4  1.0
    #> c[5]    -0.2592 -0.2568 0.19 -0.505 -0.015    107.7  1.0
    #> c[6]    -0.0045 -0.0047 0.18 -0.226  0.227    125.5  1.0
    #> c[7]     0.1317  0.1338 0.16 -0.073  0.335    166.1  1.0
    #> c[8]     0.1308  0.1261 0.16 -0.075  0.342    168.5  1.0
    #> c[9]     0.5059  0.5055 0.18  0.261  0.751    129.0  1.0
    #> c[10]    0.2660  0.2729 0.20  0.013  0.513    116.1  1.0
    #> A[1,1]   0.6418  0.6502 0.33  0.212  1.087     42.0  1.0
    #> A[2,1]   0.5447  0.5128 0.32  0.180  0.921     33.3  1.0
    #> A[3,1]   0.2975  0.2919 0.30 -0.093  0.692     60.4  1.0
    #> A[5,1]   1.4242  1.4131 0.28  1.076  1.790     80.9  1.0
    #> A[6,1]   1.2144  1.2189 0.24  0.900  1.511     83.2  1.0
    #> A[7,1]   1.0132  1.0025 0.24  0.708  1.322    100.4  1.0
    #> A[9,1]   1.3353  1.3078 0.36  0.882  1.836     35.9  1.0
    #> A[10,1]  1.4948  1.4684 0.37  1.047  1.994     25.2  1.0
    #> A[1,2]   1.3844  1.3619 0.40  0.903  1.897     47.2  1.0
    #> A[2,2]   0.7274  0.6794 0.34  0.370  1.090     73.0  1.0
    #> A[3,2]   0.5311  0.5289 0.29  0.169  0.906     76.9  1.0
    #> A[9,2]   0.8825  0.8392 0.35  0.441  1.397     45.3  1.0
    #> A[1,3]  -0.3572 -0.3469 0.33 -0.769  0.053      7.5  1.1
    #> A[2,3]  -0.6076 -0.6005 0.23 -0.918 -0.312     71.7  1.0
    #> A[3,3]  -1.1194 -1.0651 0.38 -1.667 -0.696     51.9  1.0
    #> A[4,3]  -0.8694 -0.8450 0.23 -1.190 -0.597     56.0  1.0
    #> A[7,3]  -0.5915 -0.5600 0.25 -0.939 -0.288     61.1  1.0
    #> A[8,3]  -0.9821 -0.9746 0.23 -1.256 -0.714     96.0  1.0
    #> A[9,3]  -0.9287 -0.9273 0.32 -1.336 -0.517     16.0  1.1
    #> A[10,3] -1.1977 -1.1589 0.39 -1.711 -0.743     21.0  1.0
    #> 
    #> Factor model parameters:
    #>              mean  median     sd    q10    q90 ess_bulk rhat
    #> B[1,1]     -0.717  -0.717  0.080 -0.823  -0.62    426.2  1.0
    #> B[1,2]      0.616   0.614  0.085  0.509   0.72    478.7  1.0
    #> B[1,3]      0.423   0.421  0.100  0.296   0.55    439.5  1.0
    #> T[1,1]      0.714   0.703  0.078  0.633   0.82     14.5  1.0
    #> T[2,2]      0.736   0.723  0.085  0.631   0.86     14.4  1.0
    #> T[3,3]      0.905   0.894  0.105  0.787   1.04      3.9  1.2
    #> phi[1]    124.369 118.818 29.093 89.196 164.94     33.6  1.0
    #> phi[2]    135.065 130.906 39.236 96.581 189.13     12.3  1.1
    #> phi[3]    125.415 121.503 33.953 82.171 176.66     30.0  1.1
    #> Corr[2,1]   0.109   0.115  0.137 -0.090   0.26     19.5  1.0
    #> Corr[3,1]   0.111   0.120  0.120 -0.037   0.27     14.0  1.0
    #> Corr[3,2]   0.089   0.095  0.119 -0.056   0.25     12.8  1.0
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
    #>            mean  median   sd    q10    q90 ess_bulk rhat
    #> c[1]    -0.4170 -0.4188 0.16 -0.630 -0.204    192.1  1.0
    #> c[2]    -0.5256 -0.5147 0.18 -0.757 -0.295    112.6  1.0
    #> c[3]    -0.6434 -0.6278 0.25 -0.970 -0.347     18.1  1.0
    #> c[4]    -0.3103 -0.3055 0.17 -0.529 -0.106    110.3  1.0
    #> c[5]    -0.4204 -0.3716 0.30 -0.856 -0.065     24.9  1.0
    #> c[6]     0.0271  0.0349 0.19 -0.219  0.267     96.7  1.0
    #> c[7]     0.2285  0.2191 0.20 -0.018  0.478     76.8  1.0
    #> c[8]     0.1752  0.1759 0.19 -0.066  0.415    102.2  1.0
    #> c[9]     0.7063  0.6894 0.29  0.353  1.104      6.0  1.1
    #> c[10]    0.5227  0.4902 0.33  0.130  0.960     16.6  1.1
    #> A[1,1]   0.5914  0.5814 0.20  0.352  0.831     88.8  1.0
    #> A[2,1]   0.6363  0.6274 0.32  0.240  1.010      4.6  1.2
    #> A[3,1]   1.3883  1.3365 0.44  0.835  2.007     32.5  1.0
    #> A[4,1]   0.7578  0.7252 0.28  0.431  1.139     71.6  1.0
    #> A[5,1]  -0.4476 -0.2722 0.66 -1.427  0.262      3.8  1.2
    #> A[6,1]  -0.3749 -0.3032 0.44 -1.021  0.124      3.3  1.2
    #> A[7,1]   0.3895  0.4888 0.44 -0.245  0.884      5.0  1.2
    #> A[8,1]   1.0542  1.0293 0.36  0.609  1.510     28.5  1.0
    #> A[9,1]   1.1131  1.0819 0.56  0.472  1.811     12.2  1.1
    #> A[10,1]  0.9868  1.0994 0.77 -0.102  1.912      8.1  1.2
    #> A[2,2]   0.4296  0.3976 0.29  0.092  0.830     11.0  1.1
    #> A[3,2]   0.1693  0.0936 0.47 -0.344  0.814     21.3  1.1
    #> A[4,2]  -0.1612 -0.1616 0.29 -0.545  0.209     23.0  1.0
    #> A[5,2]   2.3375  2.2633 0.73  1.403  3.388     12.9  1.1
    #> A[6,2]   1.3651  1.3758 0.33  0.931  1.799     61.8  1.0
    #> A[7,2]   1.2210  1.2057 0.36  0.770  1.737     47.1  1.0
    #> A[8,2]  -0.0027 -0.0099 0.33 -0.418  0.438     24.3  1.1
    #> A[9,2]   1.1227  1.0129 0.58  0.453  1.966      6.2  1.2
    #> A[10,2]  1.7876  1.7147 0.53  1.151  2.529      3.4  1.2
    #> 
    #> Factor model parameters:
    #>             mean median     sd    q10      q90 ess_bulk rhat
    #> B[1,1]    -0.205 -0.215 0.1608 -0.396  0.00019     27.4  1.0
    #> B[1,2]    -0.441 -0.435 0.1199 -0.594 -0.29127     87.8  1.0
    #> T[1,1]     0.520  0.508 0.1006  0.370  0.62799     30.8  1.0
    #> T[2,2]     0.521  0.513 0.0898  0.423  0.64680     28.1  1.0
    #> phi[1]     0.050  0.049 0.0091  0.039  0.06223     12.9  1.1
    #> phi[2]     0.049  0.047 0.0093  0.037  0.06030     27.6  1.1
    #> Corr[2,1]  0.486  0.485 0.2413  0.218  0.89649      7.4  1.2
    #> 
    #> ess_bulk is the bulk effective sample size; rhat is the potential
    #> scale reduction factor on split chains (Rhat = 1 at convergence).
    #> Use summary() for the full set of statistics (incl. ess_tail).

Use [`dic()`](https://ErickChacon.github.io/spifa/reference/dic.md) and
create a data.frame:

``` r

dplyr::bind_rows(
  "M1: 3 factors" = dic(fit1, burnin = 100, thin = 2),
  "M2: 2 factors" = dic(fit2, burnin = 100, thin = 2),
  .id = "model")
```

    #> # A tibble: 2 × 4
    #>   model         mean_deviance p_eff   dic
    #>   <chr>                 <dbl> <dbl> <dbl>
    #> 1 M1: 3 factors          834.  152.  986.
    #> 2 M2: 2 factors          869.  148. 1017.

Assuming both models have converged, the 3-factor model has the lower
DIC, favouring the discrimination structure used above over an arbitrary
2-factor alternative.
