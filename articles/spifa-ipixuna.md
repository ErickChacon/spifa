# Spatial Item Factor Analysis

## Introduction

In this vignette, we show how to use the **spifa** package to fit a
spatial multidimensional 2 parameter logistic model, as used in item
response theory. Let $`Y_{ij}`$ be the binary response for item $`j`$ in
individual $`i`$. The model can be defined using an auxiliary variable
$`Z_{ij}`$ such that

\$\$\begin{align} {Y}\_{ij} & = \left\lbrace \begin{array}\[2\]{cc} 1 &
\text{if} ~ {Z}\_{ij} \> 0\\ 0 & \text{otherwise} \end{array} \right.\\
{Z}\_{ij} & = c_j + a_j^\intercal\theta_i + \epsilon\_{ij}, \~~
\epsilon\_{ij} \sim {N}(0, 1) \end{align}\$\$

where the latent abilities $`\theta_i`$ are, in the spatial case,
explained by a multivariate Gaussian process rather than assumed
independent across individuals.

## Load required packages and data

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

library(ggplot2)
library(spifa)
library(sf)
```

    #> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE

``` r

data(ipixuna)
```

`ipixuna` has one row per (simulated) household, with a ten-item binary
`items` response matrix, a covariate `wealth`, and a spatial coordinate
`geometry`:

``` r

ipixuna |>
  dplyr::select(-geometry, -items) |>
  head(10) |>
  knitr::kable(digits = 2)
```

|  id | wealth | geometry                    |
|----:|-------:|:----------------------------|
|   1 |  -0.57 | POINT (-71.69689 -7.052244) |
|   2 |  -0.92 | POINT (-71.68695 -7.039144) |
|   3 |   1.23 | POINT (-71.68732 -7.047609) |
|   4 |   0.31 | POINT (-71.69005 -7.047915) |
|   5 |  -0.06 | POINT (-71.68464 -7.053396) |
|   6 |   2.05 | POINT (-71.69039 -7.042424) |
|   7 |   2.31 | POINT (-71.69486 -7.050021) |
|   8 |   0.17 | POINT (-71.6863 -7.047552)  |
|   9 |  -0.22 | POINT (-71.69208 -7.053133) |
|  10 |  -0.92 | POINT (-71.68778 -7.043805) |

The true generating parameters used to simulate the data (including the
discrimination structure) are stored as an attribute, which is
convenient both for setting sensible restrictions and for checking
recovery below:

``` r

parameters <- attr(ipixuna, "parameters")
str(parameters)
```

    #> List of 6
    #>  $ easiness      : num [1:10] -0.68 -0.51 -0.34 -0.17 0 0.17 0.34 0.51 0.68 0.85
    #>  $ discrimination: num [1:10, 1:3] 0.771 0.651 0.596 0 1.285 ...
    #>  $ abilities     : num [1:100, 1:3] 1.414 0.826 -1.613 -0.909 0.844 ...
    #>  $ effect        : num [1:3] -0.7 0.7 0.3
    #>  $ resid_params  :List of 2
    #>   ..$ sd  : num [1:3] 0.447 0.548 0.447
    #>   ..$ corr: num [1:3] -0.7 0.1 -0.1
    #>  $ mgp_params    :List of 2
    #>   ..$ sd : num [1:3] 0.632 0.548 0.707
    #>   ..$ phi: num [1:3] 100 150 250

## Fit a SPIFA model for the Ipixuna data

We restrict the discrimination matrix to the structure used to simulate
the data, and let the three latent factors follow independent Gaussian
processes (`W = diag(3)`) over the household coordinates. `niter` is
kept small here so the vignette builds quickly; in practice use a few
thousand iterations and check convergence (see
[`?plot_trace`](https://ErickChacon.github.io/spifa/reference/plot_trace.md)).

``` r

L_a <- (parameters$discrimination != 0) * 1
nfactors <- ncol(parameters$discrimination)

system.time(
  samples <- spifa(
    items ~ wealth, data = ipixuna, nfactors = nfactors,
    niter = 1000, thin = 1, standardize = FALSE,
    constraints = list(discrimination = L_a, loading = diag(nfactors),
      sd = parameters$resid_params$sd),
    priors = list(
      loading = list(initial = 0.6, mean = 0.6, sd = 0.4),
      range = list(initial = 200, mean = 200, sd = 0.4)))
  )
```

    #>    user  system elapsed 
    #>   3.755   3.288   1.797

``` r

samples
```

    #> Item factor analysis model: spifa_pred
    #> Formula: items ~ wealth
    #> Dimensions: 100 respondents, 10 items, 3 latent factors, 3 spatial processes
    #> MCMC: 1 chain, iter = 1000, thin = 1, samples = 1000
    #> 
    #> Item model parameters:
    #>           mean median   sd     q10    q90 ess_bulk rhat
    #> c[1]    -0.681 -0.652 0.40 -1.1601 -0.226     17.7  1.1
    #> c[2]    -0.617 -0.601 0.24 -0.9260 -0.325     27.8  1.0
    #> c[3]    -0.699 -0.687 0.31 -1.0356 -0.344      8.6  1.1
    #> c[4]    -0.398 -0.377 0.26 -0.7417 -0.098     48.8  1.0
    #> c[5]    -0.280 -0.245 0.35 -0.7170  0.124     19.8  1.1
    #> c[6]     0.082  0.088 0.27 -0.2447  0.420     27.7  1.0
    #> c[7]     0.242  0.217 0.26 -0.0696  0.605     10.3  1.1
    #> c[8]     0.139  0.119 0.30 -0.2140  0.529     34.6  1.0
    #> c[9]     0.705  0.692 0.33  0.2955  1.139     25.0  1.0
    #> c[10]    0.409  0.396 0.32  0.0061  0.845     20.0  1.1
    #> A[1,1]   0.504  0.502 0.31  0.1125  0.902     58.4  1.0
    #> A[2,1]   0.453  0.442 0.28  0.1050  0.811     86.5  1.0
    #> A[3,1]   0.271  0.273 0.28 -0.0848  0.637     80.5  1.0
    #> A[4,1]   0.000  0.000 0.00  0.0000  0.000       NA   NA
    #> A[5,1]   1.927  1.864 0.53  1.3294  2.667     38.0  1.0
    #> A[6,1]   1.491  1.446 0.44  0.9731  2.086     51.0  1.0
    #> A[7,1]   1.371  1.349 0.40  0.8660  1.933     15.2  1.1
    #> A[8,1]   0.000  0.000 0.00  0.0000  0.000       NA   NA
    #> A[9,1]   1.579  1.581 0.39  1.0982  2.074     43.6  1.1
    #> A[10,1]  1.748  1.754 0.50  1.1086  2.417     31.2  1.0
    #> A[1,2]   1.750  1.770 0.48  1.1254  2.322     38.2  1.1
    #> A[2,2]   0.884  0.860 0.31  0.5091  1.321     43.1  1.0
    #> A[3,2]   0.818  0.805 0.36  0.3655  1.295     44.1  1.1
    #> A[4,2]   0.000  0.000 0.00  0.0000  0.000       NA   NA
    #> A[5,2]   0.000  0.000 0.00  0.0000  0.000       NA   NA
    #> A[6,2]   0.000  0.000 0.00  0.0000  0.000       NA   NA
    #> A[7,2]   0.000  0.000 0.00  0.0000  0.000       NA   NA
    #> A[8,2]   0.000  0.000 0.00  0.0000  0.000       NA   NA
    #> A[9,2]   1.253  1.223 0.44  0.7055  1.844      5.6  1.2
    #> A[10,2]  0.000  0.000 0.00  0.0000  0.000       NA   NA
    #> A[1,3]  -0.109 -0.117 0.47 -0.7023  0.514      4.2  1.2
    #> A[2,3]   0.733  0.714 0.35  0.3315  1.151      8.2  1.1
    #> A[3,3]   1.325  1.303 0.31  0.9380  1.723     98.8  1.0
    #> A[4,3]   1.231  1.200 0.39  0.7591  1.763     48.1  1.0
    #> A[5,3]   0.000  0.000 0.00  0.0000  0.000       NA   NA
    #> A[6,3]   0.000  0.000 0.00  0.0000  0.000       NA   NA
    #> A[7,3]   0.901  0.893 0.33  0.4715  1.322     94.2  1.0
    #> A[8,3]   1.532  1.455 0.50  0.9408  2.256     27.6  1.0
    #> A[9,3]   1.100  1.105 0.52  0.4109  1.744     36.6  1.0
    #> A[10,3]  1.543  1.474 0.49  0.9720  2.199     52.3  1.0
    #> 
    #> Factor model parameters:
    #>              mean  median     sd     q10     q90 ess_bulk rhat
    #> B[1,1]     -0.318  -0.317  0.089  -0.434  -0.208    112.7  1.0
    #> B[1,2]      0.062   0.065  0.100  -0.068   0.190    124.5  1.0
    #> B[1,3]     -0.163  -0.165  0.080  -0.264  -0.062    330.1  1.0
    #> T[1,1]      0.766   0.748  0.117   0.619   0.930      5.6  1.1
    #> T[2,2]      0.629   0.620  0.102   0.513   0.775     12.1  1.1
    #> T[3,3]      0.561   0.558  0.038   0.518   0.608     32.2  1.0
    #> phi[1]    194.277 203.266 39.121 124.353 238.774      4.8  1.3
    #> phi[2]    197.863 188.553 47.421 140.099 256.462     47.6  1.0
    #> phi[3]    215.582 205.249 77.435 132.314 325.377     21.1  1.0
    #> Corr[2,1]  -0.091  -0.068  0.221  -0.402   0.184     19.8  1.1
    #> Corr[3,1]   0.076   0.078  0.101  -0.085   0.203     31.0  1.0
    #> Corr[3,2]   0.111   0.057  0.180  -0.055   0.409      6.6  1.2
    #> 
    #> ess_bulk is the bulk effective sample size; rhat is the potential
    #> scale reduction factor on split chains (Rhat = 1 at convergence).
    #> Use summary() for the full set of statistics (incl. ess_tail).

``` r

attr(samples, "fit_args")$model_type
```

    #> [1] "spifa_pred"

## Summarise samples

``` r

summary(samples, burnin = 100, select = "c")
```

    #> # A tibble: 10 × 12
    #>    variable    mean  median    sd   mad    q2.5      q10     q90   q97.5 ess_bulk ess_tail  rhat
    #>    <chr>      <dbl>   <dbl> <dbl> <dbl>   <dbl>    <dbl>   <dbl>   <dbl>    <dbl>    <dbl> <dbl>
    #>  1 c[1]     -0.720  -0.668  0.381 0.311 -1.70   -1.19    -0.294  -0.0830    21.8      33.5  1.09
    #>  2 c[2]     -0.621  -0.600  0.237 0.217 -1.20   -0.932   -0.336  -0.224     23.1      27.5  1.05
    #>  3 c[3]     -0.694  -0.681  0.323 0.276 -1.66   -1.04    -0.319  -0.150      8.26     25.8  1.14
    #>  4 c[4]     -0.384  -0.363  0.265 0.239 -0.988  -0.718   -0.0907  0.114     43.8      62.3  1.03
    #>  5 c[5]     -0.305  -0.272  0.352 0.336 -1.00   -0.742    0.103   0.277      6.19     14.4  1.15
    #>  6 c[6]      0.0699  0.0740 0.271 0.251 -0.522  -0.257    0.420   0.589     11.7      64.9  1.08
    #>  7 c[7]      0.244   0.218  0.268 0.253 -0.249  -0.0747   0.617   0.797      4.53     39.6  1.17
    #>  8 c[8]      0.166   0.141  0.286 0.291 -0.340  -0.189    0.552   0.768     33.7      83.4  1.01
    #>  9 c[9]      0.699   0.690  0.332 0.319  0.0458  0.280    1.13    1.36      19.4      80.7  1.04
    #> 10 c[10]     0.418   0.411  0.326 0.325 -0.206   0.00342  0.848   1.07       5.98     60.4  1.15

## Visualize results

Traceplots of the easiness (`c`) and discrimination (`a`) parameters:

``` r

plot_trace(samples, select = "c")
```

![](spifa-ipixuna_files/figure-html/unnamed-chunk-6-1.png)

``` r

plot_trace(samples, select = "A")
```

![](spifa-ipixuna_files/figure-html/unnamed-chunk-6-2.png)

Traceplots of the residual correlation and of the Gaussian process
hyperparameters (loading matrix `T` and range `phi`):

``` r

plot_trace(samples, select = "Corr")
```

![](spifa-ipixuna_files/figure-html/unnamed-chunk-7-1.png)

``` r

plot_trace(samples, select = "T")
```

![](spifa-ipixuna_files/figure-html/unnamed-chunk-7-2.png)

``` r

plot_trace(samples, select = "phi")
```

![](spifa-ipixuna_files/figure-html/unnamed-chunk-7-3.png)

Credible intervals for the discrimination parameters, against the true
simulated values:

``` r

plot_interval(samples, select = "A", burnin = 100,
              reference = parameters$discrimination)
```

![](spifa-ipixuna_files/figure-html/unnamed-chunk-8-1.png)
