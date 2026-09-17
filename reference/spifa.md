# Bayesian Spatial Item Factor Analysis

Fits exploratory, confirmatory, and spatial item factor analysis (IFA)
models for binary responses using full Bayesian inference. The model
represents each binary response as a thresholded continuous auxiliary
variable explained by `nfactors` latent abilities, optionally extended
with linear predictors and/or a multivariate Gaussian process to capture
spatial dependence in the latent factors (see
[`vignette( "spifa")`](https://ErickChacon.github.io/spifa/articles/spifa.md)
for a full worked example). Inference is done via Gibbs sampling with
adaptive Metropolis-Hastings updates for the spatial range and
correlation parameters.

## Usage

``` r
spifa(
  formula,
  data,
  nfactors,
  ngp = nfactors,
  niter = 100,
  thin = 1,
  burnin = 0,
  standardize = TRUE,
  constraints = list(discrimination = NULL, loading = NULL, sd = rep(1, nfactors)),
  priors = list(easiness = list(initial = NULL, mean = NULL, sd = NULL), discrimination =
    list(initial = NULL, mean = NULL, sd = NULL), effect = list(initial = NULL, mean =
    NULL, sd = NULL), corr = list(initial = NULL, eta = 1.5), loading = list(initial =
    NULL, mean = NULL, sd = NULL), range = list(initial = NULL, mean = NULL, sd = NULL)),
  adaptive = list(Sigma = NULL, Sigma_corr = NULL, Sigma_loading = NULL, Sigma_range =
    NULL, scale = 1, C = 0.7, alpha = 0.8, accep_prob = 0.234),
  execute = TRUE
)
```

## Arguments

- formula:

  A two-sided formula `items ~ predictors`. The left-hand side must be a
  single symbol naming a matrix-valued column of `data` holding the
  binary item responses (see Details). The right-hand side specifies
  predictors for the latent factors (e.g. `~ x1 + x2`); use `items ~ 1`
  for no predictors.

- data:

  A data frame containing the item-response matrix column named on the
  left-hand side of `formula` and, if used, the predictor columns named
  on its right-hand side. If `data` is an
  [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object, its
  geometry
  ([`st_geometry`](https://r-spatial.github.io/sf/reference/st_geometry.html))
  is used as the spatial coordinates and spatial Gaussian processes are
  added to the model.

- nfactors:

  Number of latent factors (dimensions of the ability construct).

- ngp:

  Number of independent Gaussian processes used to build the (possibly
  restricted) multivariate Gaussian process for the latent factors.
  Defaults to `nfactors` (one GP per factor). Only relevant when `data`
  is an `sf` object; set to `0` to fit a non-spatial model even when
  `data` has a geometry column (the geometry and `ngp` are otherwise
  ignored in that case).

- niter:

  Number of MCMC iterations to run (after `burnin`) and store.

- thin:

  Thinning interval for the stored MCMC samples.

- burnin:

  Number of initial MCMC iterations to discard. These iterations still
  run (and the adaptive Metropolis-Hastings proposals for
  `loading`/`range`/`corr` still adapt through them), but they are never
  stored, so `niter` counts only the iterations that end up in the
  returned samples. `0` by default (no iterations discarded during
  fitting – the previous behaviour). Prefer this over discarding a
  prefix of the samples afterwards (e.g. via `summary(..., burnin = )`):
  iterations dropped here were never stored, so they don't cost memory
  or thinning-index arithmetic, and the adaptive proposals get to keep
  converging across the burnin/niter boundary rather than someone
  accidentally analysing them as if they were post-adaptation draws.

- standardize:

  Logical; if `TRUE` (default), the stored posterior draws are rescaled
  after fitting so the latent factors have unit variance: `theta` is
  divided by its posterior SD per factor, and `discrimination`,
  `effect`, `loading`, and the residual SD (`constraints$sd`) are
  compensated by the same factor so the fitted response probabilities
  are unchanged. This only applies to models with a spatial Gaussian
  process and/or predictor effects on the latent factors
  (`cifa_pred`/`spifa`/`spifa_pred`; ignored for `eifa`/`cifa`, where
  the residual SD is the only source of `theta`'s variance and there's
  nothing to normalize against). Those model types have a multiplicative
  scale non-identifiability between `theta` and `discrimination`/the GP
  variance/the predictor effect – an equally good fit can shrink one and
  inflate the other – so leaving `TRUE` keeps draws on an interpretable,
  comparable scale. Set to `FALSE` to keep the raw, unscaled posterior,
  e.g. when checking recovery of known simulated parameters (see
  `dev/simulated/analyze-ipixuna.R`), where an extra rescale would make
  draws harder to compare directly against the true simulated values. No
  predictors are standardized by this argument – despite the name, it
  does not touch `formula`'s right-hand side at all.

- constraints:

  Named list of constraints associated to the factor model. Accepted
  names are \`discrimination\`, \`loading\`, and \`sd\`. The
  restrictions on the discrimination parameter should be placed in the
  element \`discrimination\` with same dimensions as the discrimination
  matrix (nitems x nfactors). A value of 0 indicates that the link
  between the item and the factor is disabled and 1 indicates that it
  remains active and the coefficient associated will be estimated. If
  left unspecified (\`NULL\`), a lower-triangular restriction is applied
  by default (the item on row \`i\` can only load on factors \`1:i\`) –
  this is a standard identifiability trick for exploratory factor
  analysis, not a fully unrestricted model, and gives EIFA (see
  Details). The restrictions for the multivariate Gaussian process
  loading matrix should be placed in the element \`loading\` with
  dimensions nfactors x ngp, such as a value of 0 indicates a link
  disconnected between the factor and the (independent) GP while 1
  indicates that it remains active. The restrictions with respect to the
  standard deviation of the latent factors' residual term should be
  placed in the element \`sd\`, which should be a vector (length
  nfactors) providing the fixed values for that standard deviation
  (paired with \`priors\$corr\`, together they parameterize the residual
  covariance). If the model includes predictors or a Gaussian process,
  it is recommended to be lower than 1.

- priors:

  Named list of initial values and prior hyperparameters, one element
  per parameter group: \`easiness\`, \`discrimination\`, \`effect\`
  (predictor effect on the latent factors), \`corr\` (correlation of the
  latent factors' residual term, paired with \`constraints\$sd\`),
  \`loading\` (multivariate Gaussian process loading matrix, paired with
  \`constraints\$loading\`), and \`range\` (multivariate Gaussian
  process scale parameters). Each element (except \`corr\`) accepts
  \`initial\`, \`mean\`, and \`sd\`; \`corr\` accepts \`initial\` and
  \`eta\` (the LKJ prior shape parameter). See the parameter glossary
  above for how these names map to the fitted model's sampled output.

  For \`eifa\` specifically, \`corr\` is never estimated (see Details) –
  \`priors\$corr\$initial\` (an \`nfactors x nfactors\` correlation
  matrix, 1s on the diagonal; identity if left unspecified) instead
  fixes the correlation between the latent factors at that value for the
  whole fit, and \`priors\$corr\$eta\` has no effect at all.

- adaptive:

  Named list of hyperparameters associated with the adaptive sampling.
  The adaptive sampling is done jointly for the \`correlation\`
  parameters, \`standard deviation of the gps\` and \`scale parameter of
  the gps\`. The matrix \`Sigma\` can be provided as the full covariance
  matrix of these parameters for the proposal distribution. Otherwise,
  part of this matrix can be provided by using the elements \`Sigma\`,
  \`Sigma_corr\`, \`Sigma_loading\` and \`Sigma_range\`. Additional
  elements are \`scale\`, \`C\`, \`alpha\` and \`accep_prob\` which are
  hyperparameters of the adaptive sampling proposed in Andrieu and
  Thomas (2008).

- execute:

  Logical value to run sampler or not. TRUE by default.

## Value

An object of class `spifa`: a
[`draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
of the MCMC samples (one variable per parameter, e.g. `c[1]`, `A[1,1]`,
`Theta[1,1]`, ...), with an attribute `"fit_args"` recording the data
and options used to fit the model (needed by
[`predict.spifa`](https://ErickChacon.github.io/spifa/reference/predict.spifa.md)
and [`dic`](https://ErickChacon.github.io/spifa/reference/dic.md)). See
[`summary.spifa`](https://ErickChacon.github.io/spifa/reference/summary.spifa.md),
[`plot_trace`](https://ErickChacon.github.io/spifa/reference/plot_trace.md),
[`plot_density`](https://ErickChacon.github.io/spifa/reference/plot_density.md),
and
[`plot_interval`](https://ErickChacon.github.io/spifa/reference/plot_interval.md)
for working with it.

## Details

The type of model fitted is determined automatically from `formula` and
the class of `data`: a one-sided right-hand side (`items ~ 1`) with
`constraints$discrimination` left unspecified gives exploratory IFA
(EIFA) – a default lower-triangular restriction on `discrimination` is
applied automatically for identifiability in this case, not something
the user chooses (see `constraints` below). EIFA also never estimates
the correlation between the latent factors: it stays fixed at
`priors$corr$initial` (identity by default) for the whole fit – see
`priors` below. Supplying your own (typically theory-driven) restricted
`constraints$discrimination` instead gives confirmatory IFA (CIFA);
adding predictors to the right-hand side (e.g. `items ~ x1`) gives CIFA
with predictors; `data` being an
[`sf`](https://r-spatial.github.io/sf/reference/sf.html) object adds a
spatial Gaussian process on the latent factors, giving spatial IFA
(SPIFA), with or without predictors.

The left-hand side of `formula` must be a single symbol naming a
matrix-valued column of `data` (`nobs x nitems`, one row per respondent,
one column per binary item) — the same mechanism base R uses for
multivariate [`lm`](https://rdrr.io/r/stats/lm.html). Build it with
[`I`](https://rdrr.io/r/base/AsIs.html) (or ordinary `$<-` assignment)
so it survives as a matrix column rather than being flattened into
separate columns, e.g.:


    items <- as.matrix(dplyr::select(data, `Item 1`:`Item 10`))
    data$items <- items
    spifa(items ~ x1, data = data, nfactors = 2)

Missing values are handled differently depending on where they occur. A
missing item response (`NA` in the response matrix) does not drop that
respondent: it is treated as an unobserved auxiliary variable and
sampled natively along with everything else. A missing predictor value
(right-hand side of `formula`), by contrast, drops that respondent
entirely, the same way [`lm`](https://rdrr.io/r/stats/lm.html) and
friends do.

**Parameter glossary.** `priors`/`constraints` use descriptive names;
the fitted model's sampled output (as seen via
[`summary.spifa`](https://ErickChacon.github.io/spifa/reference/summary.spifa.md),
[`plot_trace`](https://ErickChacon.github.io/spifa/reference/plot_trace.md),
and friends) instead uses short internal names matching the underlying
model notation. The two are deliberately different vocabularies (what
you configure vs. what the sampler produced) – this table maps between
them:

|  |  |  |
|----|----|----|
| `priors`/`constraints` name | output group name | meaning |
| `easiness` | `c` | item easiness (intercept) |
| `discrimination` | `A` | item-factor discrimination (loading) matrix |
| `effect` | `B` | predictor effect on the latent factors |
| `corr` | `Corr`, `Chol` | residual correlation matrix, and its Cholesky factor |
| `sd` | *(fixed, not sampled)* | residual standard deviation |
| `loading` | `T` | multivariate Gaussian process loading matrix |
| `range` | `phi` | multivariate Gaussian process spatial range |
| *(not user-set)* | `Theta` | latent abilities |
| *(not user-set)* | `Z` | augmented latent response |

## References

Chacon-Montalvan, E. A., Parry, L., Giorgi, E., Torres, P., Orellana, J.
D. Y., Moraga, P., and Taylor, B. M. (2025). Mapping food insecurity in
the Brazilian Amazon using a spatial item factor analysis model. The
Annals of Applied Statistics, 19(4), 3438-3463.
[doi:10.1214/25-AOAS2072](https://doi.org/10.1214/25-AOAS2072)

## Author

Erick A. Chacón-Montalván

## Examples

``` r
# \donttest{
data(ipixuna)
nitems <- ncol(ipixuna$items)
nfactors <- 3

# EIFA: default discrimination constraints, non-spatial
samples_eifa <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
  niter = 20)

# discrimination constraint for cifa/spifa
A <- matrix(1, nitems, nfactors)
A[c(4, 8), 1] <- 0
A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
A[c(5, 6), 3] <- 0

# CIFA with predictors: custom discrimination constraints, non-spatial
samples_cifa <- spifa(items ~ poly(wealth, 2), data = ipixuna, nfactors = nfactors,
  ngp = 0, niter = 20, burnin = 5, constraints = list(discrimination = A))

# SPIFA: custom discrimination constraints, default spatial components
samples_spifa <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors,
  niter = 20, burnin = 5, thin = 2, constraints = list(discrimination = A))

# SPIFA with predictors and a custom loading matrix: factors 1 and 2
# share one Gaussian process, factor 3 gets its own
ngp <- 2
T <- matrix(c(1, 1, 0, 0, 0, 1), nfactors, ngp)
samples_spifa_pred <- spifa(items ~ wealth, data = ipixuna, nfactors = nfactors,
  ngp = 2, niter = 20, constraints = list(discrimination = A, loading = T))
# }
```
