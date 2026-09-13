# Continue Sampling from the Posterior of a spifa Model

Runs `niter` additional MCMC iterations, warm-started from the last
posterior draw of `object` (the easiness/discrimination/residual
correlation/predictor-effect/Gaussian process parameters), reusing the
data, priors, and constraints of the original fit. Returns only the new
draws, not concatenated with `object`'s – see Details.

## Usage

``` r
# S3 method for class 'spifa'
update(object, niter = 100, thin = 1, burnin = 0, ...)
```

## Arguments

- object:

  A fitted `spifa` object, as returned by
  [`spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md) with
  `execute = TRUE`.

- niter:

  Number of additional MCMC iterations to run and store.

- thin:

  Thinning interval for the newly stored MCMC samples.

- burnin:

  Number of initial iterations of this continuation to discard before
  storing (see
  [`spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md)'s
  `burnin`).

- ...:

  Further arguments (currently unused).

## Value

A new `spifa` object holding only the continuation draws (see
Description).

## Details

The adaptive Metropolis-Hastings proposal tuning (for the residual
correlation and, for `spifa`/`spifa_pred`, the Gaussian process
parameters) resumes from wherever `object`'s own run left it off, rather
than restarting from `object`'s original `adaptive` settings – so a
chain of [`update()`](https://rdrr.io/r/stats/update.html) calls keeps
refining its proposal instead of re-paying for adaptation each time.

[`update()`](https://rdrr.io/r/stats/update.html) always uses the same
`standardize` setting `object` itself was originally fit with (see
[`spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md)); it
isn't an argument here. `standardize`'s rescale is still computed
independently for each call, from that call's own posterior draws, so
`object` and the object returned here can end up on slightly different
absolute scales even though both represent the same continuous chain.
Combine them yourself (e.g.
[`rbind()`](https://rdrr.io/r/base/cbind.html) on their
[`as_draws_matrix`](https://mc-stan.org/posterior/reference/draws_matrix.html)
form) if you want one continuous chain; fit with `standardize = FALSE`
in the first place if you want every continuation on a genuinely
identical, unrescaled scale.

## Author

Erick A. Chacón-Montalván

## Examples

``` r
# \donttest{
data(ipixuna)
samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 0, niter = 1000)
more_samples <- update(samples, niter = 1000)
# }
```
