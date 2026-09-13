# Default Plot of spifa Posterior Samples

The [`plot`](https://rdrr.io/r/base/plot.html) default for a fitted
`spifa` object: a combined trace + density view,
[`plot_trace`](https://ErickChacon.github.io/spifa/reference/plot_trace.md)
on the left and
[`plot_density`](https://ErickChacon.github.io/spifa/reference/plot_density.md)
on the right, one row per parameter and matched row for row, via
[`wrap_plots`](https://patchwork.data-imaginist.com/reference/wrap_plots.html)
– convergence and posterior shape at a glance, right after fitting.
Defaults to the easiness (`c`) and discrimination (`A`) parameters. For
a credible-interval view, or for any other parameter block, call
[`plot_trace`](https://ErickChacon.github.io/spifa/reference/plot_trace.md)/[`plot_density`](https://ErickChacon.github.io/spifa/reference/plot_density.md)/
[`plot_interval`](https://ErickChacon.github.io/spifa/reference/plot_interval.md)
directly instead.

## Usage

``` r
# S3 method for class 'spifa'
plot(x, select = c("c", "A"), nshow = 10, ...)
```

## Arguments

- x:

  A fitted `spifa` object, as returned by
  [`spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md).

- select:

  Parameters to plot, as in
  [`plot_trace`](https://ErickChacon.github.io/spifa/reference/plot_trace.md).
  Defaults to `c("c", "A")`.

- nshow:

  As in
  [`plot_trace`](https://ErickChacon.github.io/spifa/reference/plot_trace.md),
  but applied once up front so both panels show the same parameters:
  calling
  [`plot_trace`](https://ErickChacon.github.io/spifa/reference/plot_trace.md)/
  [`plot_density`](https://ErickChacon.github.io/spifa/reference/plot_density.md)
  separately with the same `select` can otherwise each pick a different
  random subsample.

- ...:

  Further arguments passed to both
  [`plot_trace`](https://ErickChacon.github.io/spifa/reference/plot_trace.md)
  and
  [`plot_density`](https://ErickChacon.github.io/spifa/reference/plot_density.md).

## Value

A `patchwork` object (a `ggplot`-like object).

## Author

Erick A. Chacón-Montalván

## Examples

``` r
# \donttest{
data(ipixuna)
samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 0, niter = 1000)

plot(samples)

plot(samples, select = "c")

# }
```
