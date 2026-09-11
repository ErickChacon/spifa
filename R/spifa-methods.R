#' @title Print a Fitted spifa Model
#'
#' @description
#' Prints a fitted \code{spifa} object (the output of \code{\link{spifa}}):
#' model type, formula, data dimensions, MCMC settings, and a posterior
#' summary table (via \code{\link{summary.spifa}}) -- following the
#' convention of \pkg{rstan}/\pkg{rstanarm}/\pkg{brms}/\pkg{R2jags}, which
#' all show actual parameter estimates by default rather than just fit
#' metadata.
#'
#' @param x An object of class \code{spifa}, as returned by
#' \code{\link{spifa}}.
#' @param ... Further arguments passed to methods (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' data(ipixuna)
#' samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 0, niter = 20)
#' samples
#'
#' @export
print.spifa <- function (x, ...) {

  fit_args <- attr(x, "fit_args")
  predict_setup <- attr(x, "predict_setup")

  cat("Item factor analysis model: ", fit_args$model_type, "\n", sep = "")
  cat("Formula: ", deparse(predict_setup$formula), "\n", sep = "")
  cat("Dimensions: ", fit_args$nobs, " respondents, ", fit_args$nitems, " items, ",
      fit_args$nfactors, " latent factors, ", fit_args$ngp,
      " spatial processes\n", sep = "")

  if (is.list(x)) {
    cat(fit_args$niter, " iterations requested (thin = ", fit_args$thin,
        "); not executed (`execute = FALSE`)\n", sep = "")
  } else {
    cat("MCMC: ", "1 chain, iter = ", fit_args$niter, ", thin = ", fit_args$thin,
        ", samples = ", dim(x)[1], "\n\n", sep = "")
    print_summary(x, c("c", "A"), "Item model parameters:")
    print_summary(x, c("B", "T", "phi", "Corr"), "Factor model parameters:")

    cat("ess_bulk is the bulk effective sample size; rhat is the potential\n",
        "scale reduction factor on split chains (Rhat = 1 at convergence).\n",
        "Use summary() for the full set of statistics (incl. ess_tail).\n", sep = "")
  }
  invisible(x)
}

print_summary <- function (x, params, label) {
  params <- intersect(params, posterior::variables(x, with_indices = FALSE))
  sx <- summary(x, select = params) |>
    tibble::column_to_rownames("variable") |>
    subset(select = c("mean", "median", "sd", "q10", "q90", "ess_bulk", "rhat")) |>
    as.data.frame() |>
    format(digits = 2)
  cat(label, "\n", sep = "")
  print(sx)
  cat("\n")
}


#' @title Convert MCMC Samples to a Wide Tibble
#'
#' @description
#' Converts a fitted \code{spifa} object (the output of \code{\link{spifa}})
#' into a plain wide \code{\link[tibble]{tibble}}, with one column per
#' parameter and one row per (thinned, post-burnin) MCMC iteration. This is
#' a one-way export for use with \code{dplyr}/\code{ggplot2}/other tidyverse
#' tools directly — the package's own methods (\code{\link{summary.spifa}},
#' \code{\link{dic.spifa}}, \code{\link{predict.spifa}}) operate on the
#' fitted \code{spifa} object itself and do not need this conversion.
#'
#' @param x An object of class \code{spifa}, as returned by
#' \code{\link{spifa}}.
#' @param burnin Number of initial iterations to discard.
#' @param thin Thinning interval applied after discarding burn-in.
#' @param select Character vector of parameter blocks to keep (defaults to
#' all of them).
#' @param ... Further arguments passed to methods (currently unused; present
#' for consistency with the \code{\link[tibble]{as_tibble}} generic).
#'
#' @return A plain wide \code{\link[tibble]{tibble}}, with the
#' \code{"fit_args"} attribute carried over from \code{x}.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' nitems <- ncol(ipixuna$items)
#' nfactors <- 3
#'
#' # discrimination constraint: start with every item free to load on every
#' # factor, then restrict a few items per factor based on what each item is
#' # meant to measure (0 = no relationship, 1 = free parameter to estimate)
#' A <- matrix(1, nitems, nfactors)
#' A[c(4, 8), 1] <- 0
#' A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
#' A[c(5, 6), 3] <- 0
#' samples <- spifa(
#'   items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
#'   niter = 20, standardize = FALSE,
#'   constraints = list(discrimination = A))
#' samples_tib <- as_tibble(samples)
#' samples_tib
#' }
#'
#' @importFrom tibble as_tibble
#'
#' @export
as_tibble.spifa <- function (x, burnin = 0, thin = 1, select = NULL, ...) {

  fit_args <- attr(x, "fit_args")

  niter <- posterior::niterations(x)
  df <- x |>
    posterior::subset_draws(variable = select, iteration = (burnin + 1):niter) |>
    posterior::thin_draws(thin) |>
    posterior::as_draws_df() |>
    tibble::as_tibble() |>
    dplyr::select(-dplyr::starts_with("."))

  attr(df, "fit_args") <- fit_args
  return(df)
}

#' @title Convert a Fitted spifa Model to a spifa.list
#'
#' @description
#' Splits a fitted \code{spifa} object back into a named list of matrices,
#' one per parameter block (e.g. \code{c}, \code{a}, \code{theta}, ...).
#' Used internally by \code{\link{predict.spifa}} and \code{\link{dic.spifa}},
#' which need matrix-shaped samples to pass to their C++ counterparts.
#'
#' @param x An object of class \code{spifa}, as returned by
#' \code{\link{spifa}}.
#' @param ... Further arguments passed to methods (currently unused).
#'
#' @return An object of class \code{spifa.list}: a named list of matrices,
#' one per parameter block, with the \code{"fit_args"} attribute carried
#' over from \code{x}.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' nitems <- ncol(ipixuna$items)
#' nfactors <- 3
#'
#' # discrimination constraint: start with every item free to load on every
#' # factor, then restrict a few items per factor based on what each item is
#' # meant to measure (0 = no relationship, 1 = free parameter to estimate)
#' A <- matrix(1, nitems, nfactors)
#' A[c(4, 8), 1] <- 0
#' A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
#' A[c(5, 6), 3] <- 0
#' samples <- spifa(
#'   items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
#'   niter = 20, standardize = FALSE,
#'   constraints = list(discrimination = A))
#' samples_list <- as.list(samples)
#' names(samples_list)
#' }
#'
#' @export
as.list.spifa <- function (x, ...) {

  fit_args <- attr(x, "fit_args")
  blocks <- posterior::variables(x, with_indices = FALSE)
  niter <- dim(x)[1]

  extract_block <- function (block) {
    xb <- posterior::subset_draws(x, variable = block)
    matrix(xb[, 1, ], nrow = niter, ncol = dim(xb)[3],
           dimnames = list(NULL, dimnames(xb)[[3]]))
  }

  samples <- lapply(blocks, extract_block)
  names(samples) <- blocks
  class(samples) <- c("spifa.list", class(samples))
  attr(samples, "fit_args") <- fit_args
  attr(samples, "predict_setup") <- attr(x, "predict_setup")
  return(samples)
}

#' @title Posterior Summary of MCMC Samples
#'
#' @description
#' Computes posterior summary statistics for every parameter in a fitted
#' \code{spifa} object, via \code{\link[posterior]{summarise_draws}}: mean,
#' median, sd, mad, the 2.5\%/10\%/50\%/90\%/97.5\% quantiles, effective
#' sample size (bulk and tail), and \code{rhat}. \code{rhat} is computed via
#' split-chain R-hat (Vehtari et al. 2021), which splits each chain in half
#' and compares the halves -- so it remains a meaningful convergence
#' diagnostic even though \code{\link{spifa}} only ever fits a single chain.
#'
#' @param object A fitted \code{spifa} object, as returned by
#' \code{\link{spifa}}.
#' @param burnin Number of initial iterations to discard.
#' @param thin Thinning interval applied after discarding burn-in.
#' @param select Character vector of parameter blocks to summarise (defaults
#' to all of them). An error if any requested block does not exist in the
#' fitted model (e.g. \code{"T"} for a model with no spatial process).
#' @param ... Further arguments passed to methods (currently unused).
#'
#' @return A \code{\link[tibble]{tibble}} with one row per parameter; see
#' \code{\link[posterior]{summarise_draws}} for column details.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' nitems <- ncol(ipixuna$items)
#' nfactors <- 3
#'
#' # discrimination constraint: start with every item free to load on every
#' # factor, then restrict a few items per factor based on what each item is
#' # meant to measure (0 = no relationship, 1 = free parameter to estimate)
#' A <- matrix(1, nitems, nfactors)
#' A[c(4, 8), 1] <- 0
#' A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
#' A[c(5, 6), 3] <- 0
#' samples <- spifa(
#'   items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
#'   niter = 20, standardize = FALSE,
#'   constraints = list(discrimination = A))
#' summary(samples, select = "c")
#' }
#'
#' @export
summary.spifa <- function (object, burnin = 0, thin = 1, select = NULL, ...) {

  niter <- posterior::niterations(object)
  object <- object |>
    posterior::subset_draws(variable = select, iteration = (burnin+1):niter) |>
    posterior::thin_draws(thin)

  posterior::summarise_draws(object, "mean", "median", "sd", "mad",
    ~posterior::quantile2(.x, probs = c(0.025, 0.1, 0.9, 0.975)),
    "ess_bulk", "ess_tail", "rhat")
}
