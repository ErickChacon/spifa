#' @title Print spifa Posterior Samples
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
  if (length(params) == 0) return(invisible())
  sx <- summary(x, select = params) |>
    tibble::column_to_rownames("variable") |>
    subset(select = c("mean", "median", "sd", "q10", "q90", "ess_bulk", "rhat")) |>
    as.data.frame() |>
    format(digits = 2)
  cat(label, "\n", sep = "")
  print(sx)
  cat("\n")
}


#' @title Summarise spifa Posterior Samples
#'
#' @description
#' Computes posterior summary statistics for every parameter in a fitted
#' \code{spifa} object, via \code{\link[posterior]{summarise_draws}}: mean,
#' median, sd, mad, the 2.5\%/10\%/50\%/90\%/97.5\% quantiles, effective
#' sample size (bulk and tail), and \code{rhat}. \code{rhat} is computed via
#' split-chain R-hat (Vehtari et al. 2021), which splits each chain in half
#' and compares the halves -- so it remains a meaningful convergence
#' diagnostic even though \code{\link{spifa}} only ever fits a single chain.
#' Discrimination parameters (\code{A}) structurally restricted to zero
#' (via \code{constraints$discrimination}) are excluded, since they are
#' fixed by construction rather than estimated.
#'
#' @param object A fitted \code{spifa} object, as returned by
#' \code{\link{spifa}}.
#' @param burnin Number of initial iterations to discard.
#' @param thin Thinning interval applied after discarding burn-in.
#' @param select Character vector of parameter groups to summarise (defaults
#' to all of them). An error if any requested group does not exist in the
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
    drop_restricted() |>
    posterior::subset_draws(variable = select, iteration = (burnin+1):niter) |>
    posterior::thin_draws(thin)

  posterior::summarise_draws(object, "mean", "median", "sd", "mad",
    ~posterior::quantile2(.x, probs = c(0.025, 0.1, 0.9, 0.975)),
    "ess_bulk", "ess_tail", "rhat")
}
