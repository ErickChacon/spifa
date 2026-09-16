#' @title Deviance Information Criterion
#'
#' @description
#' Generic function for the Deviance Information Criterion (DIC), a
#' Bayesian measure of model fit that penalises complexity. See
#' \code{\link{dic.spifa}} for the \code{spifa} method.
#'
#' @param x A fitted model object.
#' @param ... Further arguments passed to methods.
#'
#' @export
dic <- function (x, ...) {
  UseMethod("dic", x)
}

#' @title Deviance Information Criterion for a spifa Model
#'
#' @description
#' Computes the Deviance Information Criterion (DIC) for a fitted
#' \code{spifa} model, useful for comparing candidate models (e.g. different
#' numbers of factors or different restriction structures) fitted to the
#' same data.
#'
#' @param x A fitted \code{spifa} object, as returned by \code{\link{spifa}}.
#' @param burnin Number of initial iterations to discard.
#' @param thin Thinning interval applied after discarding burn-in.
#' @param ... Further arguments passed to methods (currently unused).
#'
#' @return A one-row \code{\link[tibble]{tibble}} with columns
#' \code{mean_deviance} (posterior mean of the deviance), \code{p_eff}
#' (effective number of parameters), and \code{dic}.
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
#' dic(samples)
#' }
#'
#' @export
dic.spifa <- function (x, burnin = 0, thin = 1, ...) {

  fit_args <- attr(x, "fit_args")
  niter <- posterior::niterations(x)
  x <- posterior::subset_draws(x, iteration = (burnin + 1):niter) |>
      posterior::thin_draws(thin)

  output <- dic_cpp(
    y = fit_args$response,
    c = as_matrix(x, "c"),
    a = as_matrix(x, "A"),
    theta = as_matrix(x, "Theta"),
    n = fit_args$nobs,
    q = fit_args$nitems,
    m = fit_args$nfactors,
    L = fit_args$constrain_L
  )

  tibble::as_tibble(output)
}

as_matrix <- function (x, param) {
  posterior::as_draws_matrix(posterior::subset_draws(x, variable = param))
}

as_matrix_opt <- function (x, param) {
  if (param %in% posterior::variables(x, with_indices = FALSE)) {
    as_matrix(x, param)
  } else {
    nsamples <- posterior::ndraws(x)
    matrix(nrow = nsamples, ncol = 0)
  }
}

#' @title Predict the Latent Factors of a spifa Model
#'
#' @description
#' Predicts the latent factors of (spatial) item factor analysis
#' for new subjects/locations and/or for new predictor values, using
#' the posterior samples from a fitted \code{\link{spifa}} model.
#'
#' @details
#' If the fitted model has no spatial or predictor structure (\code{eifa} or
#' \code{cifa}), or if \code{newdata} is not supplied for a model that has
#' one, there is nothing to predict beyond the latent abilities' own
#' posterior samples already available from the fit, so those are returned
#' directly (subject to \code{burnin}/\code{thin}) instead of calling the
#' \code{C++} sampler. Otherwise, prediction for the new locations and/or
#' predictor values is delegated to the \code{C++} sampler.
#'
#' If the fitted model has predictors (\code{cifa_pred}/\code{spifa_pred}),
#' \code{newdata} must include those predictor columns, the same as
#' \code{\link[stats]{predict.lm}} and similar methods require -- this is
#' an error otherwise. There's no synthesized reference-level fallback for
#' missing predictors (e.g. an \code{sf}/\code{sfc} object holding only new
#' locations, with no predictor columns at all): for a factor predictor
#' under the no-intercept encoding \code{\link{spifa}} uses, an
#' automatically-filled all-zero row wouldn't correspond to any real
#' category, so any reference profile -- including all zeros for numeric
#' predictors -- must be supplied explicitly in \code{newdata}, matching
#' the original data's format.
#'
#' @param object A fitted \code{spifa} object, as returned by
#' \code{\link{spifa}}.
#' @param newdata New data to predict at, mirroring \code{\link{spifa}}'s
#' own \code{data} argument: an \code{\link[sf]{sf}}/\code{\link[sf]{sfc}}
#' object (for spatial, i.e. \code{spifa}/\code{spifa_pred}, models --
#' its geometry gives the new locations, and its CRS is used directly, so
#' it need not match the training data's CRS) or a plain data frame (for
#' \code{cifa_pred}), with columns matching the predictors used on the
#' right-hand side of \code{formula} when the model was fitted. Its design
#' matrix is built the same way \code{\link{spifa}} built the training one,
#' using the same terms and factor levels. For a spatial model, its geometry
#' may also be \code{POLYGON}/\code{MULTIPOLYGON} (e.g. a prediction grid): the
#' centroid of each cell is used for the spatial kernel, but the original
#' polygons are kept (see \code{Value}) so \code{\link{plot_predict}} can draw
#' them as a filled map instead of points.
#' @param burnin Number of initial (post-fitting) iterations to discard
#' before using the posterior samples for prediction.
#' @param thin Thinning interval applied to the posterior samples used for
#' prediction.
#' @param joint Logical; for spatial models (\code{spifa}/\code{spifa_pred}),
#' whether the posterior predictive draws should respect the full predictive
#' covariance across new locations and factors (\code{TRUE}), or be drawn
#' marginally/independently per location-factor combination (\code{FALSE},
#' the default, cheaper). Each draw still propagates posterior parameter
#' uncertainty (one draw per retained MCMC iteration) either way -- this
#' only controls whether, within a single draw, the values across new
#' locations/factors are jointly correlated as the model implies. Marginal
#' draws are fine for per-location summaries (e.g. means, credible
#' intervals computed independently per column); set \code{TRUE} when the
#' samples themselves will be used as input to another model or computation
#' that depends on their joint structure (e.g. a spatial contrast or
#' aggregate across new locations). Ignored for \code{cifa_pred}, whose
#' draws are already jointly correct across factors.
#' @param ... Further arguments (currently unused).
#'
#' @return A \code{\link[posterior]{draws_array}} of posterior predictive
#' samples of the latent abilities (\code{theta}) for the requested new
#' locations and/or predictor values (or, if no prediction was requested,
#' for the originally observed subjects). For a spatial model, the locations
#' used (\code{newdata}, or the training locations if \code{newdata} was
#' omitted) are attached as the \code{"newdata"} attribute, so
#' \code{\link{plot_predict}} can map the result without needing it supplied
#' again.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' \donttest{
#' library(sf)
#' data(ipixuna)
#'
#' nitems <- ncol(ipixuna$items)
#' nfactors <- 3
#' A <- matrix(1, nitems, nfactors)
#' A[c(4, 8), 1] <- 0
#' A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
#' A[c(5, 6), 3] <- 0
#'
#' # Spifa model
#' samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, niter = 5,
#'   constraints = list(discrimination = A))
#' # latent abilities for observed locations
#' predict(samples)
#' # latent abilities for new locations
#' newdata <- st_make_grid(ipixuna, n = c(3, 2), what = "centers")
#' predict(samples, newdata = newdata)
#'
#' # Spifa model with predictors
#' samples_pred <- spifa(items ~ wealth, data = ipixuna, nfactors = nfactors, niter = 5,
#'   constraints = list(discrimination = A))
#' newdata_pred <- st_sf(wealth = rnorm(6), geometry = newdata)
#' predict(samples_pred, newdata = newdata_pred)
#' }
#'
#' @export
predict.spifa <- function (object, newdata = NULL, burnin = 0, thin = 1,
                                joint = FALSE, ...) {

  fit_args <- attr(object, "fit_args")
  predict_setup <- attr(object, "predict_setup")

  # Filter to the posterior samples
  niter <- posterior::niterations(object)
  object <- posterior::subset_draws(object, iteration = (burnin + 1):niter) |>
    posterior::thin_draws(thin)

  # Prediction I: for the observed subjects/locations
  has_newcoords <- inherits(newdata, "sf") || inherits(newdata, "sfc")
  if (fit_args$model_type %in% c("eifa", "cifa") |
      (fit_args$model_type == "cifa_pred" & is.null(newdata)) |
      (fit_args$model_type == "spifa" & !has_newcoords) |
      (fit_args$model_type == "spifa_pred" & is.null(newdata))) {
    result <- posterior::as_draws_array(as_matrix(object, "Theta"))
    if (!is.null(predict_setup$coordinates)) {
      attr(result, "newdata") <- sf::st_sf(geometry = predict_setup$coordinates)
    }
    return(result)
  }

  # Prediction II: for the new subjects/locations
  coordinates <- predict_setup$coordinates
  pred_terms <- predict_setup$pred_terms

  # Construct distances for spifa/spifa_pred
  if (fit_args$model_type == "cifa_pred") {
    npred <- nrow(newdata)
    newdist <- matrix(nrow = 0, ncol = 0)
    cross_distances <- matrix(nrow = 0, ncol = 0)
  } else {
    if (!has_newcoords) {
      stop("newdata must be an sf/sfc object for a spatial model.", call. = FALSE)
    }
    newcoords <- sf::st_geometry(newdata)
    if (any(sf::st_is(newcoords, c("POLYGON", "MULTIPOLYGON")))) {
      newcoords <- sf::st_centroid(newcoords)
    }
    npred <- length(newcoords)
    newdist <- matrix(as.numeric(sf::st_distance(newcoords)), npred, npred)
    cross_distances <- matrix(as.numeric(sf::st_distance(newcoords, coordinates)),
      npred, length(coordinates))
  }

  # Contruct newpredictors
  if (fit_args$model_type %in% c("cifa_pred", "spifa_pred") &&
      !all(all.vars(pred_terms) %in% names(newdata))) {
    stop("newdata is missing the predictor column(s) required (",
      paste(all.vars(pred_terms), collapse = ", "), ").", call. = FALSE)
  } else if (fit_args$model_type == "spifa") {
    newpredictors <- matrix(nrow = npred, ncol = 0)
  } else {
    mf_new <- model.frame(pred_terms, newdata, xlev = predict_setup$xlevels)
    newpredictors <- model.matrix(pred_terms, mf_new)
  }

  # List of options to call c++ function to predict
  nsamples <- posterior::ndraws(object)
  predict_args <- list(
    samples_theta = t(as_matrix(object, "Theta")),
    samples_corr_chol = t(as_matrix(object, "Chol")),
    samples_corr = t(as_matrix(object, "Corr")),
    samples_mgp_sd = t(as_matrix_opt(object, "T")),
    samples_mgp_phi = t(as_matrix_opt(object, "phi")),
    samples_betas = t(as_matrix_opt(object, "B")),
    response = fit_args$response,
    predictors = fit_args$predictors,
    newpredictors = newpredictors,
    distances = fit_args$distances,
    newdist = newdist,
    cross_distances = cross_distances,
    nobs = fit_args$nobs,
    nitems = fit_args$nitems,
    nfactors = fit_args$nfactors,
    ngp = fit_args$ngp,
    npred = npred,
    niter = nsamples,
    burnin = 0,
    thin = 1,
    constrain_L = fit_args$constrain_L,
    constrain_T = fit_args$constrain_T,
    constrain_V_sd = fit_args$constrain_V_sd,
    model_type = fit_args$model_type,
    joint = joint
  )

  # Predict calling c++ predict_cpp
  samples <- do.call(predict_cpp, predict_args)

  result <- posterior::as_draws_array(samples$theta)
  if (inherits(newdata, "sf") || inherits(newdata, "sfc")) {
    attr(result, "newdata") <- newdata
  }
  return(result)
}
