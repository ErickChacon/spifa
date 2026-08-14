
#' @title Print a spifa.list Object
#'
#' @description
#' Concise summary of a fitted \code{spifa.list} object (the raw output of
#' \code{\link{spifa}}): model type, data dimensions, MCMC settings, and the
#' stored parameter blocks — rather than dumping every raw MCMC sample
#' matrix to the console. Convert to a tidy tibble with
#' \code{\link{as_tibble.spifa.list}} for further inspection.
#'
#' @param x An object of class \code{spifa.list}, as returned by
#' \code{\link{spifa}}.
#' @param ... Further arguments passed to methods (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' data(ipixuna)
#' ipixuna_wide$items <- as.matrix(dplyr::select(ipixuna_wide, `Item 1`:`Item 10`))
#' samples <- spifa(items ~ 1, data = ipixuna_wide, nfactors = 2, niter = 20, thin = 1)
#' samples
#'
#' @export
print.spifa.list <- function (x, ...) {

  info <- attr(x, "model_info")
  cat("<spifa model fit>\n")
  cat("Model type: ", info$model_type, "\n", sep = "")
  cat(info$nobs, " respondents, ", info$nitems, " items, ", info$nfactors,
      " factors\n", sep = "")
  if (length(x) == 0) {
    cat(info$niter, " iterations requested (thin = ", info$thin,
        "); not executed (`execute = FALSE`)\n", sep = "")
  } else {
    cat(info$niter, " iterations, thinned by ", info$thin, " (", nrow(x[[1]]),
        " stored)\n", sep = "")
    cat("Parameter blocks: ", paste(names(x), collapse = ", "), "\n", sep = "")
  }
  cat("Use as_tibble() for a tidy tibble, summary() for posterior summaries.\n")
  invisible(x)
}

#' @title Convert MCMC Samples to a Wide Tibble
#'
#' @description
#' Converts the list of MCMC sample matrices returned by \code{\link{spifa}}
#' (class \code{spifa.list}) into a single wide \code{\link[tibble]{tibble}}
#' (class \code{spifa}), with one column per parameter and one row per
#' (thinned, post-burnin) MCMC iteration.
#'
#' @param x An object of class \code{spifa.list}, as returned by
#' \code{\link{spifa}}.
#' @param burnin Number of initial iterations to discard.
#' @param thin Thinning interval applied after discarding burn-in.
#' @param select Character vector of parameter blocks to keep (defaults to
#' all of them, i.e. \code{names(x)}).
#' @param ... Further arguments passed to methods (currently unused; present
#' for consistency with the \code{\link[tibble]{as_tibble}} generic).
#'
#' @return A wide \code{\link[tibble]{tibble}} of class \code{spifa}, with
#' the \code{"model_info"} attribute carried over from \code{x}.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' parameters <- attr(ipixuna_wide, "parameters")
#' L_a <- (parameters$discrimination != 0) * 1
#' ipixuna_wide$items <- as.matrix(dplyr::select(ipixuna_wide, `Item 1`:`Item 10`))
#' samples <- spifa(
#'   items ~ 1, data = ipixuna_wide, nfactors = 2,
#'   niter = 20, thin = 1, standardize = FALSE,
#'   constraints = list(discrimination = L_a, mgp = diag(2), resid_sd = rep(0.5, 2)))
#' samples_tib <- as_tibble(samples)
#' samples_tib
#' }
#'
#' @importFrom tibble as_tibble
#'
#' @export
as_tibble.spifa.list <- function (x, burnin = 0, thin = 1,
                              select = names(x), ...) {

  samples <- x
  # save model information
  model_info <- attr(samples, "model_info")
  # select, convert to tibble and subset
  samples <- samples[select]
  params <- names(samples)
  samples <- samples %>%
    purrr::map(tibble::as_tibble) %>%
    dplyr::bind_cols()
  niter <- nrow(samples)
  samples <- samples[seq(burnin + 1, niter, thin),]
  # assign model information and class
  attr(samples, "model_info") <- model_info
  class(samples) <- c("spifa", class(samples))
  return(samples)
}

#' @title Convert Samples to a spifa.list
#'
#' @description
#' Generic function converting posterior samples back into the
#' \code{spifa.list} representation (one matrix per parameter block). See
#' \code{\link{as_list.spifa}} for the \code{spifa} method.
#'
#' @param x An object to convert.
#' @param ... Further arguments passed to methods.
#'
#' @export
as_list <- function (x, ...) {
  UseMethod("as_list", x)
}

#' @title Convert a Wide Samples Tibble Back to a spifa.list
#'
#' @description
#' Inverse of \code{\link{as_tibble.spifa.list}}: splits a wide
#' \code{spifa} tibble (one column per scalar parameter) back into a named
#' list of matrices, one per parameter block (e.g. \code{c}, \code{a},
#' \code{theta}, ...). Used internally by \code{\link{predict.spifa}} and
#' \code{\link{dic.spifa}}.
#'
#' @param x An object of class \code{spifa}, as returned by
#' \code{\link{as_tibble.spifa.list}}.
#' @param ... Further arguments passed to methods (currently unused).
#'
#' @return An object of class \code{spifa.list}: a named list of matrices,
#' one per parameter block, with the \code{"model_info"} attribute carried
#' over from \code{samples_tib}.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' parameters <- attr(ipixuna_wide, "parameters")
#' L_a <- (parameters$discrimination != 0) * 1
#' ipixuna_wide$items <- as.matrix(dplyr::select(ipixuna_wide, `Item 1`:`Item 10`))
#' samples <- spifa(
#'   items ~ 1, data = ipixuna_wide, nfactors = 2,
#'   niter = 20, thin = 1, standardize = FALSE,
#'   constraints = list(discrimination = L_a, mgp = diag(2), resid_sd = rep(0.5, 2)))
#' samples_tib <- as_tibble(samples)
#' samples_list <- as_list(samples_tib)
#' names(samples_list)
#' }
#'
#' @importFrom tibble as_tibble
#'
#' @export
as_list.spifa <- function (x, ...) {

  samples_tib <- x
  model_info <- attr(samples_tib, "model_info")
  pars <- c("c", "a", "theta", "z", "corr_chol", "corr", "mgp_sd", "mgp_phi", "betas")
  names(pars) <- c("c", "A", "Theta", "Z", "Corr_chol", "Corr", "T", "mgp_phi", "B")

  param_label <- gsub("\\[.+\\]", "", names(samples_tib))
  param_label_unique <- unique(param_label)
  param_name_unique <- pars[param_label_unique]

  subset_tib <- function (param) {
    sub_tib <- samples_tib[param_label == param]
    sub_tib <- as.matrix(sub_tib)
    return(sub_tib)
  }

  samples <- lapply(param_label_unique, subset_tib)
  names(samples) <- param_name_unique
  class(samples) <- c("spifa.list", class(samples))
  attr(samples, "model_info") <- model_info
  return(samples)
}

#' @title Gather Parameters into a Long Format Tibble
#'
#' @description
#' Reshapes a wide \code{spifa} samples tibble (one column per parameter, as
#' produced by \code{\link{as_tibble.spifa.list}}) into long format (one row
#' per iteration/parameter pair), which is the shape expected by the
#' \code{\link{gg_trace}}/\code{\link{gg_density}} family of plotting
#' helpers.
#'
#' @param samples_wide A wide samples tibble, e.g. from
#' \code{\link{as_tibble.spifa.list}}.
#' @param each If not \code{NULL}, the number of columns that make up each
#' group of parameters (e.g. the number of items), used to additionally
#' split the gathered \code{Parameters} column into \code{group}/
#' \code{Parameter} columns.
#' @param keys Names to use for the group/parameter columns when \code{each}
#' is supplied.
#'
#' @return A long-format \code{\link[tibble]{tibble}} with columns
#' \code{iteration}, \code{Parameters}, and \code{Value} (plus \code{group}/
#' the second \code{keys} element when \code{each} is supplied).
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' parameters <- attr(ipixuna_wide, "parameters")
#' L_a <- (parameters$discrimination != 0) * 1
#' ipixuna_wide$items <- as.matrix(dplyr::select(ipixuna_wide, `Item 1`:`Item 10`))
#' samples <- spifa(
#'   items ~ 1, data = ipixuna_wide, nfactors = 2,
#'   niter = 20, thin = 1, standardize = FALSE,
#'   constraints = list(discrimination = L_a, mgp = diag(2), resid_sd = rep(0.5, 2)))
#' wide <- as_tibble(samples, select = "c")
#' long <- gather.spifa(wide)
#' long
#' }
#'
#' @export
gather.spifa <- function (samples_wide, each = NULL,
                           keys = c("group", "Parameter")) {

  # Convert to long format
  samples_long <- samples_wide %>%
    tibble::as_tibble() %>%
    dplyr::mutate(iteration = 1:dplyr::n()) %>%
    tidyr::gather(Parameters, Value, -iteration, factor_key = TRUE)

  if (!is.null(each)) {

    # Auxiliary variables to group
    groups <- paste0(keys[1], rep(1:each, ncol(samples_wide)/each))
    groups <- factor(groups, unique(groups))
    var <- paste0(keys[2], rep(1:(ncol(samples_wide)/each), each = each))
    names(groups) <- levels(samples_long$Parameters)
    names(var) <- levels(samples_long$Parameters)

    # Group parameters
    samples_long <- samples_long %>%
      dplyr::mutate(groups = groups[Parameters], var = var[Parameters]) %>%
    dplyr::select(-Parameters) %>%
    tidyr::spread(var, Value)

  }

  return(samples_long)
}


#' @title Posterior Summary of MCMC Samples
#'
#' @description
#' Computes the posterior median and 80\%/95\% credible intervals (10th,
#' 2.5th, 50th, 90th, and 97.5th percentiles) for every parameter in a
#' \code{spifa} samples tibble.
#'
#' @param object A wide \code{spifa} samples tibble, e.g. from
#' \code{\link{as_tibble.spifa.list}}.
#' @param select Character vector of column names to summarise (defaults to
#' all of them).
#' @param ... Further arguments passed to methods (currently unused).
#'
#' @return A \code{\link[tibble]{tibble}} with one row per parameter and
#' columns \code{Parameters}, \code{2.5\%}, \code{10\%}, \code{50\%},
#' \code{90\%}, \code{97.5\%}.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' parameters <- attr(ipixuna_wide, "parameters")
#' L_a <- (parameters$discrimination != 0) * 1
#' ipixuna_wide$items <- as.matrix(dplyr::select(ipixuna_wide, `Item 1`:`Item 10`))
#' samples <- spifa(
#'   items ~ 1, data = ipixuna_wide, nfactors = 2,
#'   niter = 20, thin = 1, standardize = FALSE,
#'   constraints = list(discrimination = L_a, mgp = diag(2), resid_sd = rep(0.5, 2)))
#' samples_tib <- as_tibble(samples, select = "c")
#' summary(samples_tib)
#' }
#'
#' @export
summary.spifa <- function (object, select = names(object), ...) {

  df <- object[select]
  df_names <- names(df)
  quantiles <- vapply(df, quantile, numeric(5),
                       probs = c(0.025, 0.1, 0.5, 0.9, 0.975), na.rm = TRUE)
  df <- tibble::as_tibble(t(quantiles)) %>%
    dplyr::mutate(Parameters = factor(df_names, df_names))
  df <- df[c(ncol(df), 1:(ncol(df)-1))]
  return(df)
}

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
#' @param x A wide \code{spifa} samples tibble, e.g. from
#' \code{\link{as_tibble.spifa.list}}.
#' @param ... Further arguments passed to methods (currently unused).
#'
#' @return A list with elements \code{average_of_deviance},
#' \code{n_effec_params} (effective number of parameters), and \code{dic}.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' parameters <- attr(ipixuna_wide, "parameters")
#' L_a <- (parameters$discrimination != 0) * 1
#' ipixuna_wide$items <- as.matrix(dplyr::select(ipixuna_wide, `Item 1`:`Item 10`))
#' samples <- spifa(
#'   items ~ 1, data = ipixuna_wide, nfactors = 2,
#'   niter = 20, thin = 1, standardize = FALSE,
#'   constraints = list(discrimination = L_a, mgp = diag(2), resid_sd = rep(0.5, 2)))
#' samples_tib <- as_tibble(samples)
#' dic(samples_tib)
#' }
#'
#' @export
dic.spifa <- function (x, ...) {

  object <- x
  # convert to spifa.list
  samples <- as_list(object)

  # DIC calling c++ dic_cpp
  dic <- dic_cpp(y = attr(object, "model_info")$response, c = samples$c,
                 a = samples$a, theta = samples$theta,
                 n = attr(object, "model_info")$nobs,
                 q = attr(object, "model_info")$nitems,
                 m = attr(object, "model_info")$nfactors,
                 L = attr(object, "model_info")$constrain_L)

  return(dic)
}

#' @title Predict from a Fitted spifa Model
#'
#' @description
#' Predicts the latent factors (and, for spatial models, the underlying
#' spatial process) at new locations and/or for new predictor values, using
#' the posterior samples from a fitted \code{\link{spifa}} model.
#'
#' @details
#' If the fitted model has no spatial or predictor structure (\code{eifa} or
#' \code{cifa}), or if neither \code{newcoords} nor \code{newdata} is
#' supplied for a model that has one, there is nothing to predict beyond the
#' posterior means of the latent abilities already available from the fitted
#' samples, and the function currently returns early. Otherwise, prediction
#' for the new locations and/or predictor values is delegated to the
#' \code{C++} sampler.
#'
#' @param object A wide \code{spifa} samples tibble, e.g. from
#' \code{\link{as_tibble.spifa.list}}.
#' @param newdata New values of the predictors used on the right-hand side
#' of \code{formula} when the model was fitted (only needed for models
#' with predictors).
#' @param newcoords New spatial coordinates to predict at (only needed for
#' spatial, i.e. \code{spifa}/\code{spifa_pred}, models).
#' @param burnin Number of initial (post-fitting) iterations to discard
#' before using the posterior samples for prediction.
#' @param thin Thinning interval applied to the posterior samples used for
#' prediction.
#' @param se.fit Currently unused; reserved for returning prediction
#' standard errors.
#' @param ... Further arguments (currently unused).
#'
#' @return A list of posterior predictive samples/summaries for the
#' requested new locations and/or predictor values.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' parameters <- attr(ipixuna_wide, "parameters")
#' L_a <- (parameters$discrimination != 0) * 1
#' ipixuna_sf <- sf::st_as_sf(ipixuna_wide)
#' ipixuna_sf$items <- as.matrix(dplyr::select(ipixuna_wide, `Item 1`:`Item 10`))
#' samples <- spifa(
#'   items ~ 1, data = ipixuna_sf,
#'   nfactors = 2, niter = 5, thin = 1, standardize = FALSE,
#'   constraints = list(discrimination = L_a, mgp = diag(2), resid_sd = rep(0.5, 2)))
#' samples_tib <- as_tibble(samples)
#' newcoords <- sf::st_coordinates(ipixuna_wide$coords)[1:5, , drop = FALSE]
#' predict(samples_tib, newcoords = newcoords)
#' }
#'
#' @export
predict.spifa <- function (object, newdata = NULL, newcoords = NULL, burnin = 0,
                                thin = 1, se.fit = FALSE, ...) {

  # if (inherits(object, "spifa")) {
  #   object <- as.list(object)
  # } else { #} if (!inherits(object, "spifa.list")) {
  #   stop("class of object must be spifa")
  # }

  # convert to spifa.list
  object <- as_list(object)

  # Information of model inference
  info <- attr(object, "model_info")

  # Prediction I: for the observed subjects
  if (info$model_type %in% c("eifa", "cifa") |
      (info$model_type == "cifa_pred" & is.null(newdata)) |
      (info$model_type == "spifa" & is.null(newcoords)) |
      (info$model_type == "spifa_pred" & is.null(newcoords) & is.null(newdata))) {
    return("here is only the means of the latent abilities")
  }

  # Prediction II: for the new subjects or new locations

  # Distances between predictive locations
  if (is.null(newcoords)) {
    newdist <- matrix(NA)
    cross_distances <- matrix(NA)
  } else {
    npred1 <- nrow(newcoords)
    newdist <- as.matrix(dist(newcoords))
      if (inherits(info$coordinates, "sfc")) {
        info$coordinates <- sf::st_transform(info$coordinates, crs = 3857) %>%
          sf::st_coordinates()
      }
    cross_distances <- as.matrix(pdist::pdist(newcoords, as.matrix(info$coordinates)))
  }

  # New data about predictors
  if (is.null(newdata)) {
    npred2 <- NULL
  } else {
    npred2 <- nrow(newdata)
  }

  # Obtain number of subjects or locations to predict
  if (!is.null(npred2) && exists("npred1")) {
    if (npred1 == npred2) {
      npred = npred1
    } else {
      stop(sprintf("number of rows of '%s' of '%s' must be the same",
                   "newcoords", "newdata"))
    }
  } else if (exists("npred1")) {
    npred = npred1
  } else if (!is.null(npred2)) {
    npred = npred2
  }

  # Build predictors for the new locations/subjects. If the model has
  # predictors but newdata was not supplied (e.g. predicting the spatial
  # component only), fall back to zeros so the predictor contribution is
  # held at its reference level instead of crashing with mismatched
  # dimensions.
  if (!is.null(newdata)) {
    newpredictors <- newdata
  } else if (info$model_type %in% c("cifa_pred", "spifa_pred")) {
    newpredictors <- matrix(0, npred, ncol(info$predictors))
  } else {
    newpredictors <- matrix(NA)
  }

  # Information about number of posterior samples to use
  nsamples <- nrow(object$z)
  # if (is.null(burnin)) burnin <- as.integer(nsamples / 2)
  # if (is.null(thin)) thin <- as.integer( (nsamples - burnin) / 1000)

  # List of options to call c++ function to predict
  pred_list <- list(samples_theta = t(object$theta),
    samples_corr_chol = t(object$corr_chol), samples_corr = t(object$corr),
    samples_mgp_sd = t(object$mgp_sd), samples_mgp_phi = t(object$mgp_phi),
    samples_betas = t(object$betas),
    response = info$response, predictors = info$predictors, newpredictors = newpredictors,
    distances = info$distances, newdist = newdist, cross_distances = cross_distances,
    nobs = info$nobs, nitems = info$nitems, nfactors = info$nfactors, ngp = info$ngp,
    npred = npred, niter = nsamples, burnin = burnin, thin = thin,
    constrain_L = info$constrain_L, constrain_T = info$constrain_T,
    constrain_V_sd = info$constrain_V_sd,
    model_type = info$model_type
    )

  # Predict calling c++ predict_spifa_cpp
  prediction <- do.call(predict_spifa_cpp, pred_list)

  return(prediction)
  # return(pred_list)
}
