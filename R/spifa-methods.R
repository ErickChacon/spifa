
#' @title Convert Samples to Tibble
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' as_tibble(samples)
#' as_tibble(samples, "beta")
#'
#' @importFrom tibble as_tibble
#'
#' @export
as_tibble.spifa.list <- function (samples, burnin = 0, thin = 1,
                              select = names(samples)) {

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

as_list <- function (x, ...) {
  UseMethod("as_list", x)
}

#' @title Convert Tibble to spifa.list
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' as_tibble(samples)
#' as_tibble(samples, "beta")
#'
#' @importFrom tibble as_tibble
#'
#' @export
as_list.spifa <- function (samples_tib) {

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
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' wide <- as_tibble(samples, "beta")
#' (long <- gather(wide))
#' class(long)
#'
#' @export
gather.spifa <- function (samples_wide, each = NULL,
                           keys = c("group", "Parameter")) {

  # Convert to long format
  samples_long <- samples_wide %>%
    tibble::as_tibble() %>%
    dplyr::mutate(iteration = 1:n()) %>%
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


#' @title Summary of Samples
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' summary(samples)
#'
#' @export
summary.spifa <- function (samples, select = names(samples)) {

  # df <- as_tibble.spifa(samples, select = select)
  df <- samples
  df_names <- names(df)
  df <- df %>%
    map(function (x) quantile(x, c(0.025, 0.1, 0.5, 0.9, 0.975))) %>%
    reduce(rbind) %>%
    as_tibble() %>%
    mutate(Parameters = factor(df_names, df_names))
  df <- df[c(ncol(df), 1:(ncol(df)-1))]
  return(df)
}

#' @export
dic <- function (x, ...) {
  UseMethod("dic", x)
}


#' @title Deviance Information Criterio for the Spifa Model
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' summary(samples)
#'
#' @export
dic.spifa <- function (object) {

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

#' @title Prediction Spatial Multidimensional Item Response Model with Predictors
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' 
#'
#' @export
predict.spifa <- function (object, newdata = NULL, newcoords = NULL, burnin = 0,
                                thin = 1, se.fit = FALSE, what = NULL, ...) {

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
    newpredictors <- matrix(NA)
  } else {
    npred2 <- nrow(newdata)
    newpredictors <- newdata
  }

  # Obtain number of subjects or locations to predict
  if (exists("npred1") & exists("npred2")) {
    if (npred1 == npred2) {
      npred = npred1
    } else {
      stop(sprintf("number of rows of '%s' of '%s' must be the same",
                   "newcoords", "newdata"))
    }
  } else if (exists("npred1")) {
    npred = npred1
  } else if (exists("npred2")) {
    npred = npred2
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
  if (is.null(what)) {
    prediction <- do.call(predict_spifa_cpp, pred_list)
  } else {
    prediction <- do.call(predict2_spifa_cpp, pred_list)
  }

  return(prediction)
  # return(pred_list)
}
