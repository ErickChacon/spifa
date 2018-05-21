
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
  samples <- samples[select]
  params <- names(samples)
  samples <- samples %>%
    purrr::map(tibble::as_tibble) %>%
    dplyr::bind_cols()
  niter <- nrow(samples)
  samples <- samples[seq(burnin + 1, niter, thin),]
  class(samples) <- c("spifa", class(samples))
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
predict.spifa.list <- function (object, newdata = NULL, newcoords = NULL, burnin = NULL,
                                thin = NULL, se.fit = FALSE, ...) {

  # Information of inference
  info <- attr(object, "model_info")

  # Distances between predictive and observed locations
  npred <- nrow(newcoords)
  if (is.null(newcoords)) {
    cross_distances <- matrix(NA)
  } else {
    cross_distances <- as.matrix(pdist::pdist(newcoords, info$coordinates))
  }

  # New data about predictors
  if (is.null(newdata)) {
    newpredictors <- matrix(NA)
  } else {
    newpredictors <- newdata
  }

  # Information about number of posterior samples to use
  nsamples <- ncol(object$z)
  if (is.null(burnin)) burnin <- as.integer(nsamples / 2)
  if (is.null(thin)) thin <- as.integer( (nsamples - burnin) / 1000)


  # Predict calling c++ predict
  prediction <- predict_spifa_cpp(
    response = response, predictors = info$predictors, newpredictors = newpredictors,
    distances = info$distances, cross_distances = cross_distances,
    nobs = nobs, nitems = nitems, nfactors = nfactors, ngp = ngp, npred = npred,
    burnin = burnin, thin = thin,
    constrain_L = constrain_L, constrain_T = constrain_T, constrain_V_sd = constrain_V_sd,
    model_type = model_type
    )

  return(prediction)

}
