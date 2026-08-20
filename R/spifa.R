#' @title Fit an (Spatial) Item Factor Analysis Model
#'
#' @description
#' Fits exploratory, confirmatory, and spatial item factor analysis (IFA)
#' models for binary responses using full Bayesian inference. The model
#' represents each binary response as a thresholded continuous auxiliary
#' variable explained by \code{nfactors} latent abilities, optionally
#' extended with linear predictors and/or a multivariate Gaussian process to
#' capture spatial dependence in the latent factors (see the "spifa-ipixuna"
#' vignette for a full worked example). Inference is done via Gibbs sampling
#' with adaptive Metropolis-Hastings updates for the spatial range and
#' correlation parameters.
#'
#' @details
#' The type of model fitted is determined automatically from \code{formula}
#' and the class of \code{data}: a one-sided right-hand side (\code{items ~
#' 1}) with unrestricted \code{constraints$discrimination} gives exploratory
#' IFA (EIFA); the same with a restricted \code{constraints$discrimination}
#' gives confirmatory IFA (CIFA); adding predictors to the right-hand side
#' (e.g. \code{items ~ x1}) gives CIFA with predictors; \code{data} being an
#' \code{\link[sf]{sf}} object adds a spatial Gaussian process on the latent
#' factors, giving spatial IFA (SPIFA), with or without predictors.
#'
#' The left-hand side of \code{formula} must be a single symbol naming a
#' matrix-valued column of \code{data} (\code{nobs x nitems}, one row per
#' respondent, one column per binary item) — the same mechanism base R uses
#' for multivariate \code{\link[stats]{lm}}. Build it with \code{\link{I}}
#' (or ordinary \code{$<-} assignment) so it survives as a matrix column
#' rather than being flattened into separate columns, e.g.:
#' \preformatted{
#' items <- as.matrix(dplyr::select(data, `Item 1`:`Item 10`))
#' data$items <- items
#' spifa(items ~ x1, data = data, nfactors = 2)
#' }
#'
#' @param formula A two-sided formula \code{items ~ predictors}. The
#' left-hand side must be a single symbol naming a matrix-valued column of
#' \code{data} holding the binary item responses (see Details). The
#' right-hand side specifies predictors for the latent factors (e.g.
#' \code{~ x1 + x2}); use \code{items ~ 1} for no predictors.
#' @param data A data frame containing the item-response matrix column named
#' on the left-hand side of \code{formula} and, if used, the predictor
#' columns named on its right-hand side. If \code{data} is an
#' \code{\link[sf]{sf}} object, its geometry (\code{\link[sf]{st_geometry}})
#' is used as the spatial coordinates and a spatial Gaussian process is added
#' to the model; pass a plain (non-\code{sf}) data frame — e.g. via
#' \code{\link[sf]{st_set_geometry}(data, NULL)} — for a non-spatial model.
#' @param nfactors Number of latent factors (dimensions of the ability
#' construct).
#' @param ngp Number of independent Gaussian processes used to build the
#' (possibly restricted) multivariate Gaussian process for the latent
#' factors. Defaults to \code{nfactors} (one GP per factor). Only relevant
#' when \code{data} is an \code{sf} object; set to \code{0} to fit a
#' non-spatial model even when \code{data} has a geometry column (the
#' geometry and \code{ngp} are otherwise ignored in that case).
#' @param niter Number of MCMC iterations to run (after \code{burnin}) and
#' store.
#' @param thin Thinning interval for the stored MCMC samples.
#' @param burnin Number of initial MCMC iterations to discard. These
#' iterations still run (and the adaptive Metropolis-Hastings proposals for
#' \code{mgp_sd}/\code{mgp_range}/\code{resid_corr} still adapt through
#' them), but they are never stored, so \code{niter} counts only the
#' iterations that end up in the returned samples. \code{0} by default (no
#' iterations discarded during fitting -- the previous behaviour). Prefer
#' this over discarding a prefix of the samples afterwards (e.g. via
#' \code{summary(..., burnin = )}): iterations dropped here were never
#' stored, so they don't cost memory or thinning-index arithmetic, and the
#' adaptive proposals get to keep converging across the burnin/niter
#' boundary rather than someone accidentally analysing them as if they were
#' post-adaptation draws.
#' @param standardize Logical; if \code{TRUE} (default), predictors are
#' standardized before fitting.
#' @param constraints Named list of constraints associated to the factor model. Accepted
#' names are `discrimination`, `mgp`, and `resid_sd`. The restrictions on the
#' discrimination paramater should be placed in the element `discrimination` with same
#' dimensions as the discrimination matrix (nitems x nfactors). A value of 0 indicates that
#' the link betwen the item and the factor is disabled and 1 indicates that it remains
#' active and the coefficient associated will be estimated. The restrictions for the
#' multivariate Gaussian process should be placed in the element `mgp` with dimensions
#' nfactors x ngp, such as a value of 0 indicates a link disconnected between the factor
#' and the (independent) GP while 1 indicates that it remains active. The restrictions with
#' respect to the standard deviation of the latent factors' residual term should be placed
#' in the element `resid_sd`, which should be a vector (length nfactors) providing the fixed
#' values for that standard deviation (paired with `priors$resid_corr`, together they
#' parameterize the residual covariance — see `dev/design/scope.md`). If the model includes
#' predictors or a Gaussian process, it is recomended to be lower than 1.
#'
#' @param priors Named list of initial values and prior hyperparameters, one
#' element per parameter block: `easiness`, `discrimination`, `effect`
#' (predictor effect on the latent factors), `resid_corr` (correlation of the
#' latent factors' residual term, paired with `constraints$resid_sd`),
#' `mgp_sd` (multivariate Gaussian process standard deviations), and
#' `mgp_range` (multivariate Gaussian process scale parameters). Each
#' element (except `resid_corr`) accepts `initial`, `mean`, and `sd`;
#' `resid_corr` accepts `initial` and `eta` (the LKJ prior shape parameter).
#' See `dev/design/scope.md` for the mapping between these names and the
#' paper's notation (`c`, `a`, `B`, `R`).
#'
#' @param adaptive Named list of hyperparameters associated with the adaptive sampling.
#' The adaptive sampling is done jointly for the `correlation` parameters, `standard
#' deviation of the gps` and `scale parameter of the gps`. The matrix `Sigma` can be
#' provided as the full covariance matrix of these parameters for the proposal
#' distribution. Otherwise, part of this matrix can be provided by using the elements
#' `Sigma`, `Sigma_resid_corr`, `Sigma_mgp_sd` and `Sigma_mgp_range`. Additional elements are `scale`,
#' `C`, `alpha` and `accep_prob` which are hyperparameters of the adaptive sampling
#' proposed in Andrieu and Thomas (2008).
#'
#' @param execute Logical value to run sampler or not. TRUE by default.
#'
#' @return
#' An object of (informal) class \code{spifa.list}: a named list of MCMC
#' sample matrices (one entry per parameter block, e.g. \code{c}, \code{a},
#' \code{theta}, \code{corr}, \code{betas}, ...), with an attribute
#' \code{"model_info"} recording the data and options used to fit the model
#' (needed by \code{\link{predict.spifa}} and \code{\link{dic}}). Convert it
#' to a tidy \code{\link[tibble]{tibble}} with \code{\link{as_tibble.spifa}}.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' data(ipixuna)
#'
#' # true discrimination structure used to simulate ipixuna
#' parameters <- attr(ipixuna, "parameters")
#' L_a <- (parameters$discrimination != 0) * 1
#' nfactors <- ncol(parameters$discrimination)
#'
#' # confirmatory item factor analysis (non-spatial: ngp = 0; small niter
#' # for a fast example)
#' samples <- spifa(
#'   items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
#'   niter = 20, thin = 1, standardize = FALSE,
#'   constraints = list(discrimination = L_a, resid_sd = rep(0.5, nfactors)))
#' summary(samples, select = c("c", "a"))
#'
#' @export
spifa <- function(formula, data, nfactors, ngp = nfactors,
    niter = 1000, thin = 10, burnin = 0, standardize = TRUE,
    constraints = list(discrimination = NULL, mgp = NULL, resid_sd = rep(1, nfactors)),
    priors = list(
      easiness = list(initial = NULL, mean = NULL, sd = NULL),
      discrimination = list(initial = NULL, mean = NULL, sd = NULL),
      effect = list(initial = NULL, mean = NULL, sd = NULL),
      resid_corr = list(initial = NULL, eta = 1.5),
      mgp_sd = list(initial = NULL, mean = NULL, sd = NULL),
      mgp_range = list(initial = NULL, mean = NULL, sd = NULL)),
    adaptive = list(Sigma = NULL, Sigma_resid_corr = NULL, Sigma_mgp_sd = NULL, Sigma_mgp_range = NULL,
                    scale = 1, C = 0.7, alpha = 0.8, accep_prob = 0.234),
    execute = TRUE) {

  # Determine dimensions, items and predictors
  mf <- model.frame(formula, data)
  response <- model.response(mf)
  if (!is.matrix(response)) stop("The left-hand side of 'formula' must be a matrix")
  nobs <- nrow(response)
  nitems <- ncol(response)
  response <- as.numeric(response)

  if (length(attr(terms(mf), "term.labels")) == 0) {
    predictors <- NULL
  } else {
    predictors_terms <- delete.response(terms(mf))
    attr(predictors_terms, "intercept") <- 0
    predictors <- model.matrix(predictors_terms, mf)
  }

  # Coordinates: spatial structure is inferred from class(data). Pass a
  # plain (non-sf) data frame, e.g. via sf::st_set_geometry(data, NULL), or
  # ngp = 0, for a non-spatial model.
  coordinates <- if (inherits(data, "sf")) sf::st_geometry(data) else NULL
  has_gp <- !is.null(coordinates) && ngp > 0

  # Restrictions for discrimination parameters and Gaussian process loadings
  constrain_L_explo <- matrix(NA, nitems, nfactors)
  constrain_L_explo <- lower.tri(constrain_L_explo, diag = TRUE) * 1
  constrain_L <- check_param_mat(constraints, "discrimination", c(nitems, nfactors), constrain_L_explo)
  constrain_T <- check_param_mat(constraints, "mgp", c(nfactors, ngp), diag(1, nfactors, ngp))

  # Sizes
  nsigmas <- sum(constrain_T)
  n_corr <- nfactors * (nfactors - 1) / 2

  # Detect type of model to be fitted: EIFA, CIFA, CIFA_PRED, SPIFA, SPIFA_PRED
  if (has_gp) {
    if (!is.null(predictors)) {
      model_type = "spifa_pred"
      constrain_V_sd <- check_param_vec(constraints, "resid_sd", nfactors, 0.2)
    } else {
      model_type = "spifa"
      constrain_V_sd <- check_param_vec(constraints, "resid_sd", nfactors, 0.2)
    }
  } else if (!is.null(predictors)) {
    model_type = "cifa_pred"
    constrain_V_sd <- check_param_vec(constraints, "resid_sd", nfactors, 0.3)
  } else if (all(constrain_L == constrain_L_explo)) {
    model_type = "eifa"
    constrain_V_sd <- check_param_vec(constraints, "resid_sd", nfactors, 1)
  } else {
    model_type = "cifa"
    constrain_V_sd <- check_param_vec(constraints, "resid_sd", nfactors, 1)
  }

  # Optional arguments for easiness parameters (c)
  c_prior_mean <- check_param_vec(priors$easiness, "mean", nitems, 0)
  c_prior_sd <- check_param_vec(priors$easiness, "sd", nitems, 1)
  c_initial <- check_param_vec(priors$easiness, "initial", nitems,
                               rnorm(nitems, c_prior_mean, c_prior_sd))

  # Optional arguments for discrimination parameters (A)
  A_prior_mean <-
    check_param_mat2(priors$discrimination, "mean", c(nitems, nfactors), diag(1, nitems, nfactors))
  A_prior_sd <-
    check_param_mat2(priors$discrimination, "sd", c(nitems, nfactors), 1-diag(0.55, nitems, nfactors))
  A_initial <-
    check_param_mat2(priors$discrimination, "initial", c(nitems, nfactors), A_prior_mean)

  # Adaptive Metropolis-Hastings arguments for proposed covariance matrix
  adap_Sigma_R <- check_param_matdiag(adaptive, "Sigma_resid_corr", n_corr, diag(n_corr) * 0.001)
  adap_Sigma_gp_sd <-
    check_param_matdiag(adaptive, "Sigma_mgp_sd", nsigmas, diag(nsigmas) * 0.001)
  adap_Sigma_gp_phi <-
    check_param_matdiag(adaptive, "Sigma_mgp_range", ngp, diag(ngp) * 0.001)
  adap_scale <- ifelse(is.null(adaptive$scale), 1, adaptive$scale)
  adap_C <- ifelse(is.null(adaptive$C), 0.7, adaptive$C)
  adap_alpha <- ifelse(is.null(adaptive$alpha), 0.8, adaptive$alpha)
  adap_accep_prob <- ifelse(is.null(adaptive$accep_prob), 0.234, adaptive$accep_prob)

  # Create general sigma proposal in order: gp_sd, gp_phi, corr_free
  if (!has_gp) {
    if (is.null(adaptive$Sigma)) {
      adap_Sigma <- adap_Sigma_R
    } else if (sum(dim(adaptive$Sigma) == c(n_corr, n_corr)) == 2) {
      adap_Sigma <- adaptive$Sigma
    }
  } else {
    if (is.null(adaptive$Sigma)) {
      adap_Sigma <- matrix(0, nsigmas + ngp + n_corr, nsigmas + ngp + n_corr)
      adap_Sigma[seq_len(nsigmas), seq_len(nsigmas)] <- adap_Sigma_gp_sd
      adap_Sigma[nsigmas + seq_len(ngp), nsigmas + seq_len(ngp)] <- adap_Sigma_gp_phi
      adap_Sigma[nsigmas + ngp + seq_len(n_corr), nsigmas + ngp + seq_len(n_corr)] <- adap_Sigma_R
    } else if (sum(dim(adaptive$Sigma) == rep(nsigmas + ngp + n_corr, 2)) == 2) {
      adap_Sigma <- adaptive$Sigma
    }
  }

  # Optional arguments for parameter of residual correlation R
  if (is.null(priors$resid_corr$initial)) {
    R_initial <- diag(nfactors)
  } else if (sum(dim(priors$resid_corr$initial) == c(nfactors, nfactors)) == 2) {
    if (all(diag(priors$resid_corr$initial) == 1)) {
      R_initial <- priors$resid_corr$initial
    } else {
      stop("element 'initial' of 'priors$resid_corr' argument is not a correlation matrix")
    }
  } else {
    stop("element 'initial' of 'priors$resid_corr' argument must be of dimension ",
         sprintf("c(%i, %i)", nfactors, nfactors))
  }

  R_prior_eta <- ifelse(is.null(priors$resid_corr$eta), 1, priors$resid_corr$eta)

  # Optional arguments for parameter of fixed effects (Beta)
  if (is.null(predictors)) {
    B_prior_mean <- matrix(NA, 1, nfactors)
    B_prior_sd <- matrix(NA, 1, nfactors)
    B_initial <- matrix(NA, 1, nfactors)
  } else {
    npred <- ncol(predictors)
    B_prior_mean <- check_param_mat2(priors$effect, "mean", c(npred, nfactors), 0)
    B_prior_sd <- check_param_mat2(priors$effect, "sd", c(npred, nfactors), 1)
    B_initial <- check_param_mat2(priors$effect, "initial", c(npred, nfactors), B_prior_mean)
  }

  # Optional arguments for GP standard deviations and  scale parameters
  if (!has_gp) {
    sigmas_gp_mean <- rep(NA_real_, nsigmas)
    sigmas_gp_sd <- rep(NA_real_, nsigmas)
    sigmas_gp_initial <- rep(NA_real_, nsigmas)
    phi_gp_mean <- rep(NA_real_, ngp)
    phi_gp_sd <- rep(NA_real_, ngp)
    phi_gp_initial <- rep(NA_real_, ngp)
  } else {
    sigmas_gp_mean <- check_param_vec(priors$mgp_sd, "mean", nsigmas, 0.6)
    sigmas_gp_sd <- check_param_vec(priors$mgp_sd, "sd", nsigmas, 0.2)
    sigmas_gp_initial <- check_param_vec(priors$mgp_sd, "initial", nsigmas, sigmas_gp_mean)
    phi_gp_mean <- check_param_vec(priors$mgp_range, "mean", ngp, 0.05)
    phi_gp_sd <- check_param_vec(priors$mgp_range, "sd", ngp, 0.2)
    phi_gp_initial <- check_param_vec(priors$mgp_range, "initial", ngp, phi_gp_mean)
  }

  # Compute predictors and distances as matrices
  if (is.null(predictors))  predictors <- matrix(NA)
  if (!has_gp) {
    distances <- matrix(NA)
  } else {
    distances <- sf::st_transform(coordinates, crs = 3857) %>%
      sf::st_coordinates()  %>%
      dist()  %>%
      as.matrix()
  }

  # List of options to call c++ spifa function
  model_info <- list(
    response = response, predictors = predictors, distances = distances,
    nobs = nobs, nitems = nitems, nfactors = nfactors, ngp = ngp,
    niter = niter, thin = thin, burnin = burnin, standardize = standardize,
    constrain_L = constrain_L, constrain_T = constrain_T, constrain_V_sd = constrain_V_sd,
    adap_Sigma = adap_Sigma, adap_scale = adap_scale, adap_C = adap_C,
    adap_alpha = adap_alpha, adap_accep_prob = adap_accep_prob,
    c_initial = c_initial, c_prior_mean = c_prior_mean, c_prior_sd = c_prior_sd,
    A_initial = A_initial, A_prior_mean = A_prior_mean, A_prior_sd = A_prior_sd,
    R_initial = R_initial, R_prior_eta = R_prior_eta,
    B_initial = B_initial, B_prior_mean = B_prior_mean, B_prior_sd = B_prior_sd,
    sigmas_gp_initial = sigmas_gp_initial, sigmas_gp_mean = sigmas_gp_mean,
    sigmas_gp_sd = sigmas_gp_sd,
    phi_gp_initial = phi_gp_initial, phi_gp_mean = phi_gp_mean, phi_gp_sd = phi_gp_sd,
    model_type = model_type
    )

  # Execute c++ if requested
  if (execute) {
    # Execute model calling c++ spifa function
    samples <- do.call(spifa_cpp, model_info)
    # Update model_info list
    constrain_V_sd <- attr(samples, "V_sd")
    attr(samples, "V_sd") <- NULL
    model_info$constrain_V_sd <- constrain_V_sd
    model_info <- append(model_info, list(coordinates = coordinates), 2)

    # src/ifa.cpp always returns all 9 blocks (c, a, theta, z, corr_chol,
    # corr, mgp_sd, mgp_phi, betas), NA-filled for blocks that aren't part
    # of this model_type (e.g. betas for eifa/cifa, mgp_sd/mgp_phi for
    # anything non-spatial) -- those parameters were never actually sampled,
    # so drop them here rather than shipping structural NAs. This also lets
    # `bayesplot`, which errors on any NA, work on the returned draws_array
    # directly without further filtering.
    has_predictors <- model_type %in% c("cifa_pred", "spifa_pred")
    has_coordinates <- model_type %in% c("spifa", "spifa_pred")
    blocks_drop <- c(
      if (!has_predictors) "betas",
      if (!has_coordinates) c("mgp_sd", "mgp_phi"))
    samples <- samples[!names(samples) %in% blocks_drop]

    # Flatten the named list of per-block sample matrices (each block's own
    # column names, e.g. "c[1]", "A[1,1]", already set in src/ifa.cpp) into
    # a single posterior::draws_array. "spifa" is prepended to the class so
    # package methods (print/summary/as_list) dispatch first, while every
    # posterior/bayesplot function that works on class "draws_array" keeps
    # working directly on this object with no conversion call.
    flat <- do.call(cbind, samples)
    arr <- array(flat, dim = c(nrow(flat), 1, ncol(flat)),
                 dimnames = list(NULL, NULL, colnames(flat)))
    samples <- posterior::as_draws_array(arr)
  } else {
    samples <- list()
  }

  # Add model_info to the fitted object.
  attr(samples, "model_info") <- model_info
  class(samples) <- unique(c("spifa", class(samples)))

  return(samples)
}

# Internal validators used throughout spifa() to fill in defaults and check
# dimensions for the many optional prior/initial-value arguments.

check_param_vec <- function (param_list, element, dimension, default) {
  argument <- deparse(substitute(param_list))
  if (is.null(param_list[[element]])) {
    if (length(default) == 1) {
      output <- rep(default, dimension)
    } else {
      output <- default
    }
  } else if (length(param_list[[element]]) == 1) {
    output <- rep(param_list[[element]], dimension)
  } else if (length(param_list[[element]]) == dimension) {
    output <- param_list[[element]]
  } else {
    stop(sprintf("element '%s' of '%s' must be of length 1 or %i",
                 element, argument, dimension))
  }
  return(output)
}


check_param_mat <- function (param_list, element, dimensions, default) {
  # It only accepts matrices
  argument <- deparse(substitute(param_list))
  if (is.null(param_list[[element]])) {
    output <- default
  } else if (sum(dim(param_list[[element]]) == dimensions) == 2) {
    output <- param_list[[element]]
  } else {
    stop(sprintf("element '%s' of '%s' must be of dimension c(%i, %i)",
                 element, argument, dimensions[[1]], dimensions[[2]]))
  }
  return(output)
}


check_param_mat2 <- function (param_list, element, dimensions, default) {
  # It accepts matrices and scalar
  argument <- deparse(substitute(param_list))
  if (is.null(param_list[[element]])) {
    if (length(default) == 1) {
      output <- matrix(default, dimensions[1], dimensions[2])
    } else {
      output <- default
    }
  } else if (length(param_list[[element]]) == 1) {
      output <- matrix(param_list[[element]], dimensions[1], dimensions[2])
  } else if (sum(dim(param_list[[element]]) == dimensions) == 2) {
    output <- param_list[[element]]
  } else {
    stop(sprintf("element '%s' of '%s' must be of length 1 or dimension c(%i, %i)",
                 element, argument, dimensions[[1]], dimensions[[2]]))
  }
  return(output)
}


check_param_matdiag <- function (param_list, element, dimension, default) {
  # It accepts matrices, vectors and scalar: only for square matrices
  argument <- deparse(substitute(param_list))
  if (is.null(param_list[[element]])) {
    output <- default
  } else if (length(param_list[[element]]) == 1) {
    output <- diag(as.numeric(param_list[[element]]), dimension, dimension)
  } else if (length(param_list[[element]]) == dimension) {
    output <- diag(as.numeric(param_list[[element]]))
  } else if (sum(dim(param_list[[element]]) == rep(dimension, 2)) == 2) {
    output <- param_list[[element]]
  } else {
    stop(sprintf("element '%s' of argument '%s' must be of length ", element, argument),
         sprintf("1 or %i, or dimension c(%i, %i)", dimension, dimension, dimension))
  }
  return(output)
}
