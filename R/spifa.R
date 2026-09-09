#' @title Bayesian Spatial Item Factor Analysis
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
#' \strong{Parameter glossary.} \code{priors}/\code{constraints} use
#' descriptive names; the fitted model's sampled output (from
#' \code{\link{as.list.spifa}}, \code{\link{summary.spifa}}, and
#' \code{\link{as_tibble.spifa}}) instead uses short internal names
#' matching the underlying model notation. The two are deliberately
#' different vocabularies (what you configure vs. what the sampler
#' produced) -- this table maps between them:
#' \tabular{lll}{
#' \code{priors}/\code{constraints} name \tab output block name \tab meaning \cr
#' \code{easiness} \tab \code{c} \tab item easiness (intercept) \cr
#' \code{discrimination} \tab \code{A} \tab item-factor discrimination (loading) matrix \cr
#' \code{effect} \tab \code{B} \tab predictor effect on the latent factors \cr
#' \code{corr} \tab \code{Corr}, \code{Chol} \tab residual correlation matrix, and its Cholesky factor \cr
#' \code{sd} \tab \emph{(fixed, not sampled)} \tab residual standard deviation \cr
#' \code{loading} \tab \code{T} \tab multivariate Gaussian process loading matrix \cr
#' \code{range} \tab \code{phi} \tab multivariate Gaussian process spatial range \cr
#' \emph{(not user-set)} \tab \code{Theta} \tab latent abilities \cr
#' \emph{(not user-set)} \tab \code{Z} \tab augmented latent response
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
#' is used as the spatial coordinates and spatial Gaussian processes are added
#' to the model.
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
#' \code{loading}/\code{range}/\code{corr} still adapt through
#' them), but they are never stored, so \code{niter} counts only the
#' iterations that end up in the returned samples. \code{0} by default (no
#' iterations discarded during fitting -- the previous behaviour). Prefer
#' this over discarding a prefix of the samples afterwards (e.g. via
#' \code{summary(..., burnin = )}): iterations dropped here were never
#' stored, so they don't cost memory or thinning-index arithmetic, and the
#' adaptive proposals get to keep converging across the burnin/niter
#' boundary rather than someone accidentally analysing them as if they were
#' post-adaptation draws.
#' @param standardize Logical; if \code{TRUE} (default), the stored posterior
#' draws are rescaled after fitting so the latent factors have unit variance:
#' \code{theta} is divided by its posterior SD per factor, and
#' \code{discrimination}, \code{effect}, \code{loading}, and the residual SD
#' (\code{constraints$sd}) are compensated by the same factor so the
#' fitted response probabilities are unchanged. This only applies to models
#' with a spatial Gaussian process and/or predictor effects on the latent
#' factors (\code{cifa_pred}/\code{spifa}/\code{spifa_pred}; ignored for
#' \code{eifa}/\code{cifa}, where the residual SD is the only source of
#' \code{theta}'s variance and there's nothing to normalize against). Those
#' model types have a multiplicative scale non-identifiability between
#' \code{theta} and \code{discrimination}/the GP variance/the predictor
#' effect -- an equally good fit can shrink one and inflate the other --
#' so leaving \code{TRUE} keeps draws on an interpretable, comparable scale.
#' Set to \code{FALSE} to keep the raw, unscaled posterior, e.g. when
#' checking recovery of known simulated parameters (see
#' \code{dev/simulated/analyze-ipixuna.R}), where an extra rescale would
#' make draws harder to compare directly against the true simulated values.
#' No predictors are standardized by this argument -- despite the name, it
#' does not touch \code{formula}'s right-hand side at all.
#' @param constraints Named list of constraints associated to the factor model. Accepted
#' names are `discrimination`, `loading`, and `sd`. The restrictions on the
#' discrimination paramater should be placed in the element `discrimination` with same
#' dimensions as the discrimination matrix (nitems x nfactors). A value of 0 indicates that
#' the link betwen the item and the factor is disabled and 1 indicates that it remains
#' active and the coefficient associated will be estimated. The restrictions for the
#' multivariate Gaussian process loading matrix should be placed in the element `loading`
#' with dimensions nfactors x ngp, such as a value of 0 indicates a link disconnected between
#' the factor and the (independent) GP while 1 indicates that it remains active. The restrictions with
#' respect to the standard deviation of the latent factors' residual term should be placed
#' in the element `sd`, which should be a vector (length nfactors) providing the fixed
#' values for that standard deviation (paired with `priors$corr`, together they
#' parameterize the residual covariance). If the model includes
#' predictors or a Gaussian process, it is recomended to be lower than 1.
#'
#' @param priors Named list of initial values and prior hyperparameters, one
#' element per parameter block: `easiness`, `discrimination`, `effect`
#' (predictor effect on the latent factors), `corr` (correlation of the
#' latent factors' residual term, paired with `constraints$sd`),
#' `loading` (multivariate Gaussian process loading matrix, paired with
#' `constraints$loading`), and `range` (multivariate Gaussian process
#' scale parameters). Each
#' element (except `corr`) accepts `initial`, `mean`, and `sd`;
#' `corr` accepts `initial` and `eta` (the LKJ prior shape parameter).
#' See the parameter glossary above for how these names map to the
#' fitted model's sampled output.
#'
#' @param adaptive Named list of hyperparameters associated with the adaptive sampling.
#' The adaptive sampling is done jointly for the `correlation` parameters, `standard
#' deviation of the gps` and `scale parameter of the gps`. The matrix `Sigma` can be
#' provided as the full covariance matrix of these parameters for the proposal
#' distribution. Otherwise, part of this matrix can be provided by using the elements
#' `Sigma`, `Sigma_corr`, `Sigma_loading` and `Sigma_range`. Additional elements are `scale`,
#' `C`, `alpha` and `accep_prob` which are hyperparameters of the adaptive sampling
#' proposed in Andrieu and Thomas (2008).
#'
#' @param execute Logical value to run sampler or not. TRUE by default.
#'
#' @return
#' An object of (informal) class \code{spifa.list}: a named list of MCMC
#' sample matrices (one entry per parameter block, e.g. \code{c}, \code{a},
#' \code{theta}, \code{corr}, \code{betas}, ...), with an attribute
#' \code{"spifa_args"} recording the data and options used to fit the model
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
#'   constraints = list(discrimination = L_a, sd = rep(0.5, nfactors)))
#' summary(samples, select = c("c", "A"))
#'
#' @export
spifa <- function(formula, data, nfactors, ngp = nfactors,
    niter = 100, thin = 1, burnin = 0, standardize = TRUE,
    constraints = list(discrimination = NULL, loading = NULL, sd = rep(1, nfactors)),
    priors = list(
      easiness = list(initial = NULL, mean = NULL, sd = NULL),
      discrimination = list(initial = NULL, mean = NULL, sd = NULL),
      effect = list(initial = NULL, mean = NULL, sd = NULL),
      corr = list(initial = NULL, eta = 1.5),
      loading = list(initial = NULL, mean = NULL, sd = NULL),
      range = list(initial = NULL, mean = NULL, sd = NULL)),
    adaptive = list(Sigma = NULL, Sigma_corr = NULL, Sigma_loading = NULL, Sigma_range = NULL,
                    scale = 1, C = 0.7, alpha = 0.8, accep_prob = 0.234),
    execute = TRUE) {

  # Trim niter to the last stored iteration
  if (niter > 0) niter <- thin * ((niter - 1) %/% thin) + 1

  # Dimensions, items and predictors
  mf <- model.frame(formula, data)
  response <- model.response(mf)
  if (!is.matrix(response)) stop("The left-hand side of 'formula' must be a matrix")
  nobs <- nrow(response)
  nitems <- ncol(response)

  predictors_terms <- delete.response(terms(mf))
  attr(predictors_terms, "intercept") <- 0
  predictors <- model.matrix(predictors_terms, mf)
  npred <- ncol(predictors)

  # Coordinates and ngp
  if (inherits(data, "sf") && ngp > 0) {
    coordinates <- sf::st_geometry(data)
  } else {
    coordinates <- NULL
    ngp <- 0
  }

  # Restrictions for discrimination parameters and Gaussian process loadings
  constrain_L_explo <- matrix(NA, nitems, nfactors)
  constrain_L_explo <- lower.tri(constrain_L_explo, diag = TRUE) * 1
  constrain_L <- check_param_mat(constraints, "discrimination", c(nitems, nfactors), constrain_L_explo)
  constrain_T <- check_param_mat(constraints, "loading", c(nfactors, ngp), diag(1, nfactors, ngp))

  # Sizes
  nsigmas <- sum(constrain_T)
  ncorr <- nfactors * (nfactors - 1) / 2

  # Model type: EIFA, CIFA, CIFA_PRED, SPIFA, SPIFA_PRED
  if (!is.null(coordinates)) {
    if (npred > 0) {
      model_type <- "spifa_pred"
      constrain_V_sd <- check_param_vec(constraints, "sd", nfactors, 0.2)
    } else {
      model_type <- "spifa"
      constrain_V_sd <- check_param_vec(constraints, "sd", nfactors, 0.2)
    }
  } else if (npred > 0) {
    model_type <- "cifa_pred"
    constrain_V_sd <- check_param_vec(constraints, "sd", nfactors, 0.3)
  } else if (all(constrain_L == constrain_L_explo)) {
    model_type <- "eifa"
    constrain_V_sd <- check_param_vec(constraints, "sd", nfactors, 1)
  } else {
    model_type <- "cifa"
    constrain_V_sd <- check_param_vec(constraints, "sd", nfactors, 1)
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
  adap_Sigma_R <- check_param_matdiag(adaptive, "Sigma_corr", ncorr, diag(ncorr) * 0.001)
  adap_Sigma_gp_sd <-
    check_param_matdiag(adaptive, "Sigma_loading", nsigmas, diag(nsigmas) * 0.001)
  adap_Sigma_gp_phi <-
    check_param_matdiag(adaptive, "Sigma_range", ngp, diag(ngp) * 0.001)
  adap_scale <- ifelse(is.null(adaptive$scale), 1, adaptive$scale)
  adap_C <- ifelse(is.null(adaptive$C), 0.7, adaptive$C)
  adap_alpha <- ifelse(is.null(adaptive$alpha), 0.8, adaptive$alpha)
  adap_accep_prob <- ifelse(is.null(adaptive$accep_prob), 0.234, adaptive$accep_prob)

  # Create general sigma proposal in order: gp_sd, gp_phi, corr_free
  if (is.null(coordinates)) {
    if (is.null(adaptive$Sigma)) {
      adap_Sigma <- adap_Sigma_R
    } else if (sum(dim(adaptive$Sigma) == c(ncorr, ncorr)) == 2) {
      adap_Sigma <- adaptive$Sigma
    }
  } else {
    if (is.null(adaptive$Sigma)) {
      adap_Sigma <- matrix(0, nsigmas + ngp + ncorr, nsigmas + ngp + ncorr)
      adap_Sigma[seq_len(nsigmas), seq_len(nsigmas)] <- adap_Sigma_gp_sd
      adap_Sigma[nsigmas + seq_len(ngp), nsigmas + seq_len(ngp)] <- adap_Sigma_gp_phi
      adap_Sigma[nsigmas + ngp + seq_len(ncorr), nsigmas + ngp + seq_len(ncorr)] <- adap_Sigma_R
    } else if (sum(dim(adaptive$Sigma) == rep(nsigmas + ngp + ncorr, 2)) == 2) {
      adap_Sigma <- adaptive$Sigma
    }
  }

  # Optional arguments for parameter of residual correlation R
  if (is.null(priors$corr$initial)) {
    R_initial <- diag(nfactors)
  } else if (sum(dim(priors$corr$initial) == c(nfactors, nfactors)) == 2) {
    if (all(diag(priors$corr$initial) == 1)) {
      R_initial <- priors$corr$initial
    } else {
      stop("'initial' of 'priors$corr' argument is not a correlation matrix")
    }
  } else {
    stop("'initial' of 'priors$corr' argument must be of dimension ",
         sprintf("c(%i, %i)", nfactors, nfactors))
  }

  R_prior_eta <- ifelse(is.null(priors$corr$eta), 1, priors$corr$eta)

  # Optional arguments for parameter of fixed effects (Beta)
  B_prior_mean <- check_param_mat2(priors$effect, "mean", c(npred, nfactors), 0)
  B_prior_sd <- check_param_mat2(priors$effect, "sd", c(npred, nfactors), 1)
  B_initial <- check_param_mat2(priors$effect, "initial", c(npred, nfactors), B_prior_mean)

  # Optional arguments for GP standard deviations and  scale parameters
  sigmas_gp_mean <- check_param_vec(priors$loading, "mean", nsigmas, 0.6)
  sigmas_gp_sd <- check_param_vec(priors$loading, "sd", nsigmas, 0.2)
  sigmas_gp_initial <- check_param_vec(priors$loading, "initial", nsigmas, sigmas_gp_mean)
  phi_gp_mean <- check_param_vec(priors$range, "mean", ngp, 0.05)
  phi_gp_sd <- check_param_vec(priors$range, "sd", ngp, 0.2)
  phi_gp_initial <- check_param_vec(priors$range, "initial", ngp, phi_gp_mean)

  # Compute distances as a matrix
  if (is.null(coordinates)) {
    distances <- matrix(nrow = 0, ncol = 0)
  } else {
    distances <- matrix(as.numeric(sf::st_distance(coordinates)), nobs)
  }

  # List of options to call c++ spifa function
  spifa_args <- list(
    response = as.numeric(response), predictors = predictors, distances = distances,
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
    samples <- do.call(spifa_cpp, spifa_args)
    spifa_args$constrain_V_sd <- attr(samples, "V_sd")
    samples <- do.call(cbind, samples) |> posterior::as_draws_array()
  } else {
    samples <- list()
  }

  # Add attributes
  attr(samples, "spifa_args") <- spifa_args
  attr(samples, "coordinates") <- coordinates
  class(samples) <- unique(c("spifa", class(samples)))
  return(samples)
}

check_param_vec <- function (param_list, element, dimension, default) {
  # Only for vectors
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
