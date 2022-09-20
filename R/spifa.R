
#' @title Spatial Multidimensional Item Response Model with Predictors
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param constrains Named list of constrains associated to the factor model. Accepted
#' names are `A`, `W`, and `V_sd`. The restrictions in the discrimination paramater should
#' be place in the element `A` with same dimensions as the discrimination matrix (nitems x
#' nfactors). A value of 0 indicates that link betwen the item and the factor is disabled
#' and 1 indicates that it remain active and the coefficient associated will be estimated.
#' The restrictions for the gaussian processes should be place in the element `W` with
#' dimensions nfactors x ngps, such as a value of 0 indicates a link disconnected between
#' the factor and the gp while 1 indicates that it remains active. The restrictions with
#' respect to the standard deviation of the error terms inside the latent factors should
#' be placed in the element `V_sd`, which should be a vector (length nfactors) providing
#' the fixed values for the error standard deviation. If the model includes predictors or
#' gaussian processes, it is recomended to be lower than 1.
#'
#' @param adaptive Named list of hyperparameters associated with the adaptive sampling.
#' The adaptive sampling is done jointly for the `correlation` parameters, `standard
#' deviation of the gps` and `scale parameter of the gps`. The matrix `Sigma` can be
#' provided as the full covariance matrix of these parameters for the proposal
#' distribution. Otherwise, part of this matrix can be provided by using the elements
#' `Sigma`, `Sigma_R`, `Sigma_gp_sd` and `Sigma_gp_phi`. Additional elements are `scale`,
#' `C`, `alpha` and `accep_prob` which are hyperparameters of the adaptive sampling
#' proposed in Andrieu and Thomas (2008).
#'
#' @param C_opt Named list of initial values and hyperparameters for the easiness
#' parameters. The initial value is provided in the element `initial`, and the prior meand
#' and standard deviation are provided in the elements `prior_mean` and `prior_sd`
#' respectively.
#'
#' @param A_opt Same as C_opt but for the discrimination parameters.
#'
#' @param R_opt Same as C_opt but for the correlation parameters. This list only accepts
#' `initial` value and `prior_eta` associated to the LKJ prior.
#'
#' @param B_opt Same as C_opt but for the regression parameters.
#'
#' @param sigmas_gp_opt Same as C_opt but for the standard deviation of the gaussian
#' processes.
#'
#' @param phi_gp_opt Same as C_opt but for the scale parameter of the gaussian processes.
#'
#' @param execute Logical value to run sampler or not. TRUE by default.
#'
#' @return return.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#'
#' bla
#'
#' @export
spifa <- function (responses, pred_formula = NULL, coords = NULL, data = NULL,
    nfactors, ngp = nfactors, niter = 1000, thin = 10, standardize = TRUE,
    constrains = list(A = NULL, W = NULL, V_sd = rep(1, nfactors)),
    adaptive = list(Sigma = NULL, Sigma_R = NULL, Sigma_gp_sd = NULL, Sigma_gp_phi = NULL,
                    scale = 1, C = 0.7, alpha = 0.8, accep_prob = 0.234),
    c_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
    A_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
    R_opt = list(initial = NULL, prior_eta = 1.5),
    B_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
    sigmas_gp_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
    phi_gp_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
    execute = TRUE) {

  # Names to substitute and variables of data
  responses <- substitute(responses)
  coords <- substitute(coords)
  vars <- setNames(as.list(seq_along(data)), names(data))
  vars_sfc <- sapply(data, function (x) inherits(x, "sfc"))

  # Coordinates
  if (!is.null(coords)) {
    coords_pos <- eval(coords, vars, parent.frame())
    if (sum(vars_sfc) > 0 && length(coords_pos) == 1) {
      coordinates <- data[[coords_pos]]
    } else {
      coordinates <- data[, coords_pos, drop = FALSE]
    }
  } else {
    coordinates <- NULL
  }

  # Responses
  responses_pos <- eval(responses, vars, parent.frame())
  response <- data[, responses_pos, drop = FALSE]
  if (inherits(response, "sf")) response <- sf::st_set_geometry(response, NULL)
  response <- as.matrix(response)
  nobs <- nrow(response)
  nitems <- ncol(response)
  response <- as.numeric(response)

  # Predictors
  if (is.null(pred_formula)) {
    predictors <- NULL
  } else {
    pred_formula <- update(pred_formula,  ~ . - 1)
    predictors <- model.matrix(pred_formula, data)
  }

  # Restrictions for discrimination parameters (A) and Gaussian processes (W)
  constrain_L_explo <- matrix(NA, nitems, nfactors)
  constrain_L_explo <- lower.tri(constrain_L_explo, diag = TRUE) * 1
  constrain_L <- check_param_mat(constrains, "A", c(nitems, nfactors), constrain_L_explo)
  constrain_T <- check_param_mat(constrains, "W", c(nfactors, ngp), diag(1, nfactors, ngp))

  # Sizes
  nsigmas <- sum(constrain_T)
  n_corr <- nfactors * (nfactors - 1) / 2

  # Detect type of model to be fitted: EIFA, CIFA, CIFA_PRED, SPIFA, SPIFA_PRED
  if (!is.null(coordinates)) {
    if (!is.null(predictors)) {
      model_type = "spifa_pred"
      constrain_V_sd <- check_param_vec(constrains, "V_sd", nfactors, 0.2)
    } else {
      model_type = "spifa"
      constrain_V_sd <- check_param_vec(constrains, "V_sd", nfactors, 0.2)
    }
  } else if (!is.null(predictors)) {
    model_type = "cifa_pred"
    constrain_V_sd <- check_param_vec(constrains, "V_sd", nfactors, 0.3)
  } else if (all(constrain_L == constrain_L_explo)) {
    model_type = "eifa"
    constrain_V_sd <- check_param_vec(constrains, "V_sd", nfactors, 1)
  } else {
    model_type = "cifa"
    constrain_V_sd <- check_param_vec(constrains, "V_sd", nfactors, 1)
  }

  # Optional arguments for difficulty parameters (c)
  c_prior_mean <- check_param_vec(c_opt, "prior_mean", nitems, 0)
  c_prior_sd <- check_param_vec(c_opt, "prior_sd", nitems, 1)
  c_initial <- check_param_vec(c_opt, "initial", nitems,
                               rnorm(nitems, c_prior_mean, c_prior_sd))

  # Optional arguments for discrimination parameters (A)
  A_prior_mean <-
    check_param_mat2(A_opt, "prior_mean", c(nitems, nfactors), diag(1, nitems, nfactors))
  A_prior_sd <-
    check_param_mat2(A_opt, "prior_sd", c(nitems, nfactors), 1-diag(0.55, nitems, nfactors))
  A_initial <-
    check_param_mat2(A_opt, "initial", c(nitems, nfactors), A_prior_mean)

  # Adaptive Metropolis-Hastings arguments for proposed covariance matrix
  adap_Sigma_R <- check_param_matdiag(adaptive, "Sigma_R", n_corr, diag(n_corr) * 0.001)
  adap_Sigma_gp_sd <-
    check_param_matdiag(adaptive, "Sigma_gp_sd", nsigmas, diag(nsigmas) * 0.001)
  adap_Sigma_gp_phi <-
    check_param_matdiag(adaptive, "Sigma_gp_phi", ngp, diag(ngp) * 0.001)
  adap_scale <- ifelse(is.null(adaptive$scale), 1, adaptive$scale)
  adap_C <- ifelse(is.null(adaptive$C), 0.7, adaptive$C)
  adap_alpha <- ifelse(is.null(adaptive$alpha), 0.8, adaptive$alpha)
  adap_accep_prob <- ifelse(is.null(adaptive$accep_prob), 0.234, adaptive$accep_prob)

  # Create general sigma proposal in order: gp_sd, gp_phi, corr_free
  if (is.null(coordinates)) {
    if (is.null(adaptive$Sigma)) {
      adap_Sigma <- adap_Sigma_R
    } else if (sum(dim(adaptive$Sigma) == c(n_corr, n_corr)) == 2) {
      adap_Sigma <- adaptive$Sigma
    }
  } else {
    if (is.null(adaptive$Sigma)) {
      adap_Sigma <- matrix(0, nsigmas + ngp + n_corr, nsigmas + ngp + n_corr)
      adap_Sigma[1:nsigmas, 1:nsigmas] <- adap_Sigma_gp_sd
      adap_Sigma[nsigmas + 1:ngp, nsigmas + 1:ngp] <- adap_Sigma_gp_phi
      adap_Sigma[nsigmas + ngp + 1:n_corr, nsigmas + ngp + 1:n_corr] <- adap_Sigma_R
    } else if (sum(dim(adaptive$Sigma) == rep(nsigmas + ngp + n_corr, 2)) == 2) {
      adap_Sigma <- adaptive$Sigma
    }
  }

  # Optional arguments for parameter of residual correlation R
  if (is.null(R_opt$initial)) {
    R_initial <- diag(nfactors)
  } else if (sum(dim(R_opt$initial) == c(nfactors, nfactors)) == 2) {
    if (all(diag(R_opt$initial) == 1)) {
      R_initial <- R_opt$initial
    } else {
      stop("element 'initial' of 'R_opt' argument is not a correlation matrix")
    }
  } else {
    stop("element 'initial' of 'R_opt' argument must be of dimension ",
         sprintf("c(%i, %i)", nfactors, nfactors))
  }

  R_prior_eta <- ifelse(is.null(R_opt$prior_eta), 1, R_opt$prior_eta)

  # Optional arguments for parameter of fixed effects (Beta)
  if (is.null(predictors)) {
    B_prior_mean <- matrix(NA, 1, nfactors)
    B_prior_sd <- matrix(NA, 1, nfactors)
    B_initial <- matrix(NA, 1, nfactors)
  } else {
    npred <- ncol(predictors)
    B_prior_mean <- check_param_mat2(B_opt, "prior_mean", c(npred, nfactors), 0)
    B_prior_sd <- check_param_mat2(B_opt, "prior_sd", c(npred, nfactors), 1)
    B_initial <- check_param_mat2(B_opt, "initial", c(npred, nfactors), B_prior_mean)
  }

  # Optional arguments for GP standard deviations and  scale parameters
  if (is.null(coordinates)) {
    sigmas_gp_mean <- rep(NA, nsigmas)
    sigmas_gp_sd <- rep(NA, nsigmas)
    sigmas_gp_initial <- rep(NA, nsigmas)
    phi_gp_mean <- rep(NA, ngp)
    phi_gp_sd <- rep(NA, ngp)
    phi_gp_initial <- rep(NA, ngp)
  } else {
    sigmas_gp_mean <- check_param_vec(sigmas_gp_opt, "prior_mean", nsigmas, 0.6)
    sigmas_gp_sd <- check_param_vec(sigmas_gp_opt, "prior_sd", nsigmas, 0.2)
    sigmas_gp_initial <- check_param_vec(sigmas_gp_opt, "initial", nsigmas, sigmas_gp_mean)
    phi_gp_mean <- check_param_vec(phi_gp_opt, "prior_mean", ngp, 0.05)
    phi_gp_sd <- check_param_vec(phi_gp_opt, "prior_sd", ngp, 0.2)
    phi_gp_initial <- check_param_vec(phi_gp_opt, "initial", ngp, phi_gp_mean)
  }

  # Compute predictors and distances as matrices
  if (is.null(predictors))  predictors <- matrix(NA)
  if (is.null(coordinates)) {
    distances <- matrix(NA)
  } else if (inherits(coordinates, "sfc")) {
    distances <- sf::st_transform(coordinates, crs = 3857) %>%
      sf::st_coordinates()  %>%
      dist()  %>%
      as.matrix()
    # distances <- unclass(sf::st_distance(coordinates))
    # distances <- (distances + t(distances)) / 2
  } else {
    distances <- as.matrix(dist(coordinates))
  }

  # List of options to call c++ spifa function
  model_info <- list(
    response = response, predictors = predictors, distances = distances,
    nobs = nobs, nitems = nitems, nfactors = nfactors, ngp = ngp,
    niter = niter, thin = thin, standardize = standardize,
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
  } else {
    samples <- list()
  }

  # Add model_info to MCMC samples
  attr(samples, "model_info") <- model_info

  return(samples)
}
