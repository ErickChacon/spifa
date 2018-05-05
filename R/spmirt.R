

spmirt <- function (response, predictors = NULL, coordinates = NULL,
    nobs, nitems, nfactors, niter = 1000, thin = 10,
    constrains = list(A = NULL, W = NULL, V_sd = rep(1, nfactors)),
    adaptive = list(Sigma = NULL, Sigma_R = NULL, Sigma_gp_sd = NULL, Sigma_gp_phi = NULL,
                    scale = 1, C = 0.7, alpha = 0.8, accep_prob = 0.234),
    c_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
    A_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
    R_opt = list(initial = NULL, prior_eta = 1.5),
    B_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
    sigmas_gp_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
    phi_gp_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL)) {


  # Restrictions for discrimination parameters (A) and Gaussian processes (W)

  constrain_L_explo <- matrix(NA, nitems, nfactors)
  constrain_L_explo <- lower.tri(constrain_L_explo, diag = TRUE) * 1
  constrain_L <- check_param_mat(constrains, "A", c(nitems, nfactors), constrain_L_explo)
  constrain_T <- check_param_mat(constrains, "W", c(nfactors, nfactors), diag(nfactors))

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

  n_corr <- nfactors * (nfactors - 1) / 2
  adap_Sigma_R <- check_param_matdiag(adaptive, "Sigma_R", n_corr, diag(n_corr) * 0.001)
  adap_Sigma_gp_sd <-
    check_param_matdiag(adaptive, "Sigma_gp_sd", nfactors, diag(nfactors) * 0.001)
  adap_Sigma_gp_phi <-
    check_param_matdiag(adaptive, "Sigma_gp_phi", nfactors, diag(nfactors) * 0.001)
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
      adap_Sigma <- matrix(0, 2 * nfactors + n_corr, 2 * nfactors + n_corr)
      adap_Sigma[1:nfactors, 1:nfactors] <- adap_Sigma_gp_sd
      adap_Sigma[nfactors + 1:nfactors, nfactors + 1:nfactors] <- adap_Sigma_gp_phi
      adap_Sigma[2*nfactors + 1:n_corr, 2*nfactors + 1:n_corr] <- adap_Sigma_R
    } else if (sum(dim(adaptive$Sigma) == rep(2*nfactors + n_corr, 2)) == 2) {
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

  R_prior_eta <- ifelse(is.null(R_opt$prior_eta), 1.5, R_opt$prior_eta)

  # Optional arguments for parameter of fixed effects (Beta)

  if (is.null(predictors)) {
    B_prior_mean <- matrix(NA)
    B_prior_sd <- matrix(NA)
    B_initial <- matrix(NA)
  } else {
    npred <- ncol(predictors)
    B_prior_mean <- check_param_mat2(B_opt, "prior_mean", c(npred, nfactors), 0)
    B_prior_sd <- check_param_mat2(B_opt, "prior_sd", c(npred, nfactors), 1)
    B_initial <- check_param_mat2(B_opt, "initial", c(npred, nfactors), B_prior_mean)
  }

  # Optional arguments for GP standard deviations and  scale parameters

  if (is.null(coordinates)) {
    sigmas_gp_mean <- NA
    sigmas_gp_sd <- NA
    sigmas_gp_initial <- NA
    phi_gp_mean <- NA
    phi_gp_sd <- NA
    phi_gp_initial <- NA
  } else {
    sigmas_gp_mean <- check_param_vec(sigmas_gp_opt, "prior_mean", nfactors, 0.6)
    sigmas_gp_sd <- check_param_vec(sigmas_gp_opt, "prior_sd", nfactors, 0.2)
    sigmas_gp_initial <- check_param_vec(sigmas_gp_opt, "initial", nfactors, sigmas_gp_mean)
    phi_gp_mean <- check_param_vec(phi_gp_opt, "prior_mean", nfactors, 0.05)
    phi_gp_sd <- check_param_vec(phi_gp_opt, "prior_sd", nfactors, 0.2)
    phi_gp_initial <- check_param_vec(phi_gp_opt, "initial", nfactors, phi_gp_mean)
  }

  # Compute predictors and distances as matrices

  if (is.null(predictors))  predictors <- matrix(NA)
  if (is.null(coordinates))  coordinates <- matrix(NA)

  # Execute model calling c++ spmirt function
  samples <- spmirt_cpp(
    response = response, predictors = predictors, distances = coordinates,
    nobs = nobs, nitems = nitems, nfactors = nfactors, niter = niter, thin = thin,
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

  attr(samples, "model") <- model_type

  # samples <- list(
  #   response = response, predictors = predictors, distances = coordinates,
  #   nobs = nobs, nitems = nitems, nfactors = nfactors, niter = niter, thin = thin,
  #   constrain_L = constrain_L, constrain_T = constrain_T, constrain_V_sd = constrain_V_sd,
  #   adap_Sigma = adap_Sigma, adap_scale = adap_scale, adap_C = adap_C,
  #   adap_alpha = adap_alpha, adap_accep_prob = adap_accep_prob,
  #   c_initial = c_initial, c_prior_mean = c_prior_mean, c_prior_sd = c_prior_sd,
  #   A_initial = A_initial, A_prior_mean = A_prior_mean, A_prior_sd = A_prior_sd,
  #   R_initial = R_initial, R_prior_eta = R_prior_eta,
  #   B_initial = B_initial, B_prior_mean = B_prior_mean, B_prior_sd = B_prior_sd,
  #   sigmas_gp_initial = sigmas_gp_initial, sigmas_gp_mean = sigmas_gp_mean,
  #   sigmas_gp_sd = sigmas_gp_sd,
  #   phi_gp_initial = phi_gp_initial, phi_gp_mean = phi_gp_mean, phi_gp_sd = phi_gp_sd,
  #   model_type = model_type
  #   )
  #
  #

  return(samples)
}


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
