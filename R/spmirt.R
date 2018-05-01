
spmirt <- function (response, predictors = NULL, coordinates = NULL,
    nobs, nitems, nfactors, niter = 100,
    constrains = list(A = NULL, W = NULL),
    adaptive = list(Sigma_R = NULL, Sigma_gp_sd = NULL, Sigma_gp_phi = NULL,
                    scale = 1, C = 0.7, alpha = 0.8, accep_prob = 0.234),
    c_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
    A_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
    R_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
    B_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
    sigmas_gp_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
    phi_gp_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL)
    ) {

  # Optional arguments for parameter c, checking adequate sizes

  if (is.null(c_opt$prior_mean)) {
    c_prior_mean <- rep(0, nitems)
  } else {
    if (length(c_opt$prior_mean) == 1) {
      c_prior_mean <- rep(c_opt$prior_mean, nitems)
    } else if (length(c_opt$prior_mean) == nitems) {
      c_prior_mean <- c_opt$prior_mean
    } else {
      stop("element 'prior_mean' of 'c_opt' must be of length 1 or", nitems)
    }
  }

  if (is.null(c_opt$prior_sd)) {
    c_prior_sd <- 1
  } else {
    if (length(c_opt$prior_sd) == 1) {
      c_prior_sd <- rep(c_opt$prior_sd, nitems)
    } else if (length(c_opt$prior_sd) == nitems) {
      c_prior_sd <- c_opt$prior_sd
    } else {
      stop("element 'prior_sd' of 'c_opt' must be of length 1 or", nitems)
    }
  }

  if (is.null(c_opt$initial)) {
    c_initial <- rnorm(nitems, c_prior_mean, c_prior_sd)
  } else {
    if (length(c_opt$initial) == 1) {
      c_initial <- rep(c_opt$initial, nitems)
    } else if (length(c_opt$initial) == nitems) {
      c_initial <- c_opt$initial
    } else {
      stop("element 'initial' of 'c_opt' must be of length 1 or ", nitems)
    }
  }

  # Optional arguments for parameter A, checking adequate sizes

  if (is.null(A_opt$prior_mean)) {
    A_prior_mean <- matrix(0, nitems, nfactors)
    diag(A_prior_mean) <- 1
  } else {
    if (length(A_opt$prior_mean) == 1) {
      A_prior_mean <- matrix(A_opt$prior_mean, nitems, nfactors)
    } else if (sum(dim(A_opt$prior_mean) == c(nitems, nfactors))) {
      A_prior_mean <- A_opt$prior_mean
    } else {
      stop("element 'prior_mean' of 'A_opt' must be of length 1 or",
           sprintf(" of dimension c(%i, %i)", nitems, nfactors))
    }
  }

  if (is.null(A_opt$prior_sd)) {
    A_prior_sd <- matrix(1, nitems, nfactors)
    if (is.null(A_opt$prior_sd)) {
      diag(A_prior_sd) <- 0.45
    }
  } else {
    if (length(A_opt$prior_sd) == 1) {
      A_prior_sd <- matrix(A_opt$prior_sd, nitems, nfactors)
    } else if (sum(dim(A_opt$prior_mean) == c(nitems, nfactors))) {
      A_prior_sd <- A_opt$prior_sd
    } else {
      stop("element 'prior_sd' of 'A_opt' must be of length 1 or",
           sprintf(" of dimension c(%i, %i)", nitems, nfactors))
    }
  }

  if (is.null(A_opt$initial)) {
      A_initial <- A_prior_mean
  } else {
    if (length(A_opt$initial) == 1) {
      A_initial <- matrix(A_opt$initial, nitems, nfactors)
    } else if (sum(dim(A_opt$initial) == c(nitems, nfactors))) {
      A_initial <- A_opt$initial
    } else {
      stop("element 'initial' of 'A_opt' must be of length 1 or",
           sprintf(" of dimension c(%i, %i)", nitems, nfactors))
    }
  }

  # Restrictions arguments for parameter A

  if (is.null(constrains$A)) {
    A_aux <- matrix(NA, nitems, nfactors)
    constrain_L <- lower.tri(A_aux, diag = TRUE) * 1
  } else {
    if (sum(dim(constrains$A) == c(nitems, nfactors))) {
      constrain_L <- constrains$A
    } else {
      stop("element 'A' of 'constrains' must be of dimension ",
           sprintf(" c(%i, %i)", nitems, nfactors))
    }
  }

  # Restrictions arguments for Gaussian Processes W(s)

  # Adaptive covariance matrix of proposal distribution

  n_corr = nfactors * (nfactors - 1) / 2
  if (is.null(adaptive$Sigma_R)) {
    adap_Sigma_R = diag(n_corr) * 0.001
  } else {
    if (length(adaptive$Sigma_R) == 1) {
      adap_Sigma_R <- diag(adaptive$Sigma_R, n_corr, n_corr)
    } else if (length(adaptive$Sigma_R) == n_corr) {
      adap_Sigma_R = diag(as.numeric(adaptive$Sigma_R))
    } else if (sum(dim(adaptive$Sigma_R) == c(n_corr, n_corr))) {
      adap_Sigma_R = adaptive$Sigma_R
    } else {
      stop("element 'Sigma_R' of 'adaptive' argument must be of length 1, ",
           nfactors, sprintf(" or of dimension c(%i, %i)", n_corr, n_corr))
    }
  }

  if (is.null(adaptive$Sigma_gp_sd)) {
    adap_Sigma_gp_sd = diag(nfactors) * 0.001
  } else {
    if (length(adaptive$Sigma_gp_sd) == 1) {
      adap_Sigma_gp_sd <- diag(adaptive$Sigma_gp_sd, nfactors, nfactors)
    } else if (length(adaptive$Sigma_gp_sd) == nfactors) {
      adap_Sigma_gp_sd = diag(as.numeric(adaptive$Sigma_gp_sd))
    } else if (sum(dim(adaptive$Sigma_gp_sd) == c(nfactors, nfactors))) {
      adap_Sigma_gp_sd = adaptive$Sigma_gp_sd
    } else {
      stop("element 'Sigma_gp_sd' of 'adaptive' argument must be of length 1, ",
           nfactors, sprintf(" or of dimension c(%i, %i)", nfactors, nfactors))
    }
  }

  if (is.null(adaptive$Sigma_gp_phi)) {
    adap_Sigma_gp_phi = diag(nfactors) * 0.001
  } else {
    if (length(adaptive$Sigma_gp_phi) == 1) {
      adap_Sigma_gp_phi <- diag(adaptive$Sigma_gp_phi, nfactors, nfactors)
    } else if (length(adaptive$Sigma_gp_phi) == nfactors) {
      adap_Sigma_gp_phi = diag(as.numeric(adaptive$Sigma_gp_phi))
    } else if (sum(dim(adaptive$Sigma_gp_phi) == c(nfactors, nfactors))) {
      adap_Sigma_gp_phi = adaptive$Sigma_gp_phi
    } else {
      stop("element 'Sigma_gp_phi' of 'adaptive' argument must be of length 1, ",
           nfactors, sprintf(" or of dimension c(%i, %i)", nfactors, nfactors))
    }
  }

  # samples <- spmirt_cpp(
  #   response = response,
  #   nobs = nobs, nitems = nitems, nfactors = nfactors, niter = niter,
  #   constrain_L = constrain_L,
  #   c_initial = c_initial, c_prior_mean = c_prior_mean, c_prior_sd = c_prior_sd,
  #   A_initial = A_initial, A_prior_mean = A_prior_mean, A_prior_sd = A_prior_sd,
  #   theta_init = theta_init
  #   )
  # return(samples)
  # return(list(constrain_L = constrain_L, A_initial = A_initial,
  #             A_prior_mean = A_prior_mean, A_prior_sd =  A_prior_sd))
  return(list(adap_Sigma_R = adap_Sigma_R, adap_Sigma_gp_sd =  adap_Sigma_gp_sd,
              adap_Sigma_gp_phi = adap_Sigma_gp_phi))
}
