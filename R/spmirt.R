
spmirt <- function (response, predictors = NULL, coordinates = NULL,
    nobs, nitems, nfactors, niter = 100,
    constrains = list(A = NULL, W = NULL),
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
      stop("element 'prior_mean' of 'c_opt' must be of length 1 or 'nitems'")
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
      stop("element 'prior_sd' of 'c_opt' must be of length 1 or 'nitems'")
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
      stop("element 'initial' of 'c_opt' must be of length 1 or 'nitems'")
    }
  }

  # Optional arguments for parameter A, checking adequate sizes

  if (is.null(A_opt$prior_mean)) {
    A_prior_mean <- matrix(0, nitems, nfactors)
    diag(A_prior_mean) <- 1
  } else {
    if (length(A_opt$prior_mean) == 1) {
      A_prior_mean <- matrix(A_opt$prior_mean, nitems, nfactors)
    } else if (all(dim(A_opt$prior_mean) == c(nitems, nfactors))) {
      A_prior_mean <- A_opt$prior_mean
    } else {
      stop("element 'prior_mean' of 'A_opt' must be of length 1 or
           dimension c(nitems, nfactors)")
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
    } else if (all(dim(A_opt$prior_mean) == c(nitems, nfactors))) {
      A_prior_sd <- A_opt$prior_sd
    } else {
      stop("element 'prior_sd' of 'A_opt' must be of length 1 or
           dimension c(nitems, nfactors)")
    }
  }

  if (is.null(A_opt$initial)) {
      A_initial <- A_prior_mean
  } else {
    if (length(A_opt$initial) == 1) {
      A_initial <- matrix(A_opt$initial, nitems, nfactors)
    } else if (all(dim(A_opt$initial) == c(nitems, nfactors))) {
      A_initial <- A_opt$initial
    } else {
      stop("element 'initial' of 'A_opt' must be of length 1 or
           dimension c(nitems, nfactors)")
    }
  }

  # Restrictions arguments for parameter A

  if (is.null(constrains$A)) {
    A_aux <- matrix(NA, nitems, nfactors)
    constrain_L <- lower.tri(A_aux, diag = TRUE) * 1
  } else {
    if (all(dim(constrains$A) == c(nitems, nfactors))) {
      constrain_L <- constrains$A
    } else {
      stop("element 'A' of 'constrains' must be of dimension c(nitems, nfactors)")
    }
  }

  # Restrictions arguments for Gaussian Processes W(s)

  # samples <- spmirt_cpp(
  #   response = response,
  #   nobs = nobs, nitems = nitems, nfactors = nfactors, niter = niter,
  #   constrain_L = constrain_L,
  #   c_initial = c_initial, c_prior_mean = c_prior_mean, c_prior_sd = c_prior_sd,
  #   A_initial = A_initial, A_prior_mean = A_prior_mean, A_prior_sd = A_prior_sd,
  #   theta_init = theta_init
  #   )
  # return(samples)
  return(list(constrain_L = constrain_L, A_initial = A_initial,
              A_prior_mean = A_prior_mean, A_prior_sd =  A_prior_sd))
}
