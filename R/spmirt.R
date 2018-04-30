
spmirt <- function (response, nobs, nitems, nfactors,
                    L_rest,
                    niter = 100,
                    theta_init,
                    c_opt = list(initial = 1, prior_mean = 1, prior_sd = 1),
                    A_opt = list(initial = 1, prior_mean = 1, prior_sd = 1)) {

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
    if (is.null(A_opt$prior_mean)) {
      A_initial <- matrix(1, nitems, nfactors)
    } else {
      A_initial <- A_prior_mean
    }
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

  samples <- spmirt_cpp(
    response = response,
    nobs = nobs, nitems = nitems, nfactors = nfactors,
    L_rest = L_a,
    niter = niter,
    theta_init = theta_init,
    c_initial = c_initial, c_prior_mean = c_prior_mean, c_prior_sd = c_prior_sd,
    A_initial = A_initial, A_prior_mean = A_prior_mean, A_prior_sd = A_prior_sd)

  return(samples)
  # return(list(c_initial, c_prior_mean, c_prior_sd, A_initial, A_prior_mean, A_prior_sd))
}
