

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


