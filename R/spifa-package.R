

#' @title spifa: Spatial (Geostatistical) Item Factor Analysis
#'
#' @docType package
#' @name spifa-package
#'
#' @description
#' Fits item factor analysis models for binary responses using full Bayesian
#' inference via Gibbs sampling and auxiliary variables with a probit link
#' function, optionally extending the latent factors with a multivariate
#' Gaussian process to capture spatial dependence. See \code{\link{spifa}}
#' for the main model-fitting function.
#'
#' @author Erick A. Chacon-Montalvan \email{e.chaconmontalvan@lancaster.ac.uk}
#' @author Emanuele Giorgi, \email{e.giorgi@lancaster.ac.uk}
#' @author Benjamin M. Taylor, \email{b.taylor1@lancaster.ac.uk}
#'
#' @importFrom magrittr %>%
#' @importFrom stats dist model.matrix quantile rnorm setNames update
#' @importFrom grDevices rgb
"_PACKAGE"

utils::globalVariables(c(
  ".", "Parameters", "Value", "iteration", "level",
  "2.5%", "10%", "50%", "90%", "97.5%"
  ))

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @importFrom tibble as_tibble
#' @export
tibble::as_tibble
