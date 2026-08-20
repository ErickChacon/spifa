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
#' @importFrom stats delete.response dist model.frame model.matrix
#' @importFrom stats model.response quantile rnorm setNames terms update
#' @importFrom grDevices rgb
#' @useDynLib spifa, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppTN
"_PACKAGE"

#' @title Simulated Food Insecurity Survey Data for Ipixuna
#'
#' @description
#' Simulated household survey data for the rural area of Ipixuna
#' (Para state, Brazil), used to illustrate spatial item factor analysis.
#' Ten binary items, a household-level `wealth` covariate, and spatial
#' coordinates are simulated for each household so that they reproduce the
#' spatial multidimensional item response structure described in the SPIFA
#' paper (Chacon-Montalvan, Giorgi and Taylor). The data is provided as an
#' `sf` object, already in the format [spifa()] requires.
#'
#' @format An `sf` tibble with one row per household (100 rows) and columns:
#' \describe{
#'   \item{id}{Household identifier.}
#'   \item{wealth}{A simulated household-level covariate.}
#'   \item{items}{A 100x10 binary (0/1) matrix of responses to the food
#'     insecurity items.}
#'   \item{geometry}{An `sfc_POINT` geometry column with the household's
#'     spatial coordinates (longitude/latitude), concentrated towards the
#'     town centre.}
#' }
#' It has an attribute `"parameters"` holding the true generating parameters
#' used in the simulation (`easiness`, `discrimination`, `abilities`
#' (the true latent factor scores, useful for checking recovery), `effect`,
#' `resid_params` (list with `sd` and `corr`, the residual covariance of the
#' latent factors), `mgp_params` (list with `sd` and `phi`, the spatial
#' Gaussian process hyperparameters)), which is convenient for specifying
#' sensible restrictions and checking recovery of the model in examples and
#' tests.
#'
#' @source Simulated; see `data-raw/ipixuna.R` for the generating code and
#'   `attr(ipixuna, "parameters")` for the generating values.
#'
#' @name ipixuna
#' @docType data
#'
#' @examples
#' data(ipixuna)
#' head(dplyr::select(ipixuna, -geometry, -items))
#' attr(ipixuna, "parameters")
NULL

#' @title Boundary of the Ipixuna Rural Area
#'
#' @description
#' Polygon boundary of the rural area of Ipixuna (Para state, Brazil), used
#' for mapping predictions from spatial item factor analysis models fitted
#' to the [ipixuna] data.
#'
#' @format An `sf`/`sfc` polygon object.
#'
#' @name ipixuna_boundary
#' @docType data
#'
#' @examples
#' data(ipixuna_boundary)
#' plot(sf::st_geometry(ipixuna_boundary))
NULL

## variables used in NSE contexts
utils::globalVariables(c(
  ".", "Parameters", "Value", "iteration", "level",
  "variable", "median", "q2.5", "q10", "q90", "q97.5"
  ))

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @importFrom tibble as_tibble
#' @export
tibble::as_tibble
