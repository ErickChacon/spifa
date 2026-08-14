#' @title Simulated Food Insecurity Survey Data for Ipixuna
#'
#' @description
#' Simulated household survey data for the rural area of Ipixuna
#' (Para state, Brazil), used to illustrate spatial item factor analysis.
#' Ten binary items (`Item 1`-`Item 10`), a covariate `x1`, and spatial
#' coordinates are simulated for each household so that they reproduce the
#' spatial multidimensional item response structure described in the SPIFA
#' paper (Chacon-Montalvan, Giorgi and Taylor). The data are provided in two
#' equivalent shapes: `ipixuna` (item responses in long format) and
#' `ipixuna_wide` (one row per household, one column per item).
#'
#' @format `ipixuna_wide` is a tibble with one row per household (200 rows)
#' and columns:
#' \describe{
#'   \item{id}{Household identifier.}
#'   \item{x1}{A simulated household-level covariate.}
#'   \item{coords}{An `sfc_POINT` geometry column with the household's
#'     spatial coordinates (longitude/latitude).}
#'   \item{Theta1, Theta2}{The true simulated values of the two latent
#'     factors, useful for checking recovery.}
#'   \item{`Item 1`, ..., `Item 10`}{Binary (0/1) responses to the food
#'     insecurity items.}
#' }
#' It has an attribute `"parameters"` holding the true generating parameters
#' used in the simulation (`n`, `q`, `m`, `g`, `easiness`, `discrimination`,
#' `beta`, `Cov`, `A`, `cor.params`), which is convenient for specifying
#' sensible restrictions and checking recovery of the model in examples and
#' tests.
#'
#' `ipixuna` is the long-format equivalent (2000 rows, one per
#' household/item pair), with columns `id`, `x1`, `coords`, `Theta1`,
#' `Theta2`, `prob` (true response probability), `response` (the binary
#' outcome), and `response_label` (the item identifier).
#'
#' Note: `coords` uses an older `sf` coordinate-reference-system
#' representation and will print a one-time "old-style crs object" message;
#' this is harmless but can be silenced with
#' `sf::st_set_crs(ipixuna_wide$coords, sf::st_crs(ipixuna_wide$coords))`.
#'
#' @source Simulated; see `attr(ipixuna_wide, "parameters")` for the
#'   generating values.
#'
#' @name ipixuna
#' @aliases ipixuna_wide
#' @docType data
#'
#' @examples
#' data(ipixuna)
#' head(dplyr::select(ipixuna, -coords))
#' head(dplyr::select(ipixuna_wide, -coords))
#' attr(ipixuna_wide, "parameters")
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
