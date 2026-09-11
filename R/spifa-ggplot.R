#' @title Gather Parameters into a Long Format Tibble
#'
#' @description
#' Reshapes a wide samples tibble (one column per parameter, as produced by
#' \code{\link{as_tibble.spifa}}) into long format (one row per
#' iteration/parameter pair), which is the shape expected by the
#' \code{\link{gg_trace}}/\code{\link{gg_density}} family of plotting
#' helpers.
#'
#' @param samples_wide A wide samples tibble, e.g. from
#' \code{\link{as_tibble.spifa}}.
#' @param each If not \code{NULL}, the number of columns that make up each
#' group of parameters (e.g. the number of items), used to additionally
#' split the gathered \code{Parameters} column into \code{group}/
#' \code{Parameter} columns.
#' @param keys Names to use for the group/parameter columns when \code{each}
#' is supplied.
#'
#' @return A long-format \code{\link[tibble]{tibble}} with columns
#' \code{iteration}, \code{Parameters}, and \code{Value} (plus \code{group}/
#' the second \code{keys} element when \code{each} is supplied).
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' nitems <- ncol(ipixuna$items)
#' nfactors <- 3
#'
#' # discrimination constraint: start with every item free to load on every
#' # factor, then restrict a few items per factor based on what each item is
#' # meant to measure (0 = no relationship, 1 = free parameter to estimate)
#' A <- matrix(1, nitems, nfactors)
#' A[c(4, 8), 1] <- 0
#' A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
#' A[c(5, 6), 3] <- 0
#' samples <- spifa(
#'   items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
#'   niter = 20, standardize = FALSE,
#'   constraints = list(discrimination = A))
#' wide <- as_tibble(samples, select = "c")
#' long <- gather.spifa(wide)
#' long
#' }
#'
#' @export
gather.spifa <- function (samples_wide, each = NULL,
                           keys = c("group", "Parameter")) {

  # Convert to long format
  samples_long <- samples_wide |>
    tibble::as_tibble() |>
    dplyr::mutate(iteration = 1:dplyr::n()) |>
    tidyr::gather(Parameters, Value, -iteration, factor_key = TRUE)

  if (!is.null(each)) {

    # Auxiliary variables to group
    groups <- paste0(keys[1], rep(1:each, ncol(samples_wide)/each))
    groups <- factor(groups, unique(groups))
    var <- paste0(keys[2], rep(1:(ncol(samples_wide)/each), each = each))
    names(groups) <- levels(samples_long$Parameters)
    names(var) <- levels(samples_long$Parameters)

    # Group parameters
    samples_long <- samples_long |>
      dplyr::mutate(groups = groups[Parameters], var = var[Parameters]) |>
    dplyr::select(-Parameters) |>
    tidyr::spread(var, Value)

  }

  return(samples_long)
}

#' @title Traceplot of Samples
#'
#' @description
#' Draws MCMC traceplots (iteration vs. value, one line per parameter) from
#' a long-format samples tibble, as produced by \code{\link{gather.spifa}}.
#'
#' @param df A wide \code{spifa} samples tibble (e.g. from
#' \code{\link{as_tibble.spifa}}), or a subset of its columns selected
#' via \code{select} in \code{\link{as_tibble.spifa}}.
#' @param wrap Logical; if \code{TRUE}, draw one facet per parameter instead
#' of overlaying them on a single panel.
#' @param legend Legend position passed to
#' \code{\link[ggplot2]{theme}(legend.position = ...)}.
#' @param ... Further arguments passed to \code{\link[ggplot2]{geom_path}}.
#'
#' @return A \code{ggplot} object.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' parameters <- attr(ipixuna, "parameters")
#' L_a <- (parameters$discrimination != 0) * 1
#' nfactors <- ncol(parameters$discrimination)
#' samples <- spifa(
#'   items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
#'   niter = 20, thin = 1, standardize = FALSE,
#'   constraints = list(discrimination = L_a, sd = rep(0.5, nfactors)))
#' as_tibble(samples, select = "c") %>% gg_trace(wrap = TRUE, alpha = 0.6)
#' }
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export
gg_trace <- function (df, wrap = FALSE, legend = "bottom", ...) {
  df <- gather.spifa(df)
  gg <- df %>%
    ggplot(aes(iteration, Value, group = Parameters, col = Parameters)) +
      geom_path(...)
  # in splits if required
  if (wrap) {
    gg <- gg + facet_wrap(~ Parameters, ncol = 1, scales = "free")
  }
  # theme
  gg <- gg + theme(legend.position = legend)
  return(gg)
}

#' @title Densities of Samples
#'
#' @description
#' Draws posterior density plots, one per parameter, from a samples tibble.
#' Densities can be overlaid as facets (default) or stacked as ridge plots
#' (\code{ridges = TRUE}, requires the \pkg{ggridges} package).
#'
#' @param df A wide \code{spifa} samples tibble (e.g. from
#' \code{\link{as_tibble.spifa}}).
#' @param ... Further arguments passed to \code{\link[ggplot2]{geom_density}}
#' (or to \code{ggridges::geom_density_ridges} when \code{ridges = TRUE}).
#' @param ridges Logical; if \code{TRUE}, draw ridge (stacked) densities
#' instead of faceted densities.
#'
#' @return A \code{ggplot} object.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' parameters <- attr(ipixuna, "parameters")
#' L_a <- (parameters$discrimination != 0) * 1
#' nfactors <- ncol(parameters$discrimination)
#' samples <- spifa(
#'   items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
#'   niter = 20, thin = 1, standardize = FALSE,
#'   constraints = list(discrimination = L_a, sd = rep(0.5, nfactors)))
#' as_tibble(samples, select = "c") %>% gg_density()
#' }
#'
#' @export
gg_density <- function (df, ..., ridges = FALSE) {
  df <- gather.spifa(df)
  df <- df %>%
    group_by(Parameters) %>%
    mutate(median = quantile(Value, 0.5))
  if (ridges) {
    if (!requireNamespace("ggridges", quietly = TRUE)) {
      stop("Package 'ggridges' is required for ridges = TRUE. ",
           "Install it with install.packages('ggridges').")
    }
    gg <- df %>%
      ggplot(aes(Value, Parameters, group = Parameters)) +
        ggridges::geom_density_ridges(...)
  } else  {
    gg <- df %>%
      ggplot(aes(Value, fill = Parameters)) +
      geom_density(...) +
      facet_wrap(~ Parameters, scales = "free")
  }
  # theme
  gg <- gg + theme(legend.position = "none")
  return(gg)
}

#' @title 2D Densities of Samples
#'
#' @description
#' Draws a 2D contour density plot of two parameters (e.g. two components
#' of a bivariate latent factor) from a samples tibble, optionally faceted
#' by group and/or highlighting a reference point (e.g. the true simulated
#' value).
#'
#' @param samples A samples tibble containing the columns \code{var1} and
#' \code{var2} (a wide samples tibble, or the output of
#' \code{\link{gather.spifa}} when \code{each} is used).
#' @param var1,var2 Bare (unquoted) names of the two columns in
#' \code{samples} to plot on the x and y axes.
#' @param each If not \code{NULL}, the number of columns making up each
#' group of parameters, passed to \code{\link{gather.spifa}} to facet the
#' plot by group.
#' @param keys Column names used for the group/parameter split, passed to
#' \code{\link{gather.spifa}} when \code{each} is supplied.
#' @param highlight An optional reference row (e.g. true parameter values)
#' to overlay as a point.
#' @param ncol Number of facet columns to use when \code{each} is supplied.
#' @param ... Further arguments passed to
#' \code{\link[ggplot2]{stat_density_2d}}.
#'
#' @return A \code{ggplot} object.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' parameters <- attr(ipixuna, "parameters")
#' L_a <- (parameters$discrimination != 0) * 1
#' nfactors <- ncol(parameters$discrimination)
#' samples <- spifa(
#'   items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
#'   niter = 20, thin = 1, standardize = FALSE,
#'   constraints = list(discrimination = L_a, sd = rep(0.5, nfactors)))
#' samples_tib <- as_tibble(samples)
#' gg_density2d(samples_tib, `c[1]`, `c[2]`)
#' }
#'
#' @export
gg_density2d <- function (samples, var1, var2, each = NULL,
                          keys = c("group", "Parameter"), highlight = NULL,
                          ncol = NULL, ...) {

  if (!is.null(highlight)) {
    aux_samples <- samples[1,]
    aux_samples[1, ] <- highlight
  }

  if (!is.null(each)) {
    samples <- gather.spifa(samples, each, keys)
    aux_samples <- gather.spifa(aux_samples, each, keys)
  }

  gg <- ggplot(samples, aes(!!substitute(var1), !!substitute(var2))) +
    stat_density2d(aes(fill = log(after_stat(level))),
                   geom = 'polygon', col = "black", ...) +
    # scale_fill_continuous(low="green",high="red") +
    guides(alpha="none")
  # +
  #   geom_point(...)

  if (!is.null(highlight)) {
    gg <- gg + geom_point(data = aux_samples, col = 2)
  }

  if (!is.null(each)) {
    if (!is.null(ncol)) {
      gg <- gg + facet_wrap(~ groups, scales = "free", ncol = ncol)
    } else {
      gg <- gg + facet_wrap(~ groups, scales = "free")
    }
  }


  return(gg)
}



#' @title 2D Scatterplot of Samples
#'
#' @description
#' Draws a 2D scatter/path plot of two parameters (e.g. two components of a
#' bivariate latent factor) from a samples tibble, optionally faceted by
#' group and/or highlighting a reference point (e.g. the true simulated
#' value).
#'
#' @param samples A samples tibble containing the columns \code{var1} and
#' \code{var2} (a wide samples tibble, or the output of
#' \code{\link{gather.spifa}} when \code{each} is used).
#' @param var1,var2 Bare (unquoted) names of the two columns in
#' \code{samples} to plot on the x and y axes.
#' @param each If not \code{NULL}, the number of columns making up each
#' group of parameters, passed to \code{\link{gather.spifa}} to facet the
#' plot by group.
#' @param keys Column names used for the group/parameter split, passed to
#' \code{\link{gather.spifa}} when \code{each} is supplied.
#' @param highlight An optional reference row (e.g. true parameter values)
#' to overlay as a point.
#' @param ncol Number of facet columns to use when \code{each} is supplied.
#' @param points_alpha Alpha (transparency) used for the scatter points.
#'
#' @return A \code{ggplot} object.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' parameters <- attr(ipixuna, "parameters")
#' L_a <- (parameters$discrimination != 0) * 1
#' nfactors <- ncol(parameters$discrimination)
#' samples <- spifa(
#'   items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
#'   niter = 20, thin = 1, standardize = FALSE,
#'   constraints = list(discrimination = L_a, sd = rep(0.5, nfactors)))
#' samples_tib <- as_tibble(samples)
#' gg_scatter(samples_tib, `c[1]`, `c[2]`)
#' }
#'
#' @export
gg_scatter <- function (samples, var1, var2, each = NULL,
                        keys = c("group", "Parameter"), highlight = NULL,
                        ncol = NULL, points_alpha = 0.5) {

  if (!is.null(highlight)) {
    aux_samples <- samples[1,]
    aux_samples[1, ] <- highlight
  }

  if (!is.null(each)) {
    samples <- gather.spifa(samples, each, keys)
    aux_samples <- gather.spifa(aux_samples, each, keys)
  }

  gg <- ggplot(samples, aes(!!substitute(var1), !!substitute(var2))) +
    geom_point(alpha = points_alpha) +
    geom_path(alpha = 0.4, linetype = 2)

  if (!is.null(highlight)) {
    gg <- gg + geom_point(data = aux_samples, col = 2, size = 2)
  }

  if (!is.null(each)) {
    if (!is.null(ncol)) {
      gg <- gg + facet_wrap(~ groups, scales = "free", ncol = ncol)
    } else {
      gg <- gg + facet_wrap(~ groups, scales = "free")
    }
  }

  return(gg)
}






#' @title Horizontal Errorbar Plot of Samples
#'
#' @description
#' Plots posterior medians with 80\%/95\% credible interval errorbars (one
#' row per parameter), using the output of \code{\link{summary.spifa}}.
#'
#' @param df_summary A summary tibble from \code{\link{summary.spifa}}, with
#' columns \code{Parameters}, \code{2.5\%}, \code{10\%}, \code{50\%},
#' \code{90\%}, \code{97.5\%}.
#' @param sorted Logical; if \code{TRUE}, plot against the posterior median
#' on both axes (for use with faceting/sorting upstream) instead of against
#' \code{Parameters}.
#' @param colors Colors used for the 95\% and 80\% credible interval bars.
#' @param ... Further arguments passed to
#' \code{\link[ggplot2]{geom_errorbarh}}.
#'
#' @return A \code{ggplot} object.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' parameters <- attr(ipixuna, "parameters")
#' L_a <- (parameters$discrimination != 0) * 1
#' nfactors <- ncol(parameters$discrimination)
#' samples <- spifa(
#'   items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
#'   niter = 20, thin = 1, standardize = FALSE,
#'   constraints = list(discrimination = L_a, sd = rep(0.5, nfactors)))
#' gg_errorbarh(summary(samples, select = "c"))
#' }
#'
#' @export
gg_errorbarh <- function (df_summary, sorted = FALSE,
                          colors = c(rgb(1,0.5,0.1), "black"), ...) {

  if (sorted) {
    gg <- df_summary %>% ggplot(., aes(median, median))
  } else {
    gg <- df_summary %>% ggplot(., aes(median, variable))
  }

  gg <- gg +
    geom_errorbarh(aes(xmin = q2.5, xmax = q97.5, col = "95%"),
                   height = 0, ...) +
    geom_errorbarh(aes(xmin = q10, xmax = q90, col = "80%"), linewidth = 2,
                   height = 0, ...) +
    geom_point(size = 2)
  # colors
  gg <- gg + scale_colour_manual(values = colors)
  # labels
  gg <- gg + labs(colour = "Credible Intervals:", x = "Value")
  # theme
  gg <- gg + theme(legend.position = "bottom")
  return(gg)
}

#' @title Errorbar Plot of Samples
#'
#' @description
#' Plots posterior medians with 80\%/95\% credible interval errorbars (one
#' column per parameter), using the output of \code{\link{summary.spifa}}.
#'
#' @param df_summary A summary tibble from \code{\link{summary.spifa}}, with
#' columns \code{variable}, \code{q2.5}, \code{q10}, \code{median},
#' \code{q90}, \code{q97.5}.
#' @param sorted Logical; if \code{TRUE} (default), plot against the
#' posterior median on both axes (for use with faceting/sorting upstream)
#' instead of against \code{variable}.
#' @param colors Colors used for the 95\% and 80\% credible interval bars.
#' @param ... Further arguments passed to \code{\link[ggplot2]{geom_errorbar}}.
#'
#' @return A \code{ggplot} object.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' parameters <- attr(ipixuna, "parameters")
#' L_a <- (parameters$discrimination != 0) * 1
#' nfactors <- ncol(parameters$discrimination)
#' samples <- spifa(
#'   items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
#'   niter = 20, thin = 1, standardize = FALSE,
#'   constraints = list(discrimination = L_a, sd = rep(0.5, nfactors)))
#' gg_errorbar(summary(samples, select = "c"), sorted = FALSE)
#' }
#'
#' @export
gg_errorbar <- function (df_summary, sorted = TRUE,
                         colors = c(rgb(1,0.5,0.1), "black"), ...) {

  if (sorted) {
    gg <- df_summary %>% ggplot(., aes(median, median))
  } else {
    gg <- df_summary %>% ggplot(., aes(variable, median))
  }

  gg <- gg +
    geom_errorbar(aes(ymin = q2.5, ymax = q97.5, col = "95%"),
                   width = 0, ...) +
    geom_errorbar(aes(ymin = q10, ymax = q90, col = "80%"), linewidth = 2,
                   width = 0, ...) +
    geom_point(size = 2)
  # colors
  gg <- gg + scale_colour_manual(values = colors)
  # labels
  gg <- gg + labs(colour = "Credible Intervals:", x = "Value")
  # theme
  gg <- gg + theme(legend.position = "bottom")
  return(gg)
}
