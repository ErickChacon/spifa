# remove restricted discrimination parameters
drop_restricted <- function (x) {

  fit_args <- attr(x, "fit_args")
  predict_setup <- attr(x, "predict_setup")

  index <- which(fit_args$constrain_L == 0, arr.ind = TRUE)
  a_restricted <- paste0("A[", index[, 1], ",", index[, 2], "]")
  a_unrestricted <- setdiff(posterior::variables(x), a_restricted)
  x <- posterior::subset_draws(x, variable = a_unrestricted)

  attr(x, "fit_args") <- fit_args
  attr(x, "predict_setup") <- predict_setup
  return(x)
}

# plotmath parse of parameters
parse_parameter <- function (parameter) {
  to_math <- function (p) gsub("\\[(.+),(.+)\\]", "[list(\\1,\\2)]", p)
  factor(to_math(as.character(parameter)), levels = to_math(levels(parameter)))
}

# converts samples to long format
draws_long <- function (x) {
  wide <- x |> posterior::as_draws_df() |> tibble::as_tibble()
  pars <- setdiff(names(wide), c(".chain", ".iteration", ".draw"))
  wide |>
    dplyr::mutate(iteration = .iteration) |>
    tidyr::pivot_longer(dplyr::all_of(pars), names_to = "parameter", values_to = "value") |>
    dplyr::mutate(parameter = factor(parameter, levels = pars)) |>
    dplyr::select(iteration, parameter, value)
}

theme_trace <- function (legend = "none") {
  theme_minimal() +
    theme(panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linetype = "dashed"),
          axis.line = element_line(),
          strip.background = element_blank(),
          strip.text.y = element_text(angle = 0),
          legend.position = legend)
}

#' @title Traceplot of Samples
#'
#' @description
#' Draws MCMC traceplots (iteration vs. value) directly from a fitted
#' \code{spifa} model, in one of two formats: faceted (\code{facet = TRUE},
#' the default; one panel per parameter, with its own free y-scale, so a
#' slow-mixing or small-variance parameter isn't visually flattened by
#' others sharing the same axis) or overlaid (\code{facet = FALSE}; all
#' series on a single panel, for a quick glance at overall convergence).
#' Replaces the older \code{\link{gg_trace}}/\code{\link{gather.spifa}}
#' pair: it reshapes draws directly, so it works on \code{x} with no
#' intermediate conversion.
#'
#' @param x A fitted \code{spifa} model.
#' @param select Parameters to plot, passed to the \code{variable} argument
#' of \code{\link[posterior]{subset_draws}}: either a block name (e.g.
#' \code{"A"}, matching every parameter in that block) or one or more full
#' indexed names (e.g. \code{"c[1]"}, \code{paste0("A[", 1:10, ",1]")}).
#' @param facet Logical; if \code{TRUE} (default), draw one panel per
#' parameter; if \code{FALSE}, overlay every series on a single panel. See
#' Description.
#' @param burnin Number of initial iterations to discard.
#' @param thin Thinning interval applied after \code{burnin}.
#' @param nshow If \code{select} matches more than \code{nshow} parameters,
#' a random (sorted) subsample of \code{nshow} of them is shown instead of
#' all of them. Set to \code{NULL} to always show every matched parameter.
#' @param ncol Number of columns in the facet grid (\code{facet = TRUE}
#' only); defaults to a single column.
#' @param legend Legend position (\code{facet = FALSE} only): one of
#' \code{"bottom"}, \code{"right"}, \code{"none"}, or a logical. Defaults to
#' \code{NULL}, which shows the legend unless there are more than 10
#' series, since a legend that large stops being readable and can dwarf the
#' plot itself.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 0, niter = 1000)
#'
#' plot_trace(samples, select = "c", facet = TRUE)
#' plot_trace(samples, select = "c", facet = FALSE)
#'
#' # more than nshow (10) parameters: a random subsample is shown
#' plot_trace(samples, select = "A", facet = TRUE, nshow = 6)
#'
#' # explicit selection instead of a random subsample
#' plot_trace(samples, select = paste0("A[", 1:10, ",1]"), facet = TRUE)
#' }
#'
#' @export
plot_trace <- function (x, select, facet = TRUE,
                        burnin = 0, thin = 1, nshow = 10, ncol = 1,
                        legend = NULL, ...) {
  x <- drop_restricted(x)

  # create data
  niter <- posterior::niterations(x)
  df <- x |>
    posterior::subset_draws(variable = select, iteration = (burnin + 1):niter) |>
    posterior::thin_draws(thin) |>
    draws_long() |>
    dplyr::mutate(parameter = parse_parameter(parameter))

  # random sorted subsample when there are more parameters than nshow
  pars <- unique(df$parameter)
  if (!is.null(nshow) && length(pars) > nshow) {
    pars_keep <- sort(sample(pars, nshow))
    df <- dplyr::filter(df, parameter %in% pars_keep)
  }

  # two figure types
  if (!facet) {
    # remove legend if more than 10 parameters
    if (is.null(legend)) {
      legend <- if (length(unique(df$parameter)) > 10) "none" else "bottom"
    } else if (is.logical(legend)) {
      legend <- if (legend) "bottom" else "none"
    }
    gg <- ggplot(df, aes(iteration, value, group = parameter, col = parameter)) +
      geom_path(alpha = 0.6) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_colour_discrete(labels = function(x) parse(text = x)) +
      labs(x = "Iteration", y = "Value", col = "Parameter") +
      theme_trace(legend = legend)
  } else {
    gg <- ggplot(df, aes(iteration, value, col = parameter)) +
      geom_path(linewidth = 0.2) +
      facet_wrap(~ parameter, ncol = ncol, scales = "free_y", strip.position = "right",
                 labeller = label_parsed) +
      scale_x_continuous(expand = c(0, 0)) +
      labs(x = "Iteration", y = "Value") +
      theme_trace(legend = "none")
  }

  return(gg)
}

#' @title Density Plot of Samples
#'
#' @description
#' Draws posterior density curves directly from a fitted \code{spifa}
#' model, in one of two formats: overlaid (\code{facet = FALSE}, the
#' default; stacked, slightly-overlapping ridgeline plot via
#' \code{\link[ggridges]{geom_density_ridges}} -- the more insightful
#' default for comparing many parameters' shapes and locations at a glance) or
#' faceted (\code{facet = TRUE}; one panel per parameter). Shows the
#' density shape only -- for credible intervals and point estimates, see
#' \code{\link{plot_interval}}. The density itself is computed directly by
#' \code{\link[ggplot2]{geom_density}}/\code{\link[ggridges]{geom_density_ridges}}
#' from the raw draws, reshaped via \code{\link[posterior]{as_draws_df}}.
#'
#' @param x A fitted \code{spifa} model.
#' @param select Parameters to plot, as in \code{\link{plot_trace}}.
#' @param facet Logical; if \code{TRUE}, draw one panel per parameter; if
#' \code{FALSE} (default), overlay every parameter as a ridgeline plot. See
#' Description.
#' @param burnin Number of initial iterations to discard.
#' @param thin Thinning interval applied after \code{burnin}.
#' @param nshow As in \code{\link{plot_trace}}.
#' @param ncol Number of columns in the facet grid (\code{facet = TRUE}
#' only); defaults to a single column.
#' @param facet_scales The \code{scales} argument of
#' \code{\link[ggplot2]{facet_wrap}} (\code{facet = TRUE} only): one of
#' \code{"free"} (default; each panel gets its own x/y-axis), \code{"free_y"}
#' (one shared x-axis, line shown only on the bottom panel of each column,
#' matching \code{\link{plot_trace}}'s facet -- useful for comparing
#' parameters on a similar scale), \code{"free_x"}, or \code{"fixed"}.
#' Parameters with very different value ranges can render compressed/hard
#' to read with a shared x-axis.
#' @param scale Amount of vertical overlap between ridges (\code{facet =
#' FALSE} only), passed to \code{\link[ggridges]{geom_density_ridges}}.
#' Defaults to \code{1.2}.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 0, niter = 1000)
#'
#' plot_density(samples, select = "c", facet = FALSE)
#' plot_density(samples, select = "c", facet = TRUE)
#'
#' # more than nshow (10) parameters: a random subsample is shown
#' plot_density(samples, select = "A", facet = FALSE, nshow = 6)
#' }
#'
#' @export
plot_density <- function (x, select, facet = FALSE,
                           burnin = 0, thin = 1, nshow = 10, ncol = 1,
                           facet_scales = "free", scale = 1.2, ...) {
  x <- drop_restricted(x)

  # create data
  niter <- posterior::niterations(x)
  df <- x |>
    posterior::subset_draws(variable = select, iteration = (burnin + 1):niter) |>
    posterior::thin_draws(thin) |>
    draws_long() |>
    dplyr::mutate(parameter = parse_parameter(parameter))

  # random sorted subsample when there are more parameters than nshow
  pars <- unique(df$parameter)
  if (!is.null(nshow) && length(pars) > nshow) {
    pars_keep <- sort(sample(pars, nshow))
    df <- dplyr::filter(df, parameter %in% pars_keep)
  }

  # figure types
  if (!facet) {
    gg <- ggplot(df, aes(value, parameter, fill = parameter)) +
      ggridges::geom_density_ridges(alpha = 0.5, scale = scale, linewidth = 0.4,
                                     rel_min_height = 0.01) +
      scale_y_discrete(labels = function(x) parse(text = x)) +
      labs(x = "Value", y = NULL) +
      theme_trace(legend = "none")
  } else {
    gg <- ggplot(df, aes(value, fill = parameter)) +
      geom_density(alpha = 0.5, linewidth = 0.4, trim = TRUE) +
      facet_wrap(~ parameter, ncol = ncol, scales = facet_scales,
                 strip.position = "right", labeller = label_parsed) +
      labs(x = "Value", y = "Density") +
      theme_trace(legend = "none")
  }

  return(gg)
}

#' @title Interval Plot of Samples
#'
#' @description
#' Draws a caterpillar/forest plot of posterior credible intervals directly
#' from a fitted \code{spifa} model: one row per parameter, with a thin
#' line for the \code{prob_outer} interval, a thick line for the
#' \code{prob} interval, and a point at the \code{point_est}, computed
#' directly from the raw draws (quantiles/median/mean). Unlike
#' \code{\link{plot_trace}}/\code{\link{plot_density}}, this is a single
#' combined view by design -- comparing intervals side by side is the whole
#' point, so there is no faceted alternative -- but \code{sort} can reorder
#' parameters by their point estimate, which a facet can't do usefully
#' across independent panels.
#'
#' @param x A fitted \code{spifa} model.
#' @param select Parameters to plot, as in \code{\link{plot_trace}}.
#' @param horizontal Logical; if \code{FALSE} (default), parameters run
#' along the x-axis and values run along the y-axis, matching the
#' \code{ci_intervals()} convention used in the SPIFA paper's own figures;
#' if \code{TRUE}, the axes are swapped (a forest-plot layout).
#' @param burnin Number of initial iterations to discard.
#' @param thin Thinning interval applied after \code{burnin}.
#' @param nshow As in \code{\link{plot_trace}} (a random subsample when
#' \code{select} matches more than \code{nshow} parameters), but
#' \code{NULL} (show every matched parameter) by default: unlike a faceted
#' plot, a single interval plot stays readable with many more than 10 rows.
#' @param prob Width of the thick (inner) credible interval (a central
#' quantile interval). Defaults to \code{0.5}.
#' @param prob_outer Width of the thin (outer) credible interval. Defaults
#' to \code{0.9}.
#' @param point_est Either \code{"median"} (default) or \code{"mean"}.
#' @param sort Logical; if \code{TRUE}, reorder parameters by their point
#' estimate instead of their natural order. Defaults to \code{FALSE}.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 0, niter = 1000)
#'
#' plot_interval(samples, select = "c")
#' plot_interval(samples, select = "c", sort = TRUE)
#' plot_interval(samples, select = "A", horizontal = TRUE)
#' }
#'
#' @export
plot_interval <- function (x, select, horizontal = FALSE,
                            burnin = 0, thin = 1, nshow = NULL,
                            prob = 0.5, prob_outer = 0.9,
                            point_est = c("median", "mean"), sort = FALSE, ...) {
  x <- drop_restricted(x)
  point_est <- match.arg(point_est)

  # create data: one row per parameter
  niter <- posterior::niterations(x)
  df <- x |>
    posterior::subset_draws(variable = select, iteration = (burnin + 1):niter) |>
    posterior::thin_draws(thin) |>
    draws_long() |>
    dplyr::mutate(parameter = parse_parameter(parameter)) |>
    dplyr::group_by(parameter) |>
    dplyr::summarise(
      ll = stats::quantile(value, (1 - prob_outer) / 2, names = FALSE),
      l = stats::quantile(value, (1 - prob) / 2, names = FALSE),
      m = if (point_est == "median") stats::median(value) else mean(value),
      h = stats::quantile(value, 1 - (1 - prob) / 2, names = FALSE),
      hh = stats::quantile(value, 1 - (1 - prob_outer) / 2, names = FALSE),
      .groups = "drop")

  # random sorted subsample when there are more parameters than nshow
  pars <- unique(df$parameter)
  if (!is.null(nshow) && length(pars) > nshow) {
    pars_keep <- base::sort(sample(pars, nshow))
    df <- dplyr::filter(df, parameter %in% pars_keep)
  }

  # sort by point estimate instead of the natural parameter order
  if (sort) df <- dplyr::mutate(df, parameter = stats::reorder(parameter, m))

  outer_colour <- "black"
  inner_colour <- grDevices::rgb(1, 0.5, 0.1)
  # figure types
  if (horizontal) {
    gg <- ggplot(df, aes(y = parameter)) +
      geom_segment(aes(x = ll, xend = hh, yend = parameter), linewidth = 0.4,
                   colour = outer_colour) +
      geom_segment(aes(x = l, xend = h, yend = parameter), linewidth = 1.5,
                   colour = inner_colour) +
      geom_point(aes(x = m), size = 2, colour = outer_colour) +
      scale_y_discrete(labels = function(x) parse(text = x)) +
      labs(x = "Value", y = NULL) +
      theme_trace(legend = "none")
  } else {
    gg <- ggplot(df, aes(x = parameter)) +
      geom_segment(aes(y = ll, yend = hh, xend = parameter), linewidth = 0.4,
                   colour = outer_colour) +
      geom_segment(aes(y = l, yend = h, xend = parameter), linewidth = 1.5,
                   colour = inner_colour) +
      geom_point(aes(y = m), size = 2, colour = outer_colour) +
      scale_x_discrete(labels = function(x) parse(text = x)) +
      labs(x = NULL, y = "Value") +
      theme_trace(legend = "none")
  }

  return(gg)
}

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
#' as_tibble(samples, select = "c") |> gg_density()
#' }
#'
#' @export
gg_density <- function (df, ..., ridges = FALSE) {
  df <- gather.spifa(df)
  df <- df |>
    group_by(Parameters) |>
    mutate(median = quantile(Value, 0.5))
  if (ridges) {
    if (!requireNamespace("ggridges", quietly = TRUE)) {
      stop("Package 'ggridges' is required for ridges = TRUE. ",
           "Install it with install.packages('ggridges').")
    }
    gg <- df |>
      ggplot(aes(Value, Parameters, group = Parameters)) +
        ggridges::geom_density_ridges(...)
  } else  {
    gg <- df |>
      ggplot(aes(Value, fill = Parameters)) +
      geom_density(...) +
      facet_wrap(~ Parameters, scales = "free")
  }
  # theme
  gg <- gg + theme(legend.position = "none")
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
    gg <- df_summary |> ggplot(aes(median, median))
  } else {
    gg <- df_summary |> ggplot(aes(median, variable))
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
    gg <- df_summary |> ggplot(aes(median, median))
  } else {
    gg <- df_summary |> ggplot(aes(variable, median))
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
