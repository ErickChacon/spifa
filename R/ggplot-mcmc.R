
#' @title Convert Samples to Tibble
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' as_tibble(samples)
#' as_tibble(samples, "beta")
#'
#' @importFrom tibble as_tibble
#'
#' @export
as_tibble.spmirt.list <- function (samples, burnin = 0, thin = 1,
                              select = names(samples)) {
  samples <- samples[select]
  params <- names(samples)
  samples <- samples %>%
    purrr::map(tibble::as_tibble) %>%
    dplyr::bind_cols()
  niter <- nrow(samples)
  samples <- samples[seq(burnin + 1, niter, thin),]
  class(samples) <- c("spmirt", class(samples))
  return(samples)
}

#' @title Gather Parameters into a Long Format Tibble
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' wide <- as_tibble(samples, "beta")
#' (long <- gather(wide))
#' class(long)
#'
#' @export
gather.spmirt <- function (samples_wide, each = NULL,
                           keys = c("group", "Parameter")) {

  # Convert to long format
  samples_long <- samples_wide %>%
    tibble::as_tibble() %>%
    dplyr::mutate(iteration = 1:n()) %>%
    tidyr::gather(Parameters, Value, -iteration, factor_key = TRUE)

  if (!is.null(each)) {

    # Auxiliary variables to group
    groups <- paste0(keys[1], rep(1:each, ncol(samples_wide)/each))
    groups <- factor(groups, unique(groups))
    var <- paste0(keys[2], rep(1:(ncol(samples_wide)/each), each = each))
    names(groups) <- levels(samples_long$Parameters)
    names(var) <- levels(samples_long$Parameters)

    # Group parameters
    samples_long <- samples_long %>%
      dplyr::mutate(groups = groups[Parameters], var = var[Parameters]) %>%
    dplyr::select(-Parameters) %>%
    tidyr::spread(var, Value)

  }

  return(samples_long)
}

#' @title Summary of Samples
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' summary(samples)
#'
#' @export
summary.spmirt <- function (samples, select = names(samples)) {

  # df <- as_tibble.spmirt(samples, select = select)
  df <- samples
  df_names <- names(df)
  df <- df %>%
    map(function (x) quantile(x, c(0.025, 0.1, 0.5, 0.9, 0.975))) %>%
    reduce(rbind) %>%
    as_tibble() %>%
    mutate(Parameters = factor(df_names, df_names))
  df <- df[c(ncol(df), 1:(ncol(df)-1))]
  return(df)
}

#' @title Traceplot of Samples
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' as_tibble(samples, "beta") %>% gg_trace(wrap = TRUE, alpha = 0.6)
#' as_tibble(samples, "corr_chol") %>% gg_trace(alpha = 0.6)
#'
#'
#' @export
gg_trace <- function (df, wrap = FALSE, legend = "bottom", ...) {
  df <- gather.spmirt(df)
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
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' as_tibble(samples, "corr_chol") %>%
#'   gg_density_ridges(aes(fill = Parameters), scale = 2, alpha = 0.5)
#'
#' @export
gg_density <- function (df, ..., ridges = FALSE) {
  df <- gather.spmirt(df)
  if (ridges) {
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
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' as_tibble(samples, "corr_chol") %>%
#'   gg_density_ridges(aes(fill = Parameters), scale = 2, alpha = 0.5)
#'
#' @export
gg_density2d <- function (samples, var1, var2, each = NULL,
                           keys = c("group", "Parameter"), highlight = NULL) {

  if (!is.null(highlight)) {
    aux_samples <- samples[1,]
    aux_samples[1, ] <- highlight
  }

  if (!is.null(each)) {
    samples <- gather.spmirt(samples, each, keys)
    aux_samples <- gather.spmirt(aux_samples, each, keys)
  }

  gg <- ggplot(samples, aes_(substitute(var1), substitute(var2))) +
    stat_density2d(aes(fill = ..level.., alpha = ..level..),
                   geom = 'polygon', colour = 'black') +
            scale_fill_continuous(low="green",high="red") +
            guides(alpha="none") +
            geom_point(alpha = 0.5)

  if (!is.null(highlight)) {
    gg <- gg + geom_point(data = aux_samples, col = 2, size = 2)
  }

  if (!is.null(each)) {
    gg <- gg + facet_wrap(~ groups, scales = "free")
  }

  return(gg)
}



#' @title 2D Scatterplot of Samples
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#'
#' @export
gg_scatter <- function (samples, var1, var2, each = NULL,
                           keys = c("group", "Parameter"), highlight = NULL) {

  if (!is.null(highlight)) {
    aux_samples <- samples[1,]
    aux_samples[1, ] <- highlight
  }

  if (!is.null(each)) {
    samples <- gather.spmirt(samples, each, keys)
    aux_samples <- gather.spmirt(aux_samples, each, keys)
  }

  gg <- ggplot(samples, aes_(substitute(var1), substitute(var2))) +
    geom_point(alpha = 0.5) +
    geom_path(alpha = 0.4, linetype = 2)

  if (!is.null(highlight)) {
    gg <- gg + geom_point(data = aux_samples, col = 2, size = 2)
  }

  if (!is.null(each)) {
    gg <- gg + facet_wrap(~ groups, scales = "free")
  }

  return(gg)
}






#' @title Horizontal Errorbar Plot of Samples
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' summary(samples, "corr_chol") %>%
#'   mutate(param = corr_chol) %>%
#'   gg_errorbarh() +
#'   geom_point(aes(param, Parameters), col = 3)
#'
#' @export
gg_errorbarh <- function (df_summary, sorted = FALSE,
                          colors = c(rgb(1,0.5,0.1), "black"), ...) {

  if (sorted) {
    gg <- df_summary %>% ggplot(., aes(`50%`, `50%`))
  } else {
    gg <- df_summary %>% ggplot(., aes(`50%`, Parameters))
  }

  gg <- gg +
    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`, col = "95%"),
                   height = 0, ...) +
    geom_errorbarh(aes(xmin = `10%`, xmax = `90%`, col = "80%"), size = 2,
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
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' summary(samples, "corr_chol") %>%
#'   mutate(param = corr_chol) %>%
#'   gg_errorbar() +
#'   geom_point(aes(Parameters, param), col = 3)
#'
#' @export
gg_errorbar <- function (df_summary, colors = c(rgb(1,0.5,0.1), "black"), ...) {
  gg <- df_summary %>%
    ggplot(., aes(Parameters, `50%`)) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, col = "95%"),
                   width = 0, ...) +
    geom_errorbar(aes(ymin = `10%`, ymax = `90%`, col = "80%"), size = 2,
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
