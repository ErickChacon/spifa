# #' @param threshold Threshold to plot the extreme events.
StatEvents <- ggplot2::ggproto("StatEvents", ggplot2::Stat,
  required_aes = c("coord_x", "coord_y", "value"),
  compute_group = function(data, scales, cutoff = 1, width = 0.1) {
    olddata <- data
    data <- gstat::variogram(value ~ 1, ~ coord_x + coord_y, data,
                             cutoff = cutoff, width = width) %>%
      subset(select = c(dist, gamma, np))

    data <- data.frame(x = c(data$x, transition$x),
                       y = c(data$y, transition$y),
                       ymin = threshold,
                       ymax = c(data$y, transition$y))
  }
)

stat_events <- function(mapping = NULL, data = NULL, geom = "ribbon",
                        position = "identity", na.rm = FALSE, show.legend = NA,
                        inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatEvents, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


vg <- data_geo %>%
  mutate(s2 = s1) %>%
  dplyr::filter(response_label == 1) %>%
  gstat::variogram(mgp.list.mean ~ 1, ~ s1 + s2, . , cutoff = 3, width = 0.005)
ggplot(vg, aes(dist, gamma)) +
  geom_point(aes(size = np)) +
  geom_smooth() +
  expand_limits(y = 0, x = 0) +
  scale_x_continuous(limits = c(0, 3)) +
  stat_function(fun = function(x) 0.6 * (1-exp(- x/0.05)), col = 2, size = 2)

