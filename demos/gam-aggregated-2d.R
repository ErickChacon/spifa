
rm(list = ls())
setwd("../src")
getwd()

# model
library(tidyverse)
library(mgcv)
library(Matrix)
library(raster)
library(RColorBrewer)
library(patchwork)


Rcpp::sourceCpp("gam.cpp")

# simulate betas igmrf
n_knots <- 10
tao <- 3
r1 <- raster(nrows = n_knots, ncols = n_knots, xmn = 0, xmx = n_knots, ymn = 0, ymx = n_knots)
ncell1 <- ncell(r1)
adj <- adjacent(r1, cells = 1:ncell1, 4, sorted = TRUE, pairs = TRUE)
adj_tab <- as.numeric(table(adj[, 1]))
Adj <- matrix(0, ncell1, ncell1)
Adj[adj] <- 1
Q <- (diag(adj_tab) - Adj) / tao ^ 2
r1 <- setValues(r1, rimvnorm_Q_eig(Q))
beta <- values(r1)

# simulate non-linear effect
sigma_ind <- 1
n_dim <- 100
n_ind <- n_dim
n <- n_dim ^ 2
data_ind <- expand.grid(x1 = seq(-10, 10, length.out = n_dim),
                    x2 = seq(-10, 10, length.out = n_dim))
S1 <- smoothCon(s(x1, k = n_knots, bs = "ps"), data = data_ind, knots = NULL)[[1]]
X1 <- S1$X
knots1 <- S1$knots[2 + 1:n_knots]
S2 <- smoothCon(s(x2, k = n_knots, bs = "ps"), data = data_ind, knots = NULL)[[1]]
X2 <- S2$X
knots2 <- S2$knots[2 + 1:n_knots]
X2_rev <- X2[, ncol(X2):1]
X_ind <- as.matrix(t(KhatriRao(t(X2_rev), t(X1))))
data_ind$f <- as.numeric(X_ind %*% cbind(beta))
data_ind$y <- data_ind$f + rnorm(n_ind ^ 2, sd = sigma_ind)

data_knots <- expand.grid(x1 = knots1, x2 = knots2) %>%
  arrange(desc(x2), x1) %>%
  mutate(beta = beta)

gg_rw1 <- ggplot(data_knots, aes(x1, x2)) +
  geom_tile(aes(fill = beta), col = 1) +
  geom_vline(xintercept = c(-10, 10)) +
  geom_hline(yintercept = c(-10, 10)) +
  scale_fill_distiller(palette = "Spectral") +
  coord_fixed()

gg_f <- ggplot(data_ind, aes(x1, x2)) +
  geom_tile(aes(fill = f)) +
  scale_fill_distiller(palette = "Spectral") +
  coord_fixed()

gg_y <- ggplot(data_ind, aes(x1, x2)) +
  geom_tile(aes(fill = y)) +
  scale_fill_distiller(palette = "Spectral") +
  coord_fixed()

# aggregated data
n_agg <- 5
breaks1 <- sort(c(runif(n_agg - 1, -10, 10), -10, 10))
breaks2 <- sort(c(runif(n_agg - 1, -10, 10), -10, 10))
lengths <- diff(breaks1)

data_ind$cut1 <- cut(data_ind$x1, breaks1, include.lowest = TRUE)
data_ind$cut2 <- cut(data_ind$x2, breaks2, include.lowest = TRUE)

data_agg <- data_ind %>%
  group_by(cut1, cut2) %>%
  summarise_all(mean) %>%
  tidyr::extract(cut1, c("low1", "high1"), "[(\\[](.+),(.+)[\\]]", remove = F, T) %>%
  tidyr::extract(cut2, c("low2", "high2"), "[(\\[](.+),(.+)[\\]]", remove = F, T)

gg_agg <- ggplot(data_agg) +
  geom_rect(aes(xmin = low1, xmax = high1, ymin = low2, ymax = high2, fill = y)) +
  scale_fill_distiller(palette = "Spectral") +
  coord_fixed()

# monte carlo approximation of the design matrix
n_mc <- 20000
data_mc <- data.frame(x1 = runif(n_mc, -10, 10), x2 = runif(n_mc, -10, 10))
data_mc$cut1 <- cut(data_mc$x1, breaks1, include.lowest = TRUE)
data_mc$cut2 <- cut(data_mc$x2, breaks2, include.lowest = TRUE)
S1_mc <- smoothCon(s(x1, k = n_knots, bs = "ps"), data = data_mc, knots = NULL)[[1]]
X1_mc <- S1_mc$X
# knots1 <- S1$knots[2 + 1:n_knots]
S2_mc <- smoothCon(s(x2, k = n_knots, bs = "ps"), data = data_mc, knots = NULL)[[1]]
X2_mc <- S2_mc$X
X2_mc_rev <- X2_mc[, ncol(X2_mc):1]
X_mc <- as.matrix(t(KhatriRao(t(X2_mc_rev), t(X1_mc))))

Q2 <- (diag(adj_tab) - Adj)
D_mc <- chol(Q2)
data_agg_mc <- data_mc %>%
  bind_cols(as.data.frame(X_mc)) %>%
  group_by(cut1, cut2, add = TRUE) %>%
  summarise_all(mean) %>%
  ungroup(cut1, cut2)

X_agg_mc <- dplyr::select(data_agg_mc, V1:V100) %>% as.matrix()

# run the aggregates model
bla <- gamcpp(data_agg$y, X_agg_mc, D_mc, 1, 1, niter = 1000)

beta_est <- apply(bla$beta, 2, median)
data_agg$y_est <- X_agg_mc %*% cbind(beta_est)
data_ind$f_est <- X_ind %*% cbind(beta_est)

gg_agg_est <- ggplot(data_agg) +
  geom_rect(aes(xmin = low1, xmax = high1, ymin = low2, ymax = high2, fill = y_est)) +
  scale_fill_distiller(palette = "Spectral") +
  coord_fixed()

gg_y_est <- ggplot(data_ind, aes(x1, x2)) +
  geom_tile(aes(fill = f_est)) +
  scale_fill_distiller(palette = "Spectral") +
  coord_fixed()

print(gg_rw1 + gg_f + gg_y + gg_agg + gg_agg_est + gg_y_est + plot_layout(2))


# pdf("gam_2d_aggregated.pdf", width = 7)
# print(gg_rw1 + gg_f + gg_y + gg_agg + gg_agg_est + gg_y_est + plot_layout(2))
# dev.off()





# # ggplot(data_ind, aes(x, y)) +
# #   geom_point() +
# #   geom_line(aes(x, f), col = 4, lwd = 1.5) +
# #   geom_line(aes(x, y_est), data_ind, col = 2, size = 1.5) +
# #   geom_point(aes(x, y), data_agg, col = 4, size = 3) +
# #   geom_errorbarh(aes(xmin = low, xmax = high, y = y), data_agg, lwd = 2) +
# #   geom_point(aes(x, y_est), data_agg, col = 2, size = 3)
#
# bla <- ggplot(data_ind, aes(x, y)) +
#   geom_point(alpha = 0.3) +
#   geom_line(aes(x, f), col = 4, lwd = 2, linetype = 2) +
#   geom_errorbarh(aes(xmin = low, xmax = high, y = y), data_agg, col = 4, lwd = 1.5) +
#   geom_line(aes(x, y_est), data_ind, col = 2, size = 2, linetype = 2) +
#   geom_errorbarh(aes(xmin = low, xmax = high, y = y_est), data_agg, col = 2, lwd = 1.5)
# bla
#
#
# # pdf(file.path("gam_aggregate.pdf"), width = 7)
# # print(bla)
# # dev.off()
#
#
#
#
# # hist(bla$sigma2 ^ 0.5, 100)
# # hist(bla$tau2 ^ 0.5, 100)
#
# # plot(bla$sigma2 ^ 0.5, type = "l")
# # plot(bla$tau2 ^ 0.5, type = "l")
#
#
#
#
#
#
#
#
