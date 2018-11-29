
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

n_dim <- 50
n <- n_dim ^ 2
data <- expand.grid(x1 = seq(-10, 10, length.out = n_dim),
                    x2 = seq(-10, 10, length.out = n_dim))
n_knots <- 10
S1 <- smoothCon(s(x1, k = n_knots, bs = "ps"), data = data, knots = NULL)[[1]]
X1 <- S1$X
knots1 <- S1$knots[2 + 1:n_knots]
S2 <- smoothCon(s(x2, k = n_knots, bs = "ps"), data = data, knots = NULL)[[1]]
X2 <- S2$X
knots2 <- S2$knots[2 + 1:n_knots]
X2_rev <- X2[, ncol(X2):1]
X <- as.matrix(t(KhatriRao(t(X2_rev), t(X1))))

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
data_knots <- expand.grid(x1 = knots1, x2 = knots2) %>%
  arrange(desc(x2), x1) %>%
  mutate(beta = beta)

sigma <- 3
data$f <- as.numeric(X %*% cbind(beta))
data$y <- data$f + rnorm(n, sd = sigma)

gg_rw1 <- ggplot(data_knots, aes(x1, x2)) +
  geom_tile(aes(fill = beta), col = 1) +
  geom_vline(xintercept = c(-10, 10)) +
  geom_hline(yintercept = c(-10, 10)) +
  scale_fill_distiller(palette = "Spectral") +
  coord_fixed()

gg_f <- ggplot(data, aes(x1, x2)) +
  geom_tile(aes(fill = f), col = 1) +
  scale_fill_distiller(palette = "Spectral") +
  coord_fixed()

gg_y <- ggplot(data, aes(x1, x2)) +
  geom_tile(aes(fill = y), col = 1) +
  scale_fill_distiller(palette = "Spectral") +
  coord_fixed()


Q2 <- (diag(adj_tab) - Adj)
D <- chol(Q2)

bla <- gamcpp(data$y, X, D, 1, 1, niter = 1000)

beta_est <- apply(bla$beta, 2, median)
# beta_est <- rep(1, length(beta))
data$y_est <- X %*% cbind(beta_est)

gg_est <- ggplot(data, aes(x1, x2)) +
  geom_tile(aes(fill = y_est), col = 1) +
  scale_fill_distiller(palette = "Spectral") +
  coord_fixed()

# pdf("gam_2d.pdf", width = 7)
print(gg_rw1 + gg_f + gg_y + gg_est)
# dev.off()

summary(bla$sigma2 ^ 0.5)
summary(bla$tau2 ^ 0.5)

# hist(bla$sigma2 ^ 0.5, 100)
# hist(bla$tau2 ^ 0.5, 100)
#
# plot(bla$sigma2 ^ 0.5, type = "l")
# plot(bla$tau2 ^ 0.5, type = "l")
#
#
#
#
#
#
#
#
