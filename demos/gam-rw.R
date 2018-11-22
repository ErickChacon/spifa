
rm(list = ls())
setwd("../src")
getwd()

# model
library(tidyverse)
library(mgcv)

Rcpp::sourceCpp("gam.cpp")

# n <- 5000
n <- 300
data <- data.frame(x = runif(n, -10, 10))
# %>%
#   mutate(
#          f = 0.1 * x ^ 2 + 4 * sin(x),
#          y = f + rnorm(n, sd = 2)
#          )

# rw1

n_knots <- 50
sm <- smoothCon(s(x, k = n_knots, bs = "ps", m = 1), data = data, knots = NULL)[[1]]
X <- sm$X
D <- sm$D
tao <- 0.01
sigma <- 9
# beta <- cumsum(cumsum(rnorm(n_knots, tao)))
beta1 <-cumsum(rnorm(n_knots, sd = tao))
plot(beta1, type = "l")

# tao <- 1
DtD <- crossprod(D)
Q <- DtD / tao ^ 2 + 1
# Matrix::rankMatrix(DtD)
# Matrix::rankMatrix(Q)

rw1 <- rmvnorm_Q(rep(0, n_knots), Q)
rw1 <- rw1 - mean(rw1)
plot(rw1, type = "b", col = 2, lwd = 2, pch = 19)



# space

library(raster)
library(RColorBrewer)
tao <- 1
n <- 60
r1 <- raster(nrows = n, ncols = n, xmn = 0, xmx = n, ymn = 0, ymx = n)
ncell1 <- ncell(r1)
adj <- adjacent(r1, cells = 1:(n ^ 2), 4, sorted = TRUE, pairs = TRUE)
adj_tab <- as.numeric(table(adj[, 1]))
Adj <- matrix(0, ncell1, ncell1)
Adj[adj] <- 1
Q <- (diag(adj_tab) - Adj) / tao ^ 2 + 1
rw1 <- rmvnorm_Q(rep(0, ncell1), Q)
rw1 <- rw1 - mean(rw1)
r1 <- setValues(r1, rw1)
plot(r1, col = colorRampPalette(brewer.pal(11,"Spectral"))(30))
# plot(r1, col = colorRampPalette(brewer.pal(9, "Greys"))(50))
plot(rasterToPolygons(r1), add = TRUE, border = 'black', lwd = 1)
# plot(raster(Q))




# lattice::levelplot(D)
#
# data$y <- X %*% cbind(beta) + rnorm(n, sd = sigma)
# ggplot(data, aes(x, y)) +
#   geom_point()
#
# bla <- gamcpp(data$y, X, D, 1, 1, niter = 10000)
# beta_est <- apply(bla$beta, 2, median)
# data$y_est <- X %*% cbind(beta_est)
#
# ggplot(data, aes(x, y)) +
#   geom_point() +
#   geom_point(aes(x, y_est), col = 2)
#
# # hist(bla$sigma2 ^ 0.5, 100)
# # hist(bla$tau2 ^ 0.5, 100)
# #
# # plot(bla$sigma2 ^ 0.5, type = "l")
# # plot(bla$tau2 ^ 0.5, type = "l")
# #
# #
# #
# #
# #
# #
# #
# #
