
rm(list = ls())
setwd("../src")
getwd()

# model
library(tidyverse)
library(mgcv)

Rcpp::sourceCpp("gam.cpp")

# simulate betas
tau <- 0.01
n_knots <- 50
beta <- cumsum(cumsum(rnorm(n_knots, tau)))
plot(beta)

# simulate non-linear effect
n_ind <- 1000
data_ind <- data.frame(x = seq(-10, 10, length.out = n_ind))
sm <- smoothCon(s(x, k = n_knots, bs = "ps"), data = data_ind, knots = NULL)[[1]]
X_ind <- sm$X
D_ind <- sm$D
sigma_ind <- 9
data_ind$f <- X_ind %*% cbind(beta)
data_ind$y <- data_ind$f + rnorm(n_ind, sd = sigma_ind)
ggplot(data_ind, aes(x, y)) +
  geom_point() +
  geom_line(aes(x, f), col = 2, lwd = 2)

bla <- gamcpp(data$y, X, D, 1, 1, niter = 10000)
beta_est <- apply(bla$beta, 2, median)
data$y_est <- X %*% cbind(beta_est)

ggplot(data, aes(x, y)) +
  geom_point() +
  geom_point(aes(x, y_est), col = 2)

hist(bla$sigma2 ^ 0.5, 100)
hist(bla$tau2 ^ 0.5, 100)

plot(bla$sigma2 ^ 0.5, type = "l")
plot(bla$tau2 ^ 0.5, type = "l")








