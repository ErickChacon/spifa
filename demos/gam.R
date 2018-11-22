
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

n_knots <- 50
sm <- smoothCon(s(x, k = n_knots, bs = "ps"), data = data, knots = NULL)[[1]]
X <- sm$X
D <- sm$D
tao <- 0.001
sigma <- 9
beta <- cumsum(cumsum(rnorm(n_knots, tao)))
plot(beta, type = "l")

data$y <- X %*% cbind(beta) + rnorm(n, sd = sigma)
ggplot(data, aes(x, y)) +
  geom_point()

bla <- gamcpp(data$y, X, D, 1, 1, niter = 10000)
beta_est <- apply(bla$beta, 2, median)
data$y_est <- X %*% cbind(beta_est)

ggplot(data, aes(x, y)) +
  geom_point() +
  geom_point(aes(x, y_est), col = 2)

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
