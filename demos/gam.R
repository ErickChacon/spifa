
setwd("../src")
getwd()

# data
n <- 1000
# mu <- rep(seq(0, 10, length.out = n), n)
mu <- seq(-10, 8, length.out = n)
z <- rnorm(n)
Q <- tcrossprod(matrix(rnorm(n^2), n))
Q_L <- t(chol(Q))
A2 <- rbind(1, c(1:n))

# code
Rcpp::sourceCpp("gam.cpp")

# test
all.equal(solve_sympd(Q, t(A2)), solve(Q) %*% t(A2))
all.equal(solve_sympd_chol(Q_L, t(A2)), solve(Q) %*% t(A2))

A <- rbind(rep(1, n))
e <- rep(1, nrow(A))
y <- rmvnorm_rest_Q(mu, Q, A, e)
sum(y)
plot(y)
all.equal(sum(y), e)


# model
library(tidyverse)
library(mgcv)

Rcpp::sourceCpp("gam.cpp")

n <- 5000
data <- data.frame(x = runif(n, -10, 10))
# %>%
  # mutate(
  #        f = 0.1 * x ^ 2 + 4 * sin(x),
  #        y = f + rnorm(n, sd = 2)
  #        )

n_knots <- 20
sm <- smoothCon(s(x, k = n_knots, bs = "ps"), data = data, knots = NULL)[[1]]
X <- sm$X
D <- sm$D
tao <- 1
beta <- cumsum(cumsum(rnorm(n_knots)))
plot(beta, type = "l")
data$y <- X %*% cbind(beta) + rnorm(n , sd = 1)
ggplot(data, aes(x, y)) +
  geom_point()



bla <- gamcpp(data$y, X, D, 4, 1)
beta_est <- apply(bla, 1, median)
data$y_est <- X %*% cbind(beta_est)

ggplot(data, aes(x, y)) +
  geom_point() +
  geom_point(aes(x, y_est), col = 2)




