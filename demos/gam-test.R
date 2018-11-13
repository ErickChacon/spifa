

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
