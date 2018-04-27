T <- cbind(c(1, 1, 0), c(0, 1, 0), c(0, 0, 0))
X <- matrix(rnorm(9), 3)

T %*% X %*% t(T)
t(T) %*% X %*% T
