
getwd()
Rcpp::sourceCpp("../src/pdf.cpp")

Corr <- matrix(c(1, -0.9, 0, -0.9, 1, 0, 0, 0, 1), nrow = 3)
sigmas <- c(1,1,3)
D <- diag(sigmas)
(Cov <- D %*% Corr %*% D)
(Cov <- t(sigmas * t(sigmas * Corr)))
test1(sigmas, Corr)

beta <- c(-10, 0, 10)
L <- t(chol(Cov))
L_inv <- solve(L)

X <- mvtnorm::rmvnorm(n = 4, mean = beta, sigma = Cov)

dmvnorm_cholinv(t(X), matrix(rep(beta, 4), 3), L_inv)
dmvnorm_chol(t(X), matrix(rep(beta, 4), 3), L)
dmvnorm(t(X), matrix(rep(beta, 4), 3), Cov)
sum(mvtnorm::dmvnorm((X), beta, sigma = Cov, log = TRUE))

k <- 5
B <- matrix(1, k, k)
(A <- matrix(rnorm(k), k))
test1(A, B)

D %*% matrix(1, 3, 3)
test1(sigmas, matrix(1, 3, 3))

