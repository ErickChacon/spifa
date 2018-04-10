
getwd()
Rcpp::sourceCpp("../src/pdf.cpp")

Corr <- matrix(c(1, -0.9, 0, -0.9, 1, 0, 0, 0, 1), nrow = 3)
sigmas <- c(1,1,3)
D <- diag(sigmas)
Cov <- D %*% Corr %*% D
beta <- c(-10, 0, 10)
L <- t(chol(Cov))
L_inv <- solve(L)

X <- mvtnorm::rmvnorm(n = 4, mean = beta, sigma = Cov)

dmvnorm_cholinv(t(X), matrix(rep(beta, 4), 3), L_inv)
dmvnorm_chol(t(X), matrix(rep(beta, 4), 3), L)
dmvnorm(t(X), matrix(rep(beta, 4), 3), Cov)
sum(mvtnorm::dmvnorm((X), beta, sigma = Cov, log = TRUE))

