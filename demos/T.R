T <- rbind(c(1, 0, 1), c(0, 1, 0), c(0, 0, 1))
S <- matrix(0, 9, 9)
S[1:3,1:3] = 1
S[4:6,4:6] = 3
S[7:9,7:9] = 7

Rcpp::sourceCpp("../src/arma-mat.cpp")
TST(S, T)



X <- matrix(rnorm(9), 3)

T %*% X %*% t(T)
t(T) %*% X %*% T

