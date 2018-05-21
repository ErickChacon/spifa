# TODO: include this on tests
# T <- rbind(c(1, 0, 1), c(0, 1, 0), c(0, 0, 1))
T <- rbind(c(0.5, 0, 0.2), c(0, 0.8, 0), c(0, 0, 0.7))
TK <- kronecker(T, diag(3))
S <- matrix(0, 9, 9)
S[1:3,1:3] = 1
S[4:6,4:6] = 3
S[7:9,7:9] = 7

Rcpp::sourceCpp("../src/arma-mat.cpp")
TST(S, T) - TK %*% S %*% t(TK)
mean(TST(S, T) - TK %*% S %*% t(TK))

